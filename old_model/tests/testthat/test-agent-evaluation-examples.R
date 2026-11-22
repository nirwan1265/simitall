test_that("evaluation examples are well formed and consistent with the analysis catalog", {
  path <- file.path(test_path("..", ".."), "inst", "agent", "prompt_eval", "study_spec_examples_draft.jsonl")
  if (!file.exists(path)) path <- system.file("agent", "prompt_eval", "study_spec_examples_draft.jsonl", package = "simitall")
  skip_if(!nzchar(path) || !file.exists(path), "evaluation examples not found")
  ex <- .simitall_eval_examples(path)
  expect_gte(length(ex), 50)

  ids <- vapply(ex, `[[`, "", "id")
  expect_false(anyDuplicated(ids) > 0)

  cat_path <- file.path(test_path("..", ".."), "inst", "agent", "capabilities", "analysis_catalog.json")
  if (!file.exists(cat_path)) cat_path <- system.file("agent", "capabilities", "analysis_catalog.json", package = "simitall")
  catalog <- jsonlite::fromJSON(cat_path, simplifyVector = FALSE)
  subs <- unlist(lapply(catalog$families, function(f) f$subanalyses), recursive = FALSE)
  sub_ids <- vapply(subs, `[[`, "", "id")
  enums <- list()
  for (s in subs) for (p in s$parameters %||% list()) {
    if (identical(p$type, "enum") && length(p$allowed_values)) enums[[paste(s$id, p$id)]] <- unlist(p$allowed_values)
  }
  valid_workflow <- c(sub_ids, "combination", "trap", "unsupported", "genome_simulation")

  for (e in ex) {
    s <- e$expected_spec
    expect_true(all(c("id", "split", "prompt", "expected_spec", "workflow") %in% names(e)), info = e$id)
    expect_true(e$split %in% c("dev", "test"), info = e$id)
    expect_true(e$workflow %in% valid_workflow, info = e$id)
    if (!is.null(s$goal)) expect_true(s$goal %in% c("simulate", "analyze", "simulate_and_analyze", "plan"), info = e$id)
    if (!is.null(s$data_source)) expect_true(s$data_source %in% c("ask", "provided", "synthetic", "package_demo"), info = e$id)
    if (!is.null(s$primary_analysis)) expect_true(s$primary_analysis %in% sub_ids, info = e$id)
    expect_true(all(unlist(s$secondary_analyses) %in% sub_ids), info = e$id)
    expect_true(all(unlist(e$must_not) %in% sub_ids), info = e$id)
    # enum parameters must use catalog values
    for (p in names(s$parameters %||% list())) {
      key <- paste(s$primary_analysis %||% "", p)
      if (!is.null(enums[[key]])) {
        expect_true(all(unlist(s$parameters[[p]]) %in% enums[[key]]), info = paste(e$id, p))
      }
    }
    # a must-ask field cannot also be an expected value (that would leak the answer)
    for (p in unlist(e$must_ask)) expect_null(.simitall_eval_get(s, p), info = paste(e$id, p))
  }
})

test_that("the evaluation view merges legacy and catalog specs", {
  answer <- list(
    spec = list(task = "gwas", organism = "maize", trait = list(type = "binary", heritability = 0.6),
                n_individuals = 300, defaulted = list("trait.heritability"), supported = TRUE),
    catalog_study_spec = NULL
  )
  v <- .simitall_eval_view(answer)
  expect_equal(v$primary_analysis, "binary_gwas")
  expect_equal(v$organism, "maize")
  expect_equal(v$design$population_size, 300)
  expect_true("parameters.heritability" %in% v$defaulted)
  ex <- list(id = "x", split = "test", workflow = "binary_gwas",
             expected_spec = list(primary_analysis = "binary_gwas", organism = "maize"),
             must_ask = list("parameters.heritability"), must_not = list("bulk_counts_simulation"))
  sc <- .simitall_eval_score_one(ex, v)
  expect_true(sc$primary_ok)
  expect_equal(sc$reading_hit, 1)
  expect_equal(sc$defaulted_not_asked, 1)
  expect_equal(sc$invented, 0)
})
