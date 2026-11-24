# Biodes 6e49dab+: generic model terms (treatment, treatment2, block, rep, row,
# col, unit, subject, time, env, period, sequence, check) and analysis.terms in
# study_spec.json. Fixtures in fixtures/biodes_terms are real Biodes packages;
# catalog_models.tsv lists all design-catalog models at that commit.

copy_terms_fixture <- function(name) {
  src <- test_path("fixtures", "biodes_terms", name)
  dst <- file.path(tempfile("biodes_"), name)
  dir.create(dst, recursive = TRUE)
  file.copy(list.files(src, full.names = TRUE), dst)
  dst
}

test_that("analysis.terms picks the layout columns; the split plot's cell label is not the main-plot factor", {
  skip_if_not_installed("lme4")
  rules <- simitall:::.suite_rules("handoff")
  for (name in c("rcbd_greenhouse", "alpha_120", "split_field", "crd_rnaseq")) {
    d <- copy_terms_fixture(name)
    spec <- jsonlite::read_json(file.path(d, "study_spec.json"))
    layout <- read.csv(file.path(d, "field_layout.csv"), colClasses = "character")
    p <- simitall:::.suite_parse_model(spec$analysis$model, spec, layout, rules)
    expect_equal(p$col[["treatment"]], if (name == "split_field") "main" else "treatment", info = name)
  }
  d <- copy_terms_fixture("split_field")
  spec <- jsonlite::read_json(file.path(d, "study_spec.json"))
  layout <- read.csv(file.path(d, "field_layout.csv"), colClasses = "character")
  p <- simitall:::.suite_parse_model(spec$analysis$model, spec, layout, rules)
  expect_equal(unname(p$col[c("treatment", "treatment2", "block")]), c("main", "sub", "block"))
  # an explicit analysis.terms.columns map wins
  spec$analysis$terms$columns <- list(treatment = "sub", treatment2 = "main")
  p2 <- simitall:::.suite_parse_model(spec$analysis$model, spec, layout, rules)
  expect_equal(unname(p2$col[c("treatment", "treatment2")]), c("sub", "main"))
  # two terms on one column is refused
  spec$analysis$terms$columns <- list(treatment = "main", treatment2 = "main")
  expect_error(simitall:::.suite_parse_model(spec$analysis$model, spec, layout, rules), "map to layout column 'main'")
})

test_that("new Biodes packages simulate and fit end to end", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("bioana")
  for (name in c("rcbd_greenhouse", "alpha_120", "split_field")) {
    d <- copy_terms_fixture(name)
    suppressMessages(simitall_from_biodes(d, seed = 11, effect = 4, residual_sd = 1))
    r <- suppressMessages(bioana::bioana_analyze(d, quiet = TRUE, allow_fallback = TRUE))
    expect_false(is.null(r$summary$analyses$outcome$recovery), info = name)
    tab <- r$recovery$outcome
    expect_true(all(abs(tab$error) < stats::qt(0.9995, tab$df) * tab$se), info = name)
    if (name == "split_field") {
      expect_true(all(c("treatmentwet", "treatment2H4", "treatmentwet:treatment2H4") %in% r$results$outcome$term))
    }
  }
  d <- copy_terms_fixture("crd_rnaseq")
  suppressMessages(simitall_from_biodes(d, seed = 5, n_genes = 300))
  r <- suppressMessages(bioana::bioana_analyze(d, quiet = TRUE, allow_fallback = TRUE))
  expect_true(file.exists(file.path(d, "bioana", "rnaseq", "results.tsv")))
})

# Balanced synthetic layout for a generic-term model: between-unit factors
# crossed, units (and subjects within units) nested in the cells, time and
# period (or two observations) within each unit, row and col assigned at random.
synthetic_layout <- function(model, seed = 1) {
  set.seed(seed)
  v <- setdiff(all.vars(stats::as.formula(model)), "y")
  quad <- grepl("^2", model, fixed = TRUE)
  lev <- list(treatment = if (quad) c("-1", "0", "1") else paste0("T", 1:4),
              treatment2 = if (quad) c("-1", "0", "1") else paste0("S", 1:3),
              check = c("entry", "check"), sequence = c("AB", "BA"), env = c("E1", "E2"),
              rep = c("1", "2", "3"), block = c("1", "2", "3"))
  grid <- expand.grid(lev[intersect(names(lev), v)], stringsAsFactors = FALSE)
  if (!ncol(grid)) grid <- data.frame(cell = 1)
  grid$cell_ <- seq_len(nrow(grid))
  units <- merge(grid, data.frame(k_ = 1:2))
  units$unit <- paste0("U", units$cell_, "_", units$k_)
  if ("subject" %in% v) {
    units <- merge(units, data.frame(s_ = 1:(if ("unit" %in% v) 2 else 1)))
    units$subject <- paste0(units$unit, "_", units$s_)
  }
  within <- list(time = c("t1", "t2", "t3"), period = c("p1", "p2"))[intersect(c("time", "period"), v)]
  # nothing measured within the finest unit: two observations each, so (1 | unit) is not confounded with the residual
  rows <- if (length(within)) merge(units, expand.grid(within, stringsAsFactors = FALSE)) else merge(units, data.frame(obs_ = 1:2))
  if ("row" %in% v) rows$row <- as.character(sample(rep_len(1:6, nrow(rows))))
  if ("col" %in% v) rows$col <- as.character(sample(rep_len(1:6, nrow(rows))))
  rows <- rows[, intersect(c(names(lev), "unit", "subject", "time", "period", "row", "col"), names(rows)), drop = FALSE]
  rows[] <- lapply(rows, as.character)
  rows
}

write_synthetic_package <- function(id, model, dir) {
  layout <- synthetic_layout(model)
  dir.create(dir, recursive = TRUE)
  utils::write.csv(layout, file.path(dir, "field_layout.csv"), row.names = FALSE)
  meta <- cbind(layout, sample_id = sprintf("S%04d", seq_len(nrow(layout))), primary_outcome = "")
  utils::write.table(meta, file.path(dir, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  spec <- list(handoff_version = "0.4.0", design = list(chosen = id, randomization_seed = 1),
               primary_outcome = list(name = "y"), power = list(minimum_detectable_difference = 4, residual_sd = 1),
               analysis = list(model = model, terms = list(layout_columns = names(layout))))
  jsonlite::write_json(spec, file.path(dir, "study_spec.json"), auto_unbox = TRUE)
  dir
}

test_that("every Biodes catalog model simulates and fits from its generic terms", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("bioana")
  cat_models <- read.delim(test_path("fixtures", "biodes_terms", "catalog_models.tsv"), stringsAsFactors = FALSE)
  expect_equal(nrow(cat_models), 74)
  root <- tempfile("catalog_")
  failed <- character()
  for (i in seq_len(nrow(cat_models))) {
    id <- cat_models$id[i]
    d <- write_synthetic_package(id, cat_models$model[i], file.path(root, id))
    res <- tryCatch({
      suppressMessages(suppressWarnings(simitall_from_biodes(d)))
      suppressMessages(suppressWarnings(bioana::bioana_analyze(d, quiet = TRUE, allow_fallback = TRUE)))
    }, error = function(e) e)
    if (!inherits(res, "error") && !is.null(res$summary$not_analysed$outcome)) res <- simpleError(res$summary$not_analysed$outcome$note)
    if (inherits(res, "error")) {
      failed <- c(failed, paste0(id, ": ", conditionMessage(res)))
      next
    }
    expect_true(all(is.finite(res$results$outcome$estimate)), info = id)
    tab <- res$recovery$outcome
    # one data set per model: no estimate outside its 99.9% t interval (coverage across seeds is checked separately)
    if (!is.null(tab)) expect_true(all(abs(tab$error) < stats::qt(0.9995, tab$df) * tab$se), info = id)
  }
  expect_equal(failed, character())
})
