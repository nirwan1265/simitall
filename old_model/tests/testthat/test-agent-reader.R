# The trained model needs Python + torch; these tests check routing with the lexicon alone.
local_reader_options <- function(env = parent.frame()) {
  old <- options(simitall.reader_model = FALSE)
  do.call(on.exit, list(bquote(options(.(old))), add = TRUE), envir = env)
}

test_that("legacy reader keeps routing and attaches a lexicon reading for comparison", {
  local_reader_options()
  q <- "Simulate a maize NAM population with 25 founders and 200 RILs per family, then run a GWAS for flowering time."
  a <- simitall_ask(q, provider = "local")
  b <- simitall_ask(q, provider = "local", reader = "legacy")
  expect_equal(a$reader, "legacy")
  expect_equal(a$plan_source, b$plan_source)
  expect_equal(a$answer, b$answer)
  expect_equal(a$reading$reader, "lexicon")
  expect_null(a$reader_route)
})

test_that("new reader stops unsupported, impossible and clinical requests", {
  local_reader_options()
  for (q in c("Ribo-seq: estimate translation efficiency changes between drought and control in rice.",
              "selfing in dogs 5 generations",
              "Simulate a trait with heritability 1.4 on 200 lines.",
              "Use my patient's RNA-seq to decide which cancer drug to give.")) {
    a <- simitall_ask(q, provider = "local", reader = "new", include_code = TRUE)
    expect_equal(a$plan_source, "reader", info = q)
    expect_match(a$answer, "was not run", info = q)
    expect_false(grepl("```", a$answer), info = q)
    expect_true(.simitall_eval_view(a)$do_not_run, info = q)
  }
})

test_that("new reader asks which analysis for an ambiguous assay and honours the chosen goal", {
  local_reader_options()
  a <- simitall_ask("RNA-seq.", provider = "local", reader = "new")
  expect_equal(a$plan_source, "catalog")
  expect_equal(a$catalog_resolution$status, "needs_goal")
  ids <- vapply(a$catalog_resolution$family$subanalyses, `[[`, "", "id")
  expect_true(all(c("bulk_counts_simulation", "differential_expression") %in% ids))
  label <- Filter(function(s) s$id == "differential_expression", a$catalog_resolution$family$subanalyses)[[1]]$label
  b <- simitall_ask(paste("RNA-seq. Requested goal:", label), provider = "local", reader = "new")
  expect_equal(b$catalog_study_spec$primary_analysis, "differential_expression")
})

test_that("new reader lets the compiled spec run only when it agrees", {
  local_reader_options()
  q <- "Simulate a maize NAM population with 25 founders and 200 RILs per family, then run a GWAS for flowering time."
  a <- simitall_ask(q, provider = "local", reader = "new")
  expect_equal(a$reader_route, "spec")
  expect_equal(a$plan_source, "spec")
  d <- simitall_ask("Run DESeq2 on my maize counts, drought vs control, then GO enrichment of the up-regulated genes.",
                    provider = "local", reader = "new")
  expect_equal(d$plan_source, "catalog")
  expect_equal(d$catalog_study_spec$primary_analysis, "differential_expression")
  expect_equal(d$catalog_study_spec$organism, "maize")
  expect_true("gene_coordinate_annotation" %in% d$catalog_study_spec$secondary_analyses)
})

test_that("the reader option sets the default", {
  local_reader_options()
  old <- options(simitall.reader = "new"); on.exit(options(old), add = TRUE)
  expect_equal(simitall_ask("selfing in dogs 5 generations", provider = "local")$plan_source, "reader")
})

test_that("typos, reversed phrasing, dropped assays and F1/selfing contradictions are caught", {
  local_reader_options()
  q <- "selfing in poplus tree for 5 generations and then rnaseq on the F1 siblings between 2 random parents"
  a <- simitall_ask(q, provider = "local", reader = "new")
  expect_equal(a$plan_source, "reader")
  expect_equal(a$reading$organism, "poplar")
  expect_equal(a$reading$design$selfing_generations, 5L)
  expect_true(any(grepl("selfing is not possible in poplar", a$reading$flags)))
  expect_true(any(grepl("assumed 'poplus'", a$reading$flags)))
  expect_true(any(grepl("RNA-seq was mentioned without saying", a$reading$flags)))
  expect_true(any(grepl("F1 siblings come from a single cross", a$reading$flags)))

  expect_equal(.simitall_detect_organism("GWAS in sorgum for plant height")$organism, "sorghum")
  expect_equal(.simitall_detect_organism("simulate drosphila sib mating")$organism, "drosophila")
  expect_true(is.na(.simitall_detect_organism("make a new avenue of research")$organism))
  expect_true(is.na(.simitall_detect_organism("please simulate the population")$organism))
})
