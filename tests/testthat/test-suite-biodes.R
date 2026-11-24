# Biodes -> Simitall -> Bioana hand-off. Fixtures are real Biodes
# (github.com/nirwan1265/BioDesignGPT) design packages, handoff_version 0.4.0.

copy_fixture <- function(name) {
  src <- test_path("fixtures", "biodes", name)
  dst <- file.path(tempfile("biodes_"), name)
  dir.create(dst, recursive = TRUE)
  file.copy(list.files(src, full.names = TRUE), dst)
  dst
}

test_that("Simitall simulates an RCBD with subsamples, power and attrition from the spec", {
  d <- copy_fixture("rcbd_subsamples")
  suppressMessages(simitall_from_biodes(d, seed = 7))
  meta <- read.delim(file.path(d, "sample_metadata_template.tsv"))
  truth <- jsonlite::read_json(file.path(d, "sim_truth.json"))
  expect_equal(nrow(meta), 24 * 3)
  expect_true(all(c("unit_id", "subsample") %in% names(meta)))
  expect_equal(truth$effect, 2)
  expect_equal(truth$residual_sd, 3)
  expect_equal(length(truth$dropped_units), 1)            # 5% of 24 units
  expect_true(all(is.na(meta$primary_outcome[meta$unit_id %in% unlist(truth$dropped_units)])))
  expect_equal(unlist(truth$treatment_effects$treatment)[["G6"]], 2)
  expect_equal(truth$model, "y ~ treatment + block")

  # same seed, same data
  d2 <- copy_fixture("rcbd_subsamples")
  suppressMessages(simitall_from_biodes(d2, seed = 7))
  expect_identical(read.delim(file.path(d2, "sample_metadata_template.tsv")), meta)
})

test_that("Bioana fits the spec's model, averages subsamples and refuses them as replicates", {
  skip_if_not_installed("bioana")
  d <- copy_fixture("rcbd_subsamples")
  suppressMessages(simitall_from_biodes(d, seed = 3))
  r <- suppressMessages(bioana::bioana_analyze(d, quiet = TRUE, allow_fallback = TRUE))
  a <- r$summary$analyses$outcome
  expect_match(a$method, "^lm")
  expect_equal(a$n_units_fitted, 23)              # 24 units, 1 lost to attrition
  expect_equal(a$subsamples, "aggregate")
  expect_true(file.exists(file.path(d, "bioana", "outcome", "results.tsv")))
  expect_true(file.exists(file.path(d, "bioana", "outcome", "recovery.tsv")))
  expect_equal(a$recovery$n_terms, 5)
  bad <- suppressMessages(bioana::bioana_analyze(d, params = list(outcome = list(subsamples = "as_replicates")), quiet = TRUE, allow_fallback = TRUE))
  expect_match(bad$summary$not_analysed$outcome$note, "not independent replicates")
  n <- suppressMessages(bioana::bioana_analyze(d, params = list(outcome = list(subsamples = "nest")), quiet = TRUE, allow_fallback = TRUE))
  expect_match(n$summary$analyses$outcome$model_fitted, "unit_id_")
})

test_that("mixed models: split plot and alpha lattice use lmer with the spec's random terms", {
  skip_if_not_installed("bioana")
  for (nm in c("split_field", "alpha_120")) {
    d <- copy_fixture(nm)
    suppressMessages(simitall_from_biodes(d))
    r <- suppressMessages(bioana::bioana_analyze(d, quiet = TRUE, allow_fallback = TRUE))
    expect_match(r$summary$analyses$outcome$method, "lmerMod", info = nm)
    expect_true(r$summary$analyses$outcome$recovery$n_terms > 0, info = nm)
    truth <- jsonlite::read_json(file.path(d, "sim_truth.json"))
    expect_true(length(truth$random_sd) == 1, info = nm)
    expect_true(length(truth$assumptions) > 0, info = nm)  # spec has no power: defaults are declared
  }
})

test_that("count assays get a count matrix with Biodes sample ids and a DE test", {
  d <- copy_fixture("crd_rnaseq")
  suppressMessages(simitall_from_biodes(d, n_genes = 200))
  counts <- read.delim(file.path(d, "counts.tsv"), check.names = FALSE)
  meta <- read.delim(file.path(d, "sample_metadata_template.tsv"))
  expect_equal(names(counts), c("gene_id", meta$sample_id))
  expect_equal(nrow(counts), 200)
  skip_if_not_installed("bioana")
  r <- suppressMessages(bioana::bioana_analyze(d, quiet = TRUE, allow_fallback = TRUE))
  expect_true(file.exists(file.path(d, "bioana", "rnaseq", "results.tsv")))
  rec <- r$summary$analyses$rnaseq$recovery$treatmentdrought
  expect_equal(rec$true_changes, 20)
  expect_gt(rec$estimate_cor, 0.5)
})

test_that("unsupported models, wrong versions and unknown columns stop with a clear message", {
  d <- copy_fixture("rcbd_greenhouse")
  s <- jsonlite::read_json(file.path(d, "study_spec.json"))
  s$analysis$model <- "y ~ genotype + f(row, col) + (1 | row) + (1 | col)"
  jsonlite::write_json(s, file.path(d, "study_spec.json"), auto_unbox = TRUE)
  expect_error(simitall_from_biodes(d), "cannot simulate or fit yet")
  s$analysis$model <- "y ~ treatment + block"
  s$handoff_version <- "9.9.0"
  jsonlite::write_json(s, file.path(d, "study_spec.json"), auto_unbox = TRUE)
  expect_error(simitall_from_biodes(d), "handoff_version")
  s$handoff_version <- "0.4.0"
  s$analysis$model <- "y ~ treatment + greenhouse_door"
  jsonlite::write_json(s, file.path(d, "study_spec.json"), auto_unbox = TRUE)
  expect_error(simitall_from_biodes(d), "no matching column")
  e <- file.path(tempdir(), "empty_dir_for_test"); dir.create(e, showWarnings = FALSE)
  expect_error(simitall_from_biodes(e), "Not a Biodes design folder")
})

test_that("analysing before simulating says to run simitall_from_biodes() first", {
  skip_if_not_installed("bioana")
  d <- copy_fixture("rcbd_greenhouse")
  expect_error(bioana::bioana_analyze(d), "simitall_from_biodes")
})

test_that("layouts that do not match the plan are flagged", {
  d <- copy_fixture("rcbd_greenhouse")
  spec <- jsonlite::read_json(file.path(d, "study_spec.json"))
  lay <- utils::read.csv(file.path(d, "field_layout.csv"))
  # every spec factor must be a layout column (Biodes layout_from_spec keeps them all)
  spec$factors <- list(list(name = "genotype", n_levels = 6), list(name = "water", n_levels = 3))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE)
  iss <- .suite_layout_issues(list(), d)
  expect_true(any(grepl("factor 'water' is in the spec but not a column", iss)))
  # a factor column holding joined levels of several factors is not that factor
  lay$genotype <- paste0(lay$treatment, ":wet"); lay$water <- "wet"
  utils::write.csv(lay, file.path(d, "field_layout.csv"), row.names = FALSE)
  expect_true(any(grepl("holds joined levels", .suite_layout_issues(list(), d))))
  lay$genotype <- lay$treatment
  utils::write.csv(lay, file.path(d, "field_layout.csv"), row.names = FALSE)
  expect_length(.suite_layout_issues(list(), d), 0)
  expect_equal(.suite_slug("Greenhouse drought on 6 sorghum genotypes, 4 benches"), "greenhouse_drought_sorghum_genotypes_benches")
})

test_that("simitall_from_question explains how to get Biodes when it is missing", {
  skip_if(requireNamespace("studyforge", quietly = TRUE), "Biodes is available")
  expect_error(simitall_from_question("anything"), "needs Biodes")
})

test_that("genotype levels sort naturally (G2 before G10)", {
  expect_equal(.suite_natural_sort(c("G10", "G2", "G1", "G11")), c("G1", "G2", "G10", "G11"))
  expect_equal(.suite_natural_sort(c("10", "9", "100")), c("9", "10", "100"))
  expect_equal(.suite_natural_sort(c("wet", "dry")), c("dry", "wet"))
})
