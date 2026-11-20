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
  d <- copy_fixture("rcbd_subsamples")
  suppressMessages(simitall_from_biodes(d, seed = 3))
  r <- suppressMessages(bioana_fit_biodes(d))
  expect_equal(r$summary$method, "lm")
  expect_equal(r$summary$n_units_fitted, 23)              # 24 units, 1 lost to attrition
  expect_equal(r$summary$subsamples, "aggregate")
  expect_true(file.exists(file.path(d, "results.tsv")))
  expect_true(file.exists(file.path(d, "recovery.tsv")))
  expect_equal(r$recovery$summary$n_terms, 5)
  expect_error(bioana_fit_biodes(d, subsamples = "as_replicates"), "not independent replicates")
  n <- suppressMessages(bioana_fit_biodes(d, subsamples = "nest"))
  expect_match(n$summary$model_fitted, "unit_id_")
})

test_that("mixed models: split plot and alpha lattice use lmer with the spec's random terms", {
  for (nm in c("split_field", "alpha_120")) {
    d <- copy_fixture(nm)
    suppressMessages(simitall_from_biodes(d))
    r <- suppressMessages(bioana_fit_biodes(d))
    expect_match(r$summary$method, "lmerMod", info = nm)
    expect_true(r$recovery$summary$n_terms > 0, info = nm)
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
  r <- suppressMessages(bioana_fit_biodes(d))
  expect_true(file.exists(file.path(d, "de_results.tsv")))
  rec <- r$summary$counts$recovery$drought
  expect_equal(rec$true_de, 20)
  expect_gt(rec$log2fc_cor, 0.5)
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
  d <- copy_fixture("rcbd_greenhouse")
  expect_error(bioana_fit_biodes(d), "Run simitall_from_biodes")
})

test_that("layouts that do not match Biodes' plan are flagged", {
  d <- copy_fixture("rcbd_greenhouse")
  spec <- jsonlite::read_json(file.path(d, "study_spec.json"))
  spec$factors <- list(list(name = "genotype", n_levels = 6), list(name = "water", n_levels = 3))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE)
  ok <- list(candidates = c("rcbd", "row_column"), recipe = 'book <- make_design("rcbd", paste0("G", 1:6), reps = 4)')
  expect_match(.suite_layout_issues(ok, d), "factor 'water' has 3 levels")
  bad <- list(candidates = c("met_fa", "alpha"), recipe = 'make_design("factorial", factor_levels = list(A = c("a1", "a2")))')
  iss <- .suite_layout_issues(bad, d)
  expect_true(any(grepl("first-choice design is 'met_fa'", iss)))
  expect_true(any(grepl("placeholder factor levels", iss)))
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
