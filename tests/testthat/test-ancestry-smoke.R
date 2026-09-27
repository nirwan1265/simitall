test_that("GenomeAdmixR ancestry example writes known-truth outputs when available", {
  local_libs <- normalizePath(
    c(".Rlib", "../.Rlib", "../../.Rlib"),
    mustWork = FALSE
  )
  local_libs <- local_libs[dir.exists(local_libs)]
  if (length(local_libs)) .libPaths(c(local_libs, .libPaths()))
  if (!requireNamespace("GenomeAdmixR", quietly = TRUE)) {
    skip("GenomeAdmixR is not installed in the active R library")
  }
  script <- normalizePath(
    file.path("..", "..", "inst", "extdata", "examples", "ancestry_tract_smoke.R"),
    mustWork = TRUE
  )
  out_dir <- tempfile("simitall-ancestry-smoke-")
  status <- system2(file.path(R.home("bin"), "Rscript"), c(script, out_dir))
  expect_identical(status, 0L)
  expect_true(file.exists(file.path(out_dir, "final_ancestry_frequency.csv")))
  expect_true(file.exists(file.path(out_dir, "ancestry_frequency_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "final_ancestry_frequency.png")))
  summary <- utils::read.csv(file.path(out_dir, "ancestry_frequency_summary.csv"))
  expect_equal(nrow(summary), 3L)
  expect_true(all(summary$mean_frequency >= 0 & summary$mean_frequency <= 1))
})
