# biopaint_biodes(): figures from any design folder.

test_that("figures are drawn from whatever the folder holds, for different designs", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("bioana")
  g <- function(n) list(name = "genotype", n_levels = n)
  cases <- list(factorial = list(g(6), list(name = "water", n_levels = 2)), alpha = list(g(12)),
                repeated_measures = list(list(name = "drug", n_levels = 2), list(name = "time", n_levels = 3)), row_column = list(g(6)))
  for (case in names(cases)) {
    spec <- list(handoff_version = "0.4.0", factors = cases[[case]],
                 design = list(reps = 3, chosen = case, randomization_seed = 1), power = list(minimum_detectable_difference = 4, residual_sd = 2), primary_outcome = list(name = "y"))
    d <- biodes_folder(spec, case)
    suppressMessages(simitall_from_biodes(d))
    suppressMessages(bioana::bioana_analyze(d, quiet = TRUE, allow_fallback = TRUE))
    f <- suppressMessages(biopaint_biodes(d))
    expect_true(all(c("layout", "outcome", "recovery") %in% names(f)), info = case)
    expect_true(all(file.exists(f)), info = case)
  }
})
