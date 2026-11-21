# biopaint_biodes(): figures from any design folder.

test_that("figures are drawn from whatever the folder holds, for different designs", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("ggplot2")
  for (case in list(c("factorial", "y ~ treatment * treatment2 + block"), c("alpha", "y ~ treatment + rep + (1 | rep:block)"),
                    c("repeated_measures", "y ~ treatment * time + (1 | unit)"), c("row_column", "y ~ treatment + rep + (1 | rep:row) + (1 | rep:col)"))) {
    spec <- list(handoff_version = "0.4.0", factors = list(list(name = "genotype", n_levels = 6), list(name = "water", n_levels = 2), list(name = "time", n_levels = 3)),
                 design = list(reps = 3, randomization_seed = 1), power = list(minimum_detectable_difference = 4, residual_sd = 2), primary_outcome = list(name = "y"))
    book <- simitall_layout_from_spec(spec, case[1])
    d <- simitall:::.suite_write_package(book, spec, case[2], file.path(tempfile(), case[1]))
    suppressMessages(simitall_from_biodes(d))
    suppressMessages(bioana_fit_biodes(d))
    f <- suppressMessages(biopaint_biodes(d))
    expect_true(all(c("layout", "outcome", "recovery") %in% names(f)), info = case[1])
    expect_true(all(file.exists(f)), info = case[1])
  }
})
