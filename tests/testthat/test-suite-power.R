# simitall_power_sweep(): generic over designs; reads everything from the folder.

test_that("power sweep runs from a design folder and power grows with replication", {
  skip_if_not_installed("lme4")
  spec <- list(handoff_version = "0.4.0", factors = list(list(name = "genotype", n_levels = 4), list(name = "water", n_levels = 2)),
               design = list(reps = 3, chosen = "factorial", randomization_seed = 42),
               power = list(minimum_detectable_difference = 6, residual_sd = 3), primary_outcome = list(name = "rwc"))
  book <- simitall_layout_from_spec(spec, "factorial")
  d <- simitall:::.suite_write_package(book, spec, "y ~ treatment * treatment2 + block", file.path(tempfile(), "gh"))
  ps <- simitall_power_sweep(d, reps = c(2, 6), n_sim = 40)
  y <- as.data.frame(unclass(ps))
  expect_true(file.exists(file.path(d, "power_sweep.tsv")))
  expect_setequal(unique(y$term), c("treatment", "treatment2", "block", "treatment:treatment2"))
  expect_true(is.na(y$simulated_effect[y$term == "block"][1]))
  expect_false(y$simulated_effect[y$term == "treatment:treatment2"][1])
  w <- y[y$term == "treatment2", ]
  expect_gt(w$power[w$reps == 6], w$power[w$reps == 2])
  expect_lt(max(y$power[y$term == "treatment:treatment2"]), 0.2)   # no simulated interaction: about alpha
  expect_output(print(ps), "Power at alpha")
})

test_that("power sweep works for repeated-measures and crossover designs too", {
  skip_if_not_installed("lme4")
  for (case in list(c("repeated_measures", "y ~ treatment * time + (1 | unit)"),
                    c("crossover_ab_ba", "y ~ treatment + period + (1 | subject)"))) {
    spec <- list(handoff_version = "0.4.0", factors = list(list(name = "drug", n_levels = 2), list(name = "time", n_levels = 3)),
                 design = list(reps = 3, chosen = case[1], randomization_seed = 1),
                 power = list(minimum_detectable_difference = 5, residual_sd = 2), primary_outcome = list(name = "y"))
    book <- simitall_layout_from_spec(spec, case[1])
    d <- simitall:::.suite_write_package(book, spec, case[2], file.path(tempfile(), case[1]))
    y <- as.data.frame(unclass(simitall_power_sweep(d, reps = c(2, 4), n_sim = 10)))
    expect_true("treatment" %in% y$term, info = case[1])
  }
})
