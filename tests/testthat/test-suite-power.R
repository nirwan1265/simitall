# simitall_power_sweep(): generic over designs; reads everything from the folder.

test_that("power sweep runs from a design folder and power grows with replication", {
  skip_if_not_installed("bioana")
  skip_if_not_installed("lme4")
  spec <- list(handoff_version = "0.4.0", factors = list(list(name = "genotype", n_levels = 4), list(name = "water", n_levels = 2)),
               design = list(reps = 3, chosen = "factorial", randomization_seed = 42),
               power = list(minimum_detectable_difference = 6, residual_sd = 3), primary_outcome = list(name = "rwc"))
  d <- biodes_folder(spec, "gh")
  ps <- simitall_power_sweep(d, reps = c(2, 6), n_sim = 40)
  y <- as.data.frame(unclass(ps))
  expect_true(file.exists(file.path(d, "power_sweep.tsv")))
  expect_setequal(unique(y$term), c("treatment", "treatment2", "block", "treatment:treatment2"))
  expect_true(is.na(y$simulated_effect[y$term == "block"][1]))
  # interaction simulated by default (simulation_defaults.json interaction_effect_ratio)
  expect_true(y$simulated_effect[y$term == "treatment:treatment2"][1])
  w <- y[y$term == "treatment:treatment2", ]
  expect_gt(w$power[w$reps == 6], w$power[w$reps == 2])
  g <- y[y$term == "treatment", ]
  expect_gte(g$power[g$reps == 6], g$power[g$reps == 2])
  # no interaction: its row is the false-positive rate, about alpha
  y0 <- as.data.frame(unclass(simitall_power_sweep(d, reps = 6, n_sim = 60, interaction_effect = 0)))
  expect_false(y0$simulated_effect[y0$term == "treatment:treatment2"][1])
  expect_lt(y0$power[y0$term == "treatment:treatment2"], 0.2)
  expect_output(print(ps), "Power at alpha")
})

test_that("power sweep works for repeated-measures and crossover designs too", {
  skip_if_not_installed("bioana")
  skip_if_not_installed("lme4")
  for (case in list(c("repeated_measures"), c("crossover_ab_ba"))) {
    spec <- list(handoff_version = "0.4.0", factors = list(list(name = "drug", n_levels = 2), list(name = "time", n_levels = 3)),
                 design = list(reps = 3, chosen = case[1], randomization_seed = 1),
                 power = list(minimum_detectable_difference = 5, residual_sd = 2), primary_outcome = list(name = "y"))
    d <- biodes_folder(spec, case[1])
    y <- as.data.frame(unclass(simitall_power_sweep(d, reps = c(2, 4), n_sim = 10)))
    expect_true("treatment" %in% y$term, info = case[1])
  }
})

test_that("a question with no replicate count gets design.reps from the power sweep", {
  skip_if_not_installed("bioana")
  skip_if_not_installed("studyforge")
  spec <- list(handoff_version = "0.4.0", organism = "Zea mays", system = "field",
               factors = list(list(name = "water", n_levels = 2)),
               design = list(chosen = "rcbd", randomization_seed = 42),
               power = list(minimum_detectable_difference = 6, residual_sd = 3), primary_outcome = list(name = "yield"))
  expect_true(simitall:::.suite_reps_open(spec))
  expect_true(simitall:::.suite_reps_open(utils::modifyList(spec, list(design = list(reps = 4), power = list(solve_for = "reps")))))
  ch <- NULL
  utils::capture.output(ch <- suppressMessages(simitall:::.suite_choose_reps(spec, tempdir(), 42, reps_range = c(2, 4, 8), n_sim = 30)))
  expect_true(ch$reps %in% c(2, 4, 8))
  y <- as.data.frame(unclass(ch$sweep))
  ok <- y$power[y$term == "treatment"] >= 0.8
  expect_equal(ch$reps, if (any(ok)) c(2, 4, 8)[which(ok)[1]] else 8)
})

test_that("count assays get power as the share of 2-fold genes found", {
  skip_if_not_installed("bioana")
  skip_if_not_installed("studyforge")
  spec <- list(handoff_version = "0.4.0", organism = "Secale cereale", system = "growth_chamber",
               factors = list(list(name = "genotype", n_levels = 2)), assays = list("rnaseq"),
               design = list(reps = 3, chosen = "crd", randomization_seed = 42))
  d <- biodes_folder(spec, "rna")
  y <- as.data.frame(unclass(suppressMessages(simitall_power_sweep(d, reps = c(2, 8), n_sim = 8, show = FALSE, allow_fallback = TRUE))))
  de <- y[startsWith(y$term, "DE genes"), ]
  expect_equal(nrow(de), 2)
  expect_gt(de$power[de$reps == 8], de$power[de$reps == 2])
})

test_that("Biodes' [min, max] replicate range becomes a grid of counts to try", {
  expect_equal(simitall:::.suite_reps_grid(list(3, 12)), c(3L, 4L, 5L, 7L, 9L, 12L))
  expect_equal(simitall:::.suite_reps_grid(c(2, 3, 4, 6)), c(2L, 3L, 4L, 6L))
  expect_equal(simitall:::.suite_reps_grid(c(3, 4)), c(3L, 4L))
})
