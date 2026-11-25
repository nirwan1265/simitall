# Analysis goals (Biodes schema 0.6): truth planted for analysis.primary_goal, written to sim_truth.json goal_truth.

goal_folder <- function(design, factors, goals, extra = list(), reps = 3, assays = NULL, mdd = 10) {
  testthat::skip_if_not_installed("studyforge")
  spec <- list(system = "field", organism = "Zea mays", factors = factors,
               design = list(reps = reps, chosen = design, randomization_seed = 42),
               primary_outcome = list(name = "yield", type = "continuous"),
               power = list(residual_sd = 10, minimum_detectable_difference = mdd),
               analysis = c(list(goals = as.list(goals)), extra))
  if (!is.null(assays)) spec$assays <- as.list(assays)
  biodes_folder(spec, design)
}

sim_goal <- function(dir, ...) {
  suppressMessages(simitall_from_biodes(dir, seed = 1, ...))
  list(truth = jsonlite::read_json(file.path(dir, "sim_truth.json")),
       manifest = jsonlite::read_json(file.path(dir, "data_manifest.json")))
}

test_that("goals.tsv covers every goal in the Biodes registry", {
  skip_if_not_installed("studyforge")
  reg <- studyforge::load_registry("analysis_goals")$goals
  tab <- simitall:::.suite_tsv("goals")
  expect_setequal(tab$goal, vapply(reg, `[[`, "", "id"))
  tests <- stats::setNames(vapply(reg, `[[`, "", "test"), vapply(reg, `[[`, "", "id"))
  expect_equal(tab$test, unname(tests[tab$goal]))
  off <- vapply(Filter(function(g) isFALSE(g$simulated), reg), `[[`, "", "id")
  expect_true(all(tab$planted[tab$goal %in% off] == "none"))
})

test_that("time shapes: 0 at the first time point, largest change 1, peak and half-time in the factor's unit", {
  pos <- c(a = 0, b = 6, c = 24, d = 72)
  tr <- simitall:::.suite_time_shape(pos, "transient", peak = list(value = 1, unit = "d"), factor_unit = "h")
  expect_equal(unname(tr["a"]), 0)
  expect_equal(names(which.max(tr)), "c")
  expect_equal(diff(range(tr)), 1)
  sat <- simitall:::.suite_time_shape(pos, "saturating", half_time = list(value = 6, unit = "h"), factor_unit = "h")
  expect_true(all(diff(sat) > 0)); expect_equal(unname(sat["d"]), 1)
  expect_equal(unname(simitall:::.suite_time_shape(pos)), unname(pos / 72))
  del <- simitall:::.suite_time_shape(pos, "delayed", half_time = list(value = 48, unit = "h"), factor_unit = "h")
  expect_lt(del[["c"]], 0.2)
})

test_that("interaction goal: the difference of differences of the extreme levels is the goal size", {
  dir <- goal_folder("factorial", list(list(name = "genotype", n_levels = 3), list(name = "water", n_levels = 2)),
                     "interaction", mdd = 6)
  r <- sim_goal(dir)
  gt <- r$truth$goal_truth
  expect_equal(gt$goal, "interaction")
  expect_equal(gt$status, "planted")
  expect_equal(gt$size, 6)
  expect_equal(gt$true_effect, 6)
  expect_length(gt$difference_of_differences, 2)
  expect_equal(unlist(gt$columns), c(treatment = "genotype", treatment2 = "water"))
  expect_equal(sort(all.vars(str2lang(unlist(gt$model_terms)))), c("treatment", "treatment2"))
  expect_equal(r$manifest$goal$goal, "interaction")
  expect_equal(r$manifest$data[[1]]$truth$goal$field, "goal_truth")
  expect_false(any(grepl("^interaction effect =", unlist(r$truth$assumptions))))
})

test_that("time_course goal, transient within one group: the reference group stays flat, the other peaks", {
  tf <- list(name = "time", n_levels = 4, levels = list("0 h", "6 h", "24 h", "72 h"), values = list(0, 6, 24, 72), unit = "h", scale = "linear")
  dir <- goal_folder("factorial", list(list(name = "infection", n_levels = 2), tf), "time_course",
                     extra = list(time_course = list(shape = "transient", peak = list(value = 24, unit = "h"), scope = "within_group")),
                     mdd = 20, reps = 4)
  r <- sim_goal(dir)
  gt <- r$truth$goal_truth
  expect_equal(gt$goal, "time_course")
  expect_equal(gt$shape, "transient")
  expect_equal(gt$scope, "within_group")
  tp <- unlist(gt$time_pattern)
  expect_equal(unname(tp[1]), 0)
  expect_equal(names(which.max(tp)), "24 h")
  expect_equal(gt$true_effect, 20)
  cells <- unlist(r$truth$interaction_effects[[unlist(gt$model_terms)]]$cells)
  ref <- gt$group_reference
  expect_true(all(cells[grepl(ref, names(cells), fixed = TRUE)] == 0))
  expect_equal(max(cells), 20)
  # groups equal at the first time point: no main effect for the group factor
  expect_true(all(unlist(r$truth$treatment_effects[[gt$group]]) == 0))
})

test_that("trend goal on a log-spaced dose: slope per log10 unit, and genes in counts follow it", {
  dose <- list(name = "dose", n_levels = 4, levels = list("0 uM", "1 uM", "10 uM", "100 uM"), values = list(0, 1, 10, 100), unit = "uM", scale = "log")
  dir <- goal_folder("rcbd", list(dose), "trend", assays = "rnaseq", reps = 4)
  r <- sim_goal(dir, n_genes = 300)
  gt <- r$truth$goal_truth
  expect_equal(gt$goal, "trend")
  expect_equal(gt$scale, "log")
  expect_equal(unname(unlist(gt$x)), c(-1, 0, 1, 2))
  # goal_effect.size (here power.minimum_detectable_difference) is a change per unit: per log10 uM
  expect_equal(gt$slope, 10)
  expect_equal(gt$change_over_range, 30)
  expect_true(file.exists(file.path(dir, gt$features$rnaseq$file)))
  ft <- utils::read.delim(file.path(dir, gt$features$rnaseq$file))
  expect_equal(nrow(ft), 30)
  expect_equal(ft$effect, ft$change / 3, tolerance = 1e-6)
  cnt <- r$manifest$data[[which(vapply(r$manifest$data, `[[`, "", "id") == "rnaseq")]]
  expect_equal(cnt$truth$goal$file, gt$features$rnaseq$file)
})

test_that("specs without goals (schema 0.5) keep treatment_effect; registry goals Simitall cannot simulate say so", {
  dir <- goal_folder("rcbd", list(list(name = "nitrogen", n_levels = 3)), "treatment_effect")
  s <- jsonlite::read_json(file.path(dir, "study_spec.json"))
  s$analysis[c("goals", "primary_goal", "goal_terms", "goal_effect")] <- NULL
  s$handoff_version <- "0.5.3"
  jsonlite::write_json(s, file.path(dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  gt <- sim_goal(dir)$truth$goal_truth
  expect_equal(gt$goal, "treatment_effect")
  expect_equal(unlist(gt$model_terms), "treatment")
  expect_equal(unlist(gt$level_differences$treatment), c(0, 5, 10), ignore_attr = TRUE)
  s$analysis$primary_goal <- "cell_composition"
  jsonlite::write_json(s, file.path(dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  gt <- sim_goal(dir)$truth$goal_truth
  expect_equal(gt$status, "not_simulated")
})
