# Survival outcomes: PH / AFT / non-PH models, censoring, competing events, frailty, truth for Bioana's scoring.

surv_folder <- function(reps = 150, extra_po = list(), analysis = list(), design = "crd", power = list()) {
  testthat::skip_if_not_installed("studyforge")
  spec <- list(system = "greenhouse", organism = "Zea mays",
               factors = list(list(name = "seed_treatment", n_levels = 2, levels = list("control", "primed"))),
               design = list(reps = reps, chosen = design, randomization_seed = 42),
               primary_outcome = c(list(name = "days_to_germination", type = "survival"), extra_po),
               power = c(list(residual_sd = 10, minimum_detectable_difference = 10), power), analysis = analysis)
  biodes_folder(spec, design)
}
read_surv <- function(dir) {
  m <- utils::read.delim(file.path(dir, "sample_metadata_template.tsv"))
  m$g <- factor(utils::read.csv(file.path(dir, "field_layout.csv"))$seed_treatment, c("control", "primed"))
  list(data = m, truth = jsonlite::read_json(file.path(dir, "sim_truth.json")), manifest = jsonlite::read_json(file.path(dir, "data_manifest.json")))
}

test_that("weibull_ph (default): Cox HR and Weibull time ratio match the truth; medians and landmarks recorded", {
  skip_if_not_installed("survival")
  dir <- surv_folder(reps = 400, analysis = list(survival = list(landmarks = list(10, 30), rmst_tau = 40)))
  suppressMessages(simitall_from_biodes(dir, seed = 3))
  r <- read_surv(dir); om <- r$truth$outcome_model
  expect_equal(om$family, "weibull_ph")
  expect_true(om$proportional_hazards)
  expect_equal(unlist(om$ratios$treatment$hazard_ratio), c(1, 0.6), ignore_attr = TRUE)
  expect_equal(unlist(om$ratios$treatment$time_ratio)[[2]], 0.6^(-1 / 1.5), tolerance = 1e-3)
  hr <- exp(stats::coef(survival::coxph(survival::Surv(primary_outcome, event) ~ g, r$data)))
  expect_equal(unname(hr), 0.6, tolerance = 0.25)
  expect_equal(om$groups$control$median, 30, tolerance = 0.01)
  expect_named(om$groups$primed$survival_at, c("10", "30"))
  expect_equal(om$rmst_tau, 40)
  gt <- r$truth$goal_truth
  expect_equal(gt$goal, "survival")     # Biodes' goal for a survival outcome; rmst_tau stated -> RMST difference
  expect_match(gt$effect_scale, "RMST difference at tau = 40: group primed minus reference group control")
  expect_equal(gt$true_effect, om$groups$primed$rmst - om$groups$control$rmst, tolerance = 1e-6)
  expect_gt(gt$true_effect, 0)
  expect_equal(names(om$groups), c("control", "primed"))
})

test_that("a stated target HR sizes the effect; loss to follow-up is counted apart from administrative censoring", {
  dir <- surv_folder(power = list(target_hazard_ratio = 0.5, loss_fraction = 0.15), reps = 200)
  suppressMessages(simitall_from_biodes(dir, seed = 1))
  om <- read_surv(dir)$truth$outcome_model
  expect_equal(om$hazard_ratio, 0.5)
  expect_match(om$hazard_ratio_source, "power.target_hazard_ratio|goal_effect.survival")
  expect_equal(om$censoring$loss$loss_fraction, 0.15)
  expect_gt(om$censoring$loss$n, 20)
  expect_gt(om$censoring$administrative$n, 20)
  expect_equal(om$n_events + om$censoring$loss$n + om$censoring$administrative$n, 400)
})

test_that("delayed_effect and crossing are non-PH; lognormal_aft gives time ratios only", {
  dir <- surv_folder()
  for (mdl in c("delayed_effect", "crossing")) {
    suppressMessages(simitall_from_biodes(dir, seed = 1, params = list(suite_outcome = list(survival_model = mdl, change_time = 10))))
    om <- read_surv(dir)$truth$outcome_model
    expect_false(om$proportional_hazards)
    expect_equal(om$change_time, 10)
  }
  expect_equal(unlist(om$ratios$treatment$hazard_ratio_after_change)[[2]], 1 / 0.6, tolerance = 1e-3)
  suppressMessages(simitall_from_biodes(dir, seed = 1, params = list(suite_outcome = list(survival_model = "lognormal_aft"))))
  om <- read_surv(dir)$truth$outcome_model
  expect_null(om$ratios$treatment$hazard_ratio)
  expect_gt(unlist(om$ratios$treatment$time_ratio)[[2]], 1)
  expect_error(suppressMessages(simitall_from_biodes(dir, params = list(suite_outcome = list(survival_model = "gompertz")))), "survival_model must be one of")
})

test_that("competing events from event_codes: codes 0/1/2, CIF truth, manifest carries the codes", {
  dir <- surv_folder(extra_po = list(event_codes = list("0" = "censored", "1" = "germinated", "2" = "seed rot")),
                     analysis = list(survival = list(landmarks = list(20))))
  suppressMessages(simitall_from_biodes(dir, seed = 1))
  r <- read_surv(dir); om <- r$truth$outcome_model
  expect_setequal(unique(r$data$event), c(0, 1, 2))
  expect_equal(unlist(om$competing$codes), 2)
  cif <- om$groups$control$cif_event_at[["20"]]; s <- om$groups$control$survival_at[["20"]]
  expect_lt(cif, 1 - s)          # part of the drop in event-free survival is the competing event
  ent <- r$manifest$data[[1]]
  expect_equal(ent$event_column, "event")
  expect_equal(ent$event_codes[["2"]], "seed rot")
})

test_that("interval scoring rounds times up to the grid; a survival goal reports the HR", {
  dir <- surv_folder(analysis = list(survival = list(life_table = 2), goals = list("survival")))
  s <- jsonlite::read_json(file.path(dir, "study_spec.json"))
  s$analysis$primary_goal <- "survival"; s$analysis$goal_effect <- list(survival = list(scale = "hazard_ratio", size = 0.7))
  jsonlite::write_json(s, file.path(dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  suppressMessages(simitall_from_biodes(dir, seed = 1))
  r <- read_surv(dir)
  expect_true(all(r$data$primary_outcome %% 2 == 0))
  gt <- r$truth$goal_truth
  expect_equal(gt$goal, "survival")
  expect_equal(gt$true_effect, 0.7, tolerance = 1e-6)
  expect_match(gt$effect_scale, "hazard ratio")
})

test_that("Biodes 0.6.6 fields: logrank_weight picks the non-PH model, censoring text gives follow-up, strata recorded", {
  dir <- surv_folder(extra_po = list(censoring = "administrative at 21 days (end of follow-up)"),
                     analysis = list(survival = list(logrank_weight = "fh_0_1")), design = "rcbd", reps = 6)
  s <- jsonlite::read_json(file.path(dir, "study_spec.json"))
  s$analysis$survival$strata <- list("block")
  jsonlite::write_json(s, file.path(dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  suppressMessages(simitall_from_biodes(dir, seed = 1))
  r <- read_surv(dir); om <- r$truth$outcome_model
  expect_equal(om$family, "delayed_effect")
  expect_equal(om$follow_up, 21)
  expect_lte(max(r$data$primary_outcome), 21)
  # no median stated: set from the follow-up so 20% of the reference group is still without the event at day 21
  expect_match(om$median_time_source, "from follow-up 21")
  expect_equal(om$median_time, 21 * (log(2) / -log(0.2))^(1 / 1.5), tolerance = 1e-6)
  expect_length(om$strata$block$log_hazard_shift, 6)
  s$analysis$survival$logrank_weight <- "fh_1_0"
  jsonlite::write_json(s, file.path(dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  suppressMessages(simitall_from_biodes(dir, seed = 1))
  om <- read_surv(dir)$truth$outcome_model
  expect_equal(om$family, "early_effect")
  expect_equal(unlist(om$ratios$treatment$hazard_ratio_after_change), c(1, 1), ignore_attr = TRUE)
  # params win over the spec's weight
  suppressMessages(simitall_from_biodes(dir, seed = 1, params = list(suite_outcome = list(survival_model = "weibull_ph"))))
  expect_equal(read_surv(dir)$truth$outcome_model$family, "weibull_ph")
})

test_that("strata outside the model (Biodes: stated blocks < 10, no + block) get their own baseline hazard", {
  dir <- surv_folder(design = "rcbd", reps = 6)
  s <- jsonlite::read_json(file.path(dir, "study_spec.json"))
  s$analysis$model <- "y ~ treatment"; s$analysis$survival <- list(strata = list("block"))
  jsonlite::write_json(s, file.path(dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  suppressMessages(simitall_from_biodes(dir, seed = 1, params = list(suite_outcome = list(strata_log_hazard_sd = 0.5))))
  om <- read_surv(dir)$truth$outcome_model
  expect_length(om$strata$block$log_hazard_shift, 6)
  expect_equal(om$strata$block$sd, 0.5)
  expect_match(om$strata$block$note, "strata\\(\\)")
})
