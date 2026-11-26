# Covariates and regression (Bioana docs/contracts/regression.md): simulated per unit, slopes planted, truth recorded.

cov_folder <- function(model, covs, po = list(name = "y", type = "continuous"), reps = 100, ge = NULL) {
  testthat::skip_if_not_installed("studyforge")
  spec <- list(system = "field", organism = "Zea mays", factors = list(list(name = "treatment", n_levels = 2, levels = list("control", "treated"))),
               design = list(reps = reps, chosen = "crd", randomization_seed = 42), primary_outcome = po,
               power = list(residual_sd = 10, minimum_detectable_difference = 5))
  d <- biodes_folder(spec, "crd")
  s <- jsonlite::read_json(file.path(d, "study_spec.json"))
  s$analysis$model <- model; s$analysis$covariates <- covs
  s$analysis$goals <- list("regression"); s$analysis$primary_goal <- "regression"
  s$analysis$goal_terms <- list(regression = lapply(covs, `[[`, "term"))
  if (!is.null(ge)) s$analysis$goal_effect <- list(regression = ge)
  jsonlite::write_json(s, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  d
}
run_cov <- function(d, ...) {
  suppressMessages(simitall_from_biodes(d, seed = 1, ...))
  list(data = utils::read.delim(file.path(d, "sample_metadata_template.tsv")), layout = utils::read.csv(file.path(d, "field_layout.csv")),
       truth = jsonlite::read_json(file.path(d, "sim_truth.json")))
}

test_that("covariates are simulated per unit (mean, sd, correlation, log transform) and the slopes recovered by lm", {
  d <- cov_folder("y ~ covariate + covariate2",
                  list(list(term = "covariate", name = "glucose", column = "glucose", mean = 100, sd = 15),
                       list(term = "covariate2", name = "BMI", column = "bmi", transform = "log", mean = 25, sd = 4)), reps = 200)
  r <- run_cov(d, params = list(suite_covariates = list(covariate_cor = 0.5)))
  x <- r$data
  expect_equal(mean(x$glucose), 100, tolerance = 0.03)
  expect_equal(stats::sd(x$glucose), 15, tolerance = 0.1)
  expect_true(all(x$bmi > 0))
  expect_equal(stats::cor(x$glucose, log(x$bmi)), 0.5, tolerance = 0.12)
  gt <- r$truth$goal_truth
  expect_equal(gt$goal, "regression")
  expect_setequal(names(gt$true_slopes), c("covariate", "covariate2"))
  fit <- stats::lm(x$primary_outcome ~ x$glucose + log(x$bmi))
  expect_equal(unname(stats::coef(fit)[2]), gt$true_slopes$covariate, tolerance = 0.3)
  expect_equal(unname(stats::coef(fit)[3]), gt$true_slopes$covariate2, tolerance = 0.3)
  expect_equal(r$truth$regression_truth$covariate_cor, 0.5)
  expect_gt(r$truth$regression_truth$r_squared, 0)
})

test_that("a stated slope is planted per unit; yes/no outcomes get a log odds ratio per unit", {
  d <- cov_folder("y ~ treatment + covariate", list(list(term = "covariate", name = "baseline", column = "baseline")), ge = list(size = 0.4, scale = "per_unit"))
  expect_equal(run_cov(d)$truth$goal_truth$true_slopes$covariate, 0.4)
  d <- cov_folder("y ~ covariate", list(list(term = "covariate", name = "dose", column = "dose")), po = list(name = "dead", type = "binary"),
                  reps = 300, ge = list(size = 1.1, scale = "odds_ratio_per_unit"))
  r <- run_cov(d)
  expect_equal(r$truth$outcome_model$family, "logit")
  expect_equal(r$truth$goal_truth$true_log_odds_ratios$covariate, log(1.1), tolerance = 1e-6)
  fit <- stats::glm(r$data$primary_outcome ~ r$data$dose, family = stats::binomial())
  expect_equal(unname(stats::coef(fit)[2]), log(1.1), tolerance = 0.5)
})

test_that("covariate interactions report the coefficients of the uncentred model; planted problems are recorded", {
  d <- cov_folder("y ~ covariate * covariate2", list(list(term = "covariate", column = "x1"), list(term = "covariate2", column = "x2")), reps = 300)
  r <- run_cov(d)
  ts <- r$truth$goal_truth$true_slopes
  fit <- stats::lm(primary_outcome ~ x1 * x2, data = r$data)
  expect_equal(unname(stats::coef(fit)["x1:x2"]), ts[["covariate:covariate2"]], tolerance = 0.5)
  expect_equal(unname(stats::coef(fit)["x1"]), ts$covariate, tolerance = 0.5)
  d <- cov_folder("y ~ covariate", list(list(term = "covariate", column = "x1")))
  r <- run_cov(d, params = list(suite_covariates = list(heteroscedasticity = 0.8, outliers = 3)))
  rt <- r$truth$regression_truth
  expect_length(rt$planted_outliers, 3)
  expect_equal(rt$heteroscedasticity$strength, 0.8)
  res <- stats::resid(stats::lm(primary_outcome ~ x1 + r$layout$treatment, data = r$data))
  hi <- r$data$x1 > stats::median(r$data$x1)
  expect_gt(stats::sd(res[hi]), stats::sd(res[!hi]))
})
