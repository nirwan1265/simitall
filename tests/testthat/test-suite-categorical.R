# Nominal outcomes, segregation ratios, paired yes/no, traits for correlation, raters for agreement.

cat_folder <- function(po, analysis = list(), design = "crd", reps = 150, extra = list(),
                       factors = list(list(name = "treatment", n_levels = 2, levels = list("control", "treated")))) {
  testthat::skip_if_not_installed("studyforge")
  spec <- c(list(system = "field", organism = "Zea mays", factors = factors, design = list(reps = reps, chosen = design, randomization_seed = 42),
                 primary_outcome = po, power = list(residual_sd = 10, minimum_detectable_difference = 10), analysis = analysis), extra)
  biodes_folder(spec, design)
}
edit_spec <- function(dir, f) {
  s <- jsonlite::read_json(file.path(dir, "study_spec.json")); s <- f(s)
  jsonlite::write_json(s, file.path(dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
}
run_cat <- function(dir, ...) {
  suppressMessages(simitall_from_biodes(dir, seed = 1, ...))
  list(data = utils::read.delim(file.path(dir, "sample_metadata_template.tsv"), stringsAsFactors = FALSE),
       layout = utils::read.csv(file.path(dir, "field_layout.csv")),
       truth = jsonlite::read_json(file.path(dir, "sim_truth.json")))
}

test_that("nominal outcome: multinomial logit, category probabilities per group in the truth", {
  dir <- cat_folder(list(name = "colour", type = "nominal", categories = list("white", "pink", "purple")), reps = 300)
  r <- run_cat(dir)
  om <- r$truth$outcome_model
  expect_equal(r$truth$outcome_type, "nominal")
  expect_equal(om$family, "multinomial_logit")
  expect_setequal(unique(r$data$primary_outcome), c("white", "pink", "purple"))
  expect_equal(sum(unlist(om$group_probability$treated)), 1, tolerance = 1e-3)
  expect_equal(unlist(om$group_probability$control), c(white = 1, pink = 1, purple = 1) / 3, tolerance = 1e-3)
  obs <- mean(r$data$primary_outcome[r$layout$treatment == "treated"] == "purple")
  expect_equal(obs, om$group_probability$treated$purple, tolerance = 0.08)
})

test_that("segregation: classes from analysis.expected_ratio, and a distorted ratio as the planted effect", {
  dir <- cat_folder(list(name = "colour", type = "nominal", categories = list("purple", "white")), reps = 300)
  edit_spec(dir, function(s) { s$analysis$model <- "y ~ 1"; s$analysis$expected_ratio <- list(purple = 3, white = 1); s })
  r <- run_cat(dir)
  expect_equal(unlist(r$truth$outcome_model$true_ratio), c(purple = 0.75, white = 0.25))
  expect_false(r$truth$outcome_model$distorted)
  expect_equal(mean(r$data$primary_outcome == "purple"), 0.75, tolerance = 0.05)
  r <- run_cat(dir, params = list(suite_outcome = list(segregation_distortion = -0.5)))
  tr <- unlist(r$truth$outcome_model$true_ratio)
  expect_true(r$truth$outcome_model$distorted)
  expect_equal(tr[["purple"]], 0.75 * exp(-0.5) / (0.75 * exp(-0.5) + 0.25), tolerance = 1e-3)
})

test_that("paired yes/no: logit with the unit's effect; the conditional odds ratio is recovered", {
  dir <- cat_folder(list(name = "wilted", type = "binary"), design = "rcbd", reps = 400,
                    factors = list(list(name = "time", n_levels = 2, levels = list("before", "after"))))
  edit_spec(dir, function(s) { s$analysis$terms$paired <- TRUE; s$analysis$terms$repeated_factors <- list("time"); s })
  r <- run_cat(dir)
  om <- r$truth$outcome_model
  expect_equal(om$family, "logit")
  expect_equal(unlist(om$conditional_odds_ratios$treatment), c(before = 1, after = 3))
  expect_gt(om$random_effect_sd_logit$block, 0)
  d <- r$data; d$t <- factor(r$layout$time, c("before", "after")); d$b <- r$layout$block
  # conditional OR from the discordant pairs (McNemar): before 0 -> after 1 over before 1 -> after 0
  w <- stats::reshape(data.frame(b = d$b, t = as.character(d$t), y = d$primary_outcome), idvar = "b", timevar = "t", direction = "wide")
  expect_equal(sum(w$y.before == 0 & w$y.after == 1) / sum(w$y.before == 1 & w$y.after == 0), 3, tolerance = 0.35)
  # one measurement per unit keeps the liability threshold unless asked
  dir2 <- cat_folder(list(name = "wilted", type = "binary"))
  expect_equal(run_cat(dir2)$truth$outcome_model$family, "liability_threshold")
  expect_equal(run_cat(dir2, params = list(suite_outcome = list(binary_model = "logit")))$truth$outcome_model$family, "logit")
})

test_that("secondary traits: columns written, within-design residual correlation planted and in the truth", {
  dir <- cat_folder(list(name = "yield", type = "continuous"), reps = 300,
                    extra = list(secondary_outcomes = list(list(name = "plant height", column = "height", unit = "cm", baseline_mean = 200, residual_sd = 15),
                                                           list(name = "ear length", column = "ear_length"))))
  r <- run_cat(dir, params = list(suite_outcome = list(trait_residual_cor = 0.6)))
  expect_true(all(c("height", "ear_length") %in% names(r$data)))
  R <- r$truth$secondary_traits$within_design_correlation
  expect_equal(R$height$primary_outcome, 0.6)
  res <- sapply(c("primary_outcome", "height", "ear_length"), function(v) stats::resid(stats::lm(r$data[[v]] ~ r$layout$treatment)))
  expect_equal(stats::cor(res)[1, 2], 0.6, tolerance = 0.1)
  expect_equal(stats::sd(res[, "height"]), 15, tolerance = 2)
})

test_that("raters: continuous ratings with the planted ICC2; ordinal ratings with the expected kappa", {
  dir <- cat_folder(list(name = "rust", type = "continuous"), reps = 200,
                    analysis = list(agreement = list(raters = list("rater_a", "rater_b"), scale = "continuous")))
  r <- run_cat(dir)
  ag <- r$truth$agreement
  expect_equal(ag$icc2, ag$true_value_variance / (ag$true_value_variance + ag$rater_bias_sd^2 + ag$rater_noise_sd^2))
  expect_equal(stats::cor(r$data$rater_a, r$data$rater_b), ag$icc3, tolerance = 0.06)
  dir <- cat_folder(list(name = "rust", type = "ordinal"), reps = 400,
                    analysis = list(agreement = list(raters = list("rater_a", "rater_b", "rater_c"), scale = "ordinal")))
  r <- run_cat(dir)
  ag <- r$truth$agreement
  po <- mean(r$data$rater_a == r$data$rater_b)
  pe <- sum(prop.table(table(factor(r$data$rater_a, unlist(ag$levels)))) * prop.table(table(factor(r$data$rater_b, unlist(ag$levels)))))
  expect_equal((po - pe) / (1 - pe), ag$expected_kappa, tolerance = 0.07)
  expect_gt(ag$expected_weighted_kappa_quadratic, ag$expected_kappa)
})

test_that("Biodes 0.6.7: a yes/no trait with expected_ratio draws labels; raters in long format share a true score", {
  dir <- cat_folder(list(name = "height_class", type = "binary", categories = list("tall", "dwarf")), reps = 300)
  edit_spec(dir, function(s) { s$analysis$model <- "y ~ 1"; s$analysis$expected_ratio <- list(tall = 3, dwarf = 1); s })
  r <- run_cat(dir)
  expect_equal(r$truth$outcome_type, "nominal")
  expect_equal(mean(r$data$primary_outcome == "tall"), 0.75, tolerance = 0.06)
  # long format: one row per unit and rater
  dir <- cat_folder(list(name = "rust", type = "ordinal"), reps = 100)
  lay <- utils::read.csv(file.path(dir, "field_layout.csv"))
  meta <- utils::read.delim(file.path(dir, "sample_metadata_template.tsv"))
  long_lay <- lay[rep(seq_len(nrow(lay)), each = 3), ]; long_lay$rater <- rep(c("r1", "r2", "r3"), nrow(lay))
  long_meta <- meta[rep(seq_len(nrow(meta)), each = 3), ]; long_meta$rater <- long_lay$rater
  long_meta$sample_id <- paste0(long_meta$sample_id, "_", long_meta$rater)
  utils::write.csv(long_lay, file.path(dir, "field_layout.csv"), row.names = FALSE)
  utils::write.table(long_meta, file.path(dir, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  edit_spec(dir, function(s) { s$analysis$agreement <- list(rater_column = "rater", score_column = "primary_outcome", unit_column = "plot", scale = "ordinal"); s })
  r <- run_cat(dir)
  ag <- r$truth$agreement
  expect_equal(ag$format, "long")
  expect_setequal(unlist(ag$raters), c("r1", "r2", "r3"))
  w <- stats::reshape(r$data[, c("plot", "rater", "primary_outcome")], idvar = "plot", timevar = "rater", direction = "wide")
  po <- mean(w$primary_outcome.r1 == w$primary_outcome.r2)
  expect_equal(po, sum(diag(t(matrix(unlist(ag$misclassification), length(ag$levels), byrow = TRUE)) %*% diag(unlist(ag$true_class_shares)) %*%
                         matrix(unlist(ag$misclassification), length(ag$levels), byrow = TRUE))), tolerance = 0.12)
})

test_that("a run uses Mersenne-Twister whatever the session uses, and gives the caller's generator back", {
  dir <- cat_folder(list(name = "colour", type = "nominal", categories = list("tall", "dwarf")), reps = 200)
  edit_spec(dir, function(s) { s$analysis$model <- "y ~ 1"; s$analysis$expected_ratio <- list(tall = 3, dwarf = 1); s })
  a <- run_cat(dir)$data$primary_outcome
  old <- RNGkind(); on.exit(RNGkind(old[1], old[2], old[3]), add = TRUE)
  suppressWarnings(RNGkind("Super-Duper"))          # what agricolae (Biodes layout_from_spec) leaves behind
  b <- run_cat(dir)$data$primary_outcome
  expect_identical(a, b)
  expect_equal(RNGkind()[1], "Super-Duper")
})

test_that("secondary traits keep their own type (a count trait is counts)", {
  dir <- cat_folder(list(name = "plant_height", type = "continuous"), reps = 100,
                    extra = list(secondary_outcomes = list(list(name = "tiller_number", type = "count", column = "tiller_number"))))
  r <- run_cat(dir)
  x <- r$data$tiller_number
  expect_true(all(x == round(x)) && all(x >= 0))
  expect_equal(r$truth$secondary_traits$traits$tiller_number$type, "count")
  expect_gt(stats::cor(r$data$primary_outcome, x), 0.15)        # latent 0.5, thinned by the Poisson draw
})
