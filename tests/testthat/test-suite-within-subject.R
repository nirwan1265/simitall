# Repeated measures: correlation among one unit's occasions (suite_outcome within_subject), recorded in sim_truth.json.

repeated_folder <- function(reps = 40) {
  testthat::skip_if_not_installed("studyforge")
  tf <- list(name = "time", n_levels = 4, levels = list("1 wk", "2 wk", "4 wk", "8 wk"), values = list(1, 2, 4, 8), unit = "wk", scale = "linear")
  spec <- list(system = "animal", organism = "Bos taurus", factors = list(list(name = "diet", n_levels = 2), tf),
               repeated_factors = list("time"),
               design = list(reps = reps, chosen = "randomized_complete_block_repeated", randomization_seed = 42),
               primary_outcome = list(name = "milk", type = "continuous"), power = list(residual_sd = 10, minimum_detectable_difference = 5))
  biodes_folder(spec, "randomized_complete_block_repeated")
}

goal_like_folder <- function() {
  testthat::skip_if_not_installed("studyforge")
  spec <- list(system = "field", organism = "Zea mays", factors = list(list(name = "nitrogen", n_levels = 3)),
               design = list(reps = 4, chosen = "rcbd", randomization_seed = 42),
               primary_outcome = list(name = "yield", type = "continuous"), power = list(residual_sd = 10, minimum_detectable_difference = 5))
  biodes_folder(spec, "rcbd")
}

within_cor <- function(dir) {
  m <- utils::read.delim(file.path(dir, "sample_metadata_template.tsv"))
  lay <- utils::read.csv(file.path(dir, "field_layout.csv"))
  m <- cbind(m["primary_outcome"], lay[, c("plot", "diet", "time")])
  m$r <- m$primary_outcome - stats::ave(m$primary_outcome, m$diet, m$time)
  w <- stats::reshape(m[, c("plot", "time", "r")], idvar = "plot", timevar = "time", direction = "wide")
  list(cor = stats::cor(w[, -1]), var = apply(w[, -1], 2, stats::var))
}

test_that("within_subject none keeps the random intercept only and says so", {
  dir <- repeated_folder()
  suppressMessages(simitall_from_biodes(dir, seed = 2))
  t <- jsonlite::read_json(file.path(dir, "sim_truth.json"))
  expect_equal(t$within_subject$structure, "random_intercept")
  expect_equal(t$within_subject$over, "time")
})

test_that("within_subject ar1: adjacent occasions correlate more than distant ones; rho recorded", {
  dir <- repeated_folder()
  suppressMessages(simitall_from_biodes(dir, seed = 2, params = list(suite_outcome = list(within_subject = "ar1", within_rho = 0.8))))
  t <- jsonlite::read_json(file.path(dir, "sim_truth.json"))
  expect_equal(t$within_subject$structure, "ar1")
  expect_equal(t$within_subject$rho, 0.8)
  expect_equal(t$within_subject$over, "time")
  expect_equal(t$within_subject$units, 80)
  cm <- within_cor(dir)$cor
  expect_gt(mean(cm[cbind(1:3, 2:4)]) - cm[1, 4], 0.1)
})

test_that("within_subject random_slope: variance grows over time; slope SD recorded; bad values stop", {
  dir <- repeated_folder()
  suppressMessages(simitall_from_biodes(dir, seed = 2, params = list(suite_outcome = list(within_subject = "random_slope", slope_sd_ratio = 2))))
  t <- jsonlite::read_json(file.path(dir, "sim_truth.json"))
  expect_equal(t$within_subject$structure, "random_slope")
  expect_equal(t$within_subject$slope_sd, 20)
  v <- within_cor(dir)$var
  expect_gt(v[[4]], 1.5 * v[[1]])
  expect_error(suppressMessages(simitall_from_biodes(dir, seed = 2, params = list(suite_outcome = list(within_subject = "toeplitz")))), "none, ar1 or random_slope")
})

test_that("a per-row id column in the layout does not hide the units; no repeated factor warns", {
  dir <- repeated_folder()
  lay <- utils::read.csv(file.path(dir, "field_layout.csv"))
  lay$sample_id <- sprintf("S%03d", seq_len(nrow(lay)))
  utils::write.csv(lay, file.path(dir, "field_layout.csv"), row.names = FALSE)
  suppressMessages(simitall_from_biodes(dir, seed = 2, params = list(suite_outcome = list(within_subject = "ar1"))))
  t <- jsonlite::read_json(file.path(dir, "sim_truth.json"))
  expect_equal(t$within_subject$structure, "ar1")
  expect_equal(t$within_subject$units, 80)
  expect_false("sample_id" %in% unlist(t$within_subject$unit_columns))
  dir2 <- goal_like_folder()
  expect_warning(suppressMessages(simitall_from_biodes(dir2, seed = 1, params = list(suite_outcome = list(within_subject = "ar1")))),
                 "not applied: the design has no repeated factor")
  t2 <- jsonlite::read_json(file.path(dir2, "sim_truth.json"))
  expect_equal(t2$within_subject$structure, "independent")
  expect_equal(t2$within_subject$requested, "ar1")
})
