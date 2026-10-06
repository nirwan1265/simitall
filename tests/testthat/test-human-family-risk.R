# Recessive family-risk model: the relationship effect must come from
# inheritance, and the simulation must agree with the exact expectation.

test_that("the default pedigree model is unchanged", {
  markers <- simulate_human_marker_frequencies(12L, seed = 2L)
  run <- simulate_human_pedigree_groups(markers, n_per_group = 10L, seed = 3L)
  expect_null(run$model)
  expect_false("expected_probability" %in% names(run$individuals))
  expect_equal(ncol(run$genotypes), 12L)
})

test_that("the recessive intercept is calibrated to the requested baseline", {
  q <- rep(0.05, 10L)
  a <- .simitall_recessive_intercept(0.01, q, 4)
  expect_equal(.simitall_recessive_expected(0, a, q, 4), 0.01, tolerance = 1e-8)
  # Risk rises with inbreeding even though the model has no relatedness term.
  expect_gt(.simitall_recessive_expected(1 / 64, a, q, 4), 0.01)
  expect_gt(.simitall_recessive_expected(1 / 16, a, q, 4), .simitall_recessive_expected(1 / 64, a, q, 4))
})

test_that("simulated family risk agrees with theory and rises with relatedness", {
  markers <- simulate_human_marker_frequencies(16L, seed = 4L)
  run <- simulate_human_pedigree_groups(markers, n_per_group = 400L, seed = 5L,
                                        baseline_probability = 0.02, risk_model = "recessive")
  summary <- summarize_family_risk(run)
  expect_equal(summary$group, c("unrelated", "distantly_related", "first_cousin_descendant"))
  expect_equal(summary$expected_probability[1L], 0.02, tolerance = 1e-8)
  expect_true(all(diff(summary$expected_probability) > 0))
  expect_true(all(summary$agrees_with_theory))
  expect_equal(ncol(run$genotypes), 16L + 10L)
  expect_equal(sum(run$marker_truth$role == "recessive_risk_locus"), 10L)
  path <- plot_family_risk(summary, file.path(tempdir(), "family_risk.png"))
  expect_true(file.exists(path))
})
