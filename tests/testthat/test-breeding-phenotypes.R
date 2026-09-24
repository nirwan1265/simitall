test_that("simstudy generates family, ordinal, and count traits", {
  skip_if_not_installed("simstudy", minimum_version = "0.8.1")
  skip_if_not_installed("data.table")

  records <- data.frame(
    sample = paste0("line", seq_len(120)),
    family = rep(paste0("fam", seq_len(12)), each = 10),
    genetic_score = seq(-1.5, 1.5, length.out = 120),
    stringsAsFactors = FALSE
  )
  records <- simulate_family_effects(
    records,
    family_id = "family",
    variance = 0.3,
    seed = 101
  )

  family_effect_counts <- tapply(
    records$family_effect,
    records$family,
    function(value) length(unique(value))
  )
  expect_true(all(family_effect_counts == 1L))
  expect_equal(length(unique(records$family_effect)), 12L)

  records$liability <- records$genetic_score + records$family_effect
  ordinal <- simulate_ordinal_trait(
    records,
    predictor = "liability",
    probabilities = c(0.15, 0.35, 0.35, 0.15),
    labels = c("low", "moderate", "high", "extreme"),
    seed = 102
  )
  expect_true(is.ordered(ordinal$ordinal_trait))
  expect_equal(
    levels(ordinal$ordinal_trait),
    c("low", "moderate", "high", "extreme")
  )
  expect_identical(ordinal$sample, records$sample)

  poisson <- simulate_count_trait(
    records,
    formula = "0.5 + 0.4 * genetic_score + family_effect",
    distribution = "poisson",
    seed = 103
  )
  expect_true(all(poisson$count_trait >= 0))
  expect_true(all(poisson$count_trait == floor(poisson$count_trait)))

  negative_binomial <- simulate_count_trait(
    records,
    formula = "0.5 + 0.4 * genetic_score + family_effect",
    distribution = "negative_binomial",
    dispersion = 0.7,
    trait_name = "disease_count",
    seed = 104
  )
  expect_true(all(negative_binomial$disease_count >= 0))
  expect_true(all(
    negative_binomial$disease_count ==
      floor(negative_binomial$disease_count)
  ))
})

test_that("pedtricks generates direct and maternal truth components", {
  skip_if_not_installed("pedtricks", minimum_version = "0.4.2")

  pedigree <- data.frame(
    id = c(
      "dam1", "sire1", "dam2", "sire2",
      "child1", "child2", "child3", "child4"
    ),
    sire = c(
      NA, NA, NA, NA,
      "sire1", "sire1", "sire2", "sire2"
    ),
    dam = c(
      NA, NA, NA, NA,
      "dam1", "dam1", "dam2", "dam2"
    ),
    stringsAsFactors = FALSE
  )
  out_prefix <- file.path(tempfile("simitall-parental-"), "traits")
  result <- simulate_parental_effects(
    pedigree,
    genetic_covariance = diag(c(1, 0.35)),
    environmental_covariance = diag(c(0.8, 0.25)),
    parental_genetic = c("d", "m"),
    parental_environmental = c("d", "m"),
    out_prefix = out_prefix,
    seed = 105
  )

  expect_equal(nrow(result$phenotypes), nrow(pedigree))
  expect_true(all(c(
    "a_tr1", "bv_m_tr1", "e_tr1", "e_m_tr1",
    "P_bv_m_tr1", "P_e_m_tr1", "Phen_tr1"
  ) %in% names(result$allEffects)))
  expect_true(all(file.exists(unlist(result$paths))))
})

test_that("breeding phenotype wrappers reject malformed identifiers", {
  skip_if_not_installed("simstudy", minimum_version = "0.8.1")
  skip_if_not_installed("data.table")

  expect_error(
    simulate_family_effects(
      data.frame(sample = 1:2, family = c("a", NA))
    ),
    "cannot be missing"
  )
  expect_error(
    simulate_ordinal_trait(
      data.frame(sample = c("x", "x"), genetic_score = c(0, 1))
    ),
    "unique"
  )
})
