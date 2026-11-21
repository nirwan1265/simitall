# Suggestions for what the study details leave open (inst/suite/suggestions.tsv + Biodes registries).

test_that("missing study details get suggestions matched to the question, never tied to one prompt", {
  ans <- list(question = "single cell sequencing in mouse uterus", spec = list(system = "mouse", organism = "Mus musculus"))
  sug <- simitall:::.suite_suggestions(ans)
  expect_setequal(unique(sug$field), c("factors", "experimental_unit", "primary_outcome", "power"))
  expect_true(any(grepl("pseudobulk", sug$suggestion)))
  filled <- simitall:::.suite_suggested_answers(sug, ans$spec)
  expect_equal(filled$factors, "genotype:2")
  # power is filled as a pair (difference and SD), never one without the other
  expect_true(all(c("power.minimum_detectable_difference", "power.residual_sd") %in% names(filled)))

  plant <- list(question = "drought on sorghum in the greenhouse", spec = list(system = "greenhouse",
                factors = list(list(name = "genotype", n_levels = 8)), power = list(minimum_detectable_difference = 2, residual_sd = 3)))
  s2 <- simitall:::.suite_suggestions(plant)
  expect_false("factors" %in% s2$field)
  expect_false("power" %in% s2$field)

  other <- simitall:::.suite_suggestions(list(question = "something new", spec = list()))
  expect_true("treatment:2" %in% sub('factors = "(.*)"', "\\1", other$answer[other$field == "factors"]))
})
