test_that("dependency profiles contain the expected components", {
  minimal <- check_simitall_dependencies("minimal", verbose = FALSE)
  population <- check_simitall_dependencies("population", verbose = FALSE)
  omics <- check_simitall_dependencies("omics", verbose = FALSE)
  sequencing <- check_simitall_dependencies("sequencing", verbose = FALSE)
  full <- check_simitall_dependencies("full", verbose = FALSE)

  expect_setequal(minimal$component, c("jsonlite", "Matrix", "reticulate"))
  expect_true(all(c("simplePHENOTYPES", "SimuPOP") %in% population$component))
  expect_true(all(c("Rsubread", "ChIPsim", "splatter") %in% omics$component))
  expect_true(all(c("ART", "PBSIM/PBSIM3", "QUAST") %in%
                    sequencing$component))
  expect_gt(nrow(full), nrow(population))
  expect_true(all(c(
    "component", "kind", "required_for", "available", "detected",
    "install_hint"
  ) %in% names(full)))
})

test_that("dependency installation supports a safe dry run", {
  command <- install_simitall_dependencies(
    "minimal", envname = "simitall-test", ask = FALSE, dry_run = TRUE
  )
  expect_true(file.exists(command$command))
  expect_identical(command$profile, "minimal")
  expect_identical(command$envname, "simitall-test")
  expect_true(any(grepl("install_simitall_", command$arguments)))
  expect_error(
    install_simitall_dependencies(
      "minimal", envname = "bad environment", dry_run = TRUE
    ),
    "envname"
  )
})
