test_that("bundled examples are available", {
  genome <- system.file(
    "extdata", "examples", "demo_genome.fa",
    package = "simitall"
  )
  panel <- system.file(
    "extdata", "panels", "demo_panel.fa",
    package = "simitall"
  )
  multi_panel <- system.file(
    "extdata", "panels", "demo_multichrom_panel.fa",
    package = "simitall"
  )

  expect_true(nzchar(genome))
  expect_true(nzchar(panel))
  expect_true(file.exists(genome))
  expect_true(file.exists(panel))
  expect_true(nzchar(multi_panel))
  expect_true(file.exists(multi_panel))
})
