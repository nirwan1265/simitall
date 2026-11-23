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

test_that("bundled gene-annotation slices are valid lightweight GFF3 examples", {
  files <- c(
    "maize_b73_nam5_gene_slice.gff3",
    "arabidopsis_tair10_gene_slice.gff3",
    "wheat_iwgsc_gene_slice.gff3",
    "rice_irgsp10_gene_slice.gff3",
    "sorghum_ncbi_v3_gene_slice.gff3",
    "fly_bdgp6_gene_slice.gff3"
  )
  paths <- vapply(files, function(file) {
    system.file("extdata", "annotations", file, package = "simitall")
  }, character(1))
  expect_true(all(nzchar(paths)))
  expect_true(all(file.exists(paths)))
  expect_true(all(vapply(paths, function(path) {
    rows <- readLines(path, warn = FALSE)
    sum(vapply(strsplit(rows[!grepl("^#", rows)], "\t", fixed = TRUE),
               function(x) length(x) >= 9L && identical(x[[3L]], "gene"), logical(1))) == 25L
  }, logical(1))))

  manifest <- system.file("extdata", "annotations", "annotation_manifest.tsv", package = "simitall")
  expect_true(file.exists(manifest))
  expect_equal(nrow(read.delim(manifest, check.names = FALSE)), 7L)
})


