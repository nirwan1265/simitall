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

test_that("compact gene-coordinate catalogues retain complete gene intervals", {
  manifest_path <- system.file("extdata", "gene_coordinates", "manifest.tsv", package = "simitall")
  manifest <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(nrow(manifest), 6L)

  for (i in seq_len(nrow(manifest))) {
    path <- system.file("extdata", "gene_coordinates", manifest$coordinate_file[[i]], package = "simitall")
    expect_true(file.exists(path))
    genes <- read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
    expect_equal(nrow(genes), manifest$gene_records[[i]])
    expect_true(all(nzchar(genes$seqname)))
    expect_true(all(nzchar(genes$gene_id)))
    expect_true(all(genes$start_bp <= genes$end_bp))
  }
})

test_that("GAF annotations join to compact gene coordinates", {
  gaf <- system.file("extdata", "go", "arabidopsis_tair10_go_demo.gaf.gz", package = "simitall")
  expect_true(file.exists(gaf))
  annotations <- read_go_gaf(gaf)
  expect_gt(nrow(annotations), 0L)
  expect_true(all(grepl("^GO:", annotations$go_id)))

  coordinates <- read_gene_coordinates("Arabidopsis_thaliana")
  joined <- join_gene_go(coordinates, annotations)
  expect_gt(nrow(joined), 0L)
  expect_true(all(joined$match_method == "stable_id"))
  summary <- summarize_go_terms(joined, by = "go_id_aspect")
  expect_gt(nrow(summary), 0L)
  expect_true(all(summary$gene_count <= summary$annotation_count))
})
