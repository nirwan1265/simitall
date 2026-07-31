test_that("a small random genome and truth files are generated", {
  out_dir <- tempfile("simitall-genome-")
  dir.create(out_dir)
  out_fa <- file.path(out_dir, "genome.fa")

  expect_true(simulate_genome(
    out_fa = out_fa,
    random_genome = TRUE,
    random_length = 2000,
    n_events = 2,
    seg_len = 25,
    copies = 2,
    seed = 11
  ))

  expect_true(file.exists(out_fa))
  prefix <- tools::file_path_sans_ext(out_fa)
  expect_true(file.exists(paste0(prefix, ".repeats.tsv")))
  expect_true(file.exists(paste0(prefix, ".repeats.bed")))
  expect_true(file.exists(paste0(prefix, ".repeats.gff3")))
  expect_true(file.exists(paste0(prefix, ".summary.json")))
})

test_that("a random haplotype panel is generated", {
  out_fa <- tempfile(fileext = ".fa")

  expect_true(generate_random_haplotype_panel(
    out_fa,
    n_haplotypes = 4,
    length = 1000,
    seed = 3
  ))

  headers <- readLines(out_fa)
  expect_length(grep("^>", headers), 4)
})
