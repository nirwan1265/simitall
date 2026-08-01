test_that("annotation-aware ChIP-seq writes truth, reads, and QC", {
  out_dir <- tempfile("simitall-chipseq-")
  dir.create(out_dir)
  genome <- system.file(
    "extdata", "examples", "demo_genome.fa", package = "simitall"
  )
  seqname <- sub("^>([^ ]+).*", "\\1", readLines(genome, n = 1L))
  gff <- file.path(out_dir, "targets.gff3")
  writeLines(c(
    "##gff-version 3",
    paste(seqname, "test", "promoter", 500, 650, ".", "+", ".",
          "ID=promoter1", sep = "\t"),
    paste(seqname, "test", "enhancer", 1500, 1700, ".", "+", ".",
          "ID=enhancer1", sep = "\t"),
    paste(seqname, "test", "promoter", 3000, 3150, ".", "-", ".",
          "ID=promoter2", sep = "\t"),
    paste(seqname, "test", "enhancer", 5000, 5200, ".", "+", ".",
          "ID=enhancer2", sep = "\t")
  ), gff)
  prefix <- file.path(out_dir, "tf")

  paths <- simulate_chipseq(
    genome_fa = genome,
    annotation_gff3 = gff,
    out_prefix = prefix,
    assay_type = "TF",
    n_peaks = 6,
    peak_width = 120,
    peak_width_sd = 10,
    conditions = c("control", "treated"),
    biological_replicates = 2,
    technical_replicates = 1,
    differential_binding_fraction = 0.5,
    n_reads = 300,
    input_reads = 300,
    read_length = 40,
    fragment_mean = 100,
    fragment_sd = 15,
    fragment_min = 60,
    fragment_max = 150,
    backend = "native",
    write_fastq = TRUE,
    seed = 111
  )

  expect_identical(attr(paths, "backend"), "native")
  expect_true(all(file.exists(unlist(paths))))
  peaks <- read.delim(paths$peaks, check.names = FALSE)
  metadata <- read.delim(paths$sample_metadata, check.names = FALSE)
  qc <- read.delim(paths$qc, check.names = FALSE)
  manifest <- read.delim(paths$read_manifest, check.names = FALSE)
  expect_equal(nrow(peaks), 6)
  expect_equal(nrow(metadata), 4)
  expect_equal(nrow(qc), 8)
  expect_equal(nrow(manifest), 8)
  expect_true(all(file.exists(manifest$fastq)))
  expect_true(all(qc$frip >= 0 & qc$frip <= 1))
  expect_gt(mean(qc$frip[qc$assay == "chip"]),
            mean(qc$frip[qc$assay == "input"]))

  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  figure <- plot_chipseq_results(
    chip_prefix = prefix,
    out_prefix = file.path(out_dir, "figure"),
    profile_window = 400,
    profile_bin = 25,
    max_profile_reads = 2400,
    width = 10,
    height = 7,
    dpi = 72
  )
  expect_true(file.exists(figure$paths$pdf))
  expect_true(file.exists(figure$paths$png))
  expect_true(file.exists(figure$paths$summary))
  expect_length(list.files(figure$paths$source_data, pattern = "\\.csv$"), 6)
})

test_that("the optional ChIPsim backend samples strand-specific positions", {
  skip_if_not_installed("ChIPsim")
  genome <- c(chr1 = paste(rep("ACGT", 500), collapse = ""))
  peaks <- data.frame(
    peak_id = "peak_00001", seqname = "chr1", start = 900L, end = 1100L,
    center = 1000L, width = 201L, peak_shape = "narrow",
    stringsAsFactors = FALSE
  )
  positions <- .simitall_chip_chipsim_positions(
    nreads = 100, genome = genome, peaks = peaks, peak_strength = 10,
    signal_fraction = 0.4, read_length = 36, fragment_mean = 120,
    fragment_min = 80, fragment_max = 180
  )
  expect_equal(nrow(positions), 100)
  expect_true(all(positions$strand %in% c("+", "-")))
  expect_true(all(positions$read_position >= 1))
})
