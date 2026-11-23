test_that("GWAS individuals can drive RNA-seq and eQTL simulation", {
  out_dir <- tempfile("simitall-rnaseq-")
  dir.create(out_dir)
  genome <- system.file(
    "extdata", "examples", "demo_genome.fa",
    package = "simitall"
  )
  gwas_prefix <- file.path(out_dir, "cohort")
  rna_prefix <- file.path(out_dir, "rna")

  simulate_gwas_cohort(
    genome_fa = genome,
    out_prefix = gwas_prefix,
    n_samples = 12,
    snp_rate = 0.01,
    indel_rate = 0,
    n_causal = 2,
    seed = 21
  )
  paths <- simulate_rnaseq_from_gwas(
    genotype_file = paste0(gwas_prefix, ".geno.tsv"),
    out_prefix = rna_prefix,
    n_genes = 20,
    n_cis_eqtl = 5,
    n_trans_eqtl = 3,
    cis_window_bp = 500,
    n_latent = 2,
    library_size_mean = 1e5,
    seed = 22
  )

  expect_true(all(file.exists(unlist(paths))))
  counts <- read.delim(paths$counts, check.names = FALSE)
  metadata <- read.delim(paths$sample_metadata, check.names = FALSE)
  truth <- read.delim(paths$eqtl_truth, check.names = FALSE)

  expect_equal(nrow(counts), 20)
  expect_equal(names(counts)[-1L], metadata$sample)
  expect_equal(nrow(metadata), 12)
  expect_true(all(c("condition", "batch", "condition_code") %in% names(metadata)))
  expect_true(all(c(
    "type", "gene_id", "variant_id", "effect_log2", "gxe_effect_log2"
  ) %in% names(truth)))
  expect_equal(sum(truth$type == "cis"), 5)
  expect_equal(sum(truth$type == "trans"), 3)
})


test_that("VCF genotypes and GFF3 genes interoperate", {
  out_dir <- tempfile("simitall-rnaseq-gff-")
  dir.create(out_dir)
  genome <- system.file(
    "extdata", "examples", "demo_genome.fa",
    package = "simitall"
  )
  seqname <- sub("^>", "", strsplit(readLines(genome, n = 1L), " ")[[1L]][1L])
  gff <- file.path(out_dir, "genes.gff3")
  starts <- seq(500L, 6500L, length.out = 6L)
  writeLines(c(
    "##gff-version 3",
    vapply(seq_along(starts), function(i) {
      paste(
        seqname, "test", "gene", starts[i], starts[i] + 399L,
        ".", "+", ".", paste0("ID=g", i), sep = "\t"
      )
    }, character(1L))
  ), gff)
  gwas_prefix <- file.path(out_dir, "cohort")
  simulate_gwas_cohort(
    genome_fa = genome,
    out_prefix = gwas_prefix,
    n_samples = 10,
    snp_rate = 0.01,
    indel_rate = 0,
    n_causal = 2,
    seed = 41
  )

  paths <- simulate_rnaseq_from_gwas(
    genotype_file = paste0(gwas_prefix, ".vcf"),
    out_prefix = file.path(out_dir, "rna"),
    annotation_gff3 = gff,
    n_cis_eqtl = 3,
    n_trans_eqtl = 0,
    cis_window_bp = 1000,
    n_latent = 0,
    library_size_mean = 1e5,
    seed = 42
  )

  expect_equal(nrow(read.delim(paths$gene_metadata)), 6)
  expect_equal(nrow(read.delim(paths$eqtl_truth)), 3)
})

test_that("Rsubread can emit small paired-end FASTQ files", {
  skip_if_not_installed("Rsubread")
  out_dir <- tempfile("simitall-rnaseq-reads-")
  dir.create(out_dir)
  transcript_fasta <- file.path(out_dir, "transcripts.fa")
  writeLines(c(
    ">gene1", paste(rep("ACGT", 80), collapse = ""),
    ">gene2", paste(rep("TGCA", 80), collapse = "")
  ), transcript_fasta)
  expression <- matrix(
    c(10, 20, 5, 15),
    nrow = 2,
    dimnames = list(c("gene1", "gene2"), c("sample1", "sample2"))
  )

  manifest <- simulate_rnaseq_reads(
    expression = expression,
    transcript_fasta = transcript_fasta,
    out_dir = file.path(out_dir, "fastq"),
    library_size = 20,
    read_length = 75,
    paired_end = TRUE,
    seed = 7
  )

  expect_true(all(file.exists(manifest$read1)))
  expect_true(all(file.exists(manifest$read2)))
  expect_true(file.exists(attr(manifest, "manifest_path")))
})
