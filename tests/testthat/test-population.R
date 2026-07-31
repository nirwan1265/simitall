test_that("a small GWAS cohort is generated", {
  out_dir <- tempfile("simitall-gwas-")
  dir.create(out_dir)
  genome <- system.file(
    "extdata", "examples", "demo_genome.fa",
    package = "simitall"
  )
  prefix <- file.path(out_dir, "cohort")

  expect_true(simulate_gwas_cohort(
    genome_fa = genome,
    out_prefix = prefix,
    n_samples = 8,
    snp_rate = 0.002,
    indel_rate = 0,
    n_causal = 2,
    seed = 5
  ))

  expect_true(file.exists(paste0(prefix, ".vcf")))
  expect_true(file.exists(paste0(prefix, ".geno.tsv")))
  expect_true(file.exists(paste0(prefix, ".pheno.tsv")))
})

test_that("a small F2 population is generated", {
  out_dir <- tempfile("simitall-f2-")
  dir.create(out_dir)
  panel <- system.file(
    "extdata", "panels", "demo_panel.fa",
    package = "simitall"
  )
  prefix <- file.path(out_dir, "f2")
  vcf <- paste0(prefix, ".vcf")

  expect_true(simulate_breeding(
    haplotype_fa = panel,
    out_prefix = prefix,
    scheme = "F2",
    n_offspring = 6,
    n_founders = 2,
    vcf_out = vcf,
    seed = 7
  ))

  expect_true(file.exists(paste0(prefix, ".fa")))
  expect_true(file.exists(paste0(prefix, ".meta.tsv")))
  expect_true(file.exists(vcf))
})
