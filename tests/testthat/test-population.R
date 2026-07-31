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

test_that("simplePHENOTYPES generates a quantitative trait", {
  skip_if_not_installed("simplePHENOTYPES")
  out_dir <- tempfile("simitall-phenotype-")
  dir.create(out_dir)
  genome <- system.file(
    "extdata", "examples", "demo_genome.fa",
    package = "simitall"
  )
  cohort_prefix <- file.path(out_dir, "cohort")
  trait_prefix <- file.path(out_dir, "trait")

  simulate_gwas_cohort(
    genome_fa = genome,
    out_prefix = cohort_prefix,
    n_samples = 12,
    snp_rate = 0.01,
    indel_rate = 0,
    n_causal = 2,
    seed = 71
  )
  expect_true(simulate_phenotypes(
    geno_file = paste0(cohort_prefix, ".vcf"),
    out_prefix = trait_prefix,
    h2 = 0.4,
    n_add_qtn = 2,
    n_reps = 1,
    seed = 72
  ))

  phenotype_path <- paste0(trait_prefix, ".pheno.tsv")
  expect_true(file.exists(phenotype_path))
  phenotype <- read.delim(phenotype_path, check.names = FALSE)
  expect_equal(nrow(phenotype), 12)
  expect_true(any(vapply(phenotype[-1L], is.numeric, logical(1L))))
  expect_true(file.exists(paste0(trait_prefix, ".simulation_summary.txt")))
})
