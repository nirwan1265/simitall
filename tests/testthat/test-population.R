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
  expect_true(file.exists(paste0(prefix, ".causal.tsv")))
  truth <- read.delim(paste0(prefix, ".causal.tsv"))
  phenotype <- read.delim(paste0(prefix, ".pheno.tsv"))
  expect_equal(nrow(truth), 2L)
  expect_true(all(c("marker_id", "effect", "phenotype") %in% names(truth)))
  expect_true("true_breeding_value" %in% names(phenotype))
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

test_that("monomorphic breeding outputs can write an empty VCF truth set", {
  out_dir <- tempfile("simitall-monomorphic-")
  dir.create(out_dir)
  panel <- file.path(out_dir, "identical.fa")
  writeLines(c(">parent1", "ACGTACGT", ">parent2", "ACGTACGT"), panel)
  prefix <- file.path(out_dir, "f1")
  vcf <- paste0(prefix, ".vcf")

  expect_true(simulate_breeding(
    haplotype_fa = panel,
    out_prefix = prefix,
    parents = "parent1,parent2",
    n_founders = 2,
    sequence = "F1",
    n_offspring = 4,
    vcf_out = vcf,
    seed = 8
  ))

  expect_true(file.exists(vcf))
  expect_true(any(grepl("^#CHROM", readLines(vcf))))
  expect_false(any(grepl("^chr1\\t", readLines(vcf))))
})

test_that("multi-chromosome breeding uses chromosome maps and exports truth", {
  out_dir <- tempfile("simitall-multichrom-")
  dir.create(out_dir)
  panel <- file.path(out_dir, "founders.fa")
  writeLines(c(
    ">P1|chrA", paste(rep("A", 200), collapse = ""),
    ">P1|chrB", paste(rep("C", 160), collapse = ""),
    ">P2|chrA", paste(rep("G", 200), collapse = ""),
    ">P2|chrB", paste(rep("T", 160), collapse = "")
  ), panel)
  map <- file.path(out_dir, "map.tsv")
  write.table(
    data.frame(
      chromosome = rep(c("chrA", "chrB"), each = 2),
      pos_bp = c(1, 200, 1, 160),
      cM = c(0, 300, 0, 250)
    ),
    map, sep = "\t", quote = FALSE, row.names = FALSE
  )
  prefix <- file.path(out_dir, "f2")
  vcf <- paste0(prefix, ".vcf")

  expect_true(simulate_breeding(
    haplotype_fa = panel,
    out_prefix = prefix,
    parents = "P1,P2",
    n_founders = 2,
    scheme = "F2",
    n_offspring = 24,
    recomb_map_in = map,
    sv_rate = 0.01,
    sv_maxlen = 10,
    vcf_out = vcf,
    seed = 17
  ))

  fasta_headers <- sub("^>", "", grep(
    "^>", readLines(paste0(prefix, ".fa")), value = TRUE
  ))
  expect_length(fasta_headers, 24 * 2 * 2)
  expect_true(all(c("sample1_hap1|chrA", "sample1_hap1|chrB") %in%
                    fasta_headers))

  vcf_lines <- readLines(vcf)
  expect_true(any(grepl("^##contig=<ID=chrA,length=200>", vcf_lines)))
  expect_true(any(grepl("^##contig=<ID=chrB,length=160>", vcf_lines)))
  expect_true(any(grepl("^chrA\\t", vcf_lines)))
  expect_true(any(grepl("^chrB\\t", vcf_lines)))
  expect_true(any(grepl("<(DEL|INS|DUP)>", vcf_lines)))

  ancestry <- read.delim(paste0(prefix, ".ancestry.tsv"))
  breakpoints <- read.delim(paste0(prefix, ".breakpoints.tsv"))
  realized_map <- read.delim(paste0(prefix, ".recombination_map.tsv"))
  sv_truth <- read.delim(paste0(prefix, ".sv_truth.tsv"))
  expect_setequal(unique(ancestry$chromosome), c("chrA", "chrB"))
  expect_setequal(unique(realized_map$chromosome), c("chrA", "chrB"))
  expect_true(all(c(
    "sample", "haplotype", "chromosome", "position_bp",
    "left_founder", "right_founder"
  ) %in% names(breakpoints)))
  expect_gt(nrow(breakpoints), 0)
  expect_gt(nrow(sv_truth), 0)
  expect_true(all(sv_truth$chromosome %in% c("chrA", "chrB")))

  chromosome_origins <- ancestry[ancestry$start == 1, ]
  origins <- reshape(
    chromosome_origins[c("sample", "haplotype", "chromosome", "founder")],
    idvar = c("sample", "haplotype"), timevar = "chromosome",
    direction = "wide"
  )
  expect_true(any(origins$founder.chrA != origins$founder.chrB))
})

test_that("the random panel generator can emit multiple chromosomes", {
  panel <- tempfile(fileext = ".fa")
  expect_true(generate_random_haplotype_panel(
    out_fa = panel,
    n_haplotypes = 3,
    n_chromosomes = 2,
    chromosome_lengths = c(80, 60),
    snp_rate = 0.02,
    indel_rate = 0,
    seed = 18
  ))
  headers <- sub("^>", "", grep("^>", readLines(panel), value = TRUE))
  expect_setequal(headers, as.vector(outer(
    paste0("hap", 1:3), paste0("chr", 1:2), paste, sep = "|"
  )))
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
