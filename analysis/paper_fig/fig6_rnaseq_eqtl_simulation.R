#!/usr/bin/env Rscript

# Figure 6: GWAS-linked RNA-seq and eQTL simulation
#
# This script creates the cohort, RNA-seq experiment, eQTL benchmark, combined
# publication figure, and panel-level source data using the public simitall API.
# It does not use hand-written example values in the figure.
#
# Usage:
#   Rscript analysis/paper_fig/fig6_rnaseq_eqtl_simulation.R
#   Rscript analysis/paper_fig/fig6_rnaseq_eqtl_simulation.R \
#     --out_dir analysis/results/rnaseq_eqtl --seed 42

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  index <- match(flag, args, nomatch = 0L)
  if (!index || index == length(args)) default else args[index + 1L]
}

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) {
  sub("^--file=", "", script_argument[1L])
} else {
  file.path("analysis", "paper_fig", "fig6_rnaseq_eqtl_simulation.R")
}
repo_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  mustWork = TRUE
)
setwd(repo_root)

if (requireNamespace("devtools", quietly = TRUE) &&
    file.exists(file.path(repo_root, "DESCRIPTION"))) {
  devtools::load_all(repo_root, quiet = TRUE)
} else {
  library(simitall)
}

out_dir <- get_arg("--out_dir", file.path("analysis", "results", "rnaseq_eqtl"))
seed <- as.integer(get_arg("--seed", "42"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Figure 6: GWAS-linked RNA-seq and eQTL simulation ===\n")
cat("Output directory:", normalizePath(out_dir, mustWork = TRUE), "\n")
cat("Seed:", seed, "\n\n")

demo_genome <- system.file(
  "extdata", "examples", "demo_genome.fa",
  package = "simitall"
)
if (!nzchar(demo_genome)) {
  stop("The bundled demo genome could not be located")
}

gwas_prefix <- file.path(out_dir, "cohort")
rnaseq_prefix <- file.path(out_dir, "rnaseq")
eqtl_prefix <- file.path(out_dir, "eqtl")
figure_prefix <- file.path(out_dir, "figure6_rnaseq_eqtl")

cat("1. Simulating structured GWAS cohort...\n")
simulate_gwas_cohort(
  genome_fa = demo_genome,
  out_prefix = gwas_prefix,
  n_samples = 120,
  snp_rate = 0.02,
  indel_rate = 0,
  n_pops = 2,
  fst = 0.05,
  ld_block_size = 1000,
  ld_haplotypes = 8,
  phenotype = "quantitative",
  n_causal = 10,
  effect_sd = 0.4,
  seed = seed
)

cat("2. Simulating RNA-seq for the same individuals...\n")
rna_paths <- simulate_rnaseq_from_gwas(
  genotype_file = paste0(gwas_prefix, ".geno.tsv"),
  out_prefix = rnaseq_prefix,
  sample_metadata = paste0(gwas_prefix, ".pheno.tsv"),
  n_genes = 120,
  n_cis_eqtl = 30,
  n_trans_eqtl = 10,
  cis_window_bp = 500,
  condition_levels = c("control", "treatment"),
  batch_levels = c("batch1", "batch2"),
  condition_effect_fraction = 0.15,
  gxe_fraction = 0.25,
  ase_fraction = 0.4,
  n_latent = 3,
  library_size_mean = 1e6,
  seed = seed + 1L
)

cat("3. Testing and benchmarking eQTL associations...\n")
eqtl_results <- benchmark_eqtl(
  genotype_file = paste0(gwas_prefix, ".geno.tsv"),
  expression_file = rna_paths$expression,
  out_prefix = eqtl_prefix,
  sample_metadata = rna_paths$sample_metadata,
  truth_file = rna_paths$eqtl_truth,
  pair_mode = "truth_neighborhood",
  cis_window_bp = 500,
  covariates = c("condition", "batch", "pop"),
  test_gxe = TRUE,
  fdr_threshold = 0.05
)

cat("4. Building publication figure and source-data tables...\n")
figure_results <- plot_rnaseq_eqtl_results(
  rnaseq_prefix = rnaseq_prefix,
  eqtl_prefix = eqtl_prefix,
  out_prefix = figure_prefix,
  top_variable_genes = 100,
  fdr_threshold = 0.05,
  width = 15,
  height = 10,
  dpi = 300
)

cat("\nCompleted successfully.\n")
cat("Figure PDF:", figure_results$paths$pdf, "\n")
cat("Figure PNG:", figure_results$paths$png, "\n")
cat("Result summary:", figure_results$paths$summary, "\n")
cat("Panel source data:", figure_results$paths$source_data, "\n")
cat("eQTL discoveries:", eqtl_results$metrics$discoveries, "\n")
