#!/usr/bin/env Rscript

# Figure 7: donor-aware single-cell RNA-seq and cell-type eQTL simulation
#
# This runner starts from the bundled genome and produces a structured GWAS
# cohort, single-cell counts, donor-cell-type pseudobulk expression, cell-type
# eQTL benchmarks, a six-panel figure, and panel-level source data.
#
# Usage:
#   Rscript analysis/paper_fig/fig7_scrnaseq_celltype_eqtl.R
#   Rscript analysis/paper_fig/fig7_scrnaseq_celltype_eqtl.R \
#     --out_dir analysis/results/scrnaseq_celltype_eqtl --seed 52 \
#     --backend native

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  index <- match(flag, args, nomatch = 0L)
  if (!index || index == length(args)) default else args[index + 1L]
}

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) {
  sub("^--file=", "", script_argument[1L])
} else {
  file.path("analysis", "paper_fig", "fig7_scrnaseq_celltype_eqtl.R")
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

out_dir <- get_arg(
  "--out_dir",
  file.path("analysis", "results", "scrnaseq_celltype_eqtl")
)
seed <- as.integer(get_arg("--seed", "52"))
backend <- get_arg("--backend", "auto")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Figure 7: donor-aware scRNA-seq and cell-type eQTLs ===\n")
cat("Output directory:", normalizePath(out_dir, mustWork = TRUE), "\n")
cat("Seed:", seed, "\n")
cat("Backend:", backend, "\n\n")

demo_genome <- system.file(
  "extdata", "examples", "demo_genome.fa",
  package = "simitall"
)
if (!nzchar(demo_genome)) {
  stop("The bundled demo genome could not be located")
}

gwas_prefix <- file.path(out_dir, "cohort")
scrna_prefix <- file.path(out_dir, "scrna")
eqtl_prefix <- file.path(out_dir, "celltype_eqtl")
figure_prefix <- file.path(out_dir, "figure7_scrnaseq_celltype_eqtl")

cat("1. Simulating structured GWAS donors...\n")
simulate_gwas_cohort(
  genome_fa = demo_genome,
  out_prefix = gwas_prefix,
  n_samples = 40,
  snp_rate = 0.03,
  indel_rate = 0,
  n_pops = 2,
  fst = 0.05,
  ld_block_size = 1000,
  ld_haplotypes = 8,
  phenotype = "quantitative",
  n_causal = 8,
  effect_sd = 0.4,
  seed = seed
)

cat("2. Simulating cells from the same GWAS donors...\n")
sc <- simulate_scrnaseq_from_gwas(
  genotype_file = paste0(gwas_prefix, ".geno.tsv"),
  out_prefix = scrna_prefix,
  sample_metadata = paste0(gwas_prefix, ".pheno.tsv"),
  n_genes = 120,
  cells_per_donor = 48,
  cell_type_proportions = c(
    T_cell = 0.35,
    B_cell = 0.25,
    Monocyte = 0.25,
    Stromal = 0.15
  ),
  trajectory_cell_types = c("T_cell", "Monocyte"),
  n_cis_eqtl = 24,
  n_trans_eqtl = 8,
  cis_window_bp = 700,
  cell_type_eqtl_fraction = 0.5,
  condition_levels = c("control", "treatment"),
  batch_levels = c("batch1", "batch2"),
  condition_effect_fraction = 0.15,
  gxe_fraction = 0.25,
  markers_per_cell_type = 6,
  marker_effect_mean = 1.8,
  pseudotime_gene_fraction = 0.12,
  cis_effect_sd = 0.9,
  trans_effect_sd = 0.6,
  donor_effect_sd = 0.15,
  residual_sd = 0.08,
  library_size_mean = 3000,
  dropout_rate = 0.05,
  doublet_rate = 0.03,
  ambient_fraction = 0.02,
  backend = backend,
  write_sce = TRUE,
  min_cells_pseudobulk = 5,
  seed = seed + 1L
)

cat("3. Benchmarking cell-type eQTLs with donor pseudobulk...\n")
benchmark <- benchmark_celltype_eqtls(
  genotype_file = paste0(gwas_prefix, ".geno.tsv"),
  pseudobulk_expression = sc$pseudobulk_expression,
  pseudobulk_metadata = sc$pseudobulk_metadata,
  truth_file = sc$eqtl_truth,
  out_prefix = eqtl_prefix,
  cis_window_bp = 700,
  covariates = c("condition", "batch", "pop"),
  test_gxe = TRUE,
  fdr_threshold = 0.05,
  min_donors = 20
)

cat("4. Building the publication figure and source-data tables...\n")
figure <- plot_scrnaseq_results(
  scrna_prefix = scrna_prefix,
  eqtl_prefix = eqtl_prefix,
  out_prefix = figure_prefix,
  top_variable_genes = 100,
  markers_per_cell_type = 4,
  max_plot_cells = 2000,
  width = 15,
  height = 10,
  dpi = 300,
  seed = seed + 2L
)

cat("\nCompleted successfully.\n")
cat("Backend used:", attr(sc, "backend"), "\n")
cat("Single-cell object:", sc$object_rds, "\n")
cat("Pseudobulk expression:", sc$pseudobulk_expression, "\n")
cat("eQTL metrics:", benchmark$metrics_path, "\n")
cat("Figure PDF:", figure$paths$pdf, "\n")
cat("Figure PNG:", figure$paths$png, "\n")
cat("Panel source data:", figure$paths$source_data, "\n")

