#!/usr/bin/env Rscript

# Figure 15: simple NAM quantitative-trait GWAS.
#
# This runner intentionally performs only the requested workflow: a NAM-style
# population, one synthetic drought-response trait, a GWAS, truth recovery,
# and phenotype/GWAS figures. It does not generate RNA-seq or selection data.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args, nomatch = 0L)
  if (!i || i == length(args)) default else args[i + 1L]
}
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) sub("^--file=", "", script_argument[1L]) else "analysis/paper_fig/fig15_maize_nam_quantitative_gwas.R"
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
if (!dir.exists(repo_root)) repo_root <- normalizePath(".")
setwd(repo_root)
if (requireNamespace("devtools", quietly = TRUE)) devtools::load_all(repo_root, quiet = TRUE) else library(simitall)

out_dir <- get_arg("--out_dir", file.path("analysis", "results", "maize_nam_quantitative_gwas"))
seed <- as.integer(get_arg("--seed", "1501"))
n_lines <- as.integer(get_arg("--n_lines", "100"))
heritability <- as.numeric(get_arg("--heritability", "0.60"))
input_mode <- get_arg("--input_mode", "package_demo")
provided_panel <- get_arg("--haplotype_fa", "")
provided_map <- get_arg("--recomb_map", "")
if (!is.finite(heritability) || heritability <= 0 || heritability >= 1) {
  stop("--heritability must be between zero and one")
}
if (!input_mode %in% c("package_demo", "synthetic", "provided")) {
  stop("--input_mode must be package_demo, synthetic, or provided")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(out_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

if (identical(input_mode, "package_demo")) {
  panel <- system.file("extdata", "panels", "demo_maize_nam_chr10.fa", package = "simitall")
  map <- system.file("extdata", "maps", "demo_maize_nam_chr10_map.tsv", package = "simitall")
  stopifnot(nzchar(panel), nzchar(map), file.exists(panel), file.exists(map))
} else if (identical(input_mode, "synthetic")) {
  input_dir <- file.path(out_dir, "inputs")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  panel <- file.path(input_dir, "synthetic_nam_founders_chr10.fa")
  generate_random_haplotype_panel(
    out_fa = panel, n_haplotypes = 8L, length = 4000L, n_chromosomes = 1L,
    snp_rate = 0.02, indel_rate = 0, seed = seed + 101L
  )
  panel_lines <- readLines(panel)
  panel_lines[grepl("^>", panel_lines)] <- paste0(panel_lines[grepl("^>", panel_lines)], "|chr10")
  writeLines(panel_lines, panel)
  map <- file.path(input_dir, "synthetic_nam_chr10_map.tsv")
  write.table(
    data.frame(chromosome = "chr10", pos_bp = c(1L, 1000L, 1800L, 2600L, 4000L), cM = c(0, 1, 1.1, 1.2, 2.6)),
    map, sep = "\t", row.names = FALSE, quote = FALSE
  )
} else {
  if (!nzchar(provided_panel) || !nzchar(provided_map) || !file.exists(provided_panel) || !file.exists(provided_map)) {
    stop("Provided mode requires existing --haplotype_fa and --recomb_map files")
  }
  panel <- normalizePath(provided_panel)
  map <- normalizePath(provided_map)
}

prefix <- file.path(out_dir, "nam_quantitative")
vcf <- paste0(prefix, ".vcf")
cat("1. Simulating a NAM population with ", n_lines, " lines...\n", sep = "")
simulate_breeding(
  haplotype_fa = panel, out_prefix = prefix, scheme = "NAM", n_offspring = n_lines,
  recomb_map_in = map, vcf_out = vcf,
  ancestry_out = paste0(prefix, ".ancestry.tsv"),
  breakpoints_out = paste0(prefix, ".breakpoints.tsv"), seed = seed
)
geno <- simitall:::.simitall_read_genotypes(vcf)
G <- geno$genotype
variants <- geno$variants
samples <- colnames(G)
metadata <- read.delim(paste0(prefix, ".meta.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
metadata <- metadata[match(samples, metadata$sample), , drop = FALSE]
maf <- pmin(rowMeans(G, na.rm = TRUE) / 2, 1 - rowMeans(G, na.rm = TRUE) / 2)
eligible <- which(is.finite(maf) & maf >= 0.20 & maf <= 0.45)
if (length(eligible) < 2L) stop("The NAM panel did not contain enough common markers for the programmed trait")
focal <- eligible[which.min(abs(maf[eligible] - 0.35))]
background <- setdiff(eligible, focal)
background <- sample(background, min(6L, length(background)))
effects <- c(1.25, stats::rnorm(length(background), 0, 0.16))
genetic_value <- as.numeric(t(G[c(focal, background), , drop = FALSE]) %*% effects)
residual_sd <- sqrt(stats::var(genetic_value) * (1 - heritability) / heritability)
drought_response_index <- 10 + genetic_value + stats::rnorm(n_lines, 0, residual_sd)
phenotype <- data.frame(
  sample = samples, family = metadata$family,
  drought_response_index = drought_response_index,
  genetic_value = genetic_value, focal_marker_dosage = G[focal, ],
  stringsAsFactors = FALSE
)
truth <- data.frame(
  marker_id = variants$id[c(focal, background)], seqname = variants$seqname[c(focal, background)],
  pos = variants$pos[c(focal, background)], effect = effects,
  role = c("focal_drought_response_QTL", rep("polygenic_background", length(background))),
  stringsAsFactors = FALSE
)
phenotype_file <- file.path(out_dir, "drought_response_phenotypes.tsv")
truth_file <- file.path(out_dir, "drought_response_trait_truth.tsv")
write.table(phenotype, phenotype_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(truth, truth_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("2. Running family-aware quantitative GWAS...\n")
gwas_prefix <- file.path(out_dir, "drought_response_gwas")
gwas <- analyze_gwas(
  vcf, phenotype, gwas_prefix, trait = "drought_response_index",
  fixed_effects = "family", n_pcs = 2L, min_maf = 0.05
)
benchmark <- benchmark_gwas(gwas$paths$results, truth, gwas_prefix, window_bp = 150L)
plot_gwas_results(
  benchmark$paths$results, file.path(figure_dir, "drought_response_gwas_overview"),
  individual_panels_dir = figure_dir
)
plot_trait_diagnostics(
  phenotype, "drought_response_index", figure_dir,
  file_stem = "drought_response_phenotype_distribution"
)

results <- read.delim(benchmark$paths$results, check.names = FALSE)
metrics <- read.delim(benchmark$paths$metrics, check.names = FALSE)
pass <- data.frame(
  nam_lines = n_lines, programmed_heritability = heritability,
  trait_is_quantitative = is.numeric(phenotype$drought_response_index),
  focal_qtl_markers = sum(truth$role == "focal_drought_response_QTL"),
  gwas_markers = nrow(results),
  all_core_outputs_exist = all(file.exists(c(
    vcf, phenotype_file, truth_file, gwas$paths$results, benchmark$paths$results,
    benchmark$paths$metrics, file.path(figure_dir, "drought_response_phenotype_distribution.png"),
    file.path(figure_dir, "gwas_manhattan.png"), file.path(figure_dir, "gwas_qq.png"),
    file.path(figure_dir, "drought_response_gwas_overview.png")
  ))),
  stringsAsFactors = FALSE
)
write.table(pass, file.path(out_dir, "validation_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(metrics, file.path(out_dir, "truth_recovery_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(input_mode, file.path(out_dir, "input_mode.txt"))
cat("Completed Figure 15 quantitative NAM GWAS demonstration:\n", normalizePath(out_dir), "\n", sep = "")
