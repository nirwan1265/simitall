#!/usr/bin/env Rscript

# Figure 14: synthetic NAM binary-trait GWAS on the bundled toy maize chr10.
#
# This is a mechanism demonstration. It makes synthetic case/control outcomes
# from one causal marker dosage and fits a simple logistic scan. NAM founder
# family labels are retained in the metadata but are not shown or used unless a
# user explicitly asks for structure-adjusted analysis. This is not a disease
# model or a replacement for a binary mixed-model GWAS in a real related cohort.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args, nomatch = 0L)
  if (!i || i == length(args)) default else args[i + 1L]
}
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) sub("^--file=", "", script_argument[1L]) else "analysis/paper_fig/fig14_maize_nam_binary_gwas.R"
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
if (!dir.exists(repo_root)) repo_root <- normalizePath(".")
setwd(repo_root)
if (requireNamespace("devtools", quietly = TRUE)) devtools::load_all(repo_root, quiet = TRUE) else library(simitall)

out_dir <- get_arg("--out_dir", file.path("analysis", "results", "maize_nam_binary_gwas"))
seed <- as.integer(get_arg("--seed", "1401"))
n_lines <- as.integer(get_arg("--n_lines", "100"))
target_prevalence <- as.numeric(get_arg("--prevalence", "0.30"))
input_mode <- get_arg("--input_mode", "package_demo")
provided_panel <- get_arg("--haplotype_fa", "")
provided_map <- get_arg("--recomb_map", "")
if (!is.finite(target_prevalence) || target_prevalence <= 0 || target_prevalence >= 1) {
  stop("--prevalence must be between zero and one")
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
  headers <- grepl("^>", panel_lines)
  panel_lines[headers] <- paste0(panel_lines[headers], "|chr10")
  writeLines(panel_lines, panel)
  map <- file.path(input_dir, "synthetic_nam_chr10_map.tsv")
  synthetic_map <- data.frame(
    chromosome = "chr10", pos_bp = c(1L, 1000L, 1800L, 2600L, 4000L),
    cM = c(0, 1, 1.1, 1.2, 2.6)
  )
  write.table(synthetic_map, map, sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  if (!nzchar(provided_panel) || !nzchar(provided_map) ||
      !file.exists(provided_panel) || !file.exists(provided_map)) {
    stop("Provided mode requires existing --haplotype_fa and --recomb_map files")
  }
  panel <- normalizePath(provided_panel)
  map <- normalizePath(provided_map)
}
prefix <- file.path(out_dir, "nam_binary")
vcf <- paste0(prefix, ".vcf")

cat("1. Simulating a synthetic NAM population with ", n_lines, " lines...\n", sep = "")
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
if (!length(eligible)) stop("The synthetic NAM panel did not contain an eligible causal marker")
# Choose a common marker so a 100-line binary demonstration has usable power.
causal <- eligible[which.min(abs(maf[eligible] - 0.35))]
effect <- 3.5
genetic_score <- as.numeric(G[causal, ] * effect)
intercept <- uniroot(
  function(value) mean(stats::plogis(value + genetic_score)) - target_prevalence,
  interval = c(-20, 20)
)$root
probability <- stats::plogis(intercept + genetic_score)
case <- stats::rbinom(length(samples), 1L, probability)
if (length(unique(case)) < 2L) stop("Synthetic case/control outcome has one class; change --seed or --prevalence")
phenotype <- data.frame(
  sample = samples, family = metadata$family, case = case,
  programmed_probability = probability, genetic_score = genetic_score,
  stringsAsFactors = FALSE
)
truth <- data.frame(
  marker_id = variants$id[causal], seqname = variants$seqname[causal],
  pos = variants$pos[causal], beta_log_odds = effect,
  odds_ratio = exp(effect), stringsAsFactors = FALSE
)
phenotype_file <- file.path(out_dir, "binary_phenotypes.tsv")
truth_file <- file.path(out_dir, "binary_trait_truth.tsv")
write.table(phenotype, phenotype_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(truth, truth_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("2. Fitting a simple logistic GWAS...\n")
gwas_prefix <- file.path(out_dir, "binary_gwas")
gwas <- analyze_binary_gwas(
  vcf, phenotype, gwas_prefix, trait = "case", min_maf = 0.05
)
benchmark <- benchmark_gwas(gwas$paths$results, truth, gwas_prefix, window_bp = 150L)
plot_gwas_results(
  benchmark$paths$results, file.path(figure_dir, "binary_gwas_overview"),
  individual_panels_dir = figure_dir
)
plot_trait_diagnostics(
  phenotype, "case", figure_dir, file_stem = "binary_phenotype_counts"
)

results <- read.delim(benchmark$paths$results, check.names = FALSE)

metrics <- read.delim(benchmark$paths$metrics, check.names = FALSE)
pass <- data.frame(
  nam_lines = n_lines,
  cases = sum(case), controls = sum(case == 0L), causal_markers = nrow(truth),
  gwas_markers = nrow(results),
  phenotype_is_binary = all(phenotype$case %in% c(0L, 1L)),
  all_core_outputs_exist = all(file.exists(c(
    vcf, phenotype_file, truth_file, gwas$paths$results, benchmark$paths$results,
    benchmark$paths$metrics, file.path(figure_dir, "binary_phenotype_counts.png"),
    file.path(figure_dir, "gwas_manhattan.png"), file.path(figure_dir, "gwas_qq.png"),
    file.path(figure_dir, "binary_gwas_overview.png")
  ))),
  stringsAsFactors = FALSE
)
write.table(pass, file.path(out_dir, "validation_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(metrics, file.path(out_dir, "truth_recovery_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(input_mode, file.path(out_dir, "input_mode.txt"))
cat("Completed Figure 14 demonstration:\n", normalizePath(out_dir), "\n", sep = "")
