#!/usr/bin/env Rscript

# Figure 8: annotation-aware ChIP-seq simulation
#
# This lightweight runner creates promoter, enhancer, and silencer annotations
# on the bundled genome, simulates control and treatment ChIP/input libraries,
# and exports a six-panel publication figure plus panel source data.
#
# Usage:
#   Rscript analysis/paper_fig/fig8_chipseq_simulation.R
#   Rscript analysis/paper_fig/fig8_chipseq_simulation.R \
#     --out_dir analysis/results/chipseq --seed 62 --backend native

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  index <- match(flag, args, nomatch = 0L)
  if (!index || index == length(args)) default else args[index + 1L]
}

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) {
  sub("^--file=", "", script_argument[1L])
} else {
  file.path("analysis", "paper_fig", "fig8_chipseq_simulation.R")
}
repo_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."), mustWork = TRUE
)
setwd(repo_root)

if (requireNamespace("devtools", quietly = TRUE) &&
    file.exists(file.path(repo_root, "DESCRIPTION"))) {
  devtools::load_all(repo_root, quiet = TRUE)
} else {
  library(simitall)
}

out_dir <- get_arg("--out_dir", file.path("analysis", "results", "chipseq"))
seed <- as.integer(get_arg("--seed", "62"))
backend <- get_arg("--backend", "auto")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Figure 8: annotation-aware ChIP-seq simulation ===\n")
cat("Output directory:", normalizePath(out_dir, mustWork = TRUE), "\n")
cat("Seed:", seed, "\n")
cat("Backend:", backend, "\n\n")

demo_genome <- system.file(
  "extdata", "examples", "demo_genome.fa", package = "simitall"
)
if (!nzchar(demo_genome)) stop("The bundled demo genome could not be located")
seqname <- sub("^>([^[:space:]]+).*$", "\\1", readLines(demo_genome, n = 1L))

gff_path <- file.path(out_dir, "chipseq_targets.gff3")
feature_lines <- unlist(lapply(seq_len(8L), function(i) {
  anchor <- 400L + (i - 1L) * 950L
  c(
    paste(
      seqname, "simitall", "promoter", anchor, anchor + 139L,
      ".", "+", ".", paste0("ID=promoter_", i), sep = "\t"
    ),
    paste(
      seqname, "simitall", "enhancer", anchor + 300L, anchor + 479L,
      ".", "+", ".", paste0("ID=enhancer_", i), sep = "\t"
    ),
    paste(
      seqname, "simitall", "silencer", anchor + 620L, anchor + 779L,
      ".", "-", ".", paste0("ID=silencer_", i), sep = "\t"
    )
  )
}))
writeLines(c("##gff-version 3", feature_lines), gff_path)

chip_prefix <- file.path(out_dir, "tf_chipseq")
figure_prefix <- file.path(out_dir, "figure8_chipseq")

cat("1. Simulating TF binding, matched input, and differential peaks...\n")
chip <- simulate_chipseq(
  genome_fa = demo_genome,
  annotation_gff3 = gff_path,
  out_prefix = chip_prefix,
  assay_type = "TF",
  target_features = c("promoter", "enhancer", "silencer"),
  n_peaks = 18,
  peak_width = 180,
  peak_width_sd = 35,
  peak_jitter = 30,
  enrichment_mean = 10,
  signal_fraction = 0.42,
  conditions = c("control", "treatment"),
  biological_replicates = 3,
  technical_replicates = 1,
  differential_binding_fraction = 0.45,
  differential_effect_log2 = 1.4,
  n_reads = 5000,
  input_reads = 5000,
  read_length = 50,
  fragment_mean = 160,
  fragment_sd = 25,
  fragment_min = 90,
  fragment_max = 240,
  duplicate_rate = 0.06,
  backend = backend,
  seed = seed
)

cat("2. Building publication figure and source-data tables...\n")
figure <- plot_chipseq_results(
  chip_prefix = chip_prefix,
  out_prefix = figure_prefix,
  profile_window = 750,
  profile_bin = 30,
  width = 15,
  height = 10,
  dpi = 300,
  seed = seed + 1L
)

cat("\nCompleted successfully.\n")
cat("Backend used:", attr(chip, "backend"), "\n")
cat("Peak truth:", chip$peaks, "\n")
cat("QC table:", chip$qc, "\n")
cat("Figure PDF:", figure$paths$pdf, "\n")
cat("Figure PNG:", figure$paths$png, "\n")
cat("Panel source data:", figure$paths$source_data, "\n")
