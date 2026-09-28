#!/usr/bin/env Rscript

# Figure 10: synthetic biparental BC2S7 allele-frequency demonstration
#
# This example uses two bundled synthetic founders, labels their single
# chromosome as chr10, and produces an F1 -> BC2 to parent 1 -> S7 population.
# It is an implementation example, not a simulation of named maize NAM lines.
#
# Usage:
#   Rscript analysis/paper_fig/fig10_biparental_bc2s7_allele_frequency.R
#   Rscript analysis/paper_fig/fig10_biparental_bc2s7_allele_frequency.R \
#     --out_dir analysis/results/biparental_bc2s7 --seed 72

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  index <- match(flag, args, nomatch = 0L)
  if (!index || index == length(args)) default else args[index + 1L]
}

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) {
  sub("^--file=", "", script_argument[1L])
} else {
  file.path("analysis", "paper_fig", "fig10_biparental_bc2s7_allele_frequency.R")
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
setwd(repo_root)

if (requireNamespace("devtools", quietly = TRUE) && file.exists(file.path(repo_root, "DESCRIPTION"))) {
  devtools::load_all(repo_root, quiet = TRUE)
} else {
  library(simitall)
}

out_dir <- get_arg("--out_dir", file.path("analysis", "results", "biparental_bc2s7"))
seed <- as.integer(get_arg("--seed", "72"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Figure 10: synthetic biparental BC2S7 allele frequencies ===\n")
cat("Output directory:", normalizePath(out_dir, mustWork = TRUE), "\n")
cat("Seed:", seed, "\n\n")

# Copy two synthetic founders and give the demonstration chromosome an explicit
# chr10 label. A real NAM analysis must instead supply real founder haplotypes
# and a compatible chromosome-10 recombination map.
demo_panel <- system.file("extdata", "panels", "demo_panel.fa", package = "simitall")
if (!nzchar(demo_panel)) stop("The bundled demo panel could not be located")
panel_lines <- readLines(demo_panel)
header_index <- grep("^>", panel_lines)
panel_lines[header_index] <- paste0(panel_lines[header_index], "|chr10")
panel <- file.path(out_dir, "synthetic_founders_chr10.fa")
writeLines(panel_lines, panel)

prefix <- file.path(out_dir, "bc2s7")
vcf <- paste0(prefix, ".vcf")

cat("1. Simulating F1 -> BC2 to parent 1 -> S7...\n")
simulate_breeding(
  haplotype_fa = panel,
  out_prefix = prefix,
  parents = "hap1,hap2",
  sequence = "F1,BC:P1:2,SELF:7",
  n_offspring = 200,
  vcf_out = vcf,
  ancestry_out = paste0(prefix, ".ancestry.tsv"),
  breakpoints_out = paste0(prefix, ".breakpoints.tsv"),
  recomb_map_out = paste0(prefix, ".recombination_map.tsv"),
  seed = seed
)

cat("2. Calculating alternate-allele frequency from the final VCF...\n")
vcf_lines <- readLines(vcf, warn = FALSE)
header <- strsplit(sub("^#", "", vcf_lines[grepl("^#CHROM", vcf_lines)][1L]), "\t", fixed = TRUE)[[1L]]
records <- vcf_lines[!grepl("^#", vcf_lines)]
if (!length(records)) stop("The breeding run produced no polymorphic VCF records")
fields <- strsplit(records, "\t", fixed = TRUE)

as_dosage <- function(genotype) {
  allele_text <- sub(":.*$", "", genotype)
  alleles <- strsplit(allele_text, "[/|]")[[1L]]
  if (any(alleles == ".")) return(NA_real_)
  sum(as.numeric(alleles))
}

allele_frequency <- vapply(fields, function(x) {
  dosage <- vapply(x[10:length(x)], as_dosage, numeric(1))
  mean(dosage, na.rm = TRUE) / 2
}, numeric(1))

frequency <- data.frame(
  chromosome = vapply(fields, `[[`, character(1), 1L),
  position_bp = as.integer(vapply(fields, `[[`, character(1), 2L)),
  alternate_allele_frequency = allele_frequency,
  stringsAsFactors = FALSE
)
frequency_file <- file.path(out_dir, "bc2s7_allele_frequencies.tsv")
write.table(frequency, frequency_file, sep = "\t", row.names = FALSE, quote = FALSE)

summary <- data.frame(
  population = "synthetic_biparental_BC2S7",
  founders = "hap1,hap2",
  chromosome = "chr10",
  final_lines = 200L,
  sequence = "F1,BC:P1:2,SELF:7",
  polymorphic_sites = nrow(frequency),
  mean_alternate_allele_frequency = mean(frequency$alternate_allele_frequency),
  expected_donor_allele_frequency_without_selection = 0.125,
  seed = seed,
  stringsAsFactors = FALSE
)
summary_file <- file.path(out_dir, "figure10_biparental_bc2s7.summary.tsv")
write.table(summary, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("3. Writing allele-frequency figure...\n")
figure_file <- file.path(out_dir, "figure10_biparental_bc2s7_allele_frequency.png")
png(figure_file, width = 1600, height = 1000, res = 180)
plot(
  frequency$position_bp, frequency$alternate_allele_frequency,
  pch = 16, cex = 0.75, col = grDevices::adjustcolor("#156E8A", alpha.f = 0.75),
  xlab = "Position on synthetic chr10 (bp)", ylab = "Alternate-allele frequency",
  ylim = c(0, 1),
  main = "Synthetic biparental BC2S7 population: final allele frequencies"
)
abline(
  h = summary$expected_donor_allele_frequency_without_selection,
  col = "#B4422B", lty = 2, lwd = 2
)
legend(
  "topright",
  legend = c("Observed polymorphic-site frequency", "Expected donor frequency without selection (1/8)"),
  col = c("#156E8A", "#B4422B"), pch = c(16, NA), lty = c(NA, 2), lwd = c(NA, 2), bty = "n"
)
dev.off()

cat("Completed successfully.\n")
cat("Figure:", normalizePath(figure_file), "\n")
cat("Allele frequencies:", normalizePath(frequency_file), "\n")
cat("Summary:", normalizePath(summary_file), "\n")
