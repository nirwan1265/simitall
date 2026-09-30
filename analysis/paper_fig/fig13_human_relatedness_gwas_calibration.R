#!/usr/bin/env Rscript

# Figure 13: synthetic human cohort relatedness and GWAS-calibration benchmark.
# This creates no clinical phenotype and uses no patient data. It demonstrates
# why ancestry and close relatedness can inflate association statistics.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) { i <- match(flag, args, nomatch = 0L); if (!i || i == length(args)) default else args[i + 1L] }
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) sub("^--file=", "", script_argument[1L]) else "analysis/paper_fig/fig13_human_relatedness_gwas_calibration.R"
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
if (!dir.exists(repo_root)) repo_root <- normalizePath(".")
setwd(repo_root)

out_dir <- get_arg("--out_dir", file.path("analysis", "results", "human_relatedness_gwas_calibration"))
seed <- as.integer(get_arg("--seed", "1301"))
n_per_group <- as.integer(get_arg("--n_per_group", "120"))
n_markers <- as.integer(get_arg("--n_markers", "800"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

# Two synthetic ancestry groups differ in allele frequencies. Within each
# group we create unrelated individuals, sibling pairs, and an autozygous
# inbreeding proxy. The latter increases homozygosity at heterozygous markers
# and is explicitly synthetic; it is not a reconstruction of human pedigrees.
if (n_per_group %% 3L != 0L) stop("n_per_group must be divisible by three")
n <- 2L * n_per_group
ancestry <- character(n)
cohort <- character(n)
family <- character(n)
base_af <- runif(n_markers, 0.08, 0.42)
delta <- rnorm(n_markers, 0, 0.10)
af_a <- pmin(0.49, pmax(0.01, base_af + delta))
af_b <- pmin(0.49, pmax(0.01, base_af - delta))
G <- matrix(NA_real_, nrow = n_markers, ncol = n)
cursor <- 0L
for (g in c("group_A", "group_B")) {
  n_each <- n_per_group / 3L
  af <- if (g == "group_A") af_a else af_b
  unrelated_idx <- seq.int(cursor + 1L, cursor + n_each)
  G[, unrelated_idx] <- matrix(rbinom(n_markers * n_each, 2L, rep(af, each = n_each)), nrow = n_markers)
  ancestry[unrelated_idx] <- g; cohort[unrelated_idx] <- "unrelated"
  family[unrelated_idx] <- paste0(g, "_unrelated_", seq_len(n_each))
  cursor <- cursor + n_each

  sibling_idx <- seq.int(cursor + 1L, cursor + n_each)
  parents <- matrix(rbinom((n_each / 2L) * n_markers * 2L, 2L, rep(af, each = n_each)), nrow = n_markers)
  for (pair in seq_len(n_each / 2L)) {
    father <- parents[, 2L * pair - 1L]; mother <- parents[, 2L * pair]
    gamete <- function(x) rbinom(length(x), 1L, x / 2)
    G[, sibling_idx[2L * pair - 1L]] <- gamete(father) + gamete(mother)
    G[, sibling_idx[2L * pair]] <- gamete(father) + gamete(mother)
  }
  ancestry[sibling_idx] <- g; cohort[sibling_idx] <- "related_siblings"
  family[sibling_idx] <- rep(paste0(g, "_sibling_", seq_len(n_each / 2L)), each = 2L)
  cursor <- cursor + n_each

  inbred_idx <- seq.int(cursor + 1L, cursor + n_each)
  inbred <- matrix(rbinom(n_markers * n_each, 2L, rep(af, each = n_each)), nrow = n_markers)
  heterozygous <- inbred == 1L & matrix(runif(n_markers * n_each), nrow = n_markers) < 0.45
  inbred[heterozygous] <- 2L * rbinom(sum(heterozygous), 1L, rep(af, length.out = sum(heterozygous)))
  G[, inbred_idx] <- inbred
  ancestry[inbred_idx] <- g; cohort[inbred_idx] <- "inbred_autozygous_proxy"
  family[inbred_idx] <- paste0(g, "_autozygous_", seq_len(n_each))
  cursor <- cursor + n_each
}
samples <- sprintf("human_%03d", seq_len(n))
markers <- data.frame(marker_id = sprintf("chr10_marker_%04d", seq_len(n_markers)), seqname = "chr10", pos = sort(sample.int(10000000L, n_markers)), ref = "A", alt = "G", G, check.names = FALSE, stringsAsFactors = FALSE)
names(markers)[6:ncol(markers)] <- samples
geno_file <- file.path(out_dir, "synthetic_human_chr10_genotypes.tsv")
write.table(markers, geno_file, sep = "\t", row.names = FALSE, quote = FALSE)

# No marker causes this trait. A group-level environmental shift is deliberately
# inserted so an unadjusted GWAS has confounding to discover.
trait <- 1.2 * as.integer(ancestry == "group_B") + rnorm(n, 0, 1)
pheno <- data.frame(sample = samples, ancestry = ancestry, cohort = cohort, family = family, marker_homozygosity = colMeans(G == 0L | G == 2L), synthetic_null_trait = trait, stringsAsFactors = FALSE)
write.table(pheno, file.path(out_dir, "synthetic_null_phenotype.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# PCs and kinship eigenvectors are derived from all markers. The adjusted scan
# uses both plus known family labels, providing an explicit family/kinship-aware
# fixed-effect calibration demonstration (not a replacement for an LMM).
Z <- scale(t(G), center = TRUE, scale = FALSE)
pc <- stats::prcomp(Z, rank. = 6L)$x
K <- tcrossprod(Z) / ncol(Z)
kin_eigen <- eigen(K, symmetric = TRUE)$vectors[, 1:6, drop = FALSE]
association <- function(adjusted = FALSE) {
  p <- numeric(n_markers)
  for (j in seq_len(n_markers)) {
    dat <- data.frame(y = trait, dosage = G[j, ], PC1 = pc[, 1], PC2 = pc[, 2], K1 = kin_eigen[, 1], K2 = kin_eigen[, 2], family = family)
    fit <- if (adjusted) stats::lm(y ~ dosage + PC1 + PC2 + K1 + K2 + family, data = dat) else stats::lm(y ~ dosage, data = dat)
    p[j] <- summary(fit)$coefficients["dosage", "Pr(>|t|)"]
  }
  p
}
cat("1. Running unadjusted and PC-plus-kinship-adjusted null GWAS scans...\n")
p_naive <- association(FALSE)
p_adjusted <- association(TRUE)
lambda <- function(p) stats::median(stats::qchisq(1 - p, 1L), na.rm = TRUE) / stats::qchisq(0.5, 1L)
results <- data.frame(marker_id = markers$marker_id, seqname = "chr10", pos = markers$pos, naive_p_value = p_naive, pc_kinship_p_value = p_adjusted, stringsAsFactors = FALSE)
write.table(results, file.path(out_dir, "gwas_calibration_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
metrics <- data.frame(method = c("unadjusted", "PC_plus_kinship_eigenvectors"), genomic_inflation_lambda = c(lambda(p_naive), lambda(p_adjusted)), false_positive_rate_0_05 = c(mean(p_naive < 0.05), mean(p_adjusted < 0.05)), stringsAsFactors = FALSE)
write.table(metrics, file.path(out_dir, "gwas_calibration_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

figure <- file.path(out_dir, "figure13_human_relatedness_gwas_calibration.png")
png(figure, width = 2000, height = 1500, res = 180)
par(mfrow = c(2, 2), mar = c(4.5, 4.7, 3.2, 1))
plot(pc[, 1], pc[, 2], col = ifelse(ancestry == "group_A", "#2E7D6A", "#C9653E"), pch = ifelse(cohort == "related_siblings", 16, ifelse(cohort == "inbred_autozygous_proxy", 17, 1)), xlab = "Genotype PC1", ylab = "Genotype PC2", main = "A. Ancestry and relatedness cohorts")
legend("topright", legend = c("group A", "group B"), col = c("#2E7D6A", "#C9653E"), pch = 16, bty = "n")
boxplot(marker_homozygosity ~ cohort, data = pheno, las = 2, col = c("#BED8E6", "#A7D3E0", "#E9B2A8"), ylab = "Marker homozygosity", main = "B. Unrelated, sibling, and inbred cohorts")
qq <- function(p, col, title) { obs <- -log10(sort(p)); exp <- -log10(ppoints(length(p))); plot(exp, obs, pch = 16, cex = 0.55, col = col, xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", main = title); abline(0, 1, lty = 2, col = "#444444") }
qq(p_naive, "#B4422B", paste0("C. Unadjusted null GWAS; lambda=", sprintf("%.2f", lambda(p_naive))))
qq(p_adjusted, "#315C6B", paste0("D. PC + family/kinship proxy; lambda=", sprintf("%.2f", lambda(p_adjusted))))
dev.off()

pass <- data.frame(
  samples = n, markers = n_markers, unrelated_individuals = sum(cohort == "unrelated"), related_sibling_individuals = sum(cohort == "related_siblings"), inbred_proxy_individuals = sum(cohort == "inbred_autozygous_proxy"),
  naive_lambda = lambda(p_naive), adjusted_lambda = lambda(p_adjusted),
  adjustment_reduces_inflation = lambda(p_adjusted) < lambda(p_naive),
  figure_exists = file.exists(figure),
  limitation = "Synthetic null phenotype with programmed ancestry confounding; this is a calibration illustration, not human disease inference.",
  stringsAsFactors = FALSE
)
write.table(pass, file.path(out_dir, "validation_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("Completed Figure 13 demonstration:\n", normalizePath(out_dir), "\n", sep = "")
