#!/usr/bin/env Rscript

# Figure 11: synthetic IRF6-region pedigree-liability demonstration.
# This uses allele-frequency summaries from a public reference panel to choose
# anonymous marker positions. Effect sizes, environment, and phenotype outcomes
# are deliberately synthetic and must never be interpreted as clinical risk.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args, nomatch = 0L)
  if (!i || i == length(args)) default else args[i + 1L]
}
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])), "..", ".."), mustWork = FALSE)
if (!dir.exists(root)) root <- normalizePath(".")
setwd(root)
data_dir <- get_arg("--data_dir", "data/raw/human_irf6")
out_dir <- get_arg("--out_dir", "analysis/results/human_irf6_synthetic_liability")
seed <- as.integer(get_arg("--seed", "81"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
vcf <- file.path(data_dir, "irf6_1kgp_grch38_chr1_209M_211M.vcf.gz")
gff <- file.path(data_dir, "irf6_grch38_chr1_209M_211M.gff3")
stopifnot(file.exists(vcf), file.exists(gff), nzchar(Sys.which("bcftools")))
set.seed(seed)

# Obtain biallelic SNP allele frequencies without reading 3,202 individual
# genotypes into R. These positions are markers, not asserted disease variants.
query_file <- tempfile(fileext = ".tsv")
query_command <- paste(
  "bcftools query -f",
  shQuote("%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/AF\\n"),
  shQuote(vcf), ">", shQuote(query_file)
)
if (system(query_command) != 0L) stop("bcftools query failed")
markers <- read.delim(query_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
names(markers) <- c("chromosome", "position_bp", "id", "ref", "alt", "af")
markers$af <- suppressWarnings(as.numeric(markers$af))
markers <- markers[
  is.finite(markers$af) & markers$af >= 0.10 & markers$af <= 0.40 &
    nchar(markers$ref) == 1L & nchar(markers$alt) == 1L & !grepl(",", markers$alt, fixed = TRUE),
]
if (nrow(markers) < 8L) stop("Too few common biallelic SNPs in the supplied IRF6 window")
marker_index <- unique(round(seq(1, nrow(markers), length.out = min(12L, nrow(markers)))))
markers <- markers[marker_index, ]
markers$marker_id <- paste0("IRF6_region_marker_", seq_len(nrow(markers)))
markers$synthetic_log_odds_effect <- c(0.45, rnorm(nrow(markers) - 1L, 0.12, 0.05))
write.table(markers, file.path(out_dir, "irf6_region_marker_truth.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

L <- nrow(markers)
make_unrelated <- function(n) matrix(rbinom(n * L, 2, rep(markers$af, each = n)), nrow = n, byrow = FALSE)
gamete <- function(g) rbinom(length(g), 1, g / 2)
make_child <- function(a, b) gamete(a) + gamete(b)
make_random_generation <- function(pool, n) {
  out <- matrix(NA_real_, n, L)
  for (i in seq_len(n)) {
    parents <- sample(seq_len(nrow(pool)), 2L, replace = TRUE)
    out[i, ] <- make_child(pool[parents[1L], ], pool[parents[2L], ])
  }
  out
}

# A random-mating lineage gives three comparison generations.
g0 <- make_unrelated(160L)
g1 <- make_random_generation(g0, 160L)
g2 <- make_random_generation(g1, 160L)
g3 <- make_random_generation(g2, 160L)

# Build first-cousin offspring from two siblings who each mate an unrelated
# partner. Their descendants share one grandparental pair; F = 1/16 in the
# idealized pedigree used here.
n_families <- 40L
grandparents <- make_unrelated(2L * n_families)
external <- make_unrelated(2L * n_families)
cousin_offspring <- matrix(NA_real_, n_families, L)
for (i in seq_len(n_families)) {
  grandparents_i <- grandparents[c(2L * i - 1L, 2L * i), , drop = FALSE]
  sibling_a <- make_child(grandparents_i[1L, ], grandparents_i[2L, ])
  sibling_b <- make_child(grandparents_i[1L, ], grandparents_i[2L, ])
  cousin_a <- make_child(sibling_a, external[2L * i - 1L, ])
  cousin_b <- make_child(sibling_b, external[2L * i, ])
  cousin_offspring[i, ] <- make_child(cousin_a, cousin_b)
}

summarize_group <- function(geno, group, generation, F) {
  score <- as.vector(geno %*% markers$synthetic_log_odds_effect) +
    0.75 * F + rnorm(nrow(geno), 0, 0.20)
  probability <- plogis(qlogis(0.005) + score)
  data.frame(
    sample = paste0(group, "_", seq_len(nrow(geno))), group = group,
    generation = generation, theoretical_inbreeding_F = F,
    marker_homozygosity = rowMeans(geno == 0 | geno == 2),
    synthetic_liability = score, synthetic_cleft_probability = probability,
    stringsAsFactors = FALSE
  )
}
individuals <- rbind(
  summarize_group(g0, "random_mating", 0L, 0),
  summarize_group(g1, "random_mating", 1L, 0),
  summarize_group(g2, "random_mating", 2L, 0),
  summarize_group(g3, "random_mating", 3L, 0),
  summarize_group(cousin_offspring, "first_cousin_descendant", 3L, 1 / 16)
)
write.table(individuals, file.path(out_dir, "synthetic_pedigree_liability.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
summary <- aggregate(cbind(synthetic_cleft_probability, marker_homozygosity) ~ group + generation + theoretical_inbreeding_F, individuals, mean)
summary$n <- as.integer(table(interaction(individuals$group, individuals$generation, individuals$theoretical_inbreeding_F, drop = TRUE)))
write.table(summary, file.path(out_dir, "synthetic_pedigree_liability_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

png(file.path(out_dir, "figure11_human_irf6_synthetic_liability.png"), width = 1800, height = 1300, res = 180)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.5, 1))
plot(markers$position_bp, markers$af, pch = 16, col = "#156E8A", ylim = c(0, 0.5), xlab = "Position on chr1 (GRCh38)", ylab = "1000 Genomes alternate-allele frequency", main = "A. Anonymous common markers in the IRF6 region")
abline(v = 209785617, col = "#B4422B", lty = 2, lwd = 2)
text(209785617, 0.48, "IRF6", pos = 4, col = "#B4422B")
plot(summary$generation + ifelse(summary$group == "first_cousin_descendant", 0.08, -0.08), summary$synthetic_cleft_probability, pch = 16, col = ifelse(summary$group == "random_mating", "#156E8A", "#B4422B"), xaxt = "n", xlab = "Pedigree generation", ylab = "Mean synthetic probability", main = "B. Synthetic liability by pedigree group")
axis(1, at = 0:3)
legend("topleft", legend = c("Random mating", "First-cousin descendants"), col = c("#156E8A", "#B4422B"), pch = 16, bty = "n")
boxplot(synthetic_cleft_probability ~ group, data = individuals[individuals$generation == 3L, ], col = c("#A7D3E0", "#E9B2A8"), ylab = "Synthetic probability", main = "C. Generation-3 probability distribution")
boxplot(marker_homozygosity ~ group, data = individuals[individuals$generation == 3L, ], col = c("#A7D3E0", "#E9B2A8"), ylab = "Fraction homozygous at selected markers", main = "D. Marker homozygosity")
dev.off()

metadata <- list(seed = seed, reference = "GRCh38 chr1:209,000,000-211,000,000", gene = "IRF6", marker_count = L, pedigree = "random-mating generations 0-3 plus synthetic first-cousin descendants at generation 3", limitation = "Synthetic liability effects are illustrative only; this is not a clinical predictor or a claim that selected markers cause cleft lip.")
if (requireNamespace("jsonlite", quietly = TRUE)) jsonlite::write_json(metadata, file.path(out_dir, "metadata.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Completed synthetic IRF6-region pedigree-liability demonstration:\n", normalizePath(out_dir), "\n", sep = "")
