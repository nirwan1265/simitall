#!/usr/bin/env Rscript

# Figure 12: single-chromosome maize NAM GxE, RNA-seq, GWAS, and selection.
#
# The bundled panel is synthetic and has eight founder haplotypes on chr10.
# Supply --gff3 with a compatible chromosome-10 annotation to anchor the focal
# gene to a real annotation. The default generated gene is a demonstration
# feature, not a claim about a named maize drought gene.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args, nomatch = 0L)
  if (!i || i == length(args)) default else args[i + 1L]
}
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) sub("^--file=", "", script_argument[1L]) else "analysis/paper_fig/fig12_maize_nam_chr10_gxe_multiomics.R"
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
if (!dir.exists(repo_root)) repo_root <- normalizePath(".")
setwd(repo_root)
if (requireNamespace("devtools", quietly = TRUE)) devtools::load_all(repo_root, quiet = TRUE) else library(simitall)

out_dir <- get_arg("--out_dir", file.path("analysis", "results", "maize_nam_chr10_gxe_multiomics"))
seed <- as.integer(get_arg("--seed", "1201"))
n_lines <- as.integer(get_arg("--n_lines", "300"))
user_gff3 <- get_arg("--gff3", "")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(out_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

panel <- system.file("extdata", "panels", "demo_maize_nam_chr10.fa", package = "simitall")
map <- system.file("extdata", "maps", "demo_maize_nam_chr10_map.tsv", package = "simitall")
stopifnot(nzchar(panel), nzchar(map), file.exists(panel), file.exists(map))
prefix <- file.path(out_dir, "nam_chr10")
vcf <- paste0(prefix, ".vcf")

cat("1. Simulating chr10 NAM population with ", n_lines, " lines...\n", sep = "")
simulate_breeding(
  haplotype_fa = panel, out_prefix = prefix, scheme = "NAM", n_offspring = n_lines,
  recomb_map_in = map, vcf_out = vcf,
  ancestry_out = paste0(prefix, ".ancestry.tsv"),
  breakpoints_out = paste0(prefix, ".breakpoints.tsv"), seed = seed
)
geno <- simitall:::.simitall_read_genotypes(vcf)
G <- geno$genotype
G[!is.finite(G)] <- rowMeans(G, na.rm = TRUE)[row(G)[!is.finite(G)]]
variants <- geno$variants
samples <- colnames(G)
metadata <- read.delim(paste0(prefix, ".meta.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
metadata <- metadata[match(samples, metadata$sample), , drop = FALSE]

# Resolve or create the single focal gene annotation.
focal_pos <- variants$pos[which.min(abs(variants$pos - median(variants$pos)))]
if (nzchar(user_gff3)) {
  if (!file.exists(user_gff3)) stop("--gff3 does not exist: ", user_gff3)
  annotation <- read.delim(user_gff3, comment.char = "#", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  if (ncol(annotation) < 9L) stop("--gff3 must have nine GFF3 columns")
  names(annotation)[1:9] <- c("seqname", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  genes <- annotation[annotation$type %in% c("gene", "mRNA"), , drop = FALSE]
  genes <- genes[genes$seqname == "chr10", , drop = FALSE]
  if (!nrow(genes)) stop("--gff3 has no chr10 gene or mRNA features")
  midpoint <- (as.numeric(genes$start) + as.numeric(genes$end)) / 2
  focal <- genes[which.min(abs(midpoint - focal_pos)), , drop = FALSE]
  focal_gene <- sub(".*(?:ID|Name)=([^;]+).*", "\\1", focal$attributes[1L])
  focal_start <- as.integer(focal$start[1L]); focal_end <- as.integer(focal$end[1L])
  annotation_source <- normalizePath(user_gff3)
} else {
  focal_gene <- "SIM_DROUGHT_RESPONSE_CHR10"
  focal_start <- max(1L, focal_pos - 100L); focal_end <- focal_pos + 100L
  annotation_source <- "synthetic demonstration annotation"
  gene_centers <- unique(round(seq(min(variants$pos), max(variants$pos), length.out = 12L)))
  focal_gene_index <- which.min(abs(gene_centers - focal_pos))
  gene_ids <- paste0("SIM_BACKGROUND_GENE_", sprintf("%02d", seq_along(gene_centers)))
  gene_ids[focal_gene_index] <- focal_gene
  gene_rows <- vapply(seq_along(gene_centers), function(i) {
    start <- max(1L, gene_centers[i] - 80L)
    end <- gene_centers[i] + 80L
    paste("chr10", "simitall", "gene", start, end, ".", "+", ".", paste0("ID=", gene_ids[i], ";Name=", gene_ids[i]), sep = "\t")
  }, character(1))
  writeLines(c("##gff-version 3", gene_rows), file.path(out_dir, "synthetic_chr10_focal_gene.gff3"))
}
focal_marker <- which.min(abs(variants$pos - focal_pos))
polygenic_markers <- setdiff(seq_len(nrow(G)), focal_marker)
polygenic_markers <- sample(polygenic_markers, min(6L, length(polygenic_markers)))
effects <- c(1.15, stats::rnorm(length(polygenic_markers), 0, 0.18))
truth <- data.frame(marker_id = variants$id[c(focal_marker, polygenic_markers)], seqname = variants$seqname[c(focal_marker, polygenic_markers)], pos = variants$pos[c(focal_marker, polygenic_markers)], effect = effects, role = c("focal_drought_GxE", rep("polygenic", length(polygenic_markers))), focal_gene = focal_gene, stringsAsFactors = FALSE)
write.table(truth, file.path(out_dir, "trait_truth.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# The same NAM lines are observed in both environments for the phenotype study.
cat("2. Simulating irrigated and drought phenotypes with a focal GxE locus...\n")
trait_rows <- do.call(rbind, lapply(c("irrigated", "drought"), function(environment) {
  drought <- as.integer(environment == "drought")
  genetic <- as.numeric(t(G[c(focal_marker, polygenic_markers), , drop = FALSE]) %*% effects)
  trait <- 10 + genetic + drought * (-1.1 + 1.35 * G[focal_marker, ]) + stats::rnorm(n_lines, 0, 0.8)
  data.frame(sample = samples, family = metadata$family, environment = environment, drought_index = drought, grain_yield_index = trait, focal_marker_dosage = G[focal_marker, ], stringsAsFactors = FALSE)
}))
write.table(trait_rows, file.path(out_dir, "gxe_phenotypes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
plot_trait_diagnostics(
  trait_rows, trait = "grain_yield_index", group = "environment",
  out_dir = figure_dir, file_stem = "gxe_phenotype_distribution"
)

# GWAS is performed in drought, where the programmed focal interaction is active.
drought_pheno <- trait_rows[trait_rows$environment == "drought", c("sample", "family", "grain_yield_index")]
cat("3. Running kinship-aware drought GWAS...\n")
gwas_prefix <- file.path(out_dir, "drought_gwas")
analyze_gwas(vcf, drought_pheno, gwas_prefix, trait = "grain_yield_index", fixed_effects = "family", n_pcs = 2L, min_maf = 0.02)
benchmark_gwas(paste0(gwas_prefix, ".gwas.tsv"), truth, file.path(out_dir, "drought_gwas"), window_bp = 150L)
gwas <- read.delim(paste0(gwas_prefix, ".gwas.tsv"), check.names = FALSE)
# The package plotting API writes a Manhattan and QQ diagnostic pair from the
# saved benchmark table, including the programmed causal window highlight.
gwas_diagnostics <- plot_gwas_results(
  paste0(gwas_prefix, ".benchmark.tsv"),
  file.path(figure_dir, "gwas_overview")
)
export_simitall_plot_panels(
  gwas_diagnostics, figure_dir, c("gwas_manhattan", "gwas_qq")
)

# Assign each genotype to one experimental condition for a compact RNA-seq run.
rna_metadata <- metadata[, c("sample", "family"), drop = FALSE]
rna_metadata$condition <- rep(c("irrigated", "drought"), length.out = nrow(rna_metadata))
rna_metadata$batch <- rep(c("batch1", "batch2"), length.out = nrow(rna_metadata))
rna_metadata_file <- file.path(out_dir, "rnaseq_sample_metadata_input.tsv")
write.table(rna_metadata, rna_metadata_file, sep = "\t", row.names = FALSE, quote = FALSE)
annotation_for_rna <- if (nzchar(user_gff3)) user_gff3 else file.path(out_dir, "synthetic_chr10_focal_gene.gff3")
cat("4. Generating bulk RNA-seq counts and cis/trans eQTL truth from NAM genotypes...\n")
rna_prefix <- file.path(out_dir, "nam_drought_rnaseq")
simulate_rnaseq_from_gwas(vcf, rna_prefix, annotation_gff3 = annotation_for_rna, sample_metadata = rna_metadata_file, n_genes = if (nzchar(user_gff3)) 20L else 12L, n_cis_eqtl = 5L, n_trans_eqtl = 4L, condition_effect_fraction = 0.5, gxe_fraction = 0.6, seed = seed + 1L)

# Genomic selection uses drought phenotypes and preserves NAM-family limits.
cat("5. Fitting GBLUP and selecting diverse parents...\n")
selection <- run_genomic_selection(vcf, drought_pheno, file.path(out_dir, "drought_selection"), trait = "grain_yield_index", model = "gblup", n_parents = 12L, n_crosses = 18L, diversity_penalty = 0.25, family_data = metadata[, c("sample", "family")], seed = seed + 2L)

# Four validation panels: GxE, GWAS evidence, expression separation, and selection gain.
counts <- read.delim(paste0(rna_prefix, ".counts.tsv"), check.names = FALSE, row.names = 1)
log_counts <- log2(as.matrix(counts) + 1)
pca <- stats::prcomp(t(log_counts), scale. = TRUE)
rna_meta <- read.delim(paste0(rna_prefix, ".sample_metadata.tsv"), check.names = FALSE)
pca_df <- data.frame(sample = rownames(pca$x), PC1 = pca$x[, 1L], PC2 = if (ncol(pca$x) > 1L) pca$x[, 2L] else 0, stringsAsFactors = FALSE)
pca_df <- merge(pca_df, rna_meta[, c("sample", "condition")], by = "sample", all.x = TRUE)

rna_diagnostics <- plot_rnaseq_diagnostics(
  counts, rna_meta, out_dir = figure_dir, fdr_threshold = 0.10
)
figure <- file.path(figure_dir, "maize_nam_chr10_gxe_multiomics_overview.png")
png(figure, width = 2000, height = 1500, res = 180)
par(mfrow = c(2, 2), mar = c(4.5, 4.7, 3.2, 1))
boxplot(grain_yield_index ~ environment, data = trait_rows, col = c("#9BC9B5", "#D99062"), ylab = "Synthetic grain-yield index", main = "A. Programmed genotype-by-environment design")
points(jitter(as.integer(factor(trait_rows$environment)), amount = 0.10), trait_rows$grain_yield_index, pch = 16, cex = 0.35, col = grDevices::adjustcolor("#17324D", 0.25))
plot(gwas$pos, -log10(pmax(gwas$p_value, .Machine$double.xmin)), pch = 16, col = "#315C6B", xlab = "chr10 position (bp)", ylab = expression(-log[10](p)), main = "B. Drought GWAS; focal-GxE truth in orange")
points(gwas$pos[gwas$marker_id %in% truth$marker_id], -log10(pmax(gwas$p_value[gwas$marker_id %in% truth$marker_id], .Machine$double.xmin)), pch = 16, col = "#D97706")
abline(v = variants$pos[focal_marker], lty = 2, col = "#B4422B")
plot(pca_df$PC1, pca_df$PC2, pch = 16, col = ifelse(pca_df$condition == "drought", "#D99062", "#2E7D6A"), xlab = "Expression PC1", ylab = "Expression PC2", main = "C. RNA-seq counts from the same NAM genotypes")
legend("topright", legend = c("Irrigated", "Drought"), col = c("#2E7D6A", "#D99062"), pch = 16, bty = "n")
boxplot(selection$predictions$gebv, selection$selected$gebv, names = c("All candidates", "Selected parents"), col = c("#BED8E6", "#E9B2A8"), ylab = "GBLV / predicted value", main = "D. Genomic selection differential")
dev.off()
selection_figure <- file.path(figure_dir, "genomic_selection_parent_values.png")
png(selection_figure, width = 1400, height = 900, res = 180)
boxplot(selection$predictions$gebv, selection$selected$gebv,
  names = c("All candidates", "Selected parents"),
  col = c("#BED8E6", "#E9B2A8"), ylab = "GEBV / predicted value",
  main = "Genomic-selection differential"
)
dev.off()

metrics <- read.delim(paste0(out_dir, "/drought_gwas.benchmark_metrics.tsv"), check.names = FALSE)
pass <- data.frame(
  nam_lines = n_lines,
  focal_gene = focal_gene,
  focal_gene_annotation = annotation_source,
  phenotype_rows = nrow(trait_rows),
  gwas_markers = nrow(gwas),
  rnaseq_genes = nrow(counts),
  selected_parents = nrow(selection$selected),
  figure_exists = file.exists(figure),
  all_core_outputs_exist = all(file.exists(c(vcf, paste0(prefix, ".ancestry.tsv"), paste0(prefix, ".breakpoints.tsv"), file.path(out_dir, "gxe_phenotypes.tsv"), paste0(gwas_prefix, ".gwas.tsv"), paste0(rna_prefix, ".counts.tsv"), figure, selection_figure, rna_diagnostics$paths$library_size, rna_diagnostics$paths$pca, rna_diagnostics$paths$volcano, rna_diagnostics$paths$deg, file.path(figure_dir, "gwas_manhattan.png"), file.path(figure_dir, "gwas_qq.png")))),
  stringsAsFactors = FALSE
)
write.table(pass, file.path(out_dir, "validation_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("Completed Figure 12 demonstration:\n", normalizePath(out_dir), "\n", sep = "")
