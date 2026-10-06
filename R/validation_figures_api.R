.simitall_require_figure_packages <- function(feature) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for ", feature, ".", call. = FALSE)
  }
}

.simitall_figure_data <- function(data, label) {
  if (is.character(data) && length(data) == 1L) {
    if (!file.exists(data)) stop(label, " does not exist: ", data)
    data <- utils::read.delim(data, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(data)) stop(label, " must be a data frame or TSV path")
  data
}

#' Export individual panels from a simitall result figure
#'
#' Save each component of a figure returned by [plot_gwas_results()],
#' [plot_rnaseq_eqtl_results()], [plot_chipseq_results()], or
#' [plot_scrnaseq_results()] as its own PNG and PDF. This lets a workflow keep
#' a compact multi-panel overview while also providing individual plots for a
#' report, Shiny download, or agent-generated result folder.
#'
#' @param figure_result Result object returned by a simitall plotting function.
#' @param out_dir Directory for individual figure files, normally
#'   `results/<run>/figures`.
#' @param panel_names Optional file stems. Supply one name per panel.
#' @param width,height Dimensions in inches for each panel.
#' @param dpi PNG resolution.
#'
#' @return A named list of PNG and PDF paths, invisibly.
#' @examples
#' \dontrun{
#' gwas_plot <- plot_gwas_results("results/gwas.benchmark.tsv", "results/gwas")
#' export_simitall_plot_panels(
#'   gwas_plot, "results/figures", c("gwas_manhattan", "gwas_qq")
#' )
#' }
#' @export
export_simitall_plot_panels <- function(
    figure_result,
    out_dir,
    panel_names = NULL,
    width = 7,
    height = 5,
    dpi = 300) {
  .simitall_require_figure_packages("individual figure export")
  if (!is.list(figure_result)) {
    stop("figure_result must be the list returned by a simitall plotting function")
  }
  panels <- figure_result$panels
  if (is.null(panels)) {
    panels <- Filter(Negate(is.null), list(
      manhattan = figure_result$manhattan,
      qq = figure_result$qq
    ))
  }
  if (!length(panels)) stop("No individual panels were found in figure_result")
  if (is.null(panel_names)) {
    panel_names <- names(panels)
    if (is.null(panel_names) || any(!nzchar(panel_names))) {
      panel_names <- paste0("panel_", seq_along(panels))
    }
  }
  if (length(panel_names) != length(panels) || any(!nzchar(panel_names))) {
    stop("panel_names must contain one non-empty name per panel")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- lapply(seq_along(panels), function(i) {
    prefix <- file.path(out_dir, panel_names[i])
    png <- paste0(prefix, ".png")
    pdf <- paste0(prefix, ".pdf")
    ggplot2::ggsave(png, panels[[i]], width = width, height = height, dpi = dpi)
    ggplot2::ggsave(pdf, panels[[i]], width = width, height = height)
    c(png = png, pdf = pdf)
  })
  names(paths) <- panel_names
  invisible(paths)
}

#' Plot a quantitative phenotype distribution
#'
#' Create a standalone phenotype QC figure. Binary 0/1 traits are displayed as
#' case/control counts (or groupwise case fractions), continuous traits use a
#' histogram/density plot without groups and a violin/box/jitter plot with
#' groups. This is suitable for simulated traits and real numeric phenotype
#' tables.
#'
#' @param phenotype Data frame or TSV containing the trait.
#' @param trait Numeric trait column.
#' @param out_dir Directory for the figure, normally `results/<run>/figures`.
#' @param group Optional categorical grouping column such as environment,
#'   treatment, family, or population.
#' @param file_stem Output file stem.
#' @param width,height,dpi Figure dimensions and PNG resolution.
#'
#' @return Invisibly returns the plot and PNG/PDF paths.
#' @examples
#' \dontrun{
#' plot_trait_diagnostics(
#'   "results/phenotypes.tsv", trait = "grain_yield", group = "environment",
#'   out_dir = "results/figures"
#' )
#' }
#' @export
plot_trait_diagnostics <- function(
    phenotype,
    trait,
    out_dir,
    group = NULL,
    file_stem = "phenotype_distribution",
    width = 7,
    height = 5,
    dpi = 300) {
  .simitall_require_figure_packages("phenotype diagnostics")
  data <- .simitall_figure_data(phenotype, "phenotype")
  if (!trait %in% names(data) || !is.numeric(data[[trait]])) {
    stop("trait must identify a numeric phenotype column")
  }
  data <- data[is.finite(data[[trait]]), , drop = FALSE]
  if (!nrow(data)) stop("trait contains no finite values")
  if (!is.null(group) && !group %in% names(data)) stop("group column not found: ", group)
  is_binary <- all(data[[trait]] %in% c(0, 1))
  if (is_binary && is.null(group)) {
    counts <- data.frame(
      outcome = factor(c("Control", "Case"), levels = c("Control", "Case")),
      count = c(sum(data[[trait]] == 0), sum(data[[trait]] == 1))
    )
    plot <- ggplot2::ggplot(counts, ggplot2::aes(x = outcome, y = count, fill = outcome)) +
      ggplot2::geom_col(show.legend = FALSE, width = 0.65) +
      ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.4, size = 4) +
      ggplot2::scale_fill_manual(values = c("Control" = "#315C6B", "Case" = "#D97706")) +
      ggplot2::labs(title = "Binary phenotype counts", x = NULL, y = "Samples")
  } else if (is_binary) {
    rates <- stats::aggregate(data[[trait]], list(group = data[[group]]), mean)
    names(rates) <- c("group", "case_fraction")
    plot <- ggplot2::ggplot(rates, ggplot2::aes(x = group, y = case_fraction, fill = group)) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::geom_hline(yintercept = mean(data[[trait]]), linetype = 2, color = "#B4422B") +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(title = "Binary phenotype fraction by group", x = group, y = "Case fraction")
  } else if (is.null(group)) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[trait]])) +
      ggplot2::geom_histogram(bins = 30, fill = "#2E7D6A", color = "white") +
      # `..count..` remains compatible with the older ggplot2 versions
      # commonly bundled with local R installations.
      ggplot2::geom_density(ggplot2::aes(y = ..count..), color = "#B4422B", linewidth = 0.9) +
      ggplot2::labs(title = "Phenotype distribution", x = trait, y = "Samples")
  } else {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[group]], y = .data[[trait]], fill = .data[[group]])) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
      ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
      ggplot2::geom_jitter(width = 0.10, alpha = 0.35, size = 1) +
      ggplot2::labs(title = "Phenotype distribution by group", x = group, y = trait, fill = NULL) +
      ggplot2::theme(legend.position = "none")
  }
  plot <- plot + ggplot2::theme_minimal(base_size = 11)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list(png = file.path(out_dir, paste0(file_stem, ".png")), pdf = file.path(out_dir, paste0(file_stem, ".pdf")))
  ggplot2::ggsave(paths$png, plot, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(paths$pdf, plot, width = width, height = height)
  invisible(list(plot = plot, paths = paths))
}

#' Write standalone bulk RNA-seq QC and differential-expression figures
#'
#' Produce individual library-size, expression-PCA, and volcano plots, plus a
#' DEG table. If `de_results` is omitted, a simple Welch t-test on log2 counts
#' is used as a transparent screening demonstration. For real studies, supply
#' DESeq2, edgeR, limma-voom, or another precomputed differential-expression
#' result instead.
#'
#' @param counts Gene-by-sample count matrix, data frame, or TSV path.
#' @param sample_metadata Data frame or TSV with sample IDs and conditions.
#' @param out_dir Directory for figures and the DEG TSV.
#' @param sample_id Sample-ID column in `sample_metadata`.
#' @param condition Condition column in `sample_metadata`.
#' @param de_results Optional differential-expression table with a gene ID,
#'   log2 fold-change, and p-value column.
#' @param fdr_threshold Adjusted p-value threshold for the DEG list.
#' @param width,height,dpi Figure dimensions and PNG resolution.
#'
#' @return Invisibly returns DE results and output paths.
#' @examples
#' \dontrun{
#' plot_rnaseq_diagnostics(
#'   "results/rna.counts.tsv", "results/rna.sample_metadata.tsv",
#'   out_dir = "results/figures"
#' )
#' }
#' @export
plot_rnaseq_diagnostics <- function(
    counts,
    sample_metadata,
    out_dir,
    sample_id = "sample",
    condition = "condition",
    de_results = NULL,
    fdr_threshold = 0.05,
    width = 7,
    height = 5,
    dpi = 300) {
  .simitall_require_figure_packages("RNA-seq diagnostics")
  if (is.character(counts) && length(counts) == 1L) {
    counts <- .simitall_read_expression_matrix(counts)
  }
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  if (is.null(rownames(counts)) || is.null(colnames(counts)) || nrow(counts) < 2L || ncol(counts) < 4L) {
    stop("counts must be a named gene-by-sample matrix with at least two genes and four samples")
  }
  metadata <- .simitall_figure_data(sample_metadata, "sample_metadata")
  if (!all(c(sample_id, condition) %in% names(metadata))) {
    stop("sample_metadata must contain ", sample_id, " and ", condition, " columns")
  }
  metadata <- metadata[match(colnames(counts), metadata[[sample_id]]), , drop = FALSE]
  if (anyNA(metadata[[sample_id]]) || anyNA(metadata[[condition]])) {
    stop("every count-matrix sample needs a non-missing condition")
  }
  groups <- unique(as.character(metadata[[condition]]))
  if (length(groups) != 2L) stop("RNA-seq diagnostics currently require exactly two conditions")
  log_counts <- log2(pmax(counts, 0) + 1)
  variable <- apply(log_counts, 1L, stats::var, na.rm = TRUE)
  use <- which(is.finite(variable) & variable > 0)
  if (length(use) < 2L) stop("at least two variable genes are required for RNA-seq PCA")
  pca <- stats::prcomp(t(log_counts[use, , drop = FALSE]), scale. = TRUE)
  pca_data <- data.frame(sample = colnames(counts), PC1 = pca$x[, 1L], PC2 = pca$x[, 2L], condition = metadata[[condition]], stringsAsFactors = FALSE)
  if (is.null(de_results)) {
    group_a <- metadata[[condition]] == groups[1L]
    de <- data.frame(
      gene_id = rownames(counts),
      log2_fold_change = rowMeans(log_counts[, !group_a, drop = FALSE]) - rowMeans(log_counts[, group_a, drop = FALSE]),
      p_value = vapply(seq_len(nrow(log_counts)), function(i) stats::t.test(log_counts[i, !group_a], log_counts[i, group_a])$p.value, numeric(1)),
      stringsAsFactors = FALSE
    )
  } else {
    de <- .simitall_figure_data(de_results, "de_results")
    id <- intersect(c("gene_id", "gene", "id", "feature_id"), names(de))[1L]
    fc <- intersect(c("log2_fold_change", "log2FoldChange", "logFC"), names(de))[1L]
    p <- intersect(c("p_value", "pvalue", "PValue"), names(de))[1L]
    if (anyNA(c(id, fc, p))) stop("de_results needs gene ID, log2 fold-change, and p-value columns")
    de <- data.frame(gene_id = de[[id]], log2_fold_change = as.numeric(de[[fc]]), p_value = as.numeric(de[[p]]), stringsAsFactors = FALSE)
  }
  de$q_value <- stats::p.adjust(de$p_value, method = "BH")
  de$significant <- is.finite(de$q_value) & de$q_value <= fdr_threshold
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list(
    library_size = file.path(out_dir, "rnaseq_library_size.png"),
    pca = file.path(out_dir, "rnaseq_pca.png"),
    volcano = file.path(out_dir, "rnaseq_volcano.png"),
    deg = file.path(out_dir, "rnaseq_deg.tsv")
  )
  library_data <- data.frame(sample = colnames(counts), condition = metadata[[condition]], library_size = colSums(counts), stringsAsFactors = FALSE)
  library_plot <- ggplot2::ggplot(library_data, ggplot2::aes(x = condition, y = library_size, fill = condition)) +
    ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.10, alpha = 0.5) +
    ggplot2::scale_y_log10() + ggplot2::labs(title = "RNA-seq library-size QC", x = NULL, y = "Library size") + ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "none")
  pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = condition)) +
    ggplot2::geom_point(size = 2.3, alpha = 0.8) + ggplot2::labs(title = "RNA-seq expression PCA", x = "PC1", y = "PC2", color = NULL) + ggplot2::theme_minimal(base_size = 11)
  volcano_plot <- ggplot2::ggplot(de, ggplot2::aes(x = log2_fold_change, y = -log10(pmax(p_value, .Machine$double.xmin)), color = significant)) +
    ggplot2::geom_point(alpha = 0.75) + ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "#315C6B", `TRUE` = "#D97706")) +
    ggplot2::labs(title = "RNA-seq differential-expression volcano plot", x = "Log2 fold change", y = expression(-log[10](p)), color = "FDR significant") + ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(paths$library_size, library_plot, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(paths$pca, pca_plot, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(paths$volcano, volcano_plot, width = width, height = height, dpi = dpi)
  utils::write.table(de, paths$deg, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(list(de_results = de, plots = list(library_size = library_plot, pca = pca_plot, volcano = volcano_plot), paths = paths))
}
