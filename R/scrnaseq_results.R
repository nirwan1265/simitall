#' Plot donor-aware single-cell RNA-seq simulation results
#'
#' Build a six-panel publication figure showing cells per donor, single-cell
#' PCA, marker-program recovery, pseudobulk structure, cell-type eQTL effect
#' recovery, and benchmark precision/recall. Every panel's source data are
#' exported beside PDF and PNG figures.
#'
#' @param scrna_prefix Prefix passed to [simulate_scrnaseq_from_gwas()].
#' @param eqtl_prefix Optional prefix passed to [benchmark_celltype_eqtls()].
#'   When omitted, the two eQTL panels explain that no benchmark was supplied.
#' @param out_prefix Prefix for figures, summary TSV, and source-data CSVs.
#' @param top_variable_genes Number of variable genes used for PCA.
#' @param markers_per_cell_type Number of marker genes shown per cell type.
#' @param max_plot_cells Maximum cells displayed in the single-cell PCA.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param dpi PNG resolution.
#' @param seed Seed used for plot-cell subsampling.
#'
#' @return Invisibly returns the combined figure, panel objects, source data,
#'   and output paths.
#' @examples
#' \dontrun{
#' plot_scrnaseq_results(
#'   scrna_prefix = "results/scrna/demo",
#'   eqtl_prefix = "results/scrna/demo_celltype_eqtl",
#'   out_prefix = "results/figures/scrna_demo"
#' )
#' }
#' @export
plot_scrnaseq_results <- function(
    scrna_prefix,
    eqtl_prefix = NULL,
    out_prefix,
    top_variable_genes = 500L,
    markers_per_cell_type = 4L,
    max_plot_cells = 4000L,
    width = 15,
    height = 10,
    dpi = 300,
    seed = 1L) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "ggplot2 and patchwork are required. Install them with ",
      "install.packages(c('ggplot2', 'patchwork'))."
    )
  }
  if (top_variable_genes < 2L || markers_per_cell_type < 1L ||
      max_plot_cells < 10L) {
    stop("PCA, marker, and cell plotting limits are too small")
  }
  set.seed(seed)
  object <- .simitall_read_scrna_object(scrna_prefix)
  counts <- object$counts
  cells <- object$cell_metadata
  genes <- object$gene_metadata
  marker_truth <- object$marker_truth
  if (!all(colnames(counts) == cells$cell_id)) {
    cells <- cells[match(colnames(counts), cells$cell_id), , drop = FALSE]
  }

  palette <- c(
    navy = "#153243",
    blue = "#2D6A8A",
    cyan = "#4B9DA9",
    green = "#6A994E",
    sand = "#E5B25D",
    orange = "#D9783D",
    coral = "#C8553D",
    gray = "#A8B1B7"
  )
  cell_type_colors <- stats::setNames(
    rep(
      c(
        palette[["blue"]], palette[["coral"]], palette[["green"]],
        palette[["sand"]], palette[["cyan"]], palette[["orange"]]
      ),
      length.out = length(unique(cells$cell_type))
    ),
    unique(cells$cell_type)
  )
  result_theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", color = palette[["navy"]], size = 11
      ),
      plot.subtitle = ggplot2::element_text(color = "#52616B", size = 8.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold")
    )

  cell_table <- as.data.frame(
    table(cells$donor, cells$cell_type),
    stringsAsFactors = FALSE
  )
  names(cell_table) <- c("donor", "cell_type", "n_cells")
  panel_a_data <- cell_table
  panel_a <- ggplot2::ggplot(
    panel_a_data,
    ggplot2::aes(x = cell_type, y = n_cells, fill = cell_type)
  )
  if (length(unique(cells$donor)) >= 2L) {
    panel_a <- panel_a +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.72, color = NA) +
      ggplot2::geom_boxplot(
        width = 0.16, outlier.shape = NA, fill = "white"
      ) +
      ggplot2::geom_jitter(width = 0.08, size = 0.7, alpha = 0.35)
  } else {
    panel_a <- panel_a +
      ggplot2::geom_col(width = 0.7, alpha = 0.8) +
      ggplot2::geom_text(
        ggplot2::aes(label = n_cells), vjust = -0.35, size = 3
      )
  }
  donor_label <- if (identical(object$experiment_type, "standalone")) {
    "biological samples"
  } else {
    "GWAS donors"
  }
  panel_a <- panel_a +
    ggplot2::scale_fill_manual(values = cell_type_colors) +
    ggplot2::labs(
      title = "A. Cells per donor",
      subtitle = paste(length(unique(cells$donor)), donor_label),
      x = NULL,
      y = "Cells"
    ) +
    result_theme +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
    )

  plot_cell_index <- seq_len(ncol(counts))
  if (length(plot_cell_index) > max_plot_cells) {
    plot_cell_index <- sort(sample(plot_cell_index, max_plot_cells))
  }
  gene_mean <- Matrix::rowMeans(counts)
  gene_mean_square <- Matrix::rowMeans(counts^2)
  gene_variance <- pmax(0, gene_mean_square - gene_mean^2)
  variable_gene_index <- which(is.finite(gene_variance) & gene_variance > 0)
  if (length(variable_gene_index) < 2L) {
    stop("At least two variable genes are required for single-cell PCA")
  }
  variable_gene_index <- variable_gene_index[
    order(gene_variance[variable_gene_index], decreasing = TRUE)
  ]
  variable_gene_index <- variable_gene_index[
    seq_len(min(top_variable_genes, length(variable_gene_index)))
  ]
  selected_counts <- as.matrix(
    counts[variable_gene_index, plot_cell_index, drop = FALSE]
  )
  selected_libraries <- pmax(1, Matrix::colSums(counts[, plot_cell_index, drop = FALSE]))
  selected_logcounts <- log2(
    sweep(selected_counts, 2L, selected_libraries, "/") * 1e4 + 1
  )
  pca <- stats::prcomp(t(selected_logcounts), center = TRUE, scale. = FALSE)
  explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  panel_b_data <- data.frame(
    cell_id = rownames(pca$x),
    PC1 = pca$x[, 1L],
    PC2 = pca$x[, 2L],
    cell_type = cells$cell_type[
      match(rownames(pca$x), cells$cell_id)
    ],
    donor = cells$donor[match(rownames(pca$x), cells$cell_id)],
    condition = cells$condition[match(rownames(pca$x), cells$cell_id)],
    stringsAsFactors = FALSE
  )
  panel_b <- ggplot2::ggplot(
    panel_b_data,
    ggplot2::aes(x = PC1, y = PC2, color = cell_type)
  ) +
    ggplot2::geom_point(size = 0.65, alpha = 0.55) +
    ggplot2::scale_color_manual(values = cell_type_colors) +
    ggplot2::labs(
      title = "B. Single-cell expression structure",
      subtitle = paste(nrow(panel_b_data), "cells;", length(variable_gene_index),
                       "variable genes"),
      x = sprintf("PC1 (%.1f%%)", explained[1L]),
      y = sprintf("PC2 (%.1f%%)", explained[2L]),
      color = "Cell type"
    ) +
    result_theme

  marker_truth <- marker_truth[order(
    marker_truth$cell_type,
    -marker_truth$effect_log2
  ), , drop = FALSE]
  selected_marker_rows <- unlist(lapply(
    split(seq_len(nrow(marker_truth)), marker_truth$cell_type),
    function(index) utils::head(index, markers_per_cell_type)
  ))
  marker_subset <- marker_truth[selected_marker_rows, , drop = FALSE]
  marker_index <- match(marker_subset$gene_id, rownames(counts))
  marker_counts <- as.matrix(counts[marker_index, , drop = FALSE])
  marker_logcounts <- log2(
    sweep(marker_counts, 2L, pmax(1, Matrix::colSums(counts)), "/") * 1e4 + 1
  )
  type_levels <- unique(cells$cell_type)
  marker_means <- sapply(type_levels, function(cell_type) {
    rowMeans(marker_logcounts[, cells$cell_type == cell_type, drop = FALSE])
  })
  if (is.null(dim(marker_means))) {
    marker_means <- matrix(marker_means, ncol = 1L)
    colnames(marker_means) <- type_levels
  }
  rownames(marker_means) <- marker_subset$gene_id
  panel_c_data <- as.data.frame(
    as.table(marker_means),
    stringsAsFactors = FALSE,
    responseName = "mean_log_expression"
  )
  names(panel_c_data)[1:2] <- c("gene_id", "observed_cell_type")
  panel_c_data$marker_cell_type <- marker_subset$cell_type[
    match(panel_c_data$gene_id, marker_subset$gene_id)
  ]
  panel_c_data$gene_label <- paste0(
    panel_c_data$gene_id,
    " [",
    panel_c_data$marker_cell_type,
    "]"
  )
  panel_c <- ggplot2::ggplot(
    panel_c_data,
    ggplot2::aes(
      x = observed_cell_type,
      y = gene_label,
      fill = mean_log_expression
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradientn(colors = c("#F5F2E8", palette[["sand"]], palette[["coral"]])) +
    ggplot2::labs(
      title = "C. Cell-type marker recovery",
      subtitle = "Labels show the programmed marker cell type",
      x = "Observed cell type",
      y = NULL,
      fill = "Mean log2\nCP10K"
    ) +
    result_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 6.5)
    )

  pseudobulk_expression <- .simitall_read_expression_matrix(
    paste0(scrna_prefix, ".pseudobulk.expression.tsv")
  )
  pseudobulk_metadata <- .simitall_results_table(
    paste0(scrna_prefix, ".pseudobulk.metadata.tsv"),
    "Pseudobulk metadata"
  )
  pseudobulk_metadata <- pseudobulk_metadata[
    match(colnames(pseudobulk_expression), pseudobulk_metadata$pseudobulk_id),
    , drop = FALSE
  ]
  pseudobulk_variance <- apply(pseudobulk_expression, 1L, stats::var)
  pseudobulk_genes <- order(pseudobulk_variance, decreasing = TRUE)
  pseudobulk_genes <- pseudobulk_genes[
    seq_len(min(top_variable_genes, length(pseudobulk_genes)))
  ]
  pseudobulk_pca <- stats::prcomp(
    t(pseudobulk_expression[pseudobulk_genes, , drop = FALSE]),
    center = TRUE,
    scale. = FALSE
  )
  pseudobulk_explained <- 100 * pseudobulk_pca$sdev^2 /
    sum(pseudobulk_pca$sdev^2)
  panel_d_data <- data.frame(
    pseudobulk_id = rownames(pseudobulk_pca$x),
    PC1 = pseudobulk_pca$x[, 1L],
    PC2 = pseudobulk_pca$x[, 2L],
    donor = pseudobulk_metadata$donor,
    cell_type = pseudobulk_metadata$cell_type,
    condition = pseudobulk_metadata$condition,
    stringsAsFactors = FALSE
  )
  panel_d <- ggplot2::ggplot(
    panel_d_data,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = cell_type,
      shape = condition
    )
  ) +
    ggplot2::geom_point(size = 1.8, alpha = 0.75) +
    ggplot2::scale_color_manual(values = cell_type_colors) +
    ggplot2::labs(
      title = "D. Donor-level pseudobulk structure",
      subtitle = "Donors remain the biological replication unit",
      x = sprintf("PC1 (%.1f%%)", pseudobulk_explained[1L]),
      y = sprintf("PC2 (%.1f%%)", pseudobulk_explained[2L]),
      color = "Cell type",
      shape = "Condition"
    ) +
    result_theme +
    ggplot2::guides(color = "none")

  if (!is.null(eqtl_prefix)) {
    associations <- .simitall_results_table(
      paste0(eqtl_prefix, ".celltype_eqtl_results.tsv"),
      "Cell-type eQTL results"
    )
    metrics <- .simitall_results_table(
      paste0(eqtl_prefix, ".celltype_eqtl_metrics.tsv"),
      "Cell-type eQTL metrics"
    )
    truth <- object$eqtl_truth
    truth_rows <- do.call(rbind, lapply(unique(cells$cell_type), function(cell_type) {
      active <- truth$scope == "shared" |
        (truth$scope == "cell_type_specific" & truth$cell_type == cell_type)
      data.frame(
        truth[active, c("gene_id", "variant_id", "effect_log2", "scope")],
        cell_type = cell_type,
        stringsAsFactors = FALSE
      )
    }))
    effect_data <- merge(
      truth_rows,
      associations[, c("gene_id", "variant_id", "cell_type", "beta", "q_value")],
      by = c("gene_id", "variant_id", "cell_type"),
      all.x = TRUE
    )
    complete <- is.finite(effect_data$effect_log2) & is.finite(effect_data$beta)
    effect_correlation <- if (sum(complete) >= 3L) {
      stats::cor(effect_data$effect_log2[complete], effect_data$beta[complete])
    } else NA_real_
    panel_e_data <- effect_data
    panel_e <- ggplot2::ggplot(
      panel_e_data,
      ggplot2::aes(
        x = effect_log2,
        y = beta,
        color = cell_type,
        shape = scope
      )
    ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = palette[["gray"]]) +
      ggplot2::geom_point(size = 1.5, alpha = 0.7, na.rm = TRUE) +
      ggplot2::scale_color_manual(values = cell_type_colors) +
      ggplot2::labs(
        title = "E. Cell-type eQTL effect recovery",
        subtitle = sprintf("Truth versus estimate; r = %.2f", effect_correlation),
        x = "True effect (log2)",
        y = "Estimated effect",
        color = "Cell type",
        shape = "Scope"
      ) +
      result_theme +
      ggplot2::guides(color = "none")

    metric_names <- c("precision", "recall", "empirical_fdr")
    metric_labels <- c(
      precision = "Precision",
      recall = "Recall",
      empirical_fdr = "Empirical FDR"
    )
    panel_f_data <- do.call(rbind, lapply(metric_names, function(metric) {
      data.frame(
        cell_type = metrics$cell_type,
        metric = metric_labels[[metric]],
        value = metrics[[metric]],
        stringsAsFactors = FALSE
      )
    }))
    panel_f <- ggplot2::ggplot(
      panel_f_data,
      ggplot2::aes(x = cell_type, y = value, fill = metric)
    ) +
      ggplot2::geom_col(position = "dodge", width = 0.72, na.rm = TRUE) +
      ggplot2::scale_fill_manual(values = c(
        Precision = palette[["blue"]],
        Recall = palette[["green"]],
        `Empirical FDR` = palette[["coral"]]
      )) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(
        title = "F. Donor-level eQTL benchmark",
        x = NULL,
        y = "Metric value",
        fill = NULL
      ) +
      result_theme +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
  } else {
    panel_e_data <- data.frame(message = "Run benchmark_celltype_eqtls()")
    panel_f_data <- data.frame(message = "Run benchmark_celltype_eqtls()")
    panel_e <- .simitall_empty_panel(
      "E. Cell-type eQTL effect recovery",
      "Run benchmark_celltype_eqtls() to add this panel"
    )
    panel_f <- .simitall_empty_panel(
      "F. Donor-level eQTL benchmark",
      "Run benchmark_celltype_eqtls() to add this panel"
    )
  }

  figure_title <- if (identical(object$experiment_type, "standalone")) {
    "Standalone single-cell RNA-seq experiment simulation"
  } else {
    "Donor-aware single-cell RNA-seq and cell-type eQTL simulation"
  }
  combined <- (panel_a | panel_b | panel_c) /
    (panel_d | panel_e | panel_f) +
    patchwork::plot_annotation(
      title = figure_title,
      subtitle = paste0(
        length(unique(cells$donor)), " donors | ", nrow(cells), " cells | ",
        nrow(genes), " genes | ", object$backend, " backend"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold", size = 16, color = palette[["navy"]]
        ),
        plot.subtitle = ggplot2::element_text(size = 10, color = "#52616B")
      )
    )
  source_data <- list(
    panel_a_cells_per_donor = panel_a_data,
    panel_b_single_cell_pca = panel_b_data,
    panel_c_marker_expression = panel_c_data,
    panel_d_pseudobulk_pca = panel_d_data,
    panel_e_eqtl_effects = panel_e_data,
    panel_f_eqtl_metrics = panel_f_data
  )
  source_directory <- paste0(out_prefix, "_source_data")
  dir.create(source_directory, recursive = TRUE, showWarnings = FALSE)
  for (name in names(source_data)) {
    utils::write.csv(
      source_data[[name]],
      file.path(source_directory, paste0(name, ".csv")),
      row.names = FALSE
    )
  }
  paths <- list(
    pdf = paste0(out_prefix, ".pdf"),
    png = paste0(out_prefix, ".png"),
    summary = paste0(out_prefix, ".summary.tsv"),
    source_data = source_directory
  )
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paths$pdf, combined, width = width, height = height)
  ggplot2::ggsave(
    paths$png, combined, width = width, height = height, dpi = dpi
  )
  summary <- data.frame(
    donors = length(unique(cells$donor)),
    cells = nrow(cells),
    genes = nrow(genes),
    cell_types = length(unique(cells$cell_type)),
    median_umis = stats::median(cells$observed_library_size),
    median_detected_genes = stats::median(cells$detected_genes),
    doublet_fraction = mean(cells$is_doublet),
    zero_fraction = 1 - Matrix::nnzero(counts) / length(counts),
    backend = object$backend,
    stringsAsFactors = FALSE
  )
  utils::write.table(
    summary, paths$summary, sep = "\t", row.names = FALSE, quote = FALSE
  )
  invisible(list(
    figure = combined,
    panels = list(
      panel_a, panel_b, panel_c, panel_d, panel_e, panel_f
    ),
    source_data = source_data,
    paths = paths
  ))
}
