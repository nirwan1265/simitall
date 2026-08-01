utils::globalVariables(c(
  "condition", "library_size", "PC1", "PC2", "batch", "variant_mb",
  "neglog10_p", "causal_label", "discovery", "p", "effect_log2", "type",
  "alt_probability", "observed_alt_fraction", "total_count", "metric", "value",
  "cell_type", "n_cells", "gene_label", "observed_cell_type",
  "mean_log_expression", "scope", "beta", "start", "end", "feature_type",
  "distance_bin", "reads_per_million_per_peak", "assay", "frip",
  "unique_fraction", "counts_per_million", "observed_effect_log2",
  "direction", "baseline_enrichment"
))

.simitall_results_table <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " does not exist: ", path)
  }
  result <- utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  result[] <- lapply(
    result,
    utils::type.convert,
    as.is = TRUE,
    tryLogical = FALSE
  )
  result
}

.simitall_empty_panel <- function(title, message) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0.5, y = 0.5, label = message,
      color = "#52616B", size = 4
    ) +
    ggplot2::labs(title = title) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", color = "#17324D")
    )
}

#' Summarize RNA-seq and eQTL simulation results
#'
#' Collect key dimensions, count quality metrics, causal-feature counts, and
#' eQTL recovery statistics from the files written by
#' [simulate_rnaseq_from_gwas()] and [benchmark_eqtl()]. The resulting
#' one-row table is suitable for result reports and simulation manifests.
#'
#' @param rnaseq_prefix Prefix passed to [simulate_rnaseq_from_gwas()].
#' @param eqtl_prefix Optional prefix passed to [benchmark_eqtl()].
#' @param out_tsv Optional output TSV path.
#'
#' @return A one-row data frame, invisibly when `out_tsv` is supplied.
#' @examples
#' \dontrun{
#' summarize_rnaseq_results(
#'   rnaseq_prefix = "results/rnaseq/demo",
#'   eqtl_prefix = "results/eqtl/demo",
#'   out_tsv = "results/figures/rnaseq_summary.tsv"
#' )
#' }
#' @export
summarize_rnaseq_results <- function(
    rnaseq_prefix,
    eqtl_prefix = NULL,
    out_tsv = NULL) {
  counts <- .simitall_read_expression_matrix(
    paste0(rnaseq_prefix, ".counts.tsv")
  )
  metadata <- .simitall_results_table(
    paste0(rnaseq_prefix, ".sample_metadata.tsv"),
    "RNA-seq sample metadata"
  )
  truth <- .simitall_results_table(
    paste0(rnaseq_prefix, ".eqtl_truth.tsv"),
    "eQTL truth table"
  )
  ase <- .simitall_results_table(
    paste0(rnaseq_prefix, ".ase_counts.tsv"),
    "ASE count table"
  )

  observed_library <- if ("observed_library_size" %in% names(metadata)) {
    metadata$observed_library_size
  } else {
    colSums(counts)
  }
  result <- data.frame(
    samples = ncol(counts),
    genes = nrow(counts),
    median_library_size = stats::median(observed_library),
    zero_fraction = mean(counts == 0),
    cis_eqtls = sum(truth$type == "cis"),
    trans_eqtls = sum(truth$type == "trans"),
    gxe_eqtls = sum(truth$gxe_effect_log2 != 0),
    ase_pairs = length(unique(paste(ase$gene_id, ase$variant_id))),
    ase_observations = nrow(ase),
    eqtl_tests = NA_integer_,
    discoveries = NA_integer_,
    precision = NA_real_,
    recall = NA_real_,
    empirical_fdr = NA_real_,
    stringsAsFactors = FALSE
  )

  if (!is.null(eqtl_prefix)) {
    metrics <- .simitall_results_table(
      paste0(eqtl_prefix, ".eqtl_metrics.tsv"),
      "eQTL benchmark metrics"
    )
    if (nrow(metrics)) {
      result$eqtl_tests <- metrics$tests[1L]
      result$discoveries <- metrics$discoveries[1L]
      result$precision <- metrics$precision[1L]
      result$recall <- metrics$recall[1L]
      result$empirical_fdr <- metrics$empirical_fdr[1L]
    }
  }

  if (!is.null(out_tsv)) {
    dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      result, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE
    )
    return(invisible(result))
  }
  result
}

#' Plot RNA-seq and eQTL simulation results
#'
#' Build a six-panel publication figure directly from simulated RNA-seq,
#' causal-truth, ASE, and eQTL benchmark files. Panels show library sizes,
#' expression PCA, association evidence, causal-effect recovery, ASE recovery,
#' and benchmark performance. Source data for every panel are exported beside
#' the PDF and PNG figures.
#'
#' @param rnaseq_prefix Prefix passed to [simulate_rnaseq_from_gwas()].
#' @param eqtl_prefix Prefix passed to [benchmark_eqtl()].
#' @param out_prefix Prefix for PDF, PNG, summary TSV, and source-data outputs.
#' @param top_variable_genes Number of variable genes used for expression PCA.
#' @param fdr_threshold FDR threshold used to label discoveries.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param dpi PNG resolution.
#'
#' @return Invisibly returns a list containing the combined figure, individual
#'   panels, source data, and output paths.
#' @examples
#' \dontrun{
#' plot_rnaseq_eqtl_results(
#'   rnaseq_prefix = "results/rnaseq/demo",
#'   eqtl_prefix = "results/eqtl/demo",
#'   out_prefix = "results/figures/rnaseq_eqtl"
#' )
#' }
#' @export
plot_rnaseq_eqtl_results <- function(
    rnaseq_prefix,
    eqtl_prefix,
    out_prefix,
    top_variable_genes = 500L,
    fdr_threshold = 0.05,
    width = 15,
    height = 10,
    dpi = 300) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "ggplot2 and patchwork are required for result figures. Install them ",
      "with install.packages(c('ggplot2', 'patchwork'))."
    )
  }
  if (top_variable_genes < 2L) {
    stop("top_variable_genes must be at least two")
  }
  if (fdr_threshold <= 0 || fdr_threshold > 1) {
    stop("fdr_threshold must be greater than zero and at most one")
  }

  counts <- .simitall_read_expression_matrix(
    paste0(rnaseq_prefix, ".counts.tsv")
  )
  expression <- .simitall_read_expression_matrix(
    paste0(rnaseq_prefix, ".expression.tsv")
  )
  metadata <- .simitall_results_table(
    paste0(rnaseq_prefix, ".sample_metadata.tsv"),
    "RNA-seq sample metadata"
  )
  truth <- .simitall_results_table(
    paste0(rnaseq_prefix, ".eqtl_truth.tsv"),
    "eQTL truth table"
  )
  ase_truth <- .simitall_results_table(
    paste0(rnaseq_prefix, ".ase_truth.tsv"),
    "ASE truth table"
  )
  ase_counts <- .simitall_results_table(
    paste0(rnaseq_prefix, ".ase_counts.tsv"),
    "ASE count table"
  )
  associations <- .simitall_results_table(
    paste0(eqtl_prefix, ".eqtl_results.tsv"),
    "eQTL association results"
  )
  metrics <- .simitall_results_table(
    paste0(eqtl_prefix, ".eqtl_metrics.tsv"),
    "eQTL benchmark metrics"
  )

  shared_samples <- Reduce(
    intersect,
    list(colnames(counts), colnames(expression), metadata$sample)
  )
  if (length(shared_samples) < 3L) {
    stop("At least three samples must be shared across result files")
  }
  counts <- counts[, shared_samples, drop = FALSE]
  expression <- expression[, shared_samples, drop = FALSE]
  metadata <- metadata[match(shared_samples, metadata$sample), , drop = FALSE]
  metadata$condition <- if ("condition" %in% names(metadata)) {
    as.character(metadata$condition)
  } else {
    "all samples"
  }
  metadata$batch <- if ("batch" %in% names(metadata)) {
    as.character(metadata$batch)
  } else {
    "batch 1"
  }
  metadata$library_size <- if ("observed_library_size" %in% names(metadata)) {
    metadata$observed_library_size
  } else {
    colSums(counts)
  }

  palette <- c(
    navy = "#17324D",
    blue = "#2F6690",
    teal = "#3A7D78",
    sand = "#E7C27D",
    coral = "#D95D4F",
    gray = "#A7B0B7"
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

  panel_a_data <- metadata
  panel_a <- ggplot2::ggplot(
    panel_a_data,
    ggplot2::aes(x = condition, y = library_size, fill = condition)
  ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.75, color = NA) +
    ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white") +
    ggplot2::geom_jitter(width = 0.08, size = 0.8, alpha = 0.35) +
    ggplot2::scale_y_log10(labels = function(value) format(value, scientific = TRUE)) +
    ggplot2::scale_fill_manual(values = rep(
      c(palette[["blue"]], palette[["coral"]], palette[["teal"]], palette[["sand"]]),
      length.out = length(unique(panel_a_data$condition))
    )) +
    ggplot2::labs(
      title = "A. RNA-seq library sizes",
      x = NULL,
      y = "Observed counts (log scale)"
    ) +
    result_theme +
    ggplot2::theme(legend.position = "none")

  gene_variance <- apply(expression, 1L, stats::var, na.rm = TRUE)
  variable_genes <- which(is.finite(gene_variance) & gene_variance > 0)
  if (length(variable_genes) < 2L) {
    stop("Expression matrix needs at least two variable genes for PCA")
  }
  variable_genes <- variable_genes[
    order(gene_variance[variable_genes], decreasing = TRUE)
  ]
  selected_genes <- variable_genes[
    seq_len(min(top_variable_genes, length(variable_genes)))
  ]
  pca <- stats::prcomp(t(expression[selected_genes, , drop = FALSE]), scale. = TRUE)
  explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  panel_b_data <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1L],
    PC2 = pca$x[, 2L],
    condition = metadata$condition[match(rownames(pca$x), metadata$sample)],
    batch = metadata$batch[match(rownames(pca$x), metadata$sample)],
    stringsAsFactors = FALSE
  )
  panel_b <- ggplot2::ggplot(
    panel_b_data,
    ggplot2::aes(x = PC1, y = PC2, color = condition, shape = batch)
  ) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::scale_color_manual(values = rep(
      c(palette[["blue"]], palette[["coral"]], palette[["teal"]], palette[["sand"]]),
      length.out = length(unique(panel_b_data$condition))
    )) +
    ggplot2::labs(
      title = "B. Expression structure",
      subtitle = paste(length(selected_genes), "most-variable genes"),
      x = sprintf("PC1 (%.1f%%)", explained[1L]),
      y = sprintf("PC2 (%.1f%%)", explained[2L]),
      color = "Condition",
      shape = "Batch"
    ) +
    result_theme

  associations$is_causal <- as.logical(associations$is_causal)
  associations$neglog10_p <- -log10(pmax(associations$p_value, .Machine$double.xmin))
  associations$discovery <- !is.na(associations$q_value) &
    associations$q_value <= fdr_threshold
  associations$causal_label <- ifelse(associations$is_causal, "Causal", "Background")
  associations$variant_mb <- associations$variant_pos / 1e6
  panel_c_data <- associations
  panel_c <- ggplot2::ggplot(
    panel_c_data,
    ggplot2::aes(
      x = variant_mb,
      y = neglog10_p,
      color = causal_label,
      shape = discovery
    )
  ) +
    ggplot2::geom_point(alpha = 0.65, size = 1.25, na.rm = TRUE) +
    ggplot2::facet_wrap(~variant_seqname, scales = "free_x") +
    ggplot2::scale_color_manual(values = c(
      Background = palette[["gray"]], Causal = palette[["coral"]]
    )) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
    ggplot2::labs(
      title = "C. eQTL association landscape",
      subtitle = paste0("Triangles: FDR <= ", fdr_threshold),
      x = "Variant position (Mb)",
      y = expression(-log[10](p)),
      color = "Truth",
      shape = "Discovery"
    ) +
    result_theme

  effect_data <- merge(
    truth[, c("gene_id", "variant_id", "type", "effect_log2")],
    associations[, c("gene_id", "variant_id", "beta", "q_value")],
    by = c("gene_id", "variant_id"),
    all.x = TRUE
  )
  effect_complete <- is.finite(effect_data$effect_log2) & is.finite(effect_data$beta)
  effect_correlation <- if (sum(effect_complete) >= 3L) {
    stats::cor(effect_data$effect_log2[effect_complete], effect_data$beta[effect_complete])
  } else {
    NA_real_
  }
  effect_rmse <- if (any(effect_complete)) {
    sqrt(mean((effect_data$effect_log2[effect_complete] - effect_data$beta[effect_complete])^2))
  } else {
    NA_real_
  }
  panel_d_data <- effect_data
  if (any(effect_complete)) {
    panel_d <- ggplot2::ggplot(
      panel_d_data[effect_complete, , drop = FALSE],
      ggplot2::aes(x = effect_log2, y = beta, color = type)
    ) +
      ggplot2::geom_abline(
        slope = 1, intercept = 0, linetype = "dashed", color = "#52616B"
      ) +
      ggplot2::geom_point(size = 2, alpha = 0.8) +
      ggplot2::scale_color_manual(values = c(
        cis = palette[["blue"]], trans = palette[["coral"]]
      )) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        title = "D. Causal effect recovery",
        subtitle = sprintf("Correlation = %.2f; RMSE = %.2f", effect_correlation, effect_rmse),
        x = "True effect (log2)",
        y = "Estimated effect",
        color = "eQTL type"
      ) +
      result_theme
  } else {
    panel_d <- .simitall_empty_panel(
      "D. Causal effect recovery", "No finite causal estimates"
    )
  }

  ase_summary <- data.frame()
  if (nrow(ase_counts) && nrow(ase_truth)) {
    ase_aggregate <- stats::aggregate(
      cbind(alt_count, total_count) ~ gene_id + variant_id,
      data = ase_counts,
      FUN = sum,
      na.rm = TRUE
    )
    ase_aggregate$observed_alt_fraction <- with(
      ase_aggregate,
      ifelse(total_count > 0, alt_count / total_count, NA_real_)
    )
    ase_summary <- merge(
      ase_truth,
      ase_aggregate,
      by = c("gene_id", "variant_id"),
      all.x = TRUE
    )
  }
  ase_complete <- nrow(ase_summary) && any(
    is.finite(ase_summary$alt_probability) &
      is.finite(ase_summary$observed_alt_fraction)
  )
  panel_e_data <- ase_summary
  if (ase_complete) {
    panel_e <- ggplot2::ggplot(
      ase_summary,
      ggplot2::aes(
        x = alt_probability,
        y = observed_alt_fraction,
        size = total_count
      )
    ) +
      ggplot2::geom_abline(
        slope = 1, intercept = 0, linetype = "dashed", color = "#52616B"
      ) +
      ggplot2::geom_point(
        color = palette[["teal"]], alpha = 0.75, na.rm = TRUE
      ) +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(
        title = "E. Allele-specific expression",
        x = "Expected alternate fraction",
        y = "Observed alternate fraction",
        size = "Allelic depth"
      ) +
      result_theme
  } else {
    panel_e <- .simitall_empty_panel(
      "E. Allele-specific expression", "No heterozygous ASE observations"
    )
  }

  performance_data <- data.frame(
    metric = c("Precision", "Recall", "1 - empirical FDR"),
    value = c(
      metrics$precision[1L],
      metrics$recall[1L],
      1 - metrics$empirical_fdr[1L]
    ),
    stringsAsFactors = FALSE
  )
  performance_data$metric <- factor(
    performance_data$metric,
    levels = performance_data$metric
  )
  panel_f_data <- performance_data
  panel_f <- ggplot2::ggplot(
    performance_data,
    ggplot2::aes(x = metric, y = value, fill = metric)
  ) +
    ggplot2::geom_col(width = 0.68, color = "white", na.rm = TRUE) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.finite(value), sprintf("%.2f", value), "NA")),
      vjust = -0.35,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = c(
      palette[["blue"]], palette[["teal"]], palette[["sand"]]
    )) +
    ggplot2::coord_cartesian(ylim = c(0, 1.12), clip = "off") +
    ggplot2::labs(
      title = "F. eQTL benchmark performance",
      subtitle = paste(
        metrics$discoveries[1L], "discoveries from", metrics$tests[1L], "tests"
      ),
      x = NULL,
      y = "Score"
    ) +
    result_theme +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 18, hjust = 1)
    )

  panels <- list(panel_a, panel_b, panel_c, panel_d, panel_e, panel_f)
  figure <- patchwork::wrap_plots(panels, ncol = 3, guides = "collect") +
    patchwork::plot_annotation(
      title = "GWAS-linked RNA-seq and eQTL simulation",
      subtitle = paste(
        ncol(counts), "individuals,", nrow(counts), "genes,",
        nrow(truth), "causal eQTLs"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold", size = 16, color = palette[["navy"]]
        ),
        plot.subtitle = ggplot2::element_text(size = 10, color = "#52616B")
      )
    )

  output_paths <- list(
    pdf = paste0(out_prefix, ".pdf"),
    png = paste0(out_prefix, ".png"),
    summary = paste0(out_prefix, ".summary.tsv"),
    source_data = paste0(out_prefix, "_source_data")
  )
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  dir.create(output_paths$source_data, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    output_paths$pdf, figure, width = width, height = height, units = "in"
  )
  ggplot2::ggsave(
    output_paths$png, figure, width = width, height = height,
    units = "in", dpi = dpi
  )
  summarize_rnaseq_results(rnaseq_prefix, eqtl_prefix, output_paths$summary)

  source_data <- list(
    library_sizes = panel_a_data,
    expression_pca = panel_b_data,
    eqtl_associations = panel_c_data,
    effect_recovery = panel_d_data,
    ase_recovery = panel_e_data,
    benchmark_performance = panel_f_data
  )
  for (name in names(source_data)) {
    utils::write.csv(
      source_data[[name]],
      file.path(output_paths$source_data, paste0(name, ".csv")),
      row.names = FALSE
    )
  }

  message("RNA-seq/eQTL figure written to: ", output_paths$pdf)
  invisible(list(
    figure = figure,
    panels = panels,
    source_data = source_data,
    paths = output_paths
  ))
}
