.simitall_require_package <- function(package, feature) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Package '", package, "' is required for ", feature, ". ",
      "Install the population dependency profile or run install.packages('",
      package, "').",
      call. = FALSE
    )
  }
}

.simitall_read_analysis_table <- function(data, label) {
  if (is.character(data) && length(data) == 1L) {
    if (!file.exists(data)) stop(label, " file does not exist: ", data)
    data <- utils::read.delim(
      data,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  if (!is.data.frame(data)) {
    stop(label, " must be a data frame or a tab-separated file path")
  }
  as.data.frame(data, check.names = FALSE, stringsAsFactors = FALSE)
}

.simitall_prepare_phenotypes <- function(
    phenotype,
    genotype_samples,
    sample_id,
    trait,
    required_columns = character()) {
  phenotype <- .simitall_read_analysis_table(phenotype, "Phenotype")
  columns <- unique(c(sample_id, trait, required_columns))
  missing_columns <- setdiff(columns, names(phenotype))
  if (length(missing_columns)) {
    stop(
      "Phenotype data are missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  ids <- as.character(phenotype[[sample_id]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Phenotype sample IDs must be present and unique")
  }
  phenotype[[sample_id]] <- ids
  phenotype <- phenotype[match(genotype_samples, ids), , drop = FALSE]
  phenotype[[sample_id]] <- genotype_samples
  suppressWarnings(phenotype[[trait]] <- as.numeric(phenotype[[trait]]))
  phenotype
}

.simitall_write_json <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    value,
    path,
    auto_unbox = TRUE,
    pretty = TRUE,
    na = "null"
  )
  path
}

#' Run a mixed-model GWAS with rrBLUP
#'
#' Test marker-trait associations with the mixed model implemented by
#' [rrBLUP::GWAS()]. Genotypes can be supplied as a VCF or as the marker-by-
#' sample TSV written by [simulate_gwas_cohort()]. Samples are aligned by ID,
#' population structure can be represented by principal components or fixed
#' effects, and a realized additive relationship matrix controls relatedness.
#'
#' @param genotype_file VCF, VCF.GZ, or genotype TSV path.
#' @param phenotype Data frame or phenotype TSV containing one row per sample.
#' @param out_prefix Prefix for `.gwas.tsv` and `.gwas_summary.json` outputs.
#' @param trait Numeric phenotype column to analyze.
#' @param sample_id Sample-ID column in `phenotype`.
#' @param fixed_effects Optional phenotype columns treated as categorical fixed
#'   effects by `rrBLUP`.
#' @param kinship Optional sample-by-sample kinship matrix. When `NULL`, an
#'   additive relationship matrix is estimated with [rrBLUP::A.mat()].
#' @param n_pcs Number of marker-derived principal components included as fixed
#'   effects.
#' @param min_maf Minimum minor-allele frequency passed to [rrBLUP::GWAS()].
#' @param n_cores Number of CPU cores used by `rrBLUP` where supported.
#' @param p3d Estimate variance components once before scanning markers. This is
#'   faster than marker-specific REML but can be less significant.
#' @param fdr_threshold BH-adjusted p-value threshold used to label discoveries.
#'
#' @return Invisibly returns a list containing the GWAS table, kinship matrix,
#'   summary, and output paths.
#' @examples
#' \dontrun{
#' result <- analyze_gwas(
#'   genotype_file = "results/gwas/demo.geno.tsv",
#'   phenotype = "results/gwas/demo.pheno.tsv",
#'   out_prefix = "results/gwas/demo_analysis",
#'   trait = "trait",
#'   fixed_effects = "pop",
#'   n_pcs = 2
#' )
#' }
#' @export
analyze_gwas <- function(
    genotype_file,
    phenotype,
    out_prefix,
    trait = "trait",
    sample_id = "sample",
    fixed_effects = NULL,
    kinship = NULL,
    n_pcs = 0L,
    min_maf = 0.05,
    n_cores = 1L,
    p3d = TRUE,
    fdr_threshold = 0.05) {
  .simitall_require_package("rrBLUP", "GWAS analysis")
  genotype_data <- .simitall_read_genotypes(genotype_file)
  sample_ids <- colnames(genotype_data$genotype)
  phenotype <- .simitall_prepare_phenotypes(
    phenotype,
    sample_ids,
    sample_id,
    trait,
    required_columns = fixed_effects
  )
  observed <- is.finite(phenotype[[trait]])
  if (sum(observed) < 10L) {
    stop("At least ten samples with finite phenotypes are required for GWAS")
  }
  if (!is.numeric(min_maf) || length(min_maf) != 1L ||
      !is.finite(min_maf) || min_maf < 0 || min_maf >= 0.5) {
    stop("min_maf must be one number in [0, 0.5)")
  }
  if (!is.numeric(fdr_threshold) || length(fdr_threshold) != 1L ||
      !is.finite(fdr_threshold) || fdr_threshold <= 0 ||
      fdr_threshold >= 1) {
    stop("fdr_threshold must be one number between zero and one")
  }

  observed_ids <- sample_ids[observed]
  dosage <- genotype_data$genotype[, observed, drop = FALSE]
  rr_dosage <- dosage - 1
  variants <- genotype_data$variants
  geno <- data.frame(
    marker = variants$id,
    chrom = variants$seqname,
    pos = variants$pos,
    rr_dosage,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(geno)[seq.int(4L, ncol(geno))] <- observed_ids

  pheno_columns <- unique(c(sample_id, trait, fixed_effects))
  pheno <- phenotype[observed, pheno_columns, drop = FALSE]
  names(pheno)[names(pheno) == sample_id] <- "line"
  if (is.null(kinship)) {
    kinship <- rrBLUP::A.mat(t(rr_dosage))
  } else {
    kinship <- as.matrix(kinship)
    if (is.null(rownames(kinship)) || is.null(colnames(kinship))) {
      stop("kinship must have sample IDs as row and column names")
    }
    if (any(!observed_ids %in% rownames(kinship)) ||
        any(!observed_ids %in% colnames(kinship))) {
      stop("kinship does not contain every phenotyped sample")
    }
    kinship <- kinship[observed_ids, observed_ids, drop = FALSE]
  }

  scores <- rrBLUP::GWAS(
    pheno = pheno,
    geno = geno,
    fixed = fixed_effects,
    K = kinship,
    n.PC = as.integer(n_pcs),
    min.MAF = min_maf,
    n.core = as.integer(n_cores),
    P3D = isTRUE(p3d),
    plot = FALSE
  )
  score_column <- setdiff(names(scores), c("marker", "chrom", "pos"))
  if (length(score_column) != 1L) {
    stop("rrBLUP returned an unexpected number of trait columns")
  }
  score <- as.numeric(scores[[score_column]])
  p_value <- pmin(1, 10^(-score))
  q_value <- stats::p.adjust(p_value, method = "BH")
  allele_frequency <- rowMeans(dosage, na.rm = TRUE) / 2
  maf <- pmin(allele_frequency, 1 - allele_frequency)
  result <- data.frame(
    marker_id = as.character(scores$marker),
    seqname = as.character(scores$chrom),
    pos = as.numeric(scores$pos),
    score = score,
    p_value = p_value,
    q_value = q_value,
    maf = maf[match(as.character(scores$marker), variants$id)],
    significant = !is.na(q_value) & q_value <= fdr_threshold,
    stringsAsFactors = FALSE
  )
  result <- result[order(result$p_value, na.last = TRUE), , drop = FALSE]
  rownames(result) <- NULL

  result_path <- paste0(out_prefix, ".gwas.tsv")
  summary_path <- paste0(out_prefix, ".gwas_summary.json")
  dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    result,
    result_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  summary <- list(
    backend = "rrBLUP::GWAS",
    trait = trait,
    samples = sum(observed),
    markers_tested = nrow(result),
    fixed_effects = fixed_effects %||% character(),
    principal_components = as.integer(n_pcs),
    minimum_maf = min_maf,
    p3d = isTRUE(p3d),
    fdr_threshold = fdr_threshold,
    discoveries = sum(result$significant, na.rm = TRUE)
  )
  .simitall_write_json(summary, summary_path)
  invisible(list(
    results = result,
    kinship = kinship,
    summary = summary,
    paths = list(results = result_path, summary = summary_path)
  ))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Benchmark GWAS discoveries against simulated causal variants
#'
#' Compare an association table from [analyze_gwas()] with the causal-variant
#' truth table written by [simulate_gwas_cohort()]. Recovery can be evaluated
#' by exact marker identity or by a chromosome-specific distance window.
#'
#' @param gwas_results GWAS result data frame or `.gwas.tsv` path.
#' @param truth Causal-variant data frame or `.causal.tsv` path.
#' @param out_prefix Prefix for annotated results and benchmark metrics.
#' @param fdr_threshold Discovery threshold applied to `q_value`.
#' @param window_bp Non-negative distance around each causal variant. A value of
#'   zero requires an exact marker or coordinate match.
#'
#' @return Invisibly returns annotated associations, metrics, and output paths.
#' @examples
#' \dontrun{
#' benchmark_gwas(
#'   "results/gwas/demo_analysis.gwas.tsv",
#'   "results/gwas/demo.causal.tsv",
#'   "results/gwas/demo_benchmark"
#' )
#' }
#' @export
benchmark_gwas <- function(
    gwas_results,
    truth,
    out_prefix,
    fdr_threshold = 0.05,
    window_bp = 0L) {
  results <- .simitall_read_analysis_table(gwas_results, "GWAS result")
  truth <- .simitall_read_analysis_table(truth, "GWAS truth")
  required_results <- c("marker_id", "seqname", "pos", "q_value")
  required_truth <- c("marker_id", "seqname", "pos")
  if (any(!required_results %in% names(results))) {
    stop("GWAS results need marker_id, seqname, pos, and q_value columns")
  }
  if (any(!required_truth %in% names(truth))) {
    stop("GWAS truth needs marker_id, seqname, and pos columns")
  }
  if (!nrow(truth)) stop("GWAS truth contains no causal variants")
  window_bp <- as.integer(window_bp)
  if (length(window_bp) != 1L || is.na(window_bp) || window_bp < 0L) {
    stop("window_bp must be one non-negative integer")
  }

  exact <- results$marker_id %in% truth$marker_id
  near <- vapply(seq_len(nrow(results)), function(i) {
    same_chromosome <- as.character(truth$seqname) ==
      as.character(results$seqname[i])
    any(same_chromosome &
      abs(as.numeric(truth$pos) - as.numeric(results$pos[i])) <= window_bp)
  }, logical(1L))
  results$is_causal <- exact
  results$within_causal_window <- near
  results$discovery <- !is.na(results$q_value) &
    results$q_value <= fdr_threshold

  recovered_truth <- vapply(seq_len(nrow(truth)), function(i) {
    same_chromosome <- as.character(results$seqname) ==
      as.character(truth$seqname[i])
    any(results$discovery & same_chromosome &
      abs(as.numeric(results$pos) - as.numeric(truth$pos[i])) <= window_bp)
  }, logical(1L))
  true_positive <- sum(results$discovery & near)
  false_positive <- sum(results$discovery & !near)
  false_negative <- sum(!recovered_truth)
  precision <- if (true_positive + false_positive == 0L) {
    NA_real_
  } else {
    true_positive / (true_positive + false_positive)
  }
  recall <- mean(recovered_truth)
  f1 <- if (!is.finite(precision) || precision + recall == 0) {
    NA_real_
  } else {
    2 * precision * recall / (precision + recall)
  }
  metrics <- data.frame(
    fdr_threshold = fdr_threshold,
    window_bp = window_bp,
    causal_variants = nrow(truth),
    discoveries = sum(results$discovery),
    true_positive_markers = true_positive,
    false_positive_markers = false_positive,
    unrecovered_causal_variants = false_negative,
    precision = precision,
    recall = recall,
    f1 = f1,
    stringsAsFactors = FALSE
  )

  result_path <- paste0(out_prefix, ".benchmark.tsv")
  metric_path <- paste0(out_prefix, ".benchmark_metrics.tsv")
  summary_path <- paste0(out_prefix, ".benchmark_summary.json")
  dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    results, result_path, sep = "\t", quote = FALSE, row.names = FALSE
  )
  utils::write.table(
    metrics, metric_path, sep = "\t", quote = FALSE, row.names = FALSE
  )
  .simitall_write_json(as.list(metrics[1L, ]), summary_path)
  invisible(list(
    results = results,
    truth = transform(truth, recovered = recovered_truth),
    metrics = metrics,
    paths = list(
      results = result_path,
      metrics = metric_path,
      summary = summary_path
    )
  ))
}

#' Plot GWAS Manhattan and quantile-quantile panels
#'
#' Create publication-ready Manhattan and QQ plots from [analyze_gwas()] or
#' [benchmark_gwas()] output. Benchmark annotations are used to highlight
#' causal regions when present.
#'
#' @param gwas_results GWAS or benchmark result data frame or TSV path.
#' @param out_prefix Prefix for PDF and PNG figures.
#' @param fdr_threshold Adjusted-p-value threshold shown as a horizontal line.
#' @param width,height Figure dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return Invisibly returns the combined plot, source data, and output paths.
#' @examples
#' \dontrun{
#' plot_gwas_results(
#'   "results/gwas/demo_analysis.gwas.tsv",
#'   "results/figures/demo_gwas"
#' )
#' }
#' @export
plot_gwas_results <- function(
    gwas_results,
    out_prefix,
    fdr_threshold = 0.05,
    width = 13,
    height = 6,
    dpi = 300) {
  .simitall_require_package("ggplot2", "GWAS plotting")
  .simitall_require_package("patchwork", "GWAS plotting")
  results <- .simitall_read_analysis_table(gwas_results, "GWAS result")
  required <- c("marker_id", "seqname", "pos", "p_value", "q_value")
  if (any(!required %in% names(results))) {
    stop("GWAS results need marker_id, seqname, pos, p_value, and q_value")
  }
  results <- results[is.finite(results$pos) & is.finite(results$p_value), ]
  if (!nrow(results)) stop("GWAS results contain no finite association tests")
  chromosome_levels <- unique(as.character(results$seqname))
  results$seqname <- factor(results$seqname, levels = chromosome_levels)
  chromosome_max <- tapply(results$pos, results$seqname, max, na.rm = TRUE)
  offsets <- c(0, head(cumsum(as.numeric(chromosome_max)), -1L))
  names(offsets) <- names(chromosome_max)
  results$cumulative_pos <- results$pos + offsets[as.character(results$seqname)]
  axis_positions <- tapply(
    results$cumulative_pos,
    results$seqname,
    function(value) mean(range(value, na.rm = TRUE))
  )
  results$neglog10_p <- -log10(pmax(results$p_value, .Machine$double.xmin))
  results$highlight <- if ("within_causal_window" %in% names(results)) {
    ifelse(results$within_causal_window, "Causal region", "Background")
  } else {
    "Background"
  }
  threshold_candidates <- results$p_value[
    !is.na(results$q_value) & results$q_value <= fdr_threshold
  ]
  threshold <- if (length(threshold_candidates)) {
    -log10(max(threshold_candidates))
  } else {
    -log10(0.05 / nrow(results))
  }

  manhattan <- ggplot2::ggplot(
    results,
    ggplot2::aes(x = cumulative_pos, y = neglog10_p)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = highlight),
      alpha = 0.75,
      size = 1.2
    ) +
    ggplot2::geom_hline(
      yintercept = threshold,
      linetype = 2,
      color = "#B23A48"
    ) +
    ggplot2::scale_x_continuous(
      breaks = as.numeric(axis_positions),
      labels = names(axis_positions)
    ) +
    ggplot2::scale_color_manual(
      values = c("Background" = "#315C6B", "Causal region" = "#D97706")
    ) +
    ggplot2::labs(
      title = "Genome-wide association landscape",
      x = "Chromosome",
      y = expression(-log[10](p)),
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  observed <- sort(results$p_value)
  expected <- stats::ppoints(length(observed))
  qq_data <- data.frame(
    expected = -log10(expected),
    observed = -log10(pmax(observed, .Machine$double.xmin))
  )
  qq <- ggplot2::ggplot(
    qq_data,
    ggplot2::aes(x = expected, y = observed)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "#B23A48") +
    ggplot2::geom_point(color = "#315C6B", alpha = 0.75, size = 1.3) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Quantile-quantile plot",
      x = expression(Expected~-log[10](p)),
      y = expression(Observed~-log[10](p))
    ) +
    ggplot2::theme_minimal(base_size = 11)

  combined <- manhattan + qq + patchwork::plot_layout(widths = c(2, 1))
  pdf_path <- paste0(out_prefix, ".pdf")
  png_path <- paste0(out_prefix, ".png")
  source_path <- paste0(out_prefix, ".source_data.tsv")
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(pdf_path, combined, width = width, height = height)
  ggplot2::ggsave(
    png_path,
    combined,
    width = width,
    height = height,
    dpi = dpi
  )
  utils::write.table(
    results, source_path, sep = "\t", quote = FALSE, row.names = FALSE
  )
  invisible(list(
    plot = combined,
    manhattan = manhattan,
    qq = qq,
    source_data = list(manhattan = results, qq = qq_data),
    paths = list(pdf = pdf_path, png = png_path, source_data = source_path)
  ))
}
