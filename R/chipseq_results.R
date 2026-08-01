.simitall_chip_nearest_distance <- function(seqname, position, peaks) {
  distance <- rep(NA_real_, length(position))
  for (chr in intersect(unique(seqname), unique(peaks$seqname))) {
    read_index <- which(seqname == chr)
    centers <- sort(peaks$center[peaks$seqname == chr])
    insertion <- findInterval(position[read_index], centers)
    left <- pmax(1L, insertion)
    right <- pmin(length(centers), insertion + 1L)
    left_distance <- position[read_index] - centers[left]
    right_distance <- position[read_index] - centers[right]
    use_right <- abs(right_distance) < abs(left_distance)
    distance[read_index] <- ifelse(use_right, right_distance, left_distance)
  }
  distance
}

#' Plot ChIP-seq simulation results
#'
#' Generate a six-panel publication figure containing genomic peak locations,
#' aggregate peak-centered read enrichment, FRiP, library complexity, peak
#' signal by annotation class, and differential-binding recovery. Panel source
#' tables are exported alongside PDF and PNG files.
#'
#' @param chip_prefix Prefix passed to [simulate_chipseq()].
#' @param out_prefix Prefix for figures, summary, and source-data files.
#' @param profile_window Distance on either side of peak centers in the
#'   aggregate enrichment profile.
#' @param profile_bin Width of peak-centered profile bins.
#' @param max_profile_reads Maximum reads sampled for the profile panel.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param dpi PNG resolution.
#' @param seed Seed used for read subsampling.
#'
#' @return Invisibly returns the combined figure, panels, source data, and
#'   output paths.
#' @examples
#' \dontrun{
#' plot_chipseq_results(
#'   chip_prefix = "results/chipseq/ctcf",
#'   out_prefix = "results/figures/ctcf_chipseq"
#' )
#' }
#' @export
plot_chipseq_results <- function(
    chip_prefix,
    out_prefix,
    profile_window = 1000L,
    profile_bin = 50L,
    max_profile_reads = 200000L,
    width = 15,
    height = 10,
    dpi = 300,
    seed = 1L) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("ggplot2 and patchwork are required for ChIP-seq figures")
  }
  if (profile_window < 1L || profile_bin < 1L || max_profile_reads < 100L) {
    stop("Profile plotting limits must be positive")
  }
  set.seed(seed)
  read_table <- function(suffix, label) {
    path <- paste0(chip_prefix, suffix)
    if (!file.exists(path)) stop(label, " does not exist: ", path)
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  peaks <- read_table(".truth_peaks.tsv", "Peak truth")
  reads <- read_table(".read_positions.tsv", "Read positions")
  metadata <- read_table(".sample_metadata.tsv", "Sample metadata")
  qc <- read_table(".qc.tsv", "QC table")
  peak_counts <- read_table(".peak_counts.tsv", "Peak counts")
  differential <- read_table(
    ".differential_binding_truth.tsv", "Differential-binding truth"
  )
  if (!all(c("start", "end", "seqname", "library_id", "assay") %in%
           names(reads))) {
    stop("Read-position table does not have the expected BED-like columns")
  }
  reads$position <- round((reads$start + reads$end) / 2)
  reads$condition <- metadata$condition[
    match(reads$library_id, metadata$library_id)
  ]

  palette <- c(
    navy = "#153243", blue = "#2D6A8A", cyan = "#4B9DA9",
    green = "#6A994E", sand = "#E5B25D", orange = "#D9783D",
    coral = "#C8553D", gray = "#A8B1B7"
  )
  theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", color = palette[["navy"]], size = 11
      ),
      plot.subtitle = ggplot2::element_text(color = "#52616B", size = 8.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  panel_a_data <- peaks
  panel_a <- ggplot2::ggplot(
    panel_a_data,
    ggplot2::aes(
      x = start / 1e6, xend = end / 1e6,
      y = feature_type, yend = feature_type, color = feature_type
    )
  ) +
    ggplot2::geom_segment(linewidth = 3.5, alpha = 0.75) +
    ggplot2::facet_wrap(~seqname, scales = "free_x") +
    ggplot2::labs(
      title = "A. Programmed binding landscape",
      subtitle = paste(nrow(peaks), "truth peaks"),
      x = "Genomic position (Mb)", y = NULL, color = "Target"
    ) + theme + ggplot2::guides(color = "none")

  profile_reads <- reads
  if (nrow(profile_reads) > max_profile_reads) {
    profile_reads <- profile_reads[
      sort(sample(seq_len(nrow(profile_reads)), max_profile_reads)), ,
      drop = FALSE
    ]
  }
  profile_reads$distance <- .simitall_chip_nearest_distance(
    profile_reads$seqname, profile_reads$position, peaks
  )
  profile_reads <- profile_reads[
    is.finite(profile_reads$distance) &
      abs(profile_reads$distance) <= profile_window, , drop = FALSE
  ]
  breaks <- seq(-profile_window, profile_window + profile_bin, by = profile_bin)
  profile_reads$distance_bin <- breaks[pmax(
    1L, pmin(length(breaks), findInterval(profile_reads$distance, breaks))
  )]
  profile <- stats::aggregate(
    list(reads = rep(1L, nrow(profile_reads))),
    profile_reads[c("distance_bin", "assay")], sum
  )
  assay_totals <- table(reads$assay)
  profile$reads_per_million_per_peak <- profile$reads /
    as.numeric(assay_totals[profile$assay]) * 1e6 / nrow(peaks)
  panel_b_data <- profile
  panel_b <- ggplot2::ggplot(
    panel_b_data,
    ggplot2::aes(
      x = distance_bin, y = reads_per_million_per_peak, color = assay
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = c(chip = palette[["coral"]],
                                           input = palette[["gray"]])) +
    ggplot2::labs(
      title = "B. Peak-centered enrichment",
      subtitle = "Aggregate read profile around truth peak centers",
      x = "Distance from peak center (bp)",
      y = "Reads per million per peak", color = NULL
    ) + theme

  panel_c_data <- qc[qc$assay == "chip", , drop = FALSE]
  panel_c <- ggplot2::ggplot(
    panel_c_data,
    ggplot2::aes(x = condition, y = frip, color = condition)
  ) +
    ggplot2::geom_boxplot(width = 0.45, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.08, size = 1.6) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::labs(
      title = "C. Fraction of reads in peaks",
      subtitle = "ChIP libraries only", x = NULL, y = "FRiP", color = NULL
    ) + theme + ggplot2::guides(color = "none")

  panel_d_data <- qc
  panel_d_data$unique_fraction <- panel_d_data$unique_reads /
    pmax(1, panel_d_data$total_reads)
  panel_d <- ggplot2::ggplot(
    panel_d_data,
    ggplot2::aes(x = assay, y = unique_fraction, color = condition)
  ) +
    ggplot2::geom_boxplot(width = 0.45, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.08, size = 1.4) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::labs(
      title = "D. Library complexity",
      subtitle = "Unique genomic read positions", x = NULL,
      y = "Unique read fraction", color = "Condition"
    ) + theme

  chip_counts <- peak_counts[peak_counts$assay == "chip", , drop = FALSE]
  chip_counts$feature_type <- peaks$feature_type[
    match(chip_counts$peak_id, peaks$peak_id)
  ]
  chip_qc <- qc[qc$assay == "chip", , drop = FALSE]
  library_sizes <- chip_qc$total_reads[
    match(chip_counts$library_id, chip_qc$library_id)
  ]
  chip_counts$counts_per_million <- chip_counts$count /
    pmax(1, library_sizes) * 1e6
  panel_e_data <- chip_counts
  panel_e <- ggplot2::ggplot(
    panel_e_data,
    ggplot2::aes(x = feature_type, y = log2(counts_per_million + 1),
                 fill = feature_type)
  ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.65, color = NA) +
    ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white") +
    ggplot2::labs(
      title = "E. Peak signal by target class",
      subtitle = "Read counts within programmed intervals",
      x = NULL, y = "log2 counts per million + 1", fill = NULL
    ) + theme + ggplot2::guides(fill = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

  if (nrow(differential) && length(unique(chip_counts$condition)) > 1L) {
    reference <- unique(metadata$condition)[1L]
    condition_means <- stats::aggregate(
      counts_per_million ~ peak_id + condition, chip_counts, mean
    )
    observed <- do.call(rbind, lapply(
      unique(differential$condition), function(condition) {
        ref <- condition_means[condition_means$condition == reference, ]
        alt <- condition_means[condition_means$condition == condition, ]
        merged <- merge(ref, alt, by = "peak_id", suffixes = c("_ref", "_alt"))
        merged$condition <- condition
        merged$observed_effect_log2 <- log2(
          (merged$counts_per_million_alt + 0.5) /
            (merged$counts_per_million_ref + 0.5)
        )
        merged
      }
    ))
    panel_f_data <- merge(
      differential, observed[, c("peak_id", "condition", "observed_effect_log2")],
      by = c("peak_id", "condition"), all.x = TRUE
    )
    correlation <- if (nrow(panel_f_data) >= 3L) {
      stats::cor(
        panel_f_data$effect_log2, panel_f_data$observed_effect_log2,
        use = "complete.obs"
      )
    } else NA_real_
    panel_f <- ggplot2::ggplot(
      panel_f_data,
      ggplot2::aes(x = effect_log2, y = observed_effect_log2,
                   color = direction)
    ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = palette[["gray"]]) +
      ggplot2::geom_point(size = 1.7, alpha = 0.75) +
      ggplot2::scale_color_manual(values = c(
        gain = palette[["green"]], loss = palette[["coral"]]
      )) +
      ggplot2::labs(
        title = "F. Differential-binding recovery",
        subtitle = sprintf("Truth versus simulated peak counts; r = %.2f", correlation),
        x = "True log2 effect", y = "Observed log2 fold change", color = NULL
      ) + theme
  } else {
    panel_f_data <- peaks
    panel_f <- ggplot2::ggplot(
      panel_f_data,
      ggplot2::aes(x = baseline_enrichment, fill = feature_type)
    ) +
      ggplot2::geom_histogram(bins = 25, alpha = 0.75) +
      ggplot2::labs(
        title = "F. Peak enrichment truth",
        subtitle = "Add multiple conditions for differential binding",
        x = "Baseline enrichment", y = "Peaks", fill = "Target"
      ) + theme
  }

  combined <- (panel_a | panel_b | panel_c) /
    (panel_d | panel_e | panel_f) +
    patchwork::plot_annotation(
      title = "Annotation-aware ChIP-seq simulation",
      subtitle = paste0(
        unique(peaks$assay_type)[1L], " assay | ", nrow(peaks), " peaks | ",
        length(unique(metadata$biological_id)), " biological replicates | ",
        sum(qc$total_reads), " reads"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold", size = 16, color = palette[["navy"]]
        ),
        plot.subtitle = ggplot2::element_text(size = 10, color = "#52616B")
      )
    )
  source_data <- list(
    panel_a_peak_landscape = panel_a_data,
    panel_b_peak_centered_profile = panel_b_data,
    panel_c_frip = panel_c_data,
    panel_d_library_complexity = panel_d_data,
    panel_e_peak_signal = panel_e_data,
    panel_f_differential_binding = panel_f_data
  )
  source_directory <- paste0(out_prefix, "_source_data")
  dir.create(source_directory, recursive = TRUE, showWarnings = FALSE)
  for (name in names(source_data)) {
    utils::write.csv(
      source_data[[name]], file.path(source_directory, paste0(name, ".csv")),
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
  ggplot2::ggsave(paths$png, combined, width = width, height = height, dpi = dpi)
  summary <- data.frame(
    peaks = nrow(peaks), libraries = nrow(metadata),
    biological_replicates = length(unique(metadata$biological_id)),
    total_reads = sum(qc$total_reads), mean_chip_frip = mean(panel_c_data$frip),
    mean_unique_fraction = mean(panel_d_data$unique_fraction),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    summary, paths$summary, sep = "\t", row.names = FALSE, quote = FALSE
  )
  invisible(list(
    figure = combined,
    panels = list(panel_a, panel_b, panel_c, panel_d, panel_e, panel_f),
    source_data = source_data,
    paths = paths
  ))
}
