#' Simulate an annotation-aware ChIP-seq experiment
#'
#' Generate transcription-factor, histone-mark, or nucleosome ChIP-seq truth,
#' matched input controls, read positions, and optional FASTQ files. Binding
#' sites can be selected from promoters, enhancers, silencers, genes, or any
#' other GFF3 feature, supplied directly as BED intervals, or generated across
#' the reference genome. Technical and biological replication are retained in
#' sample metadata.
#'
#' @param genome_fa Reference genome FASTA.
#' @param out_prefix Prefix for truth, FASTQ, QC, and summary files.
#' @param annotation_gff3 Optional GFF3 containing eligible binding features.
#' @param peak_bed Optional BED intervals used as candidate binding sites.
#' @param assay_type `"TF"`, `"histone"`, or `"nucleosome"`.
#' @param histone_mark Histone preset used when `assay_type = "histone"`:
#'   `"H3K4me3"`, `"H3K27ac"`, `"H3K27me3"`, or `"H3K36me3"`.
#' @param target_features GFF3 feature types eligible for binding. `NULL` uses
#'   the assay preset.
#' @param n_peaks Number of true binding peaks.
#' @param peak_width Mean peak width. `NULL` uses the assay preset.
#' @param peak_width_sd Standard deviation of peak widths.
#' @param peak_jitter Standard deviation of peak-center displacement from the
#'   selected annotation.
#' @param peak_shape `"narrow"` or `"broad"`; `NULL` uses the preset.
#' @param enrichment_mean Mean baseline peak enrichment.
#' @param enrichment_sdlog Log-scale standard deviation of enrichment.
#' @param signal_fraction Expected fraction of ChIP reads sampled from the
#'   peak component before PCR duplication and sequencing error.
#' @param conditions Experimental condition names. The first is the reference.
#' @param biological_replicates Biological replicates per condition.
#' @param technical_replicates Technical libraries per biological replicate.
#' @param differential_binding_fraction Fraction of peaks with condition-
#'   dependent binding in each non-reference condition.
#' @param differential_effect_log2 Mean absolute differential-binding effect.
#' @param differential_effect_sd Standard deviation of absolute effects.
#' @param replicate_effect_sd Peak-level replicate variation in log2 units.
#' @param n_reads ChIP reads generated per library.
#' @param input_reads Matched input-control reads generated per library.
#' @param read_length Single-end read length.
#' @param fragment_mean Mean DNA fragment length.
#' @param fragment_sd Fragment-length standard deviation. Used directly by the
#'   native backend and to describe ChIPsim output.
#' @param fragment_min Minimum fragment length.
#' @param fragment_max Maximum fragment length.
#' @param duplicate_rate Fraction of reads replaced by PCR duplicates.
#' @param sequencing_error_rate Per-base substitution probability.
#' @param backend `"auto"`, `"chipsim"`, or `"native"`. Auto uses ChIPsim
#'   for references no larger than `chipsim_max_genome` and otherwise uses the
#'   scalable native interval sampler.
#' @param chipsim_max_genome Maximum total reference length for the dense
#'   ChIPsim binding-density backend.
#' @param write_fastq Write gzip-compressed ChIP and input FASTQ files.
#' @param write_read_bed Write read-position BED-like tables.
#' @param seed Random-number seed.
#'
#' @return Invisibly returns a named list of output paths. The `backend`
#'   attribute records the backend used.
#' @examples
#' \dontrun{
#' simulate_chipseq(
#'   genome_fa = "genome.fa",
#'   annotation_gff3 = "genome.gff3",
#'   out_prefix = "results/chipseq/ctcf",
#'   assay_type = "TF",
#'   target_features = c("promoter", "enhancer"),
#'   conditions = c("control", "treated"),
#'   biological_replicates = 3,
#'   n_reads = 100000,
#'   seed = 1
#' )
#' }
#' @export
simulate_chipseq <- function(
    genome_fa,
    out_prefix,
    annotation_gff3 = NULL,
    peak_bed = NULL,
    assay_type = c("TF", "histone", "nucleosome"),
    histone_mark = "H3K27ac",
    target_features = NULL,
    n_peaks = 100L,
    peak_width = NULL,
    peak_width_sd = 50,
    peak_jitter = 50,
    peak_shape = NULL,
    enrichment_mean = 8,
    enrichment_sdlog = 0.4,
    signal_fraction = 0.3,
    conditions = "control",
    biological_replicates = 3L,
    technical_replicates = 1L,
    differential_binding_fraction = 0.1,
    differential_effect_log2 = 1,
    differential_effect_sd = 0.25,
    replicate_effect_sd = 0.15,
    n_reads = 100000L,
    input_reads = n_reads,
    read_length = 50L,
    fragment_mean = 180,
    fragment_sd = 30,
    fragment_min = 100L,
    fragment_max = 300L,
    duplicate_rate = 0.05,
    sequencing_error_rate = 0.001,
    backend = c("auto", "chipsim", "native"),
    chipsim_max_genome = 5e6,
    write_fastq = TRUE,
    write_read_bed = TRUE,
    seed = 1L) {
  assay_type <- match.arg(assay_type)
  backend <- match.arg(backend)
  preset <- .simitall_chip_preset(assay_type, histone_mark)
  if (is.null(target_features)) target_features <- preset$features
  if (is.null(peak_width)) peak_width <- preset$width
  if (is.null(peak_shape)) peak_shape <- preset$shape
  peak_shape <- match.arg(peak_shape, c("narrow", "broad"))
  for (item in c(
    signal_fraction = signal_fraction,
    differential_binding_fraction = differential_binding_fraction,
    duplicate_rate = duplicate_rate,
    sequencing_error_rate = sequencing_error_rate
  )) {
    .simitall_scrna_validate_probability(item, deparse(substitute(item)))
  }
  numeric_positive <- c(
    n_peaks = n_peaks, peak_width = peak_width,
    enrichment_mean = enrichment_mean, biological_replicates = biological_replicates,
    technical_replicates = technical_replicates, n_reads = n_reads,
    input_reads = input_reads, read_length = read_length,
    fragment_mean = fragment_mean, fragment_min = fragment_min,
    fragment_max = fragment_max
  )
  if (any(!is.finite(numeric_positive)) || any(numeric_positive < 1)) {
    stop("Peak, replicate, read, and fragment counts must be positive")
  }
  if (fragment_min > fragment_mean || fragment_mean > fragment_max ||
      read_length > fragment_min || peak_width_sd < 0 || peak_jitter < 0 ||
      enrichment_sdlog < 0 || fragment_sd <= 0 ||
      differential_effect_log2 < 0 || differential_effect_sd < 0 ||
      replicate_effect_sd < 0) {
    stop("Peak, fragment, or effect parameters are inconsistent")
  }
  if (!length(conditions) || anyNA(conditions) || any(!nzchar(conditions)) ||
      anyDuplicated(conditions)) {
    stop("conditions must contain unique non-empty names")
  }

  set.seed(seed)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  fastq_dir <- paste0(out_prefix, "_fastq")
  if (isTRUE(write_fastq)) dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)
  genome <- .simitall_chip_read_fasta(genome_fa)
  genome_size <- sum(nchar(genome))
  if (backend == "auto") {
    backend <- if (requireNamespace("ChIPsim", quietly = TRUE) &&
                   genome_size <= chipsim_max_genome) "chipsim" else "native"
  }
  if (backend == "chipsim" && genome_size > chipsim_max_genome) {
    stop("The reference exceeds chipsim_max_genome; use backend = 'native'")
  }

  features <- if (!is.null(peak_bed)) {
    .simitall_chip_bed_features(peak_bed)
  } else {
    .simitall_chip_gff_features(annotation_gff3, target_features)
  }
  peaks <- .simitall_chip_make_peaks(
    genome, features, as.integer(n_peaks), peak_width, peak_width_sd,
    peak_jitter, peak_shape, enrichment_mean, enrichment_sdlog
  )
  peaks$assay_type <- assay_type
  peaks$histone_mark <- if (assay_type == "histone") histone_mark else NA_character_

  n_differential <- min(nrow(peaks), round(
    nrow(peaks) * differential_binding_fraction
  ))
  differential_truth <- data.frame(
    peak_id = character(), condition = character(), effect_log2 = numeric(),
    direction = character(), stringsAsFactors = FALSE
  )
  if (length(conditions) > 1L && n_differential > 0L) {
    for (condition in conditions[-1L]) {
      selected <- sample(seq_len(nrow(peaks)), n_differential)
      magnitude <- pmax(0, stats::rnorm(
        n_differential, differential_effect_log2, differential_effect_sd
      ))
      effect <- magnitude * sample(c(-1, 1), n_differential, TRUE)
      differential_truth <- rbind(
        differential_truth,
        data.frame(
          peak_id = peaks$peak_id[selected], condition = condition,
          effect_log2 = effect,
          direction = ifelse(effect >= 0, "gain", "loss"),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  metadata <- .simitall_expand_experiment_design(
    sample_names = conditions,
    replicates = biological_replicates,
    technical_replicates = technical_replicates,
    replicate_mode = "mixed",
    conditions = conditions,
    batch_levels = "batch1"
  )
  metadata$assay_type <- assay_type
  metadata$histone_mark <- if (assay_type == "histone") histone_mark else NA_character_
  metadata$backend <- backend
  n_reads <- as.integer(.simitall_recycle_parameter(
    n_reads, nrow(metadata), "n_reads"
  ))
  input_reads <- as.integer(.simitall_recycle_parameter(
    input_reads, nrow(metadata), "input_reads"
  ))

  all_positions <- list()
  qc_rows <- list()
  peak_count_rows <- list()
  manifest_rows <- list()
  position_index <- qc_index <- count_index <- manifest_index <- 0L
  for (library_index in seq_len(nrow(metadata))) {
    library_id <- metadata$library_id[library_index]
    condition <- metadata$condition[library_index]
    effect <- rep(0, nrow(peaks))
    condition_truth <- differential_truth[
      differential_truth$condition == condition, , drop = FALSE
    ]
    if (nrow(condition_truth)) {
      effect[match(condition_truth$peak_id, peaks$peak_id)] <-
        condition_truth$effect_log2
    }
    replicate_effect <- stats::rnorm(nrow(peaks), 0, replicate_effect_sd)
    peak_strength <- peaks$baseline_enrichment * 2^(effect + replicate_effect)

    for (assay in c("chip", "input")) {
      read_count <- if (assay == "chip") n_reads[library_index] else {
        input_reads[library_index]
      }
      assay_signal <- if (assay == "chip") signal_fraction else 0
      if (backend == "chipsim") {
        positions <- .simitall_chip_chipsim_positions(
          read_count, genome, peaks, peak_strength, assay_signal, read_length,
          fragment_mean, fragment_min, fragment_max
        )
      } else {
        positions <- .simitall_chip_native_positions(
          read_count, genome, peaks, peak_strength, assay_signal, read_length,
          fragment_mean, fragment_sd, fragment_min, fragment_max
        )
      }
      positions <- .simitall_chip_add_duplicates(positions, duplicate_rate)
      positions$library_id <- library_id
      positions$condition <- condition
      positions$assay <- assay
      positions$read_id <- sprintf(
        "%s_%s_read_%08d", library_id, assay, seq_len(nrow(positions))
      )
      safe_id <- gsub("[^A-Za-z0-9_.-]", "_", library_id)
      fastq_path <- file.path(fastq_dir, paste0(safe_id, ".", assay, ".fastq.gz"))
      if (isTRUE(write_fastq)) {
        positions <- .simitall_chip_write_reads(
          positions, genome, fastq_path, read_length, sequencing_error_rate,
          library_id, assay
        )
      } else {
        fastq_path <- NA_character_
      }
      in_peak <- .simitall_chip_overlap(
        positions$seqname, positions$read_position, peaks
      )
      positions$in_peak <- in_peak
      unique_key <- paste(
        positions$seqname, positions$read_position, positions$strand, sep = ":"
      )
      qc_index <- qc_index + 1L
      qc_rows[[qc_index]] <- data.frame(
        library_id = library_id, condition = condition, assay = assay,
        total_reads = nrow(positions), unique_reads = length(unique(unique_key)),
        duplicate_fraction = mean(duplicated(unique_key)),
        reads_in_peaks = sum(in_peak), frip = mean(in_peak),
        mean_fragment_length = mean(positions$fragment_length, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      for (peak_index in seq_len(nrow(peaks))) {
        count_index <- count_index + 1L
        peak_count_rows[[count_index]] <- data.frame(
          peak_id = peaks$peak_id[peak_index], library_id = library_id,
          condition = condition, assay = assay,
          count = sum(
            positions$seqname == peaks$seqname[peak_index] &
              positions$read_position >= peaks$start[peak_index] &
              positions$read_position <= peaks$end[peak_index]
          ),
          stringsAsFactors = FALSE
        )
      }
      position_index <- position_index + 1L
      all_positions[[position_index]] <- positions
      manifest_index <- manifest_index + 1L
      manifest_rows[[manifest_index]] <- data.frame(
        library_id = library_id, condition = condition, assay = assay,
        fastq = fastq_path, reads = nrow(positions), stringsAsFactors = FALSE
      )
    }
  }
  positions <- do.call(rbind, all_positions)
  qc <- do.call(rbind, qc_rows)
  peak_counts <- do.call(rbind, peak_count_rows)
  manifest <- do.call(rbind, manifest_rows)

  paths <- list(
    peaks = paste0(out_prefix, ".truth_peaks.tsv"),
    peak_bed = paste0(out_prefix, ".truth_peaks.bed"),
    differential_truth = paste0(out_prefix, ".differential_binding_truth.tsv"),
    sample_metadata = paste0(out_prefix, ".sample_metadata.tsv"),
    read_manifest = paste0(out_prefix, ".read_manifest.tsv"),
    read_positions = paste0(out_prefix, ".read_positions.tsv"),
    peak_counts = paste0(out_prefix, ".peak_counts.tsv"),
    qc = paste0(out_prefix, ".qc.tsv"),
    summary = paste0(out_prefix, ".summary.json")
  )
  for (item in list(
    list(peaks, paths$peaks), list(differential_truth, paths$differential_truth),
    list(metadata, paths$sample_metadata), list(manifest, paths$read_manifest),
    list(peak_counts, paths$peak_counts), list(qc, paths$qc)
  )) {
    utils::write.table(
      item[[1L]], item[[2L]], sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
  utils::write.table(
    data.frame(
      seqname = peaks$seqname, start = peaks$start - 1L, end = peaks$end,
      name = peaks$peak_id,
      score = pmin(1000L, round(100 * peaks$baseline_enrichment)),
      strand = ".", stringsAsFactors = FALSE
    ),
    paths$peak_bed, sep = "\t", row.names = FALSE, col.names = FALSE,
    quote = FALSE
  )
  if (isTRUE(write_read_bed)) {
    .simitall_chip_write_bed(positions, paths$read_positions, read_length)
  } else {
    positions_out <- positions[, c(
      "read_id", "seqname", "read_position", "strand", "library_id",
      "condition", "assay", "is_duplicate", "in_peak"
    )]
    utils::write.table(
      positions_out, paths$read_positions, sep = "\t", row.names = FALSE,
      quote = FALSE
    )
  }
  summary <- list(
    seed = seed, backend = backend, genome_fa = normalizePath(genome_fa),
    assay_type = assay_type,
    histone_mark = if (assay_type == "histone") histone_mark else NULL,
    target_features = as.list(target_features),
    counts = list(
      genome_bp = genome_size, peaks = nrow(peaks),
      differential_peaks = length(unique(differential_truth$peak_id)),
      biological_replicates = length(unique(metadata$biological_id)),
      chip_libraries = nrow(metadata), total_reads = nrow(positions)
    ),
    parameters = list(
      signal_fraction = signal_fraction, peak_width = peak_width,
      read_length = read_length, fragment_mean = fragment_mean,
      duplicate_rate = duplicate_rate, sequencing_error_rate = sequencing_error_rate
    )
  )
  jsonlite::write_json(summary, paths$summary, pretty = TRUE, auto_unbox = TRUE)
  attr(paths, "backend") <- backend
  message(
    "ChIP-seq simulation complete: ", nrow(peaks), " peaks, ",
    nrow(metadata), " ChIP libraries plus matched inputs, and ",
    nrow(positions), " total reads (", backend, " backend)."
  )
  invisible(paths)
}
