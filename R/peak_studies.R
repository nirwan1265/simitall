# Peak-based assays (ChIP-seq and ATAC-seq): simulation, analysis, and
# scoring against simulation truth. simulate_chipseq() lives in chipseq_api.R;
# this file adds ATAC-seq, a window-based peak caller, a differential test,
# and truth-recovery scores so both assays can be benchmarked end to end.

# A small synthetic genome with a GFF3 of genes, promoters, and enhancers, so
# peak simulations can target regulatory features without external data.
.simitall_simulate_peak_genome <- function(
    out_fa, out_gff3, n_chromosomes = 2L, chromosome_length = 200000L,
    n_genes = 120L, seed = 1L, chromosome_names = NULL) {
  stopifnot(n_chromosomes >= 1, chromosome_length >= 20000, n_genes >= 4)
  set.seed(seed)
  dir.create(dirname(out_fa), recursive = TRUE, showWarnings = FALSE)
  chroms <- if (is.null(chromosome_names)) paste0("chr", seq_len(n_chromosomes)) else chromosome_names
  if (length(chroms) != n_chromosomes || any(!nzchar(chroms)) || anyDuplicated(chroms)) {
    stop("chromosome_names must contain one unique non-empty name per chromosome")
  }
  fasta <- unlist(lapply(chroms, function(chr) {
    bases <- sample(c("A", "C", "G", "T"), chromosome_length, TRUE, prob = c(0.3, 0.2, 0.2, 0.3))
    c(paste0(">", chr), vapply(split(bases, ceiling(seq_along(bases) / 80)), paste, character(1), collapse = ""))
  }))
  writeLines(fasta, out_fa)
  per_chr <- as.vector(table(factor(sample(chroms, n_genes, TRUE), levels = chroms)))
  rows <- list()
  for (k in seq_along(chroms)) {
    if (!per_chr[k]) next
    # Evenly spaced gene slots keep genes, promoters, and enhancers apart.
    slots <- floor(seq(5000, chromosome_length - 8000, length.out = per_chr[k]))
    for (j in seq_along(slots)) {
      gid <- sprintf("gene_%s_%03d", chroms[k], j)
      strand <- sample(c("+", "-"), 1L)
      g_start <- slots[j] + 1000L; g_end <- g_start + sample(1500:3000, 1L)
      p_start <- if (strand == "+") g_start - 1000L else g_end + 1L
      e_start <- slots[j] - 2500L
      rows[[length(rows) + 1L]] <- c(chroms[k], "simitall", "gene", g_start, g_end, ".", strand, ".", paste0("ID=", gid))
      rows[[length(rows) + 1L]] <- c(chroms[k], "simitall", "promoter", p_start, p_start + 999L, ".", strand, ".",
                                     paste0("ID=", gid, "_promoter;Parent=", gid))
      rows[[length(rows) + 1L]] <- c(chroms[k], "simitall", "enhancer", e_start, e_start + 499L, ".", ".", ".",
                                     paste0("ID=", gid, "_enhancer"))
    }
  }
  gff <- do.call(rbind, rows)
  writeLines(c("##gff-version 3", apply(gff, 1L, paste, collapse = "\t")), out_gff3)
  invisible(list(fasta = out_fa, gff3 = out_gff3))
}

#' Simulate an ATAC-seq experiment with known accessible regions
#'
#' Places accessible peaks on regulatory features (promoters and enhancers by
#' default), samples Tn5 fragments from a two-part size distribution
#' (nucleosome-free and mono-nucleosome), applies Tn5 cut-site offsets, and
#' adds condition-dependent differential accessibility with recorded truth.
#' Output files follow simulate_chipseq() so the same analysis applies.
#'
#' @param genome_fa Reference FASTA.
#' @param out_prefix Output prefix.
#' @param annotation_gff3 Optional GFF3 with target features.
#' @param target_features Feature types eligible for peaks.
#' @param n_peaks Number of true accessible regions.
#' @param peak_width Mean peak width (bp).
#' @param conditions Condition names; the first is the reference.
#' @param biological_replicates Replicates per condition.
#' @param n_reads Read pairs per library (simulated at reduced depth for speed).
#' @param signal_fraction Expected fraction of reads in true peaks (FRiP).
#' @param nfr_fraction Fraction of fragments that are nucleosome-free.
#' @param differential_fraction Fraction of peaks with differential accessibility.
#' @param differential_effect_log2 Mean absolute log2 accessibility change.
#' @param replicate_effect_sd Peak-level replicate noise (log2).
#' @param read_length Read length (bp).
#' @param duplicate_rate Fraction of reads replaced by PCR duplicates.
#' @param seed Random seed.
#' @return Invisibly, a named list of output paths.
#' @export
simulate_atacseq <- function(
    genome_fa, out_prefix, annotation_gff3 = NULL,
    target_features = c("promoter", "enhancer"), n_peaks = 100L, peak_width = 300L,
    conditions = c("control", "treatment"), biological_replicates = 3L,
    n_reads = 40000L, signal_fraction = 0.3, nfr_fraction = 0.5,
    differential_fraction = 0.2, differential_effect_log2 = 1.5,
    replicate_effect_sd = 0.15, read_length = 36L, duplicate_rate = 0.03, seed = 1L) {
  for (p in c(signal_fraction, nfr_fraction, differential_fraction, duplicate_rate)) {
    if (!is.numeric(p) || p < 0 || p > 1) stop("Fractions must be between 0 and 1")
  }
  if (n_peaks < 1 || biological_replicates < 1 || n_reads < 100) stop("Peaks, replicates, and reads must be positive")
  if (!length(conditions) || anyDuplicated(conditions)) stop("conditions must be unique")
  set.seed(seed)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  genome <- .simitall_chip_read_fasta(genome_fa)
  features <- .simitall_chip_gff_features(annotation_gff3, target_features)
  peaks <- .simitall_chip_make_peaks(genome, features, as.integer(n_peaks), peak_width, peak_width / 5,
                                     25, "narrow", 8, 0.4)
  peaks$assay_type <- "ATAC"

  n_diff <- min(nrow(peaks), round(nrow(peaks) * differential_fraction))
  truth <- data.frame(peak_id = character(), condition = character(), effect_log2 = numeric(),
                      direction = character(), stringsAsFactors = FALSE)
  if (length(conditions) > 1L && n_diff > 0L) {
    for (condition in conditions[-1L]) {
      selected <- sample(seq_len(nrow(peaks)), n_diff)
      effect <- pmax(0.3, stats::rnorm(n_diff, differential_effect_log2, 0.25)) * sample(c(-1, 1), n_diff, TRUE)
      truth <- rbind(truth, data.frame(peak_id = peaks$peak_id[selected], condition = condition,
                                       effect_log2 = effect, direction = ifelse(effect > 0, "gain", "loss"),
                                       stringsAsFactors = FALSE))
    }
  }
  metadata <- expand.grid(replicate = seq_len(biological_replicates), condition = conditions,
                          stringsAsFactors = FALSE)[, c("condition", "replicate")]
  metadata$library_id <- paste0(metadata$condition, "_rep", metadata$replicate)

  read_rows <- list(); qc_rows <- list()
  for (i in seq_len(nrow(metadata))) {
    effect <- rep(0, nrow(peaks))
    ct <- truth[truth$condition == metadata$condition[i], , drop = FALSE]
    if (nrow(ct)) effect[match(ct$peak_id, peaks$peak_id)] <- ct$effect_log2
    strength <- peaks$baseline_enrichment * 2^(effect + stats::rnorm(nrow(peaks), 0, replicate_effect_sd))
    n_nfr <- stats::rbinom(1L, n_reads, nfr_fraction)
    # Two overlapping fragment classes: nucleosome-free (short) and
    # mono-nucleosome (~200 bp). Di-nucleosome fragments are not simulated.
    parts <- list(
      .simitall_chip_native_positions(n_nfr, genome, peaks, strength, signal_fraction, read_length, 70, 25, max(read_length, 30L), 150L),
      .simitall_chip_native_positions(n_reads - n_nfr, genome, peaks, strength, signal_fraction, read_length, 200, 30, 120L, 320L)
    )
    pos <- do.call(rbind, parts[vapply(parts, function(x) !is.null(x) && nrow(x) > 0L, logical(1))])
    pos <- .simitall_chip_add_duplicates(pos, duplicate_rate)
    # Tn5 inserts with a 9-bp duplication: shift + strand by +4 and - strand by -5.
    len <- nchar(genome)[pos$seqname]
    pos$read_position <- pmin(len, pmax(1L, pos$read_position + ifelse(pos$strand == "+", 4L, -5L)))
    in_peak <- .simitall_chip_overlap(pos$seqname, pos$read_position, peaks)
    read_rows[[i]] <- data.frame(
      read_id = sprintf("%s_read_%07d", metadata$library_id[i], seq_len(nrow(pos))),
      seqname = pos$seqname, read_position = pos$read_position, strand = pos$strand,
      fragment_length = pos$fragment_length, library_id = metadata$library_id[i],
      condition = metadata$condition[i], assay = "atac", is_duplicate = pos$is_duplicate,
      in_peak = in_peak, stringsAsFactors = FALSE)
    qc_rows[[i]] <- data.frame(
      library_id = metadata$library_id[i], condition = metadata$condition[i], total_reads = nrow(pos),
      duplicate_fraction = mean(pos$is_duplicate), frip = mean(in_peak),
      nucleosome_free_fraction = mean(pos$fragment_length < 120), stringsAsFactors = FALSE)
  }
  paths <- list(
    peaks = paste0(out_prefix, ".truth_peaks.tsv"),
    differential_truth = paste0(out_prefix, ".differential_accessibility_truth.tsv"),
    sample_metadata = paste0(out_prefix, ".sample_metadata.tsv"),
    read_positions = paste0(out_prefix, ".read_positions.tsv"),
    qc = paste0(out_prefix, ".qc.tsv")
  )
  w <- function(x, p) utils::write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE)
  w(peaks, paths$peaks); w(truth, paths$differential_truth); w(metadata, paths$sample_metadata)
  w(do.call(rbind, read_rows), paths$read_positions); w(do.call(rbind, qc_rows), paths$qc)
  message("ATAC-seq simulation complete: ", nrow(peaks), " accessible regions, ", nrow(metadata), " libraries.")
  invisible(paths)
}

