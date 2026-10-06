# Peak-based assays (ChIP-seq and ATAC-seq): simulation, analysis, and
# scoring against simulation truth. simulate_chipseq() lives in chipseq_api.R;
# this file adds ATAC-seq, a window-based peak caller, a differential test,
# and truth-recovery scores so both assays can be benchmarked end to end.

# A small synthetic genome with a GFF3 of genes, promoters, and enhancers, so
# peak simulations can target regulatory features without external data.
.simitall_simulate_peak_genome <- function(
    out_fa, out_gff3, n_chromosomes = 2L, chromosome_length = 200000L,
    n_genes = 120L, seed = 1L) {
  stopifnot(n_chromosomes >= 1, chromosome_length >= 20000, n_genes >= 4)
  set.seed(seed)
  dir.create(dirname(out_fa), recursive = TRUE, showWarnings = FALSE)
  chroms <- paste0("chr", seq_len(n_chromosomes))
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

.simitall_peak_reads <- function(reads) {
  if (is.character(reads)) reads <- utils::read.delim(reads, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("seqname", "read_position", "library_id", "condition", "assay")
  if (any(!need %in% names(reads))) stop("Read table needs columns: ", paste(need, collapse = ", "))
  reads
}

#' Call peaks with a simple window-based Poisson test
#'
#' Counts reads in fixed windows, tests each window against a background rate
#' (scaled control reads when a control assay is given, otherwise a trimmed
#' genome-wide mean), applies Benjamini-Hochberg correction, and merges
#' adjacent significant windows. This is a transparent baseline caller, not a
#' replacement for MACS.
#'
#' @param reads Read-position table or path (as written by the simulators).
#' @param genome_fa Reference FASTA (for chromosome lengths).
#' @param signal_assay Assay label of signal reads (`"chip"` or `"atac"`).
#' @param control_assay Optional control assay label (e.g. `"input"`).
#' @param window Window width in bp.
#' @param fdr BH-adjusted threshold.
#' @param min_count Minimum reads in a window.
#' @return A data frame of called peaks.
#' @export
call_peaks_windowed <- function(reads, genome_fa, signal_assay = "chip", control_assay = NULL,
                                window = 200L, fdr = 0.05, min_count = 5L) {
  reads <- .simitall_peak_reads(reads)
  if ("is_duplicate" %in% names(reads)) reads <- reads[!as.logical(reads$is_duplicate), , drop = FALSE]
  genome_lengths <- nchar(.simitall_chip_read_fasta(genome_fa))
  sig <- reads[reads$assay == signal_assay, , drop = FALSE]
  ctl <- if (!is.null(control_assay)) reads[reads$assay == control_assay, , drop = FALSE] else NULL
  if (!nrow(sig)) stop("No reads for signal assay '", signal_assay, "'")
  bins <- lapply(names(genome_lengths), function(chr) {
    nb <- ceiling(genome_lengths[[chr]] / window)
    s <- tabulate(floor((sig$read_position[sig$seqname == chr] - 1) / window) + 1L, nb)
    cc <- if (!is.null(ctl)) tabulate(floor((ctl$read_position[ctl$seqname == chr] - 1) / window) + 1L, nb) else NULL
    data.frame(seqname = chr, bin = seq_len(nb), signal = s, control = if (is.null(cc)) NA_real_ else cc)
  })
  bins <- do.call(rbind, bins)
  # Background: trimmed mean of signal windows (peaks excluded by the trim).
  background <- mean(bins$signal[bins$signal <= stats::quantile(bins$signal, 0.9)])
  lambda <- rep(max(background, 0.5), nrow(bins))
  if (!is.null(ctl) && nrow(ctl)) {
    scale <- nrow(sig) / nrow(ctl)
    local <- stats::filter(bins$control, rep(1 / 11, 11), sides = 2)
    local[is.na(local)] <- bins$control[is.na(local)]
    lambda <- pmax(lambda, scale * pmax(bins$control, local))
  }
  p <- stats::ppois(bins$signal - 1, lambda, lower.tail = FALSE)
  bins$q <- stats::p.adjust(p, "BH")
  hit <- bins$q <= fdr & bins$signal >= min_count
  if (!any(hit)) {
    return(data.frame(peak_id = character(), seqname = character(), start = integer(), end = integer(),
                      max_count = integer(), min_q = numeric()))
  }
  b <- bins[hit, , drop = FALSE]
  run <- cumsum(c(TRUE, diff(b$bin) != 1L | b$seqname[-1L] != b$seqname[-nrow(b)]))
  out <- do.call(rbind, lapply(split(b, run), function(x) {
    data.frame(seqname = x$seqname[1L], start = as.integer((min(x$bin) - 1L) * window + 1L),
               end = as.integer(max(x$bin) * window), max_count = max(x$signal), min_q = min(x$q),
               stringsAsFactors = FALSE)
  }))
  out$peak_id <- sprintf("called_%05d", seq_len(nrow(out)))
  rownames(out) <- NULL
  out[, c("peak_id", "seqname", "start", "end", "max_count", "min_q")]
}

.simitall_peak_overlaps <- function(a, b) {
  vapply(seq_len(nrow(a)), function(i) {
    any(b$seqname == a$seqname[i] & b$start <= a$end[i] & b$end >= a$start[i])
  }, logical(1))
}

#' Score called peaks against true peaks
#'
#' @param called Called peaks (seqname, start, end).
#' @param truth True peaks (seqname, start, end), or a path.
#' @return One-row data frame with recall, precision, and F1.
#' @export
score_peak_calls <- function(called, truth) {
  if (is.character(truth)) truth <- utils::read.delim(truth, stringsAsFactors = FALSE)
  recovered <- if (nrow(called)) sum(.simitall_peak_overlaps(truth, called)) else 0L
  true_calls <- if (nrow(called)) sum(.simitall_peak_overlaps(called, truth)) else 0L
  recall <- recovered / nrow(truth)
  precision <- if (nrow(called)) true_calls / nrow(called) else NA_real_
  data.frame(analysis = "peak calling", true_peaks = nrow(truth), called_peaks = nrow(called),
             recovered_peaks = recovered, false_peaks = nrow(called) - true_calls,
             recall = recall, precision = precision,
             f1 = if (is.na(precision) || recall + precision == 0) 0 else 2 * recall * precision / (recall + precision))
}

#' Test peaks for differential binding or accessibility between two conditions
#'
#' Counts reads per peak and library, then tests the second condition against
#' the first with a per-peak t-test on log-CPM, edgeR (exact test), or DESeq2.
#'
#' @param reads Read-position table or path.
#' @param peaks Peaks to test (seqname, start, end, peak_id).
#' @param assay Assay label of the reads to count.
#' @param method `"ttest"`, `"edger"`, or `"deseq2"`.
#' @param fdr BH-adjusted threshold for calling a differential peak.
#' @return Data frame with peak_id, log2_fold_change, p_value, q_value, significant.
#' @export
test_differential_peaks <- function(reads, peaks, assay = "chip", method = c("ttest", "edger", "deseq2"), fdr = 0.05) {
  method <- match.arg(method)
  reads <- .simitall_peak_reads(reads)
  reads <- reads[reads$assay == assay, , drop = FALSE]
  if (!nrow(peaks)) stop("No peaks to test")
  libraries <- unique(reads$library_id)
  condition <- reads$condition[match(libraries, reads$library_id)]
  levels <- unique(condition)
  if (length(levels) < 2L) stop("Differential testing needs two conditions")
  counts <- matrix(0L, nrow(peaks), length(libraries), dimnames = list(peaks$peak_id, libraries))
  for (chr in unique(peaks$seqname)) {
    pk <- which(peaks$seqname == chr)
    pk <- pk[order(peaks$start[pk])]
    r <- reads[reads$seqname == chr, , drop = FALSE]
    idx <- findInterval(r$read_position, peaks$start[pk])
    ok <- idx > 0L
    ok[ok] <- r$read_position[ok] <= peaks$end[pk][idx[ok]]
    if (any(ok)) {
      tab <- table(factor(peaks$peak_id[pk][idx[ok]], levels = peaks$peak_id), factor(r$library_id[ok], levels = libraries))
      counts <- counts + matrix(as.integer(tab), nrow(peaks))
    }
  }
  group <- factor(condition, levels = levels[1:2])
  keep <- !is.na(group)
  counts <- counts[, keep, drop = FALSE]; group <- droplevels(group[keep])
  if (min(table(group)) < 2L) stop("Each condition needs at least two replicates")
  result <- switch(method,
    ttest = {
      lib <- colSums(counts)
      logcpm <- log2(t(t(counts + 0.5) / (lib + 1)) * 1e6)
      a <- group == levels(group)[1L]
      stats <- t(apply(logcpm, 1L, function(x) {
        if (stats::sd(x[a]) + stats::sd(x[!a]) == 0) return(c(mean(x[!a]) - mean(x[a]), 1))
        c(mean(x[!a]) - mean(x[a]), stats::t.test(x[!a], x[a])$p.value)
      }))
      data.frame(log2_fold_change = stats[, 1L], p_value = stats[, 2L])
    },
    edger = {
      if (!requireNamespace("edgeR", quietly = TRUE)) stop("method 'edger' needs the edgeR package")
      y <- edgeR::DGEList(counts, group = group)
      y <- edgeR::calcNormFactors(y)
      y <- edgeR::estimateDisp(y, stats::model.matrix(~group))
      et <- edgeR::exactTest(y)$table
      data.frame(log2_fold_change = et$logFC, p_value = et$PValue)
    },
    deseq2 = {
      if (!requireNamespace("DESeq2", quietly = TRUE)) stop("method 'deseq2' needs the DESeq2 package")
      dds <- DESeq2::DESeqDataSetFromMatrix(counts, data.frame(group = group), ~group)
      dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
      res <- as.data.frame(DESeq2::results(dds))
      data.frame(log2_fold_change = res$log2FoldChange, p_value = ifelse(is.na(res$pvalue), 1, res$pvalue))
    })
  result <- cbind(data.frame(peak_id = peaks$peak_id, stringsAsFactors = FALSE), result)
  result$q_value <- stats::p.adjust(result$p_value, "BH")
  result$significant <- result$q_value <= fdr
  result$method <- method
  result
}

#' Score differential peaks against the simulated truth
#'
#' A called peak is truly differential when it overlaps a true peak with a
#' differential effect. Recall is over true differential peaks; precision is
#' over called peaks declared significant.
#'
#' @param differential Output of test_differential_peaks().
#' @param called Called peaks used in the test.
#' @param truth_peaks True peaks (table or path).
#' @param differential_truth True differential effects (table or path).
#' @return One-row data frame with recall, precision, and F1.
#' @export
score_differential_peaks <- function(differential, called, truth_peaks, differential_truth) {
  if (is.character(truth_peaks)) truth_peaks <- utils::read.delim(truth_peaks, stringsAsFactors = FALSE)
  if (is.character(differential_truth)) differential_truth <- utils::read.delim(differential_truth, stringsAsFactors = FALSE)
  true_diff <- truth_peaks[truth_peaks$peak_id %in% differential_truth$peak_id, , drop = FALSE]
  significant <- called[called$peak_id %in% differential$peak_id[differential$significant], , drop = FALSE]
  recovered <- if (nrow(significant) && nrow(true_diff)) sum(.simitall_peak_overlaps(true_diff, significant)) else 0L
  correct <- if (nrow(significant) && nrow(true_diff)) sum(.simitall_peak_overlaps(significant, true_diff)) else 0L
  recall <- if (nrow(true_diff)) recovered / nrow(true_diff) else NA_real_
  precision <- if (nrow(significant)) correct / nrow(significant) else NA_real_
  data.frame(analysis = "differential peaks", method = differential$method[1L], true_differential = nrow(true_diff),
             called_differential = nrow(significant), recovered = recovered, false_differential = nrow(significant) - correct,
             recall = recall, precision = precision,
             f1 = if (is.na(recall) || is.na(precision) || recall + precision == 0) 0 else 2 * recall * precision / (recall + precision))
}
