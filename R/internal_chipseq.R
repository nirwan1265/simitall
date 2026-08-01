.simitall_chip_read_fasta <- function(path) {
  if (!file.exists(path)) stop("Genome FASTA does not exist: ", path)
  lines <- readLines(path, warn = FALSE)
  headers <- which(startsWith(lines, ">"))
  if (!length(headers)) stop("Genome FASTA has no records")
  ends <- c(headers[-1L] - 1L, length(lines))
  ids <- sub("^>([^[:space:]]+).*$", "\\1", lines[headers])
  sequences <- vapply(seq_along(headers), function(i) {
    paste(lines[(headers[i] + 1L):ends[i]], collapse = "")
  }, character(1L))
  sequences <- toupper(gsub("[[:space:]]", "", sequences))
  if (any(!nzchar(sequences))) stop("Genome FASTA contains an empty record")
  stats::setNames(sequences, make.unique(ids))
}

.simitall_chip_gff_features <- function(path, target_features) {
  if (is.null(path)) return(NULL)
  if (!file.exists(path)) stop("Annotation GFF3 does not exist: ", path)
  gff <- utils::read.delim(
    path, header = FALSE, comment.char = "#", sep = "\t", quote = "",
    stringsAsFactors = FALSE, fill = TRUE
  )
  if (ncol(gff) < 9L || !nrow(gff)) stop("GFF3 contains no feature records")
  names(gff)[1:9] <- c(
    "seqname", "source", "type", "start", "end", "score", "strand",
    "phase", "attributes"
  )
  gff <- gff[gff$type %in% target_features, , drop = FALSE]
  if (!nrow(gff)) return(NULL)
  ids <- sub(".*(?:^|;)ID=([^;]+).*", "\\1", gff$attributes, perl = TRUE)
  missing <- ids == gff$attributes | !nzchar(ids)
  ids[missing] <- paste0(gff$type[missing], "_", which(missing))
  data.frame(
    feature_id = make.unique(ids),
    seqname = as.character(gff$seqname),
    start = as.integer(gff$start),
    end = as.integer(gff$end),
    feature_type = as.character(gff$type),
    strand = as.character(gff$strand),
    stringsAsFactors = FALSE
  )
}

.simitall_chip_bed_features <- function(path) {
  if (is.null(path)) return(NULL)
  if (!file.exists(path)) stop("Peak BED does not exist: ", path)
  bed <- utils::read.delim(
    path, header = FALSE, comment.char = "#", sep = "\t", quote = "",
    stringsAsFactors = FALSE, fill = TRUE
  )
  if (ncol(bed) < 3L || !nrow(bed)) stop("Peak BED has no valid intervals")
  data.frame(
    feature_id = if (ncol(bed) >= 4L) as.character(bed[[4L]]) else {
      paste0("bed_peak_", seq_len(nrow(bed)))
    },
    seqname = as.character(bed[[1L]]),
    start = as.integer(bed[[2L]]) + 1L,
    end = as.integer(bed[[3L]]),
    feature_type = "user_peak",
    strand = if (ncol(bed) >= 6L) as.character(bed[[6L]]) else ".",
    stringsAsFactors = FALSE
  )
}

.simitall_chip_preset <- function(assay_type, histone_mark) {
  if (assay_type == "TF") {
    return(list(features = c("promoter", "enhancer"), width = 200L,
                shape = "narrow"))
  }
  if (assay_type == "nucleosome") {
    return(list(features = c("promoter", "gene"), width = 147L,
                shape = "narrow"))
  }
  presets <- list(
    H3K4me3 = list(features = c("promoter", "transcription_start_site"),
                   width = 1000L, shape = "broad"),
    H3K27ac = list(features = c("promoter", "enhancer"),
                   width = 1500L, shape = "broad"),
    H3K27me3 = list(features = c("silencer", "gene"),
                    width = 10000L, shape = "broad"),
    H3K36me3 = list(features = c("gene", "CDS"),
                    width = 5000L, shape = "broad")
  )
  if (!histone_mark %in% names(presets)) {
    stop("Unsupported histone_mark. Use H3K4me3, H3K27ac, H3K27me3, or H3K36me3")
  }
  presets[[histone_mark]]
}

.simitall_chip_make_peaks <- function(
    genome, features, n_peaks, peak_width, peak_width_sd, peak_jitter,
    peak_shape, enrichment_mean, enrichment_sdlog) {
  lengths <- nchar(genome)
  if (is.null(features) || !nrow(features)) {
    seqname <- sample(names(genome), n_peaks, replace = TRUE, prob = lengths)
    centers <- vapply(seq_len(n_peaks), function(i) {
      sample.int(lengths[[seqname[i]]], 1L)
    }, integer(1L))
    features <- data.frame(
      feature_id = paste0("random_target_", seq_len(n_peaks)),
      seqname = seqname, start = centers, end = centers,
      feature_type = "random", strand = ".", stringsAsFactors = FALSE
    )
    selected <- seq_len(n_peaks)
  } else {
    valid <- features$seqname %in% names(genome)
    features <- features[valid, , drop = FALSE]
    if (!nrow(features)) stop("No target features match the FASTA contig names")
    selected <- sample(seq_len(nrow(features)), n_peaks, replace = n_peaks > nrow(features))
  }
  targets <- features[selected, , drop = FALSE]
  centers <- round((targets$start + targets$end) / 2) +
    round(stats::rnorm(n_peaks, 0, peak_jitter))
  widths <- pmax(1L, round(stats::rnorm(n_peaks, peak_width, peak_width_sd)))
  chr_lengths <- lengths[targets$seqname]
  centers <- pmax(1L, pmin(centers, chr_lengths))
  starts <- pmax(1L, centers - floor(widths / 2))
  ends <- pmin(chr_lengths, starts + widths - 1L)
  starts <- pmax(1L, ends - widths + 1L)
  data.frame(
    peak_id = sprintf("peak_%05d", seq_len(n_peaks)),
    seqname = targets$seqname,
    start = as.integer(starts),
    end = as.integer(ends),
    center = as.integer(round((starts + ends) / 2)),
    width = as.integer(ends - starts + 1L),
    feature_id = targets$feature_id,
    feature_type = targets$feature_type,
    peak_shape = peak_shape,
    baseline_enrichment = stats::rlnorm(
      n_peaks, log(enrichment_mean) - 0.5 * enrichment_sdlog^2,
      enrichment_sdlog
    ),
    stringsAsFactors = FALSE
  )
}

.simitall_chip_truncated_fragments <- function(n, mean, sd, min, max) {
  values <- round(stats::rnorm(n, mean, sd))
  invalid <- values < min | values > max
  while (any(invalid)) {
    values[invalid] <- round(stats::rnorm(sum(invalid), mean, sd))
    invalid <- values < min | values > max
  }
  as.integer(values)
}

.simitall_chip_native_positions <- function(
    nreads, genome, peaks, peak_strength, signal_fraction, read_length,
    fragment_mean, fragment_sd, fragment_min, fragment_max) {
  n_signal <- if (nrow(peaks)) stats::rbinom(1L, nreads, signal_fraction) else 0L
  n_background <- nreads - n_signal
  lengths <- nchar(genome)
  result <- vector("list", 2L)
  if (n_signal > 0L) {
    probability <- pmax(peak_strength * peaks$width, 1e-8)
    peak_index <- sample(seq_len(nrow(peaks)), n_signal, TRUE, probability)
    selected <- peaks[peak_index, , drop = FALSE]
    centers <- ifelse(
      selected$peak_shape == "narrow",
      round(stats::rnorm(n_signal, selected$center, pmax(10, selected$width / 5))),
      round(stats::runif(n_signal, selected$start, selected$end))
    )
    centers <- pmax(1L, pmin(centers, lengths[selected$seqname]))
    result[[1L]] <- data.frame(
      seqname = selected$seqname, center = centers,
      source_peak_id = selected$peak_id, stringsAsFactors = FALSE
    )
  }
  if (n_background > 0L) {
    seqname <- sample(names(genome), n_background, TRUE, prob = lengths)
    centers <- vapply(seq_len(n_background), function(i) {
      sample.int(lengths[[seqname[i]]], 1L)
    }, integer(1L))
    result[[2L]] <- data.frame(
      seqname = seqname, center = centers, source_peak_id = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  positions <- do.call(rbind, result[lengths(result) > 0L])
  positions <- positions[sample(seq_len(nrow(positions))), , drop = FALSE]
  positions$fragment_length <- .simitall_chip_truncated_fragments(
    nrow(positions), fragment_mean, fragment_sd, fragment_min, fragment_max
  )
  positions$strand <- sample(c("+", "-"), nrow(positions), TRUE)
  chr_length <- lengths[positions$seqname]
  positions$read_position <- ifelse(
    positions$strand == "+",
    positions$center - floor(positions$fragment_length / 2),
    positions$center + ceiling(positions$fragment_length / 2)
  )
  plus <- positions$strand == "+"
  positions$read_position[plus] <- pmax(
    1L, pmin(positions$read_position[plus], chr_length[plus] - read_length + 1L)
  )
  positions$read_position[!plus] <- pmax(
    read_length, pmin(positions$read_position[!plus], chr_length[!plus])
  )
  positions
}

.simitall_chip_chipsim_positions <- function(
    nreads, genome, peaks, peak_strength, signal_fraction, read_length,
    fragment_mean, fragment_min, fragment_max) {
  if (!requireNamespace("ChIPsim", quietly = TRUE)) {
    stop("backend = 'chipsim' requires BiocManager::install('ChIPsim')")
  }
  lengths <- nchar(genome)
  allocation <- as.vector(stats::rmultinom(1L, nreads, lengths / sum(lengths)))
  rows <- list()
  row_index <- 0L
  fragment <- function(x, minLength, maxLength, meanLength, bind, ...) {
    spread <- max(1, (maxLength - minLength) / 4)
    stats::dnorm(x, meanLength, spread)
  }
  for (i in seq_along(genome)) {
    if (!allocation[i]) next
    seqname <- names(genome)[i]
    chr_length <- lengths[i]
    chr_peaks <- peaks[peaks$seqname == seqname, , drop = FALSE]
    signal <- numeric(chr_length)
    if (nrow(chr_peaks)) {
      strengths <- peak_strength[match(chr_peaks$peak_id, peaks$peak_id)]
      for (j in seq_len(nrow(chr_peaks))) {
        index <- chr_peaks$start[j]:chr_peaks$end[j]
        if (chr_peaks$peak_shape[j] == "narrow") {
          sigma <- max(10, chr_peaks$width[j] / 5)
          signal[index] <- signal[index] + strengths[j] *
            exp(-0.5 * ((index - chr_peaks$center[j]) / sigma)^2)
        } else {
          signal[index] <- signal[index] + strengths[j]
        }
      }
    }
    scale <- if (sum(signal) > 0 && signal_fraction > 0) {
      signal_fraction * chr_length / ((1 - signal_fraction) * sum(signal))
    } else 0
    binding_density <- 1 + scale * signal
    read_density <- ChIPsim::bindDens2readDens(
      binding_density, fragment = fragment, bind = 1L,
      minLength = fragment_min, maxLength = fragment_max,
      meanLength = fragment_mean
    )
    sampled <- ChIPsim::sampleReads(read_density, allocation[i])
    for (strand in c("fwd", "rev")) {
      position <- as.integer(sampled[[strand]])
      if (!length(position)) next
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        seqname = seqname,
        center = position,
        source_peak_id = NA_character_,
        fragment_length = NA_integer_,
        strand = if (strand == "fwd") "+" else "-",
        read_position = position,
        stringsAsFactors = FALSE
      )
    }
  }
  positions <- do.call(rbind, rows)
  chr_length <- lengths[positions$seqname]
  plus <- positions$strand == "+"
  positions$read_position[plus] <- pmax(
    1L, pmin(positions$read_position[plus], chr_length[plus] - read_length + 1L)
  )
  positions$read_position[!plus] <- pmax(
    read_length, pmin(positions$read_position[!plus], chr_length[!plus])
  )
  positions
}

.simitall_chip_overlap <- function(seqname, position, peaks) {
  result <- rep(FALSE, length(position))
  for (chr in intersect(unique(seqname), unique(peaks$seqname))) {
    read_index <- which(seqname == chr)
    intervals <- peaks[peaks$seqname == chr, , drop = FALSE]
    result[read_index] <- vapply(position[read_index], function(x) {
      any(x >= intervals$start & x <= intervals$end)
    }, logical(1L))
  }
  result
}

.simitall_chip_add_duplicates <- function(positions, duplicate_rate) {
  positions$is_duplicate <- FALSE
  n_duplicate <- min(nrow(positions) - 1L, round(nrow(positions) * duplicate_rate))
  if (n_duplicate > 0L) {
    target <- sample(seq_len(nrow(positions)), n_duplicate)
    source <- sample(setdiff(seq_len(nrow(positions)), target), n_duplicate, TRUE)
    columns <- c(
      "seqname", "center", "source_peak_id", "fragment_length", "strand",
      "read_position"
    )
    positions[target, columns] <- positions[source, columns]
    positions$is_duplicate[target] <- TRUE
  }
  positions
}

.simitall_chip_mutate_sequence <- function(sequence, error_rate) {
  if (error_rate <= 0) return(sequence)
  bases <- strsplit(sequence, "", fixed = TRUE)[[1L]]
  mutate <- stats::runif(length(bases)) < error_rate
  if (any(mutate)) {
    alphabet <- c("A", "C", "G", "T")
    bases[mutate] <- vapply(bases[mutate], function(base) {
      sample(setdiff(alphabet, base), 1L)
    }, character(1L))
  }
  paste0(bases, collapse = "")
}

.simitall_chip_write_reads <- function(
    positions, genome, path, read_length, error_rate, library_id, assay) {
  read_ids <- sprintf("%s_%s_read_%08d", library_id, assay, seq_len(nrow(positions)))
  positions$read_id <- read_ids
  sequences <- vapply(seq_len(nrow(positions)), function(i) {
    sequence <- genome[[positions$seqname[i]]]
    if (positions$strand[i] == "+") {
      read <- substr(
        sequence, positions$read_position[i],
        positions$read_position[i] + read_length - 1L
      )
    } else {
      start <- positions$read_position[i] - read_length + 1L
      read <- substr(sequence, start, positions$read_position[i])
      read <- paste(rev(strsplit(chartr("ACGTN", "TGCAN", read), "")[[1L]]),
                    collapse = "")
    }
    .simitall_chip_mutate_sequence(read, error_rate)
  }, character(1L))
  connection <- gzfile(path, "wt")
  on.exit(close(connection), add = TRUE)
  quality <- paste(rep("I", read_length), collapse = "")
  lines <- as.vector(rbind(
    paste0("@", read_ids), sequences, "+", rep(quality, length(read_ids))
  ))
  writeLines(lines, connection)
  positions
}

.simitall_chip_write_bed <- function(positions, path, read_length) {
  start0 <- ifelse(
    positions$strand == "+", positions$read_position - 1L,
    positions$read_position - read_length
  )
  bed <- data.frame(
    seqname = positions$seqname,
    start = pmax(0L, start0),
    end = pmax(0L, start0) + read_length,
    read_id = positions$read_id,
    score = 0L,
    strand = positions$strand,
    library_id = positions$library_id,
    assay = positions$assay,
    is_duplicate = positions$is_duplicate,
    stringsAsFactors = FALSE
  )
  utils::write.table(
    bed, path, sep = "\t", row.names = FALSE, quote = FALSE
  )
  bed
}
