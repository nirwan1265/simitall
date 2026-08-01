.simitall_recycle_parameter <- function(value, n, name) {
  if (length(value) == 1L) value <- rep(value, n)
  if (length(value) != n) {
    stop(name, " must have length one or match the number of samples")
  }
  value
}

.simitall_expand_experiment_design <- function(
    sample_names,
    replicates,
    technical_replicates,
    replicate_mode,
    conditions,
    batch_levels,
    design = NULL) {
  if (!is.null(design)) {
    design <- .simitall_read_table_or_data(design, "Experiment design")
    if (!"sample" %in% names(design)) {
      stop("design must contain a sample column")
    }
    if (!"condition" %in% names(design)) design$condition <- "control"
    if (!"biological_replicates" %in% names(design)) {
      design$biological_replicates <- 1L
    }
    if (!"technical_replicates" %in% names(design)) {
      design$technical_replicates <- 1L
    }
    if (anyDuplicated(design$sample)) {
      stop("Each sample must appear once in the aggregated design")
    }
    sample_names <- as.character(design$sample)
    conditions <- as.character(design$condition)
    biological_counts <- as.integer(design$biological_replicates)
    technical_counts <- as.integer(design$technical_replicates)
    fixed_batch <- if ("batch" %in% names(design)) {
      as.character(design$batch)
    } else {
      rep(NA_character_, nrow(design))
    }
    replicate_mode <- "custom"
  } else {
    if (!length(sample_names) || anyNA(sample_names) ||
        any(!nzchar(sample_names)) || anyDuplicated(sample_names)) {
      stop("sample_names must contain unique non-empty labels")
    }
    sample_names <- as.character(sample_names)
    n_samples <- length(sample_names)
    replicates <- as.integer(.simitall_recycle_parameter(
      replicates, n_samples, "replicates"
    ))
    technical_replicates <- as.integer(.simitall_recycle_parameter(
      technical_replicates, n_samples, "technical_replicates"
    ))
    conditions <- as.character(.simitall_recycle_parameter(
      conditions, n_samples, "conditions"
    ))
    if (replicate_mode == "technical") {
      biological_counts <- rep(1L, n_samples)
      technical_counts <- replicates
    } else if (replicate_mode == "biological") {
      biological_counts <- replicates
      technical_counts <- rep(1L, n_samples)
    } else {
      biological_counts <- replicates
      technical_counts <- technical_replicates
    }
    fixed_batch <- rep(NA_character_, n_samples)
  }

  if (any(!is.finite(biological_counts)) || any(biological_counts < 1L) ||
      any(!is.finite(technical_counts)) || any(technical_counts < 1L)) {
    stop("Biological and technical replicate counts must be positive integers")
  }
  if (anyNA(conditions) || any(!nzchar(conditions))) {
    stop("conditions cannot be empty or missing")
  }
  if (!length(batch_levels) || anyNA(batch_levels) || any(!nzchar(batch_levels))) {
    stop("batch_levels must contain non-empty labels")
  }

  rows <- list()
  row_number <- 0L
  for (i in seq_along(sample_names)) {
    safe_sample <- gsub("[^A-Za-z0-9_.-]", "_", sample_names[i])
    for (biological in seq_len(biological_counts[i])) {
      biological_id <- sprintf("%s_bio%02d", safe_sample, biological)
      for (technical in seq_len(technical_counts[i])) {
        row_number <- row_number + 1L
        rows[[row_number]] <- data.frame(
          library_id = sprintf("%s_tech%02d", biological_id, technical),
          sample = sample_names[i],
          biological_id = biological_id,
          biological_replicate = biological,
          technical_replicate = technical,
          condition = conditions[i],
          batch = fixed_batch[i],
          replicate_mode = replicate_mode,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  metadata <- do.call(rbind, rows)
  if (anyNA(metadata$batch)) {
    for (condition in unique(metadata$condition)) {
      index <- which(metadata$condition == condition & is.na(metadata$batch))
      if (length(index)) {
        metadata$batch[index] <- sample(
          rep(batch_levels, length.out = length(index)),
          length(index),
          replace = FALSE
        )
      }
    }
  }
  if (anyDuplicated(metadata$library_id)) {
    metadata$library_id <- make.unique(metadata$library_id, sep = "_")
  }
  rownames(metadata) <- metadata$library_id
  metadata
}

.simitall_standalone_genes <- function(annotation_gff3, n_genes) {
  if (!is.null(annotation_gff3)) {
    if (!file.exists(annotation_gff3)) {
      stop("Annotation GFF3 does not exist: ", annotation_gff3)
    }
    genes <- .simitall_read_genes(annotation_gff3)
    if (!is.null(n_genes) && nrow(genes) > n_genes) {
      genes <- genes[sort(sample(seq_len(nrow(genes)), n_genes)), , drop = FALSE]
    }
  } else {
    if (is.null(n_genes)) n_genes <- 1000L
    n_genes <- as.integer(n_genes)
    if (!is.finite(n_genes) || n_genes < 2L) {
      stop("n_genes must be at least two")
    }
    tss <- seq.int(1000L, by = 1000L, length.out = n_genes)
    genes <- data.frame(
      gene_id = sprintf("gene_%05d", seq_len(n_genes)),
      gene_name = sprintf("gene_%05d", seq_len(n_genes)),
      seqname = "synthetic_chr1",
      start = pmax(1L, tss - 499L),
      end = tss + 500L,
      strand = rep(c("+", "-"), length.out = n_genes),
      tss = tss,
      stringsAsFactors = FALSE
    )
  }
  if (nrow(genes) < 2L) stop("At least two genes are required")
  rownames(genes) <- genes$gene_id
  genes
}

.simitall_make_experiment_cells <- function(
    library_metadata,
    cells_per_replicate,
    cell_type_proportions,
    trajectory_cell_types) {
  n_libraries <- nrow(library_metadata)
  cells_per_replicate <- as.integer(.simitall_recycle_parameter(
    cells_per_replicate, n_libraries, "cells_per_replicate"
  ))
  if (any(!is.finite(cells_per_replicate)) || any(cells_per_replicate < 1L)) {
    stop("cells_per_replicate must contain positive integers")
  }

  rows <- vector("list", n_libraries)
  cell_number <- 0L
  for (i in seq_len(n_libraries)) {
    allocation <- .simitall_scrna_allocate_types(
      cells_per_replicate[i], cell_type_proportions
    )
    cell_types <- rep(names(cell_type_proportions), allocation)
    cell_types <- sample(cell_types, length(cell_types), replace = FALSE)
    cell_ids <- sprintf("cell_%07d", cell_number + seq_along(cell_types))
    cell_number <- cell_number + length(cell_types)
    rows[[i]] <- data.frame(
      cell_id = cell_ids,
      donor = library_metadata$biological_id[i],
      library_id = library_metadata$library_id[i],
      sample = library_metadata$sample[i],
      biological_id = library_metadata$biological_id[i],
      biological_replicate = library_metadata$biological_replicate[i],
      technical_replicate = library_metadata$technical_replicate[i],
      condition = library_metadata$condition[i],
      batch = library_metadata$batch[i],
      replicate_mode = library_metadata$replicate_mode[i],
      cell_type = cell_types,
      stringsAsFactors = FALSE
    )
  }
  cells <- do.call(rbind, rows)
  cells$pseudotime <- NA_real_
  trajectory <- cells$cell_type %in% trajectory_cell_types
  cells$pseudotime[trajectory] <- stats::runif(sum(trajectory))
  rownames(cells) <- cells$cell_id
  cells
}

.simitall_validate_effect_scales <- function(values) {
  if (any(!is.finite(values)) || any(values < 0)) {
    stop("Effect and noise standard deviations must be non-negative")
  }
  invisible(values)
}

