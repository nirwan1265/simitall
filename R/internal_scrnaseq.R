.simitall_scrna_validate_probability <- function(value, name) {
  if (length(value) != 1L || !is.finite(value) || value < 0 || value > 1) {
    stop(name, " must be a single number between zero and one")
  }
  invisible(value)
}

.simitall_scrna_allocate_types <- function(n, probabilities) {
  n <- as.integer(n)
  allocation <- floor(n * probabilities)
  remainder <- n - sum(allocation)
  if (remainder > 0L) {
    residual <- n * probabilities - allocation
    order_index <- order(residual, runif(length(residual)), decreasing = TRUE)
    allocation[order_index[seq_len(remainder)]] <-
      allocation[order_index[seq_len(remainder)]] + 1L
  }
  if (n >= length(probabilities)) {
    missing <- which(allocation == 0L)
    for (index in missing) {
      donor <- which.max(allocation)
      if (allocation[donor] > 1L) {
        allocation[donor] <- allocation[donor] - 1L
        allocation[index] <- 1L
      }
    }
  }
  allocation
}

.simitall_scrna_make_cells <- function(
    donor_metadata,
    cells_per_donor,
    cell_type_proportions,
    trajectory_cell_types) {
  n_donors <- nrow(donor_metadata)
  if (length(cells_per_donor) == 1L) {
    cells_per_donor <- rep(cells_per_donor, n_donors)
  }
  if (length(cells_per_donor) != n_donors ||
      any(!is.finite(cells_per_donor)) || any(cells_per_donor < 1L)) {
    stop("cells_per_donor must be positive and have length one or n_donors")
  }
  cells_per_donor <- as.integer(round(cells_per_donor))

  rows <- vector("list", n_donors)
  cell_number <- 0L
  for (i in seq_len(n_donors)) {
    allocation <- .simitall_scrna_allocate_types(
      cells_per_donor[i], cell_type_proportions
    )
    cell_types <- rep(names(cell_type_proportions), allocation)
    cell_types <- sample(cell_types, length(cell_types), replace = FALSE)
    ids <- sprintf("cell_%07d", cell_number + seq_along(cell_types))
    cell_number <- cell_number + length(cell_types)
    rows[[i]] <- data.frame(
      cell_id = ids,
      donor = donor_metadata$sample[i],
      cell_type = cell_types,
      stringsAsFactors = FALSE
    )
  }
  cells <- do.call(rbind, rows)
  donor_columns <- setdiff(names(donor_metadata), "sample")
  if (length(donor_columns)) {
    cells <- cbind(
      cells,
      donor_metadata[
        match(cells$donor, donor_metadata$sample), donor_columns, drop = FALSE
      ]
    )
  }
  cells$pseudotime <- NA_real_
  trajectory_index <- cells$cell_type %in% trajectory_cell_types
  cells$pseudotime[trajectory_index] <- stats::runif(sum(trajectory_index))
  cells
}

.simitall_scrna_marker_truth <- function(
    genes,
    cell_types,
    markers_per_cell_type,
    marker_effect_mean,
    marker_effect_sd) {
  markers_per_cell_type <- min(
    as.integer(markers_per_cell_type),
    max(1L, floor(nrow(genes) / length(cell_types)))
  )
  available <- sample(seq_len(nrow(genes)))
  rows <- vector("list", length(cell_types))
  cursor <- 1L
  for (i in seq_along(cell_types)) {
    end <- min(length(available), cursor + markers_per_cell_type - 1L)
    selected <- available[seq.int(cursor, end)]
    cursor <- end + 1L
    effects <- pmax(
      0.25,
      stats::rnorm(length(selected), marker_effect_mean, marker_effect_sd)
    )
    rows[[i]] <- data.frame(
      gene_id = genes$gene_id[selected],
      gene_index = selected,
      cell_type = cell_types[i],
      effect_log2 = effects,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

.simitall_scrna_scaffold <- function(
    backend,
    n_genes,
    n_cells,
    baseline,
    library_sizes,
    seed) {
  requested_backend <- backend
  if (backend == "auto") {
    backend <- if (requireNamespace("splatter", quietly = TRUE) &&
                   requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      "splatter"
    } else {
      "native"
    }
  }

  if (backend == "splatter") {
    if (!requireNamespace("splatter", quietly = TRUE) ||
        !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop(
        "backend = 'splatter' requires Bioconductor packages splatter and ",
        "SummarizedExperiment"
      )
    }
    splat <- tryCatch(
      splatter::splatSimulate(
        nGenes = n_genes,
        batchCells = n_cells,
        seed = as.integer(seed),
        dropout.type = "none",
        verbose = FALSE
      ),
      error = identity
    )
    if (inherits(splat, "error")) {
      if (requested_backend == "splatter") {
        stop("Splatter simulation failed: ", conditionMessage(splat))
      }
      warning(
        "Splatter simulation failed; using the native backend instead: ",
        conditionMessage(splat),
        call. = FALSE
      )
      backend <- "native"
    } else {
      assay_names <- SummarizedExperiment::assayNames(splat)
      assay_name <- if ("CellMeans" %in% assay_names) {
        "CellMeans"
      } else if ("BaseCellMeans" %in% assay_names) {
        "BaseCellMeans"
      } else {
        "counts"
      }
      scaffold <- as.matrix(SummarizedExperiment::assay(splat, assay_name))
      scaffold <- pmax(scaffold, 1e-8)
      scaffold <- sweep(scaffold, 2L, pmax(colSums(scaffold), 1), "/")
      scaffold <- sweep(scaffold, 2L, library_sizes, "*")
      return(list(expected = scaffold, backend = backend))
    }
  }

  gene_abundance <- pmax(2^baseline, 1e-8)
  relative_abundance <- gene_abundance / sum(gene_abundance)
  scaffold <- tcrossprod(relative_abundance, library_sizes)
  list(expected = scaffold, backend = backend)
}

.simitall_scrna_write_sparse <- function(counts, genes, cells, out_prefix) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required for sparse single-cell output")
  }
  matrix_path <- paste0(out_prefix, ".matrix.mtx")
  feature_path <- paste0(out_prefix, ".features.tsv")
  barcode_path <- paste0(out_prefix, ".barcodes.tsv")
  Matrix::writeMM(Matrix::Matrix(counts, sparse = TRUE), matrix_path)
  utils::write.table(
    genes[, c("gene_id", "gene_name"), drop = FALSE],
    feature_path,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  writeLines(cells$cell_id, barcode_path)
  list(matrix = matrix_path, features = feature_path, barcodes = barcode_path)
}

.simitall_read_scrna_object <- function(prefix) {
  path <- if (grepl("\\.scrnaseq\\.rds$", prefix)) {
    prefix
  } else {
    paste0(prefix, ".scrnaseq.rds")
  }
  if (!file.exists(path)) {
    stop("Single-cell result object does not exist: ", path)
  }
  object <- readRDS(path)
  required <- c("counts", "cell_metadata", "gene_metadata")
  if (!all(required %in% names(object))) {
    stop("Single-cell result object is missing required components")
  }
  object
}

.simitall_read_table_or_data <- function(value, label) {
  if (is.character(value) && length(value) == 1L) {
    if (!file.exists(value)) stop(label, " does not exist: ", value)
    value <- utils::read.delim(
      value,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      colClasses = "character"
    )
    value[] <- lapply(
      value,
      utils::type.convert,
      as.is = TRUE,
      tryLogical = FALSE
    )
  }
  as.data.frame(value, stringsAsFactors = FALSE)
}
