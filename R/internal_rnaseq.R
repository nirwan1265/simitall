.simitall_read_genotypes <- function(path) {
  if (!file.exists(path)) {
    stop("Genotype file does not exist: ", path)
  }

  is_vcf <- grepl("\\.vcf(?:\\.gz)?$", path, ignore.case = TRUE)
  if (is_vcf) {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
      gzfile(path, open = "rt")
    } else {
      file(path, open = "rt")
    }
    on.exit(close(con), add = TRUE)
    lines <- readLines(con, warn = FALSE)
    header_index <- grep("^#CHROM\\t", lines)
    if (length(header_index) != 1L) {
      stop("VCF must contain exactly one #CHROM header line")
    }

    header <- strsplit(sub("^#", "", lines[header_index]), "\\t")[[1L]]
    sample_ids <- header[-seq_len(9L)]
    body <- if (header_index < length(lines)) {
      lines[seq.int(header_index + 1L, length(lines))]
    } else {
      character()
    }
    body <- body[nzchar(body) & !startsWith(body, "#")]
    if (!length(body)) {
      stop("VCF contains no variant records")
    }

    fields <- strsplit(body, "\\t", fixed = FALSE)
    if (any(lengths(fields) < 9L + length(sample_ids))) {
      stop("VCF record has fewer columns than the #CHROM header")
    }

    variants <- data.frame(
      seqname = vapply(fields, `[[`, character(1L), 1L),
      pos = as.integer(vapply(fields, `[[`, character(1L), 2L)),
      id = vapply(fields, `[[`, character(1L), 3L),
      ref = vapply(fields, `[[`, character(1L), 4L),
      alt = vapply(fields, `[[`, character(1L), 5L),
      stringsAsFactors = FALSE
    )
    missing_id <- is.na(variants$id) | variants$id %in% c("", ".")
    variants$id[missing_id] <- paste0(
      variants$seqname[missing_id], ":", variants$pos[missing_id]
    )

    parse_gt <- function(value) {
      gt <- strsplit(value, ":", fixed = TRUE)[[1L]][1L]
      if (!nzchar(gt) || grepl("\\.", gt)) {
        return(NA_real_)
      }
      alleles <- strsplit(gt, "[/|]")[[1L]]
      alleles <- suppressWarnings(as.numeric(alleles))
      if (anyNA(alleles)) NA_real_ else sum(alleles > 0)
    }

    genotype <- matrix(
      NA_real_,
      nrow = length(fields),
      ncol = length(sample_ids),
      dimnames = list(variants$id, sample_ids)
    )
    for (i in seq_along(fields)) {
      genotype[i, ] <- vapply(
        fields[[i]][seq.int(10L, 9L + length(sample_ids))],
        parse_gt,
        numeric(1L)
      )
    }
  } else {
    table <- utils::read.delim(
      path,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (!nrow(table)) {
      stop("Genotype table contains no variants")
    }

    lower_names <- tolower(names(table))
    find_column <- function(candidates, required = FALSE) {
      hit <- match(candidates, lower_names, nomatch = 0L)
      hit <- hit[hit > 0L]
      if (!length(hit)) {
        if (required) {
          stop("Genotype table is missing column: ", candidates[1L])
        }
        return(NULL)
      }
      names(table)[hit[1L]]
    }

    id_col <- find_column(c("id", "variant_id", "marker"), required = TRUE)
    pos_col <- find_column(c("pos", "position", "position_bp"), required = TRUE)
    seq_col <- find_column(c("seqname", "chrom", "chr", "#chrom"))
    ref_col <- find_column("ref")
    alt_col <- find_column("alt")
    metadata_names <- unique(stats::na.omit(c(
      id_col, pos_col, seq_col, ref_col, alt_col,
      names(table)[lower_names %in% c("type", "len", "length", "qual", "filter", "info")]
    )))
    sample_ids <- setdiff(names(table), metadata_names)
    if (!length(sample_ids)) {
      stop("Genotype table contains no sample columns")
    }

    genotype <- as.matrix(table[, sample_ids, drop = FALSE])
    suppressWarnings(storage.mode(genotype) <- "numeric")
    if (all(is.na(genotype))) {
      stop("Sample columns in the genotype table are not numeric")
    }
    variants <- data.frame(
      seqname = if (is.null(seq_col)) "chr1" else as.character(table[[seq_col]]),
      pos = as.integer(table[[pos_col]]),
      id = as.character(table[[id_col]]),
      ref = if (is.null(ref_col)) "N" else as.character(table[[ref_col]]),
      alt = if (is.null(alt_col)) "N" else as.character(table[[alt_col]]),
      stringsAsFactors = FALSE
    )
    rownames(genotype) <- variants$id
    colnames(genotype) <- sample_ids
  }

  if (anyDuplicated(variants$id)) {
    variants$id <- make.unique(variants$id, sep = "_")
    rownames(genotype) <- variants$id
  }
  if (anyDuplicated(colnames(genotype))) {
    stop("Sample IDs must be unique")
  }
  if (anyNA(variants$pos)) {
    stop("Variant positions must be numeric and non-missing")
  }

  list(variants = variants, genotype = genotype)
}

.simitall_read_expression_matrix <- function(expression) {
  if (is.character(expression) && length(expression) == 1L) {
    if (!file.exists(expression)) {
      stop("Expression file does not exist: ", expression)
    }
    expression <- utils::read.delim(
      expression,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (is.data.frame(expression)) {
    id_candidates <- which(tolower(names(expression)) %in% c(
      "gene_id", "transcript_id", "feature_id", "id"
    ))
    if (length(id_candidates)) {
      ids <- as.character(expression[[id_candidates[1L]]])
      expression <- expression[, -id_candidates[1L], drop = FALSE]
    } else if (!is.null(rownames(expression))) {
      ids <- rownames(expression)
    } else {
      stop("Expression data frame needs a gene_id column or row names")
    }
    expression <- as.matrix(expression)
    rownames(expression) <- ids
  } else {
    expression <- as.matrix(expression)
  }

  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("Expression matrix must have gene and sample names")
  }
  suppressWarnings(storage.mode(expression) <- "numeric")
  if (all(is.na(expression))) {
    stop("Expression matrix is not numeric")
  }
  expression
}

.simitall_read_sample_metadata <- function(sample_metadata, sample_ids) {
  if (is.null(sample_metadata)) {
    return(data.frame(sample = sample_ids, stringsAsFactors = FALSE))
  }
  if (is.character(sample_metadata) && length(sample_metadata) == 1L) {
    if (!file.exists(sample_metadata)) {
      stop("Sample metadata file does not exist: ", sample_metadata)
    }
    sample_metadata <- utils::read.delim(
      sample_metadata,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  sample_metadata <- as.data.frame(sample_metadata, stringsAsFactors = FALSE)
  id_column <- names(sample_metadata)[
    tolower(names(sample_metadata)) %in% c("sample", "sample_id", "individual", "id")
  ][1L]
  if (!length(id_column) || is.na(id_column)) {
    stop("Sample metadata needs a sample or sample_id column")
  }
  names(sample_metadata)[names(sample_metadata) == id_column] <- "sample"
  if (anyDuplicated(sample_metadata$sample)) {
    stop("Sample metadata contains duplicate sample IDs")
  }
  missing_samples <- setdiff(sample_ids, sample_metadata$sample)
  if (length(missing_samples)) {
    stop(
      "Sample metadata is missing genotype samples: ",
      paste(utils::head(missing_samples, 5L), collapse = ", ")
    )
  }
  sample_metadata[match(sample_ids, sample_metadata$sample), , drop = FALSE]
}

.simitall_parse_gff_attributes <- function(attributes, keys) {
  vapply(attributes, function(value) {
    fields <- strsplit(value, ";", fixed = TRUE)[[1L]]
    for (key in keys) {
      hit <- grep(paste0("^", key, "="), fields, value = TRUE)
      if (length(hit)) {
        return(utils::URLdecode(sub(paste0("^", key, "="), "", hit[1L])))
      }
    }
    NA_character_
  }, character(1L))
}

.simitall_read_genes <- function(annotation_gff3) {
  lines <- readLines(annotation_gff3, warn = FALSE)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (!length(lines)) {
    stop("GFF3 contains no feature records")
  }
  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[lengths(fields) >= 9L]
  gff <- data.frame(
    seqname = vapply(fields, `[[`, character(1L), 1L),
    type = vapply(fields, `[[`, character(1L), 3L),
    start = as.integer(vapply(fields, `[[`, character(1L), 4L)),
    end = as.integer(vapply(fields, `[[`, character(1L), 5L)),
    strand = vapply(fields, `[[`, character(1L), 7L),
    attributes = vapply(fields, `[[`, character(1L), 9L),
    stringsAsFactors = FALSE
  )

  genes <- gff[tolower(gff$type) == "gene", , drop = FALSE]
  if (nrow(genes)) {
    gene_id <- .simitall_parse_gff_attributes(
      genes$attributes, c("ID", "gene_id", "Name")
    )
    gene_name <- .simitall_parse_gff_attributes(
      genes$attributes, c("Name", "gene_name", "ID")
    )
  } else {
    cds <- gff[tolower(gff$type) %in% c("cds", "exon", "mrna", "transcript"), , drop = FALSE]
    if (!nrow(cds)) {
      stop("GFF3 has no gene, CDS, exon, mRNA, or transcript features")
    }
    parent <- .simitall_parse_gff_attributes(
      cds$attributes, c("Parent", "gene_id", "ID")
    )
    parent[is.na(parent)] <- paste0("feature_", which(is.na(parent)))
    split_rows <- split(seq_len(nrow(cds)), parent)
    genes <- do.call(rbind, lapply(split_rows, function(index) {
      data.frame(
        seqname = cds$seqname[index[1L]],
        type = "gene",
        start = min(cds$start[index]),
        end = max(cds$end[index]),
        strand = cds$strand[index[1L]],
        attributes = "",
        stringsAsFactors = FALSE
      )
    }))
    gene_id <- names(split_rows)
    gene_name <- gene_id
  }

  gene_id[is.na(gene_id) | !nzchar(gene_id)] <- paste0(
    "gene_", which(is.na(gene_id) | !nzchar(gene_id))
  )
  gene_id <- make.unique(gene_id, sep = "_")
  gene_name[is.na(gene_name) | !nzchar(gene_name)] <- gene_id[
    is.na(gene_name) | !nzchar(gene_name)
  ]
  tss <- ifelse(genes$strand == "-", genes$end, genes$start)

  data.frame(
    gene_id = gene_id,
    gene_name = gene_name,
    seqname = genes$seqname,
    start = genes$start,
    end = genes$end,
    strand = genes$strand,
    tss = as.integer(tss),
    stringsAsFactors = FALSE
  )
}

.simitall_make_synthetic_genes <- function(variants, n_genes) {
  n_genes <- as.integer(n_genes)
  if (is.na(n_genes) || n_genes < 1L) {
    stop("n_genes must be at least 1 when annotation_gff3 is omitted")
  }

  chromosomes <- unique(variants$seqname)
  variant_counts <- table(factor(variants$seqname, levels = chromosomes))
  allocation <- pmax(1L, round(n_genes * as.numeric(variant_counts) / sum(variant_counts)))
  while (sum(allocation) > n_genes) {
    eligible <- which(allocation > 1L)
    if (!length(eligible)) break
    allocation[eligible[which.max(allocation[eligible])]] <-
      allocation[eligible[which.max(allocation[eligible])]] - 1L
  }
  while (sum(allocation) < n_genes) {
    allocation[which.max(variant_counts)] <- allocation[which.max(variant_counts)] + 1L
  }

  rows <- vector("list", length(chromosomes))
  gene_number <- 0L
  for (i in seq_along(chromosomes)) {
    chromosome <- chromosomes[i]
    positions <- variants$pos[variants$seqname == chromosome]
    centers <- round(seq(
      min(positions), max(positions),
      length.out = allocation[i] + 2L
    ))[seq_len(allocation[i]) + 1L]
    ids <- sprintf("gene_%05d", gene_number + seq_len(allocation[i]))
    gene_number <- gene_number + allocation[i]
    rows[[i]] <- data.frame(
      gene_id = ids,
      gene_name = ids,
      seqname = chromosome,
      start = pmax(1L, centers - 499L),
      end = centers + 500L,
      strand = rep(c("+", "-"), length.out = allocation[i]),
      tss = centers,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

.simitall_write_feature_matrix <- function(matrix, path, id_name = "gene_id") {
  table <- data.frame(
    id = rownames(matrix),
    matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(table)[1L] <- id_name
  utils::write.table(
    table,
    path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

.simitall_impute_genotypes <- function(genotype) {
  for (i in seq_len(nrow(genotype))) {
    missing <- is.na(genotype[i, ])
    if (any(missing)) {
      replacement <- mean(genotype[i, ], na.rm = TRUE)
      if (!is.finite(replacement)) replacement <- 0
      genotype[i, missing] <- replacement
    }
  }
  genotype
}

.simitall_select_eqtls <- function(
    genes,
    variants,
    genotype,
    n_cis_eqtl,
    n_trans_eqtl,
    cis_window_bp,
    cis_effect_sd,
    trans_effect_sd,
    gxe_fraction,
    gxe_effect_sd,
    treatment_label) {
  variable_variant <- apply(genotype, 1L, function(value) {
    stats::var(value, na.rm = TRUE) > 0
  })
  usable_variants <- which(variable_variant & !is.na(variable_variant))
  if (!length(usable_variants)) {
    stop("No polymorphic variants are available for eQTL simulation")
  }

  cis_candidates <- lapply(seq_len(nrow(genes)), function(gene_index) {
    usable_variants[
      variants$seqname[usable_variants] == genes$seqname[gene_index] &
        abs(variants$pos[usable_variants] - genes$tss[gene_index]) <= cis_window_bp
    ]
  })
  eligible_cis_genes <- which(lengths(cis_candidates) > 0L)
  requested_cis <- min(as.integer(n_cis_eqtl), nrow(genes))
  selected_cis_genes <- if (requested_cis > 0L && length(eligible_cis_genes)) {
    sample(eligible_cis_genes, min(requested_cis, length(eligible_cis_genes)))
  } else {
    integer()
  }

  cis <- if (length(selected_cis_genes)) {
    selected_variants <- vapply(
      selected_cis_genes,
      function(index) sample(cis_candidates[[index]], 1L),
      integer(1L)
    )
    data.frame(
      type = "cis",
      gene_index = selected_cis_genes,
      variant_index = selected_variants,
      effect_log2 = stats::rnorm(length(selected_cis_genes), 0, cis_effect_sd),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      type = character(), gene_index = integer(), variant_index = integer(),
      effect_log2 = numeric(), stringsAsFactors = FALSE
    )
  }

  requested_trans <- min(as.integer(n_trans_eqtl), nrow(genes) * length(usable_variants))
  trans <- data.frame(
    type = character(), gene_index = integer(), variant_index = integer(),
    effect_log2 = numeric(), stringsAsFactors = FALSE
  )
  if (requested_trans > 0L) {
    candidates <- expand.grid(
      gene_index = seq_len(nrow(genes)),
      variant_index = usable_variants,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    is_trans <- variants$seqname[candidates$variant_index] !=
      genes$seqname[candidates$gene_index] |
      abs(variants$pos[candidates$variant_index] - genes$tss[candidates$gene_index]) >
        cis_window_bp
    candidates <- candidates[is_trans, , drop = FALSE]
    if (nrow(candidates)) {
      candidates <- candidates[
        sample(seq_len(nrow(candidates)), min(requested_trans, nrow(candidates))),
        , drop = FALSE
      ]
      trans <- data.frame(
        type = "trans",
        gene_index = candidates$gene_index,
        variant_index = candidates$variant_index,
        effect_log2 = stats::rnorm(nrow(candidates), 0, trans_effect_sd),
        stringsAsFactors = FALSE
      )
    }
  }

  eqtls <- rbind(cis, trans)
  if (!nrow(eqtls)) {
    stop(
      "No eQTLs could be placed. Increase cis_window_bp or provide variants ",
      "across more genomic positions."
    )
  }
  n_gxe <- min(nrow(eqtls), round(nrow(eqtls) * gxe_fraction))
  eqtls$gxe_effect_log2 <- 0
  if (n_gxe > 0L) {
    gxe_index <- sample(seq_len(nrow(eqtls)), n_gxe)
    eqtls$gxe_effect_log2[gxe_index] <- stats::rnorm(n_gxe, 0, gxe_effect_sd)
  }
  eqtls$eqtl_id <- sprintf("eqtl_%05d", seq_len(nrow(eqtls)))
  eqtls$gene_id <- genes$gene_id[eqtls$gene_index]
  eqtls$variant_id <- variants$id[eqtls$variant_index]
  eqtls$seqname <- genes$seqname[eqtls$gene_index]
  eqtls$variant_seqname <- variants$seqname[eqtls$variant_index]
  eqtls$variant_pos <- variants$pos[eqtls$variant_index]
  eqtls$gene_tss <- genes$tss[eqtls$gene_index]
  eqtls$distance_bp <- ifelse(
    eqtls$seqname == eqtls$variant_seqname,
    eqtls$variant_pos - eqtls$gene_tss,
    NA_integer_
  )
  eqtls$gxe_condition <- ifelse(
    eqtls$gxe_effect_log2 != 0, treatment_label, NA_character_
  )
  eqtls$causal <- TRUE
  eqtls[, c(
    "eqtl_id", "type", "gene_id", "variant_id", "seqname",
    "variant_seqname", "variant_pos", "gene_tss", "distance_bp",
    "effect_log2", "gxe_effect_log2", "gxe_condition", "causal",
    "gene_index", "variant_index"
  )]
}

.simitall_fit_eqtl <- function(y, genotype, covariate_matrix, interaction = NULL) {
  x <- cbind(covariate_matrix, genotype = genotype)
  if (!is.null(interaction)) {
    x <- cbind(x, genotype_by_condition = genotype * interaction)
  }
  complete <- is.finite(y) & is.finite(genotype) & apply(is.finite(x), 1L, all)
  x <- x[complete, , drop = FALSE]
  y <- y[complete]
  if (length(y) <= ncol(x) || stats::var(genotype[complete]) == 0) {
    return(c(beta = NA, se = NA, p_value = NA, gxe_beta = NA, gxe_p_value = NA))
  }

  cross_inverse <- tryCatch(
    solve(crossprod(x)),
    error = function(error) NULL
  )
  if (is.null(cross_inverse)) {
    return(c(beta = NA, se = NA, p_value = NA, gxe_beta = NA, gxe_p_value = NA))
  }
  coefficients <- as.numeric(cross_inverse %*% crossprod(x, y))
  names(coefficients) <- colnames(x)
  residuals <- y - as.numeric(x %*% coefficients)
  degrees_freedom <- length(y) - ncol(x)
  sigma_squared <- sum(residuals^2) / degrees_freedom
  standard_errors <- sqrt(diag(cross_inverse) * sigma_squared)
  names(standard_errors) <- colnames(x)
  statistics <- coefficients / standard_errors
  p_values <- 2 * stats::pt(abs(statistics), df = degrees_freedom, lower.tail = FALSE)

  c(
    beta = unname(coefficients["genotype"]),
    se = unname(standard_errors["genotype"]),
    p_value = unname(p_values["genotype"]),
    gxe_beta = if ("genotype_by_condition" %in% names(coefficients)) {
      unname(coefficients["genotype_by_condition"])
    } else NA_real_,
    gxe_p_value = if ("genotype_by_condition" %in% names(p_values)) {
      unname(p_values["genotype_by_condition"])
    } else NA_real_
  )
}
