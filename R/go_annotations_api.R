#' Read a compact bundled gene-coordinate catalogue
#'
#' @param organism Optional organism name from the bundled catalogue, such as
#'   `"Zea_mays"` or `"Arabidopsis_thaliana"`.
#' @param file Optional path to a coordinate TSV or TSV.GZ. Supply either
#'   `organism` or `file`.
#' @return A data frame with one row per gene interval.
#' @export
read_gene_coordinates <- function(organism = NULL, file = NULL) {
  if (!is.null(file) && !is.null(organism)) {
    stop("Supply either organism or file, not both.", call. = FALSE)
  }
  if (is.null(file)) {
    if (is.null(organism) || !nzchar(organism)) {
      stop("Supply an organism from the bundled manifest or a coordinate file.", call. = FALSE)
    }
    manifest_path <- .simitall_extdata("gene_coordinates", "manifest.tsv")
    manifest <- utils::read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
    row <- manifest[tolower(manifest$organism) == tolower(organism), , drop = FALSE]
    if (nrow(row) != 1L) {
      stop("No unique bundled coordinate catalogue for organism: ", organism,
           ". See inst/extdata/gene_coordinates/manifest.tsv.", call. = FALSE)
    }
    file <- .simitall_extdata("gene_coordinates", row$coordinate_file[[1L]])
  }
  if (!is.character(file) || length(file) != 1L || !file.exists(file)) {
    stop("Coordinate file does not exist: ", file, call. = FALSE)
  }
  utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
}

.simitall_go_columns <- c(
  "gene_db", "go_object_id", "go_symbol", "qualifier", "go_id",
  "references", "evidence_code", "with_from", "aspect", "object_name",
  "synonyms", "object_type", "taxon", "annotation_date", "assigned_by",
  "annotation_extension", "gene_product_form_id"
)

.simitall_normalize_gene_id <- function(x) {
  x <- trimws(as.character(x))
  # GFF3 and GAF commonly prepend an identifier namespace (for example,
  # `gene:`, `AGI_LocusCode:`, or `FB:`). Coordinate catalogues retain the
  # unprefixed stable identifier for portable joins.
  sub("^[^:]+:", "", x, perl = TRUE)
}

#' Read Gene Ontology annotations in GAF format
#'
#' Reads a local `.gaf` or `.gaf.gz` file. GAF annotation files are external,
#' versioned resources; the result retains evidence codes and source metadata so
#' that downstream enrichment or reporting does not hide annotation provenance.
#'
#' @param file Path to a GAF or gzipped GAF file.
#' @param include_not Whether to retain annotations qualified with `NOT`.
#' @return A data frame in GAF 2.x column order with comment lines removed.
#' @export
read_go_gaf <- function(file, include_not = FALSE) {
  if (!is.character(file) || length(file) != 1L || !file.exists(file)) {
    stop("GAF file does not exist: ", file, call. = FALSE)
  }
  connection <- if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file, "rt") else file(file, "rt")
  on.exit(close(connection), add = TRUE)
  annotations <- utils::read.delim(
    connection, header = FALSE, sep = "\t", quote = "", comment.char = "!",
    fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE
  )
  if (ncol(annotations) < length(.simitall_go_columns)) {
    stop("GAF file has fewer than 17 required columns.", call. = FALSE)
  }
  annotations <- annotations[, seq_along(.simitall_go_columns), drop = FALSE]
  names(annotations) <- .simitall_go_columns
  if (!isTRUE(include_not)) {
    annotations <- annotations[!grepl("(^|\\|)NOT(\\||$)", annotations$qualifier), , drop = FALSE]
  }
  annotations
}

#' Join gene coordinates to Gene Ontology annotations
#'
#' @param gene_coordinates Data frame returned by [read_gene_coordinates()].
#' @param go_annotations Data frame returned by [read_go_gaf()].
#' @param match How to match IDs: `"auto"` tries stable gene ID then an exact
#'   symbol (the safe default); `"gene_id"` uses stable IDs only; `"gene_name"`
#'   uses symbols and listed synonyms and should be used only with a reviewed
#'   organism-specific ID mapping.
#' @return A data frame with matched gene intervals, GO terms, evidence codes,
#'   and `match_method`. Rows with no GO mapping are omitted.
#' @export
join_gene_go <- function(gene_coordinates, go_annotations,
                         match = c("auto", "gene_id", "gene_name")) {
  match <- match.arg(match)
  required_genes <- c("gene_id", "gene_name")
  required_go <- c("go_object_id", "go_symbol", "synonyms", "go_id", "evidence_code")
  if (!is.data.frame(gene_coordinates) || !all(required_genes %in% names(gene_coordinates))) {
    stop("gene_coordinates must contain gene_id and gene_name columns.", call. = FALSE)
  }
  if (!is.data.frame(go_annotations) || !all(required_go %in% names(go_annotations))) {
    stop("go_annotations must be the output of read_go_gaf().", call. = FALSE)
  }

  maps <- list()
  if (match %in% c("auto", "gene_id")) {
    maps[[length(maps) + 1L]] <- data.frame(
      key = .simitall_normalize_gene_id(go_annotations$go_object_id),
      go_row = seq_len(nrow(go_annotations)), match_method = "stable_id",
      stringsAsFactors = FALSE
    )
  }
  if (match %in% c("auto", "gene_name")) {
    maps[[length(maps) + 1L]] <- data.frame(
      key = .simitall_normalize_gene_id(go_annotations$go_symbol),
      go_row = seq_len(nrow(go_annotations)), match_method = "symbol",
      stringsAsFactors = FALSE
    )
    synonym_rows <- strsplit(go_annotations$synonyms, "|", fixed = TRUE)
    maps[[length(maps) + 1L]] <- data.frame(
      key = .simitall_normalize_gene_id(unlist(synonym_rows, use.names = FALSE)),
      go_row = rep.int(seq_len(nrow(go_annotations)), lengths(synonym_rows)),
      match_method = "synonym", stringsAsFactors = FALSE
    )
  }
  lookup <- do.call(rbind, maps)
  lookup <- lookup[nzchar(lookup$key), , drop = FALSE]
  lookup <- lookup[!duplicated(lookup[, c("key", "go_row", "match_method")]), , drop = FALSE]

  gene_ids <- data.frame(
    gene_row = seq_len(nrow(gene_coordinates)),
    key = .simitall_normalize_gene_id(gene_coordinates$gene_id),
    stringsAsFactors = FALSE
  )
  gene_ids <- gene_ids[nzchar(gene_ids$key), , drop = FALSE]
  stable_lookup <- lookup[lookup$match_method == "stable_id", , drop = FALSE]
  stable_hits <- if (match %in% c("auto", "gene_id")) {
    merge(gene_ids, stable_lookup, by = "key", all = FALSE, sort = FALSE)
  } else {
    data.frame(key = character(), gene_row = integer(), go_row = integer(),
               match_method = character(), stringsAsFactors = FALSE)
  }

  # A matching stable ID is decisive. Only fall back to symbol/synonym for a
  # gene with no stable-ID annotations, preventing a reused gene symbol from
  # attaching unrelated GO records to an otherwise well-mapped gene.
  use_names <- if (identical(match, "gene_name")) {
    seq_len(nrow(gene_coordinates))
  } else if (identical(match, "auto")) {
    setdiff(seq_len(nrow(gene_coordinates)), unique(stable_hits$gene_row))
  } else {
    integer()
  }
  name_hits <- stable_hits[FALSE, , drop = FALSE]
  if (length(use_names)) {
    gene_names <- data.frame(
      gene_row = use_names,
      key = .simitall_normalize_gene_id(gene_coordinates$gene_name[use_names]),
      stringsAsFactors = FALSE
    )
    gene_names <- gene_names[nzchar(gene_names$key), , drop = FALSE]
    # `auto` deliberately excludes loose synonyms: symbols are frequently
    # reused across organisms and historical records. Users who have reviewed
    # an organism-specific identifier map can opt into synonym matching with
    # `match = "gene_name"`.
    name_methods <- if (identical(match, "auto")) "symbol" else c("symbol", "synonym")
    name_lookup <- lookup[lookup$match_method %in% name_methods, , drop = FALSE]
    name_hits <- merge(gene_names, name_lookup, by = "key", all = FALSE, sort = FALSE)
  }
  hits <- rbind(stable_hits, name_hits)
  if (!nrow(hits)) {
    return(cbind(gene_coordinates[FALSE, , drop = FALSE], go_annotations[FALSE, , drop = FALSE], match_method = character()))
  }
  # A stable ID wins over a symbol/synonym when a GAF row can match more than once.
  rank <- match(hits$match_method, c("stable_id", "symbol", "synonym"))
  hits <- hits[order(hits$gene_row, hits$go_row, rank), , drop = FALSE]
  hits <- hits[!duplicated(hits[, c("gene_row", "go_row")]), , drop = FALSE]
  out <- cbind(
    gene_coordinates[hits$gene_row, , drop = FALSE],
    go_annotations[hits$go_row, , drop = FALSE],
    match_method = hits$match_method
  )
  rownames(out) <- NULL
  out
}

#' Summarize GO annotations after a gene join
#'
#' @param gene_go Output from [join_gene_go()].
#' @param by Grouping level: GO term only, or GO term plus ontology aspect.
#' @return A descending data frame of annotated-gene and annotation counts.
#' @export
summarize_go_terms <- function(gene_go, by = c("go_id", "go_id_aspect")) {
  by <- match.arg(by)
  if (!is.data.frame(gene_go) || !all(c("go_id", "gene_id", "aspect") %in% names(gene_go))) {
    stop("gene_go must be the output of join_gene_go().", call. = FALSE)
  }
  grouping <- if (identical(by, "go_id")) list(go_id = gene_go$go_id) else list(go_id = gene_go$go_id, aspect = gene_go$aspect)
  annotation_count <- stats::aggregate(gene_go$gene_id, grouping, length)
  names(annotation_count)[ncol(annotation_count)] <- "annotation_count"
  gene_count <- stats::aggregate(gene_go$gene_id, grouping, function(x) length(unique(x)))
  names(gene_count)[ncol(gene_count)] <- "gene_count"
  out <- merge(annotation_count, gene_count, by = names(grouping), sort = FALSE)
  out[order(out$gene_count, out$annotation_count, decreasing = TRUE), , drop = FALSE]
}
