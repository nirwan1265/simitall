#' Simulate donor-aware single-cell RNA-seq from GWAS genotypes
#'
#' Generate single-cell UMI counts for the same donors and genotypes used in a
#' GWAS cohort. The causal model combines cell-type marker programs, shared and
#' cell-type-specific cis/trans eQTLs, condition effects, genotype-by-condition
#' interactions, donor effects, trajectories, batches, latent factors, and
#' technical sampling. The output includes sparse 10x-style files, explicit
#' truth tables, donor-by-cell-type pseudobulk data, and an optional
#' `SingleCellExperiment` object.
#'
#' @param genotype_file Genotype TSV or VCF/VCF.GZ. Output from
#'   [simulate_gwas_cohort()] can be used directly.
#' @param out_prefix Prefix for single-cell outputs.
#' @param annotation_gff3 Optional gene annotation GFF3. Synthetic genes are
#'   placed along variant coordinates when omitted.
#' @param sample_metadata Optional data frame or TSV with a `sample` column.
#'   Donor-level population, family, phenotype, condition, and batch columns
#'   are copied into cell metadata.
#' @param n_genes Number of genes. With a GFF3, `NULL` uses every gene;
#'   otherwise the default is 1000.
#' @param cells_per_donor Number of cells per donor, as one value or one value
#'   per donor.
#' @param cell_type_proportions Named numeric vector of cell-type proportions.
#' @param trajectory_cell_types Cell types assigned continuous pseudotime.
#'   `NULL` uses the first cell type; `character(0)` disables trajectories.
#' @param ploidy Maximum genotype dosage.
#' @param n_cis_eqtl Number of causal cis-eQTL genes.
#' @param n_trans_eqtl Number of causal trans-eQTL pairs.
#' @param cis_window_bp Maximum cis distance from a gene TSS.
#' @param cell_type_eqtl_fraction Fraction of eQTLs active in only one cell
#'   type. Remaining eQTLs are shared across cell types.
#' @param condition_levels Donor condition labels used when absent from
#'   `sample_metadata`.
#' @param batch_levels Donor batch labels used when absent from metadata.
#' @param condition_effect_fraction Fraction of genes affected by condition.
#' @param gxe_fraction Fraction of eQTLs with genotype-by-condition effects.
#' @param markers_per_cell_type Number of non-overlapping marker genes assigned
#'   to each cell type.
#' @param marker_effect_mean Mean marker effect in log2 units.
#' @param marker_effect_sd Marker-effect standard deviation.
#' @param pseudotime_gene_fraction Fraction of genes with trajectory effects.
#' @param pseudotime_effect_sd Standard deviation of trajectory effects.
#' @param n_latent Number of cell-level latent factors.
#' @param baseline_mean_log2 Mean baseline gene abundance on the log2 scale.
#' @param baseline_sd_log2 Standard deviation of baseline gene abundance.
#' @param cis_effect_sd Standard deviation of cis-eQTL effects.
#' @param trans_effect_sd Standard deviation of trans-eQTL effects.
#' @param condition_effect_sd Standard deviation of condition effects.
#' @param gxe_effect_sd Standard deviation of genotype-by-condition effects.
#' @param donor_effect_sd Standard deviation of gene-by-donor effects.
#' @param batch_effect_sd Standard deviation of batch effects.
#' @param latent_effect_sd Standard deviation of latent-factor loadings.
#' @param residual_sd Standard deviation of cell-level log2 residual noise.
#' @param library_size_mean Mean UMI library size per cell.
#' @param library_size_sdlog Log-scale standard deviation of cell libraries.
#' @param dispersion_mean Mean negative-binomial dispersion where
#'   `variance = mean + dispersion * mean^2`.
#' @param dispersion_shape Gamma shape controlling gene dispersion variation.
#' @param dropout_rate Additional expression-dependent dropout probability.
#' @param doublet_rate Fraction of cells mixed with a second cell profile.
#' @param ambient_fraction Fraction of each expected library drawn from the
#'   global ambient expression pool.
#' @param backend Technical scaffold backend: `"auto"` uses Splatter when it
#'   is installed and otherwise uses the native gamma/negative-binomial model;
#'   `"splatter"` and `"native"` force a backend.
#' @param write_dense_tsv Also write a dense gene-by-cell count TSV. Sparse
#'   Matrix Market output and an RDS object are always written.
#' @param write_sce Write a `SingleCellExperiment` RDS when its Bioconductor
#'   package is installed.
#' @param min_cells_pseudobulk Minimum cells required for a donor-cell-type
#'   pseudobulk library.
#' @param seed Random-number seed.
#'
#' @return Invisibly returns a named list of output paths. The `backend`
#'   attribute records the backend actually used.
#' @examples
#' \dontrun{
#' sc <- simulate_scrnaseq_from_gwas(
#'   genotype_file = "results/gwas/demo.geno.tsv",
#'   out_prefix = "results/scrna/demo",
#'   sample_metadata = "results/gwas/demo.pheno.tsv",
#'   n_genes = 500,
#'   cells_per_donor = 100,
#'   cell_type_proportions = c(T_cell = 0.35, B_cell = 0.25,
#'                             Monocyte = 0.25, Stromal = 0.15),
#'   n_cis_eqtl = 100,
#'   n_trans_eqtl = 25,
#'   cell_type_eqtl_fraction = 0.5,
#'   seed = 1
#' )
#' }
#' @export
simulate_scrnaseq_from_gwas <- function(
    genotype_file,
    out_prefix,
    annotation_gff3 = NULL,
    sample_metadata = NULL,
    n_genes = NULL,
    cells_per_donor = 100L,
    cell_type_proportions = c(
      T_cell = 0.35,
      B_cell = 0.25,
      Monocyte = 0.25,
      Stromal = 0.15
    ),
    trajectory_cell_types = NULL,
    ploidy = 2L,
    n_cis_eqtl = 100L,
    n_trans_eqtl = 25L,
    cis_window_bp = 1e6,
    cell_type_eqtl_fraction = 0.4,
    condition_levels = c("control", "treatment"),
    batch_levels = c("batch1", "batch2"),
    condition_effect_fraction = 0.1,
    gxe_fraction = 0.2,
    markers_per_cell_type = NULL,
    marker_effect_mean = 1.5,
    marker_effect_sd = 0.3,
    pseudotime_gene_fraction = 0.1,
    pseudotime_effect_sd = 0.5,
    n_latent = 2L,
    baseline_mean_log2 = 5,
    baseline_sd_log2 = 1,
    cis_effect_sd = 0.6,
    trans_effect_sd = 0.35,
    condition_effect_sd = 0.6,
    gxe_effect_sd = 0.4,
    donor_effect_sd = 0.2,
    batch_effect_sd = 0.2,
    latent_effect_sd = 0.25,
    residual_sd = 0.1,
    library_size_mean = 5000,
    library_size_sdlog = 0.35,
    dispersion_mean = 0.1,
    dispersion_shape = 2,
    dropout_rate = 0.05,
    doublet_rate = 0.03,
    ambient_fraction = 0.02,
    backend = c("auto", "splatter", "native"),
    write_dense_tsv = FALSE,
    write_sce = TRUE,
    min_cells_pseudobulk = 5L,
    seed = 1L) {
  backend <- match.arg(backend)
  for (name in c(
    "cell_type_eqtl_fraction", "condition_effect_fraction", "gxe_fraction",
    "pseudotime_gene_fraction", "dropout_rate", "doublet_rate",
    "ambient_fraction"
  )) {
    .simitall_scrna_validate_probability(get(name), name)
  }
  if (length(ploidy) != 1L || !is.finite(ploidy) || ploidy < 1L ||
      ploidy != as.integer(ploidy)) {
    stop("ploidy must be a positive integer")
  }
  if (n_cis_eqtl < 0L || n_trans_eqtl < 0L || cis_window_bp < 0) {
    stop("eQTL counts and cis_window_bp must be non-negative")
  }
  if (library_size_mean <= 0 || library_size_sdlog < 0 ||
      dispersion_mean <= 0 || dispersion_shape <= 0) {
    stop("Library-size and dispersion parameters are invalid")
  }
  effect_scales <- c(
    marker_effect_sd = marker_effect_sd,
    pseudotime_effect_sd = pseudotime_effect_sd,
    baseline_sd_log2 = baseline_sd_log2,
    cis_effect_sd = cis_effect_sd,
    trans_effect_sd = trans_effect_sd,
    condition_effect_sd = condition_effect_sd,
    gxe_effect_sd = gxe_effect_sd,
    donor_effect_sd = donor_effect_sd,
    batch_effect_sd = batch_effect_sd,
    latent_effect_sd = latent_effect_sd,
    residual_sd = residual_sd
  )
  if (any(!is.finite(effect_scales)) || any(effect_scales < 0)) {
    stop("Effect and noise standard deviations must be non-negative")
  }
  if (!length(condition_levels) || !length(batch_levels)) {
    stop("condition_levels and batch_levels cannot be empty")
  }
  if (length(n_latent) != 1L || !is.finite(n_latent) || n_latent < 0L) {
    stop("n_latent must be a non-negative integer")
  }
  if (!length(cell_type_proportions) || any(!is.finite(cell_type_proportions)) ||
      any(cell_type_proportions <= 0)) {
    stop("cell_type_proportions must contain positive finite values")
  }
  if (is.null(names(cell_type_proportions)) ||
      any(!nzchar(names(cell_type_proportions))) ||
      anyDuplicated(names(cell_type_proportions))) {
    stop("cell_type_proportions must have unique non-empty names")
  }
  cell_type_proportions <- cell_type_proportions / sum(cell_type_proportions)
  if (is.null(trajectory_cell_types)) {
    trajectory_cell_types <- names(cell_type_proportions)[1L]
  }
  if (!all(trajectory_cell_types %in% names(cell_type_proportions))) {
    stop("trajectory_cell_types must be names in cell_type_proportions")
  }

  set.seed(seed)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  genotype_data <- .simitall_read_genotypes(genotype_file)
  variants <- genotype_data$variants
  genotype_observed <- genotype_data$genotype
  donor_ids <- colnames(genotype_observed)
  if (length(donor_ids) < 4L) {
    stop("At least four GWAS donors are required")
  }
  genotype <- .simitall_impute_genotypes(genotype_observed)
  if (max(genotype, na.rm = TRUE) > ploidy) {
    stop("Observed genotype dosage exceeds ploidy = ", ploidy)
  }

  donor_metadata <- .simitall_read_sample_metadata(sample_metadata, donor_ids)
  if (!"condition" %in% names(donor_metadata)) {
    donor_metadata$condition <- sample(
      rep(condition_levels, length.out = length(donor_ids)),
      length(donor_ids),
      replace = FALSE
    )
  }
  if (!"batch" %in% names(donor_metadata)) {
    donor_metadata$batch <- NA_character_
    for (condition in unique(donor_metadata$condition)) {
      index <- which(donor_metadata$condition == condition)
      donor_metadata$batch[index] <- sample(
        rep(batch_levels, length.out = length(index)),
        length(index),
        replace = FALSE
      )
    }
  }
  if (anyNA(donor_metadata$condition) || anyNA(donor_metadata$batch)) {
    stop("condition and batch labels cannot be missing")
  }
  donor_metadata$condition <- as.character(donor_metadata$condition)
  donor_metadata$batch <- as.character(donor_metadata$batch)
  observed_conditions <- unique(donor_metadata$condition)
  observed_batches <- unique(donor_metadata$batch)
  reference_condition <- observed_conditions[1L]
  treatment_label <- if (length(observed_conditions) > 1L) {
    observed_conditions[2L]
  } else {
    reference_condition
  }
  if (length(observed_conditions) < 2L) gxe_fraction <- 0

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
    genes <- .simitall_make_synthetic_genes(variants, n_genes)
  }
  if (!length(intersect(unique(variants$seqname), unique(genes$seqname))) &&
      length(unique(variants$seqname)) == 1L &&
      length(unique(genes$seqname)) == 1L) {
    variants$seqname <- unique(genes$seqname)
  }
  rownames(genes) <- genes$gene_id
  n_genes <- nrow(genes)
  if (n_genes < length(cell_type_proportions)) {
    stop("n_genes must be at least the number of cell types")
  }
  if (is.null(markers_per_cell_type)) {
    markers_per_cell_type <- max(3L, round(n_genes * 0.04))
  }
  if (length(markers_per_cell_type) != 1L ||
      !is.finite(markers_per_cell_type) || markers_per_cell_type < 1L) {
    stop("markers_per_cell_type must be at least one")
  }

  cells <- .simitall_scrna_make_cells(
    donor_metadata,
    cells_per_donor,
    cell_type_proportions,
    trajectory_cell_types
  )
  rownames(cells) <- cells$cell_id
  n_cells <- nrow(cells)
  donor_index <- match(cells$donor, donor_ids)
  treatment_indicator <- as.numeric(cells$condition != reference_condition)

  eqtls <- .simitall_select_eqtls(
    genes = genes,
    variants = variants,
    genotype = genotype,
    n_cis_eqtl = n_cis_eqtl,
    n_trans_eqtl = n_trans_eqtl,
    cis_window_bp = cis_window_bp,
    cis_effect_sd = cis_effect_sd,
    trans_effect_sd = trans_effect_sd,
    gxe_fraction = gxe_fraction,
    gxe_effect_sd = gxe_effect_sd,
    treatment_label = treatment_label
  )
  eqtls <- simulate_celltype_eqtls(
    eqtl_truth = eqtls,
    cell_types = names(cell_type_proportions),
    cell_type_specific_fraction = cell_type_eqtl_fraction
  )

  baseline <- stats::rnorm(n_genes, baseline_mean_log2, baseline_sd_log2)
  library_sizes <- pmax(100L, round(stats::rlnorm(
    n_cells,
    meanlog = log(library_size_mean) - 0.5 * library_size_sdlog^2,
    sdlog = library_size_sdlog
  )))
  scaffold <- .simitall_scrna_scaffold(
    backend,
    n_genes,
    n_cells,
    baseline,
    library_sizes,
    seed + 101L
  )
  backend_used <- scaffold$backend
  expected <- scaffold$expected
  baseline <- log2(pmax(rowMeans(expected), 1e-8))
  delta <- matrix(0, nrow = n_genes, ncol = n_cells)

  centered_genotype <- sweep(genotype, 1L, rowMeans(genotype), "-")
  for (i in seq_len(nrow(eqtls))) {
    active <- if (eqtls$scope[i] == "shared") {
      rep(TRUE, n_cells)
    } else {
      cells$cell_type == eqtls$cell_type[i]
    }
    dosage <- centered_genotype[eqtls$variant_index[i], donor_index]
    gene_index <- eqtls$gene_index[i]
    delta[gene_index, active] <- delta[gene_index, active] +
      eqtls$effect_log2[i] * dosage[active]
    if (eqtls$gxe_effect_log2[i] != 0) {
      delta[gene_index, active] <- delta[gene_index, active] +
        eqtls$gxe_effect_log2[i] * dosage[active] * treatment_indicator[active]
    }
  }

  marker_truth <- .simitall_scrna_marker_truth(
    genes,
    names(cell_type_proportions),
    markers_per_cell_type,
    marker_effect_mean,
    marker_effect_sd
  )
  for (i in seq_len(nrow(marker_truth))) {
    active <- cells$cell_type == marker_truth$cell_type[i]
    delta[marker_truth$gene_index[i], active] <-
      delta[marker_truth$gene_index[i], active] + marker_truth$effect_log2[i]
  }

  condition_truth <- data.frame(
    gene_id = character(),
    condition = character(),
    effect_log2 = numeric(),
    stringsAsFactors = FALSE
  )
  n_condition_genes <- min(n_genes, round(n_genes * condition_effect_fraction))
  if (n_condition_genes > 0L && length(observed_conditions) > 1L) {
    selected <- sample(seq_len(n_genes), n_condition_genes)
    for (condition in observed_conditions[-1L]) {
      effects <- stats::rnorm(n_condition_genes, 0, condition_effect_sd)
      active <- cells$condition == condition
      delta[selected, active] <- sweep(
        delta[selected, active, drop = FALSE], 1L, effects, "+"
      )
      condition_truth <- rbind(
        condition_truth,
        data.frame(
          gene_id = genes$gene_id[selected],
          condition = condition,
          effect_log2 = effects,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  batch_truth <- data.frame(
    gene_id = character(),
    batch = character(),
    effect_log2 = numeric(),
    stringsAsFactors = FALSE
  )
  if (length(observed_batches) > 1L) {
    for (batch in observed_batches[-1L]) {
      effects <- stats::rnorm(n_genes, 0, batch_effect_sd)
      active <- cells$batch == batch
      delta[, active] <- sweep(delta[, active, drop = FALSE], 1L, effects, "+")
      batch_truth <- rbind(
        batch_truth,
        data.frame(
          gene_id = genes$gene_id,
          batch = batch,
          effect_log2 = effects,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  pseudotime_truth <- data.frame(
    gene_id = character(),
    effect_log2 = numeric(),
    trajectory_cell_types = character(),
    stringsAsFactors = FALSE
  )
  n_pseudotime_genes <- min(n_genes, round(n_genes * pseudotime_gene_fraction))
  trajectory_index <- which(!is.na(cells$pseudotime))
  if (n_pseudotime_genes > 0L && length(trajectory_index)) {
    selected <- sample(seq_len(n_genes), n_pseudotime_genes)
    effects <- stats::rnorm(n_pseudotime_genes, 0, pseudotime_effect_sd)
    centered_time <- cells$pseudotime[trajectory_index] - 0.5
    delta[selected, trajectory_index] <- delta[selected, trajectory_index, drop = FALSE] +
      tcrossprod(effects, centered_time)
    pseudotime_truth <- data.frame(
      gene_id = genes$gene_id[selected],
      effect_log2 = effects,
      trajectory_cell_types = paste(trajectory_cell_types, collapse = ","),
      stringsAsFactors = FALSE
    )
  }

  donor_effects <- matrix(
    stats::rnorm(n_genes * length(donor_ids), 0, donor_effect_sd),
    nrow = n_genes,
    dimnames = list(genes$gene_id, donor_ids)
  )
  delta <- delta + donor_effects[, donor_index, drop = FALSE]

  n_latent <- max(0L, as.integer(n_latent))
  latent_scores <- matrix(
    numeric(),
    nrow = n_cells,
    ncol = 0L,
    dimnames = list(cells$cell_id, character())
  )
  latent_loadings <- matrix(
    numeric(),
    nrow = n_genes,
    ncol = 0L,
    dimnames = list(genes$gene_id, character())
  )
  if (n_latent > 0L) {
    latent_names <- paste0("latent_", seq_len(n_latent))
    latent_scores <- matrix(
      stats::rnorm(n_cells * n_latent),
      nrow = n_cells,
      dimnames = list(cells$cell_id, latent_names)
    )
    latent_loadings <- matrix(
      stats::rnorm(n_genes * n_latent, 0, latent_effect_sd),
      nrow = n_genes,
      dimnames = list(genes$gene_id, latent_names)
    )
    delta <- delta + latent_loadings %*% t(latent_scores)
  }
  if (residual_sd > 0) {
    delta <- delta + matrix(
      stats::rnorm(n_genes * n_cells, 0, residual_sd),
      nrow = n_genes
    )
  }
  delta <- pmin(20, pmax(-20, delta))
  expected <- expected * 2^delta
  expected <- sweep(expected, 2L, pmax(colSums(expected), 1e-8), "/")
  expected <- sweep(expected, 2L, library_sizes, "*")

  if (ambient_fraction > 0) {
    ambient_profile <- rowMeans(expected)
    ambient_profile <- ambient_profile / sum(ambient_profile)
    expected <- (1 - ambient_fraction) * expected +
      ambient_fraction * tcrossprod(ambient_profile, library_sizes)
  }

  cells$is_doublet <- FALSE
  cells$doublet_partner_cell_type <- NA_character_
  n_doublets <- min(n_cells, round(n_cells * doublet_rate))
  if (n_doublets > 0L) {
    doublets <- sample(seq_len(n_cells), n_doublets)
    partners <- sample(seq_len(n_cells), n_doublets, replace = TRUE)
    same <- partners == doublets
    while (any(same) && n_cells > 1L) {
      partners[same] <- sample(seq_len(n_cells), sum(same), replace = TRUE)
      same <- partners == doublets
    }
    expected[, doublets] <- 0.5 * (
      expected[, doublets, drop = FALSE] + expected[, partners, drop = FALSE]
    )
    cells$is_doublet[doublets] <- TRUE
    cells$doublet_partner_cell_type[doublets] <- cells$cell_type[partners]
  }

  dispersion <- stats::rgamma(
    n_genes,
    shape = dispersion_shape,
    rate = dispersion_shape / dispersion_mean
  )
  dispersion <- pmax(dispersion, 1e-8)
  counts <- matrix(
    stats::rnbinom(
      length(expected),
      mu = as.vector(expected),
      size = rep(1 / dispersion, times = n_cells)
    ),
    nrow = n_genes,
    dimnames = list(genes$gene_id, cells$cell_id)
  )
  technical_dropout_fraction <- 0
  if (dropout_rate > 0) {
    scale <- stats::median(expected[expected > 0])
    weights <- exp(-expected / max(scale, 1e-8))
    probability <- pmin(1, dropout_rate * weights / mean(weights))
    drop <- matrix(stats::runif(length(counts)), nrow = n_genes) < probability
    technical_dropout_fraction <- mean(drop & counts > 0)
    counts[drop] <- 0
  }
  observed_library_sizes <- colSums(counts)
  cells$target_library_size <- library_sizes
  cells$observed_library_size <- observed_library_sizes
  cells$detected_genes <- colSums(counts > 0)
  genes$baseline_log2 <- baseline
  genes$dispersion <- dispersion
  marker_output <- marker_truth[, setdiff(names(marker_truth), "gene_index")]
  eqtl_output <- eqtls[, setdiff(names(eqtls), c("gene_index", "variant_index"))]

  paths <- list(
    cell_metadata = paste0(out_prefix, ".cell_metadata.tsv"),
    donor_metadata = paste0(out_prefix, ".donor_metadata.tsv"),
    gene_metadata = paste0(out_prefix, ".gene_metadata.tsv"),
    eqtl_truth = paste0(out_prefix, ".celltype_eqtl_truth.tsv"),
    marker_truth = paste0(out_prefix, ".marker_truth.tsv"),
    condition_truth = paste0(out_prefix, ".condition_truth.tsv"),
    batch_truth = paste0(out_prefix, ".batch_truth.tsv"),
    pseudotime_truth = paste0(out_prefix, ".pseudotime_truth.tsv"),
    donor_effects = paste0(out_prefix, ".donor_effects.tsv"),
    latent_scores = paste0(out_prefix, ".latent_scores.tsv"),
    latent_loadings = paste0(out_prefix, ".latent_loadings.tsv"),
    object_rds = paste0(out_prefix, ".scrnaseq.rds"),
    summary = paste0(out_prefix, ".summary.json")
  )
  sparse_paths <- .simitall_scrna_write_sparse(counts, genes, cells, out_prefix)
  paths <- c(paths, sparse_paths)
  if (isTRUE(write_dense_tsv)) {
    paths$counts <- paste0(out_prefix, ".counts.tsv")
    .simitall_write_feature_matrix(counts, paths$counts)
  }
  utils::write.table(
    cells, paths$cell_metadata, sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    donor_metadata,
    paths$donor_metadata,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  for (item in list(
    list(genes, paths$gene_metadata),
    list(eqtl_output, paths$eqtl_truth),
    list(marker_output, paths$marker_truth),
    list(condition_truth, paths$condition_truth),
    list(batch_truth, paths$batch_truth),
    list(pseudotime_truth, paths$pseudotime_truth)
  )) {
    utils::write.table(
      item[[1L]], item[[2L]], sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
  .simitall_write_feature_matrix(donor_effects, paths$donor_effects)
  .simitall_write_feature_matrix(t(latent_scores), paths$latent_scores, "factor")
  .simitall_write_feature_matrix(latent_loadings, paths$latent_loadings)

  sparse_counts <- Matrix::Matrix(counts, sparse = TRUE)
  object <- structure(
    list(
      counts = sparse_counts,
      cell_metadata = cells,
      donor_metadata = donor_metadata,
      gene_metadata = genes,
      eqtl_truth = eqtl_output,
      marker_truth = marker_output,
      condition_truth = condition_truth,
      batch_truth = batch_truth,
      pseudotime_truth = pseudotime_truth,
      backend = backend_used
    ),
    class = c("simitall_scrnaseq", "list")
  )
  saveRDS(object, paths$object_rds, compress = "xz")

  if (isTRUE(write_sce) &&
      requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = sparse_counts),
      rowData = genes,
      colData = cells,
      metadata = list(
        backend = backend_used,
        eqtl_truth = eqtl_output,
        marker_truth = marker_output
      )
    )
    paths$sce_rds <- paste0(out_prefix, ".sce.rds")
    saveRDS(sce, paths$sce_rds, compress = "xz")
  } else if (isTRUE(write_sce)) {
    message(
      "SingleCellExperiment is not installed; the portable simitall RDS and ",
      "Matrix Market files were still written."
    )
  }

  pseudobulk <- aggregate_scrnaseq_pseudobulk(
    counts = counts,
    cell_metadata = cells,
    out_prefix = paste0(out_prefix, ".pseudobulk"),
    group_by = c("donor", "cell_type"),
    min_cells = min_cells_pseudobulk
  )
  paths$pseudobulk_counts <- pseudobulk$counts
  paths$pseudobulk_expression <- pseudobulk$expression
  paths$pseudobulk_metadata <- pseudobulk$metadata

  summary <- list(
    seed = seed,
    backend = backend_used,
    genotype_file = normalizePath(genotype_file, mustWork = TRUE),
    annotation_gff3 = if (is.null(annotation_gff3)) NULL else {
      normalizePath(annotation_gff3, mustWork = TRUE)
    },
    counts = list(
      donors = length(donor_ids),
      cells = n_cells,
      genes = n_genes,
      cell_types = length(cell_type_proportions),
      cis_eqtls = sum(eqtls$type == "cis"),
      trans_eqtls = sum(eqtls$type == "trans"),
      cell_type_specific_eqtls = sum(eqtls$scope == "cell_type_specific"),
      marker_pairs = nrow(marker_truth),
      doublets = sum(cells$is_doublet),
      zero_fraction = mean(counts == 0),
      technical_dropout_fraction = technical_dropout_fraction
    ),
    parameters = list(
      cells_per_donor = cells_per_donor,
      cell_type_proportions = as.list(cell_type_proportions),
      trajectory_cell_types = trajectory_cell_types,
      cis_window_bp = cis_window_bp,
      cell_type_eqtl_fraction = cell_type_eqtl_fraction,
      condition_effect_fraction = condition_effect_fraction,
      gxe_fraction = gxe_fraction,
      dropout_rate = dropout_rate,
      doublet_rate = doublet_rate,
      ambient_fraction = ambient_fraction,
      library_size_mean = library_size_mean,
      dispersion_mean = dispersion_mean,
      ploidy = ploidy
    )
  )
  jsonlite::write_json(summary, paths$summary, pretty = TRUE, auto_unbox = TRUE)
  attr(paths, "backend") <- backend_used
  message(
    "Single-cell simulation complete: ", length(donor_ids), " donors, ",
    n_cells, " cells, ", n_genes, " genes, and ", nrow(eqtls),
    " causal eQTLs (", backend_used, " backend)."
  )
  invisible(paths)
}

#' Assign shared and cell-type-specific activity to causal eQTLs
#'
#' Convert a bulk or generic causal eQTL truth table into a single-cell truth
#' table by marking each pair as shared across all cell types or active in one
#' selected cell type. This helper is used internally by
#' [simulate_scrnaseq_from_gwas()] and can also adapt user-supplied eQTL panels.
#'
#' @param eqtl_truth eQTL truth data frame or TSV. It must contain `gene_id`,
#'   `variant_id`, and `effect_log2` columns; other columns are preserved.
#' @param cell_types Character vector of cell-type labels.
#' @param cell_type_specific_fraction Fraction of pairs assigned to one cell
#'   type. Remaining pairs are labelled `"shared"` and active in all types.
#' @param out_tsv Optional output TSV path.
#' @param seed Optional random-number seed. `NULL` uses the current RNG stream.
#'
#' @return The augmented eQTL truth data frame, invisibly when written.
#' @examples
#' truth <- data.frame(
#'   gene_id = paste0("gene", 1:4),
#'   variant_id = paste0("variant", 1:4),
#'   effect_log2 = c(0.5, -0.4, 0.7, -0.2)
#' )
#' cell_truth <- simulate_celltype_eqtls(
#'   truth,
#'   cell_types = c("T_cell", "B_cell", "Monocyte"),
#'   cell_type_specific_fraction = 0.5,
#'   seed = 1
#' )
#' @export
simulate_celltype_eqtls <- function(
    eqtl_truth,
    cell_types,
    cell_type_specific_fraction = 0.4,
    out_tsv = NULL,
    seed = NULL) {
  .simitall_scrna_validate_probability(
    cell_type_specific_fraction,
    "cell_type_specific_fraction"
  )
  if (!length(cell_types) || anyNA(cell_types) || any(!nzchar(cell_types)) ||
      anyDuplicated(cell_types)) {
    stop("cell_types must contain unique non-empty labels")
  }
  truth <- .simitall_read_table_or_data(eqtl_truth, "eQTL truth")
  required <- c("gene_id", "variant_id", "effect_log2")
  if (!all(required %in% names(truth))) {
    stop(
      "eQTL truth is missing: ",
      paste(setdiff(required, names(truth)), collapse = ", ")
    )
  }
  if (!nrow(truth)) stop("eQTL truth contains no causal pairs")
  if (!is.null(seed)) set.seed(seed)
  truth$scope <- "shared"
  truth$cell_type <- "all"
  n_specific <- min(
    nrow(truth),
    round(nrow(truth) * cell_type_specific_fraction)
  )
  if (n_specific > 0L) {
    specific <- sample(seq_len(nrow(truth)), n_specific)
    truth$scope[specific] <- "cell_type_specific"
    truth$cell_type[specific] <- sample(
      cell_types, n_specific, replace = TRUE
    )
  }
  if (!is.null(out_tsv)) {
    dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      truth, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE
    )
    return(invisible(truth))
  }
  truth
}

#' Aggregate single-cell counts into donor-level pseudobulk libraries
#'
#' Sum cells within donor, cell type, condition, or other metadata groups and
#' write raw counts plus log2 counts-per-million expression. Pseudobulk keeps
#' donors, rather than cells, as the biological replication unit for eQTL and
#' differential-expression analyses.
#'
#' @param counts Gene-by-cell matrix, dense count TSV, single-cell result
#'   prefix, or `.scrnaseq.rds` path.
#' @param cell_metadata Cell metadata data frame or TSV. May be omitted when
#'   `counts` identifies a single-cell result object.
#' @param out_prefix Output prefix.
#' @param group_by Cell metadata columns defining each pseudobulk library.
#' @param min_cells Minimum cells required in a retained library.
#'
#' @return Invisibly returns paths to count, expression, and metadata TSVs.
#' @examples
#' \dontrun{
#' aggregate_scrnaseq_pseudobulk(
#'   counts = "results/scrna/demo",
#'   out_prefix = "results/scrna/demo.pseudobulk",
#'   group_by = c("donor", "cell_type")
#' )
#' }
#' @export
aggregate_scrnaseq_pseudobulk <- function(
    counts,
    cell_metadata = NULL,
    out_prefix,
    group_by = c("donor", "cell_type"),
    min_cells = 5L) {
  if (length(min_cells) != 1L || !is.finite(min_cells) || min_cells < 1L) {
    stop("min_cells must be at least one")
  }
  if (is.character(counts) && length(counts) == 1L &&
      (grepl("\\.scrnaseq\\.rds$", counts) ||
       file.exists(paste0(counts, ".scrnaseq.rds")))) {
    object <- .simitall_read_scrna_object(counts)
    counts <- object$counts
    if (is.null(cell_metadata)) cell_metadata <- object$cell_metadata
  } else if (is.character(counts) && length(counts) == 1L) {
    counts <- .simitall_read_expression_matrix(counts)
  } else {
    counts <- as.matrix(counts)
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("counts must have gene and cell names")
  }
  if (is.null(cell_metadata)) {
    stop("cell_metadata is required unless counts identifies a result object")
  }
  cell_metadata <- .simitall_read_table_or_data(cell_metadata, "Cell metadata")
  if (!"cell_id" %in% names(cell_metadata)) {
    stop("cell_metadata must contain a cell_id column")
  }
  if (!all(group_by %in% names(cell_metadata))) {
    stop(
      "cell_metadata is missing grouping columns: ",
      paste(setdiff(group_by, names(cell_metadata)), collapse = ", ")
    )
  }
  shared <- intersect(colnames(counts), cell_metadata$cell_id)
  if (!length(shared)) stop("No cell IDs are shared by counts and metadata")
  counts <- counts[, shared, drop = FALSE]
  cell_metadata <- cell_metadata[
    match(shared, cell_metadata$cell_id), , drop = FALSE
  ]
  group_id <- do.call(
    paste,
    c(cell_metadata[, group_by, drop = FALSE], sep = "__")
  )
  group_levels <- unique(group_id)
  group_counts <- table(factor(group_id, levels = group_levels))
  keep_levels <- group_levels[group_counts >= min_cells]
  keep <- group_id %in% keep_levels
  if (!any(keep)) {
    stop("No pseudobulk libraries meet min_cells = ", min_cells)
  }
  group_factor <- factor(group_id[keep], levels = keep_levels)
  design <- Matrix::sparseMatrix(
    i = seq_along(group_factor),
    j = as.integer(group_factor),
    x = 1,
    dims = c(length(group_factor), length(keep_levels))
  )
  pseudobulk_counts <- as.matrix(counts[, keep, drop = FALSE] %*% design)
  colnames(pseudobulk_counts) <- keep_levels
  library_sizes <- colSums(pseudobulk_counts)
  library_sizes[library_sizes < 1] <- 1
  expression <- log2(
    sweep(pseudobulk_counts, 2L, library_sizes, "/") * 1e6 + 1
  )

  first_index <- match(keep_levels, group_id)
  metadata_columns <- setdiff(
    names(cell_metadata),
    c(
      "cell_id", "pseudotime", "target_library_size",
      "observed_library_size", "detected_genes", "is_doublet",
      "doublet_partner_cell_type"
    )
  )
  pseudobulk_metadata <- cell_metadata[
    first_index, metadata_columns, drop = FALSE
  ]
  pseudobulk_metadata <- data.frame(
    pseudobulk_id = keep_levels,
    pseudobulk_metadata,
    n_cells = as.integer(group_counts[keep_levels]),
    library_size = as.numeric(library_sizes[keep_levels]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  paths <- list(
    counts = paste0(out_prefix, ".counts.tsv"),
    expression = paste0(out_prefix, ".expression.tsv"),
    metadata = paste0(out_prefix, ".metadata.tsv")
  )
  .simitall_write_feature_matrix(pseudobulk_counts, paths$counts)
  .simitall_write_feature_matrix(expression, paths$expression)
  utils::write.table(
    pseudobulk_metadata,
    paths$metadata,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  invisible(paths)
}

#' Benchmark cell-type-specific eQTL recovery
#'
#' Test eQTLs separately in each cell type using donor-level pseudobulk
#' expression. Shared eQTLs are causal in every cell type; cell-type-specific
#' eQTLs are causal only in their assigned cell type. This avoids treating
#' individual cells as independent biological replicates.
#'
#' @param genotype_file Genotype TSV or VCF/VCF.GZ.
#' @param pseudobulk_expression Gene-by-pseudobulk log2-expression matrix or
#'   TSV from [aggregate_scrnaseq_pseudobulk()].
#' @param pseudobulk_metadata Metadata data frame or TSV with `pseudobulk_id`,
#'   `donor`, and `cell_type` columns.
#' @param truth_file Cell-type eQTL truth table from
#'   [simulate_scrnaseq_from_gwas()].
#' @param out_prefix Output prefix.
#' @param cis_window_bp Cis candidate window.
#' @param covariates Donor metadata columns included in regression models.
#' @param test_gxe Test genotype-by-condition interactions.
#' @param fdr_threshold BH-adjusted significance threshold.
#' @param min_donors Minimum donors required per cell type.
#' @param max_tests Maximum association tests per cell type.
#'
#' @return Invisibly returns combined result and metric paths and tables.
#' @examples
#' \dontrun{
#' benchmark_celltype_eqtls(
#'   genotype_file = "results/gwas/demo.geno.tsv",
#'   pseudobulk_expression = sc$pseudobulk_expression,
#'   pseudobulk_metadata = sc$pseudobulk_metadata,
#'   truth_file = sc$eqtl_truth,
#'   out_prefix = "results/scrna/demo_celltype_eqtl"
#' )
#' }
#' @export
benchmark_celltype_eqtls <- function(
    genotype_file,
    pseudobulk_expression,
    pseudobulk_metadata,
    truth_file,
    out_prefix,
    cis_window_bp = 1e6,
    covariates = c("condition", "batch", "pop"),
    test_gxe = TRUE,
    fdr_threshold = 0.05,
    min_donors = 8L,
    max_tests = 1e6) {
  expression <- .simitall_read_expression_matrix(pseudobulk_expression)
  metadata <- .simitall_read_table_or_data(
    pseudobulk_metadata, "Pseudobulk metadata"
  )
  truth <- .simitall_read_table_or_data(truth_file, "eQTL truth")
  required_metadata <- c("pseudobulk_id", "donor", "cell_type")
  if (!all(required_metadata %in% names(metadata))) {
    stop(
      "Pseudobulk metadata is missing: ",
      paste(setdiff(required_metadata, names(metadata)), collapse = ", ")
    )
  }
  required_truth <- c(
    "gene_id", "variant_id", "seqname", "gene_tss", "scope", "cell_type"
  )
  if (!all(required_truth %in% names(truth))) {
    stop(
      "eQTL truth is missing: ",
      paste(setdiff(required_truth, names(truth)), collapse = ", ")
    )
  }
  if (fdr_threshold <= 0 || fdr_threshold > 1) {
    stop("fdr_threshold must be greater than zero and at most one")
  }
  cell_types <- unique(metadata$cell_type)
  results_list <- list()
  metrics_list <- list()
  result_number <- 0L

  for (cell_type in cell_types) {
    cell_metadata <- metadata[metadata$cell_type == cell_type, , drop = FALSE]
    shared_libraries <- intersect(cell_metadata$pseudobulk_id, colnames(expression))
    cell_metadata <- cell_metadata[
      match(shared_libraries, cell_metadata$pseudobulk_id), , drop = FALSE
    ]
    if (length(unique(cell_metadata$donor)) < min_donors) {
      warning(
        "Skipping ", cell_type, ": fewer than ", min_donors, " donors",
        call. = FALSE
      )
      next
    }
    truth_for_type <- truth[
      truth$scope == "shared" |
        (truth$scope == "cell_type_specific" & truth$cell_type == cell_type),
      , drop = FALSE
    ]
    if (!nrow(truth_for_type)) next
    cell_expression <- expression[, shared_libraries, drop = FALSE]
    colnames(cell_expression) <- cell_metadata$donor
    analysis_metadata <- cell_metadata
    analysis_metadata$sample <- analysis_metadata$donor
    analysis_metadata <- analysis_metadata[
      !duplicated(analysis_metadata$sample), , drop = FALSE
    ]
    cell_expression <- cell_expression[
      , analysis_metadata$sample, drop = FALSE
    ]
    usable_covariates <- intersect(covariates, names(analysis_metadata))
    usable_covariates <- usable_covariates[vapply(
      analysis_metadata[, usable_covariates, drop = FALSE],
      function(value) length(unique(value[!is.na(value)])) > 1L,
      logical(1L)
    )]
    temporary_prefix <- tempfile(paste0("simitall-", cell_type, "-"))
    fit <- benchmark_eqtl(
      genotype_file = genotype_file,
      expression_file = cell_expression,
      out_prefix = temporary_prefix,
      sample_metadata = analysis_metadata,
      truth_file = truth_for_type,
      pair_mode = "truth_neighborhood",
      cis_window_bp = cis_window_bp,
      covariates = usable_covariates,
      test_gxe = test_gxe,
      expression_scale = "log2",
      fdr_threshold = fdr_threshold,
      max_tests = max_tests
    )
    unlink(paste0(temporary_prefix, c(".eqtl_results.tsv", ".eqtl_metrics.tsv")))
    result_number <- result_number + 1L
    fit$results$cell_type <- cell_type
    fit$metrics$cell_type <- cell_type
    results_list[[result_number]] <- fit$results
    metrics_list[[result_number]] <- fit$metrics
  }
  if (!length(results_list)) {
    stop("No cell type had enough donors and causal truth for benchmarking")
  }
  results <- do.call(rbind, results_list)
  metrics <- do.call(rbind, metrics_list)
  overall <- data.frame(
    tests = sum(metrics$tests),
    truth_pairs_in_tests = sum(metrics$truth_pairs_in_tests),
    discoveries = sum(metrics$discoveries),
    true_positives = sum(metrics$true_positives),
    false_positives = sum(metrics$false_positives),
    false_negatives = sum(metrics$false_negatives),
    precision = if (sum(metrics$discoveries)) {
      sum(metrics$true_positives) / sum(metrics$discoveries)
    } else NA_real_,
    recall = if (sum(metrics$truth_pairs_in_tests)) {
      sum(metrics$true_positives) / sum(metrics$truth_pairs_in_tests)
    } else NA_real_,
    empirical_fdr = if (sum(metrics$discoveries)) {
      sum(metrics$false_positives) / sum(metrics$discoveries)
    } else NA_real_,
    fdr_threshold = fdr_threshold,
    cell_type = "Overall",
    stringsAsFactors = FALSE
  )
  metrics <- rbind(metrics, overall)

  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  paths <- list(
    results = paste0(out_prefix, ".celltype_eqtl_results.tsv"),
    metrics = paste0(out_prefix, ".celltype_eqtl_metrics.tsv")
  )
  utils::write.table(
    results, paths$results, sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    metrics, paths$metrics, sep = "\t", row.names = FALSE, quote = FALSE
  )
  message(
    "Cell-type eQTL benchmark complete: ", nrow(results), " tests across ",
    length(results_list), " cell types."
  )
  invisible(list(
    results_path = paths$results,
    metrics_path = paths$metrics,
    results = results,
    metrics = metrics
  ))
}
