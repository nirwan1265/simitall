#' Simulate a standalone bulk RNA-seq experiment
#'
#' Generate gene-level RNA-seq counts without a GWAS genotype file. Samples can
#' have technical replicates that share one biological profile, independent
#' biological replicates, or mixed biological and technical replication. The
#' model includes condition, sample, biological-replicate, batch, latent, and
#' technical-library effects with explicit truth outputs.
#'
#' @param out_prefix Prefix for count, metadata, truth, and summary files.
#' @param sample_names Names of the underlying samples or experimental groups.
#' @param replicates Number of replicates per sample. In `"technical"` mode
#'   this is the number of technical libraries; in `"biological"` and
#'   `"mixed"` modes it is the number of biological replicates.
#' @param replicate_mode `"technical"`, `"biological"`, or `"mixed"`.
#' @param technical_replicates Number of technical libraries per biological
#'   replicate in `"mixed"` mode.
#' @param conditions Condition label per sample, or one shared label.
#' @param batch_levels Labels used to assign libraries to batches.
#' @param design Optional aggregated data frame or TSV overriding the arguments
#'   above. It requires `sample` and may contain `condition`, `batch`,
#'   `biological_replicates`, and `technical_replicates` columns.
#' @param annotation_gff3 Optional GFF3 supplying gene coordinates.
#' @param n_genes Number of genes. With a GFF3, `NULL` uses every gene;
#'   otherwise the default is 1000.
#' @param condition_effect_fraction Fraction of genes affected by each
#'   non-reference condition.
#' @param baseline_mean_log2 Mean baseline log2 abundance.
#' @param baseline_sd_log2 Standard deviation of baseline abundance.
#' @param sample_effect_sd Standard deviation of source-sample effects.
#' @param biological_replicate_sd Standard deviation of biological-replicate
#'   effects. These are shared by all technical libraries from that replicate.
#' @param condition_effect_sd Standard deviation of condition effects.
#' @param batch_effect_sd Standard deviation of batch effects.
#' @param technical_effect_sd Standard deviation of library-specific technical
#'   effects.
#' @param n_latent Number of latent biological factors.
#' @param latent_effect_sd Standard deviation of latent-factor loadings.
#' @param residual_sd Additional library-level log2 noise.
#' @param library_size_mean Mean count library size.
#' @param library_size_sdlog Log-scale standard deviation of library sizes.
#' @param dispersion_mean Mean negative-binomial dispersion.
#' @param dispersion_shape Gamma shape controlling dispersion variation.
#' @param simulate_reads Also generate FASTQ files using
#'   [simulate_rnaseq_reads()].
#' @param transcript_fasta Transcript FASTA required for FASTQ simulation.
#' @param transcript_gene_map Optional transcript-to-gene map.
#' @param read_library_size Reads or read pairs per library.
#' @param read_length Read length.
#' @param paired_end Generate paired-end reads.
#' @param seed Random-number seed.
#'
#' @return Invisibly returns a named list of output paths.
#' @examples
#' \dontrun{
#' simulate_rnaseq_experiment(
#'   out_prefix = "results/rnaseq/sample_A",
#'   sample_names = "sample_A",
#'   replicates = 3,
#'   replicate_mode = "technical",
#'   n_genes = 500,
#'   seed = 1
#' )
#'
#' simulate_rnaseq_experiment(
#'   out_prefix = "results/rnaseq/treatment",
#'   sample_names = c("control", "treated"),
#'   replicates = 3,
#'   replicate_mode = "biological",
#'   conditions = c("control", "treated"),
#'   seed = 2
#' )
#' }
#' @export
simulate_rnaseq_experiment <- function(
    out_prefix,
    sample_names = "sample_A",
    replicates = 3L,
    replicate_mode = c("technical", "biological", "mixed"),
    technical_replicates = 2L,
    conditions = "control",
    batch_levels = "batch1",
    design = NULL,
    annotation_gff3 = NULL,
    n_genes = NULL,
    condition_effect_fraction = 0.1,
    baseline_mean_log2 = 6,
    baseline_sd_log2 = 1,
    sample_effect_sd = 0.25,
    biological_replicate_sd = 0.2,
    condition_effect_sd = 0.6,
    batch_effect_sd = 0.15,
    technical_effect_sd = 0.05,
    n_latent = 2L,
    latent_effect_sd = 0.2,
    residual_sd = 0.05,
    library_size_mean = 5e6,
    library_size_sdlog = 0.25,
    dispersion_mean = 0.1,
    dispersion_shape = 2,
    simulate_reads = FALSE,
    transcript_fasta = NULL,
    transcript_gene_map = NULL,
    read_library_size = 1e5,
    read_length = 100L,
    paired_end = TRUE,
    seed = 1L) {
  replicate_mode <- match.arg(replicate_mode)
  .simitall_scrna_validate_probability(
    condition_effect_fraction, "condition_effect_fraction"
  )
  .simitall_validate_effect_scales(c(
    baseline_sd_log2 = baseline_sd_log2,
    sample_effect_sd = sample_effect_sd,
    biological_replicate_sd = biological_replicate_sd,
    condition_effect_sd = condition_effect_sd,
    batch_effect_sd = batch_effect_sd,
    technical_effect_sd = technical_effect_sd,
    latent_effect_sd = latent_effect_sd,
    residual_sd = residual_sd
  ))
  if (library_size_mean <= 0 || library_size_sdlog < 0 ||
      dispersion_mean <= 0 || dispersion_shape <= 0) {
    stop("Library-size and dispersion parameters are invalid")
  }
  if (length(n_latent) != 1L || !is.finite(n_latent) || n_latent < 0L) {
    stop("n_latent must be a non-negative integer")
  }

  set.seed(seed)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  metadata <- .simitall_expand_experiment_design(
    sample_names = sample_names,
    replicates = replicates,
    technical_replicates = technical_replicates,
    replicate_mode = replicate_mode,
    conditions = conditions,
    batch_levels = batch_levels,
    design = design
  )
  genes <- .simitall_standalone_genes(annotation_gff3, n_genes)
  n_genes <- nrow(genes)
  biological_ids <- unique(metadata$biological_id)
  source_samples <- unique(metadata$sample)
  conditions_used <- unique(metadata$condition)
  batches_used <- unique(metadata$batch)

  baseline <- stats::rnorm(n_genes, baseline_mean_log2, baseline_sd_log2)
  biological_eta <- matrix(
    baseline,
    nrow = n_genes,
    ncol = length(biological_ids),
    dimnames = list(genes$gene_id, biological_ids)
  )
  biological_metadata <- metadata[
    match(biological_ids, metadata$biological_id),
    , drop = FALSE
  ]
  sample_effects <- matrix(
    stats::rnorm(n_genes * length(source_samples), 0, sample_effect_sd),
    nrow = n_genes,
    dimnames = list(genes$gene_id, source_samples)
  )
  biological_effects <- matrix(
    stats::rnorm(
      n_genes * length(biological_ids), 0, biological_replicate_sd
    ),
    nrow = n_genes,
    dimnames = list(genes$gene_id, biological_ids)
  )
  biological_eta <- biological_eta +
    sample_effects[, biological_metadata$sample, drop = FALSE] +
    biological_effects

  condition_truth <- data.frame(
    gene_id = character(),
    condition = character(),
    effect_log2 = numeric(),
    stringsAsFactors = FALSE
  )
  n_condition_genes <- min(
    n_genes, round(n_genes * condition_effect_fraction)
  )
  if (n_condition_genes > 0L && length(conditions_used) > 1L) {
    selected <- sample(seq_len(n_genes), n_condition_genes)
    for (condition in conditions_used[-1L]) {
      effects <- stats::rnorm(n_condition_genes, 0, condition_effect_sd)
      active <- biological_metadata$condition == condition
      biological_eta[selected, active] <- sweep(
        biological_eta[selected, active, drop = FALSE], 1L, effects, "+"
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

  n_latent <- as.integer(n_latent)
  latent_scores <- matrix(
    numeric(), nrow = length(biological_ids), ncol = 0L,
    dimnames = list(biological_ids, character())
  )
  latent_loadings <- matrix(
    numeric(), nrow = n_genes, ncol = 0L,
    dimnames = list(genes$gene_id, character())
  )
  if (n_latent > 0L) {
    latent_names <- paste0("latent_", seq_len(n_latent))
    latent_scores <- matrix(
      stats::rnorm(length(biological_ids) * n_latent),
      nrow = length(biological_ids),
      dimnames = list(biological_ids, latent_names)
    )
    latent_loadings <- matrix(
      stats::rnorm(n_genes * n_latent, 0, latent_effect_sd),
      nrow = n_genes,
      dimnames = list(genes$gene_id, latent_names)
    )
    biological_eta <- biological_eta + latent_loadings %*% t(latent_scores)
  }

  library_eta <- biological_eta[, metadata$biological_id, drop = FALSE]
  colnames(library_eta) <- metadata$library_id
  batch_truth <- data.frame(
    gene_id = character(), batch = character(), effect_log2 = numeric(),
    stringsAsFactors = FALSE
  )
  if (length(batches_used) > 1L) {
    for (batch in batches_used[-1L]) {
      effects <- stats::rnorm(n_genes, 0, batch_effect_sd)
      active <- metadata$batch == batch
      library_eta[, active] <- sweep(
        library_eta[, active, drop = FALSE], 1L, effects, "+"
      )
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
  technical_effects <- matrix(
    stats::rnorm(n_genes * nrow(metadata), 0, technical_effect_sd),
    nrow = n_genes,
    dimnames = list(genes$gene_id, metadata$library_id)
  )
  library_eta <- library_eta + technical_effects
  if (residual_sd > 0) {
    library_eta <- library_eta + matrix(
      stats::rnorm(n_genes * nrow(metadata), 0, residual_sd),
      nrow = n_genes
    )
  }
  # Clamp extreme values without dropping matrix dimensions for one-sample runs.
  library_eta[library_eta < -20] <- -20
  library_eta[library_eta > 30] <- 30
  abundance <- 2^library_eta
  relative_abundance <- sweep(abundance, 2L, colSums(abundance), "/")
  library_sizes <- pmax(1L, round(stats::rlnorm(
    nrow(metadata),
    meanlog = log(library_size_mean) - 0.5 * library_size_sdlog^2,
    sdlog = library_size_sdlog
  )))
  expected_counts <- sweep(relative_abundance, 2L, library_sizes, "*")
  dispersion <- pmax(
    stats::rgamma(
      n_genes,
      shape = dispersion_shape,
      rate = dispersion_shape / dispersion_mean
    ),
    1e-8
  )
  counts <- matrix(
    stats::rnbinom(
      length(expected_counts),
      mu = as.vector(expected_counts),
      size = rep(1 / dispersion, times = nrow(metadata))
    ),
    nrow = n_genes,
    dimnames = dimnames(expected_counts)
  )
  observed_library_sizes <- pmax(1, colSums(counts))
  expression <- log2(
    sweep(counts, 2L, observed_library_sizes, "/") * 1e6 + 1
  )
  metadata$target_library_size <- library_sizes
  metadata$observed_library_size <- observed_library_sizes
  genes$baseline_log2 <- baseline
  genes$dispersion <- dispersion

  paths <- list(
    counts = paste0(out_prefix, ".counts.tsv"),
    expression = paste0(out_prefix, ".expression.tsv"),
    expected_counts = paste0(out_prefix, ".expected_counts.tsv"),
    biological_expression = paste0(out_prefix, ".biological_expression.tsv"),
    sample_metadata = paste0(out_prefix, ".sample_metadata.tsv"),
    gene_metadata = paste0(out_prefix, ".gene_metadata.tsv"),
    condition_truth = paste0(out_prefix, ".condition_truth.tsv"),
    batch_truth = paste0(out_prefix, ".batch_truth.tsv"),
    sample_effects = paste0(out_prefix, ".sample_effects.tsv"),
    biological_effects = paste0(out_prefix, ".biological_effects.tsv"),
    technical_effects = paste0(out_prefix, ".technical_effects.tsv"),
    latent_scores = paste0(out_prefix, ".latent_scores.tsv"),
    latent_loadings = paste0(out_prefix, ".latent_loadings.tsv"),
    summary = paste0(out_prefix, ".summary.json")
  )
  .simitall_write_feature_matrix(counts, paths$counts)
  .simitall_write_feature_matrix(expression, paths$expression)
  .simitall_write_feature_matrix(expected_counts, paths$expected_counts)
  .simitall_write_feature_matrix(
    biological_eta, paths$biological_expression
  )
  .simitall_write_feature_matrix(sample_effects, paths$sample_effects)
  .simitall_write_feature_matrix(
    biological_effects, paths$biological_effects
  )
  .simitall_write_feature_matrix(
    technical_effects, paths$technical_effects
  )
  .simitall_write_feature_matrix(
    t(latent_scores), paths$latent_scores, "factor"
  )
  .simitall_write_feature_matrix(
    latent_loadings, paths$latent_loadings
  )
  utils::write.table(
    metadata, paths$sample_metadata, sep = "\t", row.names = FALSE,
    quote = FALSE
  )
  utils::write.table(
    genes, paths$gene_metadata, sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    condition_truth, paths$condition_truth, sep = "\t", row.names = FALSE,
    quote = FALSE
  )
  utils::write.table(
    batch_truth, paths$batch_truth, sep = "\t", row.names = FALSE,
    quote = FALSE
  )
  summary <- list(
    seed = seed,
    replicate_mode = if (is.null(design)) replicate_mode else "custom",
    annotation_gff3 = if (is.null(annotation_gff3)) NULL else {
      normalizePath(annotation_gff3, mustWork = TRUE)
    },
    counts = list(
      source_samples = length(source_samples),
      biological_replicates = length(biological_ids),
      libraries = nrow(metadata),
      genes = n_genes,
      conditions = length(conditions_used),
      batches = length(batches_used)
    ),
    parameters = list(
      condition_effect_fraction = condition_effect_fraction,
      library_size_mean = library_size_mean,
      dispersion_mean = dispersion_mean,
      n_latent = n_latent
    )
  )
  jsonlite::write_json(summary, paths$summary, pretty = TRUE, auto_unbox = TRUE)

  if (isTRUE(simulate_reads)) {
    if (is.null(transcript_fasta)) {
      stop("transcript_fasta is required when simulate_reads = TRUE")
    }
    manifest <- simulate_rnaseq_reads(
      expression = counts,
      transcript_fasta = transcript_fasta,
      out_dir = paste0(out_prefix, "_fastq"),
      transcript_gene_map = transcript_gene_map,
      library_size = read_library_size,
      read_length = read_length,
      paired_end = paired_end,
      seed = seed
    )
    paths$read_manifest <- attr(manifest, "manifest_path")
  }
  message(
    "Standalone RNA-seq simulation complete: ", length(biological_ids),
    " biological replicates, ", nrow(metadata), " libraries, and ",
    n_genes, " genes."
  )
  invisible(paths)
}

#' Simulate a standalone single-cell RNA-seq experiment
#'
#' Generate single-cell UMI counts without GWAS genotypes. Technical
#' replicates share a biological profile but receive independent cells and
#' technical effects; biological replicates receive independent biological
#' effects; mixed designs include both levels. Outputs are compatible with
#' [aggregate_scrnaseq_pseudobulk()] and [plot_scrnaseq_results()].
#'
#' @param out_prefix Prefix for sparse counts, metadata, truth, and summaries.
#' @param sample_names Names of underlying samples or experimental groups.
#' @param replicates Replicate count per sample. Interpretation follows
#'   `replicate_mode` as in [simulate_rnaseq_experiment()].
#' @param replicate_mode `"technical"`, `"biological"`, or `"mixed"`.
#' @param technical_replicates Technical captures per biological replicate in
#'   mixed mode.
#' @param conditions Condition label per sample, or one shared label.
#' @param batch_levels Batch labels assigned to capture libraries.
#' @param design Optional aggregated experiment design data frame or TSV.
#' @param annotation_gff3 Optional GFF3 supplying genes.
#' @param n_genes Number of genes when annotation is absent, or an optional
#'   subset size with annotation.
#' @param cells_per_replicate Cells generated for each capture library.
#' @param cell_type_proportions Named cell-type proportions.
#' @param trajectory_cell_types Cell types assigned pseudotime.
#' @param condition_effect_fraction Fraction of condition-responsive genes.
#' @param markers_per_cell_type Number of marker genes per cell type.
#' @param marker_effect_mean Mean marker effect in log2 units.
#' @param marker_effect_sd Marker-effect standard deviation.
#' @param pseudotime_gene_fraction Fraction of trajectory-responsive genes.
#' @param pseudotime_effect_sd Standard deviation of pseudotime effects.
#' @param baseline_mean_log2 Mean baseline log2 abundance.
#' @param baseline_sd_log2 Baseline standard deviation.
#' @param sample_effect_sd Source-sample effect standard deviation.
#' @param biological_replicate_sd Biological-replicate effect standard
#'   deviation, shared across technical captures.
#' @param condition_effect_sd Condition-effect standard deviation.
#' @param batch_effect_sd Batch-effect standard deviation.
#' @param technical_effect_sd Capture-library technical-effect standard
#'   deviation.
#' @param n_latent Number of cell-level latent factors.
#' @param latent_effect_sd Latent-loading standard deviation.
#' @param residual_sd Cell-level residual standard deviation.
#' @param library_size_mean Mean UMIs per cell.
#' @param library_size_sdlog Log-scale library-size standard deviation.
#' @param dispersion_mean Mean negative-binomial dispersion.
#' @param dispersion_shape Gamma shape controlling dispersion variation.
#' @param dropout_rate Additional expression-dependent dropout probability.
#' @param doublet_rate Fraction of doublet cells.
#' @param ambient_fraction Ambient-RNA fraction.
#' @param backend `"auto"`, `"splatter"`, or `"native"` technical scaffold.
#' @param write_dense_tsv Also write a dense gene-by-cell TSV.
#' @param write_sce Write a `SingleCellExperiment` when available.
#' @param min_cells_pseudobulk Minimum cells in pseudobulk libraries.
#' @param seed Random-number seed.
#'
#' @return Invisibly returns a named list of output paths. The `backend`
#'   attribute records the backend used.
#' @examples
#' \dontrun{
#' simulate_scrnaseq_experiment(
#'   out_prefix = "results/scrna/sample_A",
#'   sample_names = "sample_A",
#'   replicates = 3,
#'   replicate_mode = "technical",
#'   cells_per_replicate = 500,
#'   cell_type_proportions = c(T_cell = 0.4, B_cell = 0.3, Monocyte = 0.3),
#'   seed = 1
#' )
#' }
#' @export
simulate_scrnaseq_experiment <- function(
    out_prefix,
    sample_names = "sample_A",
    replicates = 3L,
    replicate_mode = c("technical", "biological", "mixed"),
    technical_replicates = 2L,
    conditions = "control",
    batch_levels = "batch1",
    design = NULL,
    annotation_gff3 = NULL,
    n_genes = NULL,
    cells_per_replicate = 500L,
    cell_type_proportions = c(
      T_cell = 0.4, B_cell = 0.3, Monocyte = 0.3
    ),
    trajectory_cell_types = NULL,
    condition_effect_fraction = 0.1,
    markers_per_cell_type = NULL,
    marker_effect_mean = 1.5,
    marker_effect_sd = 0.3,
    pseudotime_gene_fraction = 0.1,
    pseudotime_effect_sd = 0.5,
    baseline_mean_log2 = 5,
    baseline_sd_log2 = 1,
    sample_effect_sd = 0.25,
    biological_replicate_sd = 0.2,
    condition_effect_sd = 0.6,
    batch_effect_sd = 0.2,
    technical_effect_sd = 0.08,
    n_latent = 2L,
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
  replicate_mode <- match.arg(replicate_mode)
  backend <- match.arg(backend)
  for (name in c(
    "condition_effect_fraction", "pseudotime_gene_fraction",
    "dropout_rate", "doublet_rate", "ambient_fraction"
  )) {
    .simitall_scrna_validate_probability(get(name), name)
  }
  .simitall_validate_effect_scales(c(
    marker_effect_sd = marker_effect_sd,
    pseudotime_effect_sd = pseudotime_effect_sd,
    baseline_sd_log2 = baseline_sd_log2,
    sample_effect_sd = sample_effect_sd,
    biological_replicate_sd = biological_replicate_sd,
    condition_effect_sd = condition_effect_sd,
    batch_effect_sd = batch_effect_sd,
    technical_effect_sd = technical_effect_sd,
    latent_effect_sd = latent_effect_sd,
    residual_sd = residual_sd
  ))
  if (library_size_mean <= 0 || library_size_sdlog < 0 ||
      dispersion_mean <= 0 || dispersion_shape <= 0) {
    stop("Library-size and dispersion parameters are invalid")
  }
  if (!length(cell_type_proportions) || any(!is.finite(cell_type_proportions)) ||
      any(cell_type_proportions <= 0) || is.null(names(cell_type_proportions)) ||
      any(!nzchar(names(cell_type_proportions))) ||
      anyDuplicated(names(cell_type_proportions))) {
    stop("cell_type_proportions must be positive with unique non-empty names")
  }
  cell_type_proportions <- cell_type_proportions / sum(cell_type_proportions)
  if (is.null(trajectory_cell_types)) {
    trajectory_cell_types <- names(cell_type_proportions)[1L]
  }
  if (!all(trajectory_cell_types %in% names(cell_type_proportions))) {
    stop("trajectory_cell_types must occur in cell_type_proportions")
  }
  if (length(n_latent) != 1L || !is.finite(n_latent) || n_latent < 0L) {
    stop("n_latent must be a non-negative integer")
  }

  set.seed(seed)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  library_metadata <- .simitall_expand_experiment_design(
    sample_names = sample_names,
    replicates = replicates,
    technical_replicates = technical_replicates,
    replicate_mode = replicate_mode,
    conditions = conditions,
    batch_levels = batch_levels,
    design = design
  )
  genes <- .simitall_standalone_genes(annotation_gff3, n_genes)
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
  cells <- .simitall_make_experiment_cells(
    library_metadata,
    cells_per_replicate,
    cell_type_proportions,
    trajectory_cell_types
  )
  n_cells <- nrow(cells)
  biological_ids <- unique(library_metadata$biological_id)
  source_samples <- unique(library_metadata$sample)
  conditions_used <- unique(library_metadata$condition)
  batches_used <- unique(library_metadata$batch)

  baseline <- stats::rnorm(n_genes, baseline_mean_log2, baseline_sd_log2)
  library_sizes <- pmax(100L, round(stats::rlnorm(
    n_cells,
    meanlog = log(library_size_mean) - 0.5 * library_size_sdlog^2,
    sdlog = library_size_sdlog
  )))
  scaffold <- .simitall_scrna_scaffold(
    backend, n_genes, n_cells, baseline, library_sizes, seed + 101L
  )
  backend_used <- scaffold$backend
  expected <- scaffold$expected
  baseline <- log2(pmax(rowMeans(expected), 1e-8))
  delta <- matrix(0, nrow = n_genes, ncol = n_cells)

  sample_effects <- matrix(
    stats::rnorm(n_genes * length(source_samples), 0, sample_effect_sd),
    nrow = n_genes,
    dimnames = list(genes$gene_id, source_samples)
  )
  biological_effects <- matrix(
    stats::rnorm(
      n_genes * length(biological_ids), 0, biological_replicate_sd
    ),
    nrow = n_genes,
    dimnames = list(genes$gene_id, biological_ids)
  )
  technical_effects <- matrix(
    stats::rnorm(
      n_genes * nrow(library_metadata), 0, technical_effect_sd
    ),
    nrow = n_genes,
    dimnames = list(genes$gene_id, library_metadata$library_id)
  )
  delta <- delta + sample_effects[, cells$sample, drop = FALSE]
  delta <- delta + biological_effects[, cells$biological_id, drop = FALSE]
  delta <- delta + technical_effects[, cells$library_id, drop = FALSE]

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
    gene_id = character(), condition = character(), effect_log2 = numeric(),
    stringsAsFactors = FALSE
  )
  n_condition_genes <- min(
    n_genes, round(n_genes * condition_effect_fraction)
  )
  if (n_condition_genes > 0L && length(conditions_used) > 1L) {
    selected <- sample(seq_len(n_genes), n_condition_genes)
    for (condition in conditions_used[-1L]) {
      effects <- stats::rnorm(n_condition_genes, 0, condition_effect_sd)
      active <- cells$condition == condition
      delta[selected, active] <- sweep(
        delta[selected, active, drop = FALSE], 1L, effects, "+"
      )
      condition_truth <- rbind(
        condition_truth,
        data.frame(
          gene_id = genes$gene_id[selected], condition = condition,
          effect_log2 = effects, stringsAsFactors = FALSE
        )
      )
    }
  }
  batch_truth <- data.frame(
    gene_id = character(), batch = character(), effect_log2 = numeric(),
    stringsAsFactors = FALSE
  )
  if (length(batches_used) > 1L) {
    for (batch in batches_used[-1L]) {
      effects <- stats::rnorm(n_genes, 0, batch_effect_sd)
      active <- cells$batch == batch
      delta[, active] <- sweep(
        delta[, active, drop = FALSE], 1L, effects, "+"
      )
      batch_truth <- rbind(
        batch_truth,
        data.frame(
          gene_id = genes$gene_id, batch = batch,
          effect_log2 = effects, stringsAsFactors = FALSE
        )
      )
    }
  }

  pseudotime_truth <- data.frame(
    gene_id = character(), effect_log2 = numeric(),
    trajectory_cell_types = character(), stringsAsFactors = FALSE
  )
  n_pseudotime_genes <- min(
    n_genes, round(n_genes * pseudotime_gene_fraction)
  )
  trajectory <- which(!is.na(cells$pseudotime))
  if (n_pseudotime_genes > 0L && length(trajectory)) {
    selected <- sample(seq_len(n_genes), n_pseudotime_genes)
    effects <- stats::rnorm(n_pseudotime_genes, 0, pseudotime_effect_sd)
    delta[selected, trajectory] <- delta[selected, trajectory, drop = FALSE] +
      tcrossprod(effects, cells$pseudotime[trajectory] - 0.5)
    pseudotime_truth <- data.frame(
      gene_id = genes$gene_id[selected],
      effect_log2 = effects,
      trajectory_cell_types = paste(trajectory_cell_types, collapse = ","),
      stringsAsFactors = FALSE
    )
  }

  n_latent <- as.integer(n_latent)
  latent_scores <- matrix(
    numeric(), nrow = n_cells, ncol = 0L,
    dimnames = list(cells$cell_id, character())
  )
  latent_loadings <- matrix(
    numeric(), nrow = n_genes, ncol = 0L,
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
      stats::rnorm(n_genes * n_cells, 0, residual_sd), nrow = n_genes
    )
  }
  # Clamp extreme values without dropping matrix dimensions for one-cell runs.
  delta[delta < -20] <- -20
  delta[delta > 20] <- 20
  expected <- expected * 2^delta
  expected <- sweep(expected, 2L, pmax(colSums(expected), 1e-8), "/")
  expected <- sweep(expected, 2L, library_sizes, "*")
  if (ambient_fraction > 0) {
    ambient <- rowMeans(expected)
    ambient <- ambient / sum(ambient)
    expected <- (1 - ambient_fraction) * expected +
      ambient_fraction * tcrossprod(ambient, library_sizes)
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
  dispersion <- pmax(
    stats::rgamma(
      n_genes,
      shape = dispersion_shape,
      rate = dispersion_shape / dispersion_mean
    ),
    1e-8
  )
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
  cells$target_library_size <- library_sizes
  cells$observed_library_size <- colSums(counts)
  cells$detected_genes <- colSums(counts > 0)
  genes$baseline_log2 <- baseline
  genes$dispersion <- dispersion
  marker_output <- marker_truth[, setdiff(names(marker_truth), "gene_index")]
  eqtl_truth <- data.frame(
    gene_id = character(), variant_id = character(), effect_log2 = numeric(),
    scope = character(), cell_type = character(), stringsAsFactors = FALSE
  )

  paths <- list(
    cell_metadata = paste0(out_prefix, ".cell_metadata.tsv"),
    replicate_metadata = paste0(out_prefix, ".replicate_metadata.tsv"),
    gene_metadata = paste0(out_prefix, ".gene_metadata.tsv"),
    marker_truth = paste0(out_prefix, ".marker_truth.tsv"),
    condition_truth = paste0(out_prefix, ".condition_truth.tsv"),
    batch_truth = paste0(out_prefix, ".batch_truth.tsv"),
    pseudotime_truth = paste0(out_prefix, ".pseudotime_truth.tsv"),
    sample_effects = paste0(out_prefix, ".sample_effects.tsv"),
    biological_effects = paste0(out_prefix, ".biological_effects.tsv"),
    technical_effects = paste0(out_prefix, ".technical_effects.tsv"),
    latent_scores = paste0(out_prefix, ".latent_scores.tsv"),
    latent_loadings = paste0(out_prefix, ".latent_loadings.tsv"),
    object_rds = paste0(out_prefix, ".scrnaseq.rds"),
    summary = paste0(out_prefix, ".summary.json")
  )
  paths <- c(
    paths,
    .simitall_scrna_write_sparse(counts, genes, cells, out_prefix)
  )
  if (isTRUE(write_dense_tsv)) {
    paths$counts <- paste0(out_prefix, ".counts.tsv")
    .simitall_write_feature_matrix(counts, paths$counts)
  }
  for (item in list(
    list(cells, paths$cell_metadata),
    list(library_metadata, paths$replicate_metadata),
    list(genes, paths$gene_metadata),
    list(marker_output, paths$marker_truth),
    list(condition_truth, paths$condition_truth),
    list(batch_truth, paths$batch_truth),
    list(pseudotime_truth, paths$pseudotime_truth)
  )) {
    utils::write.table(
      item[[1L]], item[[2L]], sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
  .simitall_write_feature_matrix(sample_effects, paths$sample_effects)
  .simitall_write_feature_matrix(
    biological_effects, paths$biological_effects
  )
  .simitall_write_feature_matrix(
    technical_effects, paths$technical_effects
  )
  .simitall_write_feature_matrix(
    t(latent_scores), paths$latent_scores, "factor"
  )
  .simitall_write_feature_matrix(
    latent_loadings, paths$latent_loadings
  )

  sparse_counts <- Matrix::Matrix(counts, sparse = TRUE)
  donor_metadata <- library_metadata[
    match(biological_ids, library_metadata$biological_id),
    c("biological_id", "sample", "biological_replicate", "condition"),
    drop = FALSE
  ]
  names(donor_metadata)[1L] <- "donor"
  object <- structure(
    list(
      counts = sparse_counts,
      cell_metadata = cells,
      donor_metadata = donor_metadata,
      replicate_metadata = library_metadata,
      gene_metadata = genes,
      eqtl_truth = eqtl_truth,
      marker_truth = marker_output,
      condition_truth = condition_truth,
      batch_truth = batch_truth,
      pseudotime_truth = pseudotime_truth,
      backend = backend_used,
      experiment_type = "standalone"
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
        replicate_metadata = library_metadata,
        marker_truth = marker_output
      )
    )
    paths$sce_rds <- paste0(out_prefix, ".sce.rds")
    saveRDS(sce, paths$sce_rds, compress = "xz")
  } else if (isTRUE(write_sce)) {
    message(
      "SingleCellExperiment is not installed; portable RDS and sparse files ",
      "were still written."
    )
  }

  donor_pseudobulk <- aggregate_scrnaseq_pseudobulk(
    counts,
    cells,
    paste0(out_prefix, ".pseudobulk"),
    group_by = c("donor", "cell_type"),
    min_cells = min_cells_pseudobulk
  )
  library_pseudobulk <- aggregate_scrnaseq_pseudobulk(
    counts,
    cells,
    paste0(out_prefix, ".library_pseudobulk"),
    group_by = c("library_id", "cell_type"),
    min_cells = min_cells_pseudobulk
  )
  paths$pseudobulk_counts <- donor_pseudobulk$counts
  paths$pseudobulk_expression <- donor_pseudobulk$expression
  paths$pseudobulk_metadata <- donor_pseudobulk$metadata
  paths$library_pseudobulk_counts <- library_pseudobulk$counts
  paths$library_pseudobulk_expression <- library_pseudobulk$expression
  paths$library_pseudobulk_metadata <- library_pseudobulk$metadata

  summary <- list(
    seed = seed,
    backend = backend_used,
    replicate_mode = if (is.null(design)) replicate_mode else "custom",
    counts = list(
      source_samples = length(source_samples),
      biological_replicates = length(biological_ids),
      capture_libraries = nrow(library_metadata),
      cells = n_cells,
      genes = n_genes,
      cell_types = length(cell_type_proportions),
      doublets = sum(cells$is_doublet),
      zero_fraction = mean(counts == 0),
      technical_dropout_fraction = technical_dropout_fraction
    ),
    parameters = list(
      cells_per_replicate = cells_per_replicate,
      cell_type_proportions = as.list(cell_type_proportions),
      condition_effect_fraction = condition_effect_fraction,
      dropout_rate = dropout_rate,
      doublet_rate = doublet_rate,
      ambient_fraction = ambient_fraction
    )
  )
  jsonlite::write_json(summary, paths$summary, pretty = TRUE, auto_unbox = TRUE)
  attr(paths, "backend") <- backend_used
  message(
    "Standalone single-cell simulation complete: ", length(biological_ids),
    " biological replicates, ", nrow(library_metadata), " capture libraries, ",
    n_cells, " cells, and ", n_genes, " genes (", backend_used, " backend)."
  )
  invisible(paths)
}
