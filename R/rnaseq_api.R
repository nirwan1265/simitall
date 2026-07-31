#' Simulate RNA-seq from GWAS individuals
#'
#' Generate gene-level RNA-seq counts for the same samples and genotypes used
#' in a GWAS cohort. Expression is assembled from baseline abundance, cis- and
#' trans-eQTL effects, condition effects, genotype-by-condition interactions,
#' batch effects, latent confounders, and residual noise. The function writes
#' analysis matrices and explicit causal truth sets, including optional
#' allele-specific expression (ASE) counts.
#'
#' @param genotype_file Genotype TSV or VCF/VCF.GZ. The TSV written by
#'   [simulate_gwas_cohort()] is accepted directly.
#' @param out_prefix Prefix for all RNA-seq and truth output files.
#' @param annotation_gff3 Optional GFF3 containing gene coordinates. When
#'   omitted, synthetic genes are distributed across the variant coordinates.
#' @param sample_metadata Optional data frame or TSV with a `sample` column.
#'   Existing `condition`, `batch`, and population columns are preserved.
#' @param n_genes Number of genes to simulate. With a GFF3, `NULL` uses every
#'   gene; otherwise the default is 1000.
#' @param ploidy Number of homologous chromosome copies represented by genotype
#'   dosages. This is used to identify heterozygotes for ASE simulation.
#' @param n_cis_eqtl Number of genes assigned one causal cis-eQTL.
#' @param n_trans_eqtl Number of causal trans-eQTL pairs.
#' @param cis_window_bp Maximum distance in base pairs between a cis variant
#'   and a gene TSS.
#' @param condition_levels Character vector of condition labels. Used only when
#'   `sample_metadata` does not already contain `condition`.
#' @param batch_levels Character vector of batch labels. Used only when
#'   `sample_metadata` does not already contain `batch`.
#' @param condition_effect_fraction Fraction of genes with a nonzero condition
#'   effect.
#' @param gxe_fraction Fraction of eQTLs with a genotype-by-condition effect.
#' @param ase_fraction Fraction of cis-eQTL pairs used for ASE simulation. Set
#'   to zero to disable ASE output.
#' @param n_latent Number of unobserved sample-level factors.
#' @param baseline_mean_log2 Mean baseline log2 abundance.
#' @param baseline_sd_log2 Standard deviation of baseline log2 abundance.
#' @param cis_effect_sd Standard deviation of cis-eQTL effects in log2 units.
#' @param trans_effect_sd Standard deviation of trans-eQTL effects in log2
#'   units.
#' @param condition_effect_sd Standard deviation of condition effects in log2
#'   units.
#' @param gxe_effect_sd Standard deviation of genotype-by-condition effects in
#'   log2 units.
#' @param batch_effect_sd Standard deviation of batch effects in log2 units.
#' @param latent_effect_sd Standard deviation of latent-factor loadings.
#' @param residual_sd Standard deviation of gene-by-sample residual noise on
#'   the log2 abundance scale.
#' @param library_size_mean Mean gene-level count library size per sample.
#' @param library_size_sdlog Log-scale standard deviation of library sizes.
#' @param dispersion_mean Mean negative-binomial dispersion where
#'   `variance = mean + dispersion * mean^2`.
#' @param dispersion_shape Gamma shape used to vary dispersion among genes.
#' @param simulate_reads Also generate FASTQ reads with
#'   [simulate_rnaseq_reads()].
#' @param transcript_fasta Transcript FASTA required when
#'   `simulate_reads = TRUE`.
#' @param transcript_gene_map Optional transcript-to-gene mapping passed to
#'   [simulate_rnaseq_reads()].
#' @param read_library_size Number of reads or read pairs generated per sample.
#' @param read_length Simulated read length.
#' @param paired_end Generate paired-end reads.
#' @param seed Random-number seed.
#'
#' @return Invisibly returns a named list of output paths.
#' @examples
#' \dontrun{
#' simulate_rnaseq_from_gwas(
#'   genotype_file = "results/gwas/demo.geno.tsv",
#'   out_prefix = "results/rnaseq/demo",
#'   annotation_gff3 = "results/genome.gff3",
#'   n_cis_eqtl = 200,
#'   n_trans_eqtl = 50,
#'   condition_levels = c("control", "drought"),
#'   batch_levels = c("batch1", "batch2"),
#'   seed = 1
#' )
#' }
#' @export
simulate_rnaseq_from_gwas <- function(
    genotype_file,
    out_prefix,
    annotation_gff3 = NULL,
    sample_metadata = NULL,
    n_genes = NULL,
    ploidy = 2L,
    n_cis_eqtl = 100L,
    n_trans_eqtl = 25L,
    cis_window_bp = 1e6,
    condition_levels = c("control", "treatment"),
    batch_levels = c("batch1", "batch2"),
    condition_effect_fraction = 0.1,
    gxe_fraction = 0.2,
    ase_fraction = 0.25,
    n_latent = 3L,
    baseline_mean_log2 = 6,
    baseline_sd_log2 = 1,
    cis_effect_sd = 0.6,
    trans_effect_sd = 0.35,
    condition_effect_sd = 0.6,
    gxe_effect_sd = 0.4,
    batch_effect_sd = 0.2,
    latent_effect_sd = 0.25,
    residual_sd = 0.15,
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
  probability_parameters <- c(
    condition_effect_fraction = condition_effect_fraction,
    gxe_fraction = gxe_fraction,
    ase_fraction = ase_fraction
  )
  if (any(!is.finite(probability_parameters)) ||
      any(probability_parameters < 0 | probability_parameters > 1)) {
    stop("Effect fractions must be between zero and one")
  }
  if (n_cis_eqtl < 0L || n_trans_eqtl < 0L) {
    stop("n_cis_eqtl and n_trans_eqtl must be non-negative")
  }
  if (length(ploidy) != 1L || is.na(ploidy) || ploidy < 1L ||
      ploidy != as.integer(ploidy)) {
    stop("ploidy must be a positive integer")
  }
  ploidy <- as.integer(ploidy)
  if (!is.finite(cis_window_bp) || cis_window_bp < 0) {
    stop("cis_window_bp must be non-negative")
  }
  if (library_size_mean <= 0 || dispersion_mean <= 0 || dispersion_shape <= 0) {
    stop("Library size and dispersion parameters must be positive")
  }
  if (!length(condition_levels) || !length(batch_levels)) {
    stop("condition_levels and batch_levels cannot be empty")
  }

  set.seed(seed)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  genotype_data <- .simitall_read_genotypes(genotype_file)
  variants <- genotype_data$variants
  genotype_observed <- genotype_data$genotype
  sample_ids <- colnames(genotype_observed)
  if (length(sample_ids) < 4L) {
    stop("At least four GWAS individuals are required for RNA-seq simulation")
  }
  genotype <- .simitall_impute_genotypes(genotype_observed)
  if (max(genotype, na.rm = TRUE) > ploidy) {
    stop("Observed genotype dosage exceeds ploidy = ", ploidy)
  }

  metadata <- .simitall_read_sample_metadata(sample_metadata, sample_ids)
  if (!"condition" %in% names(metadata)) {
    metadata$condition <- sample(
      rep(condition_levels, length.out = length(sample_ids)),
      length(sample_ids),
      replace = FALSE
    )
  }
  if (!"batch" %in% names(metadata)) {
    metadata$batch <- NA_character_
    for (condition in unique(metadata$condition)) {
      index <- which(metadata$condition == condition)
      metadata$batch[index] <- sample(
        rep(batch_levels, length.out = length(index)),
        length(index),
        replace = FALSE
      )
    }
  }
  if (anyNA(metadata$condition) || anyNA(metadata$batch)) {
    stop("condition and batch labels cannot be missing")
  }
  observed_conditions <- unique(as.character(metadata$condition))
  observed_batches <- unique(as.character(metadata$batch))
  condition_levels_used <- c(
    condition_levels[condition_levels %in% observed_conditions],
    setdiff(observed_conditions, condition_levels)
  )
  batch_levels_used <- c(
    batch_levels[batch_levels %in% observed_batches],
    setdiff(observed_batches, batch_levels)
  )
  metadata$condition <- factor(metadata$condition, levels = condition_levels_used)
  metadata$batch <- factor(metadata$batch, levels = batch_levels_used)
  reference_condition <- condition_levels_used[1L]
  treatment_label <- if (length(condition_levels_used) > 1L) {
    condition_levels_used[2L]
  } else {
    reference_condition
  }
  treatment_indicator <- as.numeric(metadata$condition != reference_condition)
  if (length(condition_levels_used) < 2L) {
    gxe_fraction <- 0
  }

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
    message(
      "Mapped the single genotype contig to annotation contig '",
      unique(genes$seqname), "'."
    )
  }
  rownames(genes) <- genes$gene_id

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

  n_samples <- length(sample_ids)
  n_simulated_genes <- nrow(genes)
  baseline <- stats::rnorm(
    n_simulated_genes, baseline_mean_log2, baseline_sd_log2
  )
  eta <- matrix(
    baseline,
    nrow = n_simulated_genes,
    ncol = n_samples,
    dimnames = list(genes$gene_id, sample_ids)
  )
  centered_genotype <- sweep(genotype, 1L, rowMeans(genotype), "-")
  for (i in seq_len(nrow(eqtls))) {
    gene_index <- eqtls$gene_index[i]
    dosage <- centered_genotype[eqtls$variant_index[i], ]
    eta[gene_index, ] <- eta[gene_index, ] + eqtls$effect_log2[i] * dosage
    if (eqtls$gxe_effect_log2[i] != 0) {
      eta[gene_index, ] <- eta[gene_index, ] +
        eqtls$gxe_effect_log2[i] * dosage * treatment_indicator
    }
  }

  n_condition_genes <- min(
    n_simulated_genes,
    round(n_simulated_genes * condition_effect_fraction)
  )
  condition_truth <- data.frame(
    gene_id = character(), condition = character(), effect_log2 = numeric(),
    stringsAsFactors = FALSE
  )
  if (n_condition_genes > 0L && length(condition_levels_used) > 1L) {
    condition_genes <- sample(seq_len(n_simulated_genes), n_condition_genes)
    for (condition in condition_levels_used[-1L]) {
      effects <- stats::rnorm(n_condition_genes, 0, condition_effect_sd)
      sample_index <- metadata$condition == condition
      eta[condition_genes, sample_index] <- sweep(
        eta[condition_genes, sample_index, drop = FALSE], 1L, effects, "+"
      )
      condition_truth <- rbind(
        condition_truth,
        data.frame(
          gene_id = genes$gene_id[condition_genes],
          condition = condition,
          effect_log2 = effects,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  batch_truth <- data.frame(
    gene_id = character(), batch = character(), effect_log2 = numeric(),
    stringsAsFactors = FALSE
  )
  if (length(batch_levels_used) > 1L) {
    for (batch in batch_levels_used[-1L]) {
      effects <- stats::rnorm(n_simulated_genes, 0, batch_effect_sd)
      sample_index <- metadata$batch == batch
      eta[, sample_index] <- sweep(
        eta[, sample_index, drop = FALSE], 1L, effects, "+"
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

  n_latent <- max(0L, as.integer(n_latent))
  latent_scores <- matrix(
    numeric(), nrow = n_samples, ncol = 0L,
    dimnames = list(sample_ids, character())
  )
  latent_loadings <- matrix(
    numeric(), nrow = n_simulated_genes, ncol = 0L,
    dimnames = list(genes$gene_id, character())
  )
  if (n_latent > 0L) {
    latent_names <- paste0("latent_", seq_len(n_latent))
    latent_scores <- matrix(
      stats::rnorm(n_samples * n_latent),
      nrow = n_samples,
      dimnames = list(sample_ids, latent_names)
    )
    latent_loadings <- matrix(
      stats::rnorm(n_simulated_genes * n_latent, 0, latent_effect_sd),
      nrow = n_simulated_genes,
      dimnames = list(genes$gene_id, latent_names)
    )
    eta <- eta + latent_loadings %*% t(latent_scores)
  }
  eta <- eta + matrix(
    stats::rnorm(n_simulated_genes * n_samples, 0, residual_sd),
    nrow = n_simulated_genes
  )
  eta[eta < -20] <- -20
  eta[eta > 30] <- 30

  abundance <- 2^eta
  relative_abundance <- sweep(abundance, 2L, colSums(abundance), "/")
  library_sizes <- round(stats::rlnorm(
    n_samples,
    meanlog = log(library_size_mean) - 0.5 * library_size_sdlog^2,
    sdlog = library_size_sdlog
  ))
  expected_counts <- sweep(relative_abundance, 2L, library_sizes, "*")
  dispersion <- stats::rgamma(
    n_simulated_genes,
    shape = dispersion_shape,
    rate = dispersion_shape / dispersion_mean
  )
  dispersion <- pmax(dispersion, 1e-8)
  counts <- matrix(
    stats::rnbinom(
      length(expected_counts),
      mu = as.vector(expected_counts),
      size = rep(1 / dispersion, times = n_samples)
    ),
    nrow = n_simulated_genes,
    dimnames = dimnames(expected_counts)
  )
  observed_library_sizes <- pmax(1, colSums(counts))
  log2_cpm <- log2(sweep(counts, 2L, observed_library_sizes, "/") * 1e6 + 1)

  metadata$condition <- as.character(metadata$condition)
  metadata$batch <- as.character(metadata$batch)
  metadata$condition_code <- match(metadata$condition, condition_levels_used) - 1L
  metadata$batch_code <- match(metadata$batch, batch_levels_used) - 1L
  metadata$target_library_size <- library_sizes
  metadata$observed_library_size <- observed_library_sizes
  genes$baseline_log2 <- baseline
  genes$dispersion <- dispersion

  paths <- list(
    counts = paste0(out_prefix, ".counts.tsv"),
    expression = paste0(out_prefix, ".expression.tsv"),
    expected_counts = paste0(out_prefix, ".expected_counts.tsv"),
    sample_metadata = paste0(out_prefix, ".sample_metadata.tsv"),
    gene_metadata = paste0(out_prefix, ".gene_metadata.tsv"),
    eqtl_truth = paste0(out_prefix, ".eqtl_truth.tsv"),
    condition_truth = paste0(out_prefix, ".condition_truth.tsv"),
    batch_truth = paste0(out_prefix, ".batch_truth.tsv"),
    latent_scores = paste0(out_prefix, ".latent_scores.tsv"),
    latent_loadings = paste0(out_prefix, ".latent_loadings.tsv"),
    ase_truth = paste0(out_prefix, ".ase_truth.tsv"),
    ase_counts = paste0(out_prefix, ".ase_counts.tsv"),
    summary = paste0(out_prefix, ".summary.json")
  )
  .simitall_write_feature_matrix(counts, paths$counts)
  .simitall_write_feature_matrix(log2_cpm, paths$expression)
  .simitall_write_feature_matrix(expected_counts, paths$expected_counts)
  utils::write.table(
    metadata, paths$sample_metadata, sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    genes, paths$gene_metadata, sep = "\t", row.names = FALSE, quote = FALSE
  )
  eqtl_output <- eqtls[, setdiff(names(eqtls), c("gene_index", "variant_index"))]
  utils::write.table(
    eqtl_output, paths$eqtl_truth, sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    condition_truth, paths$condition_truth,
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    batch_truth, paths$batch_truth,
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  .simitall_write_feature_matrix(t(latent_scores), paths$latent_scores, "factor")
  .simitall_write_feature_matrix(latent_loadings, paths$latent_loadings)

  cis_eqtls <- eqtls[eqtls$type == "cis", , drop = FALSE]
  n_ase <- min(nrow(cis_eqtls), round(nrow(cis_eqtls) * ase_fraction))
  ase_truth <- data.frame(
    gene_id = character(), variant_id = character(), alt_probability = numeric(),
    allelic_log2_odds = numeric(), stringsAsFactors = FALSE
  )
  ase_counts <- data.frame(
    sample = character(), gene_id = character(), variant_id = character(),
    total_count = integer(), ref_count = integer(), alt_count = integer(),
    alt_fraction = numeric(), stringsAsFactors = FALSE
  )
  if (n_ase > 0L) {
    selected_ase <- cis_eqtls[sample(seq_len(nrow(cis_eqtls)), n_ase), , drop = FALSE]
    direction <- sample(c(-1, 1), n_ase, replace = TRUE)
    allelic_log2_odds <- direction * pmax(0.25, abs(selected_ase$effect_log2))
    alt_probability <- stats::plogis(allelic_log2_odds * log(2))
    ase_truth <- data.frame(
      gene_id = selected_ase$gene_id,
      variant_id = selected_ase$variant_id,
      alt_probability = alt_probability,
      allelic_log2_odds = allelic_log2_odds,
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(selected_ase))) {
      dosage <- genotype_observed[selected_ase$variant_index[i], ]
      heterozygous <- which(dosage > 0 & dosage < ploidy & !is.na(dosage))
      if (!length(heterozygous)) next
      gene_counts <- counts[selected_ase$gene_index[i], heterozygous]
      allele_depth <- stats::rbinom(length(heterozygous), gene_counts, 0.25)
      alt_count <- stats::rbinom(
        length(heterozygous), allele_depth, alt_probability[i]
      )
      ase_counts <- rbind(
        ase_counts,
        data.frame(
          sample = sample_ids[heterozygous],
          gene_id = selected_ase$gene_id[i],
          variant_id = selected_ase$variant_id[i],
          total_count = allele_depth,
          ref_count = allele_depth - alt_count,
          alt_count = alt_count,
          alt_fraction = ifelse(allele_depth > 0, alt_count / allele_depth, NA_real_),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  utils::write.table(
    ase_truth, paths$ase_truth, sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    ase_counts, paths$ase_counts, sep = "\t", row.names = FALSE, quote = FALSE
  )

  summary <- list(
    seed = seed,
    genotype_file = normalizePath(genotype_file, mustWork = TRUE),
    annotation_gff3 = if (is.null(annotation_gff3)) NULL else {
      normalizePath(annotation_gff3, mustWork = TRUE)
    },
    counts = list(
      samples = n_samples,
      genes = n_simulated_genes,
      variants = nrow(variants),
      cis_eqtls = sum(eqtls$type == "cis"),
      trans_eqtls = sum(eqtls$type == "trans"),
      gxe_eqtls = sum(eqtls$gxe_effect_log2 != 0),
      condition_genes = length(unique(condition_truth$gene_id)),
      ase_pairs = nrow(ase_truth),
      ase_observations = nrow(ase_counts),
      latent_factors = n_latent
    ),
    parameters = list(
      cis_window_bp = cis_window_bp,
      conditions = condition_levels_used,
      batches = batch_levels_used,
      baseline_mean_log2 = baseline_mean_log2,
      baseline_sd_log2 = baseline_sd_log2,
      cis_effect_sd = cis_effect_sd,
      trans_effect_sd = trans_effect_sd,
      condition_effect_sd = condition_effect_sd,
      gxe_effect_sd = gxe_effect_sd,
      batch_effect_sd = batch_effect_sd,
      latent_effect_sd = latent_effect_sd,
      residual_sd = residual_sd,
      library_size_mean = library_size_mean,
      library_size_sdlog = library_size_sdlog,
      dispersion_mean = dispersion_mean,
      dispersion_shape = dispersion_shape,
      ploidy = ploidy
    )
  )
  jsonlite::write_json(summary, paths$summary, pretty = TRUE, auto_unbox = TRUE)

  if (isTRUE(simulate_reads)) {
    if (is.null(transcript_fasta)) {
      stop("transcript_fasta is required when simulate_reads = TRUE")
    }
    read_dir <- paste0(out_prefix, "_fastq")
    read_manifest <- simulate_rnaseq_reads(
      expression = counts,
      transcript_fasta = transcript_fasta,
      out_dir = read_dir,
      transcript_gene_map = transcript_gene_map,
      library_size = read_library_size,
      read_length = read_length,
      paired_end = paired_end,
      seed = seed
    )
    paths$read_manifest <- attr(read_manifest, "manifest_path")
    paths$read_directory <- read_dir
  }

  message(
    "RNA-seq simulation complete: ", n_samples, " samples, ",
    n_simulated_genes, " genes, and ", nrow(eqtls), " causal eQTLs."
  )
  invisible(paths)
}

#' Simulate RNA-seq FASTQ reads
#'
#' Convert a gene- or transcript-expression matrix into sample-level FASTQ
#' files with `Rsubread::simReads()`. A transcript-to-gene table can distribute
#' gene-level abundance equally among annotated isoforms.
#'
#' @param expression Expression matrix, data frame, or TSV path. Features are
#'   rows and samples are columns.
#' @param transcript_fasta Transcript sequences in FASTA or gzipped FASTA.
#' @param out_dir Output directory for FASTQ files and the read manifest.
#' @param transcript_gene_map Optional two-column data frame or TSV containing
#'   `transcript_id` and `gene_id`. Without it, transcript and feature IDs must
#'   match exactly.
#' @param library_size Number of reads or read pairs per sample. Supply one
#'   value or one value per sample.
#' @param read_length Read length, at most 250 bp. Built-in Rsubread error
#'   profiles are available for 75- and 100-bp reads.
#' @param paired_end Generate paired-end reads.
#' @param fragment_length_min Minimum RNA fragment length.
#' @param fragment_length_max Maximum RNA fragment length.
#' @param fragment_length_mean Mean RNA fragment length.
#' @param fragment_length_sd Standard deviation of RNA fragment lengths.
#' @param simulate_sequencing_error Simulate substitution sequencing errors.
#' @param quality_reference Optional file of Phred+33 quality strings, required
#'   for error simulation at read lengths other than 75 or 100 bp.
#' @param truth_in_read_names Include true transcript locations in read names.
#' @param seed Random-number seed.
#'
#' @return Invisibly returns the read manifest data frame. Its
#'   `manifest_path` attribute gives the written TSV path.
#' @examples
#' \dontrun{
#' simulate_rnaseq_reads(
#'   expression = "results/rnaseq/demo.counts.tsv",
#'   transcript_fasta = "results/genome.cds.fa",
#'   out_dir = "results/rnaseq/fastq",
#'   library_size = 100000,
#'   paired_end = TRUE
#' )
#' }
#' @export
simulate_rnaseq_reads <- function(
    expression,
    transcript_fasta,
    out_dir,
    transcript_gene_map = NULL,
    library_size = 1e5,
    read_length = 100L,
    paired_end = TRUE,
    fragment_length_min = max(100L, read_length),
    fragment_length_max = 500L,
    fragment_length_mean = 180,
    fragment_length_sd = 40,
    simulate_sequencing_error = TRUE,
    quality_reference = NULL,
    truth_in_read_names = TRUE,
    seed = 1L) {
  if (!requireNamespace("Rsubread", quietly = TRUE)) {
    stop(
      "Rsubread is required for FASTQ simulation. Install it with ",
      "BiocManager::install('Rsubread')."
    )
  }
  if (!file.exists(transcript_fasta)) {
    stop("Transcript FASTA does not exist: ", transcript_fasta)
  }
  if (read_length < 1L || read_length > 250L) {
    stop("read_length must be between 1 and 250")
  }
  if (simulate_sequencing_error &&
      !read_length %in% c(75L, 100L) && is.null(quality_reference)) {
    stop(
      "quality_reference is required for sequencing-error simulation unless ",
      "read_length is 75 or 100"
    )
  }

  expression <- .simitall_read_expression_matrix(expression)
  transcripts <- Rsubread::scanFasta(transcript_fasta, quiet = TRUE)
  if (is.null(transcript_gene_map)) {
    mapping <- data.frame(
      transcript_id = transcripts$TranscriptID,
      gene_id = transcripts$TranscriptID,
      stringsAsFactors = FALSE
    )
  } else {
    if (is.character(transcript_gene_map) && length(transcript_gene_map) == 1L) {
      transcript_gene_map <- utils::read.delim(
        transcript_gene_map,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
    mapping <- as.data.frame(transcript_gene_map, stringsAsFactors = FALSE)
    required <- c("transcript_id", "gene_id")
    if (!all(required %in% names(mapping))) {
      stop("transcript_gene_map needs transcript_id and gene_id columns")
    }
    mapping <- mapping[match(transcripts$TranscriptID, mapping$transcript_id), , drop = FALSE]
    mapping$transcript_id <- transcripts$TranscriptID
  }
  feature_index <- match(mapping$gene_id, rownames(expression))
  if (!any(!is.na(feature_index))) {
    stop("No transcript IDs or mapped gene IDs match the expression features")
  }
  isoform_count <- table(mapping$gene_id)

  n_samples <- ncol(expression)
  if (length(library_size) == 1L) {
    library_size <- rep(library_size, n_samples)
  }
  if (length(library_size) != n_samples || any(library_size <= 0)) {
    stop("library_size must be positive and have length one or the sample count")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  safe_names <- make.unique(gsub("[^A-Za-z0-9_.-]", "_", colnames(expression)))
  manifest_rows <- vector("list", n_samples)

  for (i in seq_len(n_samples)) {
    levels <- numeric(nrow(transcripts))
    matched <- which(!is.na(feature_index))
    levels[matched] <- expression[feature_index[matched], i] /
      as.numeric(isoform_count[mapping$gene_id[matched]])
    levels[!is.finite(levels) | levels < 0] <- 0
    if (sum(levels) <= 0) {
      stop("Sample has no positive expression values: ", colnames(expression)[i])
    }

    output_prefix <- file.path(out_dir, safe_names[i])
    set.seed(seed + i - 1L)
    read_truth <- Rsubread::simReads(
      transcript.file = transcript_fasta,
      expression.levels = levels,
      output.prefix = output_prefix,
      library.size = round(library_size[i]),
      read.length = as.integer(read_length),
      truth.in.read.names = truth_in_read_names,
      simulate.sequencing.error = simulate_sequencing_error,
      quality.reference = quality_reference,
      paired.end = paired_end,
      fragment.length.min = as.integer(fragment_length_min),
      fragment.length.max = as.integer(fragment_length_max),
      fragment.length.mean = fragment_length_mean,
      fragment.length.sd = fragment_length_sd
    )
    truth_path <- paste0(output_prefix, ".read_truth.tsv")
    utils::write.table(
      read_truth, truth_path, sep = "\t", row.names = FALSE, quote = FALSE
    )
    manifest_rows[[i]] <- data.frame(
      sample = colnames(expression)[i],
      read1 = paste0(output_prefix, "_R1.fastq.gz"),
      read2 = if (paired_end) paste0(output_prefix, "_R2.fastq.gz") else NA_character_,
      read_truth = truth_path,
      requested_library_size = round(library_size[i]),
      stringsAsFactors = FALSE
    )
  }

  manifest <- do.call(rbind, manifest_rows)
  manifest_path <- file.path(out_dir, "read_manifest.tsv")
  utils::write.table(
    manifest, manifest_path, sep = "\t", row.names = FALSE, quote = FALSE
  )
  attr(manifest, "manifest_path") <- manifest_path
  invisible(manifest)
}

#' Benchmark eQTL recovery
#'
#' Test genotype-expression associations and compare discoveries with a known
#' eQTL truth set. Candidate pairs can be restricted to neighborhoods around
#' simulated causal genes, all cis pairs, or all variant-gene pairs. Batch and
#' condition covariates are included when available, and a two-condition model
#' can also test genotype-by-condition interactions.
#'
#' @param genotype_file Genotype TSV or VCF/VCF.GZ.
#' @param expression_file Expression TSV or matrix with genes in rows and
#'   samples in columns. The `.expression.tsv` output from
#'   [simulate_rnaseq_from_gwas()] is ready to use.
#' @param out_prefix Prefix for association results and benchmark metrics.
#' @param sample_metadata Optional data frame or TSV with sample covariates.
#' @param truth_file Optional eQTL truth TSV written by
#'   [simulate_rnaseq_from_gwas()].
#' @param gene_metadata Optional data frame or TSV with `gene_id`, `seqname`,
#'   and `tss`, required for `pair_mode = "cis"`.
#' @param pair_mode Candidate-pair strategy: `"truth_neighborhood"` tests cis
#'   neighborhoods for truth genes plus causal trans pairs, `"cis"` tests all
#'   local pairs, and `"all"` tests every pair.
#' @param cis_window_bp Cis testing window around each TSS.
#' @param covariates Metadata columns to include in each regression.
#' @param test_gxe Test genotype-by-condition interaction when exactly two
#'   conditions are present.
#' @param expression_scale `"log2"`, `"counts"`, or `"auto"`. Count data are
#'   converted to log2 counts per million.
#' @param fdr_threshold Adjusted-p-value threshold used for benchmark metrics.
#' @param max_tests Safety limit on the number of variant-gene tests.
#'
#' @return Invisibly returns a list with result and metric paths and tables.
#' @examples
#' \dontrun{
#' benchmark_eqtl(
#'   genotype_file = "results/gwas/demo.geno.tsv",
#'   expression_file = "results/rnaseq/demo.expression.tsv",
#'   sample_metadata = "results/rnaseq/demo.sample_metadata.tsv",
#'   truth_file = "results/rnaseq/demo.eqtl_truth.tsv",
#'   out_prefix = "results/eqtl/demo"
#' )
#' }
#' @export
benchmark_eqtl <- function(
    genotype_file,
    expression_file,
    out_prefix,
    sample_metadata = NULL,
    truth_file = NULL,
    gene_metadata = NULL,
    pair_mode = c("truth_neighborhood", "cis", "all"),
    cis_window_bp = 1e6,
    covariates = c("condition", "batch"),
    test_gxe = TRUE,
    expression_scale = c("auto", "log2", "counts"),
    fdr_threshold = 0.05,
    max_tests = 1e6) {
  pair_mode <- match.arg(pair_mode)
  expression_scale <- match.arg(expression_scale)
  genotype_data <- .simitall_read_genotypes(genotype_file)
  variants <- genotype_data$variants
  genotype <- genotype_data$genotype
  expression <- .simitall_read_expression_matrix(expression_file)
  sample_ids <- intersect(colnames(genotype), colnames(expression))
  if (length(sample_ids) < 4L) {
    stop("At least four shared genotype and expression samples are required")
  }
  genotype <- genotype[, sample_ids, drop = FALSE]
  expression <- expression[, sample_ids, drop = FALSE]
  metadata <- .simitall_read_sample_metadata(sample_metadata, sample_ids)

  if (expression_scale == "auto") {
    finite_values <- expression[is.finite(expression)]
    looks_like_counts <- length(finite_values) &&
      all(finite_values >= 0) && all(abs(finite_values - round(finite_values)) < 1e-8) &&
      max(finite_values) > 50
    expression_scale <- if (looks_like_counts) "counts" else "log2"
  }
  if (expression_scale == "counts") {
    library_sizes <- pmax(1, colSums(expression, na.rm = TRUE))
    expression <- log2(sweep(expression, 2L, library_sizes, "/") * 1e6 + 1)
  }

  read_table <- function(value, label) {
    if (is.null(value)) return(NULL)
    if (is.character(value) && length(value) == 1L) {
      if (!file.exists(value)) stop(label, " does not exist: ", value)
      value <- utils::read.delim(
        value, header = TRUE, sep = "\t", check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
    as.data.frame(value, stringsAsFactors = FALSE)
  }
  truth <- read_table(truth_file, "Truth file")
  gene_table <- read_table(gene_metadata, "Gene metadata")

  if (pair_mode == "truth_neighborhood" && is.null(truth)) {
    stop("truth_file is required for pair_mode = 'truth_neighborhood'")
  }
  if (pair_mode == "cis" && is.null(gene_table)) {
    stop("gene_metadata is required for pair_mode = 'cis'")
  }

  make_cis_pairs <- function(gene_rows) {
    rows <- vector("list", nrow(gene_rows))
    for (i in seq_len(nrow(gene_rows))) {
      variant_index <- which(
        variants$seqname == gene_rows$seqname[i] &
          abs(variants$pos - gene_rows$tss[i]) <= cis_window_bp
      )
      rows[[i]] <- if (length(variant_index)) {
        data.frame(
          gene_id = gene_rows$gene_id[i],
          variant_id = variants$id[variant_index],
          stringsAsFactors = FALSE
        )
      } else NULL
    }
    rows <- rows[!vapply(rows, is.null, logical(1L))]
    if (length(rows)) do.call(rbind, rows) else data.frame(
      gene_id = character(), variant_id = character(), stringsAsFactors = FALSE
    )
  }

  if (pair_mode == "truth_neighborhood") {
    required_truth <- c("gene_id", "variant_id", "seqname", "gene_tss")
    if (!all(required_truth %in% names(truth))) {
      stop("truth_file is missing required columns: ", paste(required_truth, collapse = ", "))
    }
    truth_genes <- unique(truth[, c("gene_id", "seqname", "gene_tss")])
    names(truth_genes)[names(truth_genes) == "gene_tss"] <- "tss"
    pairs <- make_cis_pairs(truth_genes)
    pairs <- unique(rbind(pairs, truth[, c("gene_id", "variant_id")]))
  } else if (pair_mode == "cis") {
    required_genes <- c("gene_id", "seqname", "tss")
    if (!all(required_genes %in% names(gene_table))) {
      stop("gene_metadata is missing required columns: ", paste(required_genes, collapse = ", "))
    }
    gene_table <- gene_table[gene_table$gene_id %in% rownames(expression), , drop = FALSE]
    pairs <- make_cis_pairs(gene_table)
  } else {
    pairs <- expand.grid(
      gene_id = rownames(expression),
      variant_id = variants$id,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }
  pairs <- pairs[
    pairs$gene_id %in% rownames(expression) &
      pairs$variant_id %in% variants$id,
    , drop = FALSE
  ]
  pairs <- unique(pairs)
  if (!nrow(pairs)) {
    stop("No testable variant-gene pairs were found")
  }
  if (nrow(pairs) > max_tests) {
    stop(
      "Candidate set contains ", nrow(pairs), " tests, exceeding max_tests = ",
      max_tests, ". Reduce the cis window or raise max_tests explicitly."
    )
  }

  usable_covariates <- intersect(covariates, names(metadata))
  if (length(usable_covariates)) {
    covariate_data <- metadata[, usable_covariates, drop = FALSE]
    for (name in names(covariate_data)) {
      if (is.character(covariate_data[[name]])) {
        preferred_code <- paste0(name, "_code")
        if (preferred_code %in% names(metadata)) {
          level_order <- unique(as.character(metadata[[name]])[
            order(metadata[[preferred_code]])
          ])
          covariate_data[[name]] <- factor(
            covariate_data[[name]], levels = level_order
          )
        } else {
          covariate_data[[name]] <- factor(covariate_data[[name]])
        }
      }
    }
    covariate_matrix <- stats::model.matrix(~ ., data = covariate_data)
  } else {
    covariate_matrix <- matrix(1, nrow = length(sample_ids), ncol = 1L)
    colnames(covariate_matrix) <- "(Intercept)"
  }
  interaction <- NULL
  if (isTRUE(test_gxe) && "condition" %in% names(metadata) &&
      length(unique(metadata$condition)) == 2L) {
    if ("condition_code" %in% names(metadata)) {
      interaction <- as.numeric(metadata$condition_code > min(metadata$condition_code))
    } else {
      condition_values <- unique(as.character(metadata$condition))
      interaction <- as.numeric(metadata$condition == condition_values[2L])
    }
  }

  gene_index <- match(pairs$gene_id, rownames(expression))
  variant_index <- match(pairs$variant_id, variants$id)
  fits <- matrix(NA_real_, nrow = nrow(pairs), ncol = 5L)
  colnames(fits) <- c("beta", "se", "p_value", "gxe_beta", "gxe_p_value")
  for (i in seq_len(nrow(pairs))) {
    fits[i, ] <- .simitall_fit_eqtl(
      y = expression[gene_index[i], ],
      genotype = genotype[variant_index[i], ],
      covariate_matrix = covariate_matrix,
      interaction = interaction
    )
  }
  results <- cbind(pairs, as.data.frame(fits, stringsAsFactors = FALSE))
  results$q_value <- stats::p.adjust(results$p_value, method = "BH")
  results$gxe_q_value <- stats::p.adjust(results$gxe_p_value, method = "BH")
  results$variant_seqname <- variants$seqname[variant_index]
  results$variant_pos <- variants$pos[variant_index]

  truth_keys <- character()
  if (!is.null(truth)) {
    truth_keys <- unique(paste(truth$gene_id, truth$variant_id, sep = "::"))
  }
  result_keys <- paste(results$gene_id, results$variant_id, sep = "::")
  results$is_causal <- result_keys %in% truth_keys
  discovered <- !is.na(results$q_value) & results$q_value <= fdr_threshold
  true_positive <- sum(discovered & results$is_causal)
  false_positive <- sum(discovered & !results$is_causal)
  truth_in_tests <- sum(results$is_causal)
  metrics <- data.frame(
    tests = nrow(results),
    truth_pairs_in_tests = truth_in_tests,
    discoveries = sum(discovered),
    true_positives = true_positive,
    false_positives = false_positive,
    false_negatives = truth_in_tests - true_positive,
    precision = if (sum(discovered)) true_positive / sum(discovered) else NA_real_,
    recall = if (truth_in_tests) true_positive / truth_in_tests else NA_real_,
    empirical_fdr = if (sum(discovered)) false_positive / sum(discovered) else NA_real_,
    fdr_threshold = fdr_threshold,
    stringsAsFactors = FALSE
  )

  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  result_path <- paste0(out_prefix, ".eqtl_results.tsv")
  metric_path <- paste0(out_prefix, ".eqtl_metrics.tsv")
  utils::write.table(
    results, result_path, sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    metrics, metric_path, sep = "\t", row.names = FALSE, quote = FALSE
  )
  message(
    "eQTL benchmark complete: ", nrow(results), " tests and ",
    sum(discovered), " discoveries at FDR <= ", fdr_threshold, "."
  )
  invisible(list(
    results_path = result_path,
    metrics_path = metric_path,
    results = results,
    metrics = metrics
  ))
}
