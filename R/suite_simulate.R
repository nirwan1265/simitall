#' Simulate data for a Biodes design package
#'
#' Reads a Biodes design folder (`study_spec.json`, `field_layout.csv`,
#' `sample_metadata_template.tsv`) and simulates the primary outcome for every
#' experimental unit from the spec's `analysis.model`: treatment effects of the
#' size the study was powered for, a variance for each blocking factor and
#' random term, subsamples nested in their unit, and attrition. For count
#' assays (RNA-seq) it also simulates a gene x sample count matrix with
#' [simulate_rnaseq_experiment()], using the treatment as condition and the
#' block as batch. The true values used are written to `sim_truth.json`, so
#' [bioana_fit_biodes()] can check that it recovers them.
#'
#' Values the spec does not give come from `inst/suite/simulation_defaults.json`
#' and are listed under `assumptions` in `sim_truth.json`.
#'
#' @param biodes_dir Folder written by Biodes `write_design_package()`.
#' @param outdir Where to write the outputs (default: the same folder).
#' @param seed Random seed. Default: the design's `randomization_seed`.
#' @param effect Largest treatment difference. Default:
#'   `power.minimum_detectable_difference` from the spec.
#' @param residual_sd Residual SD of the unit means. Default:
#'   `power.residual_sd` from the spec.
#' @param n_genes Genes in the count matrix (count assays only).
#' @return Invisibly, the paths written.
#' @examples
#' \dontrun{
#' simitall_from_biodes("studyforge_run")
#' bioana_fit_biodes("studyforge_run")
#' }
#' @export
simitall_from_biodes <- function(biodes_dir, outdir = biodes_dir, seed = NULL, effect = NULL,
                                 residual_sd = NULL, n_genes = NULL) {
  h <- .suite_read_biodes(biodes_dir)
  d <- .suite_rules("simulation_defaults")
  rules <- .suite_rules("handoff")
  spec <- h$spec
  p <- h$parsed
  assumptions <- character()
  pick <- function(arg, spec_value, default, what) {
    if (!is.null(arg)) return(arg)
    if (!is.null(spec_value)) return(spec_value)
    assumptions <<- c(assumptions, sprintf("%s = %s (default: %s)", what, format(default$value), default$note))
    default$value
  }
  seed <- seed %||% spec$design$randomization_seed %||% d$seed
  set.seed(seed)
  rsd <- as.numeric(pick(residual_sd, spec$power$residual_sd, d$residual_sd, "residual SD"))
  eff_default <- list(value = d$effect_in_residual_sd$value * rsd, note = d$effect_in_residual_sd$note)
  eff <- as.numeric(pick(effect, spec$power$minimum_detectable_difference, eff_default, "treatment effect"))
  mu <- as.numeric(pick(NULL, spec$primary_outcome$baseline_mean, d$intercept, "intercept"))
  att <- as.numeric(pick(NULL, spec$power$attrition, d$attrition, "attrition"))
  n_sub <- as.integer(spec$subsamples_per_unit %||% 1L)
  if (is.na(n_sub) || n_sub < 1L) n_sub <- 1L

  meta <- h$meta
  if (nrow(meta) != nrow(h$layout)) stop("sample_metadata_template.tsv and field_layout.csv have different row counts.", call. = FALSE)
  frame <- .suite_frame(meta, p)
  n <- nrow(frame)
  y <- rep(mu, n)

  # treatment main effects: levels spread evenly from 0 to the full effect
  treat_vars <- intersect(names(p$role)[p$role == "treatment"], p$fixed_labels)
  treatment_effects <- list()
  for (v in treat_vars) {
    lev <- levels(frame[[v]])
    e <- stats::setNames(if (length(lev) > 1) seq(0, eff, length.out = length(lev)) else 0, lev)
    treatment_effects[[v]] <- as.list(e)
    y <- y + e[as.character(frame[[v]])]
  }
  interactions <- p$fixed_labels[grepl(":", p$fixed_labels)]
  # fixed nuisance terms (block, rep): random draws, fitted as fixed effects
  block_vars <- intersect(names(p$role)[p$role == "block"], p$fixed_labels)
  block_effects <- list()
  for (v in block_vars) {
    lev <- levels(frame[[v]])
    e <- stats::setNames(stats::rnorm(length(lev), 0, d$block_sd_ratio$value * rsd), lev)
    block_effects[[v]] <- as.list(e)
    y <- y + e[as.character(frame[[v]])]
  }
  random_sd <- list()
  for (r in p$random) {
    g <- interaction(frame[r$vars], drop = TRUE)
    s <- d$random_sd_ratio$value * rsd
    y <- y + stats::rnorm(nlevels(g), 0, s)[as.integer(g)]
    random_sd[[r$label]] <- s
  }
  y <- y + stats::rnorm(n, 0, rsd)

  # attrition: whole units lost
  dropped <- character()
  n_drop <- round(att * n)
  if (n_drop > 0) {
    dropped <- sample(meta$sample_id, n_drop)
    y[meta$sample_id %in% dropped] <- NA
  }

  # subsamples nested in the unit (leaves of one plant): not independent replicates
  out <- meta
  if (n_sub > 1) {
    idx <- rep(seq_len(n), each = n_sub)
    out <- meta[idx, , drop = FALSE]
    out$unit_id <- meta$sample_id[idx]
    out$subsample <- rep(seq_len(n_sub), times = n)
    out$sample_id <- paste0(out$unit_id, "_s", out$subsample)
    y_obs <- y[idx] + stats::rnorm(length(idx), 0, d$subsample_sd_ratio$value * rsd)
  } else {
    y_obs <- y
  }
  has_outcome <- !is.null(spec$primary_outcome) || !.suite_is_count_assay(spec, rules)
  out$primary_outcome <- if (has_outcome) round(y_obs, 6) else NA

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  paths <- c()
  mpath <- file.path(outdir, "sample_metadata_template.tsv")
  utils::write.table(out, mpath, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  paths <- c(paths, metadata = mpath)

  truth <- list(
    simulator = "simitall", simitall_version = as.character(utils::packageVersion("simitall")),
    handoff_version = spec$handoff_version, design = spec$design$chosen, model = h$model,
    outcome = spec$primary_outcome$name %||% "primary_outcome", seed = seed,
    intercept = mu, residual_sd = rsd, effect = eff,
    level_effect_pattern = d$level_effect_pattern$value,
    treatment_effects = treatment_effects,
    interaction_effects = if (length(interactions)) stats::setNames(as.list(rep(d$interaction_effect$value, length(interactions))), interactions) else list(),
    block_effects = block_effects, random_sd = random_sd,
    subsamples_per_unit = n_sub, subsample_sd = if (n_sub > 1) d$subsample_sd_ratio$value * rsd else 0,
    dropped_units = dropped, assumptions = assumptions,
    columns = as.list(p$col), outcome_simulated = has_outcome
  )

  if (.suite_is_count_assay(spec, rules)) {
    cnt <- .suite_simulate_counts(meta, frame, p, d, n_genes, seed, outdir)
    paths <- c(paths, counts = cnt$path)
    truth$counts <- cnt$truth
  }
  paths <- c(paths, truth = .suite_write_json(truth, file.path(outdir, "sim_truth.json")))
  message(sprintf("simitall: simulated %d units%s%s for design '%s' (%s).", n,
                  if (n_sub > 1) sprintf(" x %d subsamples", n_sub) else "",
                  if (!is.null(truth$counts)) sprintf(", %d genes", truth$counts$n_genes) else "",
                  spec$design$chosen %||% "?", h$model))
  invisible(paths)
}

.suite_is_count_assay <- function(spec, rules = .suite_rules("handoff")) {
  a <- tolower(unlist(c(spec$assays, spec$system)))
  any(a %in% tolower(unlist(rules$count_assays)))
}

# Count matrix through simulate_rnaseq_experiment(): one biological replicate
# per experimental unit; condition = first treatment term, batch = first block term.
.suite_simulate_counts <- function(meta, frame, p, d, n_genes, seed, outdir) {
  treat <- names(p$role)[p$role == "treatment"][1]
  blk <- names(p$role)[p$role == "block"][1]
  if (is.na(treat)) stop("Count simulation needs a treatment term in analysis.model.", call. = FALSE)
  cond <- frame[[treat]]
  design <- data.frame(sample = meta$sample_id, condition = as.character(cond),
                       batch = if (!is.na(blk)) paste0("batch_", frame[[blk]]) else "batch1",
                       biological_replicates = 1L, technical_replicates = 1L, stringsAsFactors = FALSE)
  # the reference condition (first level) must appear first
  design <- design[order(as.integer(cond)), , drop = FALSE]
  cd <- d$counts
  tmp <- tempfile("simitall_counts_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  res <- suppressMessages(simulate_rnaseq_experiment(
    out_prefix = file.path(tmp, "x"), replicate_mode = "biological", design = design,
    n_genes = n_genes %||% cd$n_genes, condition_effect_fraction = cd$de_fraction,
    condition_effect_sd = cd$de_log2fc_sd, batch_effect_sd = cd$batch_effect_sd,
    library_size_mean = cd$library_size_mean, seed = seed))
  counts <- utils::read.delim(res$counts, check.names = FALSE)
  sm <- utils::read.delim(res$sample_metadata, stringsAsFactors = FALSE)
  names(counts)[-1] <- sm$sample[match(names(counts)[-1], sm$library_id)]
  counts <- counts[, c("gene_id", meta$sample_id)]
  path <- file.path(outdir, "counts.tsv")
  utils::write.table(counts, path, sep = "\t", row.names = FALSE, quote = FALSE)
  ct <- utils::read.delim(res$condition_truth, stringsAsFactors = FALSE)
  list(path = path, truth = list(
    n_genes = nrow(counts), condition_term = treat, batch_term = if (is.na(blk)) NULL else blk,
    reference = levels(cond)[1], de_fraction = cd$de_fraction, de_log2fc_sd = cd$de_log2fc_sd,
    batch_effect_sd = cd$batch_effect_sd,
    de_genes = lapply(split(ct, ct$condition), function(x) stats::setNames(as.list(x$effect_log2), x$gene_id))))
}
