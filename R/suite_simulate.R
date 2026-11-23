#' Simulate data for a Biodes design package
#'
#' Reads a Biodes design folder (`study_spec.json`, `field_layout.csv`,
#' `sample_metadata_template.tsv`) and creates data for every experimental
#' unit of the layout. What is simulated is looked up in
#' `inst/suite/simulators.tsv` from the spec's fields (see [simitall_routes()]):
#'
#' * **Outcome** (`primary_outcome.type`): the primary outcome from
#'   `analysis.model` (treatment effects of the size the study was powered for,
#'   a variance for each blocking factor and random term, subsamples nested in
#'   their unit, attrition), as a continuous, binary (liability threshold),
#'   count (negative binomial) or ordinal score.
#' * **Genetics** (`assays` gwas/qtl_mapping/eqtl/..., `population.mating`,
#'   `study_type` breeding_genetics): genotypes for the levels of the genotype
#'   factor, from a diversity panel or a cross ([simulate_breeding()]), and
#'   genetic values from a phenotype model in `inst/suite/phenotype_models.tsv`
#'   (additive, dominance, epistasis, GxE) at the requested heritability. These
#'   genetic values are the genotype effects in the outcome.
#' * **Assays** (`assays`): bulk RNA-seq, eQTL, single-cell RNA-seq, ChIP-seq
#'   and ATAC-seq, one library per unit named by its `sample_id`.
#'
#' Parameters come from `inst/suite/simulator_parameters.tsv` (see
#' [simitall_simulator_parameters()]); a spec field overrides a default and
#' `params` overrides both. Every default used is listed under `assumptions`
#' in `sim_truth.json`, with the true values of every effect.
#'
#' @param biodes_dir Folder written by Biodes `write_design_package()`.
#' @param outdir Where to write the outputs (default: the same folder).
#' @param seed Random seed. Default: the design's `randomization_seed`.
#' @param effect Largest treatment difference. Default:
#'   `power.minimum_detectable_difference` from the spec.
#' @param residual_sd Residual SD of the unit means. Default:
#'   `power.residual_sd` from the spec.
#' @param n_genes Genes in count matrices (shortcut for the `n_genes` parameter
#'   of the RNA-seq, eQTL and single-cell simulators).
#' @param interaction_effect Size of the interaction between treatment factors
#'   (e.g. treatment x time). Default: `power.interaction_effect` from the spec,
#'   else `interaction_effect_ratio` x the main effect
#'   (`simulation_defaults.json`). 0 for no interaction. Genotype x treatment
#'   interactions in genetics studies come from the phenotype model instead.
#' @param params Parameter overrides, a list named by simulator, e.g.
#'   `list(suite_genetics = list(heritability = 0.8, phenotype_model = "gxe"),
#'   simulate_chipseq = list(n_peaks = 500))`.
#' @return Invisibly, the paths written.
#' @examples
#' \dontrun{
#' simitall_routes("studyforge_run")
#' simitall_from_biodes("studyforge_run")
#' }
#' @export
simitall_from_biodes <- function(biodes_dir, outdir = biodes_dir, seed = NULL, effect = NULL,
                                 residual_sd = NULL, n_genes = NULL, interaction_effect = NULL,
                                 params = NULL) {
  h <- .suite_read_biodes(biodes_dir)
  d <- .suite_rules("simulation_defaults")
  spec <- h$spec
  p <- h$parsed
  if (!is.null(n_genes)) {
    for (s in c("simulate_rnaseq_experiment", "simulate_rnaseq_from_gwas", "simulate_scrnaseq_experiment")) {
      params[[s]]$n_genes <- params[[s]]$n_genes %||% n_genes
    }
  }
  routes <- .suite_match_routes(spec)
  assumptions <- character()
  note <- function(x) assumptions <<- c(assumptions, x)
  pick <- function(arg, spec_value, default, what) {
    if (!is.null(arg)) return(arg)
    if (!is.null(spec_value)) return(spec_value)
    note(sprintf("%s = %s (default: %s)", what, format(default$value), default$note))
    default$value
  }
  out_route <- routes[routes$layer == "outcome", , drop = FALSE][1, ]
  if (out_route$status != "ready") {
    stop(sprintf("primary_outcome.type '%s' has no simulator yet (%s in simulators.tsv): %s",
                 out_route$matched, out_route$status, out_route$note), call. = FALSE)
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
  # factors in the layout that the model does not hold (Biodes analysis.terms.uncovered_factors):
  # tissue measured within each unit, a third factor, ...
  uncovered <- setdiff(intersect(unlist(spec$analysis$terms$uncovered_factors), names(meta)), c(names(p$role), unlist(p$col)))
  for (v in uncovered) frame[[v]] <- factor(as.character(meta[[v]]), levels = .suite_natural_sort(unique(as.character(meta[[v]]))))
  n <- nrow(frame)
  y <- rep(mu, n)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  paths <- c()

  # ---- genetics layer: genotypes of the lines and their genetic values ----
  gen_route <- routes[routes$layer == "genetics" & routes$status == "ready", , drop = FALSE]
  geno_term <- .suite_genotype_term(spec, p)
  cond_term <- .suite_condition_term(p, geno_term)
  genetics <- NULL
  if (nrow(gen_route) && is.na(geno_term)) {
    # e.g. BSA pools, populations, evolve-and-resequence: no lines to give genotypes to.
    # Simulate the rest and list genetics (and eQTL, which needs it) as not simulated.
    note(sprintf("genetics (%s = %s) not simulated: analysis.model has no genotype term (a treatment term whose analysis.terms meaning matches genetics.json genotype_pattern)",
                 gen_route$field[1], gen_route$matched[1]))
    no_geno <- routes$layer == "genetics" | routes$route == "assay_eqtl"
    routes$status[no_geno] <- "needs_genotype_factor"
    routes$note[no_geno] <- "No genotype factor in the design (lines as a treatment factor), so there are no lines to give genotypes to."
    gen_route <- gen_route[0, , drop = FALSE]
  }
  if (nrow(gen_route)) {
    lines <- levels(frame[[geno_term]])
    env_levels <- if (!is.na(cond_term) && cond_term != geno_term) levels(frame[[cond_term]]) else character()
    geno <- .suite_simulate_genotypes(gen_route[1, ], lines, spec, params, outdir, seed, note)
    gv <- .suite_genetic_values(geno$dosage, geno$variants, spec, params, rsd, env_levels, seed, note)
    genetics <- list(dosage = geno$dosage, variants = geno$variants, source = geno$info,
                     g = gv$g, gxe = gv$gxe, qtn = gv$qtn, info = gv$info)
    vcf <- .suite_write_vcf(geno$dosage, geno$variants, file.path(outdir, "genotypes.vcf"), "simitall genotypes of the design's lines")
    gdir <- file.path(outdir, "genetic_truth")
    gtab <- data.frame(line = lines, genetic_value = unname(gv$g[lines]), stringsAsFactors = FALSE)
    if (!is.null(gv$gxe)) for (l in colnames(gv$gxe)) gtab[[paste0("gxe_", l)]] <- unname(gv$gxe[lines, l])
    gpath <- file.path(gdir, "genetic_values.tsv")
    utils::write.table(gtab, gpath, sep = "\t", row.names = FALSE, quote = FALSE)
    qpath <- file.path(gdir, "qtn_truth.tsv")
    utils::write.table(gv$qtn, qpath, sep = "\t", row.names = FALSE, quote = FALSE)
    paths <- c(paths, genotypes = vcf, genetic_values = gpath, qtn_truth = qpath, geno$files)
    if (!is.null(spec$power$minimum_detectable_difference)) {
      note("genotype effects come from the phenotype model and heritability, not from power.minimum_detectable_difference")
    }
    if (is.null(params$suite_genetics$heritability)) note(sprintf("heritability = %s (default: simulator_parameters.tsv)", gv$info$heritability))
    if (is.null(params$suite_genetics$phenotype_model)) note(sprintf("phenotype model = %s (default: simulator_parameters.tsv)", gv$info$phenotype_model))
  }

  # ---- outcome: treatment main effects ----
  treat_vars <- intersect(names(p$role)[p$role == "treatment"], p$fixed_labels)
  if (!is.null(genetics) && !geno_term %in% treat_vars) treat_vars <- c(treat_vars, geno_term)
  treatment_effects <- numeric_effects <- list()
  for (v in treat_vars) {
    if (v %in% p$numeric) {
      # numeric factor (response surface): linear over its range, no curvature
      rng <- range(frame[[v]], na.rm = TRUE)
      slope <- if (diff(rng) > 0) eff / diff(rng) else 0
      numeric_effects[[v]] <- slope
      y <- y + slope * frame[[v]]
      next
    }
    lev <- levels(frame[[v]])
    e <- if (!is.null(genetics) && identical(v, geno_term)) genetics$g[lev]
         else stats::setNames(if (length(lev) > 1) seq(0, eff, length.out = length(lev)) else 0, lev)
    treatment_effects[[v]] <- as.list(e)
    y <- y + e[as.character(frame[[v]])]
  }
  # genotype x treatment from the phenotype model (whether or not the model fits it)
  if (!is.null(genetics$gxe)) {
    y <- y + genetics$gxe[cbind(as.character(frame[[geno_term]]), as.character(frame[[cond_term]]))]
  }
  # interactions among treatment factors: linear-by-linear, the first factor's
  # first-to-last effect changes by `inter` from the second factor's first to last level
  interactions <- p$fixed_labels[grepl(":", p$fixed_labels)]
  inter_default <- list(value = d$interaction_effect_ratio$value * eff, note = d$interaction_effect_ratio$note)
  is_gxe <- function(term) !is.null(genetics) && geno_term %in% all.vars(str2lang(term))
  treat_inter <- Filter(function(l) !is_gxe(l) && all(all.vars(str2lang(l)) %in% names(p$role)[p$role == "treatment"]), interactions)
  inter <- if (length(treat_inter)) as.numeric(pick(interaction_effect, spec$power$interaction_effect, inter_default, "interaction effect")) else 0
  interaction_effects <- list()
  for (term in interactions) {
    vars <- all.vars(str2lang(term))
    key <- do.call(paste, c(lapply(vars, function(v) as.character(frame[[v]])), sep = ":"))
    if (is_gxe(term)) {
      # cells hold the genotype x treatment values already added above (0 without GxE)
      add <- if (!is.null(genetics$gxe) && identical(sort(vars), sort(c(geno_term, cond_term)))) {
        genetics$gxe[cbind(as.character(frame[[geno_term]]), as.character(frame[[cond_term]]))]
      } else rep(0, n)
      cells <- tapply(add, key, `[`, 1)
      interaction_effects[[term]] <- list(size = if (is.null(genetics$gxe)) 0 else genetics$info$gxe_share,
                                          pattern = "phenotype_model_gxe", cells = as.list(cells))
      next
    }
    size <- if (term %in% treat_inter) inter else 0
    score <- lapply(vars, function(v) {
      x <- frame[[v]]
      if (is.numeric(x)) { r <- range(x, na.rm = TRUE); if (diff(r) > 0) (x - r[1]) / diff(r) else 0 * x }
      else if (nlevels(x) > 1) (as.integer(x) - 1) / (nlevels(x) - 1) else 0 * as.integer(x)
    })
    add <- size * Reduce(`*`, score)
    y <- y + add
    cells <- tapply(add, key, `[`, 1)
    interaction_effects[[term]] <- list(size = size, pattern = "linear_by_linear", cells = as.list(cells))
  }
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
  # uncovered factors: main effects spread over their levels like treatment levels
  uncovered_effects <- list()
  for (v in uncovered) {
    lev <- levels(frame[[v]])
    e <- stats::setNames(if (length(lev) > 1) seq(0, eff, length.out = length(lev)) else 0, lev)
    uncovered_effects[[v]] <- as.list(e)
    y <- y + e[as.character(frame[[v]])]
  }
  # rows of one unit (one per within-unit level, e.g. leaf and root of one plant) share that unit's effect
  unit_sd <- 0
  unit_key <- NULL
  if (length(uncovered)) {
    lay_cols <- setdiff(intersect(names(h$layout), names(meta)), uncovered)
    unit_key <- do.call(paste, c(lapply(meta[lay_cols], as.character), sep = "|"))
    if (anyDuplicated(unit_key)) {
      unit_sd <- d$random_sd_ratio$value * rsd
      u <- factor(unit_key, levels = unique(unit_key))
      y <- y + stats::rnorm(nlevels(u), 0, unit_sd)[as.integer(u)]
    } else unit_key <- NULL
  }
  y_mean <- y
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
  idx <- seq_len(n)
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
  # outcome family (simulators.tsv outcome_* rows)
  oargs <- .suite_args("suite_outcome", spec, params)$args
  fam <- .suite_outcome_family(out_route$route, y_obs, y_mean[idx], oargs)
  assay_routes <- routes[routes$layer == "assay" & routes$status == "ready", , drop = FALSE]
  has_outcome <- !is.null(spec$primary_outcome) || !nrow(assay_routes)
  out$primary_outcome <- if (has_outcome) fam$y else NA
  if (!is.null(fam$event)) out$event <- if (has_outcome) fam$event else NA

  mpath <- file.path(outdir, "sample_metadata_template.tsv")
  utils::write.table(out, mpath, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  paths <- c(metadata = mpath, paths)

  truth <- list(
    simulator = "simitall", simitall_version = as.character(utils::packageVersion("simitall")),
    handoff_version = spec$handoff_version, design = spec$design$chosen, model = h$model,
    outcome = spec$primary_outcome$name %||% "primary_outcome", seed = seed,
    outcome_type = fam$type, outcome_model = fam$model,
    intercept = mu, residual_sd = rsd, effect = eff,
    level_effect_pattern = d$level_effect_pattern$value,
    treatment_effects = treatment_effects, numeric_effects = numeric_effects,
    interaction_effects = interaction_effects,
    block_effects = block_effects, random_sd = random_sd,
    uncovered_effects = if (length(uncovered_effects)) uncovered_effects,
    within_unit_sd = if (!is.null(unit_key)) list(sd = unit_sd, note = "rows of one unit (one per level of the within-unit factor) share this effect; they are not independent units",
                                                  unit_columns = setdiff(intersect(names(h$layout), names(meta)), uncovered)),
    subsamples_per_unit = n_sub, subsample_sd = if (n_sub > 1) d$subsample_sd_ratio$value * rsd else 0,
    dropped_units = dropped, assumptions = NULL,
    columns = as.list(p$col), outcome_simulated = has_outcome,
    routes = lapply(seq_len(nrow(routes)), function(i) as.list(routes[i, c("route", "layer", "simulator", "status", "field", "matched")]))
  )
  if (!is.null(genetics)) {
    truth$genetics <- c(list(genotype_term = geno_term, environment_term = if (length(env_levels)) cond_term,
                             n_lines = length(lines), n_markers = nrow(genetics$dosage),
                             genotypes = "genotypes.vcf", genetic_values = "genetic_truth/genetic_values.tsv",
                             qtn_truth = "genetic_truth/qtn_truth.tsv"),
                        genetics$info, list(genotype_source = genetics$source))
  }

  # ---- assays ----
  if (nrow(assay_routes)) {
    owners <- c("assay_eqtl", "assay_rnaseq", "assay_single_cell", "assay_small_rnaseq", "assay_microbiome")
    counts_owner <- owners[owners %in% assay_routes$route][1]
    if ("assay_eqtl" %in% assay_routes$route && "assay_rnaseq" %in% assay_routes$route) {
      note("RNA-seq and eQTL both requested: counts.tsv comes from the eQTL simulator (it includes the treatment effects)")
      assay_routes <- assay_routes[assay_routes$route != "assay_rnaseq", , drop = FALSE]
    }
    ctx <- list(spec = spec, params = params, meta = meta, frame = frame, p = p, outdir = outdir, seed = seed, uncovered = uncovered,
                cond_term = cond_term, geno_term = geno_term, genetics = genetics)
    truth$assays <- list()
    for (i in seq_len(nrow(assay_routes))) {
      r <- assay_routes[i, ]
      # one assay owns counts.tsv; another count assay keeps its matrix in its own folder
      own <- identical(r$route, counts_owner)
      ctx$counts_file <- if (own) file.path(outdir, "counts.tsv") else file.path(outdir, sub("^assay_", "", r$route), "unit_counts.tsv")
      res <- .suite_assay_runner(r$simulator)(ctx)
      if (!is.null(res$counts_truth)) {
        if (own) truth$counts <- res$counts_truth else res$truth$counts <- c(list(file = sub(paste0("^", outdir, "/"), "", ctx$counts_file)), res$counts_truth)
      }
      truth$assays[[r$route]] <- c(list(simulator = r$simulator),
                                   res$truth, list(parameters = res$args$args, parameter_source = as.list(res$args$source)))
      paths <- c(paths, stats::setNames(res$paths, paste0(sub("^assay_", "", r$route), ".", names(res$paths))))
    }
    if (!is.null(truth$counts)) paths <- c(paths, counts = file.path(outdir, "counts.tsv"))
  }
  treat_only <- routes[routes$status == "treatment", , drop = FALSE]
  if (nrow(treat_only)) truth$treatment_assays <- as.list(treat_only$matched)   # e.g. CRISPR lines: a treatment, not data
  skipped <- routes[!routes$status %in% c("ready", "treatment"), , drop = FALSE]
  if (nrow(skipped)) {
    truth$not_simulated <- lapply(seq_len(nrow(skipped)), function(i) list(
      route = skipped$route[i], requested = skipped$matched[i], status = skipped$status[i], note = skipped$note[i]))
    message("simitall: not simulated (", paste(sprintf("%s: %s", skipped$matched, skipped$status), collapse = "; "), ")")
  }
  truth$assumptions <- assumptions
  paths <- c(paths, truth = .suite_write_json(truth, file.path(outdir, "sim_truth.json")))
  message(sprintf("simitall: simulated %d units%s (%s outcome)%s%s for design '%s' (%s).", n,
                  if (n_sub > 1) sprintf(" x %d subsamples", n_sub) else "", fam$type,
                  if (!is.null(genetics)) sprintf(", %d lines x %d markers (%s)", length(lines), nrow(genetics$dosage), genetics$source$source) else "",
                  if (length(truth$assays)) paste0(", assays: ", paste(sub("^assay_", "", names(truth$assays)), collapse = ", ")) else "",
                  spec$design$chosen %||% "?", h$model))
  invisible(paths)
}

# Turn the Gaussian observation into the requested outcome type.
.suite_outcome_family <- function(route, y_obs, y_mean, a) {
  ok <- !is.na(y_obs)
  switch(route,
    outcome_continuous = list(type = "continuous", y = round(y_obs, 6), model = list(family = "gaussian")),
    outcome_binary = {
      thr <- stats::quantile(y_obs[ok], 1 - a$prevalence, names = FALSE)
      list(type = "binary", y = ifelse(ok, as.integer(y_obs > thr), NA_integer_),
           model = list(family = "liability_threshold", prevalence = a$prevalence, threshold = thr,
                        note = "treatment_effects are on the liability (Gaussian) scale"))
    },
    outcome_count = {
      m <- pmax(y_mean, a$min_mean)
      yc <- if (identical(a$count_distribution, "poisson")) stats::rpois(length(m), m)
            else stats::rnbinom(length(m), mu = m, size = a$dispersion)
      list(type = "count", y = ifelse(ok, yc, NA_integer_),
           model = list(family = a$count_distribution, link = "identity", dispersion = a$dispersion, min_mean = a$min_mean,
                        note = "unit mean = intercept + effects (no residual draw); counts drawn around it"))
    },
    outcome_ordinal = {
      pr <- a$category_probabilities / sum(a$category_probabilities)
      cuts <- stats::quantile(y_obs[ok], cumsum(pr)[-length(pr)], names = FALSE)
      list(type = "ordinal", y = ifelse(ok, findInterval(y_obs, cuts) + 1L, NA_integer_),
           model = list(family = "ordinal_threshold", category_probabilities = pr, thresholds = cuts,
                        note = "score 1..k; treatment_effects are on the latent (Gaussian) scale"))
    },
    outcome_survival = {
      # log-normal accelerated failure time: log T = log(median) + time_sd * z, z the standardised Gaussian outcome
      z <- (y_obs - mean(y_obs[ok])) / stats::sd(y_obs[ok])
      t <- a$median_time * exp(a$time_sd * z)
      follow <- stats::quantile(t[ok], 1 - a$censored_fraction, names = FALSE)
      ev <- as.integer(t <= follow)
      list(type = "survival", y = ifelse(ok, round(pmin(t, follow), 3), NA_real_), event = ifelse(ok, ev, NA_integer_),
           model = list(family = "lognormal_aft", median_time = a$median_time, time_sd = a$time_sd,
                        follow_up = follow, censored_fraction = a$censored_fraction,
                        scale = list(sd_outcome = stats::sd(y_obs[ok])),
                        note = "primary_outcome = time to event or end of follow-up; event 1 = event, 0 = censored. treatment_effects are on the Gaussian scale: an effect e multiplies time by exp(time_sd * e / scale$sd_outcome)"))
    },
    stop("No outcome model for route ", route, call. = FALSE))
}

# Kept for simitall_power_sweep(): does the spec ask for a counts-only assay?
.suite_is_count_assay <- function(spec, rules = NULL) {
  r <- .suite_match_routes(spec)
  any(r$route %in% c("assay_rnaseq", "assay_eqtl", "assay_single_cell"))
}
