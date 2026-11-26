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
  # the outcome's scale (Biodes 0.6.4: primary_outcome.unit, baseline_mean, power.residual_sd, power.cv in %)
  po <- spec$primary_outcome %||% list()
  sd_spec <- spec$power$residual_sd %||% if (!is.null(spec$power$cv) && !is.null(po$baseline_mean)) as.numeric(spec$power$cv) / 100 * as.numeric(po$baseline_mean)
  if (is.null(spec$power$residual_sd) && !is.null(sd_spec)) note(sprintf("residual SD = %s (power.cv %s%% x baseline mean %s)", format(sd_spec), format(spec$power$cv), format(po$baseline_mean)))
  sd_asked <- grepl("^not stated", spec$assumptions[["power.residual_sd"]] %||% "")
  # an outcome with a unit (or an SD Biodes asked for) but no spread: SD 10 / mean 100 would be read as that unit,
  # so the outcome is simulated on a standardized scale (residual SD 1, mean 0, effects in SDs) and says so.
  # Count outcomes draw counts around the mean and keep their defaults.
  standardized <- is.null(residual_sd) && is.null(sd_spec) && (!is.null(po$unit) || sd_asked) && out_route$route != "outcome_count"
  outcome_scale <- list(unit = po$unit %||% "none stated (arbitrary units)", standardized = standardized)
  if (standardized) {
    rsd <- 1; mu <- 0
    stated <- spec$power$minimum_detectable_difference
    eff <- as.numeric(effect %||% d$effect_in_residual_sd$value)
    outcome_scale <- c(outcome_scale, list(residual_sd_known = FALSE, residual_sd = 1, intercept = 0,
      stated_effect = if (!is.null(stated)) sprintf("%s %s", format(stated), po$unit %||% ""),
      note = sprintf("the residual SD%s was not stated%s: the outcome is standardized (residual SD 1, mean 0) and every effect is in residual SDs; the treatment effect is %s SD. Power depends on the real SD.",
                     if (!is.null(po$unit)) paste0(" in ", po$unit) else "", if (sd_asked) " (Biodes asked for it)" else "", format(eff))))
    outcome_scale <- outcome_scale[!vapply(outcome_scale, is.null, logical(1))]
    note(outcome_scale$note)
  } else {
    rsd <- as.numeric(pick(residual_sd, sd_spec, d$residual_sd, "residual SD"))
    eff_default <- list(value = d$effect_in_residual_sd$value * rsd, note = d$effect_in_residual_sd$note)
    eff <- as.numeric(pick(effect, spec$power$minimum_detectable_difference, eff_default, "treatment effect"))
    mu <- as.numeric(pick(NULL, po$baseline_mean, d$intercept, "intercept"))
    outcome_scale <- c(outcome_scale, list(residual_sd_known = !is.null(residual_sd) || !is.null(sd_spec), residual_sd = rsd, intercept = mu,
      residual_sd_source = if (!is.null(residual_sd)) "argument" else if (!is.null(spec$power$residual_sd)) "power.residual_sd" else if (!is.null(sd_spec)) "power.cv x primary_outcome.baseline_mean" else "default"))
  }
  att <- as.numeric(pick(NULL, spec$power$attrition, d$attrition, "attrition"))
  n_sub <- as.integer(spec$subsamples_per_unit %||% 1L)
  if (is.na(n_sub) || n_sub < 1L) n_sub <- 1L

  meta <- h$meta
  if (nrow(meta) != nrow(h$layout)) stop("sample_metadata_template.tsv and field_layout.csv have different row counts.", call. = FALSE)
  frame <- .suite_frame(meta, p)
  # factors in the layout that the model does not hold (Biodes analysis.terms.uncovered_factors):
  # tissue measured within each unit, a third factor, ...
  # occasion factors measured on the same units (Biodes repeated_factors, e.g. the same mice at 24/48/72 h) count as
  # within-unit too when the model does not hold them
  repeated <- intersect(unlist(spec$analysis$terms$repeated_factors %||% spec$repeated_factors), names(meta))
  uncovered <- setdiff(intersect(unique(c(unlist(spec$analysis$terms$uncovered_factors), repeated)), names(meta)), c(names(p$role), unlist(p$col)))
  spec_levels <- function(v) { f <- Filter(function(z) identical(z$name, v), spec$factors); if (length(f)) as.character(unlist(f[[1]]$levels)) }
  for (v in uncovered) frame[[v]] <- factor(as.character(meta[[v]]), levels = .suite_level_order(as.character(meta[[v]]), spec_levels(v)))
  # the analysis goal (Biodes schema 0.6 analysis.primary_goal) decides which truth is planted; see suite_goals.R
  gl <- .suite_goal_setup(spec, p, frame, eff, note, standardized)
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
  # time_course: the time factor follows the goal's shape; within_group: only the non-reference groups change
  tc_var <- tc_group <- tc_levels <- tc_pattern <- tr_var <- trend <- goal_pair <- NULL
  if (gl$goal == "time_course" && is.factor(frame[[gl$vars[1]]])) {
    tc_var <- gl$vars[1]
    tpos <- .suite_time_positions(tc_var, frame, p, spec)
    if (is.null(tpos)) {
      tpos <- stats::setNames(seq_len(nlevels(frame[[tc_var]])) - 1, levels(frame[[tc_var]]))
      note(sprintf("time_course: levels of %s carry no sampling times; the shape follows level order", tc_var))
    }
    tcol <- if (tc_var %in% names(p$col)) p$col[[tc_var]] else tc_var
    funit <- unlist(Filter(function(z) identical(z$name, tcol), spec$factors %||% list())[[1]]$unit %||% NULL)
    if (is.null(gl$shape)) note("time_course shape, peak and half-time not stated: monotone (change grows with time)")
    tc_pattern <- gl$size * .suite_time_shape(tpos, gl$shape, gl$peak, gl$half_time, funit, note)
    if (identical(gl$scope, "within_group")) {
      # Biodes 0.6.2 group_factor (a layout column) and group (the level that changes, as the question words it)
      gf <- gl$group_factor
      cand <- if (!is.null(gf)) c(names(p$col)[vapply(p$col, function(z) identical(z, gf), logical(1))], gf)
      cand <- intersect(cand, setdiff(names(frame), tc_var))
      tc_group <- if (length(cand)) cand[1] else setdiff(treat_vars, c(tc_var, if (!is.null(genetics)) geno_term))[1]
      if (is.na(tc_group)) { tc_group <- NULL; note("time_course scope within_group: no group factor besides time; every unit follows the pattern") }
      else {
        lev <- levels(frame[[tc_group]])
        hit <- if (!is.null(gl$group)) { e <- lev[tolower(lev) == tolower(gl$group)]; if (length(e)) e else lev[grepl(tolower(gl$group), tolower(lev), fixed = TRUE)] }
        tc_levels <- if (length(hit)) hit else lev[-1]
        if (!is.null(gl$group) && !length(hit)) note(sprintf("time_course within_group: group '%s' is not a level of %s (%s); the non-reference levels change", gl$group, tc_group, paste(lev, collapse = ", ")))
        note(sprintf("time_course within_group: %s levels are equal at the first time point; only %s change over time", tc_group, paste(tc_levels, collapse = ", ")))
      }
    }
  } else if (gl$goal == "trend") {
    tr_var <- gl$vars[1]
    trend <- .suite_trend_x(tr_var, frame, p, spec, note)
    if (is.null(trend)) { note(sprintf("trend: %s has no numeric values; planted as treatment_effect", tr_var)); tr_var <- NULL }
    else {
      trend$slope <- if (gl$size_stated) gl$size else if (diff(range(trend$x)) > 0) gl$size / diff(range(trend$x)) else 0
      if (!gl$size_stated) note(sprintf("trend slope = effect %s over the range of %s%s", format(gl$size), tr_var,
        if (identical(trend$scale, "log")) sprintf(" (log%s of value + %s)", format(trend$log_base), format(trend$log_offset)) else ""))
    }
  } else if (gl$goal == "interaction") goal_pair <- gl$vars
  for (v in treat_vars) {
    if (identical(v, tr_var)) {
      e <- trend$slope * (trend$x - min(trend$x))
      if (is.numeric(frame[[v]])) { numeric_effects[[v]] <- trend$slope; y <- y + e; next }
      treatment_effects[[v]] <- as.list(e)
      y <- y + e[as.character(frame[[v]])]
      next
    }
    if (v %in% p$numeric) {
      # numeric factor (response surface): linear over its range, no curvature
      rng <- range(frame[[v]], na.rm = TRUE)
      slope <- if (diff(rng) > 0) eff / diff(rng) else 0
      numeric_effects[[v]] <- slope
      y <- y + slope * frame[[v]]
      next
    }
    lev <- levels(frame[[v]])
    tpos <- .suite_time_positions(v, frame, p, spec)
    e <- if (!is.null(genetics) && identical(v, geno_term)) genetics$g[lev]
         else if (identical(v, tc_var)) (if (is.null(tc_group)) tc_pattern else 0 * tc_pattern)
         else if (identical(v, tc_group)) stats::setNames(rep(0, length(lev)), lev)
         else if (!is.null(tpos)) eff * (tpos - min(tpos)) / diff(range(tpos))   # sampling times: effect grows with real time
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
  goal_inter <- if (!is.null(goal_pair)) goal_pair else if (!is.null(tc_group)) c(tc_group, tc_var)
  other_inter <- Filter(function(l) is.null(goal_inter) || !setequal(all.vars(str2lang(l)), goal_inter), treat_inter)
  inter <- if (length(other_inter)) as.numeric(pick(interaction_effect, spec$power$interaction_effect, inter_default, "interaction effect")) else 0
  interaction_effects <- list()
  same <- function(a, b) length(a) == length(b) && setequal(a, b)
  goal_label <- NULL
  if (!is.null(goal_inter) && !any(vapply(interactions, function(l) same(all.vars(str2lang(l)), goal_inter), logical(1)))) {
    interactions <- c(interactions, paste(goal_inter, collapse = ":"))
    note(sprintf("goal %s: %s is not in analysis.model; planted anyway", gl$goal, paste(goal_inter, collapse = ":")))
  }
  for (term in interactions) {
    vars <- all.vars(str2lang(term))
    if (!is.null(goal_inter) && same(vars, goal_inter)) goal_label <- term
    key <- do.call(paste, c(lapply(vars, function(v) as.character(frame[[v]])), sep = ":"))
    if (!is.null(tc_group) && same(vars, goal_inter)) {
      ch <- frame[[tc_group]] %in% tc_levels
      add <- unname(tc_pattern[as.character(frame[[tc_var]])]) * ch
      y <- y + add
      interaction_effects[[term]] <- list(size = gl$size, pattern = "time_course_within_group", cells = as.list(tapply(add, key, `[`, 1)))
      next
    }
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
    size <- if (!is.null(goal_pair) && same(vars, goal_pair)) gl$size else if (term %in% treat_inter) inter else 0
    score <- lapply(vars, function(v) {
      x <- frame[[v]]
      tpos <- .suite_time_positions(v, frame, p, spec)
      if (!is.null(tpos)) return(unname(((tpos - min(tpos)) / diff(range(tpos)))[as.character(x)]))
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
    tpos <- .suite_time_positions(v, frame, p, spec)
    e <- if (identical(v, tc_var)) (if (is.null(tc_group)) tc_pattern else 0 * tc_pattern)
      else if (identical(v, tc_group)) stats::setNames(rep(0, length(lev)), lev)
      else if (identical(v, tr_var)) trend$slope * (trend$x - min(trend$x))
      else if (!is.null(tpos)) eff * (tpos - min(tpos)) / diff(range(tpos)) else
      stats::setNames(if (length(lev) > 1) seq(0, eff, length.out = length(lev)) else 0, lev)
    uncovered_effects[[v]] <- as.list(e)
    y <- y + e[as.character(frame[[v]])]
  }
  # rows of one unit (one per within-unit level, e.g. leaf and root of one plant) share that unit's effect
  unit_sd <- 0
  unit_key <- NULL
  # a repeated factor the model holds still needs the unit effect when the model has no random term for the unit
  unit_random <- any(vapply(p$random, function(r) any(p$role[r$vars] %in% "unit", na.rm = TRUE), logical(1)))
  # paired design (Biodes analysis.terms.paired: the repeated occasion is the only factor, each block is one
  # animal): the block term already is the animal effect, no extra unit effect
  paired <- isTRUE(spec$analysis$terms$paired)
  rep_in_model <- length(intersect(repeated, unlist(p$col))) > 0 && !unit_random && !paired
  if (length(uncovered) || rep_in_model) {
    lay_cols <- setdiff(intersect(names(h$layout), names(meta)), c(uncovered, repeated))
    unit_key <- do.call(paste, c(lapply(meta[lay_cols], as.character), sep = "|"))
    if (anyDuplicated(unit_key)) {
      unit_sd <- d$random_sd_ratio$value * rsd
      u <- factor(unit_key, levels = unique(unit_key))
      y <- y + stats::rnorm(nlevels(u), 0, unit_sd)[as.integer(u)]
    } else unit_key <- NULL
  }
  # repeated measures: correlation among one unit's occasions (suite_outcome within_subject; suite_within.R)
  ws <- .suite_within_series(meta, names(h$layout), repeated, spec)
  wsub <- .suite_within_subject(ws, n, rsd, .suite_args("suite_outcome", spec, params)$args)
  if (!is.null(wsub$truth$note) && !identical(wsub$truth$structure, "random_intercept")) note(wsub$truth$note)
  y <- y + wsub$add
  y_mean <- y
  y <- y + wsub$eps

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
    intercept = mu, residual_sd = rsd, effect = eff, outcome_scale = outcome_scale,
    level_effect_pattern = d$level_effect_pattern$value,
    treatment_effects = treatment_effects, numeric_effects = numeric_effects,
    interaction_effects = interaction_effects,
    block_effects = block_effects, random_sd = random_sd,
    uncovered_effects = if (length(uncovered_effects)) uncovered_effects,
    within_subject = wsub$truth,
    within_unit_sd = if (!is.null(unit_key)) list(sd = unit_sd, note = "rows of one unit (one per level of the within-unit factor) share this effect; they are not independent units",
                                                  unit_columns = setdiff(intersect(names(h$layout), names(meta)), c(uncovered, repeated))),
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
  # the goal in count assays: a share of features follows the goal's pattern (interaction, time course, trend)
  gpat <- .suite_goal_pattern(gl, frame, p, spec, tc_var, tc_group, tc_pattern, tr_var, trend, tc_levels)
  goal_features <- list()
  if (!is.null(gpat) && length(truth$assays)) {
    fa <- .suite_args("suite_goal_features", spec, params)
    gterm <- goal_label %||% tc_var %||% tr_var
    per_unit <- if (!is.null(tr_var)) diff(range(trend$x))
    # goal_effect.size on log2 (Biodes goal_effect.scale for counts) sizes every goal feature; a trend size is per unit
    if (!is.null(gl$feature_size)) { fa$args$log2fc <- abs(gl$feature_size) * (per_unit %||% 1); fa$args$log2fc_sd <- 0
      fa$source[c("log2fc", "log2fc_sd")] <- "spec analysis.goal_effect size (scale log2)" }
    for (route in intersect(c("assay_eqtl", "assay_rnaseq", "assay_single_cell", "assay_small_rnaseq", "assay_microbiome"), names(truth$assays))) {
      ct <- truth$assays[[route]]$counts %||% if (!is.null(truth$counts) && is.null(truth$assays[[route]]$counts)) c(list(file = "counts.tsv"), truth$counts)
      if (is.null(ct$file) || !file.exists(file.path(outdir, ct$file))) next
      gt <- .suite_goal_features(file.path(outdir, ct$file), gpat, meta$sample_id, names(unlist(unname(ct$de_genes))), fa$args, gterm, per_unit, seed)
      if (is.null(gt)) next
      id <- sub("^assay_", "", route)
      dir.create(file.path(outdir, "truth"), showWarnings = FALSE)
      f <- file.path("truth", paste0(id, "_goal.tsv"))
      utils::write.table(gt, file.path(outdir, f), sep = "\t", row.names = FALSE, quote = FALSE)
      goal_features[[id]] <- list(file = f, counts = ct$file, n_features = nrow(gt), term = gterm,
        effect_scale = if (!is.null(per_unit)) sprintf("log2 fold change per %s (column change: over the factor's range)",
                         if (identical(trend$scale, "log")) sprintf("log%s unit of (%s + %s)", format(trend$log_base), trend$unit %||% "x", format(trend$log_offset)) else trend$unit %||% "unit")
                       else if (gl$goal == "interaction") "log2 difference of differences, last levels vs references" else "log2 largest change between time points",
        parameters = fa$args, parameter_source = as.list(fa$source))
    }
  }
  truth$goal_truth <- .suite_goal_extra(gl, truth, frame, p, list(
    treat_vars = setdiff(treat_vars, if (!is.null(genetics)) geno_term), tc_var = tc_var, tc_group = tc_group, tc_levels = tc_levels,
    tc_pattern = tc_pattern, trend = trend, tr_var = tr_var, goal_label = goal_label, fam = fam, genetics = genetics,
    has_outcome = has_outcome, assay_routes = assay_routes$route, features = goal_features))
  truth$assumptions <- assumptions
  paths <- c(paths, truth = .suite_write_json(truth, file.path(outdir, "sim_truth.json")))
  paths <- c(paths, manifest = .suite_write_manifest(outdir, truth, has_outcome, fam$type, cond_term))
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

# Sampling times of a time factor's levels, so effects follow real time (0, 6, 24, 48 h are not evenly spaced).
# From the spec factor's `values` (same order as its `levels`) when given, else the number in each level label
# ("24 h", "day_10", "T48"). NULL when the variable is not a time column or the levels carry no numbers.
.suite_time_positions <- function(v, frame, p, spec) {
  col <- if (v %in% names(p$col)) p$col[[v]] else v
  time_cols <- c(.suite_rules("model_terms")$terms$time$columns, "hour", "hours", "week", "weeks", "month", "sampling_time")
  if (!tolower(col) %in% time_cols) return(NULL)
  x <- frame[[v]]
  if (!is.factor(x) || nlevels(x) < 2) return(NULL)
  lev <- levels(x)
  f <- Filter(function(z) identical(z$name, col), spec$factors %||% list())
  vals <- if (length(f)) suppressWarnings(as.numeric(unlist(f[[1]]$values))) else numeric()
  sl <- if (length(f)) as.character(unlist(f[[1]]$levels)) else character()
  num <- if (length(vals) && length(vals) == length(sl) && all(lev %in% sl)) vals[match(lev, sl)] else
    suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+(\\.[0-9]+)?).*$", "\\1", lev)))
  if (anyNA(num) || length(unique(num)) < 2) return(NULL)
  stats::setNames(num, lev)
}
