# Analysis goals (Biodes schema 0.6: analysis.goals, primary_goal, goal_terms, goal_effect).
# Simitall plants the truth the primary goal asks about and writes it to sim_truth.json -> goal_truth.
# What each goal plants: inst/suite/goals.tsv (mirrors Biodes inst/registry/analysis_goals.json).
# Specs without goals (schema 0.5 and older) get treatment_effect, today's behaviour.

# The primary goal and what it needs: the frame variables of its terms, its effect size and, for a time
# course, the shape. `eff` is the effect Simitall already uses (power.minimum_detectable_difference or default).
.suite_goal_setup <- function(spec, p, frame, eff, note) {
  tab <- .suite_tsv("goals")
  goal <- spec$analysis$primary_goal %||% unlist(spec$analysis$goals)[1] %||% "treatment_effect"
  row <- tab[tab$goal == goal, , drop = FALSE]
  if (!nrow(row)) {
    note(sprintf("primary_goal '%s' is not in inst/suite/goals.tsv; truth planted for treatment_effect", goal))
    goal <- "treatment_effect"; row <- tab[tab$goal == goal, , drop = FALSE]
  }
  ge <- spec$analysis$goal_effect[[goal]] %||% list()
  # goal_effect.scale (Biodes 0.6.1): identity/outcome = the primary outcome's scale; log2 = count features
  # (the engine's scale). A size on log2 sizes the goal's features; the outcome then uses the treatment effect.
  gscale <- ge$scale %||% "outcome"
  on_outcome <- gscale %in% c("outcome", "identity")
  stated <- !is.null(ge$size) && on_outcome
  size <- if (stated) as.numeric(ge$size) else eff
  feature_size <- if (!is.null(ge$size) && identical(gscale, "log2")) as.numeric(ge$size)
  if (!is.null(ge$size) && !on_outcome && is.null(feature_size))
    note(sprintf("goal effect size %s on scale %s: no %s-scale data simulated for it; the outcome uses the treatment effect", format(ge$size), gscale, gscale))
  if (!stated && !is.null(spec$analysis$primary_goal) && on_outcome)
    note(sprintf("goal effect size for %s = %s (analysis.goal_effect.%s.size not stated: the treatment effect)", goal, format(size), goal))
  # a stated peak means a rise and return (transient); a stated half-time alone, a saturating rise
  shape <- ge$shape %||% if (!is.null(ge$peak)) "transient" else if (!is.null(ge$half_time)) "saturating"
  if (is.null(ge$shape) && !is.null(shape)) note(sprintf("time_course shape not stated: %s (a %s is given)", shape, if (shape == "transient") "peak" else "half-time"))
  # goal terms are the model's generic terms (treatment, treatment2, time, env); find each one's frame variable
  terms <- unlist(spec$analysis$goal_terms[[goal]])
  var_of <- function(t) {
    if (t %in% names(frame)) return(t)
    col <- unlist(spec$analysis$terms[[t]])[1]
    if (!is.null(col) && col %in% names(frame)) return(col)
    hit <- names(p$col)[vapply(p$col, function(z) identical(z, col), logical(1))]
    if (length(hit) && hit[1] %in% names(frame)) hit[1] else NA_character_
  }
  vars <- vapply(terms, var_of, "")
  if (anyNA(vars)) {
    note(sprintf("goal %s: term %s not in the layout; truth planted for treatment_effect", goal, paste(terms[is.na(vars)], collapse = ", ")))
    goal <- "treatment_effect"; row <- tab[tab$goal == goal, , drop = FALSE]; terms <- vars <- character()
  }
  list(goal = goal, test = row$test, planted = row$planted, truth_in = row$truth, note = row$note,
       definition = ge$definition, size = size, size_stated = stated, scale = gscale, feature_size = feature_size,
       terms = unname(terms), vars = unname(vars),
       shape = shape, peak = ge$peak, half_time = ge$half_time, scope = ge$scope, group = ge$group, group_factor = ge$group_factor,
       per_unit = ge$per_unit)
}

# hours in a time unit (Biodes .sf_unit_name names); NA when unknown
.suite_unit_hours <- function(u) {
  u <- tolower(u %||% "")
  h <- c(s = 1 / 3600, sec = 1 / 3600, second = 1 / 3600, seconds = 1 / 3600, min = 1 / 60, minute = 1 / 60, minutes = 1 / 60,
         h = 1, hr = 1, hour = 1, hours = 1, d = 24, day = 24, days = 24, wk = 168, week = 168, weeks = 168,
         mo = 730, month = 730, months = 730, y = 8760, yr = 8760, year = 8760, years = 8760)
  if (u %in% names(h)) h[[u]] else NA_real_
}

# Pattern over time points, 0 at the first time point and a largest change between time points of 1.
# pos: sampling times (same unit as the factor); shape: monotone (linear), saturating (half-time),
# transient (rises to a peak and returns), delayed (logistic switch at the half-time).
.suite_time_shape <- function(pos, shape = NULL, peak = NULL, half_time = NULL, factor_unit = NULL, note = function(x) NULL) {
  shape <- shape %||% "monotone"
  s <- pos - min(pos); span <- max(s)
  at <- function(x, what, default) {
    if (is.null(x$value)) { note(sprintf("time_course %s not stated: %s%s", what, format(signif(default + min(pos), 3)), if (is.null(factor_unit)) "" else paste0(" ", factor_unit))); return(default) }
    v <- as.numeric(x$value)
    fh <- .suite_unit_hours(factor_unit); xh <- .suite_unit_hours(x$unit)
    if (!is.na(fh) && !is.na(xh)) v <- v * xh / fh else if (!is.null(x$unit) && !identical(x$unit, factor_unit))
      note(sprintf("time_course %s unit '%s' taken as the time factor's unit", what, x$unit))
    v - min(pos)
  }
  f <- switch(shape,
    monotone = s,
    saturating = { h <- at(half_time %||% list(), "half_time", span / 4); 1 - exp(-log(2) * s / max(h, 1e-9)) },
    transient = { pk <- at(peak %||% list(), "peak", stats::median(s)); pk <- max(pk, 1e-9); (s / pk) * exp(1 - s / pk) },
    delayed = { h <- at(half_time %||% list(), "half_time", span / 2); w <- span / 10; 1 / (1 + exp(-(s - h) / max(w, 1e-9))) },
    { note(sprintf("time_course shape '%s' unknown: monotone", shape)); s })
  f <- f - f[which.min(pos)]
  r <- diff(range(f))
  stats::setNames(if (r > 0) f / r else f, names(pos))
}

# Numeric values of a factor's levels for a trend (spec factor values, else numbers in the labels), on the
# factor's scale: for scale = log, log_base(x + log_offset) (Biodes 0.6.2, agreed with Bioana: base 10, offset
# the smallest positive level when a level is 0; the same predictor Bioana fits).
.suite_trend_x <- function(v, frame, p, spec, note) {
  col <- if (v %in% names(p$col)) p$col[[v]] else v
  f <- Filter(function(z) identical(z$name, col), spec$factors %||% list())
  f <- if (length(f)) f[[1]] else list()
  x <- frame[[v]]
  if (is.numeric(x)) { lev <- NULL; num <- x } else {
    lev <- levels(x)
    vals <- suppressWarnings(as.numeric(unlist(f$values))); sl <- as.character(unlist(f$levels))
    num <- if (length(vals) && length(vals) == length(sl) && all(lev %in% sl)) vals[match(lev, sl)] else
      suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+(\\.[0-9]+)?).*$", "\\1", lev)))
    if (anyNA(num) || length(unique(num)) < 2) return(NULL)
    names(num) <- lev
  }
  raw <- num
  scale <- f$scale %||% "linear"
  base <- offset <- NULL
  if (identical(scale, "log")) {
    pos <- num[num > 0]
    if (!length(pos)) return(NULL)
    base <- as.numeric(f$log_base %||% 10)
    offset <- as.numeric(f$log_offset %||% if (any(num <= 0)) min(pos) else 0)
    if (is.null(f$log_offset) && offset > 0) note(sprintf("trend on log scale: log_offset not given; %s (the smallest positive level)", format(offset)))
    if (any(num + offset <= 0)) return(NULL)
    num <- log(num + offset, base)
  }
  list(x = num, raw = raw, scale = scale, unit = f$unit, log_base = base, log_offset = offset)
}

# goal_truth block for sim_truth.json, built after planting
.suite_goal_truth <- function(g, tr, frame, p, extra = list()) {
  lev_diff <- function(e) { e <- unlist(e); as.list(e - e[1]) }
  out <- list(goal = g$goal, terms = as.list(g$terms), variables = as.list(g$vars), test = g$test,
              definition = g$definition %||% NULL, scale = g$scale, size = g$size,
              size_source = if (g$size_stated) "analysis.goal_effect" else "treatment effect (power.minimum_detectable_difference or default)",
              feature_size = g$feature_size)
  out <- c(out, extra)
  out$planted <- g$planted
  out$note <- g$note
  out[!vapply(out, is.null, logical(1))]
}

# The goal-specific part of goal_truth: true_effect, what it is measured on, and model_terms (the model terms
# that test it; the power sweep chooses replicates for these).
.suite_goal_extra <- function(g, truth, frame, p, x) {
  lev_diff <- function(e) { e <- unlist(e); e - e[1] }
  status <- "planted"; ex <- list()
  outcome_goal <- g$planted == "outcome"
  if (g$planted == "none") {
    status <- "not_simulated"
    ex$reason <- "Biodes registry simulated: false"
  } else if (g$goal == "treatment_effect" || (g$goal == "trend" && is.null(x$tr_var))) {
    tv <- intersect(x$treat_vars, names(truth$treatment_effects))
    d <- lapply(truth$treatment_effects[tv], function(e) as.list(lev_diff(e)))
    ex$terms <- as.list(tv); ex$variables <- as.list(tv)
    ex$level_differences <- d
    ex$true_effect <- if (length(tv)) max(abs(unlist(d))) else 0
    ex$effect_scale <- "outcome: each level minus the reference (first) level"
    ex$model_terms <- as.list(tv)
  } else if (g$goal == "interaction") {
    a <- g$vars[1]; b <- g$vars[2]
    cells <- unlist(truth$interaction_effects[[x$goal_label]]$cells)
    la <- levels(frame[[a]]); lb <- levels(frame[[b]])
    # cell keys follow the model term's variable order
    ord <- all.vars(str2lang(x$goal_label))
    cell <- function(i, j) { k <- if (identical(ord[1], a)) paste(i, j, sep = ":") else paste(j, i, sep = ":"); if (k %in% names(cells)) cells[[k]] else NA_real_ }
    dod <- list()
    for (i in la[-1]) for (j in lb[-1]) dod[[length(dod) + 1]] <- list(level_a = i, level_b = j,
      value = cell(i, j) - cell(i, lb[1]) - cell(la[1], j) + cell(la[1], lb[1]))
    ex$reference = list(a = la[1], b = lb[1])
    ex$difference_of_differences <- dod
    ex$true_effect <- dod[[length(dod)]]$value
    ex$effect_scale <- sprintf("outcome: (%s %s - %s %s) - (%s %s - %s %s), last levels vs references", a, la[length(la)], a, la[1], b, lb[length(lb)], b, lb[1])
    ex$pattern <- truth$interaction_effects[[x$goal_label]]$pattern
    if (identical(ex$pattern, "phenotype_model_gxe")) ex$note_gxe <- "genotype x environment from the phenotype model (gxe_share); not the goal size"
    ex$model_terms <- list(x$goal_label)
  } else if (g$goal == "time_course" && !is.null(x$tc_var)) {
    ex$time_pattern <- as.list(x$tc_pattern)
    ex$true_effect <- diff(range(x$tc_pattern))
    ex$effect_scale <- "outcome: largest change between time points"
    ex$shape <- g$shape %||% "monotone"
    ex$peak <- g$peak; ex$half_time <- g$half_time
    ex$scope <- if (!is.null(x$tc_group)) "within_group" else "all_units"
    if (!is.null(x$tc_group)) {
      ex$group <- x$tc_group
      ex$groups_changing <- as.list(x$tc_levels)
      ex$group_reference <- levels(frame[[x$tc_group]])[1]
      ex$model_terms <- list(x$goal_label)
    } else ex$model_terms <- list(x$tc_var)
  } else if (g$goal == "trend") {
    tr <- x$trend
    ex$slope <- tr$slope
    ex$true_effect <- tr$slope
    ex$effect_scale <- if (identical(tr$scale, "log")) sprintf("outcome per log%s unit of (%s + %s)", format(tr$log_base), tr$unit %||% "x", format(tr$log_offset)) else
      sprintf("outcome per %s", tr$unit %||% "unit")
    ex$values <- as.list(tr$raw); ex$x <- as.list(tr$x); ex$factor_scale <- tr$scale
    ex$log_base <- tr$log_base; ex$log_offset <- tr$log_offset
    # Biodes 0.6.3 per_unit: what a stated size is per ("kg N/ha", "log10 uM")
    if (!is.null(g$per_unit)) { ex$per_unit <- g$per_unit; ex$effect_scale <- sprintf("outcome per %s%s", g$per_unit,
      if (identical(tr$scale, "log")) sprintf(" (log%s of value + %s)", format(tr$log_base), format(tr$log_offset)) else "") }
    ex$change_over_range <- tr$slope * diff(range(tr$x))
    ex$model_terms <- list(x$tr_var)
  } else if (g$goal == "survival") {
    m <- x$fam$model
    tv <- intersect(x$treat_vars, names(truth$treatment_effects))
    if (identical(m$family, "lognormal_aft")) {
      ex$time_ratios <- lapply(truth$treatment_effects[tv], function(e) as.list(exp(m$time_sd * lev_diff(e) / m$scale$sd_outcome)))
      ex$true_effect <- if (length(tv)) unlist(ex$time_ratios[[1]])[length(ex$time_ratios[[1]])] else 1
      ex$effect_scale <- "time ratio (accelerated failure time), last level vs reference; hazards are not proportional"
    } else {
      status <- "not_planted"; ex$reason <- "primary_outcome.type is not survival"
    }
    ex$model_terms <- as.list(tv)
  } else if (g$goal == "variance_components") {
    ex$true_effect <- list(residual_sd = truth$residual_sd, random_sd = truth$random_sd,
                           within_unit_sd = truth$within_unit_sd$sd, heritability = truth$genetics$heritability)
    ex$effect_scale <- "standard deviations on the outcome scale"
  } else if (g$goal == "ranking") {
    v <- if (!is.null(x$genetics)) unlist(x$genetics$g) else unlist(truth$treatment_effects[[x$treat_vars[1]]])
    if (length(v)) {
      k <- max(1L, round(0.1 * length(v)))
      ex$entries <- as.list(v); ex$top_fraction <- 0.1
      ex$true_effect <- mean(sort(v, decreasing = TRUE)[seq_len(k)]) - mean(v)
      ex$effect_scale <- sprintf("outcome: mean of the top %d of %d entries minus the mean of all", k, length(v))
    } else { status <- "not_planted"; ex$reason <- "no entries" }
  } else if (g$planted == "genetics") {
    if (is.null(x$genetics)) { status <- "not_simulated"; ex$reason <- "no genotype term in analysis.model: no lines to give genotypes to" }
    else {
      ex$truth_files <- list(qtn_truth = "genetic_truth/qtn_truth.tsv", genetic_values = "genetic_truth/genetic_values.tsv", genotypes = "genotypes.vcf")
      ex$true_effect <- list(n_qtn = nrow(x$genetics$qtn), heritability = truth$genetics$heritability)
      ex$effect_scale <- "QTN effects per allele on the outcome scale (qtn_truth.tsv)"
    }
  } else if (g$planted == "assay") {
    if (!"assay_microbiome" %in% x$assay_routes) { status <- "not_simulated"; ex$reason <- "no microbiome assay simulated" }
    else {
      ex$truth_files <- list(features = "truth/microbiome.tsv")
      ex$effect_scale <- "log fold change of shifted taxa between condition levels (truth/microbiome.tsv)"
    }
  }
  if (outcome_goal && !isTRUE(x$has_outcome) && length(x$features) && !is.null(ex$true_effect)) {
    # assay-only design: the goal is in the features; the outcome-scale values are kept apart so they are not read as the data's
    ex$outcome_truth <- list(true_effect = ex$true_effect, effect_scale = ex$effect_scale, note = "the primary outcome is not written for this design")
    ex$true_effect <- "per feature: features[].file, column effect"
    ex$effect_scale <- unique(vapply(x$features, `[[`, "", "effect_scale"))[1]
  }
  if (outcome_goal && !isTRUE(x$has_outcome)) ex$outcome_note <- if (length(x$features)) "the primary outcome is not written (assay-only design): the goal is in the count assays listed under features" else
    "the primary outcome is not written (assay-only design) and no count assay took the goal: the goal truth here is not observable in the data"
  if (length(x$features)) ex$features <- x$features
  rest <- setdiff(x$assay_routes, paste0("assay_", names(x$features)))
  if (outcome_goal && length(rest) && g$goal %in% c("interaction", "time_course", "trend"))
    ex$features_note <- sprintf("%s: features carry condition-level differences only; this goal is planted on the primary outcome%s",
      paste(sub("^assay_", "", rest), collapse = ", "), if (length(x$features)) " and the count assays listed under features" else "")
  base <- .suite_goal_truth(g, truth, frame, p)
  base[names(ex)] <- ex
  if (length(base$variables)) base$columns <- as.list(vapply(unlist(base$variables), function(v) if (v %in% names(p$col)) as.character(p$col[[v]]) else v, ""))
  base$status <- status
  base[!vapply(base, is.null, logical(1))]
}

# Per-unit pattern of the goal (one value per frame row), 0 at the reference and 1 at the full goal effect:
# interaction = linear-by-linear score of the two terms (difference of differences of the extremes is 1),
# time_course = the time shape (largest change 1; only non-reference groups when within_group),
# trend = position over the numeric factor's range. NULL for goals with no per-unit pattern.
.suite_goal_pattern <- function(gl, frame, p, spec, tc_var = NULL, tc_group = NULL, tc_pattern = NULL, tr_var = NULL, trend = NULL, tc_levels = NULL) {
  score <- function(v) {
    x <- frame[[v]]
    tpos <- .suite_time_positions(v, frame, p, spec)
    if (!is.null(tpos)) return(unname(((tpos - min(tpos)) / diff(range(tpos)))[as.character(x)]))
    if (is.numeric(x)) { r <- range(x, na.rm = TRUE); return(if (diff(r) > 0) (x - r[1]) / diff(r) else 0 * x) }
    if (nlevels(x) > 1) (as.integer(x) - 1) / (nlevels(x) - 1) else 0 * as.integer(x)
  }
  if (gl$goal == "interaction" && length(gl$vars) == 2) return(score(gl$vars[1]) * score(gl$vars[2]))
  if (!is.null(tc_var)) {
    s <- tc_pattern / diff(range(tc_pattern))
    z <- unname(s[as.character(frame[[tc_var]])])
    if (!is.null(tc_group)) z <- z * (frame[[tc_group]] %in% tc_levels)
    return(z)
  }
  if (!is.null(tr_var)) {
    x <- if (is.numeric(frame[[tr_var]])) trend$x else unname(trend$x[as.character(frame[[tr_var]])])
    return((x - min(trend$x)) / diff(range(trend$x)))
  }
  NULL
}

# Plant the goal in a count matrix (features x units, first column the feature id): a share of features gets
# log2 fold change lfc_f * pattern(unit), applied by binomial thinning (each feature scaled so its largest
# factor is 1), so the count noise stays that of the engine. Features already given a condition difference by
# the engine are avoided when possible. Returns the truth table (feature_id, term, effect, change).
.suite_goal_features <- function(file, pattern, sample_ids, avoid, a, term, per_unit = NULL, seed = 1) {
  x <- utils::read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- as.character(x[[1]]); m <- as.matrix(x[, -1, drop = FALSE])
  s <- match(colnames(m), sample_ids)
  if (anyNA(s) || !nrow(m)) return(NULL)
  set.seed(seed + 7919)
  pool <- setdiff(seq_along(ids), which(ids %in% avoid))
  if (length(pool) < 10) pool <- seq_along(ids)
  k <- max(1L, round(a$fraction * length(ids)))
  pick <- pool[sample.int(length(pool), min(k, length(pool)))]
  lfc <- sample(c(-1, 1), length(pick), replace = TRUE) * (a$log2fc + abs(stats::rnorm(length(pick), 0, a$log2fc_sd)))
  z <- pattern[s]
  for (j in seq_along(pick)) {
    r <- 2^(lfc[j] * z); r <- r / max(r)
    m[pick[j], ] <- stats::rbinom(ncol(m), round(m[pick[j], ]), r)
  }
  x[, -1] <- m
  utils::write.table(x, file, sep = "\t", row.names = FALSE, quote = FALSE)
  data.frame(feature_id = ids[pick], term = term, effect = if (is.null(per_unit)) lfc else lfc / per_unit, change = lfc, stringsAsFactors = FALSE)
}
