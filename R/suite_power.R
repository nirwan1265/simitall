#' Power over replication for any Biodes design
#'
#' Takes a design folder (or the result of [simitall_from_question()]) and,
#' for each value of `design.reps`, rebuilds the layout from the study details
#' (Biodes `studyforge::layout_from_spec()`), simulates `n_sim` experiments
#' ([simitall_from_biodes()]) and fits the design's model (`bioana::bioana_analyze()`; needs Bioana).
#' Power is the share of experiments in which a model term's ANOVA test is
#' significant. Nothing here depends on the question: the design, factors,
#' model, effect size and residual SD all come from the folder's
#' `study_spec.json`, and defaults from `inst/suite/power_defaults.json`.
#'
#' Interactions between treatment factors are simulated with
#' `interaction_effect` (default from the spec or `simulation_defaults.json`);
#' terms with no simulated effect (e.g. `interaction_effect = 0`) show the
#' false-positive rate instead of power; the `simulated_effect` column says which.
#'
#' @param x A design folder written by [simitall_from_question()] or Biodes
#'   `write_design_package()`, or a `simitall_suite_run`.
#' @param reps Values of `design.reps` to try.
#' @param n_sim Simulated experiments per value.
#' @param alpha Significance level.
#' @param effect,residual_sd,interaction_effect Override the spec's power values.
#' @param seed First simulation seed.
#' @param show Show the power figure in the plot window (default: in
#'   interactive sessions); `plot(ps)` shows it again.
#' @param max_minutes Stop before running when the sweep is estimated to take
#'   longer (one fit is timed first).
#' @param allow_fallback Passed to `bioana::bioana_analyze()`: `FALSE` (Bioana's
#'   default) fits each data type only with its citable package (DESeq2, ...),
#'   the same one the real data will get; `TRUE` lets Bioana use its own
#'   method when that package is not installed.
#' @return Invisibly, a data frame (reps, n_units, term, factor, role,
#'   simulated_effect (NA for blocking terms), goal (the spec's primary analysis goal),
#'   goal_term (TRUE for the terms that test it), power, mc_se, n_fits); also written to `power_sweep.tsv`,
#'   `power_sweep.json` and `figures/power.png` in the folder.
#' @export
simitall_power_sweep <- function(x, reps = NULL, n_sim = NULL, alpha = NULL, effect = NULL, residual_sd = NULL, seed = 1,
                                 max_minutes = NULL, interaction_effect = NULL, show = interactive(), allow_fallback = FALSE) {
  dir <- if (inherits(x, "simitall_suite_run")) x$dir else x
  if (is.null(dir) || !file.exists(file.path(dir, "study_spec.json"))) {
    stop("simitall_power_sweep() needs a design folder with study_spec.json (from simitall_from_question() or Biodes write_design_package()).", call. = FALSE)
  }
  if (!requireNamespace("bioana", quietly = TRUE)) {
    stop("simitall_power_sweep() fits each simulated experiment with Bioana: remotes::install_github('nirwan1265/bioana')", call. = FALSE)
  }
  if (!requireNamespace("studyforge", quietly = TRUE)) {
    stop("simitall_power_sweep() lays out each replicate count with Biodes (R package 'studyforge'); load or install it first.", call. = FALSE)
  }
  d <- .suite_rules("power_defaults")
  reps <- reps %||% unlist(d$reps$value)
  n_sim <- n_sim %||% d$n_sim$value
  alpha <- alpha %||% d$alpha$value
  max_minutes <- max_minutes %||% d$max_minutes$value
  checked <- FALSE
  spec <- jsonlite::read_json(file.path(dir, "study_spec.json"), simplifyVector = FALSE)
  design <- spec$design$chosen %||% stop("study_spec.json has no design.chosen.", call. = FALSE)
  model <- spec$analysis$model %||% stop("study_spec.json has no analysis.model.", call. = FALSE)
  # count assays (counts.tsv): power is also the share of genes with a true >= 2-fold change found (Bioana DE, padj < alpha)
  counts <- .suite_is_count_assay(spec)
  n_genes <- if (counts) d$n_genes_counts$value %||% 300 else d$n_genes$value
  meaning <- spec$layout$factor_map %||% spec$analysis$terms %||% list()
  root <- tempfile("simitall_power_")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  rows <- list()
  for (r in reps) {
    s <- spec
    s$design$reps <- r
    book <- tryCatch(studyforge::layout_from_spec(s, seed = spec$design$randomization_seed %||% 42), error = function(e) e)
    if (inherits(book, "error")) {
      warning(sprintf("reps = %d: %s", r, .suite_layout_missing(conditionMessage(book))), call. = FALSE)
      next
    }
    base <- file.path(root, paste0("reps_", r))
    studyforge::write_design_package(book, s, base)
    sim <- file.path(base, "sim")
    dir.create(sim)
    file.copy(file.path(base, c("study_spec.json", "field_layout.csv")), sim)
    sig <- list(); found <- list(); found_goal <- list()
    n_units <- NA
    if (!exists("refused", inherits = FALSE)) refused <- character()
    for (k in seq_len(n_sim)) {
      t0 <- Sys.time()
      res <- tryCatch({
        suppressMessages(suppressWarnings(simitall_from_biodes(base, outdir = sim, seed = seed + k - 1, effect = effect,
                                                               residual_sd = residual_sd, n_genes = n_genes,
                                                               interaction_effect = interaction_effect)))
        suppressMessages(suppressWarnings(bioana::bioana_analyze(sim, quiet = TRUE, allow_fallback = allow_fallback)))
      }, error = function(e) NULL)
      # an analysis Bioana refused (e.g. no citable engine installed): say why once, the sweep goes on without it
      for (id in names(res$summary$not_analysed)) {
        why <- res$summary$not_analysed[[id]]
        key <- paste(id, why$note)
        if (!key %in% refused) { refused <- c(refused, key); warning(sprintf("Bioana did not analyse '%s' (%s): %s", id, why$status, why$note), call. = FALSE) }
      }
      if (!checked) {
        # one fit timed: scale by the number of units at the largest reps value
        checked <- TRUE
        secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        est <- secs * n_sim * sum(reps / r) / 60
        if (est > max_minutes) {
          stop(sprintf("The sweep would take about %.0f minutes (one simulated fit of %d rows took %.1f s; %d experiments x %d reps values). Limit: max_minutes = %s (power_defaults.json). Use fewer reps values or a smaller n_sim, or raise max_minutes.",
                       est, nrow(book), secs, n_sim, length(reps), max_minutes), call. = FALSE)
        }
      }
      if (is.null(res)) next
      an <- res$summary$analyses
      n_units <- an$outcome$n_units_fitted %||% unlist(lapply(an, `[[`, "n_samples"))[1] %||% n_units
      # Bioana's recovery of the goal's features (manifest truth.goal; summary goal_recovery$tpr), when it reports it
      for (a in Filter(function(z) !is.null(z$goal_recovery$tpr), an)) {
        k <- paste0("goal features: ", a$goal_recovery$goal %||% "goal", " (", a$data_type, ")")
        found_goal[[k]] <- c(found_goal[[k]], as.numeric(a$goal_recovery$tpr))
      }
      for (a in Filter(function(z) identical(z$data_type, "feature_counts"), an)) {
        for (cr in a$recovery %||% list()) {
          if (!is.null(cr$contrast) && !is.null(cr$tpr_large) && !is.na(cr$tpr_large)) found[[cr$contrast]] <- c(found[[cr$contrast]], cr$tpr_large)
        }
      }
      for (a in an$outcome$anova %||% list()) {
        pv <- a[["Pr(>F)"]]
        if (is.null(a$term) || a$term == "Residuals" || is.null(pv) || is.na(pv)) next
        sig[[a$term]] <- c(sig[[a$term]], pv < alpha)
      }
    }
    truth <- tryCatch(jsonlite::read_json(file.path(sim, "sim_truth.json")), error = function(e) list())
    # the primary goal's model terms (sim_truth.json goal_truth.model_terms): the replicate choice is made for these
    goal_sets <- lapply(unlist(truth$goal_truth$model_terms), function(z) sort(all.vars(str2lang(z))))
    is_goal <- function(term) any(vapply(goal_sets, function(z) identical(z, sort(all.vars(str2lang(term)))), logical(1)))
    p <- .suite_parse_model(model, s, book)
    for (term in names(sig)) {
      vars <- all.vars(str2lang(term))
      role <- if (all(p$role[vars] == "treatment")) "treatment" else "nuisance"
      eff <- if (role != "treatment") NA else if (grepl(":", term)) isTRUE(abs(as.numeric(truth$interaction_effects[[term]]$size %||% 0)) > 0) else
        term %in% names(truth$treatment_effects) || term %in% names(truth$numeric_effects)
      rows[[length(rows) + 1]] <- data.frame(
        reps = r, n_units = n_units, term = term,
        factor = paste(vapply(vars, function(v) as.character(unlist(meaning[[v]])[1] %||% v), ""), collapse = " x "),
        role = role, simulated_effect = eff, goal = truth$goal_truth$goal %||% NA_character_, goal_term = is_goal(term), power = mean(sig[[term]]),
        mc_se = sqrt(mean(sig[[term]]) * (1 - mean(sig[[term]])) / length(sig[[term]])), n_fits = length(sig[[term]]),
        stringsAsFactors = FALSE)
    }
    for (k in names(found_goal)) {
      f <- found_goal[[k]]
      rows[[length(rows) + 1]] <- data.frame(
        reps = r, n_units = n_units, term = k, factor = "share of the goal's features found (truth.goal in data_manifest.json)",
        role = "treatment", simulated_effect = TRUE, goal = truth$goal_truth$goal %||% NA_character_, goal_term = TRUE, power = mean(f),
        mc_se = if (length(f) > 1) stats::sd(f) / sqrt(length(f)) else NA_real_, n_fits = length(f), stringsAsFactors = FALSE)
    }
    for (ct in names(found)) {
      f <- found[[ct]]
      rows[[length(rows) + 1]] <- data.frame(
        reps = r, n_units = n_units, term = paste0("DE genes: ", ct), factor = "share of genes with a true 2-fold change found (counts.tsv)",
        role = "treatment", simulated_effect = TRUE, goal = truth$goal_truth$goal %||% NA_character_,
        goal_term = is.null(truth$goal_truth$features) && (truth$goal_truth$goal %||% "treatment_effect") == "treatment_effect", power = mean(f),
        mc_se = if (length(f) > 1) stats::sd(f) / sqrt(length(f)) else NA_real_, n_fits = length(f), stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) stop("No reps value could be laid out and fitted; see the warnings.", call. = FALSE)
  out <- do.call(rbind, rows)
  target <- d$target_power$value
  real <- out$simulated_effect %in% TRUE
  enough <- lapply(split(out[real, ], out$term[real]), function(z) {
    hit <- z$reps[z$power >= target]
    if (length(hit)) min(hit) else NA
  })
  utils::write.table(out, file.path(dir, "power_sweep.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  .suite_write_json(list(design = design, model = model, alpha = alpha, n_sim = n_sim, reps = reps,
                         effect = effect %||% spec$power$minimum_detectable_difference,
                         residual_sd = residual_sd %||% spec$power$residual_sd,
                         interaction_effect = interaction_effect %||% spec$power$interaction_effect %||% "simulation_defaults.json interaction_effect_ratio x effect",
                         target_power = target, smallest_reps_reaching_target = enough,
                         note = "Terms with simulated_effect = false show the false-positive rate, not power."),
                    file.path(dir, "power_sweep.json"))
  if (requireNamespace("ggplot2", quietly = TRUE)) tryCatch(suppressMessages(biopaint_biodes(dir, which = "power", show = show)), error = function(e) NULL)
  structure(out, class = c("simitall_power_sweep", "data.frame"), target = target, enough = enough, alpha = alpha, n_sim = n_sim, dir = dir)
}

#' @export
print.simitall_power_sweep <- function(x, ...) {
  y <- as.data.frame(unclass(x))
  wide <- stats::reshape(y[, c("reps", "n_units", "term", "power")], idvar = c("reps", "n_units"), timevar = "term", direction = "wide")
  names(wide) <- sub("^power\\.", "", names(wide))
  labels <- unique(y[, c("term", "factor", "simulated_effect")])
  cat(sprintf("Power at alpha = %s over %d simulated experiments per row (Monte Carlo SE up to %.2f)\n", attr(x, "alpha"), attr(x, "n_sim"),
              max(y$mc_se[y$simulated_effect %in% TRUE], 0)))
  for (i in seq_len(nrow(labels))) {
    cat(sprintf("  %-22s %s%s\n", labels$term[i], labels$factor[i],
                if (is.na(labels$simulated_effect[i])) "  (blocking/nuisance term)" else if (!labels$simulated_effect[i]) "  (no simulated effect: false-positive rate)" else ""))
  }
  print(format(wide, digits = 2), row.names = FALSE)
  en <- attr(x, "enough")
  for (t in names(en)) {
    cat(sprintf("  %s: %s\n", t, if (is.na(en[[t]])) sprintf("power %s not reached in this range", attr(x, "target")) else
      sprintf("power >= %s from reps = %d", attr(x, "target"), en[[t]])))
  }
  invisible(x)
}
