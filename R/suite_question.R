#' From a plain-language question to simulated and analysed data
#'
#' One call across the suite: Biodes (package `studyforge`,
#' github.com/nirwan1265/BioDesignGPT) reads the question, chooses the design
#' and lists what it still needs to know; Simitall builds the layout from those
#' study details (Biodes `studyforge::layout_from_spec()`, every factor kept)
#' and writes the design folder, then simulates data
#' under that exact design ([simitall_from_biodes()]); Bioana fits the model the
#' design calls for and checks it recovers the truth (`bioana::bioana_analyze()`, when Bioana is installed).
#' When the study details lack something the layout needs (a factor's number of
#' levels, the number of environments), nothing is simulated and the missing
#' details are listed.
#'
#' Biodes' open questions are shown before anything is simulated. Answer them
#' with `answers` (Biodes field paths, as in `studyforge::studyforge_fix()`),
#' e.g. `answers = list(factors = "genotype:6, water:3", design.reps = 4, subsamples_per_unit = 3,
#' power.minimum_detectable_difference = 2, power.residual_sd = 3)`. Unanswered
#' values fall back to the defaults Simitall declares in `sim_truth.json`.
#'
#' @param question The experiment, in plain language.
#' @param answers Named list of answers / corrections for Biodes (field paths).
#' @param outdir Design folder to write. Default: `~/biodes_runs/<question words>`.
#' @param engine Biodes reader: "rules" (default), "hybrid", "llm" or "auto".
#' @param simulate,analyze Run Simitall and Bioana after writing the design.
#' @param confirm Ask before simulating when Biodes still has open questions
#'   (default: in interactive sessions).
#' @param seed Simulation seed (default: the design's randomization seed).
#' @param params Simulator parameter overrides passed to [simitall_from_biodes()],
#'   e.g. `list(suite_genetics = list(heritability = 0.7))`.
#' @param force Simulate even when the layout does not match the plan (a study
#'   factor missing from the model, level counts that differ). Off by default:
#'   that would simulate a different experiment from the one asked.
#' @param design Catalog design to lay out instead of Biodes' first candidate.
#' @param layout_seed Randomization seed for the layout.
#' @param figures Draw the figures ([biopaint_biodes()]) into `<outdir>/figures`
#'   and, in interactive sessions, show them in the plot window (`plot(r)`
#'   shows them again).
#' @param use_suggestions Fill each item the study details leave open with its
#'   first suggestion (Biodes' traits and units, `inst/suite/suggestions.tsv`)
#'   and run. Every filled value is printed and kept in `$filled_from_suggestions`.
#'   Your own `answers` always win.
#' @param allow_fallback Passed to `bioana::bioana_analyze()` and the power
#'   sweep: `FALSE` (default) analyses each data type only with its citable
#'   package; `TRUE` allows Bioana's own method when that package is missing.
#' @param layout_defaults When a factor has no level count (or an augmented
#'   design has no checks), fill it with Biodes' stated default
#'   (`studyforge::fill_layout_defaults()`) and list each as an assumption.
#'   Default: the same as `use_suggestions`.
#' @param reps_range,power_n_sim For a "how many per group?" question (no
#'   `design.reps`, or `power.solve_for = "reps"`): the replicate counts to try
#'   (default: the spec's `power.reps_range`, else `power_defaults.json`) and the
#'   simulated experiments per count (default `n_sim_question`). The smallest
#'   count giving every treatment main effect (and, for count assays, the true
#'   genes with a 2-fold change) the target power is used; the sweep is kept in
#'   `power_sweep.tsv` and the choice in `power_choice.json`.
#' @return Invisibly, a `simitall_suite_run`: the Biodes answer (`$answer`),
#'   the design folder (`$dir`), open questions, and the simulation and
#'   analysis results (`$analysis$results`, `$analysis$recovery`). Printing it
#'   shows a short summary.
#' @examples
#' \dontrun{
#' devtools::load_all("~/Documents/Github/BioDesignGPT")   # Biodes
#' r <- simitall_from_question(
#'   "Greenhouse drought on 6 sorghum genotypes, 4 benches as blocks, pot is the unit, chlorophyll content",
#'   answers = list(power.minimum_detectable_difference = 2, power.residual_sd = 3))
#' }
#' @export
simitall_from_question <- function(question, answers = NULL, outdir = NULL,
                                   engine = c("rules", "hybrid", "llm", "auto"),
                                   simulate = TRUE, analyze = TRUE, confirm = interactive(), seed = NULL,
                                   force = FALSE, design = NULL, layout_seed = 42, use_suggestions = FALSE,
                                   figures = TRUE, params = NULL, reps_range = NULL, power_n_sim = NULL,
                                   layout_defaults = use_suggestions, allow_fallback = FALSE) {
  engine <- match.arg(engine)
  if (!requireNamespace("studyforge", quietly = TRUE)) {
    stop("simitall_from_question() needs Biodes (R package 'studyforge'). Load it with devtools::load_all(\"~/Documents/Github/BioDesignGPT\") or install it with remotes::install_github(\"nirwan1265/BioDesignGPT\").",
         call. = FALSE)
  }
  ans <- studyforge::studyforge_ask(question, engine = engine, fix = answers)
  sug <- .suite_suggestions(ans)
  suggested <- list()
  if (isTRUE(use_suggestions) && !is.null(sug)) {
    # fill only what the user did not answer; every filled value is reported
    suggested <- .suite_suggested_answers(sug, ans$spec)
    suggested <- suggested[setdiff(names(suggested), names(answers))]
    if (length(suggested)) ans <- studyforge::studyforge_ask(question, engine = engine, fix = c(answers, suggested))
    sug <- .suite_suggestions(ans)
  }
  outdir <- path.expand(outdir %||% file.path("~", "biodes_runs", .suite_slug(question)))

  # the layout comes from Biodes (studyforge::layout_from_spec): every factor kept as its own column
  spec_in <- ans$spec
  if (!is.null(design)) spec_in$design$chosen <- design
  # a factor with no level count (or augmented checks) the user left open: Biodes' stated defaults, each an assumption
  defaults_filled <- character()
  if (isTRUE(layout_defaults) && "fill_layout_defaults" %in% getNamespaceExports("studyforge")) {
    filled <- studyforge::fill_layout_defaults(spec_in)
    defaults_filled <- setdiff(names(filled$assumptions), names(spec_in$assumptions))
    if (length(defaults_filled)) {
      spec_in <- filled; ans$spec <- filled
      val <- function(k) {
        if (!startsWith(k, "factors.")) return(paste(format(unlist(.suite_get(spec_in, k))), collapse = ","))
        nm <- strsplit(k, ".", fixed = TRUE)[[1]][2]
        format(Filter(function(f) identical(f$name, nm), spec_in$factors)[[1]]$n_levels)
      }
      cat("Layout defaults (assumed, Biodes rules.tsv layout_default; change with answers = list(...)):\n",
          paste0("  - ", defaults_filled, " = ", vapply(defaults_filled, val, ""), "\n"), sep = "")
    }
  }
  # "how many per group?": no design.reps (or power.solve_for = reps) -> choose it with a power sweep
  power <- NULL
  if (.suite_reps_open(spec_in) && isTRUE(simulate)) {
    power <- .suite_choose_reps(spec_in, outdir, layout_seed, reps_range, power_n_sim, allow_fallback)
    if (!is.null(power$reps)) {
      spec_in$design$reps <- power$reps
      ans$spec$design$reps <- power$reps
    }
  }
  book <- tryCatch(studyforge::layout_from_spec(spec_in, seed = layout_seed), error = function(e) e)
  if (inherits(book, "error")) {
    needs <- .suite_layout_missing(conditionMessage(book))
    .suite_show_plan(ans, NULL, design)
    .suite_show_suggestions(ans, sug)
    cat("Biodes cannot lay out this study yet:\n", paste0("  - ", needs, "\n"), sep = "")
    message("Not simulated. Give the missing details with answers = list(...) and run again.")
    return(invisible(structure(list(answer = ans, dir = NULL, design = design, suggestions = sug, filled_from_suggestions = suggested,
                                    layout_needs = needs, open_questions = needs, layout_issues = character(), simulation = NULL, analysis = NULL),
                               class = "simitall_suite_run")))
  }
  design <- attr(book, "design") %||% design
  studyforge::write_design_package(book, spec_in, outdir, ans$audit)
  open <- .suite_show_plan(ans, outdir, design)
  if (length(suggested)) cat("Filled from suggestions (assumptions, change with answers = list(...)):\n",
                             paste0("  - ", names(suggested), " = ", vapply(suggested, function(v) paste(format(v), collapse = ", "), ""), "\n"), sep = "")
  .suite_show_suggestions(ans, sug)
  if (length(attr(book, "layout_notes"))) cat("Layout notes (Biodes):\n", paste0("  - ", attr(book, "layout_notes"), "\n"), sep = "")
  issues <- .suite_layout_issues(ans, outdir, design)
  if (length(issues)) {
    cat("Layout does NOT match Biodes' plan (simulating it would test a different experiment):\n",
        paste0("  ! ", issues, "\n"), sep = "")
    open <- c(open, issues)
  }

  out <- structure(list(answer = ans, dir = outdir, design = design, suggestions = sug, filled_from_suggestions = suggested, open_questions = open, layout_issues = issues, simulation = NULL,
                        analysis = NULL, power = power, layout_defaults = defaults_filled), class = "simitall_suite_run")
  if (!is.null(power$sweep)) {
    # keep the sweep with the design it chose
    utils::write.table(as.data.frame(unclass(power$sweep)), file.path(outdir, "power_sweep.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    .suite_write_json(power[c("reps", "reason", "target_power", "n_sim", "reps_tried")], file.path(outdir, "power_choice.json"))
  }
  if (!simulate) return(invisible(out))
  if (length(issues) && !isTRUE(force)) {
    message("Not simulated: the layout does not match Biodes' plan (see above). This needs a fix in Biodes ",
            "(build the layout from the spec). To test the layout as written anyway, run again with force = TRUE; ",
            "or make the layout yourself with studyforge::make_design() + write_design_package() and call simitall_from_biodes().")
    return(invisible(out))
  }
  if (length(open) && isTRUE(confirm)) {
    go <- tolower(trimws(readline("Simulate now, using declared defaults for the open questions? [y/N] ")))
    if (!go %in% c("y", "yes")) {
      message("Stopped before simulating. Answer the questions with answers = list(<field> = <value>, ...) and run again; ",
              "the design folder is in ", outdir)
      return(invisible(out))
    }
  }
  out$simulation <- simitall_from_biodes(outdir, seed = seed, params = params)
  if (analyze) {
    if (requireNamespace("bioana", quietly = TRUE)) {
      out$analysis <- bioana::bioana_analyze(outdir, quiet = TRUE, allow_fallback = allow_fallback)
    } else {
      message("simitall: analysis skipped; install Bioana to analyse the data: remotes::install_github('nirwan1265/bioana')")
    }
  }
  if (isTRUE(figures) && requireNamespace("ggplot2", quietly = TRUE)) {
    out$figures <- tryCatch(biopaint_biodes(outdir), error = function(e) {
      message("biopaint: figures skipped (", conditionMessage(e), ")")
      NULL
    })
  }
  invisible(out)
}

.suite_slug <- function(x, n = 6) {
  w <- strsplit(gsub("[^a-z0-9 ]", " ", tolower(x)), "\\s+")[[1]]
  w <- w[nchar(w) > 2 & !w %in% c("the", "and", "for", "with", "per", "are", "from", "measure")]
  paste(utils::head(w, n), collapse = "_")
}

# Print what Biodes understood and what it still needs; return the open items.
.suite_show_plan <- function(ans, outdir, design = NULL) {
  s <- ans$spec
  spec <- if (!is.null(outdir)) jsonlite::read_json(file.path(outdir, "study_spec.json"), simplifyVector = FALSE) else
    list(design = list(chosen = design), analysis = list(model = tryCatch(studyforge::lookup_design(design)$model, error = function(e) NULL)))
  fl <- function(x) if (is.null(x) || !length(x)) "-" else paste(unlist(x), collapse = ", ")
  factors <- vapply(s$factors %||% list(), function(f) sprintf("%s (%s levels)", f$name %||% "?", f$n_levels %||% "?"), "")
  # fields the plan itself fills (design.candidates) are Biodes bookkeeping, not questions for the user
  missing <- setdiff(unlist(ans$validation$missing), "design.candidates")
  open <- unique(c(if (length(missing)) paste0("missing: ", missing), unlist(ans$questions)))
  cat("Biodes plan\n",
      "  organism:          ", fl(s$organism), "\n",
      "  system:            ", fl(s$system), "\n",
      "  experimental unit: ", fl(s$experimental_unit), if (!is.null(s$subsamples_per_unit)) sprintf(" (%s subsamples each)", s$subsamples_per_unit), "\n",
      "  factors:           ", fl(factors), "\n",
      "  blocking:          ", fl(s$blocking_factors), "\n",
      "  outcome:           ", fl(s$primary_outcome$name), "\n",
      "  design:            ", fl(spec$design$chosen), " (candidates: ", fl(ans$candidates), ")\n",
      "  model:             ", fl(spec$analysis$model), "\n",
      "  power:             effect ", fl(s$power$minimum_detectable_difference), ", residual SD ", fl(s$power$residual_sd), "\n",
      "  design folder:     ", outdir %||% "-", "\n", sep = "")
  if (length(open)) cat("Open questions (answer with answers = list(...)):\n", paste0("  - ", open, "\n"), sep = "")
  open
}

# Does the written layout match the plan? Flag a different design, study factors
# whose level counts no treatment column has, and factors the model leaves out.
.suite_layout_issues <- function(ans, outdir, design = NULL) {
  # Biodes' layout keeps every factor as its own column (layout_from_spec); a spec factor that is not a
  # column, or a factor column whose levels join several factors, means the data would test another experiment
  spec <- jsonlite::read_json(file.path(outdir, "study_spec.json"), simplifyVector = FALSE)
  layout <- utils::read.csv(file.path(outdir, "field_layout.csv"), check.names = FALSE, colClasses = "character")
  issues <- character()
  fac <- Filter(nzchar, vapply(spec$factors %||% list(), function(f) f$name %||% "", ""))
  for (f in fac) {
    if (!f %in% names(layout)) {
      issues <- c(issues, sprintf("factor '%s' is in the spec but not a column of the layout.", f))
    } else if (length(fac) > 1 && any(grepl(":", layout[[f]], fixed = TRUE))) {
      issues <- c(issues, sprintf("layout column '%s' holds joined levels of several factors (e.g. '%s'), not factor '%s' itself.", f, layout[[f]][1], f))
    }
  }
  issues
}

#' @export
print.simitall_suite_run <- function(x, ...) {
  s <- x$answer$spec
  cat("Biodes -> Simitall -> Bioana\n")
  cat("  question:  ", x$answer$question, "\n", sep = "")
  cat("  folder:    ", x$dir, "\n", sep = "")
  bs <- x$analysis$summary
  fs <- bs$analyses$outcome
  if (is.null(x$simulation)) {
    cat("  status:    not simulated", if (length(x$layout_needs)) " (the layout needs more study details)" else if (length(x$layout_issues)) " (layout does not match the plan)", "\n", sep = "")
  } else if (is.null(bs)) {
    cat("  status:    simulated; not analysed (install Bioana: remotes::install_github('nirwan1265/bioana'))\n")
  } else {
    if (!is.null(fs)) cat(sprintf("  design:    %s, %s; %d units fitted with %s\n", bs$design, bs$model, fs$n_units_fitted, fs$method))
    for (id in setdiff(names(bs$analyses), "outcome")) {
      a <- bs$analyses[[id]]; sig <- unlist(a$significant)
      cat(sprintf("  %-10s %d features tested, %s significant\n", paste0(id, ":"), a$n_features_tested %||% 0L, if (length(sig)) paste(sig, collapse = "/") else "0"))
    }
    for (a in fs$anova %||% list()) {
      if (!is.null(a$term) && a$term != "Residuals" && !is.null(a[["Pr(>F)"]])) {
        cat(sprintf("  %-10s F = %.2f on %s df, p = %.3g\n", paste0(a$term, ":"), a[["F value"]], a$Df, a[["Pr(>F)"]]))
      }
    }
    if (!is.null(fs$recovery)) {
      cat(sprintf("  recovery:  %d/%d treatment effects inside their %d%% CI (bias %.3g, RMSE %.3g)\n",
                  fs$recovery$covered, fs$recovery$n_terms, round(100 * fs$recovery$nominal), fs$recovery$bias, fs$recovery$rmse))
    }
    if (length(fs$warnings)) cat("  warnings:  ", paste(fs$warnings, collapse = "; "), "\n", sep = "")
  }
  if (!is.null(x$power$reason)) cat("  reps:      ", x$power$reason, "\n", sep = "")
  if (length(x$open_questions)) cat("  open questions: ", length(x$open_questions), " (see $open_questions)\n", sep = "")
  if (!is.null(bs)) cat("  files:     bioana/ (results per data set, bioana_summary.json) in the folder\n")
  if (length(x$figures)) cat("  figures:   ", paste(basename(x$figures), collapse = ", "), " in ", dirname(x$figures[[1]]), "\n", sep = "")
  invisible(x)
}

# "layout_from_spec: no level count for time (missing: factors.n_levels)." -> what the user must give
.suite_layout_missing <- function(msg) {
  field <- regmatches(msg, regexpr("missing: [^)]+", msg))
  what <- sub("^layout_from_spec: ", "", sub(" \\(missing: [^)]+\\)\\.?$", "", msg))
  out <- if (length(field)) paste0(what, " (", field, ")") else msg
  if (grepl("design.reps", msg, fixed = TRUE)) out <- paste0(out, "; or run with simulate = TRUE to choose it by a power sweep")
  out
}

# Spec factors by part: treatment factors in spec order, plus time, check, env,
# period (factor_roles.tsv) and the main-plot factor (handoff.json main_plot_applied_to).
.suite_spec_factors <- function(spec) {
  p <- system.file("suite", "factor_roles.tsv", package = "simitall")
  if (!nzchar(p)) p <- file.path("inst", "suite", "factor_roles.tsv")
  roles <- utils::read.delim(p, colClasses = "character", quote = "")
  main_applied <- unlist(.suite_rules("handoff")$main_plot_applied_to)
  out <- list(treatment = list())
  for (x in spec$factors %||% list()) {
    nm <- x$name %||% "factor"
    n <- suppressWarnings(as.integer(x$n_levels %||% if (length(x$levels)) length(x$levels) else NA))
    item <- list(name = nm, n = n, levels = unlist(x$levels))
    hit <- roles$role[vapply(roles$pattern, function(pt) grepl(pt, tolower(nm), perl = TRUE), logical(1))][1]
    if (!is.na(hit) && is.null(out[[hit]])) {
      out[[hit]] <- item
      next
    }
    if (is.null(out$main) && any(unlist(x$applied_to) %in% main_applied)) out$main <- item
    out$treatment[[length(out$treatment) + 1]] <- item
  }
  out
}

# design.reps still to be chosen: not given, or Biodes marks the question as solving for it
.suite_reps_open <- function(spec) {
  is.null(spec$design$reps) || identical(tolower(as.character(spec$power$solve_for %||% "")), "reps")
}

# Lay the design out at the smallest count, run simitall_power_sweep() over the range and pick the
# smallest count where every treatment main effect (and the DE genes of a count assay) reaches the target.
.suite_choose_reps <- function(spec, outdir, layout_seed, reps_range = NULL, n_sim = NULL, allow_fallback = FALSE) {
  d <- .suite_rules("power_defaults")
  reps <- .suite_reps_grid(reps_range %||% spec$power$reps_range %||% d$reps$value)
  n_sim <- n_sim %||% d$n_sim_question$value %||% 50
  target <- d$target_power$value
  s <- spec; s$design$reps <- min(reps)
  book <- tryCatch(studyforge::layout_from_spec(s, seed = layout_seed), error = function(e) e)
  if (inherits(book, "error")) return(NULL)   # the layout needs something else first; reported by the caller
  tmp <- file.path(tempfile("simitall_reps_"), "design")
  on.exit(unlink(dirname(tmp), recursive = TRUE), add = TRUE)
  suppressMessages(studyforge::write_design_package(book, s, tmp))
  cat(sprintf("No replicate count given: choosing it by power (reps %s, %d simulated experiments each)...\n", paste(reps, collapse = ", "), n_sim))
  ps <- tryCatch(suppressWarnings(simitall_power_sweep(tmp, reps = reps, n_sim = n_sim, show = FALSE, allow_fallback = allow_fallback)), error = function(e) e)
  if (inherits(ps, "error") && grepl("would take about", conditionMessage(ps))) {
    n_sim <- max(10L, floor(n_sim / 5))
    ps <- tryCatch(suppressWarnings(simitall_power_sweep(tmp, reps = reps, n_sim = n_sim, show = FALSE, allow_fallback = allow_fallback)), error = function(e) e)
  }
  if (inherits(ps, "error")) {
    r <- d$reps_if_sweep_fails$value %||% 3L
    msg <- sprintf("design.reps = %d assumed: the power sweep failed (%s)", r, conditionMessage(ps))
    message("simitall: ", msg)
    return(list(reps = as.integer(r), reason = msg, target_power = target, n_sim = n_sim, reps_tried = reps, sweep = NULL))
  }
  y <- as.data.frame(unclass(ps))
  main <- y[y$simulated_effect %in% TRUE & !grepl(":", y$term), , drop = FALSE]
  if (!nrow(main)) main <- y[y$simulated_effect %in% TRUE, , drop = FALSE]
  ok <- vapply(reps, function(r) { z <- main[main$reps == r, , drop = FALSE]; nrow(z) > 0 && all(z$power >= target) }, logical(1))
  pick <- if (any(ok)) reps[which(ok)[1]] else max(reps)
  reason <- if (any(ok)) sprintf("design.reps = %d: smallest count tried where %s reach power %s", pick, paste(unique(main$term), collapse = ", "), target) else
    sprintf("design.reps = %d: power %s not reached for %s within reps %s; the largest count tried is used", pick, target, paste(unique(main$term), collapse = ", "), paste(reps, collapse = ", "))
  print(ps)
  cat(reason, "\n")
  list(reps = as.integer(pick), reason = reason, target_power = target, n_sim = n_sim, reps_tried = reps, sweep = ps)
}

# Replicate counts to try. Biodes writes power.reps_range as [min, max] (e.g. 3-12): spread six counts over it
# (geometric, as power changes fastest at small counts). A longer vector is used as given.
.suite_reps_grid <- function(x) {
  x <- as.integer(unlist(x))
  if (length(x) != 2L || x[2] - x[1] < 2L) return(sort(unique(x)))
  sort(unique(as.integer(round(exp(seq(log(x[1]), log(x[2]), length.out = 6))))))
}
