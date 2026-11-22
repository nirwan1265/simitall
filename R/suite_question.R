#' From a plain-language question to simulated and analysed data
#'
#' One call across the suite: Biodes (package `studyforge`,
#' github.com/nirwan1265/BioDesignGPT) reads the question, chooses the design
#' and lists what it still needs to know; Simitall builds the layout from those
#' study details ([simitall_layout_from_spec()], rules in
#' `inst/suite/layouts.tsv`) and writes the design folder, then simulates data
#' under that exact design ([simitall_from_biodes()]); Bioana fits the model the
#' design calls for and checks it recovers the truth ([bioana_fit_biodes()]).
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
                                   figures = TRUE) {
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

  # the layout, built here from Biodes' study details for its first-choice design (inst/suite/layouts.tsv)
  design <- design %||% ans$candidates[1]
  if (is.null(design) || is.na(design)) stop("Biodes proposed no design for this question.", call. = FALSE)
  book <- simitall_layout_from_spec(ans$spec, design, seed = layout_seed)
  if (inherits(book, "simitall_layout_needs")) {
    .suite_show_plan(ans, NULL, design)
    .suite_show_suggestions(ans, sug)
    print(book)
    message("Not simulated. Give the missing details with answers = list(...) (suggestions above), or run again with use_suggestions = TRUE, ",
            "or choose another of Biodes' candidate designs with design = (", paste(ans$candidates, collapse = ", "), ").")
    return(invisible(structure(list(answer = ans, dir = NULL, design = design, suggestions = sug, filled_from_suggestions = suggested, layout_needs = book$needs,
                                    open_questions = book$needs, layout_issues = character(), simulation = NULL, analysis = NULL),
                               class = "simitall_suite_run")))
  }
  spec_out <- ans$spec
  spec_out$layout <- list(built_by = "simitall", rule = "inst/suite/layouts.tsv", design = design,
                          factor_map = attr(book, "factor_map"), assumptions = attr(book, "assumptions"))
  studyforge::write_design_package(book, spec_out, outdir, ans$audit)
  open <- .suite_show_plan(ans, outdir, design)
  if (length(suggested)) cat("Filled from suggestions (assumptions, change with answers = list(...)):\n",
                             paste0("  - ", names(suggested), " = ", vapply(suggested, function(v) paste(format(v), collapse = ", "), ""), "\n"), sep = "")
  .suite_show_suggestions(ans, sug)
  if (length(attr(book, "assumptions"))) cat("Layout assumptions:\n", paste0("  - ", attr(book, "assumptions"), "\n"), sep = "")
  issues <- .suite_layout_issues(ans, outdir, design)
  if (length(issues)) {
    cat("Layout does NOT match Biodes' plan (simulating it would test a different experiment):\n",
        paste0("  ! ", issues, "\n"), sep = "")
    open <- c(open, issues)
  }

  out <- structure(list(answer = ans, dir = outdir, design = design, suggestions = sug, filled_from_suggestions = suggested, open_questions = open, layout_issues = issues, simulation = NULL,
                        analysis = NULL), class = "simitall_suite_run")
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
  out$simulation <- simitall_from_biodes(outdir, seed = seed)
  if (analyze) out$analysis <- bioana_fit_biodes(outdir)
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
  spec <- jsonlite::read_json(file.path(outdir, "study_spec.json"), simplifyVector = FALSE)
  layout <- utils::read.csv(file.path(outdir, "field_layout.csv"), check.names = FALSE, colClasses = "character")
  issues <- character()
  chosen <- spec$design$chosen %||% ""
  first <- design %||% ans$candidates[1] %||% chosen
  if (nzchar(first) && !identical(chosen, first)) {
    issues <- c(issues, sprintf("the planned design is '%s' but the layout written is '%s'.", first, chosen))
  }
  # compare the study's factors with the treatment columns the model will use (not rep/block columns)
  parsed <- tryCatch(.suite_parse_model(spec$analysis$model, spec, layout), error = function(e) NULL)
  tvars <- if (is.null(parsed)) character() else names(parsed$role)[parsed$role == "treatment"]
  tcols <- if (is.null(parsed)) names(layout) else unique(unname(parsed$col[tvars]))
  counts <- vapply(layout[tcols], function(x) length(unique(x)), integer(1))
  if (length(tcols) > 1) counts <- c(counts, all = length(unique(unlist(layout[tcols]))))   # diallel parents span two columns
  f <- .suite_spec_factors(spec)
  study <- c(f$treatment, if (!is.null(f$time)) list(f$time))
  n_checks <- if (!is.null(f$check) && !"check" %in% tvars) f$check$n else 0
  for (x in c(study, if (!is.null(f$check)) list(f$check))) {
    n <- x$n
    ok <- any(counts %in% c(n, n + if (identical(x, f$check)) 0 else n_checks)) ||
      (identical(x, f$check) && any(counts == sum(vapply(study, function(z) z$n, 0L), na.rm = TRUE) + n))
    if (!is.na(n) && n > 1 && !ok) {
      issues <- c(issues, sprintf("factor '%s' has %d levels in the spec, but no treatment column in the layout has %d levels (it has %s).",
                                  x$name, n, n, paste(counts, collapse = ", ")))
    }
  }
  if (!length(study) || all(is.na(vapply(study, function(z) z$n, 0L)))) {
    issues <- c(issues, sprintf("the spec gives no factor with a known number of levels, so the layout's %s treatment levels are a default, not the study's.",
                                paste(counts, collapse = " x ")))
  }
  real <- sum(vapply(study, function(z) is.na(z$n) || z$n > 1, logical(1)))
  if (!is.null(parsed) && real > length(tvars)) {
    issues <- c(issues, sprintf("the spec has %d factors (%s) but the model has %d treatment term(s); the design leaves factors out.",
                                real, paste(vapply(study, function(z) z$name, ""), collapse = ", "), length(tvars)))
  }
  issues
}

#' @export
print.simitall_suite_run <- function(x, ...) {
  s <- x$answer$spec
  cat("Biodes -> Simitall -> Bioana\n")
  cat("  question:  ", x$answer$question, "\n", sep = "")
  cat("  folder:    ", x$dir, "\n", sep = "")
  fs <- x$analysis$summary
  if (is.null(fs)) {
    cat("  status:    not simulated", if (length(x$layout_needs)) " (the layout needs more study details)" else if (length(x$layout_issues)) " (layout does not match the plan)", "\n", sep = "")
  } else {
    cat(sprintf("  design:    %s, %s; %d units fitted with %s\n", fs$design, fs$model_from_spec, fs$n_units_fitted, fs$method))
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
  if (length(x$open_questions)) cat("  open questions: ", length(x$open_questions), " (see $open_questions)\n", sep = "")
  cat("  files:     results.tsv, recovery.tsv, fit_summary.json in the folder\n")
  if (length(x$figures)) cat("  figures:   ", paste(basename(x$figures), collapse = ", "), " in ", dirname(x$figures[[1]]), "\n", sep = "")
  invisible(x)
}
