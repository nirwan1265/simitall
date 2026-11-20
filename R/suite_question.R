#' From a plain-language question to simulated and analysed data
#'
#' One call across the suite: Biodes (package `studyforge`,
#' github.com/nirwan1265/BioDesignGPT) reads the question, chooses the design
#' and lists what it still needs to know; the design folder is written; Simitall
#' simulates data under that exact design ([simitall_from_biodes()]); Bioana fits
#' the model the design calls for and checks it recovers the truth
#' ([bioana_fit_biodes()]).
#'
#' Biodes' open questions are shown before anything is simulated. Answer them
#' with `answers` (Biodes field paths, as in `studyforge::studyforge_fix()`),
#' e.g. `answers = list(design.reps = 4, subsamples_per_unit = 3,
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
#' @param force Simulate even when the layout Biodes wrote does not match its
#'   own plan (a different design, or placeholder factor levels). Off by
#'   default: that would simulate a different experiment from the one asked.
#' @return Invisibly: the Biodes answer, the design folder, and the simulation
#'   and analysis results.
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
                                   force = FALSE) {
  engine <- match.arg(engine)
  if (!requireNamespace("studyforge", quietly = TRUE)) {
    stop("simitall_from_question() needs Biodes (R package 'studyforge'). Load it with devtools::load_all(\"~/Documents/Github/BioDesignGPT\") or install it with remotes::install_github(\"nirwan1265/BioDesignGPT\").",
         call. = FALSE)
  }
  ans <- studyforge::studyforge_ask(question, engine = engine, fix = answers)
  outdir <- path.expand(outdir %||% file.path("~", "biodes_runs", .suite_slug(question)))

  # the layout code Biodes generated for its first-choice design
  env <- new.env(parent = asNamespace("studyforge"))
  eval(parse(text = ans$recipe), envir = env)
  book <- if (exists("field", envir = env, inherits = FALSE)) env$field else if (exists("book", envir = env, inherits = FALSE)) env$book else NULL
  if (is.null(book)) {
    stop(sprintf("Biodes chose '%s', which writes an allocation (trays or vials) rather than a field book; the bridge cannot simulate it yet. Biodes' code:\n%s",
                 ans$candidates[1] %||% "?", ans$recipe), call. = FALSE)
  }
  studyforge::write_design_package(book, ans$spec, outdir, ans$audit)
  open <- .suite_show_plan(ans, outdir)
  issues <- .suite_layout_issues(ans, outdir)
  if (length(issues)) {
    cat("Layout does NOT match Biodes' plan (simulating it would test a different experiment):\n",
        paste0("  ! ", issues, "\n"), sep = "")
    open <- c(open, issues)
  }

  out <- list(answer = ans, dir = outdir, open_questions = open, layout_issues = issues, simulation = NULL, analysis = NULL)
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
  invisible(out)
}

.suite_slug <- function(x, n = 6) {
  w <- strsplit(gsub("[^a-z0-9 ]", " ", tolower(x)), "\\s+")[[1]]
  w <- w[nchar(w) > 2 & !w %in% c("the", "and", "for", "with", "per", "are", "from", "measure")]
  paste(utils::head(w, n), collapse = "_")
}

# Print what Biodes understood and what it still needs; return the open items.
.suite_show_plan <- function(ans, outdir) {
  s <- ans$spec
  spec <- jsonlite::read_json(file.path(outdir, "study_spec.json"), simplifyVector = FALSE)
  fl <- function(x) if (is.null(x) || !length(x)) "-" else paste(unlist(x), collapse = ", ")
  factors <- vapply(s$factors %||% list(), function(f) sprintf("%s (%s levels)", f$name %||% "?", f$n_levels %||% "?"), "")
  missing <- unlist(ans$validation$missing)
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
      "  design folder:     ", outdir, "\n", sep = "")
  if (length(open)) cat("Open questions (answer with answers = list(...)):\n", paste0("  - ", open, "\n"), sep = "")
  open
}

# Does the written layout match what Biodes planned? Biodes' layout code is a
# template; flag a different design or factors whose level counts do not appear.
.suite_layout_issues <- function(ans, outdir) {
  spec <- jsonlite::read_json(file.path(outdir, "study_spec.json"), simplifyVector = FALSE)
  layout <- utils::read.csv(file.path(outdir, "field_layout.csv"), check.names = FALSE, colClasses = "character")
  issues <- character()
  chosen <- spec$design$chosen %||% ""
  first <- ans$candidates[1] %||% chosen
  if (nzchar(first) && !identical(chosen, first)) {
    issues <- c(issues, sprintf("Biodes' first-choice design is '%s' but its layout code produced '%s'.", first, chosen))
  }
  counts <- vapply(layout, function(x) length(unique(x)), integer(1))
  for (f in spec$factors %||% list()) {
    n <- suppressWarnings(as.integer(f$n_levels %||% NA))
    if (!is.na(n) && n > 1 && !any(counts == n)) {
      issues <- c(issues, sprintf("factor '%s' has %d levels in the spec, but no layout column has %d levels.", f$name, n, n))
    }
  }
  if (grepl("\"a1\"|\"b1\"", ans$recipe)) {
    issues <- c(issues, "the layout uses placeholder factor levels (a1/a2, b1/b2), not the study's factors.")
  }
  issues
}
