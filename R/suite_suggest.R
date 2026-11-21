# Suggestions for what the study details leave open. Sources, in order:
#   - Biodes' own registries: traits (what to measure), experimental units per
#     system, setup hazards (each with a citation);
#   - inst/suite/suggestions.tsv: what to compare, outcomes for omics, effect
#     size / SD, matched to the question, system, organism and assays.
# Nothing here is specific to one question.

.suite_suggestion_rules <- function() {
  p <- system.file("suite", "suggestions.tsv", package = "simitall")
  if (!nzchar(p)) p <- file.path("inst", "suite", "suggestions.tsv")
  utils::read.delim(p, colClasses = "character", quote = "", na.strings = character())
}

.suite_context <- function(ans) {
  s <- ans$spec
  tolower(paste(c(ans$question, unlist(s$system), unlist(s$organism), unlist(s$assays), unlist(s$tissue)), collapse = " "))
}

# Missing fields -> data frame (field, suggestion, answer, source, from)
.suite_suggestions <- function(ans) {
  s <- ans$spec
  ctx <- .suite_context(ans)
  missing <- character()
  if (!length(s$factors) || all(vapply(s$factors, function(f) is.null(f$n_levels), logical(1)))) missing <- c(missing, "factors")
  if (is.null(s$experimental_unit)) missing <- c(missing, "experimental_unit")
  if (is.null(s$primary_outcome$name)) missing <- c(missing, "primary_outcome")
  if (is.null(s$power$minimum_detectable_difference) || is.null(s$power$residual_sd)) missing <- c(missing, "power")
  rules <- .suite_suggestion_rules()
  hit <- rules[vapply(rules$context, function(p) grepl(p, ctx, perl = TRUE), logical(1)), , drop = FALSE]
  out <- list()
  add <- function(field, suggestion, answer, source, from) {
    out[[length(out) + 1]] <<- data.frame(field = field, suggestion = suggestion, answer = answer, source = source, from = from, stringsAsFactors = FALSE)
  }
  for (f in missing) {
    rows <- hit[hit$field == f, , drop = FALSE]
    # omics rows first, then the general ones (the table lists the specific contexts first)
    if (f == "primary_outcome" && !nrow(rows) && !is.null(ans$traits) && nrow(ans$traits)) {
      for (i in seq_len(min(4, nrow(ans$traits)))) {
        t <- ans$traits[i, ]
        add(f, sprintf("%s (%s)", t$trait, t$measured_by),
            sprintf("primary_outcome.name = \"%s\", primary_outcome.type = \"continuous\"", tolower(t$trait)), t$source %||% "", "Biodes traits")
      }
    }
    if (f == "experimental_unit") {
      units <- tryCatch(studyforge::load_registry("experimental_units")$units, error = function(e) list())
      sys <- unlist(s$system)
      for (u in Filter(function(u) u$system %in% sys, units)) {
        add(f, sprintf("%s: %s", u$id, u$note %||% ""), sprintf("experimental_unit = \"%s\"", u$id), "", "Biodes experimental units")
      }
    }
    for (i in seq_len(nrow(rows))) add(f, rows$suggestion[i], rows$answer[i], rows$source[i], "suggestions.tsv")
  }
  if (!length(out)) return(NULL)
  res <- do.call(rbind, out)
  res[!duplicated(paste(res$field, res$answer, res$suggestion)), , drop = FALSE]
}

.suite_show_suggestions <- function(ans, sug) {
  if (!is.null(sug) && nrow(sug)) {
    cat("Suggestions for what is missing (put any of these in answers = list(...)):\n")
    for (f in unique(sug$field)) {
      cat("  ", f, ":\n", sep = "")
      x <- sug[sug$field == f, , drop = FALSE]
      for (i in seq_len(min(4, nrow(x)))) {
        cat(sprintf("    - %s%s\n%s", x$suggestion[i], if (nzchar(x$source[i])) paste0("   [", x$source[i], "]") else "",
                    if (nzchar(x$answer[i])) paste0("        ", x$answer[i], "\n") else ""))
      }
    }
    cat("  or run again with use_suggestions = TRUE to fill each missing item with its first suggestion (listed as assumptions).\n")
  }
  hz <- ans$hazards %||% list()
  if (length(hz)) {
    cat("Setup checks from Biodes:\n")
    for (h in utils::head(hz, 3)) cat(sprintf("  - %s -> %s [%s]\n", h$hazard, h$mitigation[[1]], h$source %||% ""))
  }
  invisible(sug)
}

# For each missing field, the first suggestion that fills every still-missing
# part of it (power needs both the difference and the SD), as answers = list(...)
.suite_suggested_answers <- function(sug, spec = NULL) {
  if (is.null(sug) || !nrow(sug)) return(list())
  sug <- sug[nzchar(sug$answer), , drop = FALSE]
  parse_answer <- function(a) eval(parse(text = paste0("list(", a, ")")), envir = baseenv())
  needed <- function(field) {
    if (field != "power") return(NULL)
    c(if (is.null(spec$power$minimum_detectable_difference)) "power.minimum_detectable_difference",
      if (is.null(spec$power$residual_sd)) "power.residual_sd")
  }
  out <- list()
  for (f in unique(sug$field)) {
    opts <- lapply(sug$answer[sug$field == f], parse_answer)
    need <- needed(f)
    pick <- if (length(need)) Filter(function(o) all(need %in% names(o)), opts) else opts
    if (length(pick)) out <- c(out, pick[[1]][setdiff(names(pick[[1]]), names(out))])
  }
  out
}
