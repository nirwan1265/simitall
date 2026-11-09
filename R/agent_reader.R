# The new prompt reader inside simitall_ask() (option `simitall.reader`).
#
# reader = "legacy" (default): routing is unchanged; a fast lexicon-only
#   reading is attached as `answer$reading` so the app can show it side by side.
# reader = "new": the lexicon tables + trained model (.simitall_read_prompt)
#   decide which route may act:
#     stop     the reading says do not run (unsupported technique, impossible
#              or out-of-range value, clinical use, real data missing,
#              mid-sentence correction) -> explain and ask; nothing runs.
#     spec     the legacy compiler supports the request AND agrees with the
#              reading -> the compiled, runnable spec is used as before.
#     catalog  otherwise the catalog plans the analysis the reading found, or
#              asks which analysis was meant when the request is ambiguous.
#     fallback nothing was read -> the legacy catalog / provider path.

.simitall_reader_catalog_index <- function() {
  catalog <- .simitall_agent_analysis_catalog()
  idx <- list()
  for (family in catalog$families) for (sub in family$subanalyses) idx[[sub$id]] <- list(family = family, subanalysis = sub)
  idx
}

# "Requested goal: <catalog label>" is appended by the app's clarification card.
.simitall_reader_requested_goal <- function(question) {
  m <- regmatches(question, regexec("Requested goal:\\s*(.+?)\\s*$", question, perl = TRUE))[[1]]
  if (length(m) < 2) return(NULL)
  idx <- .simitall_reader_catalog_index()
  labels <- vapply(idx, function(x) x$subanalysis$label, "")
  hit <- names(labels)[tolower(labels) == tolower(m[2])]
  if (length(hit)) hit[1] else NULL
}

.simitall_reader_value <- function(reading, id) {
  reading$parameters[[id]] %||% reading$design[[id]]
}

# Build a catalog resolution (the object the catalog planner and the app use)
# from a reading.
.simitall_reader_resolution <- function(reading, question, data_source = "auto", input_files = character()) {
  idx <- .simitall_reader_catalog_index()
  primary <- reading$primary_analysis
  if (is.null(primary) && length(reading$candidates)) {
    cands <- intersect(reading$candidates, names(idx))
    if (!length(cands)) return(NULL)
    label <- reading$candidate_label %||% "this kind of data"
    return(list(
      status = "needs_goal", question = question, data_source = data_source,
      input_files = as.character(input_files),
      family = list(id = "reader_candidates", label = label, clarification_mode = "ask_goal",
                    subanalyses = lapply(cands, function(a) idx[[a]]$subanalysis)),
      recognized_family_ids = unique(vapply(cands, function(a) idx[[a]]$family$id, "")),
      unresolved_parameters = character()
    ))
  }
  if (is.null(primary) || is.null(idx[[primary]])) return(NULL)
  p <- idx[[primary]]
  params <- lapply(p$subanalysis$parameters %||% list(), function(parameter) {
    list(id = parameter$id, value = .simitall_reader_value(reading, parameter$id),
         required = isTRUE(parameter$required), missing_value_policy = parameter$missing_value_policy)
  })
  unresolved <- vapply(Filter(function(x) is.null(x$value) && isTRUE(x$required) &&
                                (x$missing_value_policy %||% "") %in% c("ask", "required_for_DE"), params), `[[`, "", "id")
  secondary <- lapply(intersect(unlist(reading$secondary_analyses), names(idx)), function(a) c(idx[[a]], score = NA))
  list(
    status = "recognized", question = question, data_source = data_source, input_files = as.character(input_files),
    primary = c(p, score = NA), secondary = secondary,
    recognized_family_ids = unique(c(p$family$id, vapply(secondary, function(x) x$family$id, ""))),
    parameters = params, unresolved_parameters = unresolved
  )
}

.simitall_reader_technique_label <- function(ids) {
  lex <- .simitall_lexicon()
  out <- lex$techniques$label[match(ids, lex$techniques$technique)]
  ifelse(is.na(out), ids, out)
}

.simitall_reader_analysis_label <- function(id) {
  x <- .simitall_reader_catalog_index()[[id]]
  if (is.null(x)) id else x$subanalysis$label
}

# Plain-language summary of a reading (used in plans and in the app's Design panel).
.simitall_reader_summary <- function(reading, show_checks = TRUE) {
  if (is.null(reading)) return("")
  lex <- .simitall_lexicon()
  line <- function(k, v) if (length(v) && !all(is.na(v)) && any(nzchar(as.character(v)))) paste0("- ", k, ": ", paste(unlist(v), collapse = ", ")) else NULL
  org <- reading$organism
  org_txt <- if (is.null(org)) NULL else {
    f <- lex$organisms[lex$organisms$organism == org, , drop = FALSE]
    if (nrow(f)) paste0(f$scientific_name, " (", org, "; ploidy ", f$ploidy,
                        if (nzchar(f$ploidy_note)) paste0(", ", f$ploidy_note) else "", "; selfing ",
                        f$can_self, if (nzchar(f$generation_time)) paste0("; generation ", f$generation_time) else "", ")") else org
  }
  vals <- c(reading$design, reading$parameters)
  vals <- vals[!vapply(vals, is.null, logical(1))]
  settings <- if (length(vals)) paste0("  - ", gsub("_", " ", names(vals)), ": ",
                                       vapply(vals, function(v) paste(unlist(v), collapse = ", "), ""), collapse = "\n") else NULL
  paste(c(
    paste0("How the request was read (", reading$reader %||% "lexicon", "):"),
    line("Goal", reading$goal),
    line("Organism", org_txt),
    line("Primary analysis", if (!is.null(reading$primary_analysis)) .simitall_reader_analysis_label(reading$primary_analysis)),
    line("Further requested analyses", vapply(unlist(reading$secondary_analyses), .simitall_reader_analysis_label, "")),
    line("Excluded or used as input", reading$negated),
    line("Data source", reading$data_source),
    line("Region", reading$region),
    if (!is.null(settings)) paste0("- Values read from the request:\n", settings),
    line("Not supported yet", if (length(reading$unsupported_techniques)) .simitall_reader_technique_label(reading$unsupported_techniques)),
    if (show_checks && length(reading$flags)) paste0("- Checks:\n", paste0("  - ", reading$flags, collapse = "\n"))
  ), collapse = "\n")
}

.simitall_reader_stop_plan <- function(reading, question) {
  why <- c(
    if (length(reading$unsupported_techniques))
      paste0("This uses a technique simitall does not support yet: ",
             paste(.simitall_reader_technique_label(reading$unsupported_techniques), collapse = ", "),
             ". It will not be substituted with a different assay."),
    reading$flags
  )
  paste(
    "simitall reader: this request was not run.",
    paste0("Question: ", question),
    paste0("Why it stopped:\n", paste0("- ", why, collapse = "\n")),
    .simitall_reader_summary(reading, show_checks = FALSE),
    "What to do: correct or confirm the points above (for example give real input files, a valid value, or a supported technique) and ask again.",
    sep = "\n\n"
  )
}

# Does the legacy compiler agree with the reading?
.simitall_reader_spec_agrees <- function(spec, reading) {
  if (!isTRUE(spec$supported) || is.null(reading$primary_analysis)) return(FALSE)
  legacy <- .simitall_eval_view(list(spec = spec))
  same <- function(a, b) identical(a, b) || all(c(a, b) %in% c("quantitative_gwas", "binary_gwas"))
  if (!same(legacy$primary_analysis %||% "", reading$primary_analysis)) return(FALSE)
  extra <- setdiff(legacy$secondary_analyses, c(reading$primary_analysis, unlist(reading$secondary_analyses)))
  if (length(extra)) return(FALSE)
  if (!is.null(reading$organism) && !is.null(legacy$organism) && !identical(reading$organism, legacy$organism)) return(FALSE)
  TRUE
}

# Decide the route for reader = "new".
.simitall_reader_route <- function(reading, spec, question, data_source = "auto", input_files = character()) {
  goal <- .simitall_reader_requested_goal(question)
  if (!is.null(goal)) {
    reading$primary_analysis <- goal
    reading$secondary_analyses <- setdiff(unlist(reading$secondary_analyses), goal)
    reading$do_not_run <- FALSE
  }
  if (isTRUE(reading$do_not_run)) return(list(action = "stop", reading = reading))
  if (.simitall_reader_spec_agrees(spec, reading)) return(list(action = "spec", reading = reading))
  res <- .simitall_reader_resolution(reading, question, data_source, input_files)
  if (!is.null(res)) return(list(action = "catalog", reading = reading, resolution = res))
  list(action = "fallback", reading = reading)
}

# Fill a catalog StudySpec with what the reading found.
.simitall_reader_enrich_spec <- function(study_spec, reading) {
  if (is.null(study_spec) || is.null(reading)) return(study_spec)
  for (k in c("organism", "region", "goal")) if (is.null(study_spec[[k]]) && !is.null(reading[[k]])) study_spec[[k]] <- reading[[k]]
  if (identical(study_spec$data_source, "ask") && !identical(reading$data_source, "ask")) study_spec$data_source <- reading$data_source
  for (k in names(reading$design)) if (is.null(study_spec$design[[k]])) study_spec$design[[k]] <- reading$design[[k]]
  for (k in names(reading$parameters)) if (is.null(study_spec$parameters[[k]])) study_spec$parameters[[k]] <- reading$parameters[[k]]
  if (length(reading$flags)) study_spec$assumptions <- c(study_spec$assumptions, reading$flags)
  study_spec$reader <- reading$reader
  study_spec
}
