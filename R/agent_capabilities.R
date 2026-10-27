# Internal registry loaders and generic planning resolver. The resolver can
# recognize and scope catalogued analyses, but it never turns metadata into
# executable code. Only separately validated generic primitives may do that.
.simitall_agent_capability_path <- function(file) {
  path <- system.file("agent", "capabilities", file, package = "simitall")
  if (!nzchar(path)) path <- file.path("inst", "agent", "capabilities", file)
  if (!file.exists(path)) stop("The SimItAll capability registry is not installed: ", file)
  path
}

.simitall_agent_capability_registry <- function() {
  path <- .simitall_agent_capability_path("workflows.json")
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

.simitall_agent_analysis_catalog <- function() {
  jsonlite::fromJSON(.simitall_agent_capability_path("analysis_catalog.json"), simplifyVector = FALSE)
}

.simitall_agent_input_type_registry <- function() {
  jsonlite::fromJSON(.simitall_agent_capability_path("input_types.json"), simplifyVector = FALSE)
}

.simitall_agent_default_policy <- function() {
  jsonlite::fromJSON(.simitall_agent_capability_path("defaults.json"), simplifyVector = FALSE)
}

.simitall_agent_sample_naming_contract <- function() {
  jsonlite::fromJSON(.simitall_agent_capability_path("sample_naming.json"), simplifyVector = FALSE)
}

.simitall_agent_sample_naming_text <- function() {
  contract <- .simitall_agent_sample_naming_contract()
  paste0(
    "For new data, use `", contract$canonical_sample_id,
    "` (S = unique sample, R = biological replicate, Trt = treatment, Ti = time). ",
    "Do not use a bare `T`: `Trt` and `Ti` remain unambiguous. ",
    "Also provide a metadata table with `",
    paste(contract$metadata_columns$recommended_order, collapse = "`, `"),
    "` as applicable; metadata is authoritative and must match count-matrix columns exactly."
  )
}

.simitall_agent_catalog_uses_sample_metadata <- function(subanalysis) {
  inputs <- tolower(c(subanalysis$required_inputs, subanalysis$optional_inputs))
  any(grepl("sample|replicate|metadata|subject|condition|batch", inputs))
}

.simitall_agent_catalog_goal <- function(question) {
  text <- .simitall_agent_catalog_normalise(question)
  if (grepl("simulate|generate|create", text)) return("simulate")
  if (grepl("compare|test|map|analy[sz]e|call|find|identify", text)) return("analyze")
  "plan"
}

#' Convert a catalog resolution into a schema-compatible partial StudySpec
#'
#' It intentionally retains unanswered design choices and is marked do_not_run.
#' The future final resolver will validate a completed version of this object
#' before generic execution primitives can be composed.
#'
#' @keywords internal
.simitall_agent_catalog_study_spec <- function(resolution) {
  source <- switch(resolution$data_source,
    auto = "ask",
    resolution$data_source
  )
  if (identical(resolution$status, "needs_goal")) {
    family <- resolution$family
    return(list(
      spec_version = "1.0",
      question = resolution$question,
      goal = "plan",
      organism = NULL,
      assembly = NULL,
      region = NULL,
      data_source = source,
      recognized_concepts = family$id,
      primary_analysis = NULL,
      secondary_analyses = character(),
      inputs = list(supplied_paths = resolution$input_files),
      design = list(),
      parameters = list(),
      outputs = character(),
      figures = character(),
      truth = character(),
      pass_criteria = character(),
      assumptions = character(),
      unresolved_questions = paste0(
        "Choose one ", family$label, " goal: ",
        paste(vapply(family$subanalyses, `[[`, character(1), "label"), collapse = "; "), "."
      ),
      execution_status = "needs_clarification",
      do_not_run = TRUE
    ))
  }
  if (!identical(resolution$status, "recognized")) return(NULL)

  primary <- resolution$primary$subanalysis
  parameter_values <- stats::setNames(
    lapply(resolution$parameters, `[[`, "value"),
    vapply(resolution$parameters, `[[`, character(1), "id")
  )
  unresolved <- c(
    paste0("Specify `", gsub("_", " ", resolution$unresolved_parameters, fixed = TRUE), "`."),
    if (!identical(primary$execution_status, "partial")) {
      paste0("The selected route is `", primary$execution_status,
        "`; it can be planned but has no validated generic execution route.")
    } else {
      "Generic execution primitives are not connected yet; review this partial StudySpec before running."
    }
  )
  list(
    spec_version = "1.0",
    question = resolution$question,
    goal = .simitall_agent_catalog_goal(resolution$question),
    organism = NULL,
    assembly = NULL,
    region = NULL,
    data_source = source,
    recognized_concepts = unique(c(
      resolution$primary$family$id,
      primary$id,
      vapply(resolution$secondary, function(x) x$subanalysis$id, character(1))
    )),
    primary_analysis = primary$id,
    secondary_analyses = vapply(resolution$secondary, function(x) x$subanalysis$id, character(1)),
    inputs = list(
      required = primary$required_inputs,
      optional = primary$optional_inputs,
      supplied_paths = resolution$input_files
    ),
    design = list(),
    parameters = parameter_values,
    outputs = primary$expected_outputs,
    figures = primary$recommended_figures,
    truth = if (any(grepl("simulation truth", primary$expected_outputs, fixed = TRUE))) "simulation truth" else character(),
    pass_criteria = primary$qc_and_guardrails,
    assumptions = if (.simitall_agent_catalog_uses_sample_metadata(primary)) .simitall_agent_sample_naming_text() else character(),
    unresolved_questions = unresolved,
    execution_status = if (identical(primary$execution_status, "knowledge_catalog_only")) "knowledge_only" else "needs_clarification",
    do_not_run = TRUE
  )
}

.simitall_agent_catalog_normalise <- function(x) {
  x <- tolower(enc2utf8(paste(x, collapse = " ")))
  x <- gsub("[^a-z0-9]+", " ", x)
  trimws(gsub("\\s+", " ", x))
}

.simitall_agent_catalog_phrase_match <- function(text, phrase) {
  phrase <- .simitall_agent_catalog_normalise(phrase)
  if (!nzchar(phrase)) return(FALSE)
  grepl(paste0("(^| )", phrase, "( |$)"), text, perl = TRUE)
}

.simitall_agent_catalog_is_input_format_cue <- function(cue) {
  extensions <- unlist(lapply(.simitall_agent_input_type_registry()$types, `[[`, "extensions"),
    use.names = FALSE
  )
  formats <- vapply(
    gsub("\\.gz$", "", sub("^\\.", "", extensions)),
    .simitall_agent_catalog_normalise,
    character(1)
  )
  .simitall_agent_catalog_normalise(cue) %in% formats
}

.simitall_agent_catalog_score <- function(text, cues) {
  cues <- unique(unlist(cues, use.names = FALSE))
  cues <- cues[nzchar(cues)]
  text_terms <- strsplit(text, " ", fixed = TRUE)[[1L]]
  stemmed_text_terms <- sub("s$", "", text_terms)
  scores <- vapply(cues, function(cue) {
    cue_terms <- strsplit(.simitall_agent_catalog_normalise(cue), " ", fixed = TRUE)[[1L]]
    if (.simitall_agent_catalog_phrase_match(text, cue)) {
      # Exact phrases are strongest, but a phrase may contain harmless inserted
      # words (for example "simulate bulk RNA-seq counts").
      return(length(cue_terms) + 3L)
    }
    cue_stems <- sub("s$", "", cue_terms)
    if (length(cue_terms) >= 2L && all(cue_stems %in% stemmed_text_terms)) {
      return(length(cue_terms))
    }
    if (length(cue_terms) == 1L && cue_stems %in% stemmed_text_terms) return(1L)
    0L
  }, integer(1))
  max(scores)
}

.simitall_agent_catalog_parameter_value <- function(text, parameter) {
  allowed <- parameter$allowed_values
  if (!is.null(allowed) && length(allowed)) {
    hits <- vapply(allowed, function(value) {
      .simitall_agent_catalog_phrase_match(text, gsub("_", " ", value, fixed = TRUE))
    }, logical(1))
    if (any(hits)) return(allowed[[which(hits)[1L]]])
  }
  id <- parameter$id %||% ""
  if (identical(id, "replicates")) {
    hit <- regmatches(text, regexec("([0-9]+) biological replicates?", text, perl = TRUE))[[1L]]
    if (length(hit) == 2L) return(as.integer(hit[2L]))
  }
  NULL
}

#' Resolve a question against the static analysis catalog
#'
#' This deliberately produces a partial planning specification, not runnable
#' code. It provides one generic bridge between audited aliases/catalog records
#' and the current agent while executable composition is still being built.
#'
#' @keywords internal
.simitall_agent_catalog_resolve <- function(question, data_source = "auto",
                                             input_files = character()) {
  catalog <- .simitall_agent_analysis_catalog()
  text <- .simitall_agent_catalog_normalise(question)
  families <- catalog$families
  family_scores <- vapply(families, function(family) {
    cues <- c(family$aliases, family$label, family$id)
    # A file format such as GTF is evidence for an input contract, not an
    # automatic request for gene-annotation analysis.
    cues <- cues[!vapply(cues, .simitall_agent_catalog_is_input_format_cue, logical(1))]
    .simitall_agent_catalog_score(text, cues)
  }, integer(1))
  recognized_family_ids <- vapply(families[family_scores > 0L], `[[`, character(1), "id")
  candidates <- list()
  for (index in seq_along(families)) {
    family <- families[[index]]
    for (subanalysis in family$subanalyses) {
      score <- .simitall_agent_catalog_score(
        text,
        c(subanalysis$aliases, subanalysis$label, subanalysis$id)
      )
      if (score > 0L) {
        candidates[[length(candidates) + 1L]] <- list(
          family = family,
          subanalysis = subanalysis,
          score = score
        )
      }
    }
  }
  if (!length(candidates) && !any(family_scores > 0L)) {
    return(list(status = "unrecognized", question = question))
  }

  if (length(candidates)) {
    scores <- vapply(candidates, `[[`, integer(1), "score")
    # Stable ordering makes ties reproducible: catalog order is the tie-breaker.
    candidates <- candidates[order(scores, decreasing = TRUE)]
    primary <- candidates[[1L]]
    # A secondary analysis must have its own direct alias/label match. We never
    # infer eQTL, GWAS, or other downstream work from a broad upstream request.
    secondary <- candidates[-1L]
    secondary <- Filter(function(candidate) {
      candidate$subanalysis$id != primary$subanalysis$id
    }, secondary)
    resolved_parameters <- lapply(primary$subanalysis$parameters, function(parameter) {
      list(
        id = parameter$id,
        value = .simitall_agent_catalog_parameter_value(text, parameter),
        required = isTRUE(parameter$required),
        missing_value_policy = parameter$missing_value_policy
      )
    })
    unresolved <- vapply(resolved_parameters, function(parameter) {
      is.null(parameter$value) && isTRUE(parameter$required) &&
        parameter$missing_value_policy %in% c("ask", "required_for_DE")
    }, logical(1))
    unresolved <- vapply(resolved_parameters[unresolved], `[[`, character(1), "id")
    return(list(
      status = "recognized",
      question = question,
      data_source = data_source,
      input_files = as.character(input_files),
      primary = primary,
      secondary = secondary,
      recognized_family_ids = recognized_family_ids,
      parameters = resolved_parameters,
      unresolved_parameters = unresolved
    ))
  }

  family <- families[[which.max(family_scores)]]
  list(
    status = "needs_goal",
    question = question,
    data_source = data_source,
    input_files = as.character(input_files),
    family = family,
    recognized_family_ids = family$id,
    unresolved_parameters = character()
  )
}

.simitall_agent_catalog_plan <- function(resolution, preflight = NULL) {
  if (identical(resolution$status, "needs_goal")) {
    family <- resolution$family
    choices <- vapply(family$subanalyses, `[[`, character(1), "label")
    return(paste(
      "simitall catalog-based planning result (not runnable yet).",
      paste0("Question: ", resolution$question),
      paste0("Clarification needed:\n- You asked about ", family$label,
        ", which covers several different goals.\n- Choose the intended goal: ",
        paste(choices, collapse = "; "), "."),
      "No downstream analysis has been added automatically.",
      sep = "\n\n"
    ))
  }

  primary <- resolution$primary$subanalysis
  secondary <- resolution$secondary
  declared <- unlist(lapply(resolution$parameters, function(parameter) {
    if (is.null(parameter$value)) return(character())
    paste0("- ", parameter$id, ": ", parameter$value)
  }), use.names = FALSE)
  declared <- declared[nzchar(declared)]
  required <- paste0("- ", primary$required_inputs)
  if (.simitall_agent_catalog_uses_sample_metadata(primary)) {
    required <- c(required, paste0("- ", .simitall_agent_sample_naming_text()))
  }
  optional <- if (length(primary$optional_inputs)) {
    paste0("- Optional: ", primary$optional_inputs)
  } else character()
  unresolved <- if (length(resolution$unresolved_parameters)) {
    paste0("- Specify: ", gsub("_", " ", resolution$unresolved_parameters, fixed = TRUE))
  } else "- No mandatory catalog parameter remains unresolved; implementation-specific assumptions still need disclosure."
  secondary_text <- if (!length(secondary)) {
    "- None. No unrequested downstream analysis was inferred."
  } else {
    paste0("- ", vapply(secondary, function(x) x$subanalysis$label, character(1)))
  }
  input_note <- if (!is.null(preflight) && !identical(preflight$status, "ready")) {
    paste0("\nInput status: ", preflight$action)
  } else ""
  paste(
    "simitall catalog-based planning result (not runnable yet).",
    paste0("Question: ", resolution$question),
    paste0("Interpreted design:\n- Primary analysis: ", primary$label,
      " (", resolution$primary$family$label, ").\n",
      "- Data source: ", resolution$data_source, ".\n",
      "- Explicit secondary analyses:\n", paste(secondary_text, collapse = "\n"),
      if (length(declared)) paste0("\n- Recognized parameters:\n", paste(declared, collapse = "\n")) else ""
    ),
    paste0("Required inputs:\n", paste(c(required, optional), collapse = "\n"), input_note),
    paste0("Required choices before execution:\n", paste(unresolved, collapse = "\n")),
    paste0("Expected results:\n", paste0("- ", primary$expected_outputs, collapse = "\n")),
    paste0("Required figures:\n", paste0("- ", primary$recommended_figures, collapse = "\n")),
    paste0("QC and guardrails:\n", paste0("- ", primary$qc_and_guardrails, collapse = "\n")),
    paste0("Execution status: `", primary$execution_status,
      "`. This catalog record can guide planning but cannot generate a recipe until generic primitives and validation are connected."),
    sep = "\n\n"
  )
}
