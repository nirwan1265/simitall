# Evaluation of the agent's prompt reading against reviewed/draft
# prompt -> StudySpec examples (inst/agent/prompt_eval/).
#
# The agent currently exposes two partial views of a request: the legacy
# keyword spec (`answer$spec`) and the catalog StudySpec
# (`answer$catalog_study_spec`). `.simitall_eval_view()` merges them into one
# flat view in StudySpec vocabulary so a single scorer can judge both.
# Nothing here changes agent behaviour.

.simitall_eval_examples <- function(path = NULL) {
  path <- path %||% .simitall_eval_default_path()
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  lapply(lines, jsonlite::fromJSON, simplifyVector = FALSE)
}

.simitall_eval_default_path <- function() {
  p <- system.file("agent", "prompt_eval", "study_spec_examples_draft.jsonl", package = "simitall")
  if (!nzchar(p)) p <- file.path("inst", "agent", "prompt_eval", "study_spec_examples_draft.jsonl")
  p
}

.simitall_eval_empty <- function(x) {
  is.null(x) || length(x) == 0 || (length(x) == 1 && (is.na(x) || identical(x, "") ||
    identical(x, "unspecified") || identical(x, "auto")))
}

# Flatten one agent answer into StudySpec vocabulary.
.simitall_eval_view <- function(answer) {
  # With the new reader, the reading is what the agent understood; when the
  # legacy compiler also ran, its values (and its defaults) fill the gaps.
  if (identical(answer$reader, "new") && !is.null(answer$reading)) {
    v <- answer$reading
    if (identical(answer$plan_source, "spec")) {
      lv <- .simitall_eval_view(list(spec = answer$spec))
      for (k in c("organism", "primary_analysis")) if (is.null(v[[k]])) v[[k]] <- lv[[k]]
      for (k in c("design", "parameters")) for (f in names(lv[[k]])) if (is.null(v[[k]][[f]])) v[[k]][[f]] <- lv[[k]][[f]]
      v$defaulted <- lv$defaulted
      v$do_not_run <- FALSE
    }
    if (identical(answer$plan_source, "reader")) v$do_not_run <- TRUE
    return(v)
  }
  legacy <- answer$spec %||% list()
  catalog <- answer$catalog_study_spec %||% list()
  v <- list(design = list(), parameters = list())
  defaulted <- character()

  # analyses -----------------------------------------------------------
  task_map <- c(breeding = "cross_design", population_study = "cross_design",
                single_cell = "scrnaseq_simulation", family_risk = "pedigree_inbreeding",
                peak = "chipseq_atacseq", peak_study = "chipseq_atacseq")
  analysis_map <- c(rnaseq = "bulk_counts_simulation", genomic_selection = "genomic_prediction",
                    eqtl = "bulk_eqtl", single_cell = "scrnaseq_simulation")
  trait_type <- legacy$trait$type %||% NULL
  gwas_id <- if (identical(trait_type, "binary")) "binary_gwas" else "quantitative_gwas"
  primary <- catalog$primary_analysis %||% NULL
  if (.simitall_eval_empty(primary)) {
    task <- legacy$task %||% NULL
    primary <- if (identical(task, "gwas")) gwas_id else if (!is.null(task) && task %in% names(task_map)) task_map[[task]] else NULL
  }
  secondary <- unlist(catalog$secondary_analyses %||% list())
  legacy_analyses <- unlist(legacy$analyses %||% list())
  legacy_secondary <- unname(ifelse(legacy_analyses == "gwas", gwas_id,
    ifelse(legacy_analyses %in% names(analysis_map), analysis_map[legacy_analyses], legacy_analyses)))
  secondary <- setdiff(unique(c(secondary, legacy_secondary)), primary)
  if (!is.null(primary)) v$primary_analysis <- primary
  v$secondary_analyses <- secondary

  # facts --------------------------------------------------------------
  if (!.simitall_eval_empty(legacy$organism)) v$organism <- tolower(legacy$organism)
  ds <- catalog$data_source %||% legacy$data_source %||% NULL
  if (!is.null(ds)) v$data_source <- if (identical(ds, "auto")) "ask" else ds
  if (!is.null(catalog$goal)) v$goal <- catalog$goal

  scheme <- legacy$scheme %||% legacy$population$design %||% NULL
  if (!.simitall_eval_empty(scheme)) v$design$scheme <- scheme
  put <- function(where, key, value, legacy_name) {
    if (.simitall_eval_empty(value) || is.list(value)) return()
    v[[where]][[key]] <<- value
    if (legacy_name %in% unlist(legacy$defaulted)) defaulted <<- c(defaulted, paste(where, key, sep = "."))
  }
  put("design", "founder_count", legacy$n_founders, "n_founders")
  put("design", "selfing_generations", legacy$self_generations, "self_generations")
  put("design", "backcross_generations", legacy$backcrosses, "backcrosses")
  put("design", "population_size", legacy$n_individuals %||% legacy$n_lines, "n_individuals")
  put("design", "donors", legacy$n_donors, "n_donors")
  put("parameters", "heritability", legacy$trait$heritability, "trait.heritability")
  put("parameters", "n_causal", legacy$trait$n_qtl, "trait.n_qtl")
  if (!is.null(trait_type) && !is.null(primary) && grepl("gwas", primary)) {
    v$parameters$trait_type <- if (identical(trait_type, "binary")) "binary" else "continuous"
  }
  if ("population.design" %in% unlist(legacy$defaulted) && !is.null(v$design$scheme)) defaulted <- c(defaulted, "design.scheme")
  for (p in names(catalog$parameters %||% list())) {
    val <- catalog$parameters[[p]]
    if (!.simitall_eval_empty(val) && !is.list(val)) v$parameters[[p]] <- val
  }

  # status -------------------------------------------------------------
  if (isTRUE(legacy$supported)) {
    v$execution_status <- "supported"
    v$do_not_run <- FALSE
  } else if (length(catalog)) {
    v$execution_status <- catalog$execution_status %||% NULL
    v$do_not_run <- isTRUE(catalog$do_not_run)
  } else {
    v$execution_status <- "unsupported"
    v$do_not_run <- TRUE
  }
  v$defaulted <- unique(defaulted)
  v
}

.simitall_eval_get <- function(x, path) {
  for (k in strsplit(path, ".", fixed = TRUE)[[1]]) {
    if (!is.list(x) || is.null(x[[k]])) return(NULL)
    x <- x[[k]]
  }
  x
}

.simitall_eval_flatten <- function(x, prefix = NULL) {
  out <- list()
  for (n in names(x)) {
    key <- if (is.null(prefix)) n else paste(prefix, n, sep = ".")
    if (is.list(x[[n]]) && !is.null(names(x[[n]])) && length(x[[n]])) {
      out <- c(out, .simitall_eval_flatten(x[[n]], key))
    } else {
      out[[key]] <- x[[n]]
    }
  }
  out
}

.simitall_eval_same <- function(got, want) {
  if (is.null(got)) return(FALSE)
  g <- tolower(as.character(unlist(got)))
  w <- tolower(as.character(unlist(want)))
  num <- suppressWarnings(as.numeric(c(g, w)))
  if (length(g) == 1 && length(w) == 1 && !anyNA(num)) return(isTRUE(all.equal(num[1], num[2])))
  setequal(g, w)
}

# Score one example against an agent view.
.simitall_eval_score_one <- function(example, view) {
  exp <- example$expected_spec
  ignore <- c("spec_version", "question", "recognized_concepts", "secondary_analyses", "primary_analysis",
              "execution_status", "do_not_run", "goal", "data_source")
  flat <- .simitall_eval_flatten(exp[setdiff(names(exp), ignore)])
  reading <- vapply(names(flat), function(p) .simitall_eval_same(.simitall_eval_get(view, p), flat[[p]]), logical(1))
  context <- c(
    goal = if (!is.null(exp$goal)) .simitall_eval_same(view$goal, exp$goal) else NA,
    data_source = if (!is.null(exp$data_source)) .simitall_eval_same(view$data_source, exp$data_source) else NA
  )
  primary_ok <- if (!is.null(exp$primary_analysis)) identical(view$primary_analysis %||% "", exp$primary_analysis) else NA
  secondary_ok <- setequal(unlist(view$secondary_analyses), unlist(exp$secondary_analyses))
  analyses <- c(view$primary_analysis, unlist(view$secondary_analyses))
  must_not_hit <- intersect(unlist(example$must_not), analyses)
  status_ok <- if (identical(exp$execution_status, "unsupported") || isTRUE(exp$do_not_run)) {
    isTRUE(view$do_not_run)
  } else NA
  invented <- defaulted <- character()
  for (p in unlist(example$must_ask)) {
    val <- .simitall_eval_get(view, p)
    if (identical(p, "data_source") && identical(val, "ask")) next
    if (!.simitall_eval_empty(val)) {
      if (p %in% view$defaulted) defaulted <- c(defaulted, p) else invented <- c(invented, p)
    }
  }
  data.frame(
    id = example$id, split = example$split %||% NA, workflow = example$workflow %||% NA,
    situation = example$situation %||% NA,
    reading_hit = sum(reading), reading_n = length(reading),
    context_hit = sum(context, na.rm = TRUE), context_n = sum(!is.na(context)),
    primary_ok = primary_ok, secondary_ok = secondary_ok,
    must_not_violations = length(must_not_hit), stop_ok = status_ok,
    invented = length(invented), defaulted_not_asked = length(defaulted),
    misses = paste(names(reading)[!reading], collapse = ";"),
    invented_fields = paste(c(invented, if (length(defaulted)) paste0(defaulted, "(default)")), collapse = ";"),
    got_primary = view$primary_analysis %||% NA,
    stringsAsFactors = FALSE
  )
}

# Run the local agent on every example and score it.
# `view`, when given, maps a prompt straight to a StudySpec view (for example
# `.simitall_lexicon_read_prompt`) and replaces `ask` + `.simitall_eval_view`.
.simitall_evaluate_examples <- function(examples = .simitall_eval_examples(), split = NULL,
                                        ask = function(q) simitall_ask(q, provider = "local"),
                                        view = NULL) {
  if (!is.null(split)) examples <- Filter(function(e) identical(e$split, split), examples)
  rows <- lapply(examples, function(e) {
    v <- if (is.null(view)) {
      answer <- tryCatch(suppressMessages(suppressWarnings(ask(e$prompt))), error = function(err) list())
      .simitall_eval_view(answer)
    } else {
      tryCatch(view(e$prompt), error = function(err) list(design = list(), parameters = list()))
    }
    .simitall_eval_score_one(e, v)
  })
  res <- do.call(rbind, rows)
  pct <- function(a, b) if (sum(b)) round(100 * sum(a) / sum(b), 1) else NA
  summary <- list(
    n = nrow(res),
    reading_recall = pct(res$reading_hit, res$reading_n),
    goal_and_source = pct(res$context_hit, res$context_n),
    primary_analysis = round(100 * mean(res$primary_ok, na.rm = TRUE), 1),
    secondary_exact = round(100 * mean(res$secondary_ok), 1),
    must_not_violations = sum(res$must_not_violations),
    stops_when_it_should = round(100 * mean(res$stop_ok, na.rm = TRUE), 1),
    examples_with_invented_values = sum(res$invented > 0),
    examples_defaulted_instead_of_asking = sum(res$defaulted_not_asked > 0)
  )
  list(summary = summary, by_example = res)
}
