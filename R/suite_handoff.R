# Reading a Biodes design package (the suite hand-off).
#
# Biodes (github.com/nirwan1265/BioDesignGPT) writes study_spec.json,
# field_layout.csv and sample_metadata_template.tsv. Simitall
# (simitall_from_biodes) and Bioana (bioana_fit_biodes) both read the folder
# through these helpers so they agree on the model, the columns and the units.
# Vocabulary lives in inst/suite/*.json; this file only holds the logic.

.suite_rules <- function(name) {
  p <- system.file("suite", paste0(name, ".json"), package = "simitall")
  if (!nzchar(p)) p <- file.path("inst", "suite", paste0(name, ".json"))
  jsonlite::read_json(p, simplifyVector = FALSE)
}

.suite_check_version <- function(spec, rules) {
  v <- spec$handoff_version %||% ""
  mm <- sub("^(\\d+\\.\\d+).*$", "\\1", v)
  if (!nzchar(v)) {
    stop("study_spec.json has no handoff_version: it was written by a Biodes (BioDesignGPT) checkout older than commit 30f184b. Run `git pull` in BioDesignGPT, reload it, and write the design folder again with write_design_package().",
         call. = FALSE)
  }
  if (!mm %in% unlist(rules$supported_handoff_versions)) {
    stop(sprintf(
      "study_spec.json has handoff_version '%s'; this simitall reads %s. Regenerate the package with a matching Biodes, or add the version to inst/suite/handoff.json after checking the schema.",
      v, paste(unlist(rules$supported_handoff_versions), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

.suite_read_biodes <- function(biodes_dir) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("The Biodes bridge needs the lme4 package to read and fit design models: install.packages('lme4')", call. = FALSE)
  }
  rules <- .suite_rules("handoff")
  if (!dir.exists(biodes_dir)) stop("Biodes folder not found: ", biodes_dir, call. = FALSE)
  missing <- setdiff(unlist(rules$required_files), list.files(biodes_dir))
  if (length(missing)) {
    found <- list.files(biodes_dir, pattern = "^study_spec\\.json$", recursive = TRUE)
    is_repo <- file.exists(file.path(biodes_dir, "DESCRIPTION")) &&
      any(grepl("^Package:\\s*studyforge", readLines(file.path(biodes_dir, "DESCRIPTION"), warn = FALSE)))
    hint <- if (length(found)) {
      paste0(" Design folders found inside it: ", paste(file.path(biodes_dir, dirname(utils::head(found, 5))), collapse = ", "), ".")
    } else if (is_repo) {
      paste0(" This is the Biodes source repository, not a design folder. Create one first, e.g.:\n",
             "  book <- make_design(\"rcbd\", paste0(\"G\", 1:6), reps = 4, seed = 42)\n",
             "  write_design_package(book, extract_spec_local(\"<your question>\"), \"~/biodes_runs/test\")\n",
             "then call simitall_from_biodes(\"~/biodes_runs/test\").")
    } else ""
    stop("Not a Biodes design folder (the folder write_design_package() writes); missing: ",
         paste(missing, collapse = ", "), ".", hint, call. = FALSE)
  }
  spec <- jsonlite::read_json(file.path(biodes_dir, "study_spec.json"), simplifyVector = FALSE)
  .suite_check_version(spec, rules)
  layout <- utils::read.csv(file.path(biodes_dir, "field_layout.csv"), stringsAsFactors = FALSE,
                            check.names = FALSE, colClasses = "character")
  meta <- utils::read.delim(file.path(biodes_dir, "sample_metadata_template.tsv"), stringsAsFactors = FALSE,
                            check.names = FALSE, colClasses = "character", na.strings = c("", "NA"))
  model <- spec$analysis$model
  if (is.null(model) || !nzchar(model)) {
    stop("study_spec.json has no analysis.model; regenerate the package with Biodes (handoff 0.4.0 or later).",
         call. = FALSE)
  }
  list(dir = biodes_dir, spec = spec, layout = layout, meta = meta, model = model,
       parsed = .suite_parse_model(model, spec, layout, rules))
}

# Turn analysis.model into variables, roles and layout columns.
.suite_parse_model <- function(model, spec, layout, rules = .suite_rules("handoff")) {
  bad <- Filter(function(tok) grepl(tok, model, fixed = TRUE), unlist(rules$unsupported_model_tokens))
  f <- tryCatch(stats::as.formula(model), error = function(e) NULL)
  if (length(bad) || is.null(f)) {
    stop(sprintf("Design '%s' uses model '%s', which simitall cannot simulate or fit yet%s.",
                 spec$design$chosen %||% "?", model,
                 if (length(bad)) paste0(" (", paste(bad, collapse = ", "), ")") else ""), call. = FALSE)
  }
  term_rules <- .suite_rules("model_terms")
  fixed <- .suite_nobars(f)
  fixed_labels <- attr(stats::terms(fixed), "term.labels")
  bars <- .suite_findbars(f)
  random <- lapply(bars, function(b) {
    g <- b[[3]]
    list(label = paste(deparse(g), collapse = ""), vars = all.vars(g))
  })
  vars <- setdiff(all.vars(f), all.vars(f[[2]]))
  factor_names <- vapply(spec$factors %||% list(), function(x) x$name %||% "", "")
  role <- vapply(vars, function(v) term_rules$terms[[v]]$role %||% term_rules$unknown_term_role, "")
  # variables used inside I(), e.g. I(treatment^2) in response-surface designs, are numeric, not factors
  numeric_vars <- intersect(vars, .suite_numeric_vars(fixed))
  # analysis.terms (Biodes 6e49dab+): which generic terms are layout columns, optionally term -> column
  aterms <- spec$analysis$terms %||% list()
  layout_terms <- unlist(aterms$layout_columns)
  explicit <- aterms$columns %||% list()
  fixed_treat <- intersect(vars[role == "treatment"], unlist(lapply(fixed_labels, function(l) all.vars(str2lang(l)))))
  sep <- term_rules$combined_level_separator %||% ":"
  is_cell_label <- function(column) {
    x <- layout[[column]]
    length(fixed_treat) > 1 && length(x) && all(grepl(sep, x, fixed = TRUE))
  }
  col <- stats::setNames(character(length(vars)), vars)
  for (v in vars) {
    r <- term_rules$terms[[v]]
    # Biodes f8c2317+: analysis.terms.<term> is the layout column itself (treatment = "soil")
    named <- if (is.character(aterms[[v]]) && length(aterms[[v]]) == 1L) aterms[[v]]
    cands <- c(explicit[[v]], named, if (v %in% layout_terms) v)
    if (!is.null(r$from_spec_factor) && length(factor_names) >= r$from_spec_factor) {
      cands <- c(cands, factor_names[[r$from_spec_factor]])
    }
    cands <- unique(c(cands, unlist(r$columns), v))
    ok <- cands[cands %in% names(layout)]
    skipped <- Filter(is_cell_label, ok)
    hit <- setdiff(ok, skipped)[1]
    if (is.na(hit)) {
      stop(sprintf("Model term '%s' in '%s' has no matching column in field_layout.csv (looked for: %s%s). Add the column name to inst/suite/model_terms.json, or have Biodes list it in analysis.terms.",
                   v, model, paste(cands, collapse = ", "),
                   if (length(skipped)) sprintf("; skipped %s: its levels join several factors", paste(skipped, collapse = ", ")) else ""),
           call. = FALSE)
    }
    col[[v]] <- hit
  }
  # two fixed terms read from one column would make the model rank deficient
  fixed_vars <- intersect(vars, all.vars(fixed))
  dup <- col[fixed_vars][duplicated(col[fixed_vars])]
  if (length(dup)) {
    stop(sprintf("Model '%s': terms %s all map to layout column '%s'. Biodes should write one column per model term (or give analysis.terms.columns).",
                 model, paste(fixed_vars[col[fixed_vars] %in% dup], collapse = ", "), dup[1]), call. = FALSE)
  }
  list(formula = f, response = all.vars(f[[2]]), fixed_labels = fixed_labels, random = random,
       vars = vars, col = col, role = role, numeric = numeric_vars)
}

# Variables that appear inside I(...) in a formula.
.suite_numeric_vars <- function(f) {
  out <- character()
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1]], as.name("I"))) out <<- c(out, all.vars(e)) else for (a in as.list(e)[-1]) walk(a)
    }
  }
  walk(f[[length(f)]])
  unique(out)
}

# Model frame: one factor column per model variable, named as in the model.
.suite_frame <- function(data, parsed) {
  out <- data
  for (v in parsed$vars) {
    x <- as.character(data[[parsed$col[[v]]]])
    if (v %in% parsed$numeric) {
      num <- suppressWarnings(as.numeric(x))
      if (anyNA(num[!is.na(x)])) stop(sprintf("Model uses %s as a number (inside I()), but layout column '%s' is not numeric.", v, parsed$col[[v]]), call. = FALSE)
      out[[v]] <- num
      next
    }
    lev <- .suite_natural_sort(unique(x))
    out[[v]] <- factor(x, levels = lev)
  }
  out
}

# Unit id for each metadata row: subsample rows share the unit's id.
.suite_unit_id <- function(meta) {
  if ("unit_id" %in% names(meta)) meta$unit_id else meta$sample_id
}

.suite_write_json <- function(x, path) {
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
  path
}

# lme4 moved nobars()/findbars() to the reformulas package; use it when present.
.suite_nobars <- function(f) {
  if (requireNamespace("reformulas", quietly = TRUE)) reformulas::nobars(f) else lme4::nobars(f)
}
.suite_findbars <- function(f) {
  if (requireNamespace("reformulas", quietly = TRUE)) reformulas::findbars(f) else lme4::findbars(f)
}

# G1, G2, ..., G10 (not G1, G10, G11, G2): order by text with each number compared as a number.
.suite_natural_sort <- function(x) {
  num <- suppressWarnings(as.numeric(x))
  if (!anyNA(num)) return(x[order(num)])
  prefix <- sub("\\d+$", "", x)
  suffix <- suppressWarnings(as.numeric(sub("^.*?(\\d+)$", "\\1", x)))
  suffix[!grepl("\\d$", x)] <- -Inf
  x[order(prefix, suffix, x)]
}
