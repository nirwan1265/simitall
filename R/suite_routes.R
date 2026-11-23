# Which simulators a Biodes design runs, and with which parameters.
#
# Routing is a lookup on Biodes' StudySpec fields (never on the user's words:
# Biodes already read those). inst/suite/simulators.tsv maps a field value
# (primary_outcome.type, assays, population.mating, study_type) to a simulator;
# inst/suite/simulator_parameters.tsv gives each simulator's arguments, their
# defaults and the spec field that overrides a default. Any other argument of a
# package simulator keeps that function's own default and can be set through
# `params`.

.suite_tsv <- function(name) {
  p <- system.file("suite", paste0(name, ".tsv"), package = "simitall")
  if (!nzchar(p)) p <- file.path("inst", "suite", paste0(name, ".tsv"))
  utils::read.delim(p, stringsAsFactors = FALSE, colClasses = "character", na.strings = character(), quote = "")
}

# Value at a dotted path ("population.mating") in a nested list.
.suite_get <- function(x, path) {
  for (k in strsplit(path, ".", fixed = TRUE)[[1]]) {
    if (!is.list(x) || is.null(x[[k]])) return(NULL)
    x <- x[[k]]
  }
  x
}

.suite_norm <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(as.character(unlist(x)))))

# Routes matched by a spec: one row per matched simulators.tsv row.
.suite_match_routes <- function(spec) {
  tab <- .suite_tsv("simulators")
  hit <- logical(nrow(tab))
  matched_value <- character(nrow(tab))
  assay_vals <- .suite_assay_values(spec)
  for (i in seq_len(nrow(tab))) {
    vals <- if (tab$field[i] == "assays") assay_vals else .suite_norm(.suite_get(spec, tab$field[i]))
    want <- strsplit(tab$values[i], "|", fixed = TRUE)[[1]]
    m <- intersect(vals, want)
    if (length(m)) { hit[i] <- TRUE; matched_value[i] <- paste(m, collapse = ",") }
  }
  out <- tab[hit, , drop = FALSE]
  out$matched <- matched_value[hit]
  # no outcome type given: the continuous outcome
  if (!any(out$layer == "outcome")) {
    row <- tab[tab$route == "outcome_continuous", , drop = FALSE]
    row$matched <- "(not given)"
    out <- rbind(row, out)
  }
  # several genetics rows can match (assays gwas + population.mating ril): the cross wins over the panel
  gen <- out$layer == "genetics"
  if (sum(gen) > 1) {
    keep <- if (any(out$simulator[gen] == "suite_genotype_cross")) which(gen & out$simulator == "suite_genotype_cross")[1] else which(gen)[1]
    out$matched[keep] <- paste(unique(out$matched[gen]), collapse = ",")
    out <- out[!gen | seq_len(nrow(out)) == keep, , drop = FALSE]
  }
  # eQTL needs genotypes: turn the panel on when nothing else did
  if (any(out$route == "assay_eqtl") && !any(out$layer == "genetics")) {
    row <- tab[tab$route == "genetics_panel", , drop = FALSE]
    row$matched <- "eqtl"
    out <- rbind(out, row)
  }
  # an assay listed twice under different names runs once
  out <- out[!duplicated(out$route), , drop = FALSE]
  rownames(out) <- NULL
  out[order(match(out$layer, c("outcome", "genetics", "assay"))), , drop = FALSE]
}

# The assay values routing sees: spec$assays (coarse) plus spec$assay_kinds (finer, Biodes >= 0.5.1). A coarse
# label is dropped when a finer kind replaces it and none keeps it (inst/suite/assay_kinds.tsv: amplicon with only
# amplicon_genotyping is not a microbiome).
.suite_assay_values <- function(spec) {
  coarse <- .suite_norm(spec$assays)
  kinds <- .suite_norm(spec$assay_kinds)
  if (!length(kinds)) return(coarse)
  k <- .suite_tsv("assay_kinds")
  k <- k[k$kind %in% kinds, , drop = FALSE]
  drop <- setdiff(k$coarse[k$action == "replace"], k$coarse[k$action == "keep"])
  unique(c(setdiff(coarse, drop), kinds))
}

# Parse one simulator_parameters.tsv value by its type.
.suite_parse_value <- function(value, type) {
  if (!nzchar(value)) return(NULL)
  switch(type,
         integer = as.integer(value),
         numeric = as.numeric(value),
         logical = as.logical(value),
         numeric_vector = as.numeric(strsplit(value, ",", fixed = TRUE)[[1]]),
         value)
}

# Arguments for one simulator: user `params` > spec field > table default.
# Returns list(args, source) where source says where each value came from.
.suite_args <- function(simulator, spec, params = NULL) {
  tab <- .suite_tsv("simulator_parameters")
  tab <- tab[tab$simulator == simulator, , drop = FALSE]
  args <- list(); src <- character()
  for (i in seq_len(nrow(tab))) {
    a <- tab$argument[i]
    # spec_field may list alternatives ("population.self_generations|population.generations"): the first given wins
    sv <- NULL; sf <- ""
    for (f in if (nzchar(tab$spec_field[i])) strsplit(tab$spec_field[i], "|", fixed = TRUE)[[1]] else character()) {
      sv <- .suite_get(spec, f); sf <- f
      if (!is.null(sv) && length(sv)) break
    }
    if (!is.null(sv) && length(sv)) {
      args[a] <- list(sv); src[a] <- paste0("spec ", sf)
    } else {
      v <- .suite_parse_value(tab$value[i], tab$type[i])
      if (!is.null(v)) { args[a] <- list(v); src[a] <- "default (simulator_parameters.tsv)" }
    }
  }
  user <- params[[simulator]] %||% list()
  for (a in names(user)) { args[a] <- list(user[[a]]); src[a] <- "params" }
  list(args = args, source = src)
}

#' Which simulators a Biodes design will run
#'
#' Reads a Biodes design folder and lists, without simulating anything, the
#' simulators that [simitall_from_biodes()] will run (from
#' `inst/suite/simulators.tsv`, matched on the spec's `primary_outcome.type`,
#' `assays`, `population.mating` and `study_type`), the parameters each gets
#' and where each value comes from (your `params`, a spec field, or the default
#' in `inst/suite/simulator_parameters.tsv`). Requests with no simulator yet
#' (`planned`, `not_supported`) are listed too.
#'
#' @param biodes_dir Folder written by Biodes `write_design_package()`.
#' @param params Optional overrides: a list named by simulator, e.g.
#'   `list(simulate_chipseq = list(n_peaks = 500), suite_genetics = list(heritability = 0.8))`.
#' @return Invisibly, a list with `routes` (data frame) and `args` (per simulator).
#' @examples
#' \dontrun{
#' simitall_routes("~/biodes_runs/my_design")
#' }
#' @export
simitall_routes <- function(biodes_dir, params = NULL) {
  h <- .suite_read_biodes(biodes_dir)
  routes <- .suite_match_routes(h$spec)
  sims <- unique(c(routes$simulator[routes$status == "ready"],
                   if (any(routes$layer == "genetics" & routes$status == "ready")) "suite_genetics",
                   if (any(routes$route %in% c("assay_chip_seq", "assay_atac_seq"))) "suite_peak_genome"))
  sims <- setdiff(sims, "-")
  # template-based count assays share simulate_counts_from_template()'s parameters
  if (any(grepl("_template$", sims))) sims <- unique(c(sims[!grepl("_template$", sims)], "simulate_counts_from_template"))
  sims[sims == "methylation_simulator"] <- "simulate_methylation"
  args <- stats::setNames(lapply(sims, function(s) .suite_args(s, h$spec, params)), sims)
  # outcome settings: only those of the matched outcome type
  keep <- switch(routes$route[routes$layer == "outcome"][1],
                 outcome_binary = "prevalence", outcome_count = c("count_distribution", "dispersion", "min_mean"),
                 outcome_survival = c("median_time", "time_sd", "censored_fraction"),
                 outcome_ordinal = "category_probabilities", character())
  if (!is.null(args$suite_outcome)) {
    args$suite_outcome$args <- args$suite_outcome$args[intersect(names(args$suite_outcome$args), keep)]
    args$suite_outcome$source <- args$suite_outcome$source[intersect(names(args$suite_outcome$source), keep)]
  }
  cat("simitall routes for design '", h$spec$design$chosen %||% "?", "' (", h$model, "):\n", sep = "")
  for (i in seq_len(nrow(routes))) {
    r <- routes[i, ]
    cat(sprintf("  [%s] %-22s %-30s %s (%s = %s)\n", r$layer, r$route, r$simulator, r$status, r$field, r$matched))
  }
  for (s in names(args)) {
    a <- args[[s]]
    if (!length(a$args)) next
    cat("  ", s, ":\n", sep = "")
    for (n in names(a$args)) cat(sprintf("    %-30s %-22s %s\n", n, paste(format(a$args[[n]]), collapse = ","), a$source[[n]]))
  }
  # questions worth answering (simulator_parameters.tsv ask column) whose answer is still a default
  tab <- .suite_tsv("simulator_parameters")
  q <- tab[tab$simulator %in% names(args) & nzchar(tab$ask), , drop = FALSE]
  src_of <- function(i) { v <- unname(args[[q$simulator[i]]]$source[q$argument[i]]); if (is.na(v)) "" else v }
  q <- q[!vapply(seq_len(nrow(q)), function(i) src_of(i) == "params" || startsWith(src_of(i), "spec"), logical(1)), , drop = FALSE]
  # only arguments this design uses (e.g. prevalence only for a binary outcome); no-default ones always apply
  q <- q[!nzchar(q$value) | vapply(seq_len(nrow(q)), function(i) q$argument[i] %in% names(args[[q$simulator[i]]]$args), logical(1)), , drop = FALSE]
  if (nrow(q)) {
    cat("Questions (answer with params = list(<simulator> = list(<argument> = ...)); otherwise the default is used):\n")
    for (i in seq_len(nrow(q))) cat(sprintf("  - %s [%s$%s, default %s]\n", q$ask[i], q$simulator[i], q$argument[i], if (nzchar(q$value[i])) q$value[i] else "from the data"))
  }
  invisible(list(routes = routes, args = args, questions = q[, c("simulator", "argument", "ask", "value")]))
}

#' Every parameter a simulator takes
#'
#' Lists all arguments of a simitall simulator with the function's own default,
#' the value the Biodes bridge uses (`inst/suite/simulator_parameters.tsv`) and
#' the spec field that can set it. Use it to see what can be passed through
#' `params` in [simitall_from_biodes()].
#'
#' @param simulator A simulator name as in `inst/suite/simulators.tsv`, e.g.
#'   `"simulate_chipseq"` or `"suite_genetics"`.
#' @return A data frame, one row per argument.
#' @export
simitall_simulator_parameters <- function(simulator) {
  tab <- .suite_tsv("simulator_parameters")
  tab <- tab[tab$simulator == simulator, , drop = FALSE]
  fn <- if (exists(simulator, envir = asNamespace("simitall"), inherits = FALSE) && !startsWith(simulator, "suite_")) {
    get(simulator, envir = asNamespace("simitall"))
  } else NULL
  if (is.null(fn) && !nrow(tab)) stop("Unknown simulator: ", simulator, call. = FALSE)
  fargs <- if (is.null(fn)) character() else setdiff(names(formals(fn)), "...")
  all_args <- unique(c(tab$argument, fargs))
  fdef <- vapply(all_args, function(a) {
    if (!a %in% fargs) return(NA_character_)
    txt <- paste(deparse(formals(fn)[[a]]), collapse = "")
    if (nzchar(txt)) txt else "(required)"
  }, "")
  i <- match(all_args, tab$argument)
  data.frame(argument = all_args, function_default = unname(fdef),
             bridge_value = ifelse(is.na(i), NA, tab$value[i]),
             spec_field = ifelse(is.na(i), NA, tab$spec_field[i]),
             note = ifelse(is.na(i), NA, tab$note[i]), stringsAsFactors = FALSE)
}
