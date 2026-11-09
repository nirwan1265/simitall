# Organism and technique lexicon (inst/agent/lexicon/).
#
# Deterministic recognition of *which organism* and *which technique or
# analysis* a request is about, from audited tables built out of NCBI Taxonomy,
# EFO, OBI, EDAM, the ENA/SRA vocabularies and curated strain/panel/breed lists
# (see data-raw/lexicon/ and inst/agent/lexicon/SOURCES.md).
#
# These helpers only *read* a prompt. They do not change agent behaviour yet;
# the agent still routes through the legacy keyword spec and the catalog.

.simitall_lexicon_cache <- new.env(parent = emptyenv())

.simitall_lexicon_file <- function(name) {
  p <- system.file("agent", "lexicon", name, package = "simitall")
  if (!nzchar(p)) p <- file.path("inst", "agent", "lexicon", name)
  p
}

.simitall_lexicon_read <- function(name) {
  utils::read.delim(.simitall_lexicon_file(name), sep = "\t", quote = "", comment.char = "",
                    colClasses = "character", na.strings = character(), check.names = FALSE)
}

# Normalise text and aliases the same way: "-", "_" and "/" become spaces.
.simitall_lexicon_norm <- function(x) {
  x <- gsub("[-_/]", " ", x)
  trimws(gsub("[[:space:]]+", " ", x))
}

.simitall_lexicon_prepare <- function(aliases, target_col) {
  a <- aliases
  a$key <- .simitall_lexicon_norm(a$alias)
  a$key <- ifelse(a$case_sensitive == "yes", a$key, tolower(a$key))
  a$target <- a[[target_col]]
  a <- a[nzchar(a$key) & !duplicated(paste(a$key, a$case_sensitive, a$target)), , drop = FALSE]
  # optional plural "s" after aliases that end in a letter ("eQTLs", "isoform-usage QTLs")
  a$pattern <- paste0("(?<![[:alnum:]])", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", a$key),
                      ifelse(grepl("[[:alpha:]]$", a$key), "s?", ""), "(?![[:alnum:]])")
  a[order(-nchar(a$key)), , drop = FALSE]
}

.simitall_lexicon <- function() {
  if (is.null(.simitall_lexicon_cache$lex)) {
    lex <- list(
      organisms = .simitall_lexicon_read("organisms.tsv"),
      techniques = .simitall_lexicon_read("techniques.tsv")
    )
    lex$organism_aliases <- .simitall_lexicon_prepare(.simitall_lexicon_read("organism_aliases.tsv"), "organism")
    lex$technique_aliases <- .simitall_lexicon_prepare(.simitall_lexicon_read("technique_aliases.tsv"), "technique")
    .simitall_lexicon_cache$lex <- lex
  }
  .simitall_lexicon_cache$lex
}

# Find alias hits, longest first, without overlapping spans.
.simitall_lexicon_hits <- function(text, aliases) {
  norm <- .simitall_lexicon_norm(text)
  lower <- tolower(norm)
  taken <- logical(nchar(norm))
  out <- list()
  for (i in seq_len(nrow(aliases))) {
    hay <- if (aliases$case_sensitive[i] == "yes") norm else lower
    if (!grepl(aliases$key[i], hay, fixed = TRUE)) next
    m <- gregexpr(aliases$pattern[i], hay, perl = TRUE)[[1]]
    if (m[1] < 0) next
    for (j in seq_along(m)) {
      span <- m[j]:(m[j] + attr(m, "match.length")[j] - 1)
      if (any(taken[span])) next
      taken[span] <- TRUE
      out[[length(out) + 1]] <- data.frame(target = aliases$target[i], alias = aliases$alias[i],
                                           strength = aliases$strength[i], start = m[j],
                                           stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(data.frame(target = character(), alias = character(), strength = character(), start = integer()))
  h <- do.call(rbind, out)
  h[order(h$start), , drop = FALSE]
}

# Which organism(s) does the text mention?
# Strong aliases win; weak ones ("flies", "patients", "N2") count only when no
# strong alias matched. Generic "bacteria" yields to a named host or species.
.simitall_detect_organism <- function(text) {
  lex <- .simitall_lexicon()
  h <- .simitall_lexicon_hits(text, lex$organism_aliases)
  use <- if (any(h$strength == "strong")) h[h$strength == "strong", , drop = FALSE] else h
  orgs <- unique(use$target)
  if (length(orgs) > 1) orgs <- setdiff(orgs, "bacteria")
  facts <- lex$organisms[match(orgs, lex$organisms$organism), , drop = FALSE]
  list(organism = if (length(orgs)) orgs[1] else NA_character_, all = orgs,
       hits = h, facts = facts, multiple = length(orgs) > 1)
}

.simitall_lexicon_binary_cue <- function(text) {
  grepl("case.?control|cases? and controls?|\\bbinary\\b|binomial|yes or no|yes/no|dead/alive|\\baffected\\b|\\bdisease\\b|logistic|firth",
        tolower(text))
}

# Which techniques/analyses does the text mention, and where would they route?
.simitall_detect_techniques <- function(text) {
  lex <- .simitall_lexicon()
  h <- .simitall_lexicon_hits(text, lex$technique_aliases)
  if (any(h$strength == "strong")) {
    weak_targets <- setdiff(h$target[h$strength == "weak"], h$target[h$strength == "strong"])
    h <- h[h$strength == "strong" | h$target %in% weak_targets, , drop = FALSE]
  }
  t <- lex$techniques[match(unique(h$target), lex$techniques$technique), , drop = FALSE]
  unsupported <- t$technique[t$support == "not_supported"]
  analyses <- character()
  params <- list()
  ambiguous <- character()
  for (i in seq_len(nrow(t))) {
    r <- t[i, ]
    if (nzchar(r$routes_to)) {
      analyses <- c(analyses, r$routes_to)
      if (nzchar(r$set_parameter)) {
        kv <- strsplit(r$set_parameter, "=", fixed = TRUE)[[1]]
        params[[kv[1]]] <- kv[2]
      }
    } else if (r$technique == "gwas") {
      analyses <- c(analyses, if (.simitall_lexicon_binary_cue(text)) "binary_gwas" else "quantitative_gwas")
    } else if (nzchar(r$candidates)) {
      ambiguous <- c(ambiguous, r$technique)
    }
  }
  # An assay with several candidate analyses (bulk RNA-seq, scRNA-seq, WGS) is
  # resolved only by a specific analysis word that is one of its candidates.
  unresolved <- Filter(function(a) {
    cand <- strsplit(t$candidates[t$technique == a], ";", fixed = TRUE)[[1]]
    !any(cand %in% analyses)
  }, ambiguous)
  resolved <- unique(analyses)
  # A population or phenotype set-up yields to the analysis run on it:
  # "NAM GWAS" is a GWAS on a NAM, "F2 ... map QTL" is QTL mapping.
  setup <- c("cross_design", "phenotype_outcome_model")
  if (length(resolved) > 1 && resolved[1] %in% setup && any(!resolved %in% setup)) {
    resolved <- c(resolved[!resolved %in% setup], resolved[resolved %in% setup])
  }
  list(
    techniques = t$technique,
    hits = h,
    unsupported = unsupported,
    unresolved_assays = unlist(unresolved),
    primary_analysis = if (length(unsupported) || !length(resolved)) NA_character_ else resolved[1],
    secondary_analyses = if (length(unsupported)) character() else resolved[-1],
    needs_goal = !length(unsupported) && !length(resolved),
    candidates = unique(unlist(strsplit(t$candidates[t$technique %in% ambiguous], ";", fixed = TRUE))),
    parameters = params,
    execution_status = if (length(unsupported)) "unsupported" else if (length(resolved)) "recognized" else "needs_clarification"
  )
}

# Biology checks that need no model: selfing in species that cannot self,
# diploid assumptions in polyploids, calendar time without a generation time.
.simitall_lexicon_checks <- function(text, organism = .simitall_detect_organism(text)$organism) {
  lex <- .simitall_lexicon()
  if (is.na(organism)) return(character())
  f <- lex$organisms[lex$organisms$organism == organism, , drop = FALSE]
  if (!nrow(f)) return(character())
  low <- tolower(text)
  flags <- character()
  if (grepl("\\bself(ing|ed|s)?\\b|self.?(pollinat|fertili[sz])", low) && f$can_self == "no") {
    flags <- c(flags, sprintf("selfing is not possible in %s (%s)", organism, f$mating_system))
  }
  ploidy <- suppressWarnings(as.numeric(f$ploidy))
  if (!is.na(ploidy) && ploidy > 2 && grepl("\\bdiploid\\b", low)) {
    flags <- c(flags, sprintf("%s is %s (ploidy %s); diploid assumptions conflict", organism,
                              if (nzchar(f$ploidy_note)) f$ploidy_note else "polyploid", f$ploidy))
  }
  if (grepl("\\b(years?|months?|weeks?|days?) ago\\b|over \\d+ (years|months)|\\d+ (years|months)\\b", low) &&
      !grepl("generation", low)) {
    flags <- c(flags, sprintf("calendar time given; convert to generations (%s generation time: %s)", organism,
                              if (nzchar(f$generation_time)) f$generation_time else "unknown; ask"))
  }
  flags
}
