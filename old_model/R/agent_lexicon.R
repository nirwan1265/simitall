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
    lex$settings <- .simitall_lexicon_read("settings.tsv")
    lex$settings$row <- as.character(seq_len(nrow(lex$settings)))
    lex$setting_aliases <- .simitall_lexicon_prepare(lex$settings, "row")
    lex$quantities <- .simitall_lexicon_read("quantities.tsv")
    lex$intent <- .simitall_lexicon_read("intent.tsv")
    lex$intent_routes <- .simitall_lexicon_read("intent_routes.tsv")
    cat_path <- system.file("agent", "capabilities", "analysis_catalog.json", package = "simitall")
    if (!nzchar(cat_path)) cat_path <- file.path("inst", "agent", "capabilities", "analysis_catalog.json")
    catalog <- jsonlite::fromJSON(cat_path, simplifyVector = FALSE)
    set_fields <- character()
    for (f in catalog$families) for (s in f$subanalyses) for (p in s$parameters %||% list()) {
      if (identical(p$type, "set")) set_fields <- c(set_fields, paste(s$id, paste0("parameters.", p$id)))
    }
    lex$set_fields <- set_fields
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
  assumed <- NULL
  if (!length(orgs)) {
    fz <- .simitall_fuzzy_organism(text)
    if (!is.null(fz)) {
      orgs <- fz$organism
      assumed <- fz
    }
  }
  facts <- lex$organisms[match(orgs, lex$organisms$organism), , drop = FALSE]
  list(organism = if (length(orgs)) orgs[1] else NA_character_, all = orgs,
       hits = h, facts = facts, multiple = length(orgs) > 1, assumed = assumed)
}

# Typo fallback ("poplus" -> Populus, "sorgum" -> sorghum), used only when no
# alias matched. Targets are genus names and long common names (>= 6 letters),
# which rarely collide with English words; one edit allowed (two for >= 9).
.simitall_fuzzy_organism <- function(text) {
  lex <- .simitall_lexicon()
  if (is.null(lex$fuzzy_targets)) {
    a <- lex$organism_aliases
    genus <- vapply(strsplit(a$alias[a$alias_type %in% c("scientific_name", "genus")], " "), `[`, "", 1)
    keep <- (a$alias_type %in% c("scientific_name", "genus")) |
      (a$alias_type %in% c("common_name", "informal", "canonical_id") & !grepl(" ", a$alias) & nchar(a$alias) >= 7)
    tg <- data.frame(word = tolower(c(genus, a$alias[keep & !a$alias_type %in% c("scientific_name")])),
                     organism = c(a$organism[a$alias_type %in% c("scientific_name", "genus")],
                                  a$organism[keep & !a$alias_type %in% c("scientific_name")]),
                     stringsAsFactors = FALSE)
    tg <- tg[nchar(tg$word) >= 6 & !grepl("[^a-z]", tg$word), , drop = FALSE]
    tg <- tg[!duplicated(tg$word), , drop = FALSE]
    # a word that names two organisms (a shared genus) is not a safe guess
    amb <- tapply(tg$organism, tg$word, function(x) length(unique(x)))
    lex$fuzzy_targets <- tg[tg$word %in% names(amb)[amb == 1], , drop = FALSE]
    .simitall_lexicon_cache$lex <- lex
  }
  tg <- lex$fuzzy_targets
  words <- unique(tolower(unlist(regmatches(text, gregexpr("[A-Za-z]{5,}", text)))))
  for (w in words) {
    d <- utils::adist(w, tg$word)[1, ]
    lim <- ifelse(nchar(tg$word) >= 9, 2, 1)
    ok <- which(d <= lim & d > 0 & abs(nchar(tg$word) - nchar(w)) <= lim)
    if (length(ok)) {
      best <- ok[which.min(d[ok])]
      return(list(word = w, organism = tg$organism[best], matched = tg$word[best]))
    }
  }
  NULL
}

.simitall_lexicon_binary_cue <- function(text) {
  grepl("case.?control|cases? and controls?|\\bbinary\\b|binomial|yes or no|yes/no|dead/alive|\\baffected\\b|\\bdisease\\b|logistic|firth",
        tolower(text))
}

# Which techniques/analyses does the text mention, and where would they route?
.simitall_detect_techniques <- function(text, drop = character()) {
  lex <- .simitall_lexicon()
  h <- .simitall_lexicon_hits(text, lex$technique_aliases)
  h <- h[!h$target %in% drop, , drop = FALSE]
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

# ---------------------------------------------------------------------------
# Table 3: settings (methods/design terms) and quantities; table 4: intent.

# Lower-case, digits for number words, no thousands separators.
.simitall_lexicon_numbers <- function(text) {
  x <- tolower(text)
  words <- c(one = 1, two = 2, three = 3, four = 4, five = 5, six = 6, seven = 7, eight = 8, nine = 9, ten = 10,
             eleven = 11, twelve = 12, fifteen = 15, twenty = 20, thirty = 30, forty = 40, fifty = 50, hundred = 100)
  for (w in names(words)) x <- gsub(paste0("\\b", w, "\\b"), words[[w]], x, perl = TRUE)
  while (grepl("(\\d),(\\d{3})\\b", x, perl = TRUE)) x <- gsub("(\\d),(\\d{3})\\b", "\\1\\2", x, perl = TRUE)
  x
}

# Settings for the routed analyses: primary first, then secondaries for fields
# still unset. Weak aliases count only when no strong alias set that field.
.simitall_detect_settings <- function(text, analyses) {
  lex <- .simitall_lexicon()
  analyses <- analyses[!is.na(analyses)]
  out <- list()
  if (!length(analyses)) return(out)
  h <- .simitall_lexicon_hits(text, lex$setting_aliases[lex$setting_aliases$analysis %in% analyses, , drop = FALSE])
  if (!nrow(h)) return(out)
  s <- lex$settings[as.integer(h$target), , drop = FALSE]
  s$strength <- h$strength
  for (a in analyses) {
    sa <- s[s$analysis == a, , drop = FALSE]
    for (f in unique(sa$field)) {
      if (!is.null(out[[f]])) next
      sf <- sa[sa$field == f, , drop = FALSE]
      if (any(sf$strength == "strong")) sf <- sf[sf$strength == "strong", , drop = FALSE]
      vals <- unique(sf$value)
      vals <- if (paste(a, f) %in% lex$set_fields) vals else vals[1]
      vals <- if (identical(vals, "true")) TRUE else if (identical(vals, "false")) FALSE else vals
      out[[f]] <- vals
    }
  }
  # short and long reads together make a hybrid assembly
  if ("genome_assembly" %in% analyses && grepl("illumina|short.read", tolower(text)) &&
      grepl("hifi|pacbio|nanopore|long.read", tolower(text))) out[["parameters.assembly_strategy"]] <- "hybrid"
  out
}

# Unit- or noun-anchored numbers ("h2 0.4", "30x", "10 generations", "4,000 cases").
.simitall_detect_quantities <- function(text, primary = NA_character_) {
  lex <- .simitall_lexicon()
  x <- .simitall_lexicon_numbers(text)
  q <- lex$quantities
  out <- list()
  for (g in unique(q$group)) {
    rows <- q[q$group == g, , drop = FALSE]
    applies <- vapply(strsplit(rows$analyses, ";", fixed = TRUE), function(a) "*" %in% a || (!is.na(primary) && primary %in% a), logical(1))
    rows <- rows[applies, , drop = FALSE]
    if (!nrow(rows)) next
    r <- rows[1, ]
    if (!is.null(out[[r$field]])) next
    m <- regmatches(x, regexec(r$pattern, x, perl = TRUE))[[1]]
    if (!length(m)) next
    caps <- m[-1]
    first <- caps[nzchar(caps)][1]
    val <- switch(sub(":.*", "", r$transform),
      int = as.integer(first),
      num = as.numeric(first),
      frac = { v <- as.numeric(sub("\\s*%", "", first)); if (grepl("%", first) || v > 1) v / 100 else v },
      bp = as.numeric(caps[1]) * c(bp = 1, kb = 1e3, mb = 1e6, gb = 1e9)[[caps[2]]],
      range = seq(as.integer(caps[1]), as.integer(caps[2])),
      sum = sum(as.integer(caps[nzchar(caps)])),
      count = length(regmatches(caps[1], gregexpr("\\d+", caps[1]))[[1]]),
      chr = paste0("chr", toupper(caps[1])),
      const = { v <- sub("^const:", "", r$transform); if (v %in% c("true", "false")) v == "true" else v },
      first)
    if (length(val) && !anyNA(val)) out[[r$field]] <- val
  }
  out
}

# Goal and data source from intent words.
.simitall_detect_intent <- function(text) {
  lex <- .simitall_lexicon()
  x <- tolower(text)
  hit <- vapply(lex$intent$cue, function(p) grepl(p, x, perl = TRUE), logical(1))
  sig <- lex$intent[hit & lex$intent$signal != "negation", , drop = FALSE]
  has <- function(s, strong = FALSE) any(sig$signal == s & (!strong | sig$strength == "strong"))
  list(simulate = has("simulate"), benchmark = has("benchmark"), analyze = has("analyze"),
       plan = has("plan"), provided_strong = has("provided", TRUE), provided = has("provided"),
       demo = has("demo"), only = has("only"), real_data = has("real_data"),
       clinical = has("clinical", TRUE) || (has("clinical") && has("provided")), correction = has("correction"),
       cues = sig$cue)
}

# Techniques used as *inputs* ("my GWAS hits", "300 DE genes", "eQTL results"):
# dropped when another analysis is requested.
.simitall_input_techniques <- function(text, hits) {
  if (!nrow(hits)) return(character())
  low <- tolower(.simitall_lexicon_norm(text))
  out <- character()
  for (i in seq_len(nrow(hits))) {
    after <- substr(low, hits$start[i] + nchar(hits$alias[i]), hits$start[i] + nchar(hits$alias[i]) + 25)
    if (grepl("(hits|results|genes|sumstats|summary statistics)$", tolower(hits$alias[i])) ||
        grepl("^\\s*(hits|results|output|summary statistics|sumstats|signals)\\b", after, perl = TRUE)) {
      out <- c(out, hits$target[i])
    }
  }
  unique(out)
}

# Techniques the user excludes: "no eQTL", "without DE", "don't run GWAS".
.simitall_negated_techniques <- function(text, hits) {
  if (!nrow(hits)) return(character())
  lex <- .simitall_lexicon()
  neg <- lex$intent$cue[lex$intent$signal == "negation"]
  low <- tolower(.simitall_lexicon_norm(text))
  out <- character()
  for (i in seq_len(nrow(hits))) {
    before <- substr(low, max(1, hits$start[i] - 30), hits$start[i] - 1)
    if (any(vapply(neg, function(p) grepl(p, before, perl = TRUE), logical(1)))) out <- c(out, hits$target[i])
  }
  unique(out)
}

.simitall_lexicon_quant_trait_cue <- function(text) {
  grepl("height|yield|weight|time|length|content|bmi|cholesterol|milk|protein|quantitative|continuous|flowering|days to|resistance|tolerance|growth|size|count|level|score|ldl|glucose|fat",
        tolower(text))
}

# Read a prompt into a StudySpec-shaped view using only the lexicon tables.
.simitall_lexicon_read_prompt <- function(text) {
  lex <- .simitall_lexicon()
  org <- .simitall_detect_organism(text)
  tq0 <- .simitall_detect_techniques(text)
  negated <- .simitall_negated_techniques(text, tq0$hits)
  inputs <- .simitall_input_techniques(text, tq0$hits)
  routed0 <- unique(stats::na.omit(c(tq0$primary_analysis, tq0$secondary_analyses)))
  if (length(inputs) && length(routed0) > 1) {
    t_routes <- lex$techniques$routes_to[match(inputs, lex$techniques$technique)]
    if (length(setdiff(routed0, c(t_routes, "quantitative_gwas", "binary_gwas")))) negated <- unique(c(negated, inputs))
  }
  tq <- if (length(negated)) .simitall_detect_techniques(text, drop = negated) else tq0
  it <- .simitall_detect_intent(text)

  analyses <- c(tq$primary_analysis, tq$secondary_analyses)
  analyses <- analyses[!is.na(analyses)]
  # "kinship", "PCs", "GRM" inside a GWAS request are corrections, not a separate analysis
  gwas_ids <- c("quantitative_gwas", "binary_gwas", "multitrait_longitudinal_gwas", "population_relatedness_gwas_calibration")
  ps_alias <- tolower(tq$hits$alias[tq$hits$target == "population_structure"])
  if (length(analyses) && analyses[1] %in% gwas_ids &&
      all(ps_alias %in% c("kinship", "kinship matrix", "grm", "relatedness", "pcs", "principal components"))) {
    analyses <- setdiff(analyses, "population_structure")
  }
  sim_only <- c("bulk_counts_simulation", "scrnaseq_simulation", "read_simulation", "cross_design", "phenotype_outcome_model")
  dual <- c("pedigree_inbreeding", "chipseq_atacseq", "time_series_rnaseq", "gxe_multiomics")

  goal <- if (it$simulate || it$benchmark) {
    if (!it$only && (it$benchmark || any(!analyses %in% c(sim_only, dual)))) "simulate_and_analyze" else "simulate"
  } else if (it$plan) "plan" else if (it$analyze || it$provided) "analyze" else
    if (length(analyses) && !analyses[1] %in% sim_only) "analyze" else NA_character_

  # An assay that alone is ambiguous is resolved by the goal (simulate scRNA-seq
  # -> scrnaseq_simulation). The simulation step comes first.
  route_goal <- if (identical(goal, "simulate_and_analyze")) "simulate" else goal
  if (!is.na(route_goal) && !length(tq$unsupported)) {
    r <- lex$intent_routes[lex$intent_routes$technique %in% tq$techniques & lex$intent_routes$goal == route_goal, , drop = FALSE]
    for (a in rev(unique(r$analysis))) {
      if (route_goal == "simulate" || !length(analyses)) analyses <- c(a, setdiff(analyses, a))
    }
  }
  if (it$only && length(analyses)) analyses <- analyses[1]
  primary <- if (length(analyses)) analyses[1] else NA_character_

  ds <- if (it$demo) "package_demo" else if (it$real_data && !it$provided_strong) "ask" else if (it$provided_strong) "provided" else
    if (it$provided && !(it$simulate || it$benchmark)) "provided" else
    if (it$simulate || it$benchmark) "synthetic" else "ask"

  v <- list(goal = if (is.na(goal)) NULL else goal, data_source = ds,
            organism = if (is.na(org$organism)) NULL else org$organism,
            primary_analysis = if (is.na(primary)) NULL else primary,
            secondary_analyses = analyses[-1], design = list(), parameters = list())
  vals <- c(.simitall_detect_settings(text, analyses), if (!it$correction) .simitall_detect_quantities(text, primary))
  for (f in names(vals)) {
    parts <- strsplit(f, ".", fixed = TRUE)[[1]]
    if (length(parts) == 1) {
      if (is.null(v[[f]])) v[[f]] <- vals[[f]]
    } else if (is.null(v[[parts[1]]][[parts[2]]])) {
      v[[parts[1]]][[parts[2]]] <- vals[[f]]
    }
  }
  # in donor-based analyses, counted people/samples are donors
  donor_based <- c("scrnaseq_simulation", "scrnaseq_analysis", "single_cell_eqtl", "bulk_eqtl", "splicing_qtl")
  if (!is.na(primary) && primary %in% donor_based && !is.null(v$design$population_size) && is.null(v$design$donors)) {
    v$design$donors <- v$design$population_size
    v$design$population_size <- NULL
  }
  if (identical(primary, "quantitative_gwas") && .simitall_lexicon_quant_trait_cue(text)) v$parameters$trait_type <- "continuous"

  flags <- c(.simitall_lexicon_checks(text, org$organism), .simitall_lexicon_range_checks(text))
  if (!is.null(org$assumed)) {
    flags <- c(flags, sprintf("assumed '%s' means %s (%s); correct me if not", org$assumed$word,
                              lex$organisms$scientific_name[lex$organisms$organism == org$organism][1], org$organism))
  }
  # an assay named next to other analyses but with no say on what to do with it
  if (length(analyses) && length(tq$unresolved_assays)) {
    for (a in tq$unresolved_assays) {
      lab <- lex$techniques$label[match(a, lex$techniques$technique)]
      flags <- c(flags, sprintf("%s was mentioned without saying what to do with it (for example simulate it or test differences): which analysis?", lab))
    }
    v$unresolved_assays <- tq$unresolved_assays
  }
  # F1 siblings come from one cross; selfing for generations gives S-lines, not F1s
  low_text <- tolower(text)
  if (grepl("\\bf1\\b", low_text) && grepl("self", low_text) && grepl("\\d+\\s+gen|gen\\w*\\s+of\\s+self|for\\s+\\d+\\s+gen", low_text)) {
    flags <- c(flags, "F1 siblings come from a single cross of two parents, but several generations of selfing give inbred (S) lines, not F1s: which population should the next step use?")
  }
  if (grepl("ignore (?:your|all|previous|the) (?:rules|instructions)", tolower(text), perl = TRUE)) {
    flags <- c(flags, "request asks to ignore the agent's rules")
  }
  if (it$clinical) flags <- c(flags, "clinical or individual risk decision requested: out of scope")
  if (it$real_data && !it$provided_strong) flags <- c(flags, "real data requested but not supplied: ask for the files; do not substitute synthetic data")
  if (it$correction) flags <- c(flags, "numbers changed mid-request: confirm which value to use")
  blocking <- grepl("not possible|conflict|ignore the agent|out of scope|real data requested|out of range", flags)
  v$execution_status <- if (length(tq$unsupported)) "unsupported" else if (any(blocking) || is.na(primary)) "needs_clarification" else "supported"
  v$do_not_run <- length(tq$unsupported) > 0 || any(blocking)
  v$flags <- flags
  v$negated <- negated
  v$unsupported_techniques <- tq$unsupported
  # an assay that alone is ambiguous ("RNA-seq.") -> the analyses to choose from
  if (is.na(primary) && length(tq$candidates) && !length(tq$unsupported)) {
    v$candidates <- tq$candidates
    amb <- tq$unresolved_assays[1] %||% tq$techniques[1]
    v$candidate_label <- lex$techniques$label[match(amb, lex$techniques$technique)]
  }
  v$defaulted <- character()
  v
}

# Values outside their valid range are flagged and never read as settings.
.simitall_lexicon_range_checks <- function(text) {
  x <- .simitall_lexicon_numbers(text)
  flags <- character()
  h <- regmatches(x, regexec("(?:h2|heritability)\\s*(?:of|=|:|is)?\\s*(\\d+(?:\\.\\d+)?)", x, perl = TRUE))[[1]]
  if (length(h) && as.numeric(h[2]) > 1) flags <- c(flags, sprintf("heritability %s is out of range [0, 1]", h[2]))
  p <- as.numeric(regmatches(x, gregexpr("\\d+(?:\\.\\d+)?(?=\\s*%)", x, perl = TRUE))[[1]])
  if (any(p > 100)) flags <- c(flags, sprintf("percentage %s%% is out of range", max(p)))
  if (grepl("(^|\\s)-\\s*\\d+(?:\\.\\d+)?\\s*x\\b|negative (?:coverage|replicates|samples)", x, perl = TRUE)) {
    flags <- c(flags, "negative value is out of range")
  }
  flags
}
