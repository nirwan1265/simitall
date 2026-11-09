# Optional trained prompt model (inst/python/prompt_model.py, training/).
#
# Rules first: the lexicon reader (.simitall_lexicon_read_prompt) decides, and
# the model only fills fields the lexicon left empty, reorders analyses the
# lexicon already found, or removes analyses it marks as inputs/excluded.
# Every number the model adds is a span copied from the prompt.
# Read-only; not yet wired into simitall_ask().

.simitall_model_fields <- c(
  DONORS = "design.donors", FOUNDERS = "design.founder_count", GENS = "parameters.generations",
  SELFGEN = "design.selfing_generations", BCGEN = "design.backcross_generations",
  TECHREPS = "design.technical_replicates", TIMEPTS = "design.time_points", CASES = "design.cases",
  CONTROLS = "design.controls", ENVS = "design.environments", BATCHES = "design.batches",
  NCAUSAL = "parameters.n_causal", CELLS = "parameters.cells_per_donor", NCELLTYPES = "design.cell_types",
  NCOND = "design.conditions", COV = "parameters.coverage", H2 = "parameters.heritability",
  POP = "@population_size", REPS = "@replicates"
)

# Run the Python model on one prompt (needs torch + transformers in `python`).
.simitall_model_dir <- function() {
  d <- Sys.getenv("SIMITALL_PROMPT_MODEL")
  if (nzchar(d)) return(d)
  d <- system.file("models", "simitall-tagger", package = "simitall")
  if (!nzchar(d)) d <- file.path("inst", "models", "simitall-tagger")
  d
}

.simitall_model_predict <- function(text, model_dir = .simitall_model_dir(),
                                    python = Sys.getenv("SIMITALL_PYTHON", "python3")) {
  script <- system.file("python", "prompt_model.py", package = "simitall")
  if (!nzchar(script)) script <- file.path("inst", "python", "prompt_model.py")
  if (!nzchar(model_dir) || !file.exists(file.path(model_dir, "tagger.pt"))) return(NULL)
  out <- tryCatch(suppressWarnings(system2(python, c(shQuote(script), "--model", shQuote(model_dir), "--text", shQuote(text)),
                                           stdout = TRUE, stderr = FALSE)), error = function(e) character())
  if (!length(out)) return(NULL)
  tryCatch(jsonlite::fromJSON(utils::tail(out, 1), simplifyVector = FALSE), error = function(e) NULL)
}

.simitall_model_number <- function(x) {
  v <- suppressWarnings(as.numeric(regmatches(.simitall_lexicon_numbers(x), regexpr("\\d+(?:\\.\\d+)?", .simitall_lexicon_numbers(x), perl = TRUE))))
  if (!length(v) || is.na(v)) return(NULL)
  if (v == round(v)) as.integer(v) else v
}

.simitall_model_span_analysis <- function(span, goal = NULL) {
  tq <- .simitall_detect_techniques(span)
  if (!is.na(tq$primary_analysis)) return(tq$primary_analysis)
  lex <- .simitall_lexicon()
  g <- if (identical(goal, "simulate_and_analyze")) "simulate" else goal %||% ""
  r <- lex$intent_routes$analysis[lex$intent_routes$technique %in% tq$techniques & lex$intent_routes$goal == g]
  if (length(r)) r[1] else NA_character_
}

.simitall_model_merge <- function(view, pred, text, min_conf = 0.6) {
  if (is.null(pred) || !length(pred)) return(view)
  sp <- pred$spans %||% list()
  tag <- vapply(sp, function(s) s[[1]], "")
  txt <- vapply(sp, function(s) as.character(s[[2]]), "")
  pos <- vapply(sp, function(s) as.integer(s[[3]]), 1L)
  blocked <- isTRUE(view$do_not_run) || identical(view$execution_status, "unsupported")
  if (blocked) return(view)
  toks <- tolower(unlist(regmatches(text, gregexpr("\\d+(?:[.,]\\d+)*|[A-Za-z]+(?:[-'/][A-Za-z0-9]+)*|[^\\sA-Za-z\\d]",
                                                   text, perl = TRUE))))
  before <- function(i, k = 3) toks[seq(max(1, pos[i] + 1 - k), pos[i])]
  after <- function(i, k = 3) {
    n <- length(strsplit(txt[i], " ", fixed = TRUE)[[1]])
    toks[seq(pos[i] + n + 1, min(length(toks), pos[i] + n + k))]
  }   # the lexicon decided to stop; the model adds nothing
  numbers_ok <- !any(grepl("numbers changed", unlist(view$flags)))

  # sentence-level fallbacks
  if (is.null(view$goal) && !is.null(pred$goal) && pred$goal[[2]] >= min_conf && pred$goal[[1]] != "none") {
    view$goal <- pred$goal[[1]]
  }
  analyses <- c(view$primary_analysis, unlist(view$secondary_analyses))
  if (!length(analyses) && !blocked && !is.null(pred$primary) && pred$primary[[2]] >= min_conf &&
      !pred$primary[[1]] %in% c("none", "unsupported")) {
    analyses <- pred$primary[[1]]
  }
  # roles: inputs / exclusions leave; the PRIMARY-tagged analysis moves first
  if (!blocked && length(analyses)) {
    for (i in which(tag %in% c("INPUT", "NEG"))) {
      # an input must be marked as one in the text ("my GWAS hits", "these eQTL results")
      if (tag[i] == "INPUT" && pos[i] < length(toks) &&
          !any(after(i) %in% c("hits", "results", "genes", "peaks", "output", "list", "signals", "loci", "sumstats")) &&
          !any(before(i) %in% c("my", "our", "these", "existing", "previous", "published"))) next
      a <- .simitall_model_span_analysis(txt[i], view$goal)
      if (!is.na(a) && a %in% analyses && length(analyses) > 1) analyses <- setdiff(analyses, a)
    }
    # reorder only when the model names exactly one primary analysis, and never move a
    # simulation step from the front of "simulate X, then Y" (reviewed convention) or
    # replace a specialised GWAS with a generic one
    pa <- unique(stats::na.omit(vapply(txt[tag == "PRIMARY"], .simitall_model_span_analysis, "", goal = view$goal)))
    keep_first <- (view$goal %||% "") %in% c("simulate", "simulate_and_analyze") &&
      analyses[1] %in% c("bulk_counts_simulation", "scrnaseq_simulation", "read_simulation", "pedigree_inbreeding",
                         "chipseq_atacseq", "time_series_rnaseq", "gxe_multiomics")
    special_gwas <- analyses[1] %in% c("binary_gwas", "multitrait_longitudinal_gwas", "rare_variant_and_gene_based_gwas",
                                       "population_relatedness_gwas_calibration")
    # "X, then Y": the text already gives the order, which the lexicon keeps
    sequenced <- grepl("\\bthen\\b|\\bafterwards\\b|\\bfollowed by\\b|\\bafter that\\b", tolower(text))
    if (length(pa) == 1 && pa %in% analyses && !keep_first && !sequenced && !(special_gwas && pa == "quantitative_gwas")) {
      analyses <- c(pa, setdiff(analyses, pa))
    }
  }
  if (length(analyses)) {
    view$primary_analysis <- analyses[1]
    view$secondary_analyses <- analyses[-1]
  }
  primary <- view$primary_analysis %||% NA_character_

  set_if_empty <- function(field, value) {
    parts <- strsplit(field, ".", fixed = TRUE)[[1]]
    if (is.null(value) || !is.null(view[[parts[1]]][[parts[2]]])) return()
    view[[parts[1]]][[parts[2]]] <<- value
  }
  # numbers: the first span per field; split groups joined by "and" are summed
  for (t in if (numbers_ok) intersect(names(.simitall_model_fields), unique(tag)) else character()) {
    if (t == "GENS" && !grepl("\\bgen", tolower(text))) next   # "10 rounds" is not 10 generations
    idx <- which(tag == t)
    # numbers glued to a chromosome name ("chr10", "chromosome 10") are not counts
    idx <- idx[!vapply(idx, function(i) pos[i] >= 1 && toks[pos[i]] %in% c("chr", "chromosome", "chrom"), logical(1))]
    if (!length(idx)) next
    vals <- lapply(txt[idx], .simitall_model_number)
    ok <- !vapply(vals, is.null, logical(1))
    if (!any(ok)) next
    idx <- idx[ok]; vals <- unlist(vals[ok])
    value <- vals[1]
    if (t %in% c("DONORS", "POP") && length(vals) > 1) {
      between <- toks[seq(pos[idx[1]] + 2, max(pos[idx[2]], pos[idx[1]] + 2))]
      if (any(between %in% c("and", "+", "plus")) && pos[idx[2]] - pos[idx[1]] <= 6) value <- sum(vals[1:2])
    }
    if (t == "H2" && (value < 0 || value > 1)) next
    if (t == "COV" && value <= 0) next
    field <- .simitall_model_fields[[t]]
    if (field == "@population_size") field <- if (identical(primary, "pedigree_inbreeding")) "parameters.population_size" else "design.population_size"
    if (field == "@replicates") field <- if (identical(primary, "chipseq_atacseq")) "parameters.replicates" else "design.biological_replicates_per_condition"
    set_if_empty(field, value)
  }
  # named lists
  cond <- unique(txt[tag == "COND"])
  if (length(cond) >= 2) set_if_empty("design.conditions", as.list(cond))
  ct <- unique(txt[tag == "CELLTYPE"])
  if (length(ct) >= 2) set_if_empty("design.cell_types", as.list(ct))
  if (is.null(view$organism) && any(tag == "ORG")) {
    o <- .simitall_detect_organism(txt[tag == "ORG"][1])$organism
    if (!is.na(o)) view$organism <- o
  }
  view
}

# Read a request into a StudySpec-shaped view: lexicon tables first, then the
# trained model (inst/models/simitall-tagger) when Python with torch is
# available. Falls back to the lexicon alone, silently, when it is not.
.simitall_read_prompt <- function(question, use_model = TRUE) {
  view <- .simitall_lexicon_read_prompt(question)
  view$reader <- "lexicon"
  if (isTRUE(use_model)) {
    pred <- .simitall_model_predict(question)
    if (!is.null(pred)) {
      view <- .simitall_model_merge(view, pred, question)
      view$reader <- "lexicon+model"
    }
  }
  view
}
