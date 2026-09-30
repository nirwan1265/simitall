.simitall_agent_knowledge_dir <- function(knowledge_dir = NULL) {
  if (!is.null(knowledge_dir)) {
    if (!dir.exists(knowledge_dir)) {
      stop("knowledge_dir does not exist: ", knowledge_dir, call. = FALSE)
    }
    return(normalizePath(knowledge_dir, mustWork = TRUE))
  }

  installed <- system.file("agent", "knowledge", package = "simitall")
  if (nzchar(installed) && dir.exists(installed)) return(installed)

  source_tree <- file.path("inst", "agent", "knowledge")
  if (dir.exists(source_tree)) return(normalizePath(source_tree, mustWork = TRUE))

  stop(
    "The simitall agent knowledge base was not found. Reinstall simitall or ",
      "supply knowledge_dir explicitly.",
    call. = FALSE
  )
}

.simitall_agent_aliases <- function() {
  installed <- system.file("agent", "aliases.tsv", package = "simitall")
  source_tree <- file.path("inst", "agent", "aliases.tsv")
  path <- if (nzchar(installed) && file.exists(installed)) installed else source_tree
  if (!file.exists(path)) {
    return(data.frame(alias = character(), canonical_term = character()))
  }
  aliases <- utils::read.delim(
    path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  required <- c("alias", "canonical_term")
  if (!all(required %in% names(aliases))) {
    stop("Agent alias table must contain: ", paste(required, collapse = ", "))
  }
  aliases$alias <- tolower(trimws(aliases$alias))
  aliases$canonical_term <- tolower(trimws(aliases$canonical_term))
  aliases[nzchar(aliases$alias) & nzchar(aliases$canonical_term), required,
          drop = FALSE]
}

.simitall_agent_terms <- function(x) {
  x <- tolower(enc2utf8(paste(x, collapse = " ")))
  x <- gsub("[^a-z0-9]+", " ", x)
  terms <- unlist(strsplit(trimws(x), "\\s+"), use.names = FALSE)
  # Common request phrasing otherwise overwhelms the few scientific terms that
  # should drive retrieval (for example, "I have 500 lines and want to...").
  stop_terms <- c(
    "a", "an", "and", "are", "as", "at", "be", "by", "can", "do", "does",
    "for", "from", "have", "how", "i", "if", "in", "is", "it", "of", "on",
    "or", "should", "that", "the", "this", "to", "want", "what", "with",
    "would", "you", "your", "first", "lines", "line", "plan", "simulate",
    "simulation", "500"
  )
  terms <- terms[nchar(terms) >= 2L & !terms %in% stop_terms]
  # Preserve the user's wording, then append audited canonical terms from the
  # versioned alias table. Phrase aliases are matched before retrieval.
  aliases <- .simitall_agent_aliases()
  padded_query <- paste0(" ", x, " ")
  matched <- vapply(
    aliases$alias,
    function(alias) grepl(paste0(" ", alias, " "), padded_query, fixed = TRUE),
    logical(1)
  )
  canonical <- aliases$canonical_term[matched]
  unique(c(terms, canonical))
}

.simitall_agent_has_term <- function(text, term) {
  grepl(paste0("\\b", term, "\\b"), text, perl = TRUE)
}

.simitall_agent_term_count <- function(text, term) {
  matches <- gregexpr(paste0("\\b", term, "\\b"), text, perl = TRUE)[[1L]]
  if (identical(matches[1L], -1L)) 0L else length(matches)
}

.simitall_agent_excerpt <- function(lines, terms, width = 280L) {
  if (!length(lines)) return("")
  # Metadata is useful for ranking but not as a user-facing evidence excerpt.
  if (identical(lines[1L], "---")) {
    closing <- which(lines[-1L] == "---")
    if (length(closing)) lines <- lines[(closing[1L] + 2L):length(lines)]
  }
  hits <- vapply(
    lines,
    function(line) {
      line <- tolower(line)
      any(vapply(terms, function(term) .simitall_agent_has_term(line, term), logical(1)))
    },
    logical(1)
  )
  matching_lines <- which(hits)
  index <- if (length(matching_lines)) matching_lines[1L] else 1L
  text <- trimws(lines[index])
  if (is.na(text)) text <- ""
  if (nchar(text) > width) paste0(substr(text, 1L, width - 3L), "...") else text
}

.simitall_agent_evidence_lines <- function(lines, terms, n_lines = 3L) {
  if (identical(lines[1L], "---")) {
    closing <- which(lines[-1L] == "---")
    if (length(closing)) lines <- lines[(closing[1L] + 2L):length(lines)]
  }
  lines <- trimws(lines)
  usable <- nzchar(lines) &
    !grepl("^#{1,6}\\s", lines) &
    !grepl("^```", lines) &
    !grepl("^[(){}]+$", lines) &
    !grepl("<-|=|\\(|^library\\(|^install\\.packages", lines)
  lines <- lines[usable]
  if (!length(lines)) return(character())
  scores <- vapply(lines, function(line) {
    sum(vapply(terms, function(term) .simitall_agent_has_term(tolower(line), term), logical(1)))
  }, numeric(1))
  selected <- lines[order(-scores, seq_along(lines))]
  selected <- unique(selected[scores[order(-scores, seq_along(lines))] > 0])
  if (!length(selected)) selected <- lines[seq_len(min(n_lines, length(lines)))]
  substr(selected[seq_len(min(n_lines, length(selected)))], 1L, 400L)
}

.simitall_agent_title <- function(lines, fallback) {
  yaml_title <- grep("^title:\\s*", lines, value = TRUE)
  if (length(yaml_title)) {
    title <- sub("^title:\\s*", "", yaml_title[1L])
    title <- gsub("^['\\\"]|['\\\"]$", "", trimws(title))
    if (nzchar(title)) return(title)
  }
  headings <- grep("^#{1,2}\\s+", lines, value = TRUE)
  if (length(headings)) return(trimws(sub("^#{1,2}\\s+", "", headings[1L])))
  fallback
}

.simitall_agent_operator <- function(lhs, rhs) {
  # Some metadata fields, such as retrieved query terms, are vectors rather
  # than scalars. Treat only an entirely missing value as absent.
  if (is.null(lhs) || length(lhs) == 0L || all(is.na(lhs))) rhs else lhs
}

`%||%` <- .simitall_agent_operator

.simitall_agent_documents <- function(knowledge_dir) {
  paths <- list.files(
    knowledge_dir,
    pattern = "\\.(Rmd|rmd|md|MD)$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(paths)) stop("No Markdown or R Markdown files found in knowledge_dir")
  lapply(paths, function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    list(
      path = normalizePath(path, mustWork = TRUE),
      relative_path = substring(
        normalizePath(path, mustWork = TRUE),
        nchar(normalizePath(knowledge_dir, mustWork = TRUE)) + 2L
      ),
      category = basename(dirname(path)),
      title = .simitall_agent_title(lines, tools::file_path_sans_ext(basename(path))),
      lines = lines,
      text = paste(lines, collapse = "\n")
    )
  })
}

.simitall_agent_category_priority <- function(category) {
  priorities <- c(
    validation_contracts = 14,
    workflows = 12,
    techniques = 10,
    species = 8,
    tools = 5,
    terminology = 2,
    evaluations = 1,
    papers = 0
  )
  if (!category %in% names(priorities)) return(0)
  unname(priorities[[category]])
}

.simitall_agent_intent_bonus <- function(terms, text) {
  intent_terms <- list(
    breeding = c("breeding", "nam", "magic", "ril", "nil", "f1", "f2", "backcross", "selfing", "cross"),
    gwas = c("gwas", "qtl", "association"),
    rnaseq = c("rna", "rnaseq", "eqtl", "expression", "transcript"),
    genomic_selection = c("genomic", "selection", "gblup", "rrblup")
  )
  text <- tolower(text)
  broad_bonus <- sum(vapply(intent_terms, function(keywords) {
    requested <- any(keywords %in% terms)
    documented <- any(vapply(keywords, function(term) .simitall_agent_has_term(text, term), logical(1)))
    as.integer(requested && documented) * 12L
  }, numeric(1)))
  # Named mapping-population designs are decisive constraints rather than
  # generic vocabulary. Preserve them even when a multi-topic request also
  # contains many high-frequency terms such as GWAS and RNA-seq.
  design_terms <- intersect(terms, c("nam", "magic", "ril", "nil", "f1", "f2", "dh"))
  design_bonus <- sum(vapply(design_terms, function(term) {
    as.integer(.simitall_agent_has_term(text, term)) * 25L
  }, numeric(1)))
  broad_bonus + design_bonus
}

.simitall_agent_topic_penalty <- function(terms, title) {
  title <- tolower(title)
  penalties <- 0L
  asks_single_cell <- any(c("single", "cell", "scrna", "singlecell") %in% terms)
  asks_chipseq <- any(c("chip", "chipseq", "binding", "histone") %in% terms)
  if (grepl("single.cell", title) && !asks_single_cell) penalties <- penalties + 30L
  if (grepl("chip.seq", title) && !asks_chipseq) penalties <- penalties + 30L
  penalties
}

#' Search the simitall agent knowledge base
#'
#' Perform local, transparent lexical retrieval across the versioned knowledge
#' documents stored under `inst/agent/knowledge`. This is the grounding layer
#' for [simitall_ask()] and is also useful on its own when preparing a workflow
#' or checking which documentation supports a proposed analysis. No network
#' request, language model request, or code execution occurs.
#'
#' @param query Character string describing a scientific question, method, or
#'   workflow.
#' @param n_results Maximum number of ranked documents to return.
#' @param knowledge_dir Optional alternate knowledge directory. By default,
#'   uses the copy installed with `simitall`; a source checkout is also found.
#'
#' @return A data frame with document title, category, score, excerpt, and
#'   absolute file path. Attributes `query_terms` and `knowledge_dir` record
#'   the retrieval inputs.
#' @examples
#' search_simitall_knowledge("maize NAM population recombination map")
#' @export
search_simitall_knowledge <- function(query, n_results = 6L, knowledge_dir = NULL) {
  if (!is.character(query) || length(query) != 1L || is.na(query) || !nzchar(trimws(query))) {
    stop("query must be one non-empty character string", call. = FALSE)
  }
  if (length(n_results) != 1L || is.na(n_results) || n_results < 1L) {
    stop("n_results must be a positive integer", call. = FALSE)
  }
  n_results <- as.integer(n_results)
  knowledge_dir <- .simitall_agent_knowledge_dir(knowledge_dir)
  terms <- .simitall_agent_terms(query)
  docs <- .simitall_agent_documents(knowledge_dir)

  results <- lapply(docs, function(doc) {
    body <- tolower(doc$text)
    title <- tolower(doc$title)
    filename <- tolower(doc$relative_path)
    term_scores <- vapply(terms, function(term) {
      occurrences <- .simitall_agent_term_count(body, term)
      # A review paper can contain a term dozens of times. Cap body matches so
      # package workflows and method cards outrank repeated narrative text.
      body_score <- min(3L, occurrences)
      title_bonus <- 5L * as.integer(.simitall_agent_has_term(title, term))
      path_bonus <- 3L * as.integer(.simitall_agent_has_term(filename, term))
      body_score + title_bonus + path_bonus
    }, numeric(1))
    match_score <- sum(term_scores)
    category_bonus <- if (match_score > 0) {
      .simitall_agent_category_priority(doc$category)
    } else {
      0
    }
    intent_bonus <- .simitall_agent_intent_bonus(terms, doc$text)
    topic_penalty <- .simitall_agent_topic_penalty(terms, doc$title)
    data.frame(
      title = doc$title,
      category = doc$category,
      score = match_score + category_bonus + intent_bonus - topic_penalty,
      excerpt = .simitall_agent_excerpt(doc$lines, terms),
      path = doc$path,
      stringsAsFactors = FALSE
    )
  })
  results <- do.call(rbind, results)
  # Prefer operational package knowledge over background papers when both
  # discuss the same subject. Documents with no query match remain available
  # only as a fallback for very short or unfamiliar questions.
  operational <- results$category %in% c("validation_contracts", "workflows", "techniques", "species", "tools") &
    results$score > 0
  results <- results[order(!operational, -results$score, results$path), , drop = FALSE]
  results <- results[seq_len(min(n_results, nrow(results))), , drop = FALSE]
  rownames(results) <- NULL
  attr(results, "query_terms") <- terms
  attr(results, "knowledge_dir") <- knowledge_dir
  results
}

.simitall_agent_context <- function(retrieval, max_chars = 18000L) {
  pieces <- character()
  if (!nrow(retrieval)) return("")
  terms <- attr(retrieval, "query_terms") %||% character()
  # Small local models respond more reliably to focused evidence cards than to
  # entire Rmd files, which otherwise tempt them to reproduce document text.
  chars_per_source <- max(400L, floor(max_chars / nrow(retrieval)))
  for (i in seq_len(nrow(retrieval))) {
    path <- retrieval$path[i]
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    evidence <- .simitall_agent_evidence_lines(lines, terms)
    block <- paste0(
      "SOURCE: ", path, "\n",
      "CATEGORY: ", retrieval$category[i], "\n",
      "TITLE: ", retrieval$title[i], "\n",
      "EVIDENCE:\n- ", paste(evidence, collapse = "\n- ")
    )
    block <- substr(block, 1L, chars_per_source)
    pieces <- c(pieces, block)
  }
  paste(pieces, collapse = "\n\n---\n\n")
}

#' List the safe functions available to the simitall agent
#'
#' Return the initial allowlist of high-level `simitall` operations that a
#' future execution-enabled agent may propose. The current [simitall_ask()]
#' implementation is deliberately read-only: it answers and plans but never
#' invokes these functions. Keeping the registry explicit makes later tool
#' execution auditable and avoids exposing an arbitrary shell to a model.
#'
#' @return A data frame containing tool names, domains, descriptions, and
#'   execution status.
#' @examples
#' simitall_agent_tools()
#' @export
simitall_agent_tools <- function() {
  data.frame(
    function_name = c(
      "simulate_genome", "simulate_annotations", "simulate_breeding",
      "generate_random_haplotype_panel", "simupop_api", "simulate_gwas_cohort",
      "analyze_gwas", "benchmark_gwas", "plot_gwas_results", "fit_genomic_model",
      "run_genomic_selection", "simulate_rnaseq_from_gwas",
      "simulate_rnaseq_experiment", "simulate_rnaseq_reads", "benchmark_eqtl",
      "plot_rnaseq_eqtl_results", "simulate_scrnaseq_from_gwas",
      "simulate_scrnaseq_experiment", "simulate_celltype_eqtls",
      "benchmark_celltype_eqtls", "plot_scrnaseq_results", "simulate_chipseq",
      "simulate_parental_effects", "run_grid_both", "run_unicycler_grid",
      "run_quast_grid"
    ),
    domain = c(
      "genome", "annotation", "population", "population", "population",
      "GWAS", "GWAS", "GWAS", "GWAS", "genomic selection",
      "genomic selection", "bulk RNA-seq", "bulk RNA-seq", "bulk RNA-seq",
      "bulk RNA-seq", "bulk RNA-seq", "single-cell RNA-seq",
      "single-cell RNA-seq", "single-cell RNA-seq", "single-cell RNA-seq",
      "single-cell RNA-seq", "ChIP-seq", "pedigree phenotypes", "sequencing", "assembly",
      "assembly evaluation"
    ),
    description = c(
      "Create a reference-derived or random genome with truth outputs.",
      "Add genes, regulatory elements, repeats, plasmids, and annotations.",
      "Create controlled crosses and mapping populations.",
      "Create a synthetic founder FASTA haplotype panel.",
      "Bridge to selected advanced SimuPOP population operations.",
      "Generate a genotype, phenotype, and truth-set GWAS cohort.",
      "Run a mixed-model genome-wide association scan.",
      "Compare GWAS results with known causal-variant truth.",
      "Create GWAS result figures.",
      "Fit genomic prediction models such as GBLUP or ridge regression.",
      "Select parents and iterate genomic-selection cycles.",
      "Generate bulk expression from a simulated GWAS cohort.",
      "Generate an independent bulk RNA-seq experiment.",
      "Generate RNA-seq reads from simulated expression.",
      "Compare simulated eQTL results with known eQTL truth.",
      "Create bulk RNA-seq and eQTL figures.",
      "Generate population-aware single-cell data and cell-type eQTL truth.",
      "Generate an independent single-cell RNA-seq experiment.",
      "Simulate cell-type-specific eQTL effects.",
      "Compare cell-type eQTL results with known truth.",
      "Create single-cell RNA-seq result figures.",
      "Generate binding peaks and reads from annotated regulatory elements.",
      "Simulate direct, maternal, and paternal effects from a pedigree.",
      "Create an Illumina by long-read sequencing coverage grid.",
      "Run an allowlisted hybrid-assembly grid.",
      "Evaluate generated assemblies against a known reference."
    ),
    execution_status = "planned; explicit user approval required",
    stringsAsFactors = FALSE
  )
}

.simitall_openai_response_text <- function(body) {
  if (!is.null(body$output_text) && nzchar(body$output_text)) return(body$output_text)
  output <- body$output %||% list()
  texts <- unlist(lapply(output, function(item) {
    content <- item$content %||% list()
    unlist(lapply(content, function(part) {
      if (identical(part$type, "output_text")) part$text else NULL
    }), use.names = FALSE)
  }), use.names = FALSE)
  paste(texts, collapse = "\n")
}

.simitall_openai_request <- function(question, instructions, model, api_key) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop(
      "Package 'httr2' is required for OpenAI requests. Install it with ",
      "install.packages('httr2').",
      call. = FALSE
    )
  }
  response <- httr2::request("https://api.openai.com/v1/responses") |>
    httr2::req_headers(Authorization = paste("Bearer", api_key)) |>
    httr2::req_body_json(list(
      model = model,
      instructions = instructions,
      input = question
    )) |>
    # Parse API failures below so users see the actionable OpenAI error code.
    httr2::req_error(is_error = function(response) FALSE) |>
    httr2::req_perform()
  body <- tryCatch(
    httr2::resp_body_json(response, simplifyVector = FALSE),
    error = function(error) list()
  )
  status <- httr2::resp_status(response)
  if (identical(status, 401L)) {
    stop(
      "OpenAI rejected OPENAI_API_KEY (HTTP 401). Use an API key created at ",
        "platform.openai.com, not a ChatGPT/Codex sign-in token. Confirm that ",
        "the key was loaded after restarting R with nzchar(Sys.getenv('OPENAI_API_KEY')). ",
        "Do not print or share the key.",
      call. = FALSE
    )
  }
  if (status >= 400L) {
    api_error <- body$error %||% list()
    code <- api_error$code %||% api_error$type %||% "unknown_error"
    message <- api_error$message %||% "No error message was returned by OpenAI."
    guidance <- if (identical(status, 429L)) {
      paste0(
        " Check the API Billing and Limits pages for prepaid credits or a spend ",
        "limit before retrying; a first-request 429 is commonly a quota issue."
      )
    } else ""
    stop(
      "OpenAI API request failed (HTTP ", status, ", ", code, "): ",
        message, guidance,
      call. = FALSE
    )
  }
  answer <- .simitall_openai_response_text(body)
  if (!nzchar(answer)) stop("OpenAI returned no text output", call. = FALSE)
  answer
}

.simitall_ollama_request <- function(instructions, model, endpoint) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package 'httr2' is required for Ollama requests. Install it with install.packages('httr2').")
  }
  response <- tryCatch(
    httr2::request(endpoint) |>
      httr2::req_body_json(list(
        model = model,
        prompt = instructions,
        stream = FALSE,
        options = list(temperature = 0.2)
      )) |>
      httr2::req_perform(),
    error = function(error) {
      stop(
        "Could not reach local Ollama at ", endpoint, ". Install and start Ollama, ",
          "then download the selected model (for example: ollama pull ", model, "). ",
          "No API key or paid account is required for the local Ollama backend.",
        call. = FALSE
      )
    }
  )
  body <- httr2::resp_body_json(response, simplifyVector = FALSE)
  answer <- body$response %||% ""
  if (!nzchar(answer)) stop("Local Ollama returned no text output", call. = FALSE)
  answer
}

.simitall_agent_allowed_tools_text <- function() {
  tools <- simitall_agent_tools()
  paste0(
    "- `", tools$function_name, "()`: ", tools$description,
    collapse = "\n"
  )
}

.simitall_agent_safe_fallback <- function(question, retrieval, violations) {
  paste(
    "Safety fallback: the local model response was withheld because it included",
    "unsupported code or package claims (", paste(violations, collapse = ", "), ").",
    "Use this evidence briefing instead, or refine the question.",
    "",
    .simitall_local_answer(question, retrieval),
    sep = "\n"
  )
}

.simitall_agent_validate_answer <- function(answer, question, retrieval, allow_code) {
  prohibited <- c(
    "```" = "code block",
    "install\\.packages\\s*\\(" = "package installation",
    "(?:library|require)\\s*\\(" = "package loading",
    "\\b(?:gapit|tassel|maize_gwas_data)\\b" = "unsupported external tool or data object"
  )
  if (isTRUE(allow_code)) prohibited <- prohibited[names(prohibited) != "```"]
  hits <- names(prohibited)[vapply(
    names(prohibited),
    function(pattern) grepl(pattern, answer, ignore.case = TRUE, perl = TRUE),
    logical(1)
  )]
  # Function-like text is only safe when it names an allowlisted high-level
  # simitall operation. This catches convincing-but-nonexistent suggestions
  # such as benchmark_eqtl() before they reach a user.
  calls <- regmatches(
    answer,
    gregexpr("\\b[A-Za-z][A-Za-z0-9_]*\\s*\\(", answer, perl = TRUE)
  )[[1L]]
  calls <- sub("\\s*\\($", "", calls)
  calls <- unique(calls[nzchar(calls)])
  allowed <- simitall_agent_tools()$function_name
  unsupported_calls <- setdiff(calls, allowed)
  if (length(unsupported_calls)) {
    hits <- c(hits, paste0("unsupported function: ", unsupported_calls))
  }
  if (!length(hits)) return(answer)
  labels <- unname(prohibited[intersect(hits, names(prohibited))])
  labels <- c(labels, setdiff(hits, names(prohibited)))
  .simitall_agent_safe_fallback(question, retrieval, labels)
}

.simitall_agent_source_file <- function(...) {
  relative_path <- file.path(...)
  here <- normalizePath(getwd(), mustWork = TRUE)
  for (i in 0:6) {
    candidate <- file.path(here, relative_path)
    if (file.exists(candidate)) return(candidate)
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  NULL
}

.simitall_agent_code_recipe <- function(question, preflight) {
  terms <- .simitall_agent_terms(question)
  is_nam <- "nam" %in% terms
  is_human_irf6 <- "irf6" %in% terms &&
    any(c("human", "pedigree", "liability", "cleft") %in% terms)
  is_human_inbreeding <- .simitall_agent_species(question) == "human" &&
    grepl("inbreed|inbred|consanguin|first[ -]?cousin|related[ -]?mating", tolower(question))
  is_nam_gxe_multiomics <- is_nam && grepl(
    "drought|irrigated|gxe|genotype.by.environment|rna.?seq|eqtl|gwas|genomic.?selection",
    tolower(question)
  )
  is_biparental <- grepl("biparental|bi[ -]?parental", tolower(question)) &&
    any(c("backcross", "selfing", "cross") %in% terms)
  if (!identical(preflight$status, "ready") ||
      !preflight$data_source %in% c("synthetic", "provided", "package_demo")) return("")

  if (identical(preflight$data_source, "package_demo") && is_human_irf6) {
    return(paste(
      "Bundled synthetic human IRF6-region demo:",
      "This runs a small synthetic VCF and GFF3 bundled with simitall. It does not contain real 1000 Genomes participants, does not estimate clinical risk, and must not be used to claim IRF6 variant causality.",
      "",
      "```r",
      "system2(\"Rscript\", c(",
      "  \"analysis/paper_fig/fig11_human_irf6_synthetic_liability.R\",",
      "  \"--data_dir\", \"inst/extdata/human_irf6\",",
      "  \"--out_dir\", \"results/demo_human_irf6\",",
      "  \"--seed\", \"81\"",
      "))",
      "```",
      "",
      "Expected outputs: synthetic marker truth, pedigree-liability tables, metadata JSON, and `figure11_human_irf6_synthetic_liability.png`.",
      sep = "\n"
    ))
  }

  if (identical(preflight$data_source, "package_demo") && is_nam) {
    return(paste(
      "Bundled synthetic maize NAM-style demo:",
      "This uses eight small synthetic chr10 founders and a small synthetic recombination map bundled with simitall. It is a runnable structural demo, not a named maize NAM panel or a publication-scale population.",
      "",
      "```r",
      "library(simitall)",
      "panel <- system.file(\"extdata\", \"panels\", \"demo_maize_nam_chr10.fa\", package = \"simitall\")",
      "map_file <- system.file(\"extdata\", \"maps\", \"demo_maize_nam_chr10_map.tsv\", package = \"simitall\")",
      "stopifnot(nzchar(panel), nzchar(map_file))",
      "dir.create(\"results/demo_maize_nam\", recursive = TRUE, showWarnings = FALSE)",
      "prefix <- \"results/demo_maize_nam/nam_demo\"",
      "simulate_breeding(",
      "  haplotype_fa = panel, out_prefix = prefix, scheme = \"NAM\",",
      "  founders = paste0(\"hap\", 1:8, collapse = \",\"), n_offspring = 80,",
      "  recomb_map_in = map_file, vcf_out = paste0(prefix, \".vcf\"),",
      "  ancestry_out = paste0(prefix, \".ancestry.tsv\"),",
      "  breakpoints_out = paste0(prefix, \".breakpoints.tsv\"), seed = 2026",
      ")",
      "# Validation figure: final alternate-allele frequency across the toy chr10.",
      "vcf_lines <- readLines(paste0(prefix, \".vcf\"), warn = FALSE)",
      "records <- strsplit(vcf_lines[!grepl(\"^#\", vcf_lines)], \"\\t\", fixed = TRUE)",
      "dosage <- function(gt) {",
      "  alleles <- strsplit(sub(\":.*$\", \"\", gt), \"[/|]\")[[1L]]",
      "  if (any(alleles == \".\")) return(NA_real_)",
      "  sum(as.numeric(alleles))",
      "}",
      "frequency <- data.frame(",
      "  position_bp = as.integer(vapply(records, `[[`, character(1), 2L)),",
      "  alt_frequency = vapply(records, function(x) {",
      "    values <- vapply(x[10:length(x)], dosage, numeric(1))",
      "    mean(values, na.rm = TRUE) / 2",
      "  }, numeric(1))",
      ")",
      "write.table(frequency, \"results/demo_maize_nam/nam_demo.allele_frequency.tsv\", sep = \"\\t\", row.names = FALSE, quote = FALSE)",
      "png(\"results/demo_maize_nam/nam_demo_validation.png\", width = 1500, height = 900, res = 160)",
      "plot(frequency$position_bp, frequency$alt_frequency, pch = 16, col = \"#156E8A\", ylim = c(0, 1),",
      "  xlab = \"Position on synthetic chr10 (bp)\", ylab = \"Alternate-allele frequency\",",
      "  main = \"Bundled NAM-style demo: final allele frequencies\")",
      "rect(1800, 0, 2600, 1, col = grDevices::adjustcolor(\"#E4A05D\", 0.18), border = NA)",
      "points(frequency$position_bp, frequency$alt_frequency, pch = 16, col = \"#156E8A\")",
      "legend(\"topright\", legend = \"Toy low-recombination interval\", fill = grDevices::adjustcolor(\"#E4A05D\", 0.18), bty = \"n\")",
      "dev.off()",
      "```",
      "",
      "Expected outputs: demo NAM VCF, sample/family metadata, ancestry tracts, breakpoint truth, realized recombination map, allele-frequency TSV, and `nam_demo_validation.png`.",
      sep = "\n"
    ))
  }

  if (is_human_irf6 && identical(preflight$data_source, "provided")) {
    script_path <- .simitall_agent_source_file(
      "analysis", "paper_fig", "fig11_human_irf6_synthetic_liability.R"
    )
    script_body <- if (!is.null(script_path)) {
      paste(readLines(script_path, warn = FALSE), collapse = "\n")
    } else {
      NULL
    }

    if (!is.null(script_body)) {
      return(paste(
        "Verified human IRF6-region execution recipe:",
        "The complete implementation below is read directly from the versioned simitall runner in this checkout. The VCF is used only to select anonymous common-marker frequencies. Liability effects and outcomes are synthetic; this is not a clinical risk predictor or a claim of variant causality.",
        "",
        "Run it directly from the repository root:",
        "",
        "```r",
        "system2(\"Rscript\", c(",
        "  \"analysis/paper_fig/fig11_human_irf6_synthetic_liability.R\",",
        "  \"--data_dir\", \"data/raw/human_irf6\",",
        "  \"--out_dir\", \"analysis/results/human_irf6_synthetic_liability\",",
        "  \"--seed\", \"81\"",
        "))",
        "```",
        "",
        "Full versioned R implementation:",
        "",
        "```r",
        script_body,
        "```",
        "",
        "Expected outputs: marker-truth TSV, individual synthetic-liability TSV, group summary TSV, metadata JSON, and `figure11_human_irf6_synthetic_liability.png`.",
        sep = "\n"
      ))
    }

    return("The IRF6 runner is unavailable from this installation. Run this request from a simitall source checkout containing `analysis/paper_fig/fig11_human_irf6_synthetic_liability.R` to retrieve the full verified implementation.")
  }

  if (!identical(preflight$data_source, "synthetic")) return("")

  if (is_human_inbreeding) {
    return(paste(
      "Synthetic human chromosome-10 inbreeding recipe:",
      "This is a fully synthetic marker simulation. It does not use real human variants, model a real family, or estimate clinical risk. Each round builds first-cousin-descendant pedigrees from the preceding synthetic population, so marker homozygosity is expected to increase on average.",
      "",
      "```r",
      "# Tune these values before running.",
      "out_dir <- \"results/human_chr10_inbreeding\"",
      "seed <- 2026",
      "chr_length_bp <- 10000000L   # Synthetic chr10 interval length.",
      "n_markers <- 500L            # Independent synthetic biallelic markers.",
      "n_lines <- 200L              # Individuals retained per round.",
      "n_rounds <- 3L               # Repeated first-cousin-pedigree rounds.",
      "maf_min <- 0.10; maf_max <- 0.40",
      "",
      "dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)",
      "set.seed(seed)",
      "marker_truth <- data.frame(",
      "  chromosome = \"chr10\",",
      "  position_bp = sort(sample.int(chr_length_bp, n_markers)),",
      "  alternate_allele_frequency = runif(n_markers, maf_min, maf_max)",
      ")",
      "draw_population <- function(n) {",
      "  matrix(rbinom(n * n_markers, 2, rep(marker_truth$alternate_allele_frequency, each = n)),",
      "         nrow = n, byrow = FALSE)",
      "}",
      "gamete <- function(genotype) rbinom(length(genotype), 1, genotype / 2)",
      "make_child <- function(parent_a, parent_b) gamete(parent_a) + gamete(parent_b)",
      "make_cousin_round <- function(pool, n) {",
      "  offspring <- matrix(NA_integer_, n, ncol(pool))",
      "  for (i in seq_len(n)) {",
      "    grandparents <- pool[sample.int(nrow(pool), 2, replace = TRUE), , drop = FALSE]",
      "    sibling_a <- make_child(grandparents[1, ], grandparents[2, ])",
      "    sibling_b <- make_child(grandparents[1, ], grandparents[2, ])",
      "    partners <- pool[sample.int(nrow(pool), 2, replace = TRUE), , drop = FALSE]",
      "    cousin_a <- make_child(sibling_a, partners[1, ])",
      "    cousin_b <- make_child(sibling_b, partners[2, ])",
      "    offspring[i, ] <- make_child(cousin_a, cousin_b)",
      "  }",
      "  offspring",
      "}",
      "summarize_round <- function(genotypes, round) {",
      "  data.frame(",
      "    sample = sprintf(\"round%02d_%03d\", round, seq_len(nrow(genotypes))),",
      "    round = round,",
      "    marker_homozygosity = rowMeans(genotypes == 0L | genotypes == 2L)",
      "  )",
      "}",
      "",
      "population <- draw_population(n_lines)",
      "individuals <- summarize_round(population, 0L)",
      "for (round in seq_len(n_rounds)) {",
      "  population <- make_cousin_round(population, n_lines)",
      "  individuals <- rbind(individuals, summarize_round(population, round))",
      "}",
      "round_summary <- aggregate(marker_homozygosity ~ round, individuals, mean)",
      "write.table(marker_truth, file.path(out_dir, \"marker_truth.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
      "write.table(individuals, file.path(out_dir, \"individual_homozygosity.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
      "write.table(round_summary, file.path(out_dir, \"round_summary.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
      "png(file.path(out_dir, \"chr10_inbreeding_homozygosity.png\"), width = 1500, height = 850, res = 160)",
      "par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))",
      "plot(round_summary$round, round_summary$marker_homozygosity, type = \"b\", pch = 16, col = \"#0D5F4C\",",
      "     xlab = \"Synthetic cousin-pedigree round\", ylab = \"Mean marker homozygosity\",",
      "     main = \"A. Homozygosity across rounds\")",
      "boxplot(marker_homozygosity ~ round, data = individuals, col = \"#DFF0E7\",",
      "        xlab = \"Round\", ylab = \"Marker homozygosity\",",
      "        main = \"B. Individual distributions\")",
      "dev.off()",
      "",
      "# Pass criterion for this programmed synthetic design:",
      "stopifnot(tail(round_summary$marker_homozygosity, 1) > round_summary$marker_homozygosity[1])",
      "round_summary",
      "```",
      "",
      "Tune `n_rounds`, `n_lines`, `n_markers`, allele-frequency bounds, chromosome length, and seed. The final `stopifnot()` checks only the expected direction of the synthetic mechanism: more related pedigree rounds should yield greater mean marker homozygosity than round 0.",
      sep = "\n"
    ))
  }

  if (is_nam_gxe_multiomics) {
    script_path <- .simitall_agent_source_file(
      "analysis", "paper_fig", "fig12_maize_nam_chr10_gxe_multiomics.R"
    )
    script_body <- if (!is.null(script_path)) {
      paste(readLines(script_path, warn = FALSE), collapse = "\n")
    } else {
      NULL
    }
    if (is.null(script_body)) {
      return("The NAM GxE multi-omics runner is unavailable from this installation. Run the request from a simitall source checkout containing `analysis/paper_fig/fig12_maize_nam_chr10_gxe_multiomics.R`.")
    }
    return(paste(
      "Verified NAM execution recipe: maize GxE and multi-omics workflow:",
      "The versioned runner simulates one synthetic chr10 NAM population, retains ancestry/breakpoint truth, creates irrigated and drought phenotypes with a focal GxE locus, derives RNA-seq/eQTL truth from the same lines, runs drought GWAS, and performs one GBLUP selection cycle. The bundled panel and focal gene are synthetic unless `--gff3` is supplied.",
      "",
      "Run it directly from the repository root:",
      "",
      "```r",
      "system2(\"Rscript\", c(",
      "  \"analysis/paper_fig/fig12_maize_nam_chr10_gxe_multiomics.R\",",
      "  \"--out_dir\", \"analysis/results/maize_nam_chr10_gxe_multiomics\",",
      "  \"--seed\", \"1201\"",
      "))",
      "```",
      "",
      "Expected figures: GxE phenotype distribution, regional GWAS truth overlay, GWAS Manhattan/QQ diagnostics, RNA-seq PCA, RNA-seq volcano/library-size diagnostics, and genomic-selection differential. Expected tables include breeding truth, trait/QTL truth, GWAS/eQTL results, RNA-seq counts, differential-expression results, and selected-parent/cross tables.",
      "",
      "Full versioned R implementation:",
      "",
      "```r",
      script_body,
      "```",
      sep = "\n"
    ))
  }

  if (is_biparental && !is_nam) {
    return(paste(
      "Verified biparental BC2S7 execution recipe:",
      "This is generated by simitall itself, not by the language model. It uses two synthetic founders labelled chr10, makes an F1, backcrosses twice to parent 1, self-fertilizes seven times, and writes an allele-frequency figure. It is not a simulation of named NAM varieties; use a supplied NAM founder panel and a measured chromosome-10 map for that.",
      "",
      "```r",
      "library(simitall)",
      "",
      "dir.create(\"results/biparental_bc2s7\", recursive = TRUE, showWarnings = FALSE)",
      "panel <- \"results/biparental_bc2s7/founders_chr10.fa\"",
      "generate_random_haplotype_panel(",
      "  out_fa = panel, n_haplotypes = 2, length = 2000, n_chromosomes = 1,",
      "  snp_rate = 0.02, indel_rate = 0, seed = 71",
      ")",
      "lines <- readLines(panel)",
      "lines[grepl(\"^>\", lines)] <- paste0(lines[grepl(\"^>\", lines)], \"|chr10\")",
      "writeLines(lines, panel)",
      "",
      "prefix <- \"results/biparental_bc2s7/bc2s7\"",
      "simulate_breeding(",
      "  haplotype_fa = panel, out_prefix = prefix, parents = \"hap1,hap2\",",
      "  sequence = \"F1,BC:P1:2,SELF:7\", n_offspring = 200,",
      "  vcf_out = paste0(prefix, \".vcf\"),",
      "  ancestry_out = paste0(prefix, \".ancestry.tsv\"),",
      "  breakpoints_out = paste0(prefix, \".breakpoints.tsv\"), seed = 72",
      ")",
      "",
      "# Calculate alternate-allele frequency at every final VCF site.",
      "vcf_lines <- readLines(paste0(prefix, \".vcf\"))",
      "records <- strsplit(vcf_lines[!grepl(\"^#\", vcf_lines)], \"\\t\", fixed = TRUE)",
      "dosage <- function(gt) {",
      "  a <- strsplit(sub(\":.*$\", \"\", gt), \"[/|]\")[[1L]]",
      "  if (any(a == \".\")) return(NA_real_)",
      "  sum(as.numeric(a))",
      "}",
      "freq <- data.frame(",
      "  position_bp = as.integer(vapply(records, `[[`, character(1), 2)),",
      "  alt_frequency = vapply(records, function(x) {",
      "    values <- vapply(x[10:length(x)], dosage, numeric(1))",
      "    mean(values, na.rm = TRUE) / 2",
      "  }, numeric(1))",
      ")",
      "write.table(freq, paste0(prefix, \".allele_frequencies.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
      "png(\"results/biparental_bc2s7/allele_frequency.png\", width = 1600, height = 1000, res = 180)",
      "plot(freq$position_bp, freq$alt_frequency, pch = 16, col = \"#156E8A\", ylim = c(0, 1),",
      "     xlab = \"Position on synthetic chr10 (bp)\", ylab = \"Alternate-allele frequency\",",
      "     main = \"Synthetic biparental BC2S7: final allele frequencies\")",
      "abline(h = 1 / 8, col = \"#B4422B\", lty = 2, lwd = 2)",
      "dev.off()",
      "```",
      "",
      "Expected outputs: `bc2s7.fa`, `bc2s7.meta.tsv`, `bc2s7.vcf`, ancestry and breakpoint truth TSVs, `bc2s7.allele_frequencies.tsv`, and `allele_frequency.png`. The 1/8 line is the expected mean donor-allele frequency after two backcrosses to P1 without selection; finite population and recombination make individual sites vary.",
      sep = "\n"
    ))
  }

  if (!is_nam) return("")

  paste(
    "Verified NAM execution recipe (template):",
    "This is generated by simitall itself, not by the language model. It creates a synthetic 100 Mb chr10 founder panel, makes 45-55 Mb a recombination cold spot (therefore a high-LD target region), simulates 500 NAM lines, and carries those same lines into phenotype, GWAS, RNA-seq, and genomic-selection steps.",
    "",
    "```r",
    "library(simitall)",
    "",
    "set.seed(2026)",
    "dir.create(\"results/maize_nam\", recursive = TRUE, showWarnings = FALSE)",
    "",
    "# 1. Build eight aligned founders for one 100 Mb chromosome.",
    "panel <- \"results/maize_nam/founders_chr10.fa\"",
    "generate_random_haplotype_panel(",
    "  out_fa = panel, n_haplotypes = 8, length = 100000000,",
    "  n_chromosomes = 1, snp_rate = 0.001, indel_rate = 0, seed = 101",
    ")",
    "# Label the single simulated chromosome as chr10 for the map and VCF.",
    "panel_lines <- readLines(panel)",
    "headers <- grepl(\"^>\", panel_lines)",
    "panel_lines[headers] <- paste0(panel_lines[headers], \"|chr10\")",
    "writeLines(panel_lines, panel)",
    "",
    "# 2. Make a chromosome-specific genetic map.",
    "# Lower cM/Mb from 45-55 Mb suppresses recombination and creates high LD.",
    "genome_bp <- 100000000L",
    "positions <- sort(unique(c(seq.int(1L, genome_bp, by = 1000000L), genome_bp)))",
    "midpoints <- head(positions, -1L) + diff(positions) / 2",
    "rate_cm_per_mb <- ifelse(midpoints >= 45000000 & midpoints <= 55000000, 0.05, 1.0)",
    "recomb_map <- data.frame(",
    "  chromosome = \"chr10\", pos_bp = positions,",
    "  cM = c(0, cumsum(diff(positions) * rate_cm_per_mb / 1e6))",
    ")",
    "map_file <- \"results/maize_nam/chr10_recombination_map.tsv\"",
    "write.table(recomb_map, map_file, sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "",
    "# 3. Simulate the NAM lines and retain ancestry, breakpoint, and VCF truth.",
    "nam_prefix <- \"results/maize_nam/nam500\"",
    "simulate_breeding(",
    "  haplotype_fa = panel, out_prefix = nam_prefix, scheme = \"NAM\",",
    "  founders = paste0(\"hap\", 1:8, collapse = \",\"), n_offspring = 500,",
    "  recomb_map_in = map_file, vcf_out = paste0(nam_prefix, \".vcf\"),",
    "  breakpoints_out = paste0(nam_prefix, \".breakpoints.tsv\"),",
    "  ancestry_out = paste0(nam_prefix, \".ancestry.tsv\"),",
    "  recomb_map_out = paste0(nam_prefix, \".realized_map.tsv\"), seed = 202",
    ")",
    "",
    "# 4. Simulate a quantitative trait directly from the NAM VCF.",
    "trait_prefix <- \"results/maize_nam/trait\"",
    "simulate_phenotypes(",
    "  geno_file = paste0(nam_prefix, \".vcf\"), out_prefix = trait_prefix,",
    "  heritability = 0.6, n_add_qtn = 20, n_dom_qtn = 5, n_epi_qtn = 2,",
    "  export_qtn = TRUE, seed = 303",
    ")",
    "# Add NAM family labels for family-aware GWAS and genomic-selection validation.",
    "pheno <- read.delim(paste0(trait_prefix, \".pheno.tsv\"), check.names = FALSE)",
    "names(pheno)[1] <- \"sample\"",
    "families <- read.delim(paste0(nam_prefix, \".meta.tsv\"), check.names = FALSE)",
    "pheno <- merge(pheno, families[, c(\"sample\", \"family\")], by = \"sample\", all.x = TRUE)",
    "pheno_file <- paste0(trait_prefix, \".with_family.tsv\")",
    "write.table(pheno, pheno_file, sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "trait <- setdiff(names(pheno)[vapply(pheno, is.numeric, logical(1))], character())[1]",
    "",
    "# 5. Analyze the actual NAM VCF. This is not simulate_gwas_cohort(), which creates a separate cohort.",
    "gwas <- analyze_gwas(",
    "  genotype_file = paste0(nam_prefix, \".vcf\"), phenotype = pheno_file,",
    "  out_prefix = \"results/maize_nam/gwas\", trait = trait,",
    "  fixed_effects = \"family\", n_pcs = 3, min_maf = 0.05",
    ")",
    "plot_gwas_results(gwas$paths$results, \"results/maize_nam/gwas_plot\")",
    "",
    "# 6. Simulate bulk RNA-seq from these same 500 NAM genotypes and test eQTL recovery.",
    "rna <- simulate_rnaseq_from_gwas(",
    "  genotype_file = paste0(nam_prefix, \".vcf\"),",
    "  out_prefix = \"results/maize_nam/rnaseq\", sample_metadata = pheno_file,",
    "  n_genes = 1000, n_cis_eqtl = 100, n_trans_eqtl = 25,",
    "  condition_levels = c(\"control\", \"stress\"),",
    "  batch_levels = c(\"batch1\", \"batch2\"), seed = 404",
    ")",
    "benchmark_eqtl(",
    "  genotype_file = paste0(nam_prefix, \".vcf\"),",
    "  expression_file = rna$expression, sample_metadata = rna$sample_metadata,",
    "  truth_file = rna$eqtl_truth, out_prefix = \"results/maize_nam/eqtl\"",
    ")",
    "",
    "# 7. Evaluate prediction with family-aware CV, then make a diverse parent/cross plan.",
    "cv <- cross_validate_genomic_prediction(",
    "  genotype_file = paste0(nam_prefix, \".vcf\"), phenotype = pheno_file,",
    "  trait = trait, model = \"gblup\", folds = 5, split = \"group\",",
    "  group = \"family\", seed = 505",
    ")",
    "selection <- run_genomic_selection(",
    "  genotype_file = paste0(nam_prefix, \".vcf\"), phenotype = pheno_file,",
    "  out_prefix = \"results/maize_nam/selection\", trait = trait, model = \"gblup\",",
    "  n_parents = 12, n_crosses = 30, mating = \"minimum_kinship\",",
    "  diversity_penalty = 0.25, family_data = pheno_file, family_col = \"family\",",
    "  max_per_family = 3, seed = 505",
    ")",
    "```",
    "",
    "Important: this breeding-first route currently analyzes the NAM VCF directly. `benchmark_gwas()` requires the standardized causal-truth schema written by `simulate_gwas_cohort()`, so do not call it on the simplePHENOTYPES QTN files until a converter is added.",
    sep = "\n"
  )
}

.simitall_agent_code_skeleton <- function(question, preflight) {
  terms <- .simitall_agent_terms(question)
  is_human_inbreeding <- .simitall_agent_species(question) == "human" &&
    grepl("inbreed|inbred|consanguin|first[ -]?cousin|related[ -]?mating", tolower(question))
  if (is_human_inbreeding) {
    generation_match <- regexec(
      "([0-9]+)\\s*(?:generations?|rounds?)", tolower(question), perl = TRUE
    )
    generation_value <- regmatches(tolower(question), generation_match)[[1L]]
    n_rounds <- if (length(generation_value) > 1L) generation_value[2L] else "10"
    return(paste(
      "Input-dependent human chr10 inbreeding code template:",
      "This template requires a real, build-matched chromosome-10 VCF/BCF. It reads observed allele frequencies from that input, then performs an allele-frequency-only finite-population pedigree approximation. It does not preserve phased haplotypes or LD; use a haplotype-aware pedigree engine for that extension.",
      "",
      "```r",
      "# Required: real chr10 VCF/BCF with INFO/AF and a bcftools installation.",
      "vcf <- \"PATH_TO_CHR10_VCF.gz\"",
      "out_dir <- \"results/human_chr10_inbreeding\"",
      "seed <- 2026",
      paste0("n_rounds <- ", n_rounds, "L"),
      "n_lines <- 200L",
      "maf_min <- 0.10; maf_max <- 0.40",
      "stopifnot(file.exists(vcf), nzchar(Sys.which(\"bcftools\")))",
      "dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)",
      "set.seed(seed)",
      "query_file <- tempfile(fileext = \".tsv\")",
      "command <- paste(\"bcftools query -f\", shQuote(\"%CHROM\\t%POS\\t%INFO/AF\\n\"), shQuote(vcf), \">\", shQuote(query_file))",
      "stopifnot(system(command) == 0L)",
      "markers <- read.delim(query_file, header = FALSE, sep = \"\\t\", stringsAsFactors = FALSE)",
      "names(markers) <- c(\"chromosome\", \"position_bp\", \"af\")",
      "markers$af <- suppressWarnings(as.numeric(markers$af))",
      "markers <- markers[markers$chromosome %in% c(\"10\", \"chr10\") & is.finite(markers$af) & markers$af >= maf_min & markers$af <= maf_max, ]",
      "stopifnot(nrow(markers) >= 20L)",
      "markers <- markers[sort(sample.int(nrow(markers), min(500L, nrow(markers)))), ]",
      "draw_population <- function(n) matrix(rbinom(n * nrow(markers), 2, rep(markers$af, each = n)), nrow = n, byrow = FALSE)",
      "gamete <- function(g) rbinom(length(g), 1, g / 2)",
      "child <- function(a, b) gamete(a) + gamete(b)",
      "cousin_round <- function(pool, n) {",
      "  out <- matrix(NA_integer_, n, ncol(pool))",
      "  for (i in seq_len(n)) {",
      "    gp <- pool[sample.int(nrow(pool), 2L, replace = TRUE), , drop = FALSE]",
      "    sib1 <- child(gp[1, ], gp[2, ]); sib2 <- child(gp[1, ], gp[2, ])",
      "    partners <- pool[sample.int(nrow(pool), 2L, replace = TRUE), , drop = FALSE]",
      "    out[i, ] <- child(child(sib1, partners[1, ]), child(sib2, partners[2, ]))",
      "  }; out",
      "}",
      "population <- draw_population(n_lines)",
      "summary <- data.frame(round = 0L, mean_af = mean(population) / 2, mean_homozygosity = mean(population == 0L | population == 2L))",
      "for (round in seq_len(n_rounds)) {",
      "  population <- cousin_round(population, n_lines)",
      "  summary <- rbind(summary, data.frame(round = round, mean_af = mean(population) / 2, mean_homozygosity = mean(population == 0L | population == 2L)))",
      "}",
      "write.table(markers, file.path(out_dir, \"input_marker_frequencies.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
      "write.table(summary, file.path(out_dir, \"round_summary.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
      "png(file.path(out_dir, \"allele_frequency_and_homozygosity.png\"), width = 1500, height = 750, res = 160)",
      "par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))",
      "plot(summary$round, summary$mean_af, type = \"b\", pch = 16, col = \"#156E8A\", xlab = \"Pedigree round\", ylab = \"Mean alternate-allele frequency\", main = \"A. Allele frequency: no directional expectation\")",
      "plot(summary$round, summary$mean_homozygosity, type = \"b\", pch = 16, col = \"#B4422B\", xlab = \"Pedigree round\", ylab = \"Mean marker homozygosity\", main = \"B. Expected increase in homozygosity\")",
      "dev.off()",
      "summary",
      "```",
      "",
      "Interpretation: without selection, mutation, or migration, allele frequencies have no directional expectation under inbreeding; observed changes in this finite simulation are genetic drift. The robust expected signal is increased homozygosity and reduced heterozygosity.",
      sep = "\n"
    ))
  }
  if (!"nam" %in% terms) return("")
  paste(
    "Input-dependent NAM code skeleton:",
    "Replace the two `PATH_TO_*` values with a real aligned founder panel and a compatible chr10 recombination map. This skeleton is shown even when inputs are unresolved; it becomes runnable only after those paths exist.",
    "",
    "```r",
    "library(simitall)",
    "",
    "founder_panel <- \"PATH_TO_FOUNDER_PANEL.fa\"  # >founder|chr10 records",
    "recomb_map <- \"PATH_TO_CHR10_RECOMBINATION_MAP.tsv\"  # chromosome, pos_bp, cM",
    "stopifnot(file.exists(founder_panel), file.exists(recomb_map))",
    "dir.create(\"results/study\", recursive = TRUE, showWarnings = FALSE)",
    "",
    "nam_prefix <- \"results/study/nam500\"",
    "simulate_breeding(",
    "  haplotype_fa = founder_panel, out_prefix = nam_prefix, scheme = \"NAM\",",
    "  n_offspring = 500, recomb_map_in = recomb_map,",
    "  vcf_out = paste0(nam_prefix, \".vcf\"),",
    "  ancestry_out = paste0(nam_prefix, \".ancestry.tsv\"),",
    "  breakpoints_out = paste0(nam_prefix, \".breakpoints.tsv\"), seed = 202",
    ")",
    "",
    "simulate_phenotypes(",
    "  geno_file = paste0(nam_prefix, \".vcf\"), out_prefix = \"results/study/trait\",",
    "  heritability = 0.6, n_add_qtn = 20, n_dom_qtn = 5, seed = 303",
    ")",
    "pheno <- read.delim(\"results/study/trait.pheno.tsv\", check.names = FALSE)",
    "names(pheno)[1] <- \"sample\"",
    "family <- read.delim(paste0(nam_prefix, \".meta.tsv\"), check.names = FALSE)",
    "pheno <- merge(pheno, family[, c(\"sample\", \"family\")], by = \"sample\", all.x = TRUE)",
    "write.table(pheno, \"results/study/trait.with_family.tsv\", sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "trait <- names(pheno)[vapply(pheno, is.numeric, logical(1))][1]",
    "",
    "gwas <- analyze_gwas(",
    "  genotype_file = paste0(nam_prefix, \".vcf\"), phenotype = \"results/study/trait.with_family.tsv\",",
    "  out_prefix = \"results/study/gwas\", trait = trait, fixed_effects = \"family\", n_pcs = 3",
    ")",
    "rna <- simulate_rnaseq_from_gwas(",
    "  genotype_file = paste0(nam_prefix, \".vcf\"), out_prefix = \"results/study/rnaseq\",",
    "  sample_metadata = \"results/study/trait.with_family.tsv\", n_genes = 1000,",
    "  n_cis_eqtl = 100, n_trans_eqtl = 25, seed = 404",
    ")",
    "selection <- run_genomic_selection(",
    "  genotype_file = paste0(nam_prefix, \".vcf\"), phenotype = \"results/study/trait.with_family.tsv\",",
    "  out_prefix = \"results/study/selection\", trait = trait, model = \"gblup\",",
    "  family_data = \"results/study/trait.with_family.tsv\", family_col = \"family\", seed = 505",
    ")",
    "```",
    sep = "\n"
  )
}

.simitall_agent_species <- function(question) {
  text <- tolower(question)
  if (grepl("arabidopsis|a\\. thaliana|thaliana", text)) return("arabidopsis")
  if (grepl("maize|zea mays", text)) return("maize")
  if (grepl("rice|oryza", text)) return("rice")
  if (grepl("human|grch|homo sapiens", text)) return("human")
  if (grepl("e\\. coli|ecoli|escherichia", text)) return("bacteria")
  "unspecified"
}

.simitall_agent_preflight <- function(question, input_files, data_source) {
  terms <- .simitall_agent_terms(question)
  species <- .simitall_agent_species(question)
  needs_panel <- any(c("breeding", "nam", "magic", "ril", "nil", "f1", "f2", "backcross", "selfing", "cross") %in% terms) ||
    grepl("inbreed|inbred|consanguin|first[ -]?cousin|related[ -]?mating", tolower(question))
  input_files <- input_files %||% character()
  input_files <- as.character(input_files)
  supplied <- input_files[nzchar(input_files)]
  missing_files <- supplied[!file.exists(path.expand(supplied))]
  supplied <- supplied[file.exists(path.expand(supplied))]

  status <- "ready"
  action <- ""
  if (length(missing_files)) {
    status <- "blocked"
    action <- paste0("These declared input files do not exist: ", paste(missing_files, collapse = ", "), ".")
  } else if (needs_panel && identical(data_source, "provided") && !length(supplied)) {
    status <- "needs_input"
    action <- "Provide a founder haplotype FASTA panel and a chromosome-specific recombination map."
  } else if (needs_panel && identical(species, "arabidopsis") && grepl("1001", tolower(question)) && !length(supplied)) {
    status <- "needs_input"
    action <- paste(
      "The package does not bundle an Arabidopsis 1001 Genomes haplotype panel.",
      "Provide the panel FASTA/VCF and a compatible genetic map; the agent will not substitute a generic reference."
    )
  } else if (needs_panel && identical(data_source, "auto") && !length(supplied)) {
    status <- "needs_input"
    action <- paste(
      "Choose a data source: provide a founder panel and map, or explicitly set",
      "data_source = 'synthetic' to generate a synthetic panel. The bundled demo panel is only a toy example."
    )
  }

  list(
    status = status,
    species = species,
    input_files = supplied,
    action = action,
    needs_panel = needs_panel,
    data_source = data_source
  )
}

.simitall_local_answer <- function(question, retrieval) {
  terms <- .simitall_agent_terms(question)
  has_any <- function(words) any(words %in% terms)
  is_human_inbreeding <- .simitall_agent_species(question) == "human" &&
    grepl("inbreed|inbred|consanguin|first[ -]?cousin|related[ -]?mating", tolower(question))
  is_pedigree <- has_any(c("pedigree", "parentage", "maternal", "paternal")) || is_human_inbreeding
  is_breeding <- has_any(c("breeding", "nam", "magic", "ril", "nil", "f1", "f2", "backcross", "selfing", "cross")) || is_pedigree
  is_ancestry_tools <- has_any(c("rfmix", "tractor", "genesis", "genomeadmixr"))
  # "local" commonly describes a file on disk (for example, a local GRCh38
  # VCF). Only treat it as an ancestry request when the phrase is explicit.
  is_ancestry <- has_any(c("ancestry", "admixed", "admixture")) ||
    grepl("local[[:space:]-]+ancestry", tolower(question)) || is_ancestry_tools
  is_bacterial <- has_any(c("bacteria", "bacterial", "ecoli", "coli", "escherichia"))
  is_assembly <- has_any(c("assembly", "hybrid", "unicycler", "quast", "illumina", "pacbio", "hifi", "clr"))
  is_polyploid <- has_any(c("polyploid", "tetraploid", "hexaploid", "potato", "wheat"))
  source <- function(pattern) {
    hits <- grep(pattern, paste(retrieval$category, retrieval$title), ignore.case = TRUE)
    if (!length(hits)) return("")
    paste0(" [SOURCE: ", retrieval$path[hits[1L]], "]")
  }
  stages <- character()
  required <- character()
  truth <- character()
  limits <- "This is a simulation plan, not evidence that a synthetic population reproduces real biology. Use measured maps, founder panels, and population metadata when making organism-level claims."

  if (is_human_inbreeding) {
    stages <- c(
      stages,
      paste0(
        "1. Separate allele frequency from genotype frequency: inbreeding alone does not systematically change the expected allele frequency at a neutral locus, but it reduces heterozygosity and increases homozygosity. Finite population size, selection, mutation, or migration can change allele frequencies.",
        source("Human|Pedigree")
      ),
      paste0(
        "2. Define the mating relationship for each generation (for example, repeated first-cousin unions) and simulate or analyze chromosome-10 genotypes while retaining allele-frequency and heterozygosity summaries separately.",
        source("Pedigree")
      )
    )
  }

  if (is_ancestry) {
    stages <- c(
      stages,
      paste0(
        "1. Decide whether the question requires global-structure adjustment, admixture mapping, or local-ancestry-aware SNP effects; these are different analyses.",
        source("Ancestry-Aware|Human")
      ),
      paste0(
        "2. Lock the genome build, genotype QC, cohort/sample manifest, phenotype, covariates, relatedness treatment, and reference-panel definition.",
        source("Ancestry-Aware|Human")
      ),
      paste0(
        "3. For local ancestry, require phased cohort haplotypes, suitable reference panels, and an external versioned inference method before association testing; retain tract calls and uncertainty.",
        source("Ancestry-Aware")
      )
    )
    required <- c(
      required,
      "Cohort genotypes with a fixed genome build, phenotype/covariate metadata, and relatedness handling.",
      "For local ancestry: phased haplotypes, reference-panel sample manifests, ancestry labels, and local-ancestry calls with uncertainty."
    )
    truth <- c(
      truth,
      "Global-ancestry/PC summaries, local-ancestry tract truth or calls, reference-panel provenance, and calibration metrics stratified by ancestry/group."
    )
    limits <- paste(
      limits,
      "simitall currently plans ancestry-aware workflows and can simulate structured cohorts, but it does not itself phase data, infer local ancestry, or execute RFMix/Tractor."
    )
    if (is_ancestry_tools) {
      stages <- c(stages, paste0(
        length(stages) + 1L, ". Use `GenomeAdmixR` to simulate ancestry-tract truth; use `GENESIS` for R-side structure/relatedness or admixture mapping after local ancestry is available; use external RFMix for local-ancestry inference and Tractor for ancestry-specific GWAS.",
        source("GenomeAdmixR|GENESIS|RFMix|Tractor")
      ))
    }
  }

  if (is_breeding) {
    stages <- c(
      stages,
      paste0(
        "1. Define founder haplotypes and chromosome-specific recombination maps; ",
        "use `generate_random_haplotype_panel()` only when an existing founder panel is unavailable.",
        source("Maize|Breeding")
      ),
      paste0(
        "2. Create the requested designed-cross population with `simulate_breeding()`; ",
        "retain family labels, ancestry tracts, breakpoints, and map truth.",
        source("Breeding")
      )
    )
    required <- c(required, "Founder FASTA haplotypes or a synthetic founder-panel design, genetic maps, mating design, and population size.")
    truth <- c(truth, "Per-chromosome ancestry, breakpoints, founder contributions, heterozygosity, and LD/segregation QC.")
    if (is_polyploid) {
      stages <- c(stages, paste0(
        length(stages) + 1L, ". Before simulating crosses, define ploidy, dosage encoding, chromosome/subgenome naming, and pairing assumptions. Do not apply diploid genotype-frequency expectations to a tetraploid or hexaploid design.",
        source("Potato|Wheat")
      ))
      required <- c(required, "For polyploids: ploidy, allele-dosage representation, pairing model, and whether the scenario is a simplified diploidized analysis.")
    }
  }
  if (is_pedigree) {
    stages <- c(stages, paste0(
      length(stages) + 1L, ". Validate unique individual IDs, parent IDs, unknown parents, cycles, aliases, and generations before using the pedigree. Use `simulate_parental_effects()` for direct/maternal/paternal phenotype effects, or `simupop_api()` for advanced PedigreeMating.",
      source("Pedigree")
    ))
    required <- c(required, "A pedigree table with unique id, sire, and dam columns; optional family, generation, sex, and source metadata.")
    truth <- c(truth, "Validated pedigree, cross graph, family/generation summaries, and genotype-pedigree concordance checks when markers are available.")
  }
  if (is_human_inbreeding) {
    stages <- c(stages, paste0(
      length(stages) + 1L,
      ". Choose the data mode explicitly: use a synthetic marker panel for a mechanism check, or provide build-matched real chromosome-10 genotypes and a pedigree for a data-grounded simulation.",
      source("Human|Pedigree")
    ))
    required <- c(
      required,
      "Synthetic mode: no external chromosome file is required; the recipe generates independent synthetic chr10 markers and allele frequencies.",
      "Real-data mode: a genome-build-matched chr10 reference FASTA is useful for sequence-aware simulation, but a phased VCF/BCF plus index is the essential genotype input; also provide a compatible recombination map with chromosome, position, and cM columns.",
      "For a pedigree-grounded design: a de-identified pedigree table with id, sire, dam, and generation fields, plus documented consent/governance appropriate to the source data."
    )
    truth <- c(
      truth,
      "Data-mode label, genome build, input-file checksums, marker/allele-frequency truth for synthetic mode, and realized homozygosity by pedigree round."
    )
    limits <- paste(
      limits,
      "The built-in generic chr10 inbreeding recipe is synthetic and labels chr10 only as a coordinate namespace; it does not read real GRCh38 chromosome-10 sequence or human genotypes. Use supplied, build-matched VCF/BCF and map files for a data-grounded extension."
    )
  }
  if (has_any(c("gwas", "qtl", "association"))) {
    if (is_ancestry) {
      stages <- c(stages, paste0(
        length(stages) + 1L, ". Run the pre-specified ancestry-aware association method only after the local/global ancestry inputs are resolved; report the method boundary and sensitivity to reference panels and covariates.",
        source("Ancestry-Aware")
      ))
    } else if (is_bacterial) {
      stages <- c(stages, paste0(
        length(stages) + 1L, ". Use a bacterial-GWAS design with an explicit lineage/recombination model and variant representation such as SNPs, k-mers, genes, or accessory elements. Do not treat a clonal haploid collection as a diploid human GWAS cohort.",
        source("Escherichia|Bacterial")
      ))
      required <- c(required, "Bacterial isolate panel, lineage/phylogeny or kinship representation, accessory-genome/variant encoding, phenotype, and recombination assumptions.")
      truth <- c(truth, "Lineage-aware calibration summaries, causal truth where simulated, and variant-representation provenance.")
    } else if (is_breeding) {
      stages <- c(stages, paste0(
        length(stages) + 1L, ". Simulate a phenotype from the breeding VCF with `simulate_phenotypes()`, then run `analyze_gwas()` on that same VCF with family-aware covariates.",
        source("Genome-Wide Association|Simulate and Run a GWAS")
      ))
      required <- c(required, "Trait architecture, phenotype design, and family covariates carried forward from the breeding metadata.")
      truth <- c(truth, "Breeding VCF, phenotype table, family metadata, and GWAS results. A standardized GWAS causal-truth benchmark is not yet emitted by this breeding-first route.")
    } else {
      stages <- c(stages, paste0(
        length(stages) + 1L, ". Simulate the genotype/phenotype cohort with `simulate_gwas_cohort()`, then use `analyze_gwas()` and `benchmark_gwas()` against causal truth.",
        source("Genome-Wide Association|Simulate and Run a GWAS")
      ))
      required <- c(required, "Trait architecture, causal-locus assumptions, population/family covariates, and phenotype design.")
      truth <- c(truth, "Genotypes, phenotypes, causal-variant truth, GWAS results, and recovery metrics.")
    }
  }
  if (is_assembly) {
    stages <- c(stages, paste0(
      length(stages) + 1L, ". For a bacterial hybrid-assembly benchmark, simulate/read in short and long reads, run `run_grid_both()` and `run_unicycler_grid()`, then evaluate against the known reference with `run_quast_grid()`.",
      source("Hybrid Assembly|Unicycler|QUAST|Escherichia")
    ))
    required <- c(required, "Reference genome, short-read and long-read technology settings, coverage grid, assembler mode, and known truth reference for evaluation.")
    truth <- c(truth, "Read-generation parameters, assembly FASTA/graphs/logs, and reference-aware QUAST metrics.")
  }
  if (has_any(c("rna", "rnaseq", "eqtl", "expression", "transcript"))) {
    stages <- c(stages, paste0(
      length(stages) + 1L, ". Derive expression from the same simulated donor genotypes using `simulate_rnaseq_from_gwas()` and compare eQTL results with `benchmark_eqtl()`.",
      source("Bulk RNA|GWAS-Linked")
    ))
    required <- c(required, "Conditions, biological replicates, balanced batches, and optional annotation.")
    truth <- c(truth, "Expression/eQTL truth, counts, sample metadata, and eQTL recovery metrics.")
  }
  if (has_any(c("genomic", "selection", "gblup", "rrblup"))) {
    stages <- c(stages, paste0(
      length(stages) + 1L, ". Fit and validate prediction with `fit_genomic_model()`, then run `run_genomic_selection()` with family-aware validation and diversity checks.",
      source("Genomic Selection")
    ))
    required <- c(required, "A phenotyped training set, genotyped candidates, and family labels for leakage-aware validation.")
    truth <- c(truth, "Prediction accuracy, calibration, selected-parent GEBVs, kinship, family representation, and selection differential.")
  }
  if (!length(stages)) {
    stages <- "1. Start by defining the biological unit, reference/founder inputs, and desired truth outputs before selecting a simitall workflow."
    required <- "The target organism, study design, sample count, and desired analysis endpoint."
    truth <- "Parameters, random seed, source inputs, and simulation truth files."
  }
  sources <- vapply(seq_len(nrow(retrieval)), function(i) {
    paste0("- ", retrieval$title[i], " [SOURCE: ", retrieval$path[i], "]")
  }, character(1))
  paste(
    "Offline simitall plan (no cloud model or API was used).",
    paste0("Question: ", question),
    "Interpreted design:",
    paste0("- ", question),
    "Proposed stages:", paste(stages, collapse = "\n"),
    "Required inputs:", paste0("- ", unique(required), collapse = "\n"),
    "Truth and QC outputs:", paste0("- ", unique(truth), collapse = "\n"),
    "Limitations:",
    paste0("- ", limits),
    "Grounding documents:", paste(sources, collapse = "\n"),
    sep = "\n\n"
  )
}

#' Ask the grounded simitall agent a question
#'
#' Retrieve relevant local simitall knowledge documents and answer with a
#' selected backend. The default local backend is free and offline: it returns
#' a transparent evidence briefing and candidate package functions. The
#' optional Ollama backend sends context to a model running on the user's own
#' computer, while OpenAI uses the paid API. The agent never executes R,
#' commands, downloads, or simulations.
#'
#' Install Ollama and download a model such as llama3.2 to use the local-LLM
#' backend. Set OPENAI_API_KEY in .Renviron only when explicitly choosing the
#' paid OpenAI backend. Use dry_run = TRUE to inspect a prompt without a model
#' request.
#'
#' @param question One scientific or workflow question.
#' @param provider Answer backend: local (default, free, no LLM), ollama
#'   (free local LLM), or openai (paid API).
#' @param model Model identifier. Defaults depend on provider: none for local
#'   retrieval, SIMITALL_OLLAMA_MODEL or llama3.2 for Ollama, and
#'   SIMITALL_OPENAI_MODEL or gpt-5-mini for OpenAI.
#' @param api_key OpenAI API key, used only for provider = "openai".
#' @param n_context Number of local documents to retrieve.
#' @param knowledge_dir Optional alternate knowledge directory.
#' @param ollama_endpoint URL of a local Ollama API server.
#' @param allow_code Logical; whether to permit R code in a model response.
#'   Defaults to `FALSE` so a local model cannot silently invent packages,
#'   objects, or parameter names. When `TRUE`, code remains advisory and is
#'   restricted to the explicit `simitall` function allowlist.
#' @param include_code Logical; append a deterministic, package-verified code
#'   recipe when `simitall` recognizes a supported workflow request. This is
#'   safer than asking a language model to invent executable code. Currently,
#'   recipes are supplied for synthetic NAM breeding-to-omics requests and
#'   synthetic biparental BC2S7 allele-frequency requests.
#' @param input_files Optional named or unnamed paths to user-supplied data.
#'   For example, supply a founder FASTA panel and recombination map when
#'   requesting a real breeding population. Files are checked locally and are
#'   never uploaded by the agent.
#' @param data_source Input-source policy: `"auto"` asks for real data instead
#'   of silently substituting a toy resource, `"provided"` requires
#'   `input_files`, `"synthetic"` explicitly permits a generated input, and
#'   `"package_demo"` permits only the bundled demonstration files.
#' @param dry_run Return the constructed request without calling an API.
#'
#' @return A `simitall_agent_answer` object with `answer`, ranked `sources`,
#'   `model`, and `question`. With `dry_run = TRUE`, `answer` is `NULL` and the
#'   object also contains `instructions` and `context`.
#' @examples
#' \dontrun{
#' answer <- simitall_ask(
#'   "Plan a multi-chromosome maize NAM population for GWAS and genomic selection.",
#'   provider = "local"
#' )
#' print(answer)
#' }
#' simitall_ask("What truth files should a GWAS simulation retain?", dry_run = TRUE)
#' @export
simitall_ask <- function(
    question,
    provider = c("local", "ollama", "openai"),
    model = NULL,
    api_key = Sys.getenv("OPENAI_API_KEY"),
    n_context = 6L,
    knowledge_dir = NULL,
    ollama_endpoint = Sys.getenv("SIMITALL_OLLAMA_ENDPOINT", "http://localhost:11434/api/generate"),
    allow_code = FALSE,
    include_code = FALSE,
    input_files = NULL,
    data_source = c("auto", "provided", "synthetic", "package_demo"),
    dry_run = FALSE) {
  provider <- match.arg(provider)
  data_source <- match.arg(data_source)
  if (is.null(model)) {
    model <- switch(
      provider,
      local = "none",
      ollama = Sys.getenv("SIMITALL_OLLAMA_MODEL", "llama3.2"),
      openai = Sys.getenv("SIMITALL_OPENAI_MODEL", "gpt-5-mini")
    )
  }
  retrieval <- search_simitall_knowledge(question, n_context, knowledge_dir)
  preflight <- .simitall_agent_preflight(question, input_files, data_source)
  context <- .simitall_agent_context(
    retrieval,
    max_chars = if (identical(provider, "ollama")) 8000L else 18000L
  )
  code_policy <- if (isTRUE(allow_code)) {
    paste(
      "You may include one small advisory R code block only when every simitall",
      "function and every argument is explicitly supported by the supplied context.",
      "Do not use install.packages(), library(), external packages, example data",
      "objects, shell commands, or undocumented parameter names."
    )
  } else {
    paste(
      "Do not write R code, shell commands, package-installation instructions,",
      "or code blocks. Give a conceptual, function-level plan only."
    )
  }
  instructions <- paste(
    "You are simitall, an evidence-aware genomics simulation and analysis assistant.",
    "Answer only from the supplied local knowledge context and the user's question.",
    "Treat package workflows, techniques, species notes, and tool cards as higher-priority",
    "operational evidence than review papers. Clearly distinguish simulation truth from",
    "real-world biological validation. First restate every requested design constraint",
    "under the heading 'Interpreted design'. Then use these headings: 'Proposed stages',",
    "'Required inputs', 'Truth and QC outputs', and 'Limitations'. If a requested feature",
    "is absent from the local context or allowlist, say that it is unsupported or needs",
    "implementation; never substitute an external tool. Do not invent functions, paper",
    "findings, package behavior, versions, results, files, objects, or argument names.",
    "You cannot execute R, shell commands, downloads, or simulations.",
    "For a breeding or mapping-population question that also requests downstream",
    "GWAS, RNA-seq, or genomic selection, the stage order is: obtain or create founder",
    "haplotypes and chromosome-specific recombination maps; simulate the breeding",
    "population; create the genotype/phenotype GWAS cohort; derive linked RNA-seq from",
    "those same individuals; then fit and evaluate genomic selection. Do not propose",
    "genome generation or annotation unless the user explicitly asks for it, and never",
    "treat annotation as something to add separately to every individual.",
    code_policy,
    "Only mention these simitall functions, and only when relevant:\n",
    .simitall_agent_allowed_tools_text(),
    "\nSynthesize one direct answer; do not reproduce document headings, quote source",
    "text, or answer each evidence card separately. Keep the answer concise (at most",
    "400 words), use short numbered steps, and cite",
    "every scientific or workflow claim using the relevant [SOURCE: absolute/path] from",
    "the local context.\n\nUSER QUESTION:\n", question,
    "\n\nLOCAL KNOWLEDGE CONTEXT:\n", context
  )
  result <- list(
    question = question,
    model = model,
    sources = retrieval,
    answer = NULL,
    instructions = instructions,
    context = context,
    preflight = preflight,
    execution = "disabled; planning and question answering only"
  )
  if (!isTRUE(dry_run)) {
    if (!identical(preflight$status, "ready")) {
      result$answer <- paste(
        .simitall_local_answer(question, retrieval),
        "Data-source resolution required:",
        paste0("- Status: ", preflight$status),
        paste0("- ", preflight$action),
        sep = "\n\n"
      )
    } else if (identical(provider, "local")) {
      result$answer <- .simitall_local_answer(question, retrieval)
    } else if (identical(provider, "ollama")) {
      result$answer <- .simitall_ollama_request(instructions, model, ollama_endpoint)
    } else {
      if (!nzchar(api_key)) {
        stop(
          "OPENAI_API_KEY is not set. Add it to ~/.Renviron or use dry_run = TRUE ",
            "to inspect the grounded prompt without contacting OpenAI.",
          call. = FALSE
        )
      }
      result$answer <- .simitall_openai_request(question, instructions, model, api_key)
    }
    if (!identical(provider, "local")) {
      result$answer <- .simitall_agent_validate_answer(
        result$answer, question, retrieval, allow_code
      )
    }
    if (isTRUE(include_code)) {
      recipe <- .simitall_agent_code_recipe(question, preflight)
      if (nzchar(recipe)) {
        result$answer <- paste(result$answer, recipe, sep = "\n\n")
      } else {
        skeleton <- .simitall_agent_code_skeleton(question, preflight)
        if (nzchar(skeleton)) result$answer <- paste(result$answer, skeleton, sep = "\n\n")
      }
      if (!identical(preflight$status, "ready")) {
        result$answer <- paste(
          result$answer,
          "The code skeleton is not runnable until the required input has been resolved.",
          sep = "\n\n"
        )
      }
    }
    result$instructions <- NULL
    result$context <- NULL
  }
  class(result) <- "simitall_agent_answer"
  result
}

#' @export
print.simitall_agent_answer <- function(x, ...) {
  cat("simitall agent", if (is.null(x$answer)) "request preview" else "answer", "\n")
  cat("Question: ", x$question, "\n", sep = "")
  if (!is.null(x$answer)) cat("\n", x$answer, "\n", sep = "")
  cat("\nGrounding documents:\n")
  print(x$sources[, c("category", "title", "path"), drop = FALSE], row.names = FALSE)
  invisible(x)
}

#' Launch the local simitall agent Shiny application
#'
#' Start an optional local chat interface for [simitall_ask()]. The app reads
#' `OPENAI_API_KEY` only from the local R session environment; it does not
#' connect to GitHub, upload the package repository, or run simulations. Heavy
#' simulations should remain queued through a future execution service rather
#' than being launched in an interactive Shiny worker.
#'
#' @param ... Arguments passed to [shiny::runApp()], such as `launch.browser`.
#'
#' @return Starts a Shiny application and returns invisibly when it stops.
#' @examples
#' \dontrun{
#' launch_simitall_agent(launch.browser = TRUE)
#' }
#' @export
launch_simitall_agent <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Install it with install.packages('shiny').")
  }
  app_dir <- system.file("shiny", "simitall-agent", package = "simitall")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    app_dir <- file.path("inst", "shiny", "simitall-agent")
  }
  if (!dir.exists(app_dir)) stop("The simitall Shiny app could not be located")
  shiny::runApp(app_dir, ...)
}
