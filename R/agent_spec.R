# Typed workflow specs for the simitall agent.
#
# The pipeline is: question -> spec -> validate -> resolve inputs -> compile.
# Every stage is a plain list in, plain list out, so any extractor (this
# keyword baseline today, an LLM emitting the same JSON tomorrow) can be
# swapped in and scored against the same gold specs.
#
# Supported in v1: quantitative-trait GWAS on a NAM population or an inbred
# diversity panel. Anything else returns supported = FALSE and the agent falls
# back to its existing recipe router.

# ---------------------------------------------------------------------------
# Panel catalog: what simitall knows about, and whether it ships the data.
# ---------------------------------------------------------------------------

.simitall_panel_catalog <- function() {
  list(
    maize_nam = list(
      panel_id = "maize_nam",
      organism = "maize",
      design = "NAM",
      label = "maize NAM-style population",
      bundled = TRUE,
      bundled_note = paste(
        "Bundled panel: eight small synthetic chr10 founders and a toy",
        "recombination map shipped with simitall. It is not the real maize",
        "NAM founder panel."
      ),
      missing_note = ""
    ),
    arabidopsis_1001 = list(
      panel_id = "arabidopsis_1001",
      organism = "arabidopsis",
      design = "diversity_panel",
      label = "Arabidopsis 1001 Genomes-style accession panel",
      bundled = FALSE,
      bundled_note = "",
      missing_note = paste(
        "The package does not bundle an Arabidopsis 1001 Genomes haplotype panel."
      ),
      # Synthetic stand-in assumptions. These are modelling choices, not
      # estimates from the real 1001 Genomes data.
      synthetic_defaults = list(
        n_subpops = 4L, fst = 0.10, n_chromosomes = 5L,
        chromosome_length_bp = 1000000L, n_markers = 2000L, inbred = TRUE
      )
    )
  )
}

# ---------------------------------------------------------------------------
# Stage 1: extraction (keyword baseline).
# ---------------------------------------------------------------------------

.simitall_agent_first_number <- function(text, pattern) {
  hit <- regmatches(text, regexec(pattern, text, perl = TRUE))[[1L]]
  if (length(hit) < 2L) return(NULL)
  value <- suppressWarnings(as.numeric(gsub(",", "", hit[2L], fixed = TRUE)))
  if (is.na(value)) NULL else value
}

`%|na|%` <- function(a, b) if (is.null(a) || length(a) != 1L || is.na(a)) b else a

# Synonyms from inst/agent/aliases.tsv ("corn" -> "maize", "water deficit" ->
# "drought"). Uses the same phrase matching as retrieval, and returns the
# matched alias -> canonical pairs so the plan can show what was normalised.
.simitall_agent_alias_matches <- function(question) {
  aliases <- tryCatch(.simitall_agent_aliases(), error = function(e) NULL)
  if (is.null(aliases) || !nrow(aliases)) return(data.frame(alias = character(), canonical_term = character()))
  x <- gsub("[^a-z0-9]+", " ", tolower(enc2utf8(question)))
  padded <- paste0(" ", trimws(x), " ")
  hit <- vapply(aliases$alias, function(a) grepl(paste0(" ", a, " "), padded, fixed = TRUE), logical(1))
  out <- unique(aliases[hit & aliases$alias != aliases$canonical_term, c("alias", "canonical_term"), drop = FALSE])
  rownames(out) <- NULL
  out
}

.simitall_agent_normalised_line <- function(spec) {
  n <- spec$normalised
  if (is.null(n) || !nrow(n)) return(NULL)
  paste0("- Read as: ", paste(paste(n$alias, "->", n$canonical_term), collapse = "; "))
}

.simitall_agent_extract_spec <- function(question, data_source = "auto",
                                         input_files = character()) {
  # Extraction reads the user's words plus the canonical terms of any matched
  # synonyms; the spec keeps the original question.
  normalised <- .simitall_agent_alias_matches(question)
  augmented <- paste(question, paste(unique(normalised$canonical_term), collapse = " "))
  text <- tolower(augmented)
  organism <- .simitall_agent_species(augmented)
  catalog <- .simitall_panel_catalog()
  family <- .simitall_agent_extract_family_risk(augmented, organism, data_source, input_files)
  if (!is.null(family)) {
    family$question <- question
    family$normalised <- normalised
    return(family)
  }
  peak <- .simitall_agent_extract_peak_study(augmented, organism, data_source, input_files)
  if (!is.null(peak)) {
    peak$question <- question
    peak$normalised <- normalised
    return(peak)
  }
  single_cell <- .simitall_agent_extract_single_cell(augmented, organism, data_source, input_files)
  if (!is.null(single_cell)) {
    single_cell$question <- question
    single_cell$normalised <- normalised
    return(single_cell)
  }
  breeding <- .simitall_agent_extract_breeding(augmented, organism, data_source, input_files)
  if (!is.null(breeding)) {
    breeding$question <- question
    breeding$normalised <- normalised
    return(breeding)
  }

  # Plain-language gene hunting ("find the gene for X", "which genes control
  # X", "map the genes behind X") is a GWAS request too.
  gene_hunt <- "\\bgwas\\b|association (?:study|mapping|scan)|marker[- ]trait association|\\bqtl (?:mapping|scan)|(?:find|locate|identify|discover|detect|map)\\w* (?:the |a |any )?(?:causal )?(?:genes?|loci|locus|qtls?|variants?|snps?|markers?)|which (?:genes?|loci|regions?|variants?|snps?|markers?)|(?:genes?|loci|variants?) (?:for|behind|controlling|that control|responsible for|underlying)"
  # The methods a prompt names decide which stages run. A generic gene hunt
  # with no method named means GWAS; "using RNA-seq" means RNA-seq only.
  methods <- c(
    gwas = grepl("\\bgwas\\b|association (?:study|mapping|scan)|marker[- ]trait association|\\bqtl (?:mapping|scan)", text, perl = TRUE),
    rnaseq = !grepl("single[ -]?cell|scrna", text) &&
      grepl("rna.?seq|transcriptom|gene expression|expression|eqtl|differential(?:ly)? expressed|\\bdegs?\\b", text, perl = TRUE),
    genomic_selection = grepl("genomic.?selection|genomic prediction|gblup|breeding values?|(?:select|pick|choose|identify|find)(?:ing)? (?:the )?(?:best|top|elite|right|ideal)? ?(?:parents|parent lines|lines to cross|plants to cross)|best (?:parents|lines|plants) to (?:cross|breed)|which (?:lines|plants|parents) to (?:cross|keep|select|breed)", text, perl = TRUE)
  )
  analyses <- names(methods)[methods]
  # A gene hunt adds GWAS unless RNA-seq is the named way to find the genes.
  if (grepl(gene_hunt, text, perl = TRUE) && !methods[["rnaseq"]]) analyses <- union("gwas", analyses)
  task <- if (!length(analyses)) {
    NA_character_
  } else if (identical(analyses, "gwas")) {
    "gwas"
  } else {
    "population_study"
  }

  panel_id <- NA_character_
  design <- NA_character_
  if (grepl("\\bnam\\b|nested association", text, perl = TRUE)) {
    design <- "NAM"
    # The only NAM resource is maize. Record the inference instead of hiding it.
    if (identical(organism, "unspecified")) organism <- "maize"
    if (identical(organism, "maize")) panel_id <- "maize_nam"
  } else if (grepl("\\b1001\\b", text, perl = TRUE) && identical(organism, "arabidopsis")) {
    design <- "diversity_panel"
    panel_id <- "arabidopsis_1001"
  } else if (grepl("diversity panel|association panel|accessions|natural population|wild population", text)) {
    design <- "diversity_panel"
  }
  # Maize with no population named: use the maize breeding population simitall
  # ships (NAM), and record it as an assumption the plan shows.
  design_defaulted <- FALSE
  if (is.na(design) && identical(organism, "maize")) {
    design <- "NAM"; panel_id <- "maize_nam"; design_defaulted <- TRUE
  }

  # Analyses the composable compiler cannot build yet go to the legacy router.
  # Explicit environment contrasts keep using the versioned GxE multi-omics runner.
  other_analyses <- c(
    if (grepl("single[ -]?cell|scrna", text)) "scrnaseq",
    if (grepl("chip.?seq", text)) "chipseq",
    if (grepl("gxe|genotype.by.environment|irrigated|well[ -]?watered|multi.?omics|two environments|both environments", text)) "gxe",
    if (grepl("ancestry|admix", text)) "ancestry",
    if (grepl("pedigree|cousin|consanguin", text)) "pedigree"
  )

  trait_type <- if (grepl("binary|case[ -]?control|presence[ -]?absence|bernoulli", text)) {
    "binary"
  } else {
    "quantitative"
  }
  trait_names <- c(
    "flowering time" = "flowering_time", "drought" = "drought_response",
    "plant height" = "plant_height", "height" = "plant_height",
    "yield" = "yield", "disease resistance" = "disease_resistance"
  )
  trait_hit <- names(trait_names)[vapply(names(trait_names), grepl, logical(1), x = text, fixed = TRUE)]
  trait_name <- if (length(trait_hit)) unname(trait_names[trait_hit[1L]]) else "trait"

  heritability <- .simitall_agent_first_number(
    text, "(?:heritability|h2|h\\^2)\\s*(?:of|=|:|is)?\\s*(0?\\.[0-9]+|[01](?:\\.0+)?)"
  )
  # Panel names contain numbers ("1001 Genomes"); never read them as sizes.
  count_text <- gsub("\\b1001(?:[ -]genomes?)?\\b", " ", text, perl = TRUE)
  n_individuals <- .simitall_agent_first_number(
    count_text, "([0-9][0-9,]*)[\\s-]*(?:lines?|accessions?|individuals?|samples?|offspring|plants?)\\b"
  )
  n_qtl <- .simitall_agent_first_number(text, "([0-9]+)\\s*(?:qtls?|causal (?:variants?|loci|markers?))")
  seed <- .simitall_agent_first_number(text, "seed\\s*(?:=|of|:)?\\s*([0-9]+)")

  conditions <- if (grepl("drought", text)) {
    c("irrigated", "drought")
  } else if (grepl("heat", text)) {
    c("control", "heat")
  } else if (grepl("cold|chill|freez", text)) {
    c("control", "cold")
  } else if (grepl("salt|salin", text)) {
    c("control", "salt")
  } else {
    c("control", "treatment")
  }
  n_genes <- .simitall_agent_first_number(text, "([0-9][0-9,]*)\\s*genes\\b")

  list(
    spec_version = "1",
    extractor = "keyword_v1",
    question = question,
    task = task,
    analyses = analyses,
    organism = organism,
    population = list(design = design, panel_id = panel_id),
    other_analyses = other_analyses,
    rnaseq = if ("rnaseq" %in% analyses) list(conditions = conditions, n_genes = n_genes %||% 40L),
    normalised = normalised,
    trait = list(
      type = trait_type,
      name = trait_name,
      heritability = heritability %||% 0.6,
      n_qtl = n_qtl
    ),
    n_individuals = n_individuals,
    seed = seed %||% 2026,
    data_source = data_source,
    input_files = as.character(input_files %||% character()),
    # Which values came from the question versus defaults. This is what an
    # eval scores and what the user should see.
    defaulted = c(
      if (is.null(heritability)) "trait.heritability",
      if (is.null(n_individuals)) "n_individuals",
      if (is.null(seed)) "seed",
      if (design_defaulted) "population.design"
    )
  )
}

# Family-risk questions: "what is the chance a child of first cousins gets X".
# Returns NULL when the question is not one, so GWAS extraction runs instead.
.simitall_agent_extract_family_risk <- function(question, organism, data_source, input_files) {
  text <- tolower(question)
  if (grepl("\\bgwas\\b|association (?:study|mapping|scan)|(?:find|locate|identify|map)\\w* (?:the |a )?(?:genes?|loci|qtls?)", text, perl = TRUE)) return(NULL)
  if (!organism %in% c("human", "unspecified")) return(NULL)
  relationship_hit <- grepl("cousin|consanguin|inbre|unrelated|related (?:parents|couples?|people|partners)|marr(?:y|ied|iage)|\\bkin\\b|relatives?|uncle|aunt|niece|nephew|\\bsiblings?\\b|brother|sister", text, perl = TRUE)
  risk_hit <- grepl("probabilit|risk|chance|likel|odds|disease|disorder|condition|syndrome|defect|trait", text)
  if (!relationship_hit || !risk_hit) return(NULL)

  other_analyses <- c(
    if (grepl("rna.?seq|eqtl|expression", text)) "rnaseq",
    if (grepl("ancestry|admix", text)) "ancestry"
  )
  unsupported_relationships <- c(
    sibling = "\\bsiblings?\\b|brothers?\\b|sisters?\\b",
    uncle_niece = "uncle|aunt|niece|nephew",
    half_relatives = "half[ -]",
    third_cousin = "third[ -]?cousin",
    double_cousin = "double[ -]?cousin"
  )
  unsupported <- names(unsupported_relationships)[vapply(unsupported_relationships, grepl, logical(1), x = text, perl = TRUE)]
  plain_cousin <- grepl("\\bcousins?\\b", gsub("(?:second|third|distant(?:ly)?|double)[ -]?(?:related )?cousins?", " ", text, perl = TRUE), perl = TRUE)
  groups <- c(
    if (grepl("second[ -]?cousin|distant", text)) "distantly_related",
    if (grepl("first[ -]?cousin", text) || plain_cousin) "first_cousin_descendant"
  )
  defaulted <- character()
  if (!length(groups)) {
    groups <- c("distantly_related", "first_cousin_descendant")
    defaulted <- c(defaulted, "groups")
  }
  # Relative risk needs a reference, so unrelated parents are always included.
  groups <- c("unrelated", setdiff(groups, "unrelated"))
  groups <- intersect(c("unrelated", "distantly_related", "first_cousin_descendant"), groups)

  known <- c(
    "cleft lip|cleft palate|orofacial cleft" = "cleft lip",
    "cystic fibrosis" = "cystic fibrosis", "sickle[ -]cell" = "sickle cell disease",
    "thalass" = "thalassemia", "tay[ -]?sachs" = "Tay-Sachs disease",
    "phenylketonuria|\\bpku\\b" = "phenylketonuria", "deaf|hearing loss" = "hearing loss",
    "heart defect|congenital heart" = "congenital heart defect", "birth defect" = "birth defect",
    "spinal muscular atrophy|\\bsma\\b" = "spinal muscular atrophy", "albinism" = "albinism"
  )
  hit <- names(known)[vapply(names(known), grepl, logical(1), x = text, perl = TRUE)]
  condition <- if (length(hit)) unname(known[hit[1L]]) else if (grepl("disease", text)) "disease" else "condition"

  prevalence <- NULL
  one_in <- regmatches(text, regexec("1\\s*(?:in|/|out of)\\s*([0-9][0-9,]*)", text, perl = TRUE))[[1L]]
  if (length(one_in) > 1L) prevalence <- 1 / as.numeric(gsub(",", "", one_in[2L]))
  if (is.null(prevalence)) {
    percent <- .simitall_agent_first_number(text, "([0-9]+(?:\\.[0-9]+)?)\\s*%")
    if (!is.null(percent)) prevalence <- percent / 100
  }
  if (is.null(prevalence)) {
    prevalence <- .simitall_agent_first_number(text, "(?:prevalence|baseline(?: risk| probability)?)\\s*(?:of|=|:|is)?\\s*(0?\\.[0-9]+)")
  }
  n <- .simitall_agent_first_number(
    text, "([0-9][0-9,]*)\\s*(?:families|children|kids|people|individuals|descendants|samples|couples|births)\\b"
  )
  seed <- .simitall_agent_first_number(text, "seed\\s*(?:=|of|:)?\\s*([0-9]+)")
  model <- if (grepl("dominant", text)) "dominant" else if (grepl("additive|polygenic", text)) "additive" else "recessive"

  list(
    spec_version = "1",
    extractor = "keyword_v1",
    question = question,
    task = "family_risk",
    organism = "human",
    population = list(design = "relationship_groups",
                      panel_id = if (grepl("\\birf6\\b", text, perl = TRUE)) "human_irf6" else NA_character_),
    groups = groups,
    unsupported_relationships = unsupported,
    other_analyses = other_analyses,
    condition = list(
      name = condition,
      prevalence = prevalence %||% 0.01,
      model = model,
      n_risk_loci = 10L,
      risk_allele_frequency = 0.05,
      risk_effect = 4
    ),
    n_individuals = n,
    seed = seed %||% 81,
    data_source = data_source,
    input_files = as.character(input_files %||% character()),
    defaulted = c(
      defaulted,
      if (is.null(prevalence)) "condition.prevalence",
      if (is.null(n)) "n_individuals",
      if (is.null(seed)) "seed",
      "condition.risk_loci"
    )
  )
}

.simitall_agent_validate_family_risk <- function(spec) {
  errors <- character()
  add <- function(...) errors <<- c(errors, paste0(...))
  is_whole <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)
  if (length(spec$unsupported_relationships)) {
    add("groups: ", paste(spec$unsupported_relationships, collapse = ", "),
        " not supported yet; available: unrelated parents, second cousins, first cousins")
  }
  if (!identical(spec$condition$model, "recessive")) {
    add("condition.model: '", spec$condition$model, "' inheritance is not supported; relationship comparisons use a recessive model")
  }
  prev <- spec$condition$prevalence
  if (!is.numeric(prev) || !is.finite(prev) || prev <= 0 || prev >= 0.5) {
    add("condition.prevalence: must be between 0 and 0.5; got ", format(prev))
  }
  if (!is.null(spec$n_individuals) && (!is_whole(spec$n_individuals) || spec$n_individuals < 20 || spec$n_individuals > 20000)) {
    add("n_individuals: children per group must be a whole number in [20, 20000]; got ", format(spec$n_individuals))
  }
  k <- spec$condition$n_risk_loci
  if (!is_whole(k) || k < 1 || k > 50) add("condition.n_risk_loci: must be a whole number in [1, 50]; got ", format(k))
  q <- spec$condition$risk_allele_frequency
  if (!is.numeric(q) || !is.finite(q) || q <= 0 || q >= 0.5) add("condition.risk_allele_frequency: must be between 0 and 0.5; got ", format(q))
  b <- spec$condition$risk_effect
  if (!is.numeric(b) || !is.finite(b) || b <= 0 || b > 10) add("condition.risk_effect: must be in (0, 10]; got ", format(b))
  if (!is_whole(spec$seed) || spec$seed < 1) add("seed: must be a positive whole number")
  if (!spec$data_source %in% c("auto", "synthetic", "provided", "package_demo")) {
    add("data_source: '", spec$data_source, "' is not one of auto, synthetic, provided, package_demo")
  }
  missing <- spec$input_files[nzchar(spec$input_files) & !file.exists(path.expand(spec$input_files))]
  if (length(missing)) add("input_files: these files do not exist: ", paste(missing, collapse = ", "))
  errors
}

# Marker source: a named bundled region, the user's VCF, or simulated
# anonymous markers. Risk loci are always simulated, never real variants.
.simitall_agent_resolve_family_risk <- function(spec) {
  files <- spec$input_files[nzchar(spec$input_files)]
  if (identical(spec$data_source, "provided")) {
    vcf <- files[grepl("\\.(vcf|bcf)(\\.gz)?$", files, ignore.case = TRUE)][1L]
    if (is.na(vcf)) {
      return(list(status = "needs_input", mode = NA_character_, files = files,
                  action = "Provide a human VCF with INFO/AF allele frequencies, or choose 'Generate synthetic inputs'."))
    }
    return(list(status = "ready", mode = "provided", files = files, vcf = normalizePath(path.expand(vcf)),
                note = "Neutral marker frequencies come from your VCF; disease risk loci are simulated."))
  }
  if (identical(spec$population$panel_id, "human_irf6") || identical(spec$data_source, "package_demo")) {
    return(list(status = "ready", mode = "bundled", files = files,
                note = "Neutral marker frequencies come from the small synthetic IRF6-region VCF bundled with simitall; disease risk loci are simulated, not real IRF6 variants."))
  }
  list(status = "ready", mode = "simulate", files = files,
       note = "No data set was named, so anonymous marker frequencies and disease risk loci are both simulated.")
}

.simitall_agent_compile_family_risk <- function(spec) {
  d <- .simitall_agent_spec_defaults(spec)
  cond <- spec$condition
  condition_id <- gsub("[^a-z0-9]+", "_", tolower(cond$name))
  out_dir <- file.path("results", "agent_family_risk", condition_id)
  marker_block <- switch(
    spec$resolution$mode,
    simulate = "markers <- simulate_human_marker_frequencies(n_markers = 24L, seed = params$seed)",
    bundled = c(
      "vcf <- system.file(\"extdata\", \"human_irf6\", \"demo_human_irf6_chr1.vcf.gz\", package = \"simitall\")",
      "stopifnot(nzchar(vcf), file.exists(vcf))",
      "markers <- read_human_marker_frequencies(vcf, maf_min = 0.10, maf_max = 0.40)"
    ),
    provided = c(
      paste0("vcf <- ", encodeString(spec$resolution$vcf, quote = "\"")),
      "markers <- read_human_marker_frequencies(vcf, maf_min = 0.10, maf_max = 0.40)"
    )
  )
  code <- c(
    "if (!\"simitall\" %in% loadedNamespaces()) library(simitall)",
    "params <- list(",
    paste0("  seed = ", as.integer(spec$seed), ","),
    paste0("  n_per_group = ", as.integer(spec$n_individuals %||% d$n_individuals), "L,"),
    paste0("  baseline_probability = ", signif(cond$prevalence, 6), ","),
    paste0("  n_risk_loci = ", as.integer(cond$n_risk_loci), "L,"),
    paste0("  risk_allele_frequency = ", cond$risk_allele_frequency, ","),
    paste0("  risk_effect = ", cond$risk_effect),
    ")",
    paste0("groups <- c(", paste(sprintf("\"%s\"", spec$groups), collapse = ", "), ")"),
    paste0("out_dir <- ", encodeString(out_dir, quote = "\"")),
    "figure_dir <- file.path(out_dir, \"figures\")",
    "dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)",
    "",
    "cat(\"1. Preparing marker frequencies...\\n\")",
    marker_block,
    "",
    "cat(\"2. Simulating families for each relationship group...\\n\")",
    "run <- simulate_human_pedigree_groups(",
    "  markers, groups = groups, n_per_group = params$n_per_group, seed = params$seed,",
    "  baseline_probability = params$baseline_probability, risk_model = \"recessive\",",
    "  n_risk_loci = params$n_risk_loci, risk_allele_frequency = params$risk_allele_frequency,",
    "  risk_effect = params$risk_effect",
    ")",
    "write.table(run$pedigree, file.path(out_dir, \"pedigree_truth.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "write.table(run$individuals, file.path(out_dir, \"children.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "write.table(run$marker_truth, file.path(out_dir, \"marker_and_risk_locus_truth.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "",
    "cat(\"3. Comparing simulated risk with the exact model expectation...\\n\")",
    "summary <- summarize_family_risk(run)",
    "write.table(summary, file.path(out_dir, \"family_risk_summary.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "plot_human_pedigree_diagnostics(run$individuals, figure_dir)",
    paste0("plot_family_risk(summary, file.path(figure_dir, \"family_risk_by_relationship.png\"), condition = ",
           encodeString(cond$name, quote = "\""), ")"),
    "stopifnot(all(summary$agrees_with_theory))  # simulation agrees with the exact expectation",
    "print(summary)"
  )
  paste(
    "Spec-compiled family-risk recipe:",
    "",
    "```r",
    paste(code, collapse = "\n"),
    "```",
    sep = "\n"
  )
}

.simitall_agent_family_risk_plan <- function(spec) {
  d <- .simitall_agent_spec_defaults(spec)
  cond <- spec$condition
  was_default <- function(field) if (field %in% spec$defaulted) " (default)" else ""
  n <- as.integer(spec$n_individuals %||% d$n_individuals)
  pretty <- c(unrelated = "unrelated parents", distantly_related = "second cousins (F = 1/64)",
              first_cousin_descendant = "first cousins (F = 1/16)")
  group_text <- paste(unname(pretty[spec$groups]), collapse = ", ")
  mode <- spec$resolution$mode %||% NA_character_
  mode_key <- if (is.na(mode)) "" else mode
  source_text <- switch(mode_key, simulate = "simulated anonymous markers",
                        bundled = "bundled IRF6-region marker frequencies",
                        provided = "your VCF's marker frequencies", "not resolved yet")
  pct <- function(x) paste0(format(signif(100 * x, 3)), "%")
  design <- c(
    "- Task: probability of a condition by parents' relationship",
    paste0("- Condition: ", cond$name, " (a label: the risk loci are simulated, not real variants)"),
    paste0("- Groups: ", group_text, if ("groups" %in% spec$defaulted) " (default comparison)" else ""),
    paste0("- Children per group: ", n, was_default("n_individuals")),
    paste0("- Risk for children of unrelated parents: ", pct(cond$prevalence), was_default("condition.prevalence")),
    paste0("- Model: recessive, ", cond$n_risk_loci, " risk loci, risk-allele frequency ", cond$risk_allele_frequency,
           ", +", cond$risk_effect, " log-odds per locus with two risk copies", was_default("condition.risk_loci")),
    paste0("- Data source: ", spec$data_source, " -> ", source_text),
    paste0("- Seed: ", spec$seed, was_default("seed")),
    paste0("- Parsed by: ", spec$extractor),
    .simitall_agent_normalised_line(spec)
  )
  stages <- c(
    paste0("1. Prepare neutral marker frequencies from ", source_text, ", plus ", cond$n_risk_loci, " simulated recessive risk loci."),
    paste0("2. Build ", n, " families per group (", group_text, ") and pass every allele from parents to children by Mendelian inheritance."),
    "3. Score each child from the number of risk loci where they inherited two risk copies. The model has no relatedness term, so any group difference comes from inheritance.",
    "4. Compare each group's simulated risk with the exact theoretical expectation, then plot risk and marker homozygosity by group."
  )
  inputs <- if (length(spec$errors)) {
    paste0("- Fix before running: ", spec$errors)
  } else if (identical(spec$resolution$status, "needs_input")) {
    paste0("- ", spec$resolution$action)
  } else if (identical(mode, "provided")) {
    paste0("- ", spec$resolution$files)
  } else {
    "- None: everything needed is simulated or bundled."
  }
  outputs <- c(
    "- Pedigree truth for every family, each child's genotypes, risk-locus status, and probability.",
    "- family_risk_summary: expected vs simulated risk per group, relative risk vs unrelated parents, and a pass/fail truth check.",
    "- Figures: risk by relationship, marker homozygosity, probability distribution."
  )
  limitations <- c(
    if (!is.null(spec$resolution$note)) paste0("- ", spec$resolution$note),
    "- Risk-locus count, allele frequency, and effect size are modelling assumptions. The relative risks are a consequence of those assumptions, not estimates for any real condition or family.",
    "- This is a simulation for method testing and teaching. It is not a clinical risk estimate or genetic counselling."
  )
  paste(
    "simitall spec-based plan (the plan and code below come from the same validated spec).",
    paste0("Question: ", spec$question),
    paste0("Interpreted design:\n", paste(design, collapse = "\n")),
    paste0("Proposed stages:\n", paste(stages, collapse = "\n")),
    paste0("Required inputs:\n", paste(inputs, collapse = "\n")),
    paste0("Truth and QC outputs:\n", paste(outputs, collapse = "\n")),
    paste0("Limitations:\n", paste(limitations, collapse = "\n")),
    sep = "\n\n"
  )
}

# ---------------------------------------------------------------------------
# Peak studies (ChIP-seq, ATAC-seq): simulate, call peaks, test differential
# binding or accessibility, and score both against the simulated truth.
# ---------------------------------------------------------------------------

.simitall_agent_conditions <- function(text) {
  if (grepl("drought", text)) c("control", "drought")
  else if (grepl("heat", text)) c("control", "heat")
  else if (grepl("cold|chill|freez", text)) c("control", "cold")
  else if (grepl("salt|salin", text)) c("control", "salt")
  else c("control", "treatment")
}

.simitall_agent_extract_peak_study <- function(question, organism, data_source, input_files) {
  text <- tolower(question)
  atac <- grepl("atac|chromatin accessib|open chromatin|accessible (?:chromatin|regions?)|accessibility", text, perl = TRUE)
  chip <- grepl("chip.?seq|chromatin immunoprecip|binding sites?|histone (?:marks?|modifications?)|\\bh3k[0-9]+(?:me[0-9]|ac)\\b|transcription factor binding|tf binding", text, perl = TRUE)
  if (!atac && !chip) return(NULL)
  mark <- regmatches(text, regexpr("\\bh3k[0-9]+(?:me[0-9]|ac)\\b", text, perl = TRUE))
  supported_marks <- c(h3k4me3 = "H3K4me3", h3k27ac = "H3K27ac", h3k27me3 = "H3K27me3", h3k36me3 = "H3K36me3")
  number <- function(pattern) .simitall_agent_first_number(text, pattern)
  n_peaks <- number("([0-9][0-9,]*)\\s*(?:true |binding |accessible )?(?:peaks?|binding sites?|regions?)\\b")
  replicates <- number("([0-9]+)\\s*(?:biological )?replicates?")
  reads <- .simitall_agent_first_number(text, "([0-9][0-9,]*)\\s*reads?\\b")
  seed <- number("seed\\s*(?:=|of|:)?\\s*([0-9]+)")
  histone <- chip && (length(mark) > 0L || grepl("histone", text))
  list(
    spec_version = "1", extractor = "keyword_v1", question = question,
    task = "peak_study",
    assay = if (atac && chip) "both" else if (atac) "atacseq" else "chipseq",
    target = if (!chip) NA_character_ else if (histone) "histone" else "TF",
    histone_mark = if (length(mark)) {
      if (mark %in% names(supported_marks)) unname(supported_marks[mark]) else sub("^h3k", "H3K", mark)
    } else if (histone) "H3K27ac" else NA_character_,
    organism = organism,
    population = list(design = "peak_study", panel_id = NA_character_),
    other_analyses = c(if (grepl("single[ -]?cell|scatac|sc-atac", text)) "single_cell"),
    conditions = .simitall_agent_conditions(text),
    n_peaks = n_peaks %||% 100L,
    replicates = replicates %||% 3L,
    n_reads = reads %||% 40000L,
    seed = seed %||% 2026,
    data_source = data_source,
    input_files = as.character(input_files %||% character()),
    defaulted = c(if (is.null(n_peaks)) "n_peaks", if (is.null(replicates)) "replicates",
                  if (is.null(reads)) "n_reads", if (is.null(seed)) "seed",
                  if (histone && !length(mark)) "histone_mark")
  )
}

.simitall_agent_validate_peak_study <- function(spec) {
  errors <- character()
  add <- function(...) errors <<- c(errors, paste0(...))
  is_whole <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)
  if (identical(spec$assay, "both")) add("assay: ask for ChIP-seq or ATAC-seq in separate requests")
  if (identical(spec$target, "histone") && !spec$histone_mark %in% c("H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3")) {
    add("histone_mark: '", spec$histone_mark, "' is not simulated; available: H3K4me3, H3K27ac, H3K27me3, H3K36me3")
  }
  if (!is_whole(spec$n_peaks) || spec$n_peaks < 10 || spec$n_peaks > 2000) add("n_peaks: must be a whole number in [10, 2000]; got ", format(spec$n_peaks))
  if (!is_whole(spec$replicates) || spec$replicates < 2 || spec$replicates > 10) add("replicates: must be a whole number in [2, 10] (differential tests need at least 2 per condition); got ", format(spec$replicates))
  if (!is_whole(spec$n_reads) || spec$n_reads < 5000 || spec$n_reads > 2000000) add("n_reads: reads per library must be in [5,000, 2,000,000] for a local run; got ", format(spec$n_reads))
  if (!is_whole(spec$seed) || spec$seed < 1) add("seed: must be a positive whole number")
  missing <- spec$input_files[nzchar(spec$input_files) & !file.exists(path.expand(spec$input_files))]
  if (length(missing)) add("input_files: these files do not exist: ", paste(missing, collapse = ", "))
  errors
}

.simitall_agent_resolve_peak_study <- function(spec) {
  files <- spec$input_files[nzchar(spec$input_files)]
  fasta <- files[grepl("\\.(fa|fasta|fna)$", files, ignore.case = TRUE)][1L]
  gff <- files[grepl("\\.gff3?$", files, ignore.case = TRUE)][1L]
  if (identical(spec$data_source, "provided")) {
    if (is.na(fasta)) {
      return(list(status = "needs_input", mode = NA_character_, files = files,
                  action = "Provide a reference FASTA (and optionally a GFF3 with promoter/enhancer/gene features), or choose 'Generate synthetic inputs'."))
    }
    return(list(status = "ready", mode = "provided", files = files, fasta = normalizePath(path.expand(fasta)),
                gff3 = if (is.na(gff)) NA_character_ else normalizePath(path.expand(gff)),
                note = "Peaks are simulated on your reference; with a GFF3 they target its regulatory features."))
  }
  list(status = "ready", mode = "simulate", files = files,
       note = "No reference was supplied, so a small synthetic genome with genes, promoters, and enhancers is generated for the simulation.")
}

.simitall_agent_peak_defaults <- function(spec) {
  installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)
  list(peak_fdr = 0.05,
       diff_method = if (installed("DESeq2")) "deseq2" else if (installed("edgeR")) "edger" else "ttest")
}

.simitall_agent_compile_peak_study <- function(spec) {
  atac <- identical(spec$assay, "atacseq")
  label <- if (atac) "atac" else "chip"
  q <- function(x) encodeString(x, quote = "\"")
  cond <- spec$conditions
  broad <- identical(spec$target, "histone")
  width <- if (broad && identical(spec$histone_mark, "H3K27me3")) 10000L else if (broad && identical(spec$histone_mark, "H3K36me3")) 5000L else if (broad) 1500L else 300L
  window <- if (broad) 1000L else 200L
  # Size the synthetic genome so peaks occupy at most ~10% of it.
  chrom_len <- as.integer(max(200000, ceiling(spec$n_peaks * width * 10 / 2 / 1000) * 1000))
  genome <- if (identical(spec$resolution$mode, "provided")) {
    paste0("genome <- list(fasta = ", q(spec$resolution$fasta), ", gff3 = ",
           if (is.na(spec$resolution$gff3)) "NULL" else q(spec$resolution$gff3), ")")
  } else {
    c("genome <- simitall:::.simitall_simulate_peak_genome(",
      "  file.path(out_dir, \"genome.fa\"), file.path(out_dir, \"genome.gff3\"),",
      paste0("  chromosome_length = ", chrom_len, "L, n_genes = ", max(120L, as.integer(spec$n_peaks)), "L, seed = params$seed"),
      ")")
  }
  simulate <- if (atac) {
    c("sim <- simulate_atacseq(",
      "  genome$fasta, file.path(out_dir, \"atac\"), annotation_gff3 = genome$gff3, n_peaks = params$n_peaks,",
      paste0("  conditions = c(", q(cond[1L]), ", ", q(cond[2L]), "), biological_replicates = params$replicates,"),
      "  n_reads = params$n_reads, differential_fraction = 0.2, differential_effect_log2 = 1.5, seed = params$seed",
      ")")
  } else {
    c("sim <- simulate_chipseq(",
      "  genome$fasta, file.path(out_dir, \"chip\"), annotation_gff3 = genome$gff3,",
      paste0("  assay_type = ", q(spec$target), if (broad) paste0(", histone_mark = ", q(spec$histone_mark)) else "", ", n_peaks = params$n_peaks,"),
      paste0("  conditions = c(", q(cond[1L]), ", ", q(cond[2L]), "), biological_replicates = params$replicates,"),
      "  n_reads = params$n_reads, differential_binding_fraction = 0.2, differential_effect_log2 = 1.5,",
      "  write_fastq = FALSE, write_read_bed = FALSE, backend = \"native\", seed = params$seed",
      ")")
  }
  kind <- if (atac) "accessibility" else "binding"
  code <- c(
    "if (!\"simitall\" %in% loadedNamespaces()) library(simitall)",
    "params <- list(",
    paste0("  seed = ", as.integer(spec$seed), ","),
    paste0("  n_peaks = ", as.integer(spec$n_peaks), "L,"),
    paste0("  replicates = ", as.integer(spec$replicates), "L,"),
    paste0("  n_reads = ", as.integer(spec$n_reads), "L"),
    ")",
    paste0("out_dir <- ", q(file.path("results", "agent_peaks", if (atac) "atacseq" else paste0("chipseq_", tolower(spec$histone_mark %|na|% "tf"))))),
    "figure_dir <- file.path(out_dir, \"figures\")",
    "dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)",
    "",
    "cat(\"1. Preparing the reference genome and regulatory annotation...\\n\")",
    genome,
    paste0("cat(\"2. Simulating ", if (atac) "ATAC-seq" else "ChIP-seq", " libraries with known ", kind, " truth...\\n\")"),
    simulate,
    "reads <- read.delim(sim$read_positions, stringsAsFactors = FALSE)",
    "",
    "cat(\"3. Calling peaks and scoring them against the true peaks...\\n\")",
    paste0("called <- call_peaks_windowed(reads, genome$fasta, signal_assay = ", q(label), ", control_assay = ",
           if (atac) "NULL" else "\"input\"", ", window = ", window, "L, fdr = ", .simitall_agent_choice(spec, "peak_fdr"), ")"),
    "write.table(called, file.path(out_dir, \"called_peaks.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "peak_metrics <- score_peak_calls(called, sim$peaks)",
    "write.table(peak_metrics, file.path(out_dir, \"peak_calling_truth_recovery_metrics.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "",
    paste0("cat(\"4. Testing differential ", kind, " and scoring it against the truth...\\n\")"),
    paste0("da <- test_differential_peaks(reads, called, assay = ", q(label), ", method = ", q(.simitall_agent_choice(spec, "diff_method")), ", fdr = 0.05)"),
    "write.table(da, file.path(out_dir, \"differential_peaks.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "diff_metrics <- score_differential_peaks(da, called, sim$peaks, sim$differential_truth)",
    "write.table(diff_metrics, file.path(out_dir, \"differential_truth_recovery_metrics.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "",
    "truth <- read.delim(sim$peaks); diff_truth <- read.delim(sim$differential_truth)",
    "true_diff <- truth[truth$peak_id %in% diff_truth$peak_id, ]",
    "is_true <- vapply(seq_len(nrow(called)), function(i) any(true_diff$seqname == called$seqname[i] & true_diff$start <= called$end[i] & true_diff$end >= called$start[i]), logical(1))",
    "png(file.path(figure_dir, \"differential_volcano.png\"), width = 1400, height = 1000, res = 170)",
    "plot(da$log2_fold_change, -log10(pmax(da$p_value, 1e-300)), pch = 19, col = ifelse(is_true, \"#B85252\", \"#9AA7A0\"),",
    paste0("  xlab = \"log2 fold change\", ylab = \"-log10 p\", main = ", q(paste0("Differential ", kind, ": red = truly differential (simulation truth)")), ", cex.main = 0.9)"),
    "abline(h = -log10(max(c(0, da$p_value[da$significant]))), lty = 2)",
    "dev.off()",
    if (atac) c(
      "png(file.path(figure_dir, \"fragment_lengths.png\"), width = 1400, height = 900, res = 170)",
      "hist(reads$fragment_length, breaks = 80, col = \"#7A9CC6\", border = NA, main = \"ATAC fragment lengths: nucleosome-free and mono-nucleosome\", xlab = \"Fragment length (bp)\", cex.main = 0.9)",
      "dev.off()"),
    "print(peak_metrics); print(diff_metrics)"
  )
  paste(paste0("Spec-compiled ", if (atac) "ATAC-seq" else "ChIP-seq", " recipe:"), "", "```r", paste(code, collapse = "\n"), "```", sep = "\n")
}

.simitall_agent_peak_study_plan <- function(spec) {
  atac <- identical(spec$assay, "atacseq")
  was_default <- function(field) if (field %in% spec$defaulted) " (default)" else ""
  mode <- spec$resolution$mode %||% NA_character_
  assay_name <- if (atac) "ATAC-seq (chromatin accessibility)" else if (identical(spec$assay, "both")) "ChIP-seq and ATAC-seq" else
    paste0("ChIP-seq, ", if (identical(spec$target, "histone")) paste0("histone mark ", spec$histone_mark, was_default("histone_mark")) else "transcription-factor binding")
  kind <- if (atac) "accessibility" else "binding"
  design <- c(
    paste0("- Assay: ", assay_name),
    paste0("- Conditions: ", spec$conditions[1L], " vs ", spec$conditions[2L], ", ", spec$replicates, " biological replicates each", was_default("replicates")),
    paste0("- True peaks: ", spec$n_peaks, was_default("n_peaks"), ", 20% of them differential (about 1.5 log2-fold change)"),
    paste0("- Reads per library: ", format(spec$n_reads, big.mark = ","), was_default("n_reads"), " (scaled down for a fast local run)"),
    paste0("- Peak-calling FDR: ", .simitall_agent_choice(spec, "peak_fdr"), "; differential test: ", .simitall_agent_choice(spec, "diff_method")),
    paste0("- Reference: ", if (identical(mode, "provided")) "your FASTA/GFF3" else "synthetic genome with genes, promoters, and enhancers"),
    paste0("- Seed: ", spec$seed, was_default("seed")),
    paste0("- Parsed by: ", spec$extractor),
    .simitall_agent_normalised_line(spec)
  )
  stages <- c(
    "1. Prepare the reference genome and its regulatory annotation.",
    if (atac) "2. Simulate ATAC libraries: accessible regions at promoters/enhancers, nucleosome-free and mono-nucleosome fragments, Tn5 cut-site offsets, PCR duplicates."
    else "2. Simulate ChIP and matched input libraries with known binding peaks and condition-dependent binding.",
    paste0("3. Call peaks with a window-based Poisson test", if (atac) " against the genome background" else " against the input control", ", then score recall/precision against the true peaks."),
    paste0("4. Test differential ", kind, " with `", .simitall_agent_choice(spec, "diff_method"), "` and score it against the true differential peaks.")
  )
  inputs <- if (length(spec$errors)) paste0("- Fix before running: ", spec$errors) else if (identical(spec$resolution$status, "needs_input")) paste0("- ", spec$resolution$action) else "- None: everything needed is simulated."
  outputs <- c(
    "- True peaks and true differential effects; simulated read positions and per-library QC (FRiP, duplicates).",
    "- Called peaks with peak_calling_truth_recovery_metrics.",
    paste0("- Differential ", kind, " results with differential_truth_recovery_metrics and a volcano plot marking truly differential peaks."),
    if (atac) "- Fragment-length histogram showing the nucleosome-free and mono-nucleosome modes."
  )
  limitations <- c(
    if (!is.null(spec$resolution$note)) paste0("- ", spec$resolution$note),
    "- Peaks are called with a transparent window-based baseline caller, not MACS; real pipelines usually use MACS2/3.",
    "- Read depth is scaled down for speed; real ATAC libraries are typically around 10 million read pairs.",
    "- Recovering the programmed peaks shows the analysis works under these simulation assumptions. It is a benchmark, not a biological discovery."
  )
  paste(
    "simitall spec-based plan (the plan and code below come from the same validated spec).",
    paste0("Question: ", spec$question),
    paste0("Interpreted design:\n", paste(design, collapse = "\n")),
    paste0("Proposed stages:\n", paste(stages, collapse = "\n")),
    paste0("Required inputs:\n", paste(inputs, collapse = "\n")),
    paste0("Truth and QC outputs:\n", paste(outputs, collapse = "\n")),
    paste0("Limitations:\n", paste(limitations, collapse = "\n")),
    sep = "\n\n"
  )
}

# ---------------------------------------------------------------------------
# Single-cell RNA-seq: simulate donors and cells, cluster (scored against the
# true cell types), and test condition DE (scored against the true condition
# genes), with pseudobulk and per-cell tests available for comparison.
# ---------------------------------------------------------------------------

.simitall_agent_cell_types <- function(text, organism) {
  if (grepl("\\broots?\\b", text, perl = TRUE)) {
    return(list(tissue = "root", types = c(epidermis = 0.3, cortex = 0.3, endodermis = 0.2, stele = 0.2)))
  }
  if (identical(organism, "human") || grepl("blood|immune|pbmc|lymph", text)) {
    return(list(tissue = "blood", types = c(T_cell = 0.35, B_cell = 0.25, Monocyte = 0.25, NK_cell = 0.15)))
  }
  list(tissue = "leaf", types = c(mesophyll = 0.4, epidermis = 0.3, vascular = 0.2, guard_cell = 0.1))
}

.simitall_agent_extract_single_cell <- function(question, organism, data_source, input_files) {
  text <- tolower(question)
  if (!grepl("single[ -]?cell|single[ -]?nucle|\\bscrna|\\bsnrna|sc-rna|cell atlas|\\b10x\\b|cell types?\\b|cell populations?", text, perl = TRUE)) return(NULL)
  number <- function(pattern) .simitall_agent_first_number(text, pattern)
  donors <- number("([0-9]+)\\s*(?:donors?|plants|samples|individuals|biological replicates)")
  cells <- number("([0-9][0-9,]*)\\s*cells?\\b")
  genes <- number("([0-9][0-9,]*)\\s*genes\\b")
  seed <- number("seed\\s*(?:=|of|:)?\\s*([0-9]+)")
  tissue <- .simitall_agent_cell_types(text, organism)
  list(
    spec_version = "1", extractor = "keyword_v1", question = question,
    task = "single_cell", organism = organism,
    population = list(design = "single_cell", panel_id = NA_character_),
    other_analyses = c(if (grepl("eqtl|\\bgwas\\b|association", text, perl = TRUE)) "genetics",
                       if (grepl("trajector|pseudotime", text)) "trajectory"),
    tissue = tissue$tissue, cell_types = tissue$types,
    conditions = .simitall_agent_conditions(text),
    n_donors = donors %||% 6L,
    cells_per_donor = cells %||% 200L,
    n_genes = genes %||% 400L,
    seed = seed %||% 2026,
    data_source = data_source,
    input_files = as.character(input_files %||% character()),
    defaulted = c(if (is.null(donors)) "n_donors", if (is.null(cells)) "cells_per_donor",
                  if (is.null(genes)) "n_genes", if (is.null(seed)) "seed")
  )
}

.simitall_agent_validate_single_cell <- function(spec) {
  errors <- character()
  add <- function(...) errors <<- c(errors, paste0(...))
  is_whole <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)
  if (!is_whole(spec$n_donors) || spec$n_donors < 4 || spec$n_donors > 24) add("n_donors: must be a whole number in [4, 24] (at least 2 donors per condition for pseudobulk tests); got ", format(spec$n_donors))
  if (!is_whole(spec$cells_per_donor) || spec$cells_per_donor < 50 || spec$cells_per_donor > 2000) add("cells_per_donor: must be a whole number in [50, 2000]; got ", format(spec$cells_per_donor))
  if (!is_whole(spec$n_genes) || spec$n_genes < 100 || spec$n_genes > 3000) add("n_genes: must be a whole number in [100, 3000]; got ", format(spec$n_genes))
  if (is_whole(spec$n_donors) && is_whole(spec$cells_per_donor) && is_whole(spec$n_genes) &&
      spec$n_donors * spec$cells_per_donor * spec$n_genes > 6e7) {
    add("size: donors x cells x genes is too large for a local run (limit 60 million); reduce one of them")
  }
  if (!is_whole(spec$seed) || spec$seed < 1) add("seed: must be a positive whole number")
  errors
}

.simitall_agent_resolve_single_cell <- function(spec) {
  list(status = "ready", mode = "simulate", files = character(),
       note = "Donor genotypes, cells, cell types, and condition effects are all simulated; no data are needed.")
}

.simitall_agent_single_cell_defaults <- function(spec) {
  installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)
  list(hvg_fraction = 0.2,
       sc_de_method = if (installed("DESeq2")) "pseudobulk_deseq2" else if (installed("edgeR")) "pseudobulk_edger" else "pseudobulk_ttest")
}

.simitall_agent_compile_single_cell <- function(spec) {
  q <- function(x) encodeString(x, quote = "\"")
  cond <- spec$conditions
  types <- spec$cell_types
  type_code <- paste0("c(", paste0(names(types), " = ", types, collapse = ", "), ")")
  code <- c(
    "if (!\"simitall\" %in% loadedNamespaces()) library(simitall)",
    "params <- list(",
    paste0("  seed = ", as.integer(spec$seed), ","),
    paste0("  n_donors = ", as.integer(spec$n_donors), "L,"),
    paste0("  cells_per_donor = ", as.integer(spec$cells_per_donor), "L,"),
    paste0("  n_genes = ", as.integer(spec$n_genes), "L"),
    ")",
    paste0("out_dir <- ", q(file.path("results", "agent_singlecell", spec$tissue))),
    "figure_dir <- file.path(out_dir, \"figures\")",
    "dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)",
    "",
    "cat(\"1. Simulating donor genotypes...\\n\")",
    "donors <- simitall:::.simitall_simulate_inbred_panel(file.path(out_dir, \"donors.vcf\"), n_accessions = params$n_donors,",
    "  n_markers = 300L, n_chromosomes = 2L, chromosome_length_bp = 1000000L, n_subpops = 1L, fst = 0.05, seed = params$seed)",
    "donor_meta <- data.frame(sample = sprintf(\"acc%04d\", seq_len(params$n_donors)),",
    paste0("  condition = rep(c(", q(cond[1L]), ", ", q(cond[2L]), "), length.out = params$n_donors), batch = \"batch1\")"),
    "",
    "cat(\"2. Simulating single cells with known cell types and condition effects...\\n\")",
    "sim <- simulate_scrnaseq_from_gwas(",
    "  donors$vcf, file.path(out_dir, \"sc\"), sample_metadata = donor_meta, n_genes = params$n_genes,",
    paste0("  cells_per_donor = params$cells_per_donor, cell_type_proportions = ", type_code, ","),
    "  markers_per_cell_type = 25L, marker_effect_mean = 2, marker_effect_sd = 0.3,",
    "  n_cis_eqtl = 20L, n_trans_eqtl = 5L, condition_effect_fraction = 0.1, write_sce = FALSE,",
    "  backend = \"native\", seed = params$seed",
    ")",
    "obj <- simitall:::.simitall_read_scrna_object(file.path(out_dir, \"sc\"))",
    "",
    "cat(\"3. Clustering cells and scoring them against the true cell types...\\n\")",
    paste0("clusters <- cluster_scrnaseq(obj, n_pcs = 10L, n_hvg = round(", .simitall_agent_choice(spec, "hvg_fraction"),
           " * params$n_genes), seed = params$seed)"),
    "write.table(clusters$clusters, file.path(out_dir, \"clusters.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "cluster_metrics <- score_scrnaseq_clusters(clusters, obj$cell_metadata)",
    "write.table(cluster_metrics, file.path(out_dir, \"clustering_truth_recovery_metrics.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "",
    "cat(\"4. Testing condition differences and scoring them against the true condition genes...\\n\")",
    paste0("de <- test_scrnaseq_conditions(obj, method = ", q(.simitall_agent_choice(spec, "sc_de_method")), ", fdr = 0.05)"),
    "write.table(de, file.path(out_dir, \"condition_de.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "de_metrics <- score_scrnaseq_conditions(de, sim$condition_truth)",
    "write.table(de_metrics, file.path(out_dir, \"condition_de_truth_recovery_metrics.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "",
    "true_type <- obj$cell_metadata$cell_type[match(rownames(clusters$pcs), obj$cell_metadata$cell_id)]",
    "palette <- c(\"#2F6B5A\", \"#D9984A\", \"#7A9CC6\", \"#B85252\", \"#8E6BB8\", \"#4FA3A5\", \"#C9A227\", \"#6B6B6B\", \"#D46A9F\", \"#3E5C8A\")",
    "png(file.path(figure_dir, \"cells_true_types_vs_clusters.png\"), width = 1800, height = 850, res = 170)",
    "par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))",
    "plot(clusters$pcs[, 1:2], pch = 16, cex = 0.4, col = palette[as.integer(factor(true_type))], main = \"True cell types (simulation truth)\", xlab = \"PC1\", ylab = \"PC2\")",
    "legend(\"topright\", legend = levels(factor(true_type)), col = palette, pch = 16, cex = 0.7, bty = \"n\")",
    "cluster_palette <- c(\"#1B1B1B\", \"#E07B39\", \"#3A86A8\", \"#9C4F96\", \"#7CA33B\", \"#C23B3B\", \"#5E5E9E\", \"#B59A2F\", \"#2E8B57\", \"#8B5A2B\")",
    "plot(clusters$pcs[, 1:2], pch = 16, cex = 0.4, col = cluster_palette[as.integer(factor(clusters$clusters$cluster))],",
    "  main = sprintf(\"Clusters found (ARI = %.2f)\", cluster_metrics$adjusted_rand_index), xlab = \"PC1\", ylab = \"PC2\")",
    "dev.off()",
    "truth <- read.delim(sim$condition_truth)",
    "is_true <- de$gene_id %in% truth$gene_id[abs(truth$effect_log2) >= 0.25]",
    "png(file.path(figure_dir, \"condition_de_volcano.png\"), width = 1400, height = 1000, res = 170)",
    "plot(de$log2_fold_change, -log10(pmax(de$p_value, 1e-300)), pch = 19, cex = 0.6, col = ifelse(is_true, \"#B85252\", \"#9AA7A0\"),",
    paste0("  xlab = \"log2 fold change\", ylab = \"-log10 p\", main = ", q(paste0(cond[2L], " vs ", cond[1L], ": red = truly responsive genes (simulation truth)")), ", cex.main = 0.85)"),
    "dev.off()",
    "print(cluster_metrics); print(de_metrics)"
  )
  paste("Spec-compiled single-cell RNA-seq recipe:", "", "```r", paste(code, collapse = "\n"), "```", sep = "\n")
}

.simitall_agent_single_cell_plan <- function(spec) {
  was_default <- function(field) if (field %in% spec$defaulted) " (default)" else ""
  types <- spec$cell_types
  method <- .simitall_agent_choice(spec, "sc_de_method")
  design <- c(
    paste0("- Assay: single-cell RNA-seq, ", spec$tissue, " tissue"),
    paste0("- Cell types: ", paste0(names(types), " (", round(100 * types), "%)", collapse = ", "),
           "; 25 marker genes per type at about 4-fold (2 log2 units)"),
    paste0("- Donors: ", spec$n_donors, was_default("n_donors"), ", split between ", spec$conditions[1L], " and ", spec$conditions[2L]),
    paste0("- Cells per donor: ", spec$cells_per_donor, was_default("cells_per_donor"), "; genes: ", spec$n_genes, was_default("n_genes")),
    paste0("- Clustering: top ", round(100 * .simitall_agent_choice(spec, "hvg_fraction")), "% most variable genes, 10 PCs, k chosen by silhouette"),
    paste0("- Condition test: ", method),
    paste0("- Seed: ", spec$seed, was_default("seed")),
    paste0("- Parsed by: ", spec$extractor),
    .simitall_agent_normalised_line(spec)
  )
  stages <- c(
    "1. Simulate donor genotypes and assign donors to the two conditions.",
    "2. Simulate cells with `simulate_scrnaseq_from_gwas()`: cell-type marker programs, condition effects on a subset of genes, donor and latent variation, doublets, ambient RNA, and dropout.",
    "3. Cluster cells (variable genes, PCA, k-means) and score the clusters against the true cell types (adjusted Rand index).",
    paste0("4. Test ", spec$conditions[2L], " vs ", spec$conditions[1L], " with `", method, "` and score the calls against the true condition genes.")
  )
  inputs <- if (length(spec$errors)) paste0("- Fix before running: ", spec$errors) else "- None: everything needed is simulated."
  outputs <- c(
    "- Cell metadata with true cell types, marker and condition truth, and the count matrix.",
    "- Clusters with clustering_truth_recovery_metrics (adjusted Rand index, NMI) and a true-type vs cluster figure.",
    "- Condition DE results with condition_de_truth_recovery_metrics and a volcano plot marking truly responsive genes."
  )
  limitations <- c(
    paste0("- ", spec$resolution$note %||% "Data are simulated."),
    if (identical(method, "cell_wilcoxon")) "- Per-cell tests treat every cell as a replicate, ignoring donor variation; expect inflated false positives (compare with a pseudobulk method)."
    else "- Pseudobulk tests use donors as replicates, so power depends on the number of donors, not the number of cells.",
    "- Clustering is a transparent baseline (k-means on PCs), not Leiden; rare cell types can merge with neighbours.",
    "- Recovering the programmed truth shows the analysis works under these simulation assumptions. It is a benchmark, not a biological discovery."
  )
  paste(
    "simitall spec-based plan (the plan and code below come from the same validated spec).",
    paste0("Question: ", spec$question),
    paste0("Interpreted design:\n", paste(design, collapse = "\n")),
    paste0("Proposed stages:\n", paste(stages, collapse = "\n")),
    paste0("Required inputs:\n", paste(inputs, collapse = "\n")),
    paste0("Truth and QC outputs:\n", paste(outputs, collapse = "\n")),
    paste0("Limitations:\n", paste(limitations, collapse = "\n")),
    sep = "\n\n"
  )
}

# ---------------------------------------------------------------------------
# Breeding populations (F2, RIL, NIL/backcross, DH, MAGIC): simulate lines
# from founders and compare their genome composition with Mendelian theory.
# ---------------------------------------------------------------------------

.simitall_agent_extract_breeding <- function(question, organism, data_source, input_files) {
  text <- tolower(question)
  # Requests that analyse a trait belong to the GWAS / population-study path.
  if (grepl("\\bgwas\\b|association|\\bqtl\\b|genomic.?selection|genomic prediction|rna.?seq|expression|eqtl|phenotype|trait|\\bnam\\b|nested association|(?:find|map|which)\\w* .*genes?", text, perl = TRUE)) return(NULL)
  is <- function(pattern) grepl(pattern, text, perl = TRUE)
  scheme <- if (is("\\bmagic\\b|multi[ -]?parent|[0-9]+[ -]way cross")) "MAGIC"
    else if (is("backcross|\\bbc[0-9]|near[ -]?isogenic|\\bnils?\\b|introgress")) "NIL"
    else if (is("doubled[ -]?haploid|\\bdhs?\\b")) "DH"
    else if (is("\\brils?\\b|recombinant inbred|single[ -]seed descent|\\bssd\\b|sib[ -]mating")) "RIL"
    else if (is("\\bf2\\b")) "F2"
    else if (is("biparental|bi-parental|cross(?:ing)? two|breeding population|inbred lines?|selfing")) "GENERIC"
    else return(NULL)
  number <- function(pattern) .simitall_agent_first_number(text, pattern)
  bc <- number("\\bbc([0-9]+)") %||% number("([0-9]+)\\s*(?:rounds? of |generations? of |cycles? of )?backcross")
  self <- number("\\bbc[0-9]+s([0-9]+)\\b") %||%
    number("([0-9]+)\\s*(?:rounds? of |generations? of |cycles? of )?(?:selfing|selfings|self(?:ing)? generations|generations of selfing)")
  filial <- number("\\bf([0-9]+)\\b")
  if (identical(scheme, "GENERIC")) scheme <- if (!is.null(bc)) "NIL" else if (!is.null(self) && self > 1) "RIL" else "F2"
  if (identical(scheme, "RIL") && is.null(self) && !is.null(filial) && filial >= 2) self <- filial - 1
  founders <- number("([0-9]+)[ -]way") %||% number("([0-9]+)\\s*(?:founders|parents)")
  lines <- number("([0-9][0-9,]*)\\s*(?:lines|plants|offspring|individuals|progeny|rils|nils|dhs|families)\\b")
  seed <- number("seed\\s*(?:=|of|:)?\\s*([0-9]+)")
  defaults <- list(RIL = c(0, 6), NIL = c(3, 3), MAGIC = c(0, 3), F2 = c(0, 0), DH = c(0, 0))[[scheme]]
  list(
    spec_version = "1", extractor = "keyword_v1", question = question,
    task = "breeding", organism = organism,
    population = list(design = scheme, panel_id = NA_character_),
    other_analyses = character(),
    scheme = scheme,
    backcrosses = if (identical(scheme, "NIL")) bc %||% defaults[1L] else 0,
    self_generations = if (scheme %in% c("RIL", "NIL", "MAGIC")) self %||% defaults[2L] else 0,
    n_founders = if (identical(scheme, "MAGIC")) founders %||% 8L else 2L,
    ril_mating = if (is("sib[ -]mating|\\bsib\\b")) "SIB" else "SSD",
    n_lines = lines %||% 150L,
    seed = seed %||% 2026,
    data_source = data_source,
    input_files = as.character(input_files %||% character()),
    defaulted = c(if (identical(scheme, "NIL") && is.null(bc)) "backcrosses",
                  if (scheme %in% c("RIL", "NIL", "MAGIC") && is.null(self)) "self_generations",
                  if (identical(scheme, "MAGIC") && is.null(founders)) "n_founders",
                  if (is.null(lines)) "n_lines", if (is.null(seed)) "seed")
  )
}

.simitall_agent_validate_breeding <- function(spec) {
  errors <- character()
  add <- function(...) errors <<- c(errors, paste0(...))
  is_whole <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)
  bc <- .simitall_agent_choice(spec, "backcrosses"); self <- .simitall_agent_choice(spec, "selfing")
  if (!is_whole(spec$n_lines) || spec$n_lines < 10 || spec$n_lines > 1000) add("n_lines: must be a whole number in [10, 1000]; got ", format(spec$n_lines))
  if (!is_whole(bc) || bc < 0 || bc > 8) add("backcrosses: must be a whole number in [0, 8]; got ", format(bc))
  if (identical(spec$scheme, "NIL") && is_whole(bc) && bc < 1) add("backcrosses: a backcross (NIL) population needs at least one backcross")
  if (!is_whole(self) || self < 0 || self > 10) add("self_generations: must be a whole number in [0, 10]; got ", format(self))
  if (identical(spec$scheme, "MAGIC") && (!is_whole(spec$n_founders) || spec$n_founders < 3 || spec$n_founders > 8)) {
    add("n_founders: MAGIC needs 3 to 8 founders; got ", format(spec$n_founders))
  }
  if (!is_whole(spec$seed) || spec$seed < 1) add("seed: must be a positive whole number")
  errors
}

.simitall_agent_resolve_breeding <- function(spec) {
  list(status = "ready", mode = "simulate", files = character(),
       note = "Founder haplotypes are simulated on 5 chromosomes with a 100 cM genetic map per chromosome; no data are needed.")
}

.simitall_agent_breeding_defaults <- function(spec) {
  list(backcrosses = spec$backcrosses, selfing = spec$self_generations)
}

.simitall_agent_breeding_label <- function(spec) {
  bc <- .simitall_agent_choice(spec, "backcrosses"); self <- .simitall_agent_choice(spec, "selfing")
  switch(spec$scheme,
    F2 = "F2: one F1 cross between two parents, then one selfing generation",
    RIL = paste0("RIL: F1 between two parents, then ", self, " generations of ", if (identical(spec$ril_mating, "SIB")) "sib mating" else "single-seed-descent selfing", " (F", self + 1, ")"),
    NIL = paste0("Backcross/NIL: ", bc, " backcross", if (bc == 1) "" else "es", " to parent 1 (the recurrent parent), then ", self, " selfing generation", if (self == 1) "" else "s", " (BC", bc, "S", self, ")"),
    DH = "Doubled haploids from F1 gametes between two parents",
    MAGIC = paste0("MAGIC: ", spec$n_founders, "-founder funnel per line, then ", self, " selfing generations"))
}

.simitall_agent_compile_breeding <- function(spec) {
  q <- function(x) encodeString(x, quote = "\"")
  bc <- .simitall_agent_choice(spec, "backcrosses"); self <- .simitall_agent_choice(spec, "selfing")
  magic <- identical(spec$scheme, "MAGIC")
  code <- c(
    "if (!\"simitall\" %in% loadedNamespaces()) library(simitall)",
    "params <- list(",
    paste0("  seed = ", as.integer(spec$seed), ","),
    paste0("  n_lines = ", as.integer(spec$n_lines), "L,"),
    paste0("  backcrosses = ", as.integer(bc), "L,"),
    paste0("  self_generations = ", as.integer(self), "L,"),
    paste0("  n_founders = ", as.integer(spec$n_founders), "L"),
    ")",
    paste0("out_dir <- ", q(file.path("results", "agent_breeding", tolower(spec$scheme)))),
    "figure_dir <- file.path(out_dir, \"figures\")",
    "dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)",
    "",
    "cat(\"1. Simulating founder haplotypes and a genetic map...\\n\")",
    "panel <- file.path(out_dir, \"founders.fa\")",
    "generate_random_haplotype_panel(out_fa = panel, n_haplotypes = max(2L, params$n_founders), length = 2000L,",
    "  n_chromosomes = 5L, snp_rate = 0.02, indel_rate = 0, seed = params$seed)",
    "map <- file.path(out_dir, \"genetic_map.tsv\")",
    "write.table(data.frame(chromosome = rep(paste0(\"chr\", 1:5), each = 2), pos_bp = rep(c(1L, 2000L), 5), cM = rep(c(0, 100), 5)),",
    "  map, sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "",
    paste0("cat(\"2. Simulating ", spec$n_lines, " ", spec$scheme, " lines...\\n\")"),
    paste0("prefix <- file.path(out_dir, ", q(tolower(spec$scheme)), ")"),
    "simulate_breeding(",
    paste0("  haplotype_fa = panel, out_prefix = prefix, scheme = ", q(spec$scheme), ", n_offspring = params$n_lines,"),
    if (magic) "  founders = paste0(\"hap\", seq_len(params$n_founders), collapse = \",\")," else "  parents = \"hap1,hap2\",",
    paste0("  self_generations = params$self_generations, backcross_generations = params$backcrosses, ril_mating = ", q(spec$ril_mating), ","),
    "  recomb_map_in = map, vcf_out = paste0(prefix, \".vcf\"), seed = params$seed",
    ")",
    "",
    "cat(\"3. Measuring genome composition and comparing it with Mendelian theory...\\n\")",
    "summary <- summarize_breeding_population(",
    paste0("  paste0(prefix, \".ancestry.tsv\"), ", q(spec$scheme), ", parent1 = \"hap1\","),
    paste0("  backcrosses = ", if (identical(spec$scheme, "NIL")) "params$backcrosses" else "0L",
           ", self_generations = ", if (spec$scheme %in% c("RIL", "NIL", "MAGIC")) "params$self_generations" else "0L",
           ", n_founders = params$n_founders, ril_mating = ", q(spec$ril_mating)),
    ")",
    "write.table(summary$lines, file.path(out_dir, \"line_genome_composition.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "write.table(summary$metrics, file.path(out_dir, \"breeding_truth_recovery_metrics.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "plot_breeding_population(paste0(prefix, \".ancestry.tsv\"), summary, figure_dir)",
    "print(summary$metrics)"
  )
  paste("Spec-compiled breeding-population recipe:", "", "```r", paste(code, collapse = "\n"), "```", sep = "\n")
}

.simitall_agent_breeding_plan <- function(spec) {
  was_default <- function(field) if (field %in% spec$defaulted) " (default)" else ""
  bc <- .simitall_agent_choice(spec, "backcrosses"); self <- .simitall_agent_choice(spec, "selfing")
  expected <- breeding_expectations(spec$scheme, bc, self, spec$n_founders, spec$ril_mating)
  pct <- function(x) if (is.na(x)) "not computed" else paste0(format(round(100 * x, 1), nsmall = 1), "%")
  design <- c(
    paste0("- Scheme: ", .simitall_agent_breeding_label(spec),
           if (any(c("backcrosses", "self_generations") %in% spec$defaulted)) " (generation counts partly default)" else ""),
    paste0("- Lines: ", spec$n_lines, was_default("n_lines"),
           if (identical(spec$scheme, "MAGIC")) paste0("; founders: ", spec$n_founders, was_default("n_founders")) else ""),
    paste0("- Expected ", if (identical(spec$scheme, "NIL")) "recurrent-parent" else if (identical(spec$scheme, "MAGIC")) "per-founder" else "parent-1",
           " genome: ", pct(expected$parent1_genome), "; expected heterozygosity: ", pct(expected$heterozygosity)),
    "- Genome: 5 simulated chromosomes, 100 cM each",
    paste0("- Seed: ", spec$seed, was_default("seed")),
    paste0("- Parsed by: ", spec$extractor),
    .simitall_agent_normalised_line(spec)
  )
  stages <- c(
    "1. Simulate founder haplotypes and a genetic map (100 cM per chromosome).",
    paste0("2. Simulate ", spec$n_lines, " lines with `simulate_breeding()` (", spec$scheme, "), recording which founder every chromosome segment came from."),
    "3. Measure each line's genome composition and heterozygosity from those ancestry tracts.",
    "4. Compare the population averages with Mendelian expectations and paint the chromosomes of 20 lines."
  )
  outputs <- c(
    "- Line genotypes (VCF), ancestry tracts, and crossover breakpoints.",
    "- line_genome_composition and breeding_truth_recovery_metrics (observed vs expected, with a pass/fail truth check).",
    "- Chromosome painting of 20 lines and genome-composition histograms with the expectations marked."
  )
  limitations <- c(
    paste0("- ", spec$resolution$note %||% "Data are simulated."),
    if (identical(spec$ril_mating, "SIB")) "- Heterozygosity under sib mating has no simple closed form here, so only genome composition is checked against theory.",
    "- Expectations are population averages; individual lines vary around them because of linkage (see the histograms).",
    "- No selection is applied: lines are a random sample of the cross, as in a mapping population."
  )
  paste(
    "simitall spec-based plan (the plan and code below come from the same validated spec).",
    paste0("Question: ", spec$question),
    paste0("Interpreted design:\n", paste(design, collapse = "\n")),
    paste0("Proposed stages:\n", paste(stages, collapse = "\n")),
    paste0("Required inputs:\n", if (length(spec$errors)) paste(paste0("- Fix before running: ", spec$errors), collapse = "\n") else "- None: everything needed is simulated."),
    paste0("Truth and QC outputs:\n", paste(outputs, collapse = "\n")),
    paste0("Limitations:\n", paste(limitations, collapse = "\n")),
    sep = "\n\n"
  )
}

# ---------------------------------------------------------------------------
# Stage 2: validation. Returns a character vector; empty means valid.
# ---------------------------------------------------------------------------

.simitall_agent_validate_spec <- function(spec) {
  errors <- character()
  add <- function(...) errors <<- c(errors, paste0(...))
  is_whole <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)

  if (identical(spec$task, "family_risk")) return(.simitall_agent_validate_family_risk(spec))
  if (identical(spec$task, "peak_study")) return(.simitall_agent_validate_peak_study(spec))
  if (identical(spec$task, "single_cell")) return(.simitall_agent_validate_single_cell(spec))
  if (identical(spec$task, "breeding")) return(.simitall_agent_validate_breeding(spec))
  if (!spec$task %in% c("gwas", "population_study")) {
    add("task: only 'gwas', 'population_study', and 'family_risk' are supported by the spec compiler")
    return(errors)
  }
  needs_trait <- identical(spec$task, "gwas") || any(c("gwas", "genomic_selection") %in% spec$analyses)
  design <- spec$population$design
  if (is.null(design) || is.na(design)) {
    add("population.design: name a design such as 'NAM' or a diversity panel (e.g. 'Arabidopsis 1001 accessions')")
  } else if (!design %in% c("NAM", "diversity_panel")) {
    add("population.design: '", design, "' is not one of NAM, diversity_panel")
  }
  if (identical(design, "NAM") && !identical(spec$organism, "maize")) {
    add("organism: NAM GWAS is only implemented for maize; got '", spec$organism, "'")
  }
  if (identical(design, "diversity_panel") && identical(spec$organism, "unspecified")) {
    add("organism: a diversity-panel GWAS needs a species")
  }
  if (needs_trait && !identical(spec$trait$type, "quantitative")) {
    add("trait.type: '", spec$trait$type, "' traits are not yet supported by the spec compiler (quantitative only)")
  }
  if ("rnaseq" %in% spec$analyses) {
    g <- spec$rnaseq$n_genes
    if (!is_whole(g) || g < 10 || g > 2000) add("rnaseq.n_genes: must be a whole number in [10, 2000]; got ", format(g))
  }
  h2 <- spec$trait$heritability
  if (!is.numeric(h2) || length(h2) != 1L || !is.finite(h2) || h2 <= 0 || h2 >= 1) {
    add("trait.heritability: must be strictly between 0 and 1; got ", format(h2))
  }
  if (!is.null(spec$n_individuals)) {
    n <- spec$n_individuals
    if (!is_whole(n) || n < 10 || n > 5000) {
      add("n_individuals: must be a whole number in [10, 5000]; got ", format(n))
    }
  }
  if (!is.null(spec$trait$n_qtl)) {
    if (identical(design, "NAM") && identical(spec$task, "gwas")) {
      # The versioned NAM runner fixes one focal QTL plus six background loci.
      add("trait.n_qtl: not configurable for the NAM runner (fixed 1 focal + 6 background QTL)")
    } else if (!is_whole(spec$trait$n_qtl) || spec$trait$n_qtl < 1 || spec$trait$n_qtl > 50) {
      add("trait.n_qtl: must be a whole number in [1, 50]; got ", format(spec$trait$n_qtl))
    }
  }
  if (!is_whole(spec$seed) || spec$seed < 1) add("seed: must be a positive whole number")
  if (!spec$data_source %in% c("auto", "synthetic", "provided", "package_demo")) {
    add("data_source: '", spec$data_source, "' is not one of auto, synthetic, provided, package_demo")
  }
  missing <- spec$input_files[nzchar(spec$input_files) & !file.exists(path.expand(spec$input_files))]
  if (length(missing)) add("input_files: these files do not exist: ", paste(missing, collapse = ", "))
  errors
}

# ---------------------------------------------------------------------------
# Stage 3: input resolution. Decide where the panel comes from.
#
#   provided      -> the user's files, required
#   synthetic     -> bundled panel if simitall ships one, otherwise simulate
#   package_demo  -> bundled panel only
#   auto          -> bundled panel if available, otherwise ask
# ---------------------------------------------------------------------------

.simitall_agent_resolve_inputs <- function(spec) {
  if (identical(spec$task, "family_risk")) return(.simitall_agent_resolve_family_risk(spec))
  if (identical(spec$task, "peak_study")) return(.simitall_agent_resolve_peak_study(spec))
  if (identical(spec$task, "single_cell")) return(.simitall_agent_resolve_single_cell(spec))
  if (identical(spec$task, "breeding")) return(.simitall_agent_resolve_breeding(spec))
  catalog <- .simitall_panel_catalog()
  entry <- if (!is.na(spec$population$panel_id)) catalog[[spec$population$panel_id]] else NULL
  bundled <- !is.null(entry) && isTRUE(entry$bundled)
  design <- spec$population$design
  files <- spec$input_files[nzchar(spec$input_files)]
  label <- entry$label %||% paste(spec$organism, "diversity panel")

  ready <- function(mode, note, extra = list()) {
    c(list(status = "ready", mode = mode, note = note, files = files), extra)
  }
  ask <- function(action) list(status = "needs_input", mode = NA_character_, action = action, files = files)

  missing_text <- if (!is.null(entry) && nzchar(entry$missing_note)) {
    entry$missing_note
  } else {
    paste0("simitall does not bundle a ", label, ".")
  }

  if (identical(spec$data_source, "provided")) {
    if (identical(design, "NAM")) {
      panel <- files[grepl("\\.(fa|fasta|fna)(\\.gz)?$", files, ignore.case = TRUE)][1L]
      map <- files[grepl("(recomb|genetic|map).*(\\.tsv|\\.csv|\\.txt)$", basename(files), ignore.case = TRUE)][1L]
      if (is.na(panel) || is.na(map)) {
        return(ask("Provide a founder haplotype FASTA and a recombination-map table (filename containing map, recomb, or genetic)."))
      }
      return(ready("provided", "Using your founder panel and recombination map.",
                   list(haplotype_fa = normalizePath(path.expand(panel)), recomb_map = normalizePath(path.expand(map)))))
    }
    vcf <- files[grepl("\\.vcf(\\.gz)?$", files, ignore.case = TRUE)][1L]
    if (is.na(vcf)) return(ask("Provide a genotype VCF (.vcf or .vcf.gz) for the diversity panel."))
    return(ready("provided", "Using your genotype VCF; the phenotype and its causal truth are simulated on top of it.",
                 list(vcf = normalizePath(path.expand(vcf)))))
  }

  if (bundled && spec$data_source %in% c("auto", "synthetic", "package_demo")) {
    return(ready("bundled", entry$bundled_note))
  }

  if (identical(spec$data_source, "package_demo")) {
    return(ask(paste(missing_text, "Choose 'Generate synthetic inputs' or provide your own files.")))
  }

  if (identical(spec$data_source, "synthetic")) {
    note <- if (identical(design, "diversity_panel")) {
      paste0(
        missing_text, " Simulating a synthetic stand-in instead: inbred accessions drawn from ",
        "structured subpopulations (Balding-Nichols allele frequencies). It mimics population ",
        "structure only; it does not reproduce real accessions, LD blocks, or allele frequencies."
      )
    } else {
      paste(missing_text, "Simulating synthetic founder haplotypes instead.")
    }
    return(ready("simulate", note))
  }

  ask(paste(missing_text, "Provide the panel files, or choose 'Generate synthetic inputs' to simulate a stand-in."))
}

# ---------------------------------------------------------------------------
# Stage 4: compilation. Produces the recipe text the agent appends.
# ---------------------------------------------------------------------------

# Defaults the compilers apply when the question does not say. Kept in one
# place so the plan text and the compiled code cannot disagree.
.simitall_agent_spec_defaults <- function(spec) {
  if (identical(spec$task, "family_risk")) return(list(n_individuals = 1000L, unit = "children per group"))
  if (identical(spec$task, "population_study") && identical(spec$population$design, "NAM")) {
    return(list(n_individuals = 100L, n_qtl = 5L, unit = "NAM lines"))
  }
  if (identical(spec$population$design, "NAM")) {
    list(n_individuals = 100L, n_qtl = NA_integer_, unit = "NAM lines")
  } else {
    panel <- .simitall_panel_catalog()[[spec$population$panel_id %||% ""]]$synthetic_defaults %||%
      list(n_subpops = 3L, fst = 0.10, n_chromosomes = 3L, chromosome_length_bp = 1000000L, n_markers = 2000L)
    c(list(n_individuals = 200L, n_qtl = 5L, unit = "accessions"), panel)
  }
}

# Plan prose built from the spec, in the same section format as the local
# answer so the Shiny panels parse it. Used whenever the spec handles a
# request, so the plan always describes the code that will actually run.
# Paper-derived method records live in inst/agent/paper_records/. Only rows a
# person has marked "reviewed" ever reach a plan; unreviewed and excluded rows
# are ignored, so an unchecked extraction can never steer the agent.
.simitall_paper_records <- function(path = NULL, status = "reviewed") {
  if (is.null(path)) {
    path <- system.file("agent", "paper_records", "method_records.tsv", package = "simitall")
    if (!nzchar(path)) path <- file.path("inst", "agent", "paper_records", "method_records.tsv")
  }
  if (!file.exists(path)) return(data.frame())
  records <- utils::read.delim(path, quote = "", colClasses = "character", check.names = FALSE,
                               na.strings = character(), comment.char = "")
  records[records$review_status %in% status, , drop = FALSE]
}

.simitall_agent_spec_domains <- function(spec) {
  if (identical(spec$task, "gwas")) return("gwas")
  if (identical(spec$task, "peak_study")) return(if (identical(spec$assay, "atacseq")) "atacseq" else character())
  if (identical(spec$task, "single_cell")) return("scrnaseq")
  if (identical(spec$task, "population_study")) {
    return(unname(c(gwas = "gwas", rnaseq = "rnaseq", genomic_selection = "genomic_selection")[spec$analyses]))
  }
  character()
}

# Reviewed guidance relevant to the analyses a spec runs, strongest first.
.simitall_agent_paper_guidance <- function(spec, records = .simitall_paper_records(), max_lines = 6L) {
  domains <- .simitall_agent_spec_domains(spec)
  if (!length(domains) || !nrow(records)) return(character())
  uses <- c("default", "guided_question", "validation_check", "design_rule")
  hits <- records[records$domain %in% domains & records$simitall_use %in% uses, , drop = FALSE]
  if (!nrow(hits)) return(character())
  hits <- hits[order(match(hits$confidence, c("high", "medium", "low"))), , drop = FALSE]
  hits <- utils::head(hits, max_lines)
  paste0("- Paper guidance (", hits$record_id, "): ", hits$guidance, " [", hits$source, ", p. ", hits$page, "]")
}

# The plan's last section is Limitations, so reviewed paper guidance is
# appended there, where the Shiny Limitations panel already shows it.
.simitall_agent_spec_plan <- function(spec) {
  plan <- .simitall_agent_spec_plan_body(spec)
  guidance <- tryCatch(.simitall_agent_paper_guidance(spec), error = function(e) character())
  if (length(guidance)) plan <- paste(plan, paste(guidance, collapse = "\n"), sep = "\n")
  plan
}

.simitall_agent_spec_plan_body <- function(spec) {
  if (identical(spec$task, "family_risk")) return(.simitall_agent_family_risk_plan(spec))
  if (identical(spec$task, "population_study")) return(.simitall_agent_population_study_plan(spec))
  if (identical(spec$task, "peak_study")) return(.simitall_agent_peak_study_plan(spec))
  if (identical(spec$task, "single_cell")) return(.simitall_agent_single_cell_plan(spec))
  if (identical(spec$task, "breeding")) return(.simitall_agent_breeding_plan(spec))
  d <- .simitall_agent_spec_defaults(spec)
  is_nam <- identical(spec$population$design, "NAM")
  was_default <- function(field) if (field %in% spec$defaulted) " (default)" else ""
  n <- as.integer(spec$n_individuals %||% d$n_individuals)
  n_qtl <- spec$trait$n_qtl %||% d$n_qtl
  h2 <- spec$trait$heritability
  # The versioned NAM runner always writes its trait as drought_response_index,
  # so the plan names that rather than the requested label.
  trait <- if (is_nam) "drought_response_index" else spec$trait$name
  mode <- spec$resolution$mode %||% NA_character_
  mode_key <- if (is.na(mode)) "" else mode
  status <- if (length(spec$errors)) "invalid" else spec$resolution$status %||% "invalid"

  population <- if (is_nam) {
    paste0("NAM (", spec$population$panel_id %||% "unnamed panel", ")",
           if ("population.design" %in% spec$defaulted) " (default: no population was named, so the bundled maize breeding population is used)" else "")
  } else {
    paste0("diversity panel (", if (is.na(spec$population$panel_id)) "unnamed" else spec$population$panel_id, ")")
  }
  source_text <- switch(
    mode_key,
    bundled = "use the panel bundled with simitall",
    simulate = "simulate a synthetic stand-in panel",
    provided = "use your supplied files",
    "not resolved yet"
  )
  design <- c(
    "- Task: quantitative-trait GWAS",
    paste0("- Organism: ", spec$organism),
    paste0("- Population: ", population),
    paste0("- Size: ", n, " ", d$unit, was_default("n_individuals")),
    paste0("- Trait: `", trait, "`, heritability ", h2, was_default("trait.heritability"),
           if (is_nam) paste0(", 1 focal + 6 background QTL (trait name and QTL layout are fixed by the NAM runner",
                              if (!identical(spec$trait$name, "trait")) paste0("; you asked about ", gsub("_", " ", spec$trait$name)) else "", ")")
           else paste0(", ", n_qtl, " causal markers", if (is.null(spec$trait$n_qtl)) " (default)" else "")),
    paste0("- Data source: ", spec$data_source, " -> ", source_text),
    paste0("- Seed: ", spec$seed, was_default("seed")),
    paste0("- Parsed by: ", spec$extractor),
    .simitall_agent_normalised_line(spec)
  )

  stages <- if (is_nam) {
    founders <- switch(mode_key, bundled = "the bundled synthetic chr10 founders",
                       simulate = "eight simulated founder haplotypes", provided = "your founder panel and map",
                       "the founder panel")
    c(
      paste0("1. Simulate ", n, " NAM lines from ", founders, " with `simulate_breeding()` (scheme NAM), keeping ancestry and breakpoint truth."),
      paste0("2. Program `", trait, "` with one focal QTL and six background loci at heritability ", h2, ", and save the QTL truth."),
      "3. Run a family-aware GWAS with `analyze_gwas()` (family fixed effect, 2 PCs, MAF >= 0.05).",
      "4. Score recovery with `benchmark_gwas()` (150 bp window) and draw Manhattan, QQ, and phenotype figures."
    )
  } else {
    first <- if (identical(mode, "provided")) {
      "1. Read your genotype VCF with the simitall genotype reader."
    } else {
      paste0("1. Simulate ", n, " inbred accessions in ", d$n_subpops, " subpopulations (Balding-Nichols, Fst ",
             d$fst, ") on ", d$n_chromosomes, " chromosomes, writing the VCF and subpopulation labels.")
    }
    c(
      first,
      paste0("2. Program `", trait, "` with ", n_qtl, " additive causal markers (MAF >= 0.10) at heritability ", h2, ", and save the causal truth."),
      "3. Run a mixed-model GWAS with `analyze_gwas()` (kinship + 3 PCs, MAF >= 0.05) to control population structure.",
      "4. Score recovery with `benchmark_gwas()` against the causal markers and draw Manhattan, QQ, and phenotype figures."
    )
  }

  inputs <- if (length(spec$errors)) {
    paste0("- Fix before running: ", spec$errors)
  } else if (identical(status, "needs_input")) {
    paste0("- ", spec$resolution$action)
  } else if (identical(mode, "provided")) {
    paste0("- ", spec$resolution$files)
  } else if (identical(mode, "bundled")) {
    "- None: the bundled panel and map ship with simitall."
  } else {
    "- None: genotypes, trait, and truth are all simulated."
  }

  outputs <- c(
    if (is_nam) "- NAM VCF, family metadata, ancestry tracts, and breakpoint truth."
    else "- Genotype VCF (simulated or supplied) and subpopulation labels for simulated panels.",
    paste0("- `", trait, "` phenotypes and causal-QTL truth tables."),
    "- GWAS results, truth-recovery metrics, Manhattan/QQ figures, and a phenotype distribution figure."
  )

  limitations <- c(
    if (!is.null(spec$resolution$note) && nzchar(spec$resolution$note)) paste0("- ", spec$resolution$note),
    if (!is_nam && !identical(mode, "provided")) "- Simulated markers are independent (no LD), so a hit is either the causal marker or a false positive; there are no linked tag SNPs.",
    "- Recovering the programmed loci shows the method works under these assumptions. It is a benchmark, not a biological discovery."
  )

  paste(
    "simitall spec-based plan (the plan and code below come from the same validated spec).",
    paste0("Question: ", spec$question),
    paste0("Interpreted design:\n", paste(design, collapse = "\n")),
    paste0("Proposed stages:\n", paste(stages, collapse = "\n")),
    paste0("Required inputs:\n", paste(inputs, collapse = "\n")),
    paste0("Truth and QC outputs:\n", paste(outputs, collapse = "\n")),
    paste0("Limitations:\n", paste(limitations, collapse = "\n")),
    sep = "\n\n"
  )
}

.simitall_agent_compile_nam_gwas <- function(spec) {
  script_path <- .simitall_agent_source_file("analysis", "paper_fig", "fig15_maize_nam_quantitative_gwas.R")
  if (is.null(script_path)) {
    return("The NAM quantitative-GWAS runner is unavailable from this installation. Run from a simitall source checkout containing `analysis/paper_fig/fig15_maize_nam_quantitative_gwas.R`.")
  }
  input_mode <- switch(spec$resolution$mode, bundled = "package_demo", simulate = "synthetic", provided = "provided")
  n_lines <- as.integer(spec$n_individuals %||% .simitall_agent_spec_defaults(spec)$n_individuals)
  args <- c(
    "  \"analysis/paper_fig/fig15_maize_nam_quantitative_gwas.R\",",
    "  \"--out_dir\", \"analysis/results/maize_nam_quantitative_gwas\",",
    paste0("  \"--n_lines\", \"", n_lines, "\","),
    paste0("  \"--heritability\", \"", format(spec$trait$heritability, nsmall = 2), "\","),
    paste0("  \"--input_mode\", \"", input_mode, "\","),
    if (identical(input_mode, "provided")) c(
      paste0("  \"--haplotype_fa\", ", encodeString(spec$resolution$haplotype_fa, quote = "\""), ","),
      paste0("  \"--recomb_map\", ", encodeString(spec$resolution$recomb_map, quote = "\""), ",")
    ),
    paste0("  \"--seed\", \"", as.integer(spec$seed), "\"")
  )
  paste(
    "Spec-compiled NAM quantitative-trait GWAS recipe:",
    "",
    "```r",
    "system2(\"Rscript\", c(",
    paste(args, collapse = "\n"),
    "))",
    "```",
    "",
    "Steps: simulate the NAM population from the founders, program one quantitative trait with known QTL, run a family-aware GWAS, and score recovery against the causal truth.",
    "",
    "Full versioned R implementation:",
    "",
    "```r",
    paste(readLines(script_path, warn = FALSE), collapse = "\n"),
    "```",
    sep = "\n"
  )
}

.simitall_agent_compile_panel_gwas <- function(spec) {
  defaults <- .simitall_agent_spec_defaults(spec)
  panel_id <- if (is.na(spec$population$panel_id)) paste0(spec$organism, "_panel") else spec$population$panel_id
  out_dir <- file.path("results", "agent_gwas", panel_id)
  n_lines <- as.integer(spec$n_individuals %||% defaults$n_individuals)
  n_qtl <- as.integer(spec$trait$n_qtl %||% defaults$n_qtl)
  trait <- spec$trait$name

  genotype_block <- if (identical(spec$resolution$mode, "provided")) {
    c(
      paste0("vcf <- ", encodeString(spec$resolution$vcf, quote = "\"")),
      "cat(\"1. Using the supplied genotype VCF...\\n\")"
    )
  } else {
    c(
      "cat(\"1. Simulating a synthetic inbred accession panel...\\n\")",
      "panel <- simitall:::.simitall_simulate_inbred_panel(",
      "  vcf_out = file.path(out_dir, \"synthetic_panel.vcf\"),",
      "  n_accessions = params$n_lines, n_markers = params$n_markers,",
      paste0("  n_chromosomes = ", defaults$n_chromosomes, "L, chromosome_length_bp = ", defaults$chromosome_length_bp, "L,"),
      paste0("  n_subpops = ", defaults$n_subpops, "L, fst = ", defaults$fst, ", seed = params$seed"),
      ")",
      "vcf <- panel$vcf"
    )
  }

  code <- c(
    "if (!\"simitall\" %in% loadedNamespaces()) library(simitall)",
    "params <- list(",
    paste0("  seed = ", as.integer(spec$seed), ","),
    paste0("  n_lines = ", n_lines, "L,"),
    paste0("  n_markers = ", as.integer(defaults$n_markers), "L,"),
    paste0("  heritability = ", spec$trait$heritability, ","),
    paste0("  n_qtl = ", n_qtl, "L"),
    ")",
    paste0("out_dir <- ", encodeString(out_dir, quote = "\"")),
    "figure_dir <- file.path(out_dir, \"figures\")",
    "dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)",
    "",
    genotype_block,
    "",
    "cat(\"2. Programming a quantitative trait with known causal markers...\\n\")",
    "geno <- simitall:::.simitall_read_genotypes(vcf)",
    "trait <- simitall:::.simitall_simulate_additive_trait(",
    "  geno$genotype, geno$variants, heritability = params$heritability,",
    paste0("  n_qtl = params$n_qtl, trait_name = \"", trait, "\", seed = params$seed + 1L"),
    ")",
    paste0("phenotype_file <- file.path(out_dir, \"", trait, "_phenotypes.tsv\")"),
    paste0("truth_file <- file.path(out_dir, \"", trait, "_qtl_truth.tsv\")"),
    "write.table(trait$phenotype, phenotype_file, sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "write.table(trait$truth, truth_file, sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "",
    "cat(\"3. Running a kinship + PC mixed-model GWAS...\\n\")",
    "gwas_prefix <- file.path(out_dir, \"gwas\")",
    "gwas <- analyze_gwas(",
    paste0("  vcf, trait$phenotype, gwas_prefix, trait = \"", trait, "\","),
    paste0("  n_pcs = ", .simitall_agent_choice(spec, "gwas_pcs"), "L, min_maf = 0.05, fdr_threshold = ", .simitall_agent_choice(spec, "gwas_fdr")),
    ")",
    "",
    "cat(\"4. Scoring recovery against the causal truth...\\n\")",
    paste0("benchmark <- benchmark_gwas(gwas$paths$results, trait$truth, gwas_prefix, fdr_threshold = ",
           .simitall_agent_choice(spec, "gwas_fdr"), ", window_bp = 0L)"),
    "plot_gwas_results(",
    "  benchmark$paths$results, file.path(figure_dir, \"gwas_overview\"),",
    "  individual_panels_dir = figure_dir",
    ")",
    paste0("plot_trait_diagnostics(trait$phenotype, \"", trait, "\", figure_dir, file_stem = \"phenotype_distribution\")"),
    "metrics <- read.delim(benchmark$paths$metrics, check.names = FALSE)",
    "write.table(metrics, file.path(out_dir, \"truth_recovery_metrics.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "print(metrics)"
  )

  paste(
    paste0("Spec-compiled diversity-panel GWAS recipe (", panel_id, "):"),
    "",
    "```r",
    paste(code, collapse = "\n"),
    "```",
    "",
    "Expected outputs: genotype VCF (simulated or supplied), phenotype and QTL-truth TSVs, GWAS results, truth-recovery metrics, Manhattan/QQ figures, and a phenotype distribution figure.",
    sep = "\n"
  )
}

# ---------------------------------------------------------------------------
# Composable population studies: the population stage always runs, then only
# the analyses the question named (GWAS, RNA-seq, genomic selection), each
# scored against its own simulation truth.
# ---------------------------------------------------------------------------

.simitall_agent_analysis_labels <- function() {
  c(gwas = "GWAS", rnaseq = "RNA-seq (differential expression + eQTL)", genomic_selection = "genomic selection")
}

.simitall_agent_compile_population_study <- function(spec) {
  d <- .simitall_agent_spec_defaults(spec)
  is_nam <- identical(spec$population$design, "NAM")
  analyses <- spec$analyses
  needs_trait <- any(c("gwas", "genomic_selection") %in% analyses)
  mode <- spec$resolution$mode
  panel_id <- if (is.na(spec$population$panel_id)) paste0(spec$organism, "_panel") else spec$population$panel_id
  out_dir <- file.path("results", "agent_study", paste(c(panel_id, analyses), collapse = "_"))
  trait <- spec$trait$name
  q <- function(x) encodeString(x, quote = "\"")
  step <- 0L
  say <- function(text) { step <<- step + 1L; paste0("cat(\"", step, ". ", text, "...\\n\")") }

  params <- c(
    paste0("  seed = ", as.integer(spec$seed)),
    paste0("  n_lines = ", as.integer(spec$n_individuals %||% d$n_individuals), "L"),
    if (needs_trait) paste0("  heritability = ", spec$trait$heritability),
    if (needs_trait) paste0("  n_qtl = ", as.integer(spec$trait$n_qtl %||% d$n_qtl), "L"),
    if ("rnaseq" %in% analyses) paste0("  n_genes = ", as.integer(spec$rnaseq$n_genes), "L")
  )
  code <- c(
    "if (!\"simitall\" %in% loadedNamespaces()) library(simitall)",
    "params <- list(",
    paste(params, collapse = ",\n"),
    ")",
    paste0("out_dir <- ", q(out_dir)),
    "figure_dir <- file.path(out_dir, \"figures\")",
    "dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)",
    ""
  )

  # Population stage (always).
  if (is_nam) {
    founders <- switch(
      mode,
      bundled = c(
        "panel <- system.file(\"extdata\", \"panels\", \"demo_maize_nam_chr10.fa\", package = \"simitall\")",
        "map <- system.file(\"extdata\", \"maps\", \"demo_maize_nam_chr10_map.tsv\", package = \"simitall\")",
        "stopifnot(nzchar(panel), nzchar(map))"
      ),
      provided = c(
        paste0("panel <- ", q(spec$resolution$haplotype_fa)),
        paste0("map <- ", q(spec$resolution$recomb_map))
      ),
      c(
        "panel <- file.path(out_dir, \"synthetic_founders.fa\")",
        "generate_random_haplotype_panel(out_fa = panel, n_haplotypes = 8L, length = 4000L, n_chromosomes = 1L,",
        "  snp_rate = 0.02, indel_rate = 0, seed = params$seed + 101L)",
        "panel_lines <- readLines(panel)",
        "panel_lines[grepl(\"^>\", panel_lines)] <- paste0(panel_lines[grepl(\"^>\", panel_lines)], \"|chr10\")",
        "writeLines(panel_lines, panel)",
        "map <- file.path(out_dir, \"synthetic_map.tsv\")",
        "write.table(data.frame(chromosome = \"chr10\", pos_bp = c(1L, 1000L, 1800L, 2600L, 4000L), cM = c(0, 1, 1.1, 1.2, 2.6)),",
        "  map, sep = \"\\t\", row.names = FALSE, quote = FALSE)"
      )
    )
    code <- c(code,
      say("Simulating the NAM population"),
      founders,
      "prefix <- file.path(out_dir, \"nam\")",
      "vcf <- paste0(prefix, \".vcf\")",
      "simulate_breeding(",
      "  haplotype_fa = panel, out_prefix = prefix, scheme = \"NAM\", n_offspring = params$n_lines,",
      "  recomb_map_in = map, vcf_out = vcf, ancestry_out = paste0(prefix, \".ancestry.tsv\"),",
      "  breakpoints_out = paste0(prefix, \".breakpoints.tsv\"), seed = params$seed",
      ")",
      "metadata <- read.delim(paste0(prefix, \".meta.tsv\"), check.names = FALSE, stringsAsFactors = FALSE)"
    )
  } else {
    pd <- .simitall_panel_catalog()[[spec$population$panel_id %||% ""]]$synthetic_defaults %||%
      list(n_subpops = 3L, fst = 0.10, n_chromosomes = 3L, chromosome_length_bp = 1000000L, n_markers = 2000L)
    population <- if (identical(mode, "provided")) {
      paste0("vcf <- ", q(spec$resolution$vcf))
    } else {
      c(
        "panel <- simitall:::.simitall_simulate_inbred_panel(",
        paste0("  vcf_out = file.path(out_dir, \"synthetic_panel.vcf\"), n_accessions = params$n_lines, n_markers = ", pd$n_markers, "L,"),
        paste0("  n_chromosomes = ", pd$n_chromosomes, "L, chromosome_length_bp = ", pd$chromosome_length_bp, "L, n_subpops = ",
               pd$n_subpops, "L, fst = ", pd$fst, ", seed = params$seed"),
        ")",
        "vcf <- panel$vcf"
      )
    }
    code <- c(code, say("Preparing the accession panel"), population, "metadata <- NULL")
  }
  code <- c(code,
    "geno <- simitall:::.simitall_read_genotypes(vcf)",
    "samples <- colnames(geno$genotype)",
    "metrics <- list()",
    ""
  )

  if (needs_trait) {
    code <- c(code,
      say(paste0("Programming the ", gsub("_", " ", trait), " trait with known causal markers")),
      "trait <- simitall:::.simitall_simulate_additive_trait(",
      "  geno$genotype, geno$variants, heritability = params$heritability, n_qtl = params$n_qtl,",
      paste0("  trait_name = ", q(trait), ", seed = params$seed + 1L"),
      ")",
      "phenotype <- trait$phenotype",
      "if (!is.null(metadata)) phenotype$family <- metadata$family[match(phenotype$sample, metadata$sample)]",
      paste0("write.table(phenotype, file.path(out_dir, ", q(paste0(trait, "_phenotypes.tsv")), "), sep = \"\\t\", row.names = FALSE, quote = FALSE)"),
      paste0("write.table(trait$truth, file.path(out_dir, ", q(paste0(trait, "_qtl_truth.tsv")), "), sep = \"\\t\", row.names = FALSE, quote = FALSE)"),
      paste0("plot_trait_diagnostics(phenotype, ", q(trait), ", figure_dir, file_stem = \"phenotype_distribution\")"),
      ""
    )
  }

  if ("gwas" %in% analyses) {
    code <- c(code,
      say("Running GWAS and scoring it against the causal markers"),
      "gwas <- analyze_gwas(",
      paste0("  vcf, phenotype, file.path(out_dir, \"gwas\"), trait = ", q(trait), ","),
      paste0(if (is_nam) "  fixed_effects = \"family\", " else "  ", "n_pcs = ", .simitall_agent_choice(spec, "gwas_pcs"),
             "L, min_maf = 0.05, fdr_threshold = ", .simitall_agent_choice(spec, "gwas_fdr")),
      ")",
      paste0("gwas_benchmark <- benchmark_gwas(gwas$paths$results, trait$truth, file.path(out_dir, \"gwas\"), fdr_threshold = ",
             .simitall_agent_choice(spec, "gwas_fdr"), ", window_bp = ", if (is_nam) "150L" else "0L", ")"),
      "plot_gwas_results(gwas_benchmark$paths$results, file.path(figure_dir, \"gwas_overview\"), individual_panels_dir = figure_dir)",
      "gwas_metrics <- read.delim(gwas_benchmark$paths$metrics, check.names = FALSE)",
      "write.table(gwas_metrics, file.path(out_dir, \"gwas_truth_recovery_metrics.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
      ""
    )
  }

  if ("rnaseq" %in% analyses) {
    cond <- spec$rnaseq$conditions
    gff <- spec$input_files[grepl("\\.gff3?$", spec$input_files, ignore.case = TRUE)][1L]
    code <- c(code,
      say(paste0("Simulating RNA-seq (", cond[1L], " vs ", cond[2L], ") from the same genotypes")),
      "rna_meta <- data.frame(sample = samples, stringsAsFactors = FALSE)",
      paste0("rna_meta$condition <- rep(c(", q(cond[1L]), ", ", q(cond[2L]), "), length.out = length(samples))"),
      "rna_meta$batch <- rep(c(\"batch1\", \"batch2\"), each = 2L, length.out = length(samples))",
      "rna_meta_file <- file.path(out_dir, \"rnaseq_design.tsv\")",
      "write.table(rna_meta, rna_meta_file, sep = \"\\t\", row.names = FALSE, quote = FALSE)",
      "rna_prefix <- file.path(out_dir, \"rnaseq\")",
      "simulate_rnaseq_from_gwas(",
      "  vcf, rna_prefix, sample_metadata = rna_meta_file,",
      if (!is.na(gff)) paste0("  annotation_gff3 = ", q(normalizePath(path.expand(gff))), ",") else "  n_genes = params$n_genes,",
      paste0("  condition_levels = c(", q(cond[1L]), ", ", q(cond[2L]), "),"),
      "  n_cis_eqtl = max(2L, round(params$n_genes / 5)), n_trans_eqtl = max(1L, round(params$n_genes / 15)),",
      "  condition_effect_fraction = 0.3, gxe_fraction = 0.25, seed = params$seed + 2L",
      ")",
      "",
      say("Finding responsive genes (differential expression) and their regulators (eQTL)"),
      "rna_diagnostics <- plot_rnaseq_diagnostics(",
      "  paste0(rna_prefix, \".counts.tsv\"), paste0(rna_prefix, \".sample_metadata.tsv\"),",
      paste0("  out_dir = figure_dir, fdr_threshold = ", .simitall_agent_choice(spec, "rnaseq_fdr")),
      ")",
      "de <- read.delim(rna_diagnostics$paths$deg, check.names = FALSE)",
      "de_truth <- read.delim(paste0(rna_prefix, \".condition_truth.tsv\"), check.names = FALSE)",
      "de_metrics <- simitall:::.simitall_score_gene_recovery(",
      "  called = de$gene_id[de$significant], truth = unique(de_truth$gene_id), tested = de$gene_id,",
      paste0("  analysis = ", q(paste0("genes responding to ", cond[2L], " (differential expression)"))),
      ")",
      "write.table(de_metrics, file.path(out_dir, \"rnaseq_de_truth_recovery_metrics.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
      "benchmark_eqtl(",
      "  vcf, paste0(rna_prefix, \".expression.tsv\"), file.path(out_dir, \"eqtl\"),",
      "  sample_metadata = paste0(rna_prefix, \".sample_metadata.tsv\"),",
      "  truth_file = paste0(rna_prefix, \".eqtl_truth.tsv\")",
      ")",
      ""
    )
  }

  if ("genomic_selection" %in% analyses) {
    code <- c(code,
      say("Predicting breeding values and selecting parents"),
      "selection <- run_genomic_selection(",
      paste0("  vcf, phenotype, file.path(out_dir, \"selection\"), trait = ", q(trait), ", model = ", q(.simitall_agent_choice(spec, "gs_model")), ","),
      "  n_parents = 12L, n_crosses = 18L, diversity_penalty = 0.25,",
      if (is_nam) "  family_data = metadata[, c(\"sample\", \"family\")], seed = params$seed + 3L" else "  seed = params$seed + 3L",
      ")",
      "gs <- merge(selection$predictions[, c(\"sample\", \"gebv\")], trait$phenotype[, c(\"sample\", \"genetic_value\")], by = \"sample\")",
      "picked <- gs$sample %in% selection$selected$sample",
      "# Held-out accuracy: refit on 80% of lines and score the other 20% against",
      "# their true genetic values. This is the fair number for comparing models.",
      "set.seed(params$seed + 4L)",
      "holdout <- sample(samples, max(2L, round(length(samples) * 0.2)))",
      paste0("holdout_fit <- fit_genomic_model(vcf, phenotype, trait = ", q(trait), ", model = ",
             q(.simitall_agent_choice(spec, "gs_model")), ", training_samples = setdiff(samples, holdout), seed = params$seed + 5L)"),
      "held <- merge(holdout_fit$predictions[holdout_fit$predictions$sample %in% holdout, c(\"sample\", \"gebv\")],",
      "  trait$phenotype[, c(\"sample\", \"genetic_value\")], by = \"sample\")",
      "gs_metrics <- data.frame(",
      paste0("  analysis = \"genomic selection\", model = ", q(.simitall_agent_choice(spec, "gs_model")), ", lines = nrow(gs), selected_parents = sum(picked),"),
      "  holdout_accuracy = stats::cor(held$gebv, held$genetic_value), holdout_lines = nrow(held),",
      "  in_sample_accuracy = stats::cor(gs$gebv, gs$genetic_value),",
      "  true_value_all_lines = mean(gs$genetic_value), true_value_selected = mean(gs$genetic_value[picked]),",
      "  note = \"holdout_accuracy: trained on 80% of lines, scored on the other 20%; in_sample is optimistic\"",
      ")",
      "write.table(gs_metrics, file.path(out_dir, \"genomic_selection_truth_recovery_metrics.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
      "png(file.path(figure_dir, \"genomic_selection_true_value.png\"), width = 1400, height = 900, res = 170)",
      "boxplot(gs$genetic_value, gs$genetic_value[picked], names = c(\"All lines\", \"Selected parents\"),",
      "  col = c(\"#BED8E6\", \"#E9B2A8\"), ylab = \"True genetic value (simulation truth)\",",
      "  main = \"Did genomic selection pick genuinely better parents?\")",
      "dev.off()",
      ""
    )
  }

  code <- c(code,
    "for (f in list.files(out_dir, \"truth_recovery_metrics\\\\.tsv$\", full.names = TRUE)) print(read.delim(f, check.names = FALSE))"
  )
  label <- paste(unname(.simitall_agent_analysis_labels()[analyses]), collapse = " + ")
  paste(
    paste0("Spec-compiled population study recipe (", label, "):"),
    "",
    "```r",
    paste(code, collapse = "\n"),
    "```",
    sep = "\n"
  )
}

.simitall_agent_population_study_plan <- function(spec) {
  d <- .simitall_agent_spec_defaults(spec)
  is_nam <- identical(spec$population$design, "NAM")
  analyses <- spec$analyses
  labels <- .simitall_agent_analysis_labels()
  needs_trait <- any(c("gwas", "genomic_selection") %in% analyses)
  was_default <- function(field) if (field %in% spec$defaulted) " (default)" else ""
  n <- as.integer(spec$n_individuals %||% d$n_individuals)
  mode <- spec$resolution$mode %||% NA_character_
  mode_key <- if (is.na(mode)) "" else mode
  source_text <- switch(mode_key,
    bundled = "the panel bundled with simitall", simulate = "a simulated stand-in panel",
    provided = "your supplied files", "not resolved yet")
  population <- if (is_nam) {
    paste0("NAM (", spec$population$panel_id %||% "unnamed", ")",
           if ("population.design" %in% spec$defaulted) " (default: no population was named, so the bundled maize breeding population is used)" else "")
  } else {
    paste0("diversity panel (", if (is.na(spec$population$panel_id)) "unnamed" else spec$population$panel_id, ")")
  }
  skipped <- setdiff(names(labels), analyses)
  design <- c(
    paste0("- Analyses requested: ", paste(unname(labels[analyses]), collapse = ", ")),
    if (length(skipped)) paste0("- Not requested, so not run: ", paste(unname(labels[skipped]), collapse = ", ")),
    paste0("- Organism: ", spec$organism),
    paste0("- Population: ", population),
    paste0("- Size: ", n, " ", d$unit, was_default("n_individuals")),
    if (needs_trait) paste0("- Trait: `", spec$trait$name, "`, heritability ", spec$trait$heritability,
                            was_default("trait.heritability"), ", ", spec$trait$n_qtl %||% d$n_qtl, " causal markers",
                            if (is.null(spec$trait$n_qtl)) " (default)" else ""),
    if ("rnaseq" %in% analyses) paste0("- RNA-seq: ", spec$rnaseq$n_genes, " genes, ", spec$rnaseq$conditions[1L], " vs ",
                                      spec$rnaseq$conditions[2L], ", samples alternate between the two conditions"),
    paste0("- Data source: ", spec$data_source, " -> ", source_text),
    paste0("- Seed: ", spec$seed, was_default("seed")),
    paste0("- Parsed by: ", spec$extractor),
    .simitall_agent_normalised_line(spec)
  )
  stages <- c(
    if (is_nam) paste0("Simulate ", n, " NAM lines from ", source_text, " with `simulate_breeding()`, keeping ancestry and breakpoint truth.")
    else paste0("Prepare ", n, " accessions from ", source_text, "."),
    if (needs_trait) paste0("Program `", spec$trait$name, "` with known causal markers and save the truth table."),
    if ("gwas" %in% analyses) paste0("Run GWAS with `analyze_gwas()` (", if (is_nam) "family fixed effect + 2 PCs" else "kinship + 3 PCs",
                                    ") and score recovery with `benchmark_gwas()`."),
    if ("rnaseq" %in% analyses) paste0("Simulate RNA-seq counts from the same genotypes with `simulate_rnaseq_from_gwas()` (",
                                      spec$rnaseq$conditions[1L], " vs ", spec$rnaseq$conditions[2L], ")."),
    if ("rnaseq" %in% analyses) paste0("Find genes that respond to ", spec$rnaseq$conditions[2L],
                                      " (differential expression) and variants that control expression (eQTL, `benchmark_eqtl()`); score both against the simulated truth."),
    if ("genomic_selection" %in% analyses) "Predict breeding values with GBLUP (`run_genomic_selection()`), select parents, and check them against the true genetic values."
  )
  stages <- paste0(seq_along(stages), ". ", stages)
  inputs <- if (length(spec$errors)) {
    paste0("- Fix before running: ", spec$errors)
  } else if (identical(spec$resolution$status, "needs_input")) {
    paste0("- ", spec$resolution$action)
  } else if (identical(mode, "provided")) {
    paste0("- ", spec$resolution$files)
  } else {
    "- None: everything needed is simulated or bundled."
  }
  outputs <- c(
    "- Population genotypes (VCF) and, for NAM, family, ancestry, and breakpoint truth.",
    if (needs_trait) "- Trait values and the causal-marker truth table.",
    if ("gwas" %in% analyses) "- GWAS results, Manhattan/QQ figures, and gwas_truth_recovery_metrics.",
    if ("rnaseq" %in% analyses) "- RNA-seq counts, PCA/volcano/library-size figures, rnaseq_de_truth_recovery_metrics, and eQTL metrics.",
    if ("genomic_selection" %in% analyses) "- Breeding values, selected parents, and genomic_selection_truth_recovery_metrics."
  )
  rnaseq_only_hunt <- "rnaseq" %in% analyses && !needs_trait &&
    grepl("(?:find|locate|identify|which|map)\\w* .*genes?", tolower(spec$question), perl = TRUE)
  limitations <- c(
    if (!is.null(spec$resolution$note) && nzchar(spec$resolution$note)) paste0("- ", spec$resolution$note),
    if (rnaseq_only_hunt) paste0("- RNA-seq alone finds genes that respond to ", spec$rnaseq$conditions[2L],
                                 " and the variants that regulate them. It does not show which gene causes the trait. Add GWAS to the question to link genes to a measured trait."),
    if ("genomic_selection" %in% analyses) "- Genomic-selection accuracy is reported two ways: held-out (trained on 80% of lines, scored on the other 20%) and in-sample, which is optimistic.",
    "- Recovering the programmed truth shows each method works under these simulation assumptions. It is a benchmark, not a biological discovery."
  )
  paste(
    "simitall spec-based plan (the plan and code below come from the same validated spec).",
    paste0("Question: ", spec$question),
    paste0("Interpreted design:\n", paste(design, collapse = "\n")),
    paste0("Proposed stages:\n", paste(stages, collapse = "\n")),
    paste0("Required inputs:\n", paste(inputs, collapse = "\n")),
    paste0("Truth and QC outputs:\n", paste(outputs, collapse = "\n")),
    paste0("Limitations:\n", paste(limitations, collapse = "\n")),
    sep = "\n\n"
  )
}

# One entry point. Returns the spec with $errors, $resolution, $supported and
# $recipe filled in. Never throws for an unsupported question.
.simitall_agent_spec_pipeline <- function(question, input_files = character(), data_source = "auto") {
  spec <- .simitall_agent_extract_spec(question, data_source, input_files)
  # Unsupported shapes (other tasks, no design, binary traits) go back to the
  # legacy recipe router; supported shapes with bad parameters fail loudly.
  spec$supported <- !length(spec$other_analyses) && (
    identical(spec$task, "family_risk") || identical(spec$task, "peak_study") || identical(spec$task, "single_cell") || identical(spec$task, "breeding") || (
      identical(spec$task, "gwas") &&
        !is.na(spec$population$design %||% NA_character_) &&
        identical(spec$trait$type, "quantitative")
    ) || (
      identical(spec$task, "population_study") &&
        !is.na(spec$population$design %||% NA_character_) &&
        (identical(spec$trait$type, "quantitative") || !any(c("gwas", "genomic_selection") %in% spec$analyses))
    )
  )
  spec$errors <- .simitall_agent_validate_spec(spec)
  spec$resolution <- NULL
  spec$recipe <- ""
  if (!spec$supported || length(spec$errors)) return(spec)
  .simitall_agent_compile_spec(spec)
}

# Resolve inputs and compile a validated spec into its recipe.
.simitall_agent_compile_spec <- function(spec) {
  spec$resolution <- .simitall_agent_resolve_inputs(spec)
  spec$recipe <- ""
  if (!identical(spec$resolution$status, "ready")) return(spec)
  spec$recipe <- if (identical(spec$task, "family_risk")) {
    .simitall_agent_compile_family_risk(spec)
  } else if (identical(spec$task, "peak_study")) {
    .simitall_agent_compile_peak_study(spec)
  } else if (identical(spec$task, "single_cell")) {
    .simitall_agent_compile_single_cell(spec)
  } else if (identical(spec$task, "breeding")) {
    .simitall_agent_compile_breeding(spec)
  } else if (identical(spec$task, "population_study")) {
    .simitall_agent_compile_population_study(spec)
  } else if (identical(spec$population$design, "NAM")) {
    .simitall_agent_compile_nam_gwas(spec)
  } else {
    .simitall_agent_compile_panel_gwas(spec)
  }
  spec
}

# ---------------------------------------------------------------------------
# Guided choices: method decisions the agent can ask about. Each has a
# default (quick mode needs no answers), options the code can really run,
# a reason, and the paper-record ids whose reviewed guidance supports it.
# ---------------------------------------------------------------------------

.simitall_agent_choice_defaults <- function(spec) {
  if (identical(spec$task, "peak_study")) return(.simitall_agent_peak_defaults(spec))
  if (identical(spec$task, "single_cell")) return(.simitall_agent_single_cell_defaults(spec))
  if (identical(spec$task, "breeding")) return(.simitall_agent_breeding_defaults(spec))
  is_nam <- identical(spec$population$design, "NAM")
  list(gwas_pcs = if (is_nam) 2L else 3L, gwas_fdr = 0.05, rnaseq_fdr = 0.05, gs_model = "gblup")
}

.simitall_agent_choice <- function(spec, id) {
  spec$choices[[id]] %||% .simitall_agent_choice_defaults(spec)[[id]]
}

.simitall_agent_spec_questions <- function(spec) {
  if (!isTRUE(spec$supported)) return(list())
  panel_gwas <- identical(spec$task, "gwas") && identical(spec$population$design, "diversity_panel")
  study <- identical(spec$task, "population_study")
  has <- function(a) study && a %in% spec$analyses
  d <- .simitall_agent_choice_defaults(spec)
  q <- function(id, question, values, labels, why, records) {
    current <- .simitall_agent_choice(spec, id)
    if (!as.character(current) %in% as.character(values)) {
      values <- c(current, values); labels <- c(paste0(current, " (current)"), labels)
    }
    list(id = id, question = question, values = as.character(values), labels = labels,
         selected = as.character(current), default = as.character(d[[id]]), why = why, records = records)
  }
  out <- list()
  if (identical(spec$task, "breeding")) {
    if (identical(spec$scheme, "NIL")) {
      out[[1L]] <- q(
        "backcrosses", "How many backcrosses to the recurrent parent?",
        c(1L, 2L, 3L, 4L), c("1 (BC1)", "2 (BC2)", "3 (BC3)", "4 (BC4)"),
        "Each backcross halves the remaining donor genome, but costs another season.",
        character())
    }
    if (spec$scheme %in% c("RIL", "NIL", "MAGIC")) {
      out[[length(out) + 1L]] <- q(
        "selfing", "How many selfing generations?",
        c(2L, 4L, 6L, 8L), c("2", "4", "6", "8"),
        "Each selfing generation halves the remaining heterozygosity, so lines become more inbred (fixed) but take longer to develop.",
        character())
    }
    return(out)
  }
  if (identical(spec$task, "single_cell")) {
    out[[1L]] <- q(
      "hvg_fraction", "What share of genes should be used as highly variable genes for clustering?",
      c(0.1, 0.2, 0.5, 1), c("10% of genes", "20% of genes", "50% of genes", "All genes"),
      "Keeping a minority of the most variable genes focuses clustering on cell-type markers; using every gene lets genome-wide noise blur the cell types.",
      c("scrnaseq_012", "scrnaseq_001"))
    out[[2L]] <- q(
      "sc_de_method", paste0("How should genes that respond to ", spec$conditions[2L], " be tested?"),
      c("pseudobulk_deseq2", "pseudobulk_edger", "pseudobulk_ttest", "cell_wilcoxon"),
      c("Pseudobulk DESeq2", "Pseudobulk edgeR", "Pseudobulk t-test", "Per-cell Wilcoxon"),
      "Pseudobulk methods use donors as replicates. Per-cell tests treat every cell as a replicate, which inflates false positives.",
      c("scrnaseq_014", "scrnaseq_015", "scrnaseq_011"))
    return(out)
  }
  if (identical(spec$task, "peak_study")) {
    kind <- if (identical(spec$assay, "atacseq")) "accessibility" else "binding"
    out[[1L]] <- q(
      "peak_fdr", "Which false-discovery threshold should call a peak?",
      c(0.01, 0.05, 0.10), c("FDR 0.01 (strict)", "FDR 0.05 (standard)", "FDR 0.10 (lenient)"),
      "Stricter thresholds call fewer false peaks but can miss weak real ones.",
      c("atacseq_007", "atacseq_004"))
    out[[2L]] <- q(
      "diff_method", paste0("Which test should find differential ", kind, "?"),
      c("deseq2", "edger", "ttest"), c("DESeq2", "edgeR", "t-test on log-CPM"),
      "DESeq2 tends to give fewer false positives; edgeR is often more sensitive; a simple t-test is a weak baseline with few replicates.",
      c("atacseq_009", "scatac_002", "rnaseq_011"))
    return(out)
  }
  if (panel_gwas || has("gwas")) {
    out[[length(out) + 1L]] <- q(
      "gwas_pcs", "How many principal components should the GWAS use to control population structure?",
      c(0L, 3L, 5L, 10L), c("0 PCs", "3 PCs", "5 PCs", "10 PCs"),
      "Too few PCs can leave structure that creates false hits; too many can absorb real signal. The best number depends on how structured the population is.",
      c("gwas_027", "gwas_004", "gwas_023", "gwas_025"))
    out[[length(out) + 1L]] <- q(
      "gwas_fdr", "Which false-discovery threshold should call a GWAS hit?",
      c(0.01, 0.05, 0.10), c("FDR 0.01 (strict)", "FDR 0.05 (standard)", "FDR 0.10 (lenient)"),
      "Stricter thresholds give fewer false hits but miss more real genes.",
      c("gwas_016", "gwas_007", "gwas_017"))
  }
  if (has("rnaseq")) {
    out[[length(out) + 1L]] <- q(
      "rnaseq_fdr", "Which false-discovery threshold should call a gene differentially expressed?",
      c(0.01, 0.05, 0.10), c("FDR 0.01 (strict)", "FDR 0.05 (standard)", "FDR 0.10 (lenient)"),
      "Published studies use different cut-offs; the choice trades missed genes against false ones.",
      c("rnaseq_012", "rnaseq_011"))
  }
  if (has("genomic_selection")) {
    out[[length(out) + 1L]] <- q(
      "gs_model", "Which genomic prediction model should pick the parents?",
      c("gblup", "rrblup", "bayes_b"), c("GBLUP", "rrBLUP", "Bayes B"),
      "No model wins for every trait; Bayesian models can help when a few genes have large effects.",
      c("genomic_selection_007", "genomic_selection_010", "genomic_selection_021", "genomic_selection_026"))
  }
  out
}

# Reviewed paper guidance for one question (empty until records are reviewed).
.simitall_agent_question_evidence <- function(question, records = .simitall_paper_records()) {
  if (!nrow(records)) return(character())
  hits <- records[records$record_id %in% question$records, , drop = FALSE]
  if (!nrow(hits)) return(character())
  paste0(hits$guidance, " [", hits$source, ", p. ", hits$page, "]")
}

# A recipe that runs the spec once per option of one question, each in its own
# folder and on identical simulated data (same seeds), then compares the truth
# scores side by side. Options that fail (e.g. a missing package) are recorded,
# not fatal.
.simitall_agent_compile_comparison <- function(spec, question_id) {
  questions <- .simitall_agent_spec_questions(spec)
  question <- Filter(function(x) identical(x$id, question_id), questions)
  if (!length(question)) stop("No question '", question_id, "' for this request")
  question <- question[[1L]]
  metric_file <- switch(question_id,
    gwas_pcs = , gwas_fdr = if (identical(spec$task, "gwas")) "truth_recovery_metrics.tsv" else "gwas_truth_recovery_metrics.tsv",
    rnaseq_fdr = "rnaseq_de_truth_recovery_metrics.tsv",
    gs_model = "genomic_selection_truth_recovery_metrics.tsv",
    peak_fdr = "peak_calling_truth_recovery_metrics.tsv",
    diff_method = "differential_truth_recovery_metrics.tsv",
    hvg_fraction = "clustering_truth_recovery_metrics.tsv",
    sc_de_method = "condition_de_truth_recovery_metrics.tsv",
    backcrosses = , selfing = "breeding_truth_recovery_metrics.tsv")
  base <- file.path("results", "agent_compare", paste0(question_id, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  blocks <- character()
  for (i in seq_along(question$values)) {
    value <- question$values[i]
    variant <- spec
    variant$choices[[question_id]] <- if (question_id %in% c("gs_model", "diff_method", "sc_de_method")) value else as.numeric(value)
    variant <- .simitall_agent_compile_spec(variant)
    code <- sub("^[^`]*```r\n", "", variant$recipe)
    code <- sub("\n```[^`]*$", "", code)
    option_dir <- file.path(base, paste0("option_", gsub("[^A-Za-z0-9.]+", "_", value)))
    code <- sub("(?m)^out_dir <- .*$", paste0("out_dir <- ", encodeString(option_dir, quote = "\"")), code, perl = TRUE)
    blocks <- c(blocks,
      paste0("cat(\"\\n=== Option ", i, " of ", length(question$values), ": ", question$labels[i], " ===\\n\")"),
      paste0("status[[", encodeString(value, quote = "\""), "]] <- tryCatch({ local({"),
      code,
      "}); \"ok\" }, error = function(e) { message(\"Option failed: \", conditionMessage(e)); conditionMessage(e) })"
    )
  }
  aggregate <- c(
    "",
    "cat(\"\\n=== Comparing options against the simulation truth ===\\n\")",
    paste0("values <- c(", paste(encodeString(question$values, quote = "\""), collapse = ", "), ")"),
    paste0("labels <- c(", paste(encodeString(question$labels, quote = "\""), collapse = ", "), ")"),
    "rows <- lapply(seq_along(values), function(i) {",
    "  dir <- file.path(base, paste0(\"option_\", gsub(\"[^A-Za-z0-9.]+\", \"_\", values[i])))",
    paste0("  f <- list.files(dir, ", encodeString(paste0(gsub(".", "\\.", metric_file, fixed = TRUE), "$"), quote = "\""), ", full.names = TRUE)[1]"),
    "  if (is.na(f)) return(data.frame(option = labels[i], status = status[[values[i]]], stringsAsFactors = FALSE))",
    "  m <- read.delim(f, check.names = FALSE)[1, , drop = FALSE]",
    "  cbind(data.frame(option = labels[i], status = \"ok\", stringsAsFactors = FALSE), m)",
    "})",
    "cols <- unique(unlist(lapply(rows, names)))",
    "comparison <- do.call(rbind, lapply(rows, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] }))",
    "write.table(comparison, file.path(base, \"option_comparison_truth_recovery_metrics.tsv\"), sep = \"\\t\", row.names = FALSE, quote = FALSE)",
    "score_cols <- intersect(c(\"holdout_accuracy\", \"adjusted_rand_index\", \"nmi\", \"f1\", \"recall\", \"precision\", \"parent1_genome\", \"heterozygosity\"), names(comparison))",
    "if (length(score_cols)) {",
    "  scores <- t(as.matrix(comparison[, score_cols, drop = FALSE])); colnames(scores) <- comparison$option",
    "  png(file.path(base, \"option_comparison.png\"), width = 1500, height = 950, res = 170)",
    "  barplot(scores, beside = TRUE, legend.text = rownames(scores), ylim = c(min(0, scores, na.rm = TRUE), 1.25),",
    "    args.legend = list(x = \"top\", horiz = TRUE, bty = \"n\"), col = c(\"#2F6B5A\", \"#D9984A\", \"#7A9CC6\", \"#B85252\")[seq_len(nrow(scores))],",
    if (question_id %in% c("backcrosses", "selfing")) {
      paste0("    ylab = \"Fraction of genome (simulation truth)\", main = ", encodeString(paste("What each option produces:", question$question), quote = "\""), ", cex.main = 0.8)")
    } else {
      paste0("    ylab = \"Score against the simulation truth\", main = ", encodeString(paste("Which option worked best?", question$question), quote = "\""), ", cex.main = 0.8)")
    },
    "  mtext(\"One simulated dataset per option (same seed): small differences can be chance.\", side = 1, line = 3, cex = 0.75)",
    "  dev.off()",
    "}",
    "print(comparison)"
  )
  code <- c(
    "if (!\"simitall\" %in% loadedNamespaces()) library(simitall)",
    paste0("base <- ", encodeString(base, quote = "\"")),
    "dir.create(base, recursive = TRUE, showWarnings = FALSE)",
    "status <- list()",
    blocks, aggregate
  )
  paste("Spec-compiled option comparison recipe:", "", "```r", paste(code, collapse = "\n"), "```", sep = "\n")
}

# Tunable values for a spec, with labels and ranges. The Shiny tuning panel
# is built from this list, so each workflow controls its own knobs.
.simitall_agent_spec_knobs <- function(spec) {
  d <- .simitall_agent_spec_defaults(spec)
  knob <- function(name, label, value, min, max, step) {
    list(name = name, label = label, value = value, min = min, max = max, step = step)
  }
  if (identical(spec$task, "family_risk")) {
    cond <- spec$condition
    return(list(
      knob("seed", "Random seed", spec$seed, 1, NA, 1),
      knob("n_individuals", "Children per group", as.integer(spec$n_individuals %||% d$n_individuals), 20, 20000, 100),
      knob("prevalence", "Risk, unrelated parents (0-0.5)", cond$prevalence, 0.0001, 0.49, 0.005),
      knob("n_risk_loci", "Recessive risk loci", as.integer(cond$n_risk_loci), 1, 50, 1),
      knob("risk_allele_frequency", "Risk-allele frequency", cond$risk_allele_frequency, 0.001, 0.49, 0.01),
      knob("risk_effect", "Log-odds per homozygous locus", cond$risk_effect, 0.5, 10, 0.5)
    ))
  }
  is_nam <- identical(spec$population$design, "NAM")
  if (identical(spec$task, "breeding")) {
    return(list(
      knob("seed", "Random seed", spec$seed, 1, NA, 1),
      knob("n_lines", "Lines", as.integer(spec$n_lines), 10, 1000, 10)
    ))
  }
  if (identical(spec$task, "single_cell")) {
    return(list(
      knob("seed", "Random seed", spec$seed, 1, NA, 1),
      knob("n_donors", "Donors (split across conditions)", as.integer(spec$n_donors), 4, 24, 2),
      knob("cells_per_donor", "Cells per donor", as.integer(spec$cells_per_donor), 50, 2000, 50),
      knob("n_genes", "Genes", as.integer(spec$n_genes), 100, 3000, 100)
    ))
  }
  if (identical(spec$task, "peak_study")) {
    return(list(
      knob("seed", "Random seed", spec$seed, 1, NA, 1),
      knob("n_peaks", "True peaks", as.integer(spec$n_peaks), 10, 2000, 10),
      knob("replicates", "Replicates per condition", as.integer(spec$replicates), 2, 10, 1),
      knob("n_reads", "Reads per library", as.integer(spec$n_reads), 5000, 2000000, 5000)
    ))
  }
  if (identical(spec$task, "population_study")) {
    needs_trait <- any(c("gwas", "genomic_selection") %in% spec$analyses)
    return(Filter(Negate(is.null), list(
      knob("seed", "Random seed", spec$seed, 1, NA, 1),
      knob("n_individuals", if (is_nam) "NAM lines" else "Accessions",
           as.integer(spec$n_individuals %||% d$n_individuals), 10, 5000, 10),
      if (needs_trait) knob("heritability", "Heritability", spec$trait$heritability, 0.05, 0.95, 0.05),
      if (needs_trait) knob("n_qtl", "Causal markers", as.integer(spec$trait$n_qtl %||% d$n_qtl), 1, 50, 1),
      if ("rnaseq" %in% spec$analyses) knob("n_genes", "Genes measured (RNA-seq)", as.integer(spec$rnaseq$n_genes), 10, 2000, 10)
    )))
  }
  Filter(Negate(is.null), list(
    knob("seed", "Random seed", spec$seed, 1, NA, 1),
    knob("n_individuals", if (is_nam) "NAM lines" else "Accessions",
         as.integer(spec$n_individuals %||% d$n_individuals), 10, 5000, 10),
    knob("heritability", "Heritability", spec$trait$heritability, 0.05, 0.95, 0.05),
    if (!is_nam) knob("n_qtl", "Causal markers", as.integer(spec$trait$n_qtl %||% d$n_qtl), 1, 50, 1)
  ))
}

# The values a spec will actually run with (question values or defaults).
.simitall_agent_spec_effective <- function(spec) {
  knobs <- .simitall_agent_spec_knobs(spec)
  values <- stats::setNames(lapply(knobs, `[[`, "value"), vapply(knobs, `[[`, character(1), "name"))
  for (question in .simitall_agent_spec_questions(spec)) values[[question$id]] <- question$selected
  values
}

.simitall_agent_spec_set <- function(spec, name, value) {
  if (name %in% c("gwas_pcs", "gwas_fdr", "rnaseq_fdr", "gs_model", "peak_fdr", "diff_method", "hvg_fraction", "sc_de_method", "backcrosses", "selfing")) {
    spec$choices[[name]] <- if (name %in% c("gs_model", "diff_method", "sc_de_method")) as.character(value) else as.numeric(value)
    return(spec)
  }
  value <- as.numeric(value)
  untag <- function(field) setdiff(spec$defaulted, field)
  switch(
    name,
    seed = { spec$seed <- value; spec$defaulted <- untag("seed") },
    n_individuals = { spec$n_individuals <- value; spec$defaulted <- untag("n_individuals") },
    heritability = { spec$trait$heritability <- value; spec$defaulted <- untag("trait.heritability") },
    n_qtl = { spec$trait$n_qtl <- value },
    n_genes = { if (identical(spec$task, "single_cell")) { spec$n_genes <- value; spec$defaulted <- untag("n_genes") } else spec$rnaseq$n_genes <- value },
    n_peaks = { spec$n_peaks <- value; spec$defaulted <- untag("n_peaks") },
    replicates = { spec$replicates <- value; spec$defaulted <- untag("replicates") },
    n_reads = { spec$n_reads <- value; spec$defaulted <- untag("n_reads") },
    n_donors = { spec$n_donors <- value; spec$defaulted <- untag("n_donors") },
    cells_per_donor = { spec$cells_per_donor <- value; spec$defaulted <- untag("cells_per_donor") },
    n_lines = { spec$n_lines <- value; spec$defaulted <- untag("n_lines") },
    prevalence = { spec$condition$prevalence <- value; spec$defaulted <- untag("condition.prevalence") },
    n_risk_loci = { spec$condition$n_risk_loci <- value; spec$defaulted <- untag("condition.risk_loci") },
    risk_allele_frequency = { spec$condition$risk_allele_frequency <- value; spec$defaulted <- untag("condition.risk_loci") },
    risk_effect = { spec$condition$risk_effect <- value; spec$defaulted <- untag("condition.risk_loci") }
  )
  spec
}

# Apply user overrides (e.g. from the Shiny tuning panel) to a spec, then
# re-validate and recompile. Only values that differ from the effective
# value are applied, so an untouched knob never rewrites anything. Invalid
# overrides come back as spec$errors with an empty recipe, never as code.
.simitall_agent_recompile_spec <- function(spec, overrides = list()) {
  current <- .simitall_agent_spec_effective(spec)
  for (name in names(current)) {
    value <- overrides[[name]]
    valid <- !is.null(value) && length(value) == 1L && !is.na(value) && nzchar(as.character(value))
    # Compare numbers as numbers ("3" == 3L) and text as text ("gblup" vs "rrblup").
    num_value <- suppressWarnings(as.numeric(value))
    num_current <- suppressWarnings(as.numeric(current[[name]]))
    same <- if (valid && !is.na(num_value) && length(num_current) == 1L && !is.na(num_current)) {
      isTRUE(all.equal(num_value, num_current))
    } else {
      identical(as.character(value), as.character(current[[name]]))
    }
    changed <- valid && !same
    if (changed) spec <- .simitall_agent_spec_set(spec, name, value)
  }
  spec$errors <- .simitall_agent_validate_spec(spec)
  spec$recipe <- ""
  if (length(spec$errors)) return(spec)
  .simitall_agent_compile_spec(spec)
}

# ---------------------------------------------------------------------------
# Simulation helpers used by compiled recipes.
# ---------------------------------------------------------------------------

# Simulate fully inbred accessions from structured subpopulations.
# Ancestral frequencies p ~ Uniform(0.05, 0.5); subpopulation frequencies are
# Balding-Nichols draws Beta(p(1-F)/F, (1-p)(1-F)/F). Accessions are homozygous
# (0/0 or 1/1), as for a selfing species. Markers are independent: no LD.
.simitall_simulate_inbred_panel <- function(
    vcf_out, n_accessions = 200L, n_markers = 2000L, n_chromosomes = 5L,
    chromosome_length_bp = 1000000L, n_subpops = 4L, fst = 0.10, seed = 1L,
    maf_min = 0.05) {
  stopifnot(n_accessions >= 2, n_markers >= 10, n_subpops >= 1,
            fst > 0, fst < 1, n_markers <= n_chromosomes * chromosome_length_bp)
  set.seed(seed)
  subpop <- rep_len(seq_len(n_subpops), n_accessions)
  subpop <- sample(subpop)
  ancestral <- stats::runif(n_markers, 0.05, 0.5)
  a <- ancestral * (1 - fst) / fst
  b <- (1 - ancestral) * (1 - fst) / fst
  sub_freq <- vapply(seq_len(n_subpops), function(k) stats::rbeta(n_markers, a, b), numeric(n_markers))
  sub_freq <- matrix(sub_freq, nrow = n_markers)
  alt <- matrix(stats::runif(n_markers * n_accessions) < sub_freq[, subpop, drop = FALSE],
                nrow = n_markers)
  dosage <- 2L * alt
  # Drop markers that are rare after sampling so the GWAS sees usable sites.
  freq <- rowMeans(dosage) / 2
  keep <- pmin(freq, 1 - freq) >= maf_min
  if (sum(keep) < 10L) stop("Too few common markers survived; increase n_markers or n_accessions")
  dosage <- dosage[keep, , drop = FALSE]

  chrom_index <- sort(sample(seq_len(n_chromosomes), nrow(dosage), replace = TRUE))
  positions <- unlist(lapply(seq_len(n_chromosomes), function(k) {
    sort(sample.int(chromosome_length_bp, sum(chrom_index == k)))
  }))
  chrom <- paste0("chr", chrom_index)
  ids <- paste0(chrom, "_", positions)
  samples <- sprintf("acc%04d", seq_len(n_accessions))
  gt <- ifelse(dosage == 2L, "1/1", "0/0")

  dir.create(dirname(vcf_out), recursive = TRUE, showWarnings = FALSE)
  header <- c(
    "##fileformat=VCFv4.2",
    "##source=simitall synthetic inbred panel (Balding-Nichols)",
    paste0("##contig=<ID=chr", seq_len(n_chromosomes), ",length=", chromosome_length_bp, ">"),
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples), collapse = "\t")
  )
  body <- paste(chrom, positions, ids, "A", "G", ".", "PASS", ".", "GT",
                apply(gt, 1L, paste, collapse = "\t"), sep = "\t")
  writeLines(c(header, body), vcf_out)

  population_file <- sub("\\.vcf$", ".population.tsv", vcf_out)
  utils::write.table(
    data.frame(sample = samples, subpopulation = paste0("pop", subpop)),
    population_file, sep = "\t", row.names = FALSE, quote = FALSE
  )
  list(vcf = vcf_out, population = population_file,
       n_markers = nrow(dosage), n_accessions = n_accessions)
}

# Score a list of called genes against a truth set: recall, precision, and
# counts, so every analysis stage reports how much of the truth it recovered.
.simitall_score_gene_recovery <- function(called, truth, tested, analysis = "analysis") {
  called <- unique(as.character(called)); truth <- intersect(unique(as.character(truth)), tested)
  hit <- intersect(called, truth)
  data.frame(
    analysis = analysis, genes_tested = length(unique(tested)), true_genes = length(truth),
    called = length(called), recovered = length(hit), missed = length(truth) - length(hit),
    false_positives = length(setdiff(called, truth)),
    recall = if (length(truth)) length(hit) / length(truth) else NA_real_,
    precision = if (length(called)) length(hit) / length(called) else NA_real_,
    stringsAsFactors = FALSE
  )
}

# Additive trait with known causal markers at a programmed heritability.
.simitall_simulate_additive_trait <- function(
    genotype, variants, heritability = 0.6, n_qtl = 5L, trait_name = "trait",
    seed = 1L, maf_min = 0.10) {
  stopifnot(heritability > 0, heritability < 1, n_qtl >= 1)
  set.seed(seed)
  freq <- rowMeans(genotype, na.rm = TRUE) / 2
  maf <- pmin(freq, 1 - freq)
  eligible <- which(is.finite(maf) & maf >= maf_min)
  if (length(eligible) < n_qtl) stop("Only ", length(eligible), " markers are eligible for ", n_qtl, " QTL")
  causal <- sort(sample(eligible, n_qtl))
  effects <- stats::rnorm(n_qtl, 0, 1)
  G <- genotype[causal, , drop = FALSE]
  G[is.na(G)] <- 0
  genetic_value <- as.numeric(crossprod(G, effects))
  if (stats::var(genetic_value) <= 0) stop("Causal markers produced no genetic variance")
  residual_sd <- sqrt(stats::var(genetic_value) * (1 - heritability) / heritability)
  phenotype <- data.frame(sample = colnames(genotype), stringsAsFactors = FALSE)
  phenotype[[trait_name]] <- genetic_value + stats::rnorm(ncol(genotype), 0, residual_sd)
  phenotype$genetic_value <- genetic_value
  truth <- data.frame(
    marker_id = variants$id[causal], seqname = variants$seqname[causal],
    pos = variants$pos[causal], effect = effects, maf = maf[causal],
    stringsAsFactors = FALSE
  )
  list(phenotype = phenotype, truth = truth)
}
