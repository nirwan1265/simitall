# Spec pipeline: question -> spec -> validate -> resolve -> compile.

test_that("maize NAM GWAS extracts a NAM spec and uses the bundled panel", {
  spec <- .simitall_agent_spec_pipeline("run a gwas for a phenotype in maize NAM", data_source = "synthetic")
  expect_true(spec$supported)
  expect_length(spec$errors, 0L)
  expect_equal(spec$organism, "maize")
  expect_equal(spec$population$design, "NAM")
  expect_equal(spec$population$panel_id, "maize_nam")
  expect_equal(spec$resolution$mode, "bundled")
  expect_match(spec$recipe, "\"--input_mode\", \"package_demo\"", fixed = TRUE)
  expect_true(all(c("trait.heritability", "n_individuals", "seed") %in% spec$defaulted))
})

test_that("auto data source uses a bundled panel when one exists", {
  spec <- .simitall_agent_spec_pipeline("run a gwas for a phenotype in maize NAM", data_source = "auto")
  expect_equal(spec$resolution$status, "ready")
  expect_equal(spec$resolution$mode, "bundled")
})

test_that("Arabidopsis 1001 is simulated under synthetic and asked for under auto", {
  q <- "run a gwas for flowering time in arabidopsis 1001 with 300 accessions, heritability 0.4"
  synthetic <- .simitall_agent_spec_pipeline(q, data_source = "synthetic")
  expect_equal(synthetic$population$design, "diversity_panel")
  expect_equal(synthetic$resolution$mode, "simulate")
  expect_equal(synthetic$n_individuals, 300)
  expect_equal(synthetic$trait$heritability, 0.4)
  expect_equal(synthetic$trait$name, "flowering_time")
  expect_match(synthetic$recipe, ".simitall_simulate_inbred_panel", fixed = TRUE)
  expect_match(.simitall_agent_spec_plan(synthetic), "does not reproduce real accessions", fixed = TRUE)

  auto <- .simitall_agent_spec_pipeline(q, data_source = "auto")
  expect_equal(auto$resolution$status, "needs_input")
  expect_match(auto$resolution$action, "does not bundle an Arabidopsis 1001 Genomes", fixed = TRUE)
})

test_that("panel names are not parsed as sample sizes", {
  spec <- .simitall_agent_extract_spec("gwas on arabidopsis 1001 genomes accessions")
  expect_null(spec$n_individuals)
})

test_that("bad parameters fail loudly instead of being ignored", {
  h2 <- .simitall_agent_spec_pipeline("run a gwas in maize NAM with heritability 1.0", data_source = "synthetic")
  expect_true(any(grepl("^trait.heritability", h2$errors)))
  expect_equal(h2$recipe, "")

  qtl <- .simitall_agent_spec_pipeline("run a gwas in maize NAM with 5 QTL", data_source = "synthetic")
  expect_true(any(grepl("n_qtl: not configurable for the NAM runner", qtl$errors, fixed = TRUE)))
})

test_that("unsupported shapes return to the legacy router", {
  expect_false(.simitall_agent_spec_pipeline("simulate RNA-seq for wheat")$supported)
  expect_false(.simitall_agent_spec_pipeline("binary trait gwas in maize NAM", data_source = "synthetic")$supported)
  expect_false(.simitall_agent_spec_pipeline(
    "maize NAM drought vs irrigated GxE with RNA-seq", data_source = "synthetic"
  )$supported)
})

test_that("provided diversity-panel VCF is used directly", {
  vcf <- tempfile(fileext = ".vcf")
  .simitall_simulate_inbred_panel(vcf, n_accessions = 30L, n_markers = 200L, seed = 3L)
  spec <- .simitall_agent_spec_pipeline(
    "gwas for flowering time in arabidopsis 1001", input_files = vcf, data_source = "provided"
  )
  expect_equal(spec$resolution$mode, "provided")
  expect_match(spec$recipe, normalizePath(vcf), fixed = TRUE)
  expect_false(grepl(".simitall_simulate_inbred_panel", spec$recipe, fixed = TRUE))
})

test_that("synthetic inbred panel round-trips through the genotype reader", {
  vcf <- tempfile(fileext = ".vcf")
  panel <- .simitall_simulate_inbred_panel(vcf, n_accessions = 60L, n_markers = 400L, n_subpops = 3L, seed = 5L)
  geno <- .simitall_read_genotypes(panel$vcf)
  expect_equal(ncol(geno$genotype), 60L)
  expect_true(all(geno$genotype %in% c(0, 2)))
  trait <- .simitall_simulate_additive_trait(geno$genotype, geno$variants, heritability = 0.5, n_qtl = 4L, seed = 6L)
  expect_equal(nrow(trait$truth), 4L)
  realized <- stats::var(trait$phenotype$genetic_value) / stats::var(trait$phenotype$trait)
  expect_gt(realized, 0.3)
  expect_lt(realized, 0.7)
})

# Validator: prose parentheses are not function calls.

validator_knowledge <- function() {
  knowledge <- file.path(tempdir(), "simitall-agent-validator")
  dir.create(file.path(knowledge, "workflows"), recursive = TRUE, showWarnings = FALSE)
  writeLines(c("# GWAS workflow", "Retain causal truth."), file.path(knowledge, "workflows", "gwas.Rmd"))
  search_simitall_knowledge("GWAS", knowledge_dir = knowledge)
}

test_that("validator accepts normal prose with parentheses", {
  answer <- paste(
    "Interpreted design: one population (NAM, 500 lines).",
    "Proposed stages: 1. Simulate the population with `simulate_breeding()`.",
    "Required inputs (FASTA and map) are bundled."
  )
  out <- .simitall_agent_validate_answer(answer, "q", validator_knowledge(), allow_code = FALSE)
  expect_false(attr(out, "simitall_validation")$fallback)
})

test_that("validator still catches invented functions and external tools in code", {
  invented <- .simitall_agent_validate_answer(
    "Then run benchmark_eqtls() on the output.", "q",
    validator_knowledge(), allow_code = FALSE
  )
  expect_true(attr(invented, "simitall_validation")$fallback)
  expect_true(any(grepl("benchmark_eqtls", attr(invented, "simitall_validation")$violations)))

  code <- .simitall_agent_validate_answer(
    "```r\nx <- paste0('a', 1)\nGAPIT(Y = y)\n```", "q",
    validator_knowledge(), allow_code = TRUE
  )
  expect_true(attr(code, "simitall_validation")$fallback)

  base_only <- .simitall_agent_validate_answer(
    "```r\nout <- file.path('results', paste0('run', 1))\nsimulate_breeding(panel, out)\n```", "q",
    data.frame(), allow_code = TRUE
  )
  expect_false(attr(base_only, "simitall_validation")$fallback)
})

# Plan prose and compiled code come from the same spec.

test_that("spec-handled answers describe the code that will run", {
  knowledge <- file.path(tempdir(), "simitall-agent-plan")
  dir.create(file.path(knowledge, "workflows"), recursive = TRUE, showWarnings = FALSE)
  writeLines(c("# GWAS workflow", "Use simulate_gwas_cohort()."), file.path(knowledge, "workflows", "gwas.Rmd"))
  answer <- simitall_ask(
    "run a gwas for flowering time in arabidopsis 1001 with 300 accessions",
    provider = "local", include_code = TRUE, data_source = "synthetic", knowledge_dir = knowledge
  )
  expect_equal(answer$plan_source, "spec")
  plan <- sub("\n\nSpec-compiled.*$", "", answer$answer)
  expect_match(plan, "Simulate 300 inbred accessions in 4 subpopulations", fixed = TRUE)
  expect_match(plan, "- Size: 300 accessions", fixed = TRUE)
  expect_false(grepl("simulate_gwas_cohort", plan, fixed = TRUE))
  expect_match(answer$answer, "n_lines = 300L", fixed = TRUE)
})

test_that("spec plans surface missing inputs and parameter errors", {
  knowledge <- file.path(tempdir(), "simitall-agent-plan")
  auto <- simitall_ask("run a gwas for flowering time in arabidopsis 1001",
                       provider = "local", include_code = TRUE, knowledge_dir = knowledge)
  expect_equal(auto$preflight$status, "needs_input")
  expect_match(auto$answer, "does not bundle an Arabidopsis 1001 Genomes", fixed = TRUE)
  expect_false(grepl("```", auto$answer, fixed = TRUE))

  bad <- simitall_ask("run a gwas in maize NAM with heritability 1.0",
                      provider = "local", include_code = TRUE, data_source = "synthetic", knowledge_dir = knowledge)
  expect_equal(bad$preflight$status, "invalid_spec")
  expect_match(bad$answer, "Fix before running: trait.heritability", fixed = TRUE)
})

# Tuning overrides go through validation and recompilation.

test_that("untouched tuning knobs leave the recipe unchanged", {
  spec <- .simitall_agent_spec_pipeline("run a gwas for a phenotype in maize NAM", data_source = "synthetic")
  same <- .simitall_agent_recompile_spec(spec, .simitall_agent_spec_effective(spec))
  expect_identical(same$recipe, spec$recipe)
  expect_match(same$recipe, "\"--n_lines\", \"100\"", fixed = TRUE)
})

test_that("tuning overrides recompile the recipe and plan", {
  spec <- .simitall_agent_spec_pipeline("run a gwas for flowering time in arabidopsis 1001", data_source = "synthetic")
  tuned <- .simitall_agent_recompile_spec(spec, list(n_individuals = 250, heritability = 0.3, n_qtl = 8, seed = 11))
  expect_length(tuned$errors, 0L)
  expect_match(tuned$recipe, "n_lines = 250L", fixed = TRUE)
  expect_match(tuned$recipe, "heritability = 0.3", fixed = TRUE)
  expect_match(tuned$recipe, "n_qtl = 8L", fixed = TRUE)
  expect_match(tuned$recipe, "seed = 11", fixed = TRUE)
  expect_false("n_individuals" %in% tuned$defaulted)
  expect_match(.simitall_agent_spec_plan(tuned), "- Size: 250 accessions\n", fixed = TRUE)
})

test_that("invalid tuning overrides are rejected instead of run", {
  spec <- .simitall_agent_spec_pipeline("run a gwas for a phenotype in maize NAM", data_source = "synthetic")
  bad <- .simitall_agent_recompile_spec(spec, list(heritability = 1.5, n_individuals = 3))
  expect_equal(bad$recipe, "")
  expect_true(any(grepl("^trait.heritability", bad$errors)))
  expect_true(any(grepl("^n_individuals", bad$errors)))
})

# Family-risk questions in plain language.

test_that("plain-language cousin questions become family-risk specs", {
  spec <- .simitall_agent_spec_pipeline("what is the chance my kids get a disease if I marry my cousin?")
  expect_true(spec$supported)
  expect_equal(spec$task, "family_risk")
  expect_equal(spec$groups, c("unrelated", "first_cousin_descendant"))
  expect_equal(spec$resolution$mode, "simulate")
  expect_match(spec$recipe, "risk_model = \"recessive\"", fixed = TRUE)
  expect_match(spec$recipe, "simulate_human_marker_frequencies", fixed = TRUE)
  expect_false(grepl("12 * F", spec$recipe, fixed = TRUE))
})

test_that("family-risk questions parse condition, prevalence, groups, and size", {
  spec <- .simitall_agent_spec_pipeline(
    "If first cousins and second cousins have children, how likely is cystic fibrosis? assume 1 in 2500, 500 families each"
  )
  expect_equal(spec$condition$name, "cystic fibrosis")
  expect_equal(spec$condition$prevalence, 1 / 2500)
  expect_equal(spec$groups, c("unrelated", "distantly_related", "first_cousin_descendant"))
  expect_equal(spec$n_individuals, 500)
  expect_match(.simitall_agent_spec_plan(spec), "Risk for children of unrelated parents: 0.04%", fixed = TRUE)
})

test_that("IRF6 family questions use the bundled marker frequencies", {
  spec <- .simitall_agent_spec_pipeline(
    "Simulate the probability of a cleft lip trait in the human IRF6 region for children of first cousins and unrelated parents"
  )
  expect_equal(spec$resolution$mode, "bundled")
  expect_match(spec$recipe, "demo_human_irf6_chr1.vcf.gz", fixed = TRUE)
})

test_that("unsupported relationships and inheritance fail loudly", {
  uncle <- .simitall_agent_spec_pipeline("risk of disease for children of an uncle and niece")
  expect_true(any(grepl("uncle_niece not supported", uncle$errors, fixed = TRUE)))
  expect_equal(uncle$recipe, "")
  dominant <- .simitall_agent_spec_pipeline("probability of a dominant disease in children of cousins")
  expect_true(any(grepl("'dominant' inheritance is not supported", dominant$errors, fixed = TRUE)))
})

test_that("family-risk knobs recompile and validate", {
  spec <- .simitall_agent_spec_pipeline("chance of disease for children of first cousins")
  knobs <- vapply(.simitall_agent_spec_knobs(spec), `[[`, character(1), "name")
  expect_true(all(c("prevalence", "n_risk_loci", "risk_allele_frequency", "risk_effect") %in% knobs))
  tuned <- .simitall_agent_recompile_spec(spec, list(prevalence = 0.002, n_individuals = 300))
  expect_match(tuned$recipe, "baseline_probability = 0.002", fixed = TRUE)
  expect_match(tuned$recipe, "n_per_group = 300L", fixed = TRUE)
  bad <- .simitall_agent_recompile_spec(spec, list(prevalence = 0.7))
  expect_true(any(grepl("^condition.prevalence", bad$errors)))
})

# Plain-language plant gene hunting.

test_that("plain-language gene-hunting questions become GWAS specs", {
  qs <- c(
    "find the gene for drought tolerance in maize NAM",
    "which genes control drought tolerance in maize? use a NAM population with 300 lines",
    "map the genes behind flowering time in arabidopsis 1001 accessions"
  )
  specs <- lapply(qs, .simitall_agent_spec_pipeline, data_source = "synthetic")
  expect_true(all(vapply(specs, function(s) identical(s$task, "gwas") && isTRUE(s$supported), logical(1))))
  expect_equal(specs[[2]]$n_individuals, 300)
  expect_equal(specs[[3]]$population$design, "diversity_panel")
})

test_that("gene hunting does not swallow family-risk questions", {
  expect_equal(.simitall_agent_spec_pipeline("what is the chance my kids get a disease if I marry my cousin?")$task, "family_risk")
  expect_false(isTRUE(.simitall_agent_spec_pipeline("simulate RNA-seq for wheat")$supported))
  maize <- .simitall_agent_spec_pipeline("simulate RNA-seq for maize")
  expect_equal(maize$analyses, "rnaseq")
  expect_true("population.design" %in% maize$defaulted)
})

test_that("the NAM plan names the runner's fixed trait honestly", {
  spec <- .simitall_agent_spec_pipeline("find the gene for yield in maize NAM", data_source = "synthetic")
  plan <- .simitall_agent_spec_plan(spec)
  expect_match(plan, "`drought_response_index`", fixed = TRUE)
  expect_match(plan, "you asked about yield", fixed = TRUE)
})

# Composable population studies: run exactly the analyses the prompt names.

test_that("RNA-seq-only gene hunts run RNA-seq and not GWAS", {
  spec <- .simitall_agent_spec_pipeline(
    "Find the gene for drought tolerance in maize using RNAseq. Use a NAM population with 300 lines.",
    data_source = "synthetic"
  )
  expect_true(spec$supported)
  expect_equal(spec$analyses, "rnaseq")
  expect_equal(spec$n_individuals, 300)
  expect_equal(spec$rnaseq$conditions, c("irrigated", "drought"))
  expect_match(spec$recipe, "simulate_rnaseq_from_gwas", fixed = TRUE)
  expect_match(spec$recipe, "benchmark_eqtl", fixed = TRUE)
  expect_false(grepl("analyze_gwas", spec$recipe, fixed = TRUE))
  expect_false(grepl("run_genomic_selection", spec$recipe, fixed = TRUE))
  plan <- .simitall_agent_spec_plan(spec)
  expect_match(plan, "Not requested, so not run: GWAS, genomic selection", fixed = TRUE)
  expect_match(plan, "It does not show which gene causes the trait", fixed = TRUE)
})

test_that("named analyses compose in order", {
  both <- .simitall_agent_spec_pipeline(
    "I want to do RNA-seq on this NAM maize population and also a GWAS for drought tolerance, 300 lines",
    data_source = "synthetic"
  )
  expect_equal(both$analyses, c("gwas", "rnaseq"))
  expect_lt(regexpr("analyze_gwas", both$recipe), regexpr("simulate_rnaseq_from_gwas", both$recipe))
  gs <- .simitall_agent_spec_pipeline(
    "maize NAM, 200 lines: GWAS for plant height and genomic selection to pick the best parents",
    data_source = "synthetic"
  )
  expect_equal(gs$analyses, c("gwas", "genomic_selection"))
  expect_match(gs$recipe, "run_genomic_selection", fixed = TRUE)
  expect_false(grepl("simulate_rnaseq_from_gwas", gs$recipe, fixed = TRUE))
  expect_match(.simitall_agent_spec_plan(gs), "Not requested, so not run: RNA-seq", fixed = TRUE)
})

test_that("GWAS-only requests keep the versioned NAM runner", {
  spec <- .simitall_agent_spec_pipeline("find the gene for drought tolerance in maize NAM", data_source = "synthetic")
  expect_equal(spec$task, "gwas")
  expect_match(spec$recipe, "fig15_maize_nam_quantitative_gwas.R", fixed = TRUE)
})

test_that("population-study knobs depend on the requested analyses", {
  rna <- .simitall_agent_spec_pipeline("RNA-seq on maize NAM for drought", data_source = "synthetic")
  knobs <- vapply(.simitall_agent_spec_knobs(rna), `[[`, character(1), "name")
  expect_true("n_genes" %in% knobs)
  expect_false("heritability" %in% knobs)
  tuned <- .simitall_agent_recompile_spec(rna, list(n_genes = 80))
  expect_match(tuned$recipe, "n_genes = 80L", fixed = TRUE)
  bad <- .simitall_agent_recompile_spec(rna, list(n_genes = 3))
  expect_true(any(grepl("^rnaseq.n_genes", bad$errors)))
})

test_that("a compiled RNA-seq study runs end to end and scores its truth", {
  skip_if_not_installed("ggplot2")
  spec <- .simitall_agent_spec_pipeline(
    "RNA-seq for heat stress in arabidopsis 1001 accessions, 30 genes, 60 accessions", data_source = "synthetic"
  )
  code <- sub("\\n```[^`]*$", "", sub("^[^`]*```r\\n", "", spec$recipe))
  old <- setwd(tempdir()); on.exit(setwd(old), add = TRUE)
  suppressMessages(capture.output(eval(parse(text = code), envir = new.env())))
  out <- file.path(tempdir(), "results", "agent_study", "arabidopsis_1001_rnaseq")
  de <- read.delim(file.path(out, "rnaseq_de_truth_recovery_metrics.tsv"))
  expect_equal(de$genes_tested, 30L)
  expect_gt(de$true_genes, 0L)
  expect_true(file.exists(file.path(out, "eqtl.eqtl_metrics.tsv")))
  expect_true(file.exists(file.path(out, "figures", "rnaseq_volcano.png")))
})

# Judge-friendly breeding prompts.

test_that("a plain breeding request runs gene finding and parent selection", {
  spec <- .simitall_agent_spec_pipeline(
    "Help me breed higher-yielding maize. Find the genes that control yield and pick the best parents to cross. Use 300 plants."
  )
  expect_true(spec$supported)
  expect_equal(spec$analyses, c("gwas", "genomic_selection"))
  expect_equal(spec$trait$name, "yield")
  expect_equal(spec$n_individuals, 300)
  expect_equal(spec$population$design, "NAM")
  expect_true("population.design" %in% spec$defaulted)
  expect_match(.simitall_agent_spec_plan(spec), "no population was named", fixed = TRUE)
  expect_match(spec$recipe, "analyze_gwas", fixed = TRUE)
  expect_match(spec$recipe, "run_genomic_selection", fixed = TRUE)
  expect_false(grepl("simulate_rnaseq_from_gwas", spec$recipe, fixed = TRUE))
})

test_that("a named population is never overridden by the maize default", {
  spec <- .simitall_agent_spec_pipeline("find the gene for drought tolerance in maize NAM", data_source = "synthetic")
  expect_false("population.design" %in% spec$defaulted)
})

# Synonyms in inst/agent/aliases.tsv reach the spec extractor.

test_that("plain-language synonyms are normalised and shown in the plan", {
  corn <- .simitall_agent_spec_pipeline("find the genes for yield in corn")
  expect_equal(corn$organism, "maize")
  expect_match(.simitall_agent_spec_plan(corn), "Read as: corn -> maize", fixed = TRUE)
  expect_equal(corn$question, "find the genes for yield in corn")

  cress <- .simitall_agent_spec_pipeline("find the genes for flowering in thale cress accessions", data_source = "synthetic")
  expect_equal(cress$organism, "arabidopsis")
  expect_equal(cress$trait$name, "flowering_time")

  water <- .simitall_agent_spec_pipeline("maize NAM: transcriptome profiling under water deficit", data_source = "synthetic")
  expect_equal(water$analyses, "rnaseq")
  expect_equal(water$rnaseq$conditions, c("irrigated", "drought"))
})

test_that("family-risk specs keep the original question after normalisation", {
  spec <- .simitall_agent_spec_pipeline("chance a child gets a disease if the parents are consanguineous")
  expect_equal(spec$task, "family_risk")
  expect_equal(spec$question, "chance a child gets a disease if the parents are consanguineous")
})

# Paper records: only reviewed rows reach a plan.

test_that("only reviewed paper records reach the plan", {
  path <- tempfile(fileext = ".tsv")
  cols <- c("record_id", "domain", "decision", "guidance", "numbers", "applies_to", "simitall_use",
            "source", "doi", "page", "confidence", "flag", "review_status", "reviewed_by")
  rows <- rbind(
    c("gwas_901", "gwas", "pcs", "Use 3-5 PCs.", "", "GWAS", "default", "Test 2020", "x", "3", "medium", "", "reviewed", "NT"),
    c("gwas_902", "gwas", "draft", "Unchecked advice.", "", "GWAS", "default", "Test 2021", "x", "1", "high", "", "unreviewed", ""),
    c("gwas_903", "gwas", "bad", "Wrong statement.", "", "GWAS", "glossary", "Test 2022", "x", "4", "low", "ERROR", "excluded", ""),
    c("rnaseq_901", "rnaseq", "reps", "Use 3 replicates.", "", "RNA-seq", "validation_check", "Test 2023", "x", "8", "high", "", "reviewed", "NT")
  )
  writeLines(c(paste(cols, collapse = "\t"), apply(rows, 1, paste, collapse = "\t")), path)
  records <- .simitall_paper_records(path)
  expect_equal(records$record_id, c("gwas_901", "rnaseq_901"))

  gwas <- .simitall_agent_spec_pipeline("find the gene for drought tolerance in maize NAM", data_source = "synthetic")
  lines <- .simitall_agent_paper_guidance(gwas, records)
  expect_length(lines, 1L)
  expect_match(lines, "Use 3-5 PCs. [Test 2020, p. 3]", fixed = TRUE)

  rna <- .simitall_agent_spec_pipeline("RNA-seq on maize NAM for drought", data_source = "synthetic")
  expect_match(.simitall_agent_paper_guidance(rna, records), "Use 3 replicates.", fixed = TRUE)
  expect_length(.simitall_agent_paper_guidance(.simitall_agent_spec_pipeline("chance of disease for children of cousins"), records), 0L)
})

test_that("with no reviewed records, plans carry no paper guidance", {
  records <- .simitall_paper_records(status = c("reviewed", "unreviewed", "excluded"))
  skip_if(!nrow(records), "paper records not found")
  skip_if(any(records$review_status == "reviewed"), "some records have been reviewed")
  expect_false(grepl("Paper guidance", .simitall_agent_spec_plan(
    .simitall_agent_spec_pipeline("find the gene for drought tolerance in maize NAM", data_source = "synthetic")), fixed = TRUE))
})

# Guided mode: questions, choices, and option comparisons.

test_that("questions appear only for workflows whose code can honour them", {
  panel <- .simitall_agent_spec_pipeline("find the genes for flowering time in arabidopsis 1001 accessions", data_source = "synthetic")
  expect_equal(vapply(.simitall_agent_spec_questions(panel), `[[`, "", "id"), c("gwas_pcs", "gwas_fdr"))
  study <- .simitall_agent_spec_pipeline("maize NAM, 200 lines: GWAS for yield, RNA-seq for drought, and pick the best parents", data_source = "synthetic")
  expect_equal(vapply(.simitall_agent_spec_questions(study), `[[`, "", "id"), c("gwas_pcs", "gwas_fdr", "rnaseq_fdr", "gs_model"))
  nam_runner <- .simitall_agent_spec_pipeline("find the gene for drought tolerance in maize NAM", data_source = "synthetic")
  expect_length(.simitall_agent_spec_questions(nam_runner), 0L)
  expect_length(.simitall_agent_spec_questions(.simitall_agent_spec_pipeline("chance of disease for children of cousins")), 0L)
})

test_that("defaults keep the previous recipes and choices recompile them", {
  panel <- .simitall_agent_spec_pipeline("find the genes for flowering time in arabidopsis 1001 accessions", data_source = "synthetic")
  expect_match(panel$recipe, "n_pcs = 3L", fixed = TRUE)
  tuned <- .simitall_agent_recompile_spec(panel, list(gwas_pcs = "10", gwas_fdr = "0.01"))
  expect_match(tuned$recipe, "n_pcs = 10L", fixed = TRUE)
  expect_match(tuned$recipe, "fdr_threshold = 0.01", fixed = TRUE)
  study <- .simitall_agent_spec_pipeline("maize NAM: GWAS for yield and pick the best parents", data_source = "synthetic")
  expect_match(study$recipe, "fixed_effects = \"family\", n_pcs = 2L", fixed = TRUE)
  expect_match(.simitall_agent_recompile_spec(study, list(gs_model = "rrblup"))$recipe, "model = \"rrblup\"", fixed = TRUE)
})

test_that("an option comparison runs every option in its own folder and aggregates", {
  panel <- .simitall_agent_spec_pipeline("find the genes for flowering time in arabidopsis 1001 accessions", data_source = "synthetic")
  recipe <- .simitall_agent_compile_comparison(panel, "gwas_pcs")
  for (k in c(0, 3, 5, 10)) expect_match(recipe, paste0("n_pcs = ", k, "L"), fixed = TRUE)
  expect_match(recipe, "option_comparison_truth_recovery_metrics.tsv", fixed = TRUE)
  expect_equal(lengths(regmatches(recipe, gregexpr("tryCatch\\(\\{ local\\(\\{", recipe))), 4L)
  expect_error(.simitall_agent_compile_comparison(panel, "gs_model"), "No question")
})

test_that("question evidence shows only reviewed records", {
  panel <- .simitall_agent_spec_pipeline("find the genes for flowering time in arabidopsis 1001 accessions", data_source = "synthetic")
  question <- .simitall_agent_spec_questions(panel)[[1L]]
  records <- data.frame(record_id = c("gwas_027", "gwas_004"), guidance = c("Use 3-5 PCs.", "Control structure."),
                        source = c("Yin 2020", "Sahito 2024"), page = c("3", "7"), review_status = c("reviewed", "unreviewed"))
  expect_equal(.simitall_agent_question_evidence(question, records[records$review_status == "reviewed", ]),
               "Use 3-5 PCs. [Yin 2020, p. 3]")
})

# ChIP-seq and ATAC-seq: plain-English prompts compile to plan and code.

test_that("ChIP and ATAC prompts become peak-study specs", {
  atac <- .simitall_agent_spec_pipeline("simulate ATAC-seq for maize roots under drought, 4 replicates")
  expect_equal(atac$task, "peak_study"); expect_equal(atac$assay, "atacseq")
  expect_equal(atac$conditions, c("control", "drought")); expect_equal(atac$replicates, 4)
  expect_match(atac$recipe, "simulate_atacseq(", fixed = TRUE)
  expect_match(atac$recipe, "control_assay = NULL", fixed = TRUE)
  chip <- .simitall_agent_spec_pipeline("find differential binding sites for a transcription factor with ChIP-seq, control vs heat")
  expect_equal(chip$target, "TF")
  expect_match(chip$recipe, "control_assay = \"input\"", fixed = TRUE)
  mark <- .simitall_agent_spec_pipeline("ChIP-seq for H3K27me3 in arabidopsis, 150 peaks")
  expect_equal(mark$histone_mark, "H3K27me3"); expect_equal(mark$n_peaks, 150)
  expect_match(mark$recipe, "window = 1000L", fixed = TRUE)
})

test_that("unsupported peak requests fail loudly or go to the legacy router", {
  bad <- .simitall_agent_spec_pipeline("ChIP-seq for H3K9me2")
  expect_true(any(grepl("H3K9me2' is not simulated", bad$errors, fixed = TRUE)))
  expect_false(.simitall_agent_spec_pipeline("single-cell ATAC-seq of maize leaves")$supported)
  expect_equal(.simitall_agent_spec_pipeline("run a gwas for a phenotype in maize NAM")$task, "gwas")
})

test_that("peak studies ask about peak FDR and the differential test", {
  atac <- .simitall_agent_spec_pipeline("ATAC-seq under drought")
  expect_equal(vapply(.simitall_agent_spec_questions(atac), `[[`, "", "id"), c("peak_fdr", "diff_method"))
  tuned <- .simitall_agent_recompile_spec(atac, list(diff_method = "ttest", peak_fdr = "0.01", n_peaks = 50))
  expect_match(tuned$recipe, "method = \"ttest\"", fixed = TRUE)
  expect_match(tuned$recipe, "fdr = 0.01)", fixed = TRUE)
  expect_match(tuned$recipe, "n_peaks = 50L", fixed = TRUE)
  expect_match(.simitall_agent_compile_comparison(atac, "diff_method"), "differential_truth_recovery_metrics", fixed = TRUE)
})

test_that("simulated ATAC peaks and differential accessibility are recovered", {
  dir <- file.path(tempdir(), "atac_unit")
  genome <- .simitall_simulate_peak_genome(file.path(dir, "g.fa"), file.path(dir, "g.gff3"), seed = 3)
  sim <- suppressMessages(simulate_atacseq(genome$fasta, file.path(dir, "atac"), annotation_gff3 = genome$gff3,
                                           n_peaks = 60, conditions = c("control", "drought"), n_reads = 30000, seed = 4))
  reads <- read.delim(sim$read_positions)
  qc <- read.delim(sim$qc)
  expect_true(all(qc$frip > 0.2 & qc$frip < 0.45))
  expect_true(all(abs(qc$nucleosome_free_fraction - 0.5) < 0.1))
  called <- call_peaks_windowed(reads, genome$fasta, signal_assay = "atac")
  expect_gt(score_peak_calls(called, sim$peaks)$f1, 0.8)
  da <- test_differential_peaks(reads, called, assay = "atac", method = "ttest")
  scored <- score_differential_peaks(da, called, sim$peaks, sim$differential_truth)
  expect_gt(scored$precision, 0.7)
})

# Single-cell RNA-seq: plain-English prompts compile to plan and code.

test_that("single-cell prompts become single-cell specs with tissue-aware cell types", {
  leaf <- .simitall_agent_spec_pipeline("single-cell RNA-seq of maize leaves under drought")
  expect_equal(leaf$task, "single_cell"); expect_true(leaf$supported)
  expect_equal(names(leaf$cell_types), c("mesophyll", "epidermis", "vascular", "guard_cell"))
  expect_equal(leaf$conditions, c("control", "drought"))
  expect_match(leaf$recipe, "simulate_scrnaseq_from_gwas(", fixed = TRUE)
  root <- .simitall_agent_spec_pipeline("scRNA-seq of arabidopsis roots, 8 plants, 300 cells per plant, control vs salt")
  expect_equal(root$tissue, "root"); expect_equal(root$n_donors, 8); expect_equal(root$cells_per_donor, 300)
  expect_false(.simitall_agent_spec_pipeline("single-cell eQTL mapping in maize")$supported)
  expect_equal(.simitall_agent_spec_pipeline("single-cell ATAC-seq of maize leaves")$task, "peak_study")
})

test_that("single-cell questions and knobs recompile the recipe", {
  leaf <- .simitall_agent_spec_pipeline("single-cell RNA-seq of maize leaves under drought")
  expect_equal(vapply(.simitall_agent_spec_questions(leaf), `[[`, "", "id"), c("hvg_fraction", "sc_de_method"))
  tuned <- .simitall_agent_recompile_spec(leaf, list(sc_de_method = "cell_wilcoxon", hvg_fraction = "0.5", n_donors = 10))
  expect_match(tuned$recipe, "method = \"cell_wilcoxon\"", fixed = TRUE)
  expect_match(tuned$recipe, "n_hvg = round(0.5 * params$n_genes)", fixed = TRUE)
  expect_match(tuned$recipe, "n_donors = 10L", fixed = TRUE)
  expect_match(.simitall_agent_spec_plan(tuned), "inflated false positives", fixed = TRUE)
  expect_true(any(grepl("^n_donors", .simitall_agent_recompile_spec(leaf, list(n_donors = 2))$errors)))
})

test_that("adjusted Rand index ignores label names", {
  expect_equal(.simitall_adjusted_rand(c("a", "a", "b", "b"), c("x", "x", "y", "y")), 1)
  expect_lt(abs(.simitall_adjusted_rand(rep(1:2, 50), rep(1:2, each = 50))), 0.05)
})

test_that("clustering recovers simulated cell types and per-cell DE inflates false positives", {
  skip_if_not_installed("Matrix")
  dir <- file.path(tempdir(), "sc_unit")
  donors <- .simitall_simulate_inbred_panel(file.path(dir, "donors.vcf"), n_accessions = 6L, n_markers = 200L,
                                            n_chromosomes = 2L, n_subpops = 1L, fst = 0.05, seed = 2L)
  meta <- data.frame(sample = sprintf("acc%04d", 1:6), condition = rep(c("control", "drought"), 3), batch = "batch1")
  sim <- suppressMessages(simulate_scrnaseq_from_gwas(
    donors$vcf, file.path(dir, "sc"), sample_metadata = meta, n_genes = 250L, cells_per_donor = 120L,
    cell_type_proportions = c(mesophyll = 0.4, epidermis = 0.3, vascular = 0.3),
    markers_per_cell_type = 25L, marker_effect_mean = 2, n_cis_eqtl = 10L, n_trans_eqtl = 2L,
    write_sce = FALSE, backend = "native", seed = 3L))
  obj <- .simitall_read_scrna_object(file.path(dir, "sc"))
  clusters <- cluster_scrnaseq(obj, n_pcs = 10L, seed = 1L)
  expect_gt(score_scrnaseq_clusters(clusters, obj$cell_metadata)$adjusted_rand_index, 0.7)
  pseudo <- score_scrnaseq_conditions(test_scrnaseq_conditions(obj, "pseudobulk_ttest"), sim$condition_truth)
  percell <- score_scrnaseq_conditions(test_scrnaseq_conditions(obj, "cell_wilcoxon"), sim$condition_truth)
  expect_gt(percell$false_positives, pseudo$false_positives)
})

# Breeding populations: plain-English prompts compile to plan and code.

test_that("breeding prompts parse scheme and generation counts", {
  bc <- .simitall_agent_spec_pipeline("Simulate a synthetic biparental cross with 2 backcrosses and 7 selfing generations")
  expect_equal(bc$task, "breeding"); expect_equal(bc$scheme, "NIL")
  expect_equal(c(bc$backcrosses, bc$self_generations), c(2, 7))
  ril <- .simitall_agent_spec_pipeline("make 200 RILs by single seed descent to F8")
  expect_equal(ril$scheme, "RIL"); expect_equal(ril$self_generations, 7); expect_equal(ril$n_lines, 200)
  nil <- .simitall_agent_spec_pipeline("develop near-isogenic lines with BC3S2")
  expect_equal(c(nil$backcrosses, nil$self_generations), c(3, 2))
  magic <- .simitall_agent_spec_pipeline("8-way MAGIC population with 300 lines")
  expect_equal(magic$scheme, "MAGIC"); expect_equal(magic$n_founders, 8)
  expect_equal(.simitall_agent_spec_pipeline("doubled haploid population from two inbred parents")$scheme, "DH")
  expect_equal(.simitall_agent_spec_pipeline("find the genes for yield in corn")$task, "gwas")
})

test_that("Mendelian expectations are right", {
  expect_equal(breeding_expectations("NIL", 2, 0)$parent1_genome, 0.875)
  expect_equal(breeding_expectations("NIL", 2, 1)$heterozygosity, 0.125)
  expect_equal(breeding_expectations("RIL", 0, 6)$heterozygosity, 1 / 64)
  expect_equal(breeding_expectations("DH")$heterozygosity, 0)
  expect_equal(breeding_expectations("MAGIC", n_founders = 8)$parent1_genome, 0.125)
  expect_true(is.na(breeding_expectations("RIL", 0, 6, ril_mating = "SIB")$heterozygosity))
})

test_that("breeding questions follow the scheme and recompile the recipe", {
  nil <- .simitall_agent_spec_pipeline("develop near-isogenic lines with BC3S2")
  expect_equal(vapply(.simitall_agent_spec_questions(nil), `[[`, "", "id"), c("backcrosses", "selfing"))
  tuned <- .simitall_agent_recompile_spec(nil, list(backcrosses = "4", n_lines = 80))
  expect_match(tuned$recipe, "backcrosses = 4L", fixed = TRUE)
  expect_match(tuned$recipe, "n_lines = 80L", fixed = TRUE)
  expect_match(.simitall_agent_spec_plan(tuned), "BC4S2", fixed = TRUE)
  expect_length(.simitall_agent_spec_questions(.simitall_agent_spec_pipeline("simulate an F2 population")), 0L)
  expect_match(.simitall_agent_compile_comparison(nil, "backcrosses"), "breeding_truth_recovery_metrics", fixed = TRUE)
})

test_that("simulated backcross and MAGIC lines match Mendelian theory", {
  dir <- file.path(tempdir(), "breeding_unit"); dir.create(dir, showWarnings = FALSE)
  panel <- file.path(dir, "founders.fa")
  invisible(capture.output(generate_random_haplotype_panel(out_fa = panel, n_haplotypes = 4L, length = 2000L, n_chromosomes = 3L,
                                                          snp_rate = 0.02, indel_rate = 0, seed = 1L)))
  map <- file.path(dir, "map.tsv")
  write.table(data.frame(chromosome = rep(paste0("chr", 1:3), each = 2), pos_bp = rep(c(1L, 2000L), 3), cM = rep(c(0, 100), 3)),
              map, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(capture.output(simulate_breeding(haplotype_fa = panel, out_prefix = file.path(dir, "nil"), scheme = "NIL",
                                             n_offspring = 80, parents = "hap1,hap2", backcross_generations = 2,
                                             self_generations = 1, recomb_map_in = map, seed = 2)))
  nil <- summarize_breeding_population(file.path(dir, "nil.ancestry.tsv"), "NIL", "hap1", backcrosses = 2, self_generations = 1)
  expect_true(nil$metrics$agrees_with_theory)
  # Regression test: each MAGIC line must come from its own funnel, so every
  # founder contributes about 1/n_founders on average.
  invisible(capture.output(simulate_breeding(haplotype_fa = panel, out_prefix = file.path(dir, "magic"), scheme = "MAGIC",
                                             n_offspring = 80, founders = "hap1,hap2,hap3,hap4", self_generations = 2,
                                             recomb_map_in = map, seed = 3)))
  magic <- summarize_breeding_population(file.path(dir, "magic.ancestry.tsv"), "MAGIC", "hap1", self_generations = 2, n_founders = 4)
  expect_lt(magic$metrics$founder_balance, 0.05)
  expect_true(magic$metrics$agrees_with_theory)
})
