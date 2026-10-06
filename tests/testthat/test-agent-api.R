test_that("local knowledge search ranks matching documents", {
  knowledge <- file.path(tempdir(), "simitall-agent-knowledge")
  dir.create(file.path(knowledge, "workflows"), recursive = TRUE)
  dir.create(file.path(knowledge, "species"), recursive = TRUE)
  writeLines(c("# Maize NAM workflow", "Use a NAM population and recombination map."),
             file.path(knowledge, "workflows", "nam.Rmd"))
  writeLines(c("# Human notes", "Human reference genome advice."),
             file.path(knowledge, "species", "human.Rmd"))

  found <- search_simitall_knowledge(
    "maize NAM recombination",
    n_results = 1,
    knowledge_dir = knowledge
  )
  expect_equal(found$title, "Maize NAM workflow")
  expect_equal(found$category, "workflows")
  expect_gt(found$score, 0)
})

test_that("versioned alias table normalizes biology terminology", {
  aliases <- simitall:::.simitall_agent_aliases()
  expect_gte(nrow(aliases), 70L)
  terms <- simitall:::.simitall_agent_terms(
    paste(
      "Use epistatic, polygenic effects in an admixed genome-wide association study",
      "with DEG transcript abundance from zero-inflated pseudo-bulk scRNA-seq",
      "and long-read HiFi CCS CLR data."
    )
  )
  expect_true(all(c(
    "epistasis", "polygenicity", "admixture", "gwas", "differential expression",
    "transcript abundance", "zero inflation", "pseudobulk", "scrnaseq",
    "long read sequencing", "high fidelity reads", "circular consensus sequencing",
    "continuous long reads"
  ) %in% terms))
})

test_that("figure terminology retrieves validation guidance", {
  found <- search_simitall_knowledge(
    "Which figures validate LD, breeding, GWAS, RNA-seq, and ChIP-seq?",
    n_results = 3
  )
  expect_true(any(basename(found$path) == "validation_figures.Rmd"))
})

test_that("validation contracts retrieve the four-part GWAS specification", {
  found <- search_simitall_knowledge(
    "GWAS required figure simulation truth results pass criteria",
    n_results = 3
  )
  expect_true(any(basename(found$path) == "phenotypes_gwas.Rmd"))
  expect_true(any(found$category == "validation_contracts"))

  contract_path <- found$path[basename(found$path) == "phenotypes_gwas.Rmd"][1L]
  contract <- readLines(contract_path, warn = FALSE)
  expect_true(all(c(
    "## Simulation Truth", "## Required Results", "## Required Figures", "## Pass Criteria"
  ) %in% contract))
})

test_that("retrieval prioritizes workflow knowledge over repeated paper text", {
  knowledge <- file.path(tempdir(), "simitall-agent-priority")
  dir.create(file.path(knowledge, "workflows"), recursive = TRUE)
  dir.create(file.path(knowledge, "papers"), recursive = TRUE)
  writeLines(
    c("# Maize NAM workflow", "Plan a maize NAM population for GWAS."),
    file.path(knowledge, "workflows", "maize_nam.Rmd")
  )
  writeLines(
    c("# Maize review", rep("maize NAM GWAS", 100L)),
    file.path(knowledge, "papers", "review.Rmd")
  )

  found <- search_simitall_knowledge(
    "maize NAM GWAS", n_results = 2, knowledge_dir = knowledge
  )
  expect_equal(found$category[1L], "workflows")
  expect_equal(found$title[1L], "Maize NAM workflow")
})

test_that("agent request previews are local and non-executing", {
  knowledge <- file.path(tempdir(), "simitall-agent-preview")
  dir.create(file.path(knowledge, "tools"), recursive = TRUE)
  writeLines(c("# GWAS", "GWAS requires genotypes and phenotypes."),
             file.path(knowledge, "tools", "gwas.Rmd"))
  preview <- simitall_ask(
    "What does GWAS require?",
    knowledge_dir = knowledge,
    dry_run = TRUE
  )
  expect_s3_class(preview, "simitall_agent_answer")
  expect_null(preview$answer)
  expect_match(preview$instructions, "cannot execute R")
  expect_match(preview$instructions, "Do not write R code", fixed = TRUE)
  expect_match(preview$instructions, "simulate_gwas_cohort", fixed = TRUE)
  # simitall_ask() never executes code itself; only the user's RUN does.
  expect_match(preview$execution, "only when the user presses RUN", fixed = TRUE)
})

test_that("unsafe local-model output falls back to the evidence briefing", {
  knowledge <- file.path(tempdir(), "simitall-agent-safety")
  dir.create(file.path(knowledge, "workflows"), recursive = TRUE)
  writeLines(c("# GWAS workflow", "Retain causal truth."),
             file.path(knowledge, "workflows", "gwas.Rmd"))
  safe <- simitall:::.simitall_agent_validate_answer(
    "```r\nlibrary(GAPIT)\n```",
    "Plan a GWAS simulation",
    search_simitall_knowledge("GWAS", knowledge_dir = knowledge),
    allow_code = FALSE
  )
  expect_match(safe, "Safety fallback", fixed = TRUE)
  expect_match(safe, "Offline simitall plan", fixed = TRUE)
})

test_that("agent tool registry is explicit and read-only", {
  tools <- simitall_agent_tools()
  expect_true("simulate_gwas_cohort" %in% tools$function_name)
  expect_true(all(grepl("approval", tools$execution_status, fixed = TRUE)))
})

test_that("NAM requests receive a deterministic code recipe only on request", {
  knowledge <- file.path(tempdir(), "simitall-agent-nam-recipe")
  dir.create(file.path(knowledge, "workflows"), recursive = TRUE)
  writeLines(c("# NAM workflow", "Use founders and a recombination map."),
             file.path(knowledge, "workflows", "nam.Rmd"))
  answer <- simitall_ask(
    "Plan a maize NAM population for GWAS",
    provider = "local",
    include_code = TRUE,
    data_source = "synthetic",
    knowledge_dir = knowledge
  )
  expect_match(answer$answer, "Spec-compiled NAM quantitative-trait GWAS recipe", fixed = TRUE)
  expect_equal(answer$spec$resolution$mode, "bundled")
  expect_match(answer$answer, "simulate_breeding", fixed = TRUE)
  expect_match(answer$answer, "fig15_maize_nam_quantitative_gwas.R", fixed = TRUE)
  expect_false(grepl("simulate_rnaseq_from_gwas", answer$answer, fixed = TRUE))
})

test_that("synthetic biparental BC2S7 requests compile to the breeding workflow", {
  answer <- simitall_ask(
    "Simulate a synthetic biparental cross with 2 backcrosses and 7 selfing generations, then show allele frequencies.",
    provider = "local",
    include_code = TRUE,
    data_source = "synthetic"
  )
  expect_equal(answer$plan_source, "spec")
  expect_match(answer$answer, "Spec-compiled breeding-population recipe", fixed = TRUE)
  expect_match(answer$answer, "BC2S7", fixed = TRUE)
  expect_match(answer$answer, "backcrosses = 2L", fixed = TRUE)
  expect_match(answer$answer, "self_generations = 7L", fixed = TRUE)
})

test_that("provided IRF6 requests compile against the supplied VCF", {
  vcf <- system.file("extdata", "human_irf6", "demo_human_irf6_chr1.vcf.gz", package = "simitall")
  recipe <- simitall:::.simitall_agent_code_recipe(
    "Plan a human IRF6 pedigree liability simulation",
    list(status = "ready", data_source = "provided", input_files = vcf)
  )
  expect_match(recipe, "simulate_human_pedigree_groups", fixed = TRUE)
  expect_match(recipe, normalizePath(vcf), fixed = TRUE)
  expect_false(grepl("fig11_human_irf6_synthetic_liability.R", recipe, fixed = TRUE))
})

test_that("Arabidopsis 1001 requests require a real panel instead of a toy substitute", {
  answer <- simitall_ask(
    "Plan an Arabidopsis 1001 Genomes breeding population for GWAS",
    provider = "local",
    include_code = TRUE
  )
  expect_equal(answer$preflight$status, "needs_input")
  expect_match(answer$answer, "does not bundle an Arabidopsis 1001 Genomes haplotype panel", fixed = TRUE)
  expect_false(grepl("Verified NAM execution recipe", answer$answer, fixed = TRUE))
})

test_that("local agent answers without an API key or network request", {
  knowledge <- file.path(tempdir(), "simitall-local-agent")
  dir.create(file.path(knowledge, "workflows"), recursive = TRUE)
  writeLines(c("# GWAS workflow", "Retain phenotype, genotype, and causal truth files."),
             file.path(knowledge, "workflows", "gwas.Rmd"))
  answer <- simitall_ask(
    "Which GWAS truth files should I keep?",
    provider = "local",
    knowledge_dir = knowledge
  )
  expect_match(answer$answer, "Offline simitall plan", fixed = TRUE)
  expect_match(answer$answer, "Proposed stages", fixed = TRUE)
})

test_that("agent context accepts multi-term retrieval metadata", {
  knowledge <- file.path(tempdir(), "simitall-agent-multi-term")
  dir.create(file.path(knowledge, "workflows"), recursive = TRUE)
  writeLines(
    c("# Maize GWAS workflow", "Retain founder, family, and phenotype metadata."),
    file.path(knowledge, "workflows", "maize_gwas.Rmd")
  )
  answer <- simitall_ask(
    "Plan a maize GWAS with founder family metadata",
    provider = "local",
    knowledge_dir = knowledge
  )
  expect_match(answer$answer, "Offline simitall plan", fixed = TRUE)
})

test_that("local ancestry questions do not fall back to ordinary GWAS claims", {
  answer <- simitall_ask(
    "Plan a human admixed local ancestry GWAS",
    provider = "local"
  )
  expect_match(answer$answer, "phased cohort haplotypes", fixed = TRUE)
  expect_match(answer$answer, "does not itself phase data", fixed = TRUE)
  expect_false(grepl("simulate_gwas_cohort()", answer$answer, fixed = TRUE))
})

test_that("ancestry tool questions name the documented tool boundary", {
  answer <- simitall_ask(
    "Which packages should I use for RFMix Tractor GENESIS local ancestry GWAS?",
    provider = "local"
  )
  expect_match(answer$answer, "GenomeAdmixR", fixed = TRUE)
  expect_match(answer$answer, "RFMix", fixed = TRUE)
  expect_match(answer$answer, "Tractor", fixed = TRUE)
})

test_that("E. coli GWAS requests use the bacterial analysis boundary", {
  answer <- simitall_ask(
    "Plan an E. coli bacterial GWAS and hybrid assembly benchmark",
    provider = "local"
  )
  expect_match(answer$answer, "bacterial-GWAS", fixed = TRUE)
  expect_false(grepl("simulate_gwas_cohort()", answer$answer, fixed = TRUE))
})

test_that("human IRF6 pedigree prompts compile reusable modules instead of a figure script", {
  question <- paste(
    "Simulate a synthetic human family study around the IRF6 region.",
    "Compare unrelated people, distantly related people, and first-cousin descendants.",
    "Create an illustrative inherited-trait probability and show marker homozygosity."
  )
  answer <- simitall_ask(question, provider = "local", include_code = TRUE, data_source = "package_demo")
  expect_match(answer$answer, "simulate_human_pedigree_groups", fixed = TRUE)
  expect_match(answer$answer, "plot_human_pedigree_diagnostics", fixed = TRUE)
  expect_false(grepl("fig11_human_irf6_synthetic_liability.R", answer$answer, fixed = TRUE))
})
