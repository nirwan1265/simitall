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
  expect_equal(preview$execution, "disabled; planning and question answering only")
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
  expect_match(answer$answer, "Verified NAM execution recipe", fixed = TRUE)
  expect_match(answer$answer, "simulate_breeding", fixed = TRUE)
  expect_match(answer$answer, "simulate_rnaseq_from_gwas", fixed = TRUE)
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
