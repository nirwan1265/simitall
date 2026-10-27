test_that("capability inventory is packaged and references real resources", {
  registry <- .simitall_agent_capability_registry()
  expect_identical(registry$planner_integration, "not_connected")
  entries <- registry$workflows
  ids <- vapply(entries, `[[`, character(1), "id")
  expect_equal(anyDuplicated(ids), 0L)
  expect_true(length(ids) >= 15L)
  exports <- getNamespaceExports("simitall")
  for (entry in entries) {
    expect_true(entry$implementation_status %in% c("partial", "planned"), info = entry$id)
    expect_true(length(entry$required_inputs) > 0L, info = entry$id)
    expect_true(all(unlist(entry$package_functions) %in% exports), info = entry$id)
    for (ref in c(entry$knowledge_sources, entry$paper_sources)) {
      expect_true(nzchar(system.file("agent", "knowledge", ref, package = "simitall")), info = ref)
    }
  }
})

test_that("registry preserves important scientific and support distinctions", {
  entries <- .simitall_agent_capability_registry()$workflows
  by_id <- setNames(entries, vapply(entries, `[[`, character(1), "id"))
  roh <- by_id$human_multigeneration_roh_gene_ranking
  expect_identical(roh$implementation_status, "planned")
  expect_true(any(grepl("phased founder", unlist(roh$required_inputs), fixed = TRUE)))
  expect_true(any(grepl("recombination map", unlist(roh$required_inputs), fixed = TRUE)))
  counts <- unlist(by_id$bulk_rnaseq_counts$required_inputs)
  expect_false(any(grepl("FASTA", counts, ignore.case = TRUE)))
  expect_true(any(grepl("transcript sequences", unlist(by_id$rnaseq_read_simulation$required_inputs), fixed = TRUE)))
  expect_identical(by_id$qtl_mapping$implementation_status, "planned")
})

test_that("analysis catalog separates concepts, goals, models, and inputs", {
  catalog <- .simitall_agent_analysis_catalog()
  families <- catalog$families
  ids <- vapply(families, `[[`, character(1), "id")
  expect_equal(anyDuplicated(ids), 0L)
  expect_true(length(families) >= 12L)
  subanalyses <- unlist(lapply(families, `[[`, "subanalyses"), recursive = FALSE)
  sub_ids <- vapply(subanalyses, `[[`, character(1), "id")
  expect_equal(anyDuplicated(sub_ids), 0L)
  expect_true(all(vapply(subanalyses, function(x) length(x$candidate_models) > 0L, logical(1))))
  expect_true(all(vapply(subanalyses, function(x) length(x$required_inputs) > 0L, logical(1))))
  expect_true(all(vapply(subanalyses, function(x) x$execution_status %in% c("partial", "planned", "knowledge_catalog_only"), logical(1))))
  by_id <- setNames(subanalyses, sub_ids)
  expect_true("maSigPro regression" %in% by_id$time_series_rnaseq$candidate_models)
  expect_true("differential_exon_usage" %in% by_id$alternative_splicing$parameters[[1L]]$allowed_values)
  expect_true("cis" %in% by_id$bulk_eqtl$parameters[[1L]]$allowed_values)
  expect_true(any(grepl("Do not conflate pedigree F", by_id$pedigree_inbreeding$qc_and_guardrails, fixed = TRUE)))
})

test_that("input and default registries keep data contracts explicit", {
  inputs <- .simitall_agent_input_type_registry()$types
  ids <- vapply(inputs, `[[`, character(1), "id")
  expect_true(all(c("vcf_bcf", "gene_annotation", "count_matrix", "sample_metadata", "pedigree") %in% ids))
  expect_true(all(vapply(inputs, function(x) length(x$validation) > 0L, logical(1))))
  defaults <- .simitall_agent_default_policy()
  expect_true(any(grepl("mating scheme", defaults$must_ask_or_require, ignore.case = TRUE)))
  expect_true(any(grepl("binary phenotype", defaults$never_infer, ignore.case = TRUE)))
})

test_that("sample naming contract keeps treatment and time unambiguous", {
  contract <- .simitall_agent_sample_naming_contract()
  expect_identical(contract$canonical_sample_id, "S001_R01_TrtControl_Ti00")
  expect_true(all(c("sample_id", "biological_replicate") %in% contract$metadata_columns$required))
  expect_true(any(grepl("bare `T`", contract$rules, fixed = TRUE)))
  expect_match(.simitall_agent_sample_naming_text(), "Trt = treatment, Ti = time", fixed = TRUE)
})

test_that("tutorial maps and reviewed examples are well formed", {
  map_path <- .simitall_agent_capability_path("tutorial_map.json")
  tutorial_map <- jsonlite::fromJSON(map_path, simplifyVector = FALSE)
  topics <- tutorial_map$topics
  expect_true(length(topics) >= 10L)
  references <- unlist(lapply(topics, function(x) unlist(x[c("tutorials", "workflows", "validation_contracts", "tools", "papers")])), use.names = FALSE)
  references <- references[!is.na(references)]
  for (ref in references) {
    expect_true(nzchar(system.file("agent", "knowledge", ref, package = "simitall")), info = ref)
  }
  examples_path <- .simitall_agent_capability_path("study_spec_examples.jsonl")
  examples <- lapply(readLines(examples_path, warn = FALSE), jsonlite::fromJSON, simplifyVector = FALSE)
  expect_true(length(examples) >= 10L)
  expect_true(all(vapply(examples, function(x) identical(x$review_status, "reviewed"), logical(1))))
  expect_true(all(vapply(examples, function(x) isTRUE(x$expected_spec$do_not_run), logical(1))))
  expect_true(all(vapply(examples, function(x) length(x$expected_spec$unresolved_questions) > 0L, logical(1))))
})

test_that("generic catalog planning preserves a requested RNA-seq analysis", {
  question <- paste(
    "Simulate bulk RNA-seq counts for control and drought, 4 biological",
    "replicates each, with batch effects and a negative-binomial model."
  )
  resolved <- .simitall_agent_catalog_resolve(question, data_source = "synthetic")
  expect_identical(resolved$status, "recognized")
  expect_identical(resolved$primary$subanalysis$id, "bulk_counts_simulation")
  expect_length(resolved$secondary, 0L)
  parameter_values <- setNames(
    lapply(resolved$parameters, `[[`, "value"),
    vapply(resolved$parameters, `[[`, character(1), "id")
  )
  expect_identical(parameter_values$count_model, "negative_binomial")
  expect_identical(parameter_values$replicates, 4L)

  answer <- simitall_ask(
    question,
    provider = "local",
    data_source = "synthetic",
    include_code = TRUE
  )
  expect_identical(answer$plan_source, "catalog")
  expect_identical(answer$catalog_study_spec$primary_analysis, "bulk_counts_simulation")
  expect_true(any(grepl("S001_R01_TrtControl_Ti00", answer$catalog_study_spec$assumptions, fixed = TRUE)))
  expect_false(grepl("eQTL", answer$answer, fixed = TRUE))
  expect_false(grepl("```", answer$answer, fixed = TRUE))
})

test_that("a broad catalog family asks for a goal instead of guessing", {
  resolved <- .simitall_agent_catalog_resolve(
    "I need help planning a bulk RNA-seq study.",
    data_source = "provided"
  )
  expect_identical(resolved$status, "needs_goal")
  expect_identical(resolved$family$id, "bulk_rnaseq")
})

test_that("an ambiguous RNA-seq request receives a partial StudySpec question", {
  answer <- simitall_ask("I want to do RNA-seq", provider = "local", data_source = "provided")
  expect_identical(answer$plan_source, "catalog")
  expect_identical(answer$catalog_resolution$status, "needs_goal")
  expect_identical(answer$catalog_study_spec$execution_status, "needs_clarification")
  expect_null(answer$catalog_study_spec$primary_analysis)
  expect_match(answer$answer, "Choose the intended goal", fixed = TRUE)
  expect_false(grepl("eQTL", answer$answer, fixed = TRUE))
})

test_that("catalog resolution treats supplied file formats as inputs, not analyses", {
  splicing <- .simitall_agent_catalog_resolve(
    "Compare exon usage between treated and control RNA-seq samples using my BAM files and GTF.",
    data_source = "provided"
  )
  expect_identical(splicing$status, "recognized")
  expect_identical(splicing$primary$subanalysis$id, "alternative_splicing")
  expect_identical(splicing$recognized_family_ids, "bulk_rnaseq")

  eqtl <- .simitall_agent_catalog_resolve(
    "Map cis eQTLs from a VCF and normalized gene expression for the same 300 donors.",
    data_source = "provided"
  )
  expect_identical(eqtl$status, "recognized")
  expect_identical(eqtl$primary$subanalysis$id, "bulk_eqtl")
})
