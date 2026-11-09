test_that("lexicon tables are consistent with each other and the catalog", {
  lex <- .simitall_lexicon()
  expect_gte(nrow(lex$organisms), 80)
  expect_false(anyDuplicated(lex$organisms$organism) > 0)
  expect_true(all(lex$organism_aliases$organism %in% lex$organisms$organism))
  expect_true(all(lex$organisms$can_self %in% c("yes", "no", "conditional")))
  expect_true(all(grepl("^[0-9]+$", lex$organisms$taxid)))
  expect_true(all(lex$technique_aliases$technique %in% lex$techniques$technique))
  expect_true(all(lex$organism_aliases$strength %in% c("strong", "weak")))

  cat_path <- file.path(test_path("..", ".."), "inst", "agent", "capabilities", "analysis_catalog.json")
  if (!file.exists(cat_path)) cat_path <- system.file("agent", "capabilities", "analysis_catalog.json", package = "simitall")
  catalog <- jsonlite::fromJSON(cat_path, simplifyVector = FALSE)
  sub_ids <- vapply(unlist(lapply(catalog$families, function(f) f$subanalyses), recursive = FALSE), `[[`, "", "id")
  routes <- lex$techniques$routes_to[nzchar(lex$techniques$routes_to)]
  cands <- unlist(strsplit(lex$techniques$candidates, ";", fixed = TRUE))
  expect_true(all(c(routes, cands) %in% sub_ids))
  # after normalisation every alias still points to exactly one target
  for (a in list(lex$organism_aliases, lex$technique_aliases)) {
    n_targets <- tapply(a$target, paste(a$key, a$case_sensitive), function(x) length(unique(x)))
    expect_true(all(n_targets == 1))
  }
})

test_that("organisms are recognised from common, Latin, strain, panel and breed names", {
  org <- function(q) .simitall_detect_organism(q)$organism
  expect_equal(org("GWAS of flowering in corn"), "maize")
  expect_equal(org("my Zea mays panel"), "maize")
  expect_equal(org("B73 x Mo17 RILs"), "maize")
  expect_equal(org("200 DGRP lines"), "drosophila")
  expect_equal(org("D. melanogaster heads"), "drosophila")
  expect_equal(org("Collaborative Cross mice"), "mouse")
  expect_equal(org("heterogeneous stock rats"), "rat")
  expect_equal(org("C. elegans N2 hermaphrodites"), "c_elegans")
  expect_equal(org("3,000 Holstein cows"), "cattle")
  expect_equal(org("Labradors with hip dysplasia"), "dog")
  expect_equal(org("E. coli K-12 reads"), "ecoli")
  expect_equal(org("my GRCh38 VCF"), "human")
  expect_equal(org("fission yeast tetrads"), "fission_yeast")
  expect_equal(org("durum wheat lines"), "durum_wheat")
  expect_equal(org("sick and healthy chickens, bacteria in the gut"), "chicken")
  # no false hits inside other words
  expect_true(is.na(org("what is the price of sequencing?")))
  expect_true(is.na(org("binomial phenotype GWAS in NAM")))
  expect_true(is.na(org("run gwas")))
})

test_that("techniques route to catalog analyses or stop when unsupported", {
  tq <- function(q) .simitall_detect_techniques(q)
  expect_true(tq("RNA-seq.")$needs_goal)
  expect_true(is.na(tq("RNA-seq.")$primary_analysis))
  expect_equal(tq("DESeq2 on my RNA-seq counts")$primary_analysis, "differential_expression")
  expect_equal(tq("ATAC-seq peaks, drought vs control")$primary_analysis, "chipseq_atacseq")
  expect_equal(tq("ATAC-seq peaks, drought vs control")$parameters$assay, "atacseq")
  expect_equal(tq("simulate NovaSeq X reads at 30x")$parameters$platform, "NovaSeq_X")
  expect_equal(tq("case-control GWAS, 4000 cases")$primary_analysis, "binary_gwas")
  expect_equal(tq("GWAS of plant height")$primary_analysis, "quantitative_gwas")
  expect_equal(tq("map cis-eQTLs in 300 donors")$primary_analysis, "bulk_eqtl")
  for (q in c("Ribo-seq translation efficiency", "Call TADs from Hi-C", "CUT&Tag for H3K27me3",
              "Visium spatial transcriptomics", "TMT proteomics", "genome-wide CRISPR knockout screen",
              "DMRs from bisulfite sequencing", "PacBio Iso-Seq isoforms")) {
    expect_equal(tq(q)$execution_status, "unsupported", info = q)
  }
})

test_that("biology checks catch impossible or contradictory designs", {
  expect_match(.simitall_lexicon_checks("selfing in dogs 5 generations"), "selfing is not possible")
  expect_length(.simitall_lexicon_checks("selfing maize for 7 generations"), 0)
  expect_length(.simitall_lexicon_checks("C. elegans hermaphrodites self for 10 generations"), 0)
  expect_match(.simitall_lexicon_checks("F2 in tetraploid alfalfa with diploid segregation ratios"), "diploid assumptions conflict")
  expect_match(.simitall_lexicon_checks("a wolf bottleneck 30 years ago"), "convert to generations")
})

test_that("settings, quantities and intent tables are valid", {
  lex <- .simitall_lexicon()
  cat_path <- file.path(test_path("..", ".."), "inst", "agent", "capabilities", "analysis_catalog.json")
  if (!file.exists(cat_path)) cat_path <- system.file("agent", "capabilities", "analysis_catalog.json", package = "simitall")
  catalog <- jsonlite::fromJSON(cat_path, simplifyVector = FALSE)
  subs <- unlist(lapply(catalog$families, function(f) f$subanalyses), recursive = FALSE)
  names(subs) <- vapply(subs, `[[`, "", "id")
  expect_true(all(lex$settings$analysis %in% names(subs)))
  for (i in seq_len(nrow(lex$settings))) {
    s <- lex$settings[i, ]
    pid <- sub("^(parameters|design)\\.", "", s$field)
    p <- Filter(function(p) identical(p$id, pid), subs[[s$analysis]]$parameters %||% list())
    expect_length(p, 1)
    if (length(p) && length(p[[1]]$allowed_values)) expect_true(s$value %in% unlist(p[[1]]$allowed_values), info = s$alias)
  }
  for (p in c(lex$quantities$pattern, lex$intent$cue)) expect_silent(grepl(p, "test", perl = TRUE))
  expect_true(all(grepl("^(design|parameters)\\.|^region$", lex$quantities$field)))
  expect_true(all(lex$intent_routes$analysis %in% names(subs)))
})

test_that("the lexicon reads a prompt into a StudySpec view", {
  v <- .simitall_lexicon_read_prompt("Simulate 200 RILs from two founders with eight generations of selfing.")
  expect_equal(v$goal, "simulate")
  expect_equal(v$data_source, "synthetic")
  expect_equal(v$primary_analysis, "cross_design")
  expect_equal(v$design[c("scheme", "founder_count", "selfing_generations", "population_size")],
               list(scheme = "RIL", founder_count = 2L, selfing_generations = 8L, population_size = 200L))

  v <- .simitall_lexicon_read_prompt("GWAS with an LMM, kinship plus 3 PCs, FDR 5%, on my rice panel of 413 accessions.")
  expect_equal(v[c("goal", "data_source", "organism", "primary_analysis")],
               list(goal = "analyze", data_source = "provided", organism = "rice", primary_analysis = "quantitative_gwas"))
  expect_setequal(v$parameters$structure_adjustment, c("PCs", "kinship"))
  expect_equal(v$parameters$n_PCs, 3L)
  expect_length(v$secondary_analyses, 0)

  v <- .simitall_lexicon_read_prompt("Simulate scRNA-seq with three cell types, 12 donors, control and treatment, then run pseudobulk DE.")
  expect_equal(v$goal, "simulate_and_analyze")
  expect_equal(v$primary_analysis, "scrnaseq_simulation")
  expect_equal(v$secondary_analyses, "scrnaseq_analysis")

  v <- .simitall_lexicon_read_prompt("Only simulate the RNA-seq counts, I'll do the analysis myself.")
  expect_equal(v$primary_analysis, "bulk_counts_simulation")
  expect_length(v$secondary_analyses, 0)

  v <- .simitall_lexicon_read_prompt("Could you simulate a drought time course in maize, 0, 2, 6, 24 and 48 hours, 3 reps each?")
  expect_equal(v$design$time_points, 5L)

  expect_true(.simitall_lexicon_read_prompt("RNA-seq.")$data_source == "ask")
  expect_null(.simitall_lexicon_read_prompt("RNA-seq.")$primary_analysis)
  # input roles and exclusions
  expect_equal(.simitall_lexicon_read_prompt("which genes are within 10 kb of my GWAS hits? rice")$primary_analysis,
               "gene_coordinate_annotation")
  expect_false("bulk_eqtl" %in% unlist(.simitall_lexicon_read_prompt(
    "Simulate RNA-seq counts for my cohort, no eQTL mapping.")[c("primary_analysis", "secondary_analyses")]))
})

test_that("guard rules stop impossible, clinical, real-data and corrected requests", {
  stops <- function(q) .simitall_lexicon_read_prompt(q)$do_not_run
  expect_true(stops("Simulate a trait with heritability 1.4 on 200 lines."))
  expect_null(.simitall_lexicon_read_prompt("Simulate a trait with heritability 1.4 on 200 lines.")$parameters$heritability)
  expect_true(stops("simulate reads at -5x coverage"))
  expect_true(stops("Simulate scRNA-seq where 120% of eQTLs are cell-type specific."))
  expect_true(stops("Use my patient's RNA-seq to decide which cancer drug to give."))
  expect_true(stops("Simulate the real maize NAM population with the actual 25 founders."))
  expect_true(stops("selfing in dogs 5 generations"))
  expect_true(stops("Simulate Ribo-seq footprints for heat vs control with matched RNA-seq."))
  expect_false(stops("Simulate 10 generations of selfing in a population of 500 Arabidopsis plants."))
  expect_null(.simitall_lexicon_read_prompt("Simulate 200 RILs, actually make it 300, from two founders.")$design$population_size)
})
