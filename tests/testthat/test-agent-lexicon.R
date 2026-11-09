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
