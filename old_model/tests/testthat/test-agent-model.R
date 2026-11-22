test_that("model merge only fills gaps, keeps guards, and never adds numbers to a stopped request", {
  q <- "Can you simulate scRNA-seq for 6 healthy and 6 sick people, 4 cell types, 300 cells per donor?"
  pred <- list(spans = list(list("DONORS", "6", 5L), list("COND", "healthy", 6L), list("DONORS", "6", 8L),
                            list("COND", "sick", 9L), list("NCELLTYPES", "4", 12L), list("CELLS", "300", 16L)),
               goal = list("simulate", 0.99), data_source = list("synthetic", 0.99),
               primary = list("scrnaseq_simulation", 0.95))
  lex <- .simitall_lexicon_read_prompt(q)
  v <- .simitall_model_merge(lex, pred, q)
  expect_equal(v$design$donors, 12L)                       # "6 healthy and 6 sick" summed
  expect_equal(v$design$conditions, list("healthy", "sick"))
  expect_equal(v$parameters$cells_per_donor, lex$parameters$cells_per_donor)  # lexicon value kept
  expect_equal(v$primary_analysis, "scrnaseq_simulation")

  stop_q <- "simulate reads at -5x coverage"
  stopped <- .simitall_lexicon_read_prompt(stop_q)
  v2 <- .simitall_model_merge(stopped, list(spans = list(list("COV", "5 x", 4L)), primary = list("read_simulation", 1)), stop_q)
  expect_null(v2$parameters$coverage)

  r <- "Simulate 10 rounds of first-cousin marriages in humans."
  v3 <- .simitall_model_merge(.simitall_lexicon_read_prompt(r), list(spans = list(list("GENS", "10", 1L))), r)
  expect_null(v3$parameters$generations)                   # rounds are not generations

  s <- "Simulate 10 generations of sib mating in rats, then a GWAS of a recessive trait."
  v4 <- .simitall_model_merge(.simitall_lexicon_read_prompt(s), list(spans = list(list("PRIMARY", "GWAS", 9L))), s)
  expect_equal(v4$primary_analysis, "pedigree_inbreeding")  # the simulation step stays first
})

test_that("model merge ignores chromosome numbers, unmarked inputs and keeps 'X, then Y' order", {
  q <- "Simulate a synthetic human chr10 pedigree with first-cousin marriages for 10 generations and compute pedigree F."
  pred <- list(spans = list(list("ORG", "human", 3L), list("POP", "10", 5L), list("GENS", "10", 12L)))
  v <- .simitall_model_merge(.simitall_lexicon_read_prompt(q), pred, q)
  expect_null(v$parameters$population_size)

  g <- "GWAS on my Arabidopsis VCF for flowering time, then annotate significant SNPs with nearby genes and GO terms."
  pred <- list(spans = list(list("INPUT", "GWAS", 0L), list("ORG", "Arabidopsis", 3L), list("PRIMARY", "GO terms", 19L)))
  v <- .simitall_model_merge(.simitall_lexicon_read_prompt(g), pred, g)
  expect_equal(v$primary_analysis, "quantitative_gwas")
  expect_true("gene_coordinate_annotation" %in% v$secondary_analyses)
})

test_that("the reader falls back to the lexicon without a model", {
  old <- Sys.getenv("SIMITALL_PROMPT_MODEL")
  Sys.setenv(SIMITALL_PROMPT_MODEL = tempfile())
  on.exit(Sys.setenv(SIMITALL_PROMPT_MODEL = old), add = TRUE)
  v <- .simitall_read_prompt("Simulate 200 RILs from two founders.")
  expect_equal(v$reader, "lexicon")
  expect_equal(v$primary_analysis, "cross_design")
})
