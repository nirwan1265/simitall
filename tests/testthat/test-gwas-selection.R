make_genomic_fixture <- function() {
  set.seed(901)
  n_samples <- 48L
  n_markers <- 90L
  samples <- paste0("line", seq_len(n_samples))
  frequencies <- stats::runif(n_markers, 0.15, 0.85)
  genotype <- matrix(
    stats::rbinom(
      n_markers * n_samples,
      size = 2L,
      prob = rep(frequencies, n_samples)
    ),
    nrow = n_markers,
    ncol = n_samples
  )
  genotype[1L, ] <- rep(0:2, length.out = n_samples)
  colnames(genotype) <- samples
  variants <- data.frame(
    id = paste0("variant", seq_len(n_markers)),
    seqname = rep(c("chr1", "chr2", "chr3"), each = n_markers / 3L),
    pos = rep(seq.int(100L, by = 100L, length.out = n_markers / 3L), 3L),
    ref = "A",
    alt = "G",
    genotype,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  causal <- c(5L, 31L, 70L)
  true_value <- colSums(genotype[causal, , drop = FALSE])
  phenotype <- data.frame(
    sample = samples,
    pop = rep(c("pop1", "pop2"), each = n_samples / 2L),
    family = rep(paste0("family", 1:8), each = 6L),
    trait = true_value + stats::rnorm(n_samples, sd = 0.8),
    true_breeding_value = true_value,
    stringsAsFactors = FALSE
  )
  directory <- tempfile("simitall-genomic-")
  dir.create(directory)
  genotype_file <- file.path(directory, "cohort.geno.tsv")
  utils::write.table(
    variants,
    genotype_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  list(
    directory = directory,
    genotype_file = genotype_file,
    phenotype = phenotype,
    truth = data.frame(
      marker_id = variants$id[causal],
      seqname = variants$seqname[causal],
      pos = variants$pos[causal],
      effect = 1,
      stringsAsFactors = FALSE
    )
  )
}

test_that("rrBLUP GWAS writes associations and truth benchmarks", {
  skip_if_not_installed("rrBLUP")
  fixture <- make_genomic_fixture()
  prefix <- file.path(fixture$directory, "association")
  result <- analyze_gwas(
    fixture$genotype_file,
    fixture$phenotype,
    prefix,
    fixed_effects = "pop",
    n_pcs = 1L
  )
  expect_equal(nrow(result$results), 90L)
  expect_true(all(c(
    "marker_id", "p_value", "q_value", "maf", "significant"
  ) %in% names(result$results)))
  expect_true(all(file.exists(unlist(result$paths))))

  benchmark <- benchmark_gwas(
    result$results,
    fixture$truth,
    file.path(fixture$directory, "benchmark")
  )
  expect_equal(benchmark$metrics$causal_variants, 3L)
  expect_true(all(file.exists(unlist(benchmark$paths))))

  if (requireNamespace("ggplot2", quietly = TRUE) &&
      requireNamespace("patchwork", quietly = TRUE)) {
    figures <- plot_gwas_results(
      benchmark$results,
      file.path(fixture$directory, "association_plot"),
      width = 7,
      height = 4,
      dpi = 72
    )
    expect_true(all(file.exists(unlist(figures$paths))))
  }
})

test_that("genomic models predict, select parents, and design crosses", {
  skip_if_not_installed("rrBLUP")
  fixture <- make_genomic_fixture()
  training <- fixture$phenotype$sample[seq_len(36L)]
  fit <- fit_genomic_model(
    fixture$genotype_file,
    fixture$phenotype,
    model = "gblup",
    training_samples = training,
    seed = 902
  )
  predictions <- predict_genomic_values(fit)
  expect_s3_class(fit, "simitall_genomic_model")
  expect_equal(nrow(predictions), 48L)
  expect_true(all(is.finite(predictions$gebv)))

  kinship <- calculate_genomic_kinship(fixture$genotype_file)
  selected <- select_parents(
    predictions,
    n_parents = 8L,
    kinship = kinship,
    diversity_penalty = 0.2,
    family_data = fixture$phenotype,
    max_per_family = 2L
  )
  expect_equal(nrow(selected), 8L)
  expect_lte(max(table(selected$family)), 2L)

  crosses <- design_parent_crosses(
    selected,
    n_crosses = 6L,
    mating = "minimum_kinship",
    kinship = kinship
  )
  expect_equal(nrow(crosses), 6L)
  expect_true(all(c("parent1", "parent2", "kinship") %in% names(crosses)))
})

test_that("marker-effect, Bayesian, and forest backends return predictions", {
  fixture <- make_genomic_fixture()
  training <- fixture$phenotype$sample[seq_len(36L)]

  skip_if_not_installed("rrBLUP")
  rr_fit <- fit_genomic_model(
    fixture$genotype_file,
    fixture$phenotype,
    model = "rrblup",
    training_samples = training
  )
  expect_true(all(is.finite(predict_genomic_values(rr_fit)$gebv)))

  if (requireNamespace("BGLR", quietly = TRUE)) {
    bayes_fit <- fit_genomic_model(
      fixture$genotype_file,
      fixture$phenotype,
      model = "bayes_ridge",
      training_samples = training,
      n_iter = 300L,
      burn_in = 100L,
      thin = 5L
    )
    expect_true(all(is.finite(predict_genomic_values(bayes_fit)$gebv)))
  }

  if (requireNamespace("ranger", quietly = TRUE)) {
    forest_fit <- fit_genomic_model(
      fixture$genotype_file,
      fixture$phenotype,
      model = "random_forest",
      training_samples = training,
      num_trees = 50L
    )
    expect_true(all(is.finite(predict_genomic_values(forest_fit)$gebv)))
  }
})

test_that("cross-validation and one-round genomic selection write outputs", {
  skip_if_not_installed("rrBLUP")
  fixture <- make_genomic_fixture()
  cv <- cross_validate_genomic_prediction(
    fixture$genotype_file,
    fixture$phenotype,
    model = "gblup",
    folds = 3L,
    true_value = "true_breeding_value",
    seed = 903
  )
  expect_equal(nrow(cv$metrics), 3L)
  expect_true(is.finite(cv$summary$truth_accuracy))

  selection <- run_genomic_selection(
    fixture$genotype_file,
    fixture$phenotype,
    file.path(fixture$directory, "selection"),
    model = "gblup",
    n_parents = 6L,
    n_crosses = 5L,
    diversity_penalty = 0.1,
    seed = 904
  )
  expect_equal(nrow(selection$selected), 6L)
  expect_equal(nrow(selection$crosses), 5L)
  expect_true(all(file.exists(unlist(selection$paths))))
})
