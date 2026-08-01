test_that("GWAS donors can drive single-cell and pseudobulk simulation", {
  out_dir <- tempfile("simitall-scrna-")
  dir.create(out_dir)
  genome <- system.file(
    "extdata", "examples", "demo_genome.fa",
    package = "simitall"
  )
  gwas_prefix <- file.path(out_dir, "cohort")
  scrna_prefix <- file.path(out_dir, "scrna")

  simulate_gwas_cohort(
    genome_fa = genome,
    out_prefix = gwas_prefix,
    n_samples = 12,
    snp_rate = 0.02,
    indel_rate = 0,
    n_causal = 2,
    seed = 71
  )
  paths <- simulate_scrnaseq_from_gwas(
    genotype_file = paste0(gwas_prefix, ".geno.tsv"),
    out_prefix = scrna_prefix,
    sample_metadata = paste0(gwas_prefix, ".pheno.tsv"),
    n_genes = 24,
    cells_per_donor = 15,
    cell_type_proportions = c(T_cell = 0.4, B_cell = 0.3, Monocyte = 0.3),
    n_cis_eqtl = 5,
    n_trans_eqtl = 2,
    cis_window_bp = 700,
    markers_per_cell_type = 3,
    n_latent = 1,
    library_size_mean = 1000,
    backend = "native",
    write_dense_tsv = TRUE,
    write_sce = TRUE,
    min_cells_pseudobulk = 3,
    seed = 72
  )

  expect_true(all(file.exists(unlist(paths))))
  expect_identical(attr(paths, "backend"), "native")
  object <- readRDS(paths$object_rds)
  expect_s3_class(object, "simitall_scrnaseq")
  expect_equal(dim(object$counts), c(24, 180))
  expect_equal(length(unique(object$cell_metadata$donor)), 12)
  expect_equal(length(unique(object$cell_metadata$cell_type)), 3)
  expect_true(all(c("shared", "cell_type_specific") %in%
                    unique(object$eqtl_truth$scope)))

  sparse <- Matrix::readMM(paths$matrix)
  expect_equal(dim(sparse), c(24, 180))
  pseudobulk <- read.delim(paths$pseudobulk_metadata, check.names = FALSE)
  expect_equal(nrow(pseudobulk), 36)
  expect_true(all(pseudobulk$n_cells >= 3))
  expect_true(all(is.finite(pseudobulk$library_size)))
  if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    expect_true(file.exists(paths$sce_rds))
    expect_s4_class(readRDS(paths$sce_rds), "SingleCellExperiment")
  }
})

test_that("cell-type eQTL benchmarking and figures use donor pseudobulk", {
  out_dir <- tempfile("simitall-scrna-benchmark-")
  dir.create(out_dir)
  genome <- system.file(
    "extdata", "examples", "demo_genome.fa",
    package = "simitall"
  )
  gwas_prefix <- file.path(out_dir, "cohort")
  scrna_prefix <- file.path(out_dir, "scrna")
  eqtl_prefix <- file.path(out_dir, "celltype_eqtl")

  simulate_gwas_cohort(
    genome_fa = genome,
    out_prefix = gwas_prefix,
    n_samples = 12,
    snp_rate = 0.02,
    indel_rate = 0,
    n_causal = 2,
    seed = 81
  )
  paths <- simulate_scrnaseq_from_gwas(
    genotype_file = paste0(gwas_prefix, ".geno.tsv"),
    out_prefix = scrna_prefix,
    sample_metadata = paste0(gwas_prefix, ".pheno.tsv"),
    n_genes = 24,
    cells_per_donor = 15,
    cell_type_proportions = c(T_cell = 0.4, B_cell = 0.3, Monocyte = 0.3),
    n_cis_eqtl = 5,
    n_trans_eqtl = 2,
    cis_window_bp = 700,
    markers_per_cell_type = 3,
    n_latent = 1,
    library_size_mean = 1000,
    backend = "native",
    write_sce = FALSE,
    min_cells_pseudobulk = 3,
    seed = 82
  )
  benchmark <- benchmark_celltype_eqtls(
    genotype_file = paste0(gwas_prefix, ".geno.tsv"),
    pseudobulk_expression = paths$pseudobulk_expression,
    pseudobulk_metadata = paths$pseudobulk_metadata,
    truth_file = paths$eqtl_truth,
    out_prefix = eqtl_prefix,
    cis_window_bp = 700,
    min_donors = 8
  )

  expect_true(file.exists(benchmark$results_path))
  expect_true(file.exists(benchmark$metrics_path))
  expect_setequal(
    unique(benchmark$metrics$cell_type),
    c("T_cell", "B_cell", "Monocyte", "Overall")
  )
  expect_true(all(c("cell_type", "is_causal", "q_value") %in%
                    names(benchmark$results)))

  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  figure <- plot_scrnaseq_results(
    scrna_prefix = scrna_prefix,
    eqtl_prefix = eqtl_prefix,
    out_prefix = file.path(out_dir, "scrna_figure"),
    top_variable_genes = 20,
    markers_per_cell_type = 2,
    max_plot_cells = 200,
    width = 10,
    height = 7,
    dpi = 72
  )
  expect_true(file.exists(figure$paths$pdf))
  expect_true(file.exists(figure$paths$png))
  expect_true(file.exists(figure$paths$summary))
  expect_length(
    list.files(figure$paths$source_data, pattern = "\\.csv$"),
    6
  )
})

test_that("Splatter is available as an optional technical backend", {
  skip_if_not_installed("splatter")
  skip_if_not_installed("SummarizedExperiment")
  scaffold <- .simitall_scrna_scaffold(
    backend = "splatter",
    n_genes = 10,
    n_cells = 12,
    baseline = rep(5, 10),
    library_sizes = rep(1000, 12),
    seed = 91
  )
  expect_identical(scaffold$backend, "splatter")
  expect_equal(dim(scaffold$expected), c(10, 12))
})

test_that("generic eQTL truth can be assigned to cell types", {
  truth <- data.frame(
    gene_id = paste0("gene", seq_len(10)),
    variant_id = paste0("variant", seq_len(10)),
    effect_log2 = seq(-0.5, 0.5, length.out = 10),
    stringsAsFactors = FALSE
  )
  result <- simulate_celltype_eqtls(
    truth,
    cell_types = c("T_cell", "B_cell"),
    cell_type_specific_fraction = 0.4,
    seed = 1
  )
  expect_equal(sum(result$scope == "cell_type_specific"), 4)
  expect_equal(sum(result$scope == "shared"), 6)
  expect_true(all(result$cell_type[result$scope == "shared"] == "all"))
  expect_true(all(result$cell_type[result$scope == "cell_type_specific"] %in%
                    c("T_cell", "B_cell")))
})
