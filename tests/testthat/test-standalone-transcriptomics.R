test_that("standalone bulk RNA-seq distinguishes technical replicates", {
  out_dir <- tempfile("simitall-standalone-rna-")
  dir.create(out_dir)
  prefix <- file.path(out_dir, "sample_a")

  paths <- simulate_rnaseq_experiment(
    out_prefix = prefix,
    sample_names = "sample_A",
    replicates = 3,
    replicate_mode = "technical",
    n_genes = 20,
    n_latent = 1,
    library_size_mean = 1e4,
    seed = 101
  )

  expect_true(all(file.exists(unlist(paths))))
  metadata <- read.delim(paths$sample_metadata, check.names = FALSE)
  counts <- read.delim(paths$counts, check.names = FALSE)
  expect_equal(nrow(metadata), 3)
  expect_equal(length(unique(metadata$biological_id)), 1)
  expect_equal(length(unique(metadata$library_id)), 3)
  expect_equal(unique(metadata$replicate_mode), "technical")
  expect_equal(dim(counts), c(20, 4))
})

test_that("custom bulk designs expand biological and technical replication", {
  out_dir <- tempfile("simitall-standalone-design-")
  dir.create(out_dir)
  design <- data.frame(
    sample = c("control", "treated"),
    condition = c("control", "treated"),
    batch = c("batch1", "batch2"),
    biological_replicates = c(2, 2),
    technical_replicates = c(2, 2),
    stringsAsFactors = FALSE
  )

  paths <- simulate_rnaseq_experiment(
    out_prefix = file.path(out_dir, "mixed"),
    design = design,
    n_genes = 20,
    condition_effect_fraction = 0.5,
    n_latent = 0,
    library_size_mean = 1e4,
    seed = 102
  )

  metadata <- read.delim(paths$sample_metadata, check.names = FALSE)
  condition_truth <- read.delim(paths$condition_truth, check.names = FALSE)
  expect_equal(nrow(metadata), 8)
  expect_equal(length(unique(metadata$biological_id)), 4)
  expect_equal(length(unique(metadata$library_id)), 8)
  expect_setequal(unique(metadata$condition), c("control", "treated"))
  expect_equal(nrow(condition_truth), 10)
})

test_that("standalone single-cell technical captures share one donor", {
  out_dir <- tempfile("simitall-standalone-scrna-")
  dir.create(out_dir)
  prefix <- file.path(out_dir, "sample_a")

  paths <- simulate_scrnaseq_experiment(
    out_prefix = prefix,
    sample_names = "sample_A",
    replicates = 3,
    replicate_mode = "technical",
    n_genes = 24,
    cells_per_replicate = 15,
    cell_type_proportions = c(T_cell = 0.4, B_cell = 0.3, Monocyte = 0.3),
    markers_per_cell_type = 3,
    n_latent = 1,
    library_size_mean = 1000,
    backend = "native",
    write_sce = TRUE,
    min_cells_pseudobulk = 3,
    seed = 103
  )

  expect_true(all(file.exists(unlist(paths))))
  object <- readRDS(paths$object_rds)
  cells <- object$cell_metadata
  expect_s3_class(object, "simitall_scrnaseq")
  expect_identical(object$experiment_type, "standalone")
  expect_equal(dim(object$counts), c(24, 45))
  expect_equal(length(unique(cells$donor)), 1)
  expect_equal(length(unique(cells$library_id)), 3)
  expect_equal(nrow(read.delim(paths$pseudobulk_metadata)), 3)
  expect_equal(nrow(read.delim(paths$library_pseudobulk_metadata)), 9)
  if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    expect_s4_class(readRDS(paths$sce_rds), "SingleCellExperiment")
  }
})
