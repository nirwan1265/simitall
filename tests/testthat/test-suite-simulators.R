# Routes from Biodes spec fields to simulators (inst/suite/simulators.tsv) and
# their parameters (inst/suite/simulator_parameters.tsv).

test_that("routes come from spec fields, not words, and unsupported requests are listed", {
  spec <- list(primary_outcome = list(type = "count"), assays = list("RNA-seq", "gwas", "lipidomics", "methylation"),
               population = list(mating = "RIL"))
  r <- simitall:::.suite_match_routes(spec)
  expect_equal(r$route[r$layer == "outcome"], "outcome_count")
  expect_equal(r$route[r$layer == "genetics"], "genetics_cross")          # the cross wins over the panel
  expect_true("assay_rnaseq" %in% r$route)
  expect_equal(r$status[r$route == "assay_methylation"], "planned")
  expect_equal(r$status[r$route == "assay_other"], "not_supported")
  # no outcome type: continuous; eQTL alone turns on the panel
  r2 <- simitall:::.suite_match_routes(list(assays = list("eqtl")))
  expect_equal(r2$route[r2$layer == "outcome"], "outcome_continuous")
  expect_true("genetics_panel" %in% r2$route)
  # every simulator named in the table exists
  tab <- simitall:::.suite_tsv("simulators")
  fns <- setdiff(tab$simulator, c("-", "suite_outcome", "suite_genotype_panel", "suite_genotype_cross"))
  expect_true(all(vapply(fns, exists, logical(1), envir = asNamespace("simitall"))))
})

test_that("parameters: params beat spec fields, spec fields beat defaults, and every argument is listed", {
  spec <- list(population = list(generations = "F6", parents = list("B73", "Mo17", "CML247", "Oh43")))
  a <- simitall:::.suite_args("suite_genotype_cross", spec, list(suite_genotype_cross = list(snp_rate = 0.05)))
  expect_equal(a$args$snp_rate, 0.05)
  expect_equal(a$source[["snp_rate"]], "params")
  expect_equal(a$source[["self_generations"]], "spec population.generations")
  expect_equal(simitall:::.suite_count(a$args$self_generations, filial = TRUE), 5L)   # F6 = F1 + 5 selfings
  expect_equal(simitall:::.suite_count(a$args$n_founders), 4L)                         # four parent names
  tab <- simitall_simulator_parameters("simulate_chipseq")
  expect_true(all(c("n_peaks", "fragment_mean", "duplicate_rate") %in% tab$argument))
  expect_equal(tab$bridge_value[tab$argument == "backend"], "native")
  tab_types <- simitall:::.suite_tsv("simulator_parameters")
  expect_true(all(tab_types$type %in% c("integer", "numeric", "logical", "character", "numeric_vector")))
})

test_that("binary, count and ordinal outcomes follow simulator_parameters.tsv", {
  for (ty in c("binary", "count", "ordinal")) {
    d <- make_design_folder("rcbd", "y ~ treatment + block", list(design_factor("genotype", 10)), reps = 4)
    s <- jsonlite::read_json(file.path(d, "study_spec.json"))
    s$primary_outcome$type <- ty
    jsonlite::write_json(s, file.path(d, "study_spec.json"), auto_unbox = TRUE)
    suppressMessages(simitall_from_biodes(d, params = list(suite_outcome = list(prevalence = 0.25))))
    y <- read.delim(file.path(d, "sample_metadata_template.tsv"))$primary_outcome
    t <- jsonlite::read_json(file.path(d, "sim_truth.json"))
    expect_equal(t$outcome_type, ty)
    if (ty == "binary") expect_equal(mean(y), 0.25)
    if (ty == "count") expect_true(all(y >= 0 & y == round(y)))
    if (ty == "ordinal") expect_setequal(unique(y), 1:4)
  }
  d <- make_design_folder("rcbd", "y ~ treatment + block", list(design_factor("genotype", 4)))
  s <- jsonlite::read_json(file.path(d, "study_spec.json")); s$primary_outcome$type <- "survival"
  jsonlite::write_json(s, file.path(d, "study_spec.json"), auto_unbox = TRUE)
  expect_error(simitall_from_biodes(d), "no simulator yet")
})

test_that("GWAS designs get line genotypes whose planted QTN a GWAS recovers", {
  skip_if_not_installed("rrBLUP")
  d <- make_design_folder("rcbd", "y ~ treatment + block", list(design_factor("accession", 150)), reps = 2,
                          extra = list(assays = list("gwas")))
  suppressMessages(simitall_from_biodes(d, params = list(suite_genetics = list(heritability = 0.8, n_add_qtn = 4))))
  t <- jsonlite::read_json(file.path(d, "sim_truth.json"))
  expect_equal(t$genetics$genotype_term, "treatment")
  expect_equal(t$genetics$n_lines, 150)
  gv <- read.delim(file.path(d, "genetic_truth", "genetic_values.tsv"))
  expect_equal(unlist(t$treatment_effects$treatment)[gv$line], stats::setNames(gv$genetic_value, gv$line), tolerance = 1e-8)
  expect_equal(stats::var(gv$genetic_value), 0.8 / 0.2 * 10^2, tolerance = 1e-6)   # h2 / (1 - h2) * residual SD^2
  m <- read.delim(file.path(d, "sample_metadata_template.tsv"))
  ph <- stats::aggregate(primary_outcome ~ treatment, m, mean); names(ph) <- c("sample", "yield")
  g <- suppressMessages(utils::capture.output(res <- analyze_gwas(file.path(d, "genotypes.vcf"), ph, file.path(d, "gwas"), trait = "yield", n_pcs = 3L, min_maf = 0.05)))
  out <- read.delim(res$paths$results)
  q <- read.delim(file.path(d, "genetic_truth", "qtn_truth.tsv"))
  top <- out$marker_id[order(out$p_value)][1:10]
  expect_gte(sum(q$marker_id %in% top), 3)
})

test_that("GxE phenotype model splits genetic variance and fills the genotype x treatment cells", {
  d <- make_design_folder("factorial", "y ~ treatment * treatment2 + block",
                          list(design_factor("genotype", 30), design_factor("water", 2, c("wet", "dry"))), reps = 2,
                          extra = list(assays = list("gwas")))
  suppressMessages(simitall_from_biodes(d, params = list(suite_genetics = list(phenotype_model = "gxe", heritability = 0.6))))
  t <- jsonlite::read_json(file.path(d, "sim_truth.json"))
  vg <- 0.6 / 0.4 * 100
  expect_equal(t$genetics$realized_main_variance, 0.7 * vg, tolerance = 1e-6)
  expect_equal(t$genetics$realized_gxe_variance, 0.3 * vg, tolerance = 1e-6)
  expect_equal(t$interaction_effects[["treatment:treatment2"]]$pattern, "phenotype_model_gxe")
  expect_equal(t$genetics$environment_term, "treatment2")
})

test_that("a RIL cross gives inbred lines with founder ancestry", {
  d <- make_design_folder("rcbd", "y ~ treatment + block", list(design_factor("ril", 20)), reps = 2,
                          extra = list(population = list(mating = "ril", generations = "F6")))
  suppressMessages(simitall_from_biodes(d, params = list(suite_genotype_cross = list(chromosome_length_bp = 4000L, n_chromosomes = 2L, snp_rate = 0.05))))
  t <- jsonlite::read_json(file.path(d, "sim_truth.json"))
  expect_equal(t$genetics$genotype_source$scheme, "RIL")
  expect_equal(t$genetics$genotype_source$self_generations, 5)
  expect_true(file.exists(file.path(d, "genetic_truth", "founder_ancestry.tsv")))
  v <- simitall:::.simitall_read_genotypes(file.path(d, "genotypes.vcf"))
  expect_lt(mean(v$genotype == 1), 0.1)                     # 5 selfings: ~3% heterozygous
  expect_setequal(colnames(v$genotype), paste0("ril_", 1:20))
})

test_that("omics assays write one library per unit, named by sample_id", {
  d <- make_design_folder("factorial", "y ~ treatment * treatment2 + block",
                          list(design_factor("genotype", 4), design_factor("water", 2, c("wet", "dry"))), reps = 2,
                          extra = list(assays = list("eqtl", "single_cell", "chip_seq", "atac_seq", "methylation")))
  suppressMessages(simitall_from_biodes(d, params = list(
    simulate_rnaseq_from_gwas = list(n_genes = 200), simulate_scrnaseq_experiment = list(n_genes = 200, cells_per_replicate = 40),
    simulate_chipseq = list(n_peaks = 40, n_reads = 8000), simulate_atacseq = list(n_peaks = 40, n_reads = 8000),
    suite_peak_genome = list(chromosome_length = 100000L, n_genes = 80L))))
  meta <- read.delim(file.path(d, "sample_metadata_template.tsv"))
  cnt <- read.delim(file.path(d, "counts.tsv"), check.names = FALSE)
  expect_identical(names(cnt)[-1], meta$sample_id)                         # eQTL owns counts.tsv
  sc <- read.delim(file.path(d, "single_cell", "unit_counts.tsv"), check.names = FALSE)
  expect_identical(names(sc)[-1], meta$sample_id)
  chip <- read.delim(file.path(d, "chipseq", "chip.sample_metadata.tsv"))
  water <- stats::setNames(as.character(meta$treatment2), meta$sample_id)
  expect_setequal(chip$sample, meta$sample_id)
  expect_true(all(chip$condition == water[chip$sample]))
  atac <- read.delim(file.path(d, "atacseq", "library_to_unit.tsv"))
  expect_setequal(stats::na.omit(atac$sample_id), meta$sample_id)
  t <- jsonlite::read_json(file.path(d, "sim_truth.json"))
  expect_setequal(names(t$assays), c("assay_eqtl", "assay_single_cell", "assay_chip_seq", "assay_atac_seq"))
  expect_equal(t$not_simulated[[1]]$route, "assay_methylation")
  # each unit's eQTL genotypes are its line's genotypes
  lines <- simitall:::.simitall_read_genotypes(file.path(d, "genotypes.vcf"))$genotype
  units <- simitall:::.simitall_read_genotypes(file.path(d, "eqtl", "unit_genotypes.vcf"))$genotype
  u <- meta$sample_id[1]
  expect_equal(unname(units[, u]), unname(lines[, as.character(meta$treatment[1])]))
})
