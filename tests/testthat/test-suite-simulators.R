# Routes from Biodes spec fields to simulators (inst/suite/simulators.tsv) and
# their parameters (inst/suite/simulator_parameters.tsv).

test_that("routes come from spec fields, not words, and unsupported requests are listed", {
  spec <- list(primary_outcome = list(type = "count"), assays = list("RNA-seq", "gwas", "lipidomics", "methylation"),
               population = list(mating = "RIL"))
  r <- simitall:::.suite_match_routes(spec)
  expect_equal(r$route[r$layer == "outcome"], "outcome_count")
  # amplicon is 16S/ITS microbiome, not genotyping; CRISPR / overexpression are treatments, not data
  r3 <- simitall:::.suite_match_routes(list(assays = list("amplicon", "crispr")))
  expect_false(any(r3$layer == "genetics"))
  expect_equal(r3$status[r3$route == "assay_microbiome"], "ready")
  expect_equal(r3$status[r3$route == "treatment_only"], "treatment")
  expect_equal(r$route[r$layer == "genetics"], "genetics_cross")          # the cross wins over the panel
  expect_true("assay_rnaseq" %in% r$route)
  expect_equal(r$status[r$route == "assay_methylation"], "ready")
  expect_equal(simitall:::.suite_match_routes(list(assays = list("medip")))$status[2], "planned")
  expect_equal(r$status[r$route == "assay_other"], "not_supported")
  # no outcome type: continuous; eQTL alone turns on the panel
  r2 <- simitall:::.suite_match_routes(list(assays = list("eqtl")))
  expect_equal(r2$route[r2$layer == "outcome"], "outcome_continuous")
  expect_true("genetics_panel" %in% r2$route)
  # every simulator named in the table exists
  tab <- simitall:::.suite_tsv("simulators")
  fns <- setdiff(tab$simulator[tab$layer == "assay"], "-")
  runs <- vapply(fns, function(f) is.function(tryCatch(simitall:::.suite_assay_runner(f), error = function(e) NULL)), logical(1))
  expect_true(all(runs))
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
                          extra = list(assays = list("eqtl", "single_cell", "chip_seq", "atac_seq", "medip")))
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
  water <- stats::setNames(as.character(meta$water), meta$sample_id)
  expect_setequal(chip$sample, meta$sample_id)
  expect_true(all(chip$condition == water[chip$sample]))
  atac <- read.delim(file.path(d, "atacseq", "library_to_unit.tsv"))
  expect_setequal(stats::na.omit(atac$sample_id), meta$sample_id)
  t <- jsonlite::read_json(file.path(d, "sim_truth.json"))
  expect_setequal(names(t$assays), c("assay_eqtl", "assay_single_cell", "assay_chip_seq", "assay_atac_seq"))
  expect_equal(t$not_simulated[[1]]$route, "assay_medip")
  # each unit's eQTL genotypes are its line's genotypes
  lines <- simitall:::.simitall_read_genotypes(file.path(d, "genotypes.vcf"))$genotype
  units <- simitall:::.simitall_read_genotypes(file.path(d, "eqtl", "unit_genotypes.vcf"))$genotype
  u <- meta$sample_id[1]
  expect_equal(unname(units[, u]), unname(lines[, as.character(meta$genotype[1])]))
})

test_that("factors the model does not hold (tissue within each plant) get effects, and a plant's rows share its effect", {
  d <- make_design_folder("rcbd", NULL, list(design_factor("genotype", 3), design_factor("water", 2), design_factor("tissue", 2)), reps = 4,
                          extra = list(assays = list("rnaseq", "crispr")))
  s <- jsonlite::read_json(file.path(d, "study_spec.json"))
  expect_equal(unlist(s$analysis$terms$uncovered_factors), "tissue")
  suppressMessages(simitall_from_biodes(d, params = list(simulate_rnaseq_experiment = list(n_genes = 200))))
  t <- jsonlite::read_json(file.path(d, "sim_truth.json"))
  expect_equal(names(t$uncovered_effects), "tissue")
  expect_gt(t$within_unit_sd$sd, 0)
  expect_null(t$not_simulated)                                  # CRISPR is a treatment, not a missing simulator
  expect_equal(unlist(t$treatment_assays), "crispr")
  meta <- read.delim(file.path(d, "sample_metadata_template.tsv"))
  expect_equal(nrow(meta), 3 * 2 * 4 * 2)                        # one row per plant and tissue
  # RNA-seq condition = water x tissue: true DE for every non-reference combination
  expect_equal(length(t$counts$de_genes), 2 * 2 - 1)
  # rows of one plant share its effect: leaf and root of a plant correlate after removing fixed effects
  fit <- stats::lm(primary_outcome ~ genotype * water + tissue + block, data = meta)
  r <- stats::residuals(fit)
  pr <- split(r, meta$plot)
  pairs <- do.call(rbind, pr[lengths(pr) == 2])
  expect_gt(stats::cor(pairs[, 1], pairs[, 2]), 0.2)
})

test_that("a factor named 'treatment' that is not the first factor maps through analysis.terms", {
  d <- make_design_folder("factorial", NULL, list(design_factor("soil", 2), design_factor("treatment", 3)), reps = 3)
  s <- jsonlite::read_json(file.path(d, "study_spec.json"))
  expect_equal(s$analysis$terms$treatment, "soil")
  expect_equal(s$analysis$terms$treatment2, "treatment")
  h <- simitall:::.suite_read_biodes(d)
  expect_equal(unname(h$parsed$col[c("treatment", "treatment2")]), c("soil", "treatment"))
  suppressMessages(simitall_from_biodes(d))
  t <- jsonlite::read_json(file.path(d, "sim_truth.json"))
  expect_length(t$treatment_effects$treatment, 2)
  expect_length(t$treatment_effects$treatment2, 3)
})

test_that("templates: each is real integer counts with a sample sheet and provenance", {
  tab <- simitall_templates()
  expect_true(all(c("small_rnaseq_human_mirna_atlas", "small_rnaseq_mouse_mirna_atlas", "microbiome_plant_root_16s", "microbiome_human_gut_hmp_16s") %in% tab$id))
  expect_true(all(nzchar(tab$source) & nzchar(tab$license)))
  for (id in tab$id) {
    t <- simitall_template(id)
    m <- t$counts %||% t$beta
    expect_identical(colnames(m), t$samples$sample)
    expect_true(t$info$default_group %in% t$samples$group)
    if (t$info$kind %in% c("counts", "methylation_counts")) expect_true(is.integer(t$counts) && min(t$counts) >= 0)
    if (t$info$kind == "methylation_counts") expect_true(all(t$methylated <= t$counts))
    if (t$info$kind == "beta") expect_true(all(t$beta >= 0 & t$beta <= 1))
  }
  expect_true(all(c("small_rnaseq_zea_mays_pmiren", "small_rnaseq_sorghum_bicolor_pmiren", "methylation_mouse_encode_wgbs") %in% tab$id))
})

test_that("template engine recovers the planted fold changes in both modes", {
  d <- data.frame(sample = paste0("u", 1:10), condition = rep(c("ctrl", "heat"), each = 5), batch = "b1")
  for (m in c("resample", "parametric")) {
    r <- simulate_counts_from_template("small_rnaseq_mouse_mirna_atlas", d, file.path(tempfile(), m), mode = m, seed = 2)
    y <- as.matrix(read.delim(r$counts, row.names = 1, check.names = FALSE)); tr <- read.delim(r$condition_truth)
    expect_identical(colnames(y), d$sample)
    cpm <- log2(t(t(y) / colSums(y)) * 1e6 + 1)
    obs <- rowMeans(cpm[, 6:10]) - rowMeans(cpm[, 1:5])
    expect_gt(stats::cor(obs[tr$gene_id], tr$effect_log2), 0.6)
    expect_lt(stats::median(abs(obs[setdiff(rownames(y), tr$gene_id)])), 0.2)
  }
  if (TRUE) expect_error(simulate_counts_from_template("small_rnaseq_human_mirna_atlas", data.frame(sample = paste0("u", 1:30), condition = "a"), tempfile(), mode = "resample"), "one template sample per unit")
})

test_that("small RNA-seq and microbiome designs pick the template by organism", {
  run <- function(organism, system, assay) {
    d <- make_design_folder("crd", NULL, list(design_factor("treatment", 2)), reps = 4,
                            extra = list(organism = organism, system = system, assays = list(assay)))
    suppressMessages(simitall_from_biodes(d))
    list(t = jsonlite::read_json(file.path(d, "sim_truth.json")), meta = read.delim(file.path(d, "sample_metadata_template.tsv")),
         counts = read.delim(file.path(d, "counts.tsv"), check.names = FALSE))
  }
  h <- run("Homo sapiens", "clinical", "small_rnaseq")
  expect_equal(h$t$assays$assay_small_rnaseq$template_choice$template, "small_rnaseq_human_mirna_atlas")
  expect_identical(names(h$counts)[-1], h$meta$sample_id)
  expect_true(any(grepl("^hsa-", h$counts$gene_id)))
  p <- run("Sorghum bicolor", "greenhouse", "small_rnaseq")
  expect_equal(p$t$assays$assay_small_rnaseq$template_choice$template, "small_rnaseq_sorghum_bicolor_pmiren")
  q <- run("Oryza sativa", "field", "small_rnaseq")                    # no rice template: nearest plant one
  expect_match(q$t$assays$assay_small_rnaseq$template_choice$why, "same kingdom")
  y <- run("Saccharomyces cerevisiae", "cell_culture", "small_rnaseq") # no kingdom match: native negative binomial
  expect_equal(y$t$assays$assay_small_rnaseq$template_choice$mode, "native_negative_binomial")
  m <- run("Zea mays", "field", "amplicon")
  expect_equal(m$t$assays$assay_microbiome$template_choice$template, "microbiome_plant_root_16s")
  expect_match(m$t$assays$assay_microbiome$template_choice$why, "same kingdom")
  g <- run("Mus musculus", "mouse", "16s")
  expect_equal(g$t$assays$assay_microbiome$template_choice$template, "microbiome_human_gut_hmp_16s")
  expect_gt(length(g$t$counts$de_genes[[1]]), 0)
})

test_that("methylation engines plant DMRs that the data show", {
  d <- data.frame(sample = paste0("u", 1:8), condition = rep(c("ctrl", "trt"), each = 4), batch = "b1")
  check <- function(r) {
    M <- as.matrix(read.delim(r$methylated, row.names = 1, check.names = FALSE)); N <- as.matrix(read.delim(r$coverage, row.names = 1, check.names = FALSE))
    expect_true(all(M <= N)); expect_identical(colnames(M), d$sample)
    st <- read.delim(r$site_truth); dt <- read.delim(r$dmr_truth)
    lvl <- function(i, s) colSums(M[i, s, drop = FALSE]) / pmax(colSums(N[i, s, drop = FALSE]), 1)
    obs <- vapply(dt$dmr_id, function(id) { i <- which(st$dmr_id == id); mean(lvl(i, 5:8)) - mean(lvl(i, 1:4)) }, 0)
    expect_lt(stats::median(abs(obs - dt$delta)), 0.06)                 # measured change matches the planted one
    expect_gt(mean(sign(obs) == sign(dt$delta)), 0.9)
    expect_true(all(sign(dt$delta) == ifelse(dt$base < 0.5, 1, -1)))     # gain in unmethylated regions, loss otherwise
    dt
  }
  p <- simulate_methylation(d, file.path(tempfile(), "p"), engine = "native", kingdom = "plant", n_sites = 6000, n_dmrs = 20)
  expect_setequal(p$truth$contexts, c("CG", "CHG", "CHH"))
  check(p)
  a <- simulate_methylation(d, file.path(tempfile(), "a"), engine = "native", kingdom = "animal", n_sites = 6000, n_dmrs = 20)
  expect_equal(a$truth$contexts, "CG")
  t <- simulate_methylation(d, file.path(tempfile(), "t"), engine = "template", template = "methylation_mouse_encode_wgbs", n_dmrs = 20)
  expect_equal(t$truth$mode, "parametric")                               # 8 units, 3 template samples
  check(t)
  r <- simulate_methylation(d[1:3, ], file.path(tempfile(), "r"), engine = "template", template = "methylation_mouse_encode_wgbs", n_dmrs = 20)
  expect_equal(r$truth$mode, "resample")
  ar <- simulate_methylation(d, file.path(tempfile(), "ar"), template = "methylation_array_human_tcga_tgct", n_dmrs = 50)
  B <- as.matrix(read.delim(ar$beta, row.names = 1, check.names = FALSE)); tr <- read.delim(ar$dmr_truth)
  expect_true(all(B >= 0 & B <= 1))
  expect_gt(stats::cor(rowMeans(B[tr$feature_id, 5:8]) - rowMeans(B[tr$feature_id, 1:4]), tr$delta), 0.8)
})

test_that("methylation designs: plants get the native CG/CHG/CHH engine, animals the WGBS template, and the questions are listed", {
  run <- function(organism, system, assays = list("methylation")) {
    d <- make_design_folder("crd", NULL, list(design_factor("treatment", 2)), reps = 3,
                            extra = list(organism = organism, system = system, assays = assays))
    r <- utils::capture.output(rt <- simitall_routes(d))
    suppressMessages(simitall_from_biodes(d, params = list(simulate_methylation = list(n_sites = 4000, n_dmrs = 10))))
    list(t = jsonlite::read_json(file.path(d, "sim_truth.json"))$assays$assay_methylation, q = rt$questions, d = d)
  }
  m <- run("Zea mays", "field")
  expect_equal(m$t$engine, "native"); expect_equal(m$t$kingdom, "plant")
  expect_true("dmr_delta" %in% m$q$argument)
  expect_false("prevalence" %in% m$q$argument)                          # continuous outcome: no prevalence question
  h <- run("Homo sapiens", "clinical")
  expect_equal(h$t$template, "methylation_mouse_encode_wgbs")
  meta <- read.delim(file.path(h$d, "sample_metadata_template.tsv"))
  cov <- read.delim(file.path(h$d, "methylation", "methylation.coverage.tsv.gz"), check.names = FALSE)
  expect_identical(names(cov)[-1], meta$sample_id)
  arr <- run("Homo sapiens", "clinical", list("methylation_array"))
  expect_equal(arr$t$engine, "array")
})

test_that("genetics with no genotype factor (BSA pools) simulates the rest and lists genetics as not simulated", {
  d <- make_design_folder("crd", "y ~ treatment", list(design_factor("phenotype_group", 2)), reps = 3,
                          extra = list(assays = list("bsa_seq", "rnaseq")))
  suppressMessages(simitall_from_biodes(d, params = list(simulate_rnaseq_experiment = list(n_genes = 100))))
  t <- jsonlite::read_json(file.path(d, "sim_truth.json"))
  expect_null(t$genetics)
  expect_true("assay_rnaseq" %in% names(t$assays))
  ns <- vapply(t$not_simulated, function(z) z$status, "")
  expect_true("needs_genotype_factor" %in% ns)
  expect_false(file.exists(file.path(d, "genotypes.vcf")))
})
