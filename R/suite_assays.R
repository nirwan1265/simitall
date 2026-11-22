# Assay layer of the Biodes bridge: omics data for every experimental unit of
# the layout. Each library is named by the unit's sample_id, so the data line
# up with sample_metadata_template.tsv. Simulators and their parameters come
# from inst/suite/simulators.tsv and inst/suite/simulator_parameters.tsv.

# Condition term for assays: the first treatment term that is not the genotype
# and not numeric; else the genotype term; else none.
.suite_condition_term <- function(p, geno_term = NA_character_) {
  treat <- setdiff(names(p$role)[p$role == "treatment"], p$numeric)
  other <- setdiff(treat, geno_term)
  if (length(other)) other[1] else if (length(treat)) treat[1] else NA_character_
}

# One row per unit: sample (= sample_id), condition, batch; reference level first.
.suite_unit_design <- function(meta, frame, p, cond_term, extra = character()) {
  blk <- names(p$role)[p$role == "block"][1]
  cond <- if (is.na(cond_term)) factor(rep("all", nrow(meta))) else frame[[cond_term]]
  # uncovered factors (e.g. tissue) also shape the data: the assay condition is their combination
  extra <- setdiff(extra, cond_term)
  if (length(extra)) {
    parts <- c(if (!is.na(cond_term)) list(frame[[cond_term]]), lapply(extra, function(v) frame[[v]]))
    cond <- interaction(parts, sep = ".", lex.order = TRUE, drop = TRUE)
  }
  d <- data.frame(sample = meta$sample_id, condition = as.character(cond),
                  batch = if (!is.na(blk)) paste0("batch_", frame[[blk]]) else "batch1",
                  biological_replicates = 1L, technical_replicates = 1L, stringsAsFactors = FALSE)
  d <- d[order(as.integer(cond)), , drop = FALSE]
  attr(d, "levels") <- levels(droplevels(cond))
  attr(d, "batch_term") <- if (is.na(blk)) NULL else blk
  d
}

.suite_de_truth <- function(path) {
  ct <- utils::read.delim(path, stringsAsFactors = FALSE)
  if (!nrow(ct)) return(list())
  lapply(split(ct, ct$condition), function(x) stats::setNames(as.list(x$effect_log2), x$gene_id))
}

# counts.tsv for the assay that owns it, else <assay>/pseudobulk_counts.tsv (ctx$counts_file).
.suite_write_counts <- function(counts, ctx) {
  path <- ctx$counts_file %||% file.path(ctx$outdir, "counts.tsv")
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(counts, path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}

.suite_call <- function(fn, fixed, args) {
  args <- args[setdiff(names(args), names(fixed))]
  suppressMessages(do.call(fn, c(fixed, args)))
}

# Bulk RNA-seq: one biological replicate per unit.
.suite_assay_rnaseq <- function(ctx) {
  d <- .suite_unit_design(ctx$meta, ctx$frame, ctx$p, ctx$cond_term, ctx$uncovered)
  if (is.na(ctx$cond_term) && !length(ctx$uncovered)) stop("RNA-seq simulation needs a treatment term in analysis.model.", call. = FALSE)
  a <- .suite_args("simulate_rnaseq_experiment", ctx$spec, ctx$params)
  tmp <- tempfile("simitall_counts_"); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  res <- .suite_call(simulate_rnaseq_experiment,
                     list(out_prefix = file.path(tmp, "x"), replicate_mode = "biological", design = d, seed = ctx$seed), a$args)
  counts <- utils::read.delim(res$counts, check.names = FALSE)
  sm <- utils::read.delim(res$sample_metadata, stringsAsFactors = FALSE)
  names(counts)[-1] <- sm$sample[match(names(counts)[-1], sm$library_id)]
  counts <- counts[, c("gene_id", ctx$meta$sample_id)]
  path <- .suite_write_counts(counts, ctx)
  list(paths = c(counts = path), args = a,
       counts_truth = list(n_genes = nrow(counts), condition_term = ctx$cond_term, batch_term = attr(d, "batch_term"),
                           reference = attr(d, "levels")[1], de_fraction = a$args$condition_effect_fraction,
                           de_log2fc_sd = a$args$condition_effect_sd, batch_effect_sd = a$args$batch_effect_sd,
                           de_genes = .suite_de_truth(res$condition_truth)))
}

# eQTL: each unit carries its line's genotypes; cis/trans eQTL plus treatment effects.
.suite_assay_eqtl <- function(ctx) {
  if (is.null(ctx$genetics)) stop("eQTL simulation needs the genetics layer (a genotype factor in the design).", call. = FALSE)
  d <- .suite_unit_design(ctx$meta, ctx$frame, ctx$p, ctx$cond_term, ctx$uncovered)
  line_of <- as.character(ctx$frame[[ctx$geno_term]])[match(d$sample, ctx$meta$sample_id)]
  dos <- ctx$genetics$dosage[, line_of, drop = FALSE]
  colnames(dos) <- d$sample
  edir <- file.path(ctx$outdir, "eqtl"); dir.create(edir, showWarnings = FALSE)
  uvcf <- .suite_write_vcf(dos, ctx$genetics$variants, file.path(edir, "unit_genotypes.vcf"), "simitall unit genotypes (each unit = its line)")
  a <- .suite_args("simulate_rnaseq_from_gwas", ctx$spec, ctx$params)
  sm <- d[, c("sample", "condition", "batch")]
  sm$line <- line_of
  res <- .suite_call(simulate_rnaseq_from_gwas,
                     list(genotype_file = uvcf, out_prefix = file.path(edir, "eqtl"), sample_metadata = sm,
                          condition_levels = attr(d, "levels"), batch_levels = unique(d$batch), seed = ctx$seed), a$args)
  counts <- utils::read.delim(res$counts, check.names = FALSE)
  counts <- counts[, c("gene_id", ctx$meta$sample_id)]
  path <- .suite_write_counts(counts, ctx)
  et <- utils::read.delim(res$eqtl_truth, stringsAsFactors = FALSE)
  list(paths = c(counts = path, eqtl_truth = res$eqtl_truth, unit_genotypes = uvcf), args = a,
       counts_truth = list(n_genes = nrow(counts), condition_term = ctx$cond_term, batch_term = attr(d, "batch_term"),
                           reference = attr(d, "levels")[1], de_genes = .suite_de_truth(res$condition_truth)),
       truth = list(n_eqtl = nrow(et), n_cis = sum(et$type == "cis"), n_trans = sum(et$type == "trans"),
                    eqtl_truth = "eqtl/eqtl.eqtl_truth.tsv"))
}

# Single-cell RNA-seq: cells per unit; counts.tsv = pseudobulk (sum over cells) per unit.
.suite_assay_single_cell <- function(ctx) {
  d <- .suite_unit_design(ctx$meta, ctx$frame, ctx$p, ctx$cond_term, ctx$uncovered)
  a <- .suite_args("simulate_scrnaseq_experiment", ctx$spec, ctx$params)
  sdir <- file.path(ctx$outdir, "single_cell"); dir.create(sdir, showWarnings = FALSE)
  res <- .suite_call(simulate_scrnaseq_experiment,
                     list(out_prefix = file.path(sdir, "sc"), replicate_mode = "biological", design = d, seed = ctx$seed), a$args)
  pb <- utils::read.delim(res$pseudobulk_counts, check.names = FALSE)
  pm <- utils::read.delim(res$pseudobulk_metadata, stringsAsFactors = FALSE)
  unit <- pm$sample[match(names(pb)[-1], pm$pseudobulk_id)]
  m <- as.matrix(pb[, -1, drop = FALSE])
  summed <- t(rowsum(t(m), unit))
  counts <- data.frame(gene_id = pb$gene_id, summed[, ctx$meta$sample_id, drop = FALSE], check.names = FALSE)
  path <- .suite_write_counts(counts, ctx)
  cells <- utils::read.delim(res$cell_metadata, stringsAsFactors = FALSE)
  list(paths = c(counts = path, cell_metadata = res$cell_metadata, matrix = res$matrix, marker_truth = res$marker_truth), args = a,
       counts_truth = list(n_genes = nrow(counts), condition_term = ctx$cond_term, batch_term = attr(d, "batch_term"),
                           reference = attr(d, "levels")[1], pseudobulk = "sum of cells per unit",
                           de_genes = .suite_de_truth(res$condition_truth)),
       truth = list(n_cells = nrow(cells), cell_types = as.list(table(cells$cell_type)),
                    marker_truth = "single_cell/sc.marker_truth.tsv", cell_metadata = "single_cell/sc.cell_metadata.tsv"))
}

.suite_peak_genome <- function(ctx, dir) {
  a <- .suite_args("suite_peak_genome", ctx$spec, ctx$params)$args
  .simitall_simulate_peak_genome(file.path(dir, "genome.fa"), file.path(dir, "genome.gff3"),
                                 n_chromosomes = 1L, chromosome_names = "chr1",
                                 chromosome_length = as.integer(a$chromosome_length), n_genes = as.integer(a$n_genes),
                                 seed = ctx$seed)
}

# ChIP-seq: one ChIP library (plus input) per unit.
.suite_assay_chip_seq <- function(ctx) {
  d <- .suite_unit_design(ctx$meta, ctx$frame, ctx$p, ctx$cond_term, ctx$uncovered)
  cdir <- file.path(ctx$outdir, "chipseq"); dir.create(cdir, showWarnings = FALSE)
  g <- .suite_peak_genome(ctx, cdir)
  a <- .suite_args("simulate_chipseq", ctx$spec, ctx$params)
  res <- .suite_call(simulate_chipseq,
                     list(genome_fa = g$fasta, out_prefix = file.path(cdir, "chip"), annotation_gff3 = g$gff3,
                          conditions = attr(d, "levels"), design = d, seed = ctx$seed), a$args)
  pk <- utils::read.delim(res$peaks, stringsAsFactors = FALSE)
  dt <- utils::read.delim(res$differential_truth, stringsAsFactors = FALSE)
  list(paths = unlist(res), args = a,
       truth = list(n_peaks = nrow(pk), n_differential = nrow(dt), condition_term = ctx$cond_term,
                    reference = attr(d, "levels")[1], peaks = "chipseq/chip.truth_peaks.tsv",
                    differential = "chipseq/chip.differential_binding_truth.tsv",
                    library_to_unit = "chipseq/chip.sample_metadata.tsv (column sample = sample_id)"))
}

# ATAC-seq: simulate_atacseq() takes replicates per condition; libraries are mapped to units in order.
.suite_assay_atac_seq <- function(ctx) {
  d <- .suite_unit_design(ctx$meta, ctx$frame, ctx$p, ctx$cond_term, ctx$uncovered)
  adir <- file.path(ctx$outdir, "atacseq"); dir.create(adir, showWarnings = FALSE)
  g <- .suite_peak_genome(ctx, adir)
  a <- .suite_args("simulate_atacseq", ctx$spec, ctx$params)
  lev <- attr(d, "levels")
  per <- table(factor(d$condition, levels = lev))
  res <- .suite_call(simulate_atacseq,
                     list(genome_fa = g$fasta, out_prefix = file.path(adir, "atac"), annotation_gff3 = g$gff3,
                          conditions = lev, biological_replicates = as.integer(max(per)), seed = ctx$seed), a$args)
  lib <- utils::read.delim(res$sample_metadata, stringsAsFactors = FALSE)
  lib$sample_id <- NA_character_
  for (l in lev) {
    u <- d$sample[d$condition == l]
    i <- which(lib$condition == l)
    lib$sample_id[i[seq_along(u)]] <- u
  }
  mp <- file.path(adir, "library_to_unit.tsv")
  utils::write.table(lib, mp, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  pk <- utils::read.delim(res$peaks, stringsAsFactors = FALSE)
  dt <- utils::read.delim(res$differential_truth, stringsAsFactors = FALSE)
  list(paths = c(unlist(res), library_to_unit = mp), args = a,
       truth = list(n_peaks = nrow(pk), n_differential = nrow(dt), condition_term = ctx$cond_term,
                    reference = lev[1], peaks = "atacseq/atac.truth_peaks.tsv",
                    differential = "atacseq/atac.differential_accessibility_truth.tsv",
                    library_to_unit = "atacseq/library_to_unit.tsv",
                    unused_libraries = sum(is.na(lib$sample_id))))
}

.suite_assay_runner <- function(simulator) {
  switch(simulator,
         simulate_rnaseq_experiment = .suite_assay_rnaseq,
         simulate_rnaseq_from_gwas = .suite_assay_eqtl,
         simulate_scrnaseq_experiment = .suite_assay_single_cell,
         simulate_chipseq = .suite_assay_chip_seq,
         simulate_atacseq = .suite_assay_atac_seq,
         small_rnaseq_template = function(ctx) .suite_assay_template(ctx, "small_rnaseq", c("small_rnaseq", "mirna_seq", "srna_seq")),
         microbiome_template = function(ctx) .suite_assay_template(ctx, "microbiome", c("amplicon", "16s", "its", "amplicon_16s", "metagenomics", "microbiome")),
         stop("No runner for simulator ", simulator, call. = FALSE))
}

# Template for a count assay: same species, else same kingdom (inst/suite/templates.json), else NULL.
.suite_pick_template <- function(assay_values, spec, user = NULL) {
  tab <- simitall_templates()
  if (!is.null(user) && nzchar(user)) {
    if (!user %in% tab$id) stop("Unknown template '", user, "'. See simitall_templates().", call. = FALSE)
    return(list(id = user, why = "chosen in params"))
  }
  ok <- vapply(strsplit(tab$assays, "|", fixed = TRUE), function(a) any(assay_values %in% a), logical(1))
  tab <- tab[ok, , drop = FALSE]
  if (!nrow(tab)) return(NULL)
  rules <- .suite_rules("templates")
  org <- tolower(trimws(paste(unlist(spec$organism), collapse = " ")))
  same <- tab[nzchar(org) & tolower(tab$organism) == org, , drop = FALSE]
  if (nrow(same)) return(list(id = same$id[1], why = sprintf("same species (%s)", spec$organism)))
  genus <- sub(" .*", "", org)
  kingdom <- names(Filter(function(g) genus %in% g, rules$kingdom_by_genus))
  if (!length(kingdom)) kingdom <- names(Filter(function(s) any(.suite_norm(spec$system) %in% s), rules$kingdom_by_system))
  if (length(kingdom)) {
    k <- tab[tab$kingdom == kingdom[1], , drop = FALSE]
    if (nrow(k)) return(list(id = k$id[1], why = sprintf("nearest template: same kingdom (%s); no %s template", kingdom[1], if (nzchar(org)) spec$organism else "species-specific")))
  }
  NULL
}

# Count assays simulated from a reference template (small RNA-seq, 16S/ITS microbiome), with a
# native negative-binomial fallback where templates.json names one.
.suite_assay_template <- function(ctx, kind, assay_values) {
  d <- .suite_unit_design(ctx$meta, ctx$frame, ctx$p, ctx$cond_term, ctx$uncovered)
  a <- .suite_args("simulate_counts_from_template", ctx$spec, ctx$params)
  pick <- .suite_pick_template(assay_values, ctx$spec, a$args$template)
  odir <- file.path(ctx$outdir, kind); dir.create(odir, showWarnings = FALSE)
  if (is.null(pick)) {
    fb <- .suite_rules("templates")$native_fallback[[kind]]
    if (is.null(fb)) stop(sprintf("No %s template fits organism '%s' and there is no native fallback.", kind, ctx$spec$organism %||% "?"), call. = FALSE)
    na <- .suite_args(fb, ctx$spec, ctx$params)
    na$args$n_genes <- a$args$n_features_native %||% na$args$n_genes
    res <- .suite_call(get(fb, envir = asNamespace("simitall")),
                       list(out_prefix = file.path(odir, kind), replicate_mode = "biological", design = d, seed = ctx$seed),
                       na$args[setdiff(names(na$args), "n_features_native")])
    counts <- utils::read.delim(res$counts, check.names = FALSE)
    sm <- utils::read.delim(res$sample_metadata, stringsAsFactors = FALSE)
    names(counts)[-1] <- sm$sample[match(names(counts)[-1], sm$library_id)]
    tinfo <- list(template = NULL, mode = "native_negative_binomial",
                  why = sprintf("no %s template for organism '%s' (simitall_templates()); native negative binomial", kind, ctx$spec$organism %||% "?"))
    ct <- res$condition_truth
  } else {
    targs <- a$args[setdiff(names(a$args), c("template", "n_features_native"))]
    res <- .suite_call(simulate_counts_from_template,
                       list(template = pick$id, design = d, out_prefix = file.path(odir, kind), seed = ctx$seed), targs)
    counts <- utils::read.delim(res$counts, check.names = FALSE)
    tinfo <- c(res$truth[setdiff(names(res$truth), "template_sample_per_unit")], list(why = pick$why))
    ct <- res$condition_truth
  }
  counts <- counts[, c("gene_id", ctx$meta$sample_id)]
  path <- .suite_write_counts(counts, ctx)
  list(paths = c(counts = path, condition_truth = ct), args = a,
       counts_truth = list(n_genes = nrow(counts), condition_term = ctx$cond_term, batch_term = attr(d, "batch_term"),
                           reference = attr(d, "levels")[1], feature_type = kind, de_genes = .suite_de_truth(ct)),
       truth = c(list(template_choice = tinfo), list(condition_truth = sub(paste0("^", ctx$outdir, "/"), "", ct))))
}
