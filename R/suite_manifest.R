# data_manifest.json: what simulated data the folder holds, for Bioana.
#
# Bioana (github.com/nirwan1265/bioana) reads it to know which analysis to run
# on which files: one entry per data set with its data_type (a row of Bioana's
# inst/rules/analyses.tsv), its files, and its truth in Bioana's standard
# format (features: feature_id, condition, effect; regions: region_id, chr,
# start, end, condition, effect), written under truth/. An experimenter
# writes the same manifest for measured data (bioana::bioana_write_manifest()).

.suite_write_manifest <- function(outdir, truth, has_outcome, outcome_type, cond_term) {
  rel <- function(x) sub(paste0("^", normalizePath(outdir, mustWork = FALSE), "/"), "", normalizePath(x, mustWork = FALSE))
  tdir <- file.path(outdir, "truth"); dir.create(tdir, showWarnings = FALSE)
  wtruth <- function(df, name) {
    p <- file.path(tdir, paste0(name, ".tsv"))
    utils::write.table(df, p, sep = "\t", row.names = FALSE, quote = FALSE)
    file.path("truth", paste0(name, ".tsv"))
  }
  from_de <- function(ct) {
    rows <- lapply(names(ct$de_genes), function(cond) {
      g <- unlist(ct$de_genes[[cond]])
      if (!length(g)) return(NULL)
      data.frame(feature_id = names(g), condition = cond, effect = as.numeric(g), stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    if (is.null(out)) data.frame(feature_id = character(), condition = character(), effect = numeric()) else out
  }
  # the analysis goal's features (sim_truth.json goal_truth.features): feature_id, term, effect, change
  goal_entry <- function(gf) if (!is.null(gf)) list(goal = truth$goal_truth$goal, file = gf$file, term = gf$term, effect_scale = gf$effect_scale)
  rd <- function(f) utils::read.delim(file.path(outdir, f), stringsAsFactors = FALSE, check.names = FALSE)
  data <- list()
  add <- function(e) data[[length(data) + 1]] <<- e
  if (isTRUE(has_outcome)) {
    e <- list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"),
              column = "primary_outcome", outcome_type = outcome_type, truth = list(kind = "outcome", file = "sim_truth.json",
              goal = if (identical(truth$goal_truth$planted, "outcome")) list(goal = truth$goal_truth$goal, field = "goal_truth",
                                                                              terms = truth$goal_truth$model_terms, test = truth$goal_truth$test)))
    if (identical(outcome_type, "survival")) e$event_column <- "event"
    add(e)
  }
  a <- truth$assays %||% list()
  # count matrices: counts.tsv belongs to one assay; other count assays keep <assay>/unit_counts.tsv
  for (route in intersect(c("assay_eqtl", "assay_rnaseq", "assay_single_cell", "assay_small_rnaseq", "assay_microbiome"), names(a))) {
    id <- sub("^assay_", "", route)
    ct <- a[[route]]$counts %||% if (!is.null(truth$counts)) c(list(file = "counts.tsv"), truth$counts)
    if (is.null(ct)) next
    add(list(id = id, data_type = "feature_counts", assay = id, files = list(counts = ct$file %||% "counts.tsv"),
             note = if (route == "assay_single_cell") "pseudobulk: cells summed per unit",
             truth = list(kind = "features", file = wtruth(from_de(ct), id), term = ct$condition_term %||% cond_term, large = 1,
                          goal = goal_entry(truth$goal_truth$features[[id]]))))
  }
  if (!is.null(a$assay_chip_seq)) {
    x <- a$assay_chip_seq
    dt <- rd(x$differential)
    add(list(id = "chip_seq", data_type = "peak_counts_long", assay = "chip_seq",
             files = list(counts = "chipseq/chip.peak_counts.tsv", libraries = "chipseq/chip.sample_metadata.tsv"),
             truth = list(kind = "features", term = x$condition_term %||% cond_term, large = 1,
                          file = wtruth(data.frame(feature_id = dt$peak_id, condition = dt$condition, effect = dt$effect_log2), "chip_seq"))))
  }
  if (!is.null(a$assay_atac_seq)) {
    x <- a$assay_atac_seq
    pk <- rd("atacseq/atac.truth_peaks.tsv"); dt <- rd("atacseq/atac.differential_accessibility_truth.tsv")
    i <- match(dt$peak_id, pk$peak_id)
    add(list(id = "atac_seq", data_type = "read_positions", assay = "atac_seq",
             files = list(reads = "atacseq/atac.read_positions.tsv", libraries = "atacseq/library_to_unit.tsv"),
             truth = list(kind = "regions", term = x$condition_term %||% cond_term, large = 1,
                          file = wtruth(data.frame(region_id = dt$peak_id, chr = pk$seqname[i], start = pk$start[i], end = pk$end[i],
                                                   condition = dt$condition, effect = dt$effect_log2), "atac_seq"))))
  }
  if (!is.null(a$assay_methylation)) {
    x <- a$assay_methylation; f <- x$files
    dt <- rd(x$dmr_truth %||% "methylation/methylation.dmr_truth.tsv")
    if (!is.null(f$beta)) {
      add(list(id = "methylation", data_type = "methylation_beta", assay = "methylation_array", files = list(beta = f$beta),
               truth = list(kind = "features", term = cond_term,
                            file = wtruth(data.frame(feature_id = dt$feature_id, condition = dt$condition, effect = dt$delta), "methylation"))))
    } else {
      add(list(id = "methylation", data_type = "bisulfite_counts", assay = "methylation",
               files = list(sites = f$sites, methylated = f$methylated, coverage = f$coverage),
               truth = list(kind = "regions", term = cond_term,
                            file = wtruth(data.frame(region_id = dt$dmr_id, chr = dt$chr, start = dt$start, end = dt$end,
                                                     condition = dt$condition, effect = dt$delta), "methylation"))))
    }
  }
  if (!is.null(a$assay_m6a)) {
    f <- a$assay_m6a$files
    pt <- rd(f$peak_truth)
    cc <- grep("^log2_ratio_change_", names(pt), value = TRUE)
    long <- do.call(rbind, lapply(cc, function(k) data.frame(feature_id = pt$window_id, condition = sub("^log2_ratio_change_", "", k),
                                                             effect = pt[[k]], stringsAsFactors = FALSE)))
    long <- if (is.null(long)) data.frame(feature_id = character(), condition = character(), effect = numeric()) else long[long$effect != 0, , drop = FALSE]
    add(list(id = "m6a", data_type = "ip_input_counts", assay = "m6a_seq", files = list(ip = f$ip, input = f$input, windows = f$windows),
             truth = list(kind = "features", term = cond_term, file = wtruth(long, "m6a"))))
  }
  if (!is.null(truth$genetics)) {
    add(list(id = "genetics", data_type = "genotypes", assay = "genotypes", files = list(vcf = truth$genetics$genotypes),
             genotype_term = truth$genetics$genotype_term,
             truth = list(kind = "qtn", file = truth$genetics$qtn_truth)))
  }
  m <- list(manifest_version = "0.1", written_by = paste("simitall", as.character(utils::packageVersion("simitall"))),
            note = "Data sets in this folder for Bioana (github.com/nirwan1265/bioana): data_type is a row of its inst/rules/analyses.tsv; truth/ holds the simulated truth in Bioana's format.",
            goal = if (!is.null(truth$goal_truth)) list(goal = truth$goal_truth$goal, test = truth$goal_truth$test, status = truth$goal_truth$status,
                                                         truth = "sim_truth.json goal_truth"),
            data = data)
  p <- file.path(outdir, "data_manifest.json")
  jsonlite::write_json(m, p, auto_unbox = TRUE, pretty = TRUE, null = "null")
  p
}
