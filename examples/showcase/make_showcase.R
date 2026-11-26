# Simitall showcase: simulate the three Biodes showcase designs (seed 2026) and document what was planted.
#
#   Rscript examples/showcase/make_showcase.R            (from the repository root)
#
# Inputs: the Biodes folders examples/showcase/<name>/ (study_spec.json, field_layout.csv,
# sample_metadata_template.tsv) from BIODES_SHOWCASE, else ../BioDesignGPT/examples/showcase, else GitHub.
# Outputs, per folder: the Biodes files, Simitall's data and truth, summary.png and README.md.
# The same folder names continue in Bioana (analysis scored against this truth) and Biopaint (figures).
#
# Options for 03_germination_survival (environment variables):
#   SHOWCASE_INTERACTION = 0 (default) | model   no genotype x priming interaction, or the model's default-size one
#   SHOWCASE_PRIMING_HR  = (unset: the 0.6 typed in the question) | e.g. 1.67   hazard ratio of the priming effect

suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else library(simitall)
  library(ggplot2)
})
`%||%` <- function(x, y) if (is.null(x)) y else x

out_root <- file.path("examples", "showcase")
biodes_url <- "https://github.com/nirwan1265/BioDesignGPT/tree/main/examples/showcase"
bioana_url <- "https://github.com/nirwan1265/Bioana/tree/main/examples/showcase"
raw_url <- "https://raw.githubusercontent.com/nirwan1265/BioDesignGPT/main/examples/showcase"
src <- Sys.getenv("BIODES_SHOWCASE", file.path("..", "BioDesignGPT", "examples", "showcase"))
files <- c("study_spec.json", "field_layout.csv", "sample_metadata_template.tsv")

examples <- list(
  "01_gwas_sorghum" = list(params = list(suite_genetics = list(heritability = 0.6)),
                           params_text = "heritability 0.6 (`suite_genetics`)"),
  "02_rnaseq_drought" = list(params = list(simulate_rnaseq_experiment = list(n_genes = 5000)),
                             params_text = "5,000 genes (`simulate_rnaseq_experiment`)"),
  "03_germination_survival" = local({
    inter <- Sys.getenv("SHOWCASE_INTERACTION", "0"); hr <- Sys.getenv("SHOWCASE_PRIMING_HR", "")
    list(params = if (nzchar(hr)) list(suite_outcome = list(hazard_ratio = as.numeric(hr))) else list(),
         interaction_effect = if (identical(inter, "model")) NULL else as.numeric(inter),
         params_text = paste0("default settings except ",
                              if (identical(inter, "model")) "the model's genotype x priming interaction at its default size" else "no genotype x priming interaction",
                              if (nzchar(hr)) sprintf(" and a priming hazard ratio of %s instead of the 0.6 in the question", hr) else "",
                              " (see Choices below)"))
  })
)

fmt <- function(x, d = 2) formatC(as.numeric(x), digits = d, format = "f", big.mark = ",")
theme_show <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

for (nm in names(examples)) {
  ex <- examples[[nm]]
  d <- file.path(out_root, nm)
  unlink(d, recursive = TRUE); dir.create(d, recursive = TRUE)
  for (f in files) {
    from <- file.path(src, nm, f)
    if (file.exists(from)) file.copy(from, file.path(d, f), overwrite = TRUE)
    else utils::download.file(file.path(raw_url, nm, f), file.path(d, f), quiet = TRUE)
  }
  question <- if (file.exists(file.path(src, nm, "question.txt"))) paste(readLines(file.path(src, nm, "question.txt"), warn = FALSE), collapse = " ")
              else tryCatch(paste(readLines(file.path(raw_url, nm, "question.txt"), warn = FALSE), collapse = " "), error = function(e) "")

  simitall_from_biodes(d, seed = 2026, params = ex$params, interaction_effect = ex$interaction_effect)
  for (e in list.dirs(d, recursive = FALSE)) if (!length(list.files(e))) unlink(e, recursive = TRUE)   # empty truth/

  # keep the repository light: the genotype VCF is gzipped (Bioana reads .vcf.gz) and the manifest follows
  if (file.exists(file.path(d, "genotypes.vcf"))) {
    con <- gzfile(file.path(d, "genotypes.vcf.gz"), "w"); writeLines(readLines(file.path(d, "genotypes.vcf")), con); close(con)
    unlink(file.path(d, "genotypes.vcf"))
    for (j in c("data_manifest.json", "sim_truth.json")) {
      x <- readLines(file.path(d, j)); writeLines(gsub("genotypes.vcf\"", "genotypes.vcf.gz\"", x, fixed = TRUE), file.path(d, j))
    }
  }

  spec <- jsonlite::read_json(file.path(d, "study_spec.json"))
  truth <- jsonlite::read_json(file.path(d, "sim_truth.json"))
  meta <- utils::read.delim(file.path(d, "sample_metadata_template.tsv"), check.names = FALSE)
  lay <- utils::read.csv(file.path(d, "field_layout.csv"), check.names = FALSE)
  planted <- character(); fig_caption <- ""

  if (nm == "01_gwas_sorghum") {
    g <- truth$genetics
    q <- utils::read.delim(file.path(d, g$qtn_truth))
    gv <- utils::read.delim(file.path(d, g$genetic_values))
    # every marker as a tick on one genome axis, the QTNs as lollipops
    vl <- readLines(gzfile(file.path(d, "genotypes.vcf.gz")))
    vl <- strsplit(vl[!startsWith(vl, "#")], "\t", fixed = TRUE)
    vcf <- data.frame(chr = vapply(vl, `[`, "", 1), pos = as.numeric(vapply(vl, `[`, "", 2)))
    chrs <- unique(vcf$chr[order(as.integer(sub("\\D+", "", vcf$chr)))])
    len <- tapply(vcf$pos, vcf$chr, max)[chrs]; off <- stats::setNames(c(0, cumsum(as.numeric(len)))[seq_along(chrs)], chrs)
    vcf$x <- (vcf$pos + off[vcf$chr]) / 1e6; q$x <- (q$pos + off[q$seqname]) / 1e6
    mids <- (off + len / 2) / 1e6
    p <- ggplot(q, aes(x, effect)) +
      geom_rect(data = data.frame(xmin = off / 1e6, xmax = (off + len) / 1e6, k = seq_along(chrs) %% 2), inherit.aes = FALSE,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = factor(k)), alpha = 0.35, show.legend = FALSE) +
      scale_fill_manual(values = c("grey88", "white")) +
      geom_rug(data = vcf, aes(x = x), inherit.aes = FALSE, sides = "b", alpha = 0.15, length = grid::unit(0.03, "npc")) +
      geom_hline(yintercept = 0, colour = "grey60") +
      geom_segment(aes(xend = x, yend = 0), colour = "#2b6c8f") + geom_point(aes(size = maf), colour = "#2b6c8f") +
      scale_x_continuous(breaks = mids, labels = chrs, expand = c(0.01, 0)) + scale_size_continuous(range = c(1.5, 4), name = "MAF") +
      labs(title = sprintf("The %d planted QTNs (heritability %s); %s markers along the bottom", nrow(q), g$heritability, fmt(g$n_markers, 0)),
           x = NULL, y = "additive effect (cm per allele)") + theme_show + theme(panel.grid.major.x = element_blank())
    ggsave(file.path(d, "summary.png"), p, width = 9, height = 3.6, dpi = 150)
    fig_caption <- "The planted QTNs: genome position, additive effect per allele and minor allele frequency; the ticks at the bottom are the markers in `genotypes.vcf.gz`."
    top <- q[order(-abs(q$effect)), ][1:3, ]
    planted <- c(
      sprintf("- **Genotypes:** %d inbred accessions x %s markers from a %s with %d subpopulations (F_ST %s), in `genotypes.vcf.gz`.",
              g$n_lines, fmt(g$n_markers, 0), sub("_", " ", g$genotype_source$source), g$genotype_source$n_subpops, g$genotype_source$fst),
      sprintf("- **Genetic architecture:** %d additive QTNs, effects from a %s distribution; heritability %s (genetic variance %s cm² against a plot residual SD of %s cm). Largest: %s.",
              g$n_add_qtn, g$effect_distribution, g$heritability, fmt(g$genetic_variance, 0), fmt(truth$residual_sd, 0),
              paste(sprintf("%s (%+.1f cm, MAF %.2f)", top$marker_id, top$effect, top$maf), collapse = ", ")),
      sprintf("- **Plots:** plant height = %s cm + the accession's genetic value + replicate and incomplete-block effects (block SD %s cm) + residual (SD %s cm), %d plots.",
              fmt(truth$intercept, 0), fmt(unlist(truth$random_sd)[1], 1), fmt(truth$residual_sd, 0), nrow(meta)),
      "- **Truth files:** `genetic_truth/qtn_truth.tsv` (marker, position, effect, MAF), `genetic_truth/genetic_values.tsv` (each accession's true value), `genetic_truth/population_structure.tsv`.")
  }

  if (nm == "02_rnaseq_drought") {
    a <- truth$assays$assay_rnaseq
    de <- utils::read.delim(file.path(d, "truth", "rnaseq.tsv"))
    gi <- utils::read.delim(file.path(d, truth$goal_truth$features$rnaseq$file))
    cm <- utils::read.delim(file.path(d, "counts.tsv"), check.names = FALSE)
    m <- as.matrix(cm[, -1]); rownames(m) <- cm[[1]]
    lc <- log2(sweep(m + 0.5, 2, colSums(m), "/") * 1e6)
    smp <- meta[match(colnames(m), meta$sample_id), ]
    wl <- spec$factors[[1]]$levels; gl <- spec$factors[[2]]$levels
    cell <- function(w, gn) rowMeans(lc[, smp$water == w & smp$genotype == gn, drop = FALSE])
    fc_wt <- cell(wl[[2]], gl[[1]]) - cell(wl[[1]], gl[[1]]); fc_mu <- cell(wl[[2]], gl[[2]]) - cell(wl[[1]], gl[[2]])
    cls <- ifelse(rownames(m) %in% gi$feature_id, "interaction (planted)", ifelse(rownames(m) %in% de$feature_id[abs(de$effect) > 0], "drought effect only", "no planted effect"))
    df <- data.frame(wt = fc_wt, mu = fc_mu, cls = factor(cls, c("no planted effect", "drought effect only", "interaction (planted)")))
    p <- ggplot(df[order(df$cls), ], aes(wt, mu, colour = cls)) + geom_abline(colour = "grey50", linetype = 2) +
      geom_point(size = 0.8, alpha = 0.7) + scale_colour_manual(values = c("grey80", "#2b6c8f", "#c0392b"), name = NULL) +
      labs(title = "Drought response, wild type vs mutant",
           x = sprintf("log2 FC drought, %s", gl[[1]]), y = sprintf("log2 FC drought, %s", gl[[2]])) +
      coord_equal() + theme_show + theme(legend.position = "bottom")
    ggsave(file.path(d, "summary.png"), p, width = 6.2, height = 6.4, dpi = 150)
    fig_caption <- "Each point is a gene: its observed log2 fold change under drought in each genotype. On the dashed line, drought acts the same in both genotypes; the planted interaction genes (red) sit off it."
    planted <- c(
      sprintf("- **Counts:** %s genes x %d samples (%d plants per genotype x water cell), library size about %s M reads, negative binomial counts (`counts.tsv`).",
              fmt(nrow(m), 0), ncol(m), as.integer(spec$design$reps), fmt(a$parameters$library_size_mean / 1e6, 0)),
      sprintf("- **Drought effect:** %d genes (%s%%) change between well-watered and drought in both genotypes, log2 fold changes with SD %s (`truth/rnaseq.tsv`).",
              sum(abs(de$effect) > 0), fmt(100 * a$parameters$condition_effect_fraction, 0), a$parameters$condition_effect_sd),
      sprintf("- **The goal, genotype x drought interaction:** %d genes (%s%%) whose drought response differs in the mutant. Each has a difference of differences (drought minus well-watered in the mutant, minus the same in wild type) of %s to %s log2 units, either sign (`truth/rnaseq_goal.tsv`: `feature_id, term, effect`).",
              nrow(gi), fmt(100 * nrow(gi) / nrow(m), 0), fmt(min(abs(gi$effect)), 1), fmt(max(abs(gi$effect)), 1)),
      sprintf("- **Blocks:** a block effect per replicate (SD %s on the log scale); no hidden factors (`n_latent = 0`).", a$parameters$batch_effect_sd),
      "- Every other gene has no planted effect: these give the false-positive rate.")
  }

  if (nm == "03_germination_survival") {
    om <- truth$outcome_model
    meta$grp <- paste(meta$genotype %||% lay$genotype[match(meta$plot, lay$plot)], meta$treatment, sep = " | ")
    if (!requireNamespace("survival", quietly = TRUE)) stop("the survival package is needed for this figure")
    sf <- survival::survfit(survival::Surv(primary_outcome, event == 1) ~ grp, data = meta)
    sm <- summary(sf, times = sort(unique(c(0, meta$primary_outcome))), extend = TRUE)
    df <- data.frame(time = sm$time, germinated = 1 - sm$surv, grp = sub("^grp=", "", sm$strata))
    tl <- do.call(rbind, lapply(names(om$groups), function(gn) data.frame(time = unlist(om$landmarks),
           germinated = 1 - unlist(om$groups[[gn]]$survival_at) - (if (!is.null(om$competing)) 0 else 0), grp = sub(":", " | ", gn))))
    p <- ggplot(df, aes(time, germinated, colour = grp)) + geom_step(linewidth = 0.8) +
      geom_vline(xintercept = om$change_time, linetype = 3, colour = "grey40") +
      annotate("text", x = om$change_time, y = 0.98, label = sprintf("priming acts from day %.0f", om$change_time), hjust = -0.05, size = 3.2, colour = "grey30") +
      scale_colour_manual(values = c("#9ecae1", "#2b6c8f", "#fdae6b", "#c0392b"), name = "genotype | treatment") +
      labs(title = "Germination over 21 days (1 - Kaplan-Meier; rot and censoring as censored)", x = "day", y = "share germinated") +
      coord_cartesian(ylim = c(0, 1)) + theme_show + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
    ggsave(file.path(d, "summary.png"), p, width = 7.4, height = 4.8, dpi = 150)
    fig_caption <- "Observed germination curves per genotype and priming; the dotted line is where the planted priming effect starts."
    n <- nrow(meta); ev <- table(factor(meta$event, 0:2))
    gtab <- vapply(names(om$groups), function(gn) { x <- om$groups[[gn]]
      sprintf("| %s | %s | %s | %s | %s |", sub(":", " | ", gn), fmt(x$median, 1), fmt(x$cif_event_at[[2]] %||% NA, 2), fmt(x$survival_at[[2]], 2), fmt(x$rmst, 1)) }, "")
    planted <- c(
      sprintf("- **Seeds:** %d seeds = %d tray x genotype x priming cells x %d seeds; each seed has its own time, scored every %s days up to day %s.",
              n, nrow(lay), as.integer(spec$subsamples_per_unit %||% 10), om$life_table_width, om$follow_up),
      sprintf("- **Model:** Weibull hazard (shape %s) with a delayed effect: hazard ratio 1 until day %s, then %s for primed vs not primed seeds (%s; the delay is `logrank_weight: fh_0_1`). The second line gets the same hazard ratio against the first%s.",
              om$weibull_shape, fmt(om$change_time, 1), om$hazard_ratio,
              if (grepl("^params", om$hazard_ratio_source)) "set with the SHOWCASE_PRIMING_HR option instead of the question's 0.6"
              else sprintf("the HR stated in the question, and in the spec as `power.target_hazard_ratio`; %s", if (om$hazard_ratio < 1) "below 1, so primed seeds germinate more slowly once the effect starts" else "above 1, so primed seeds germinate faster once the effect starts"),
              if (identical(ex$interaction_effect, 0)) "; no genotype x priming interaction" else ", plus the model's genotype x priming interaction"),
      sprintf("- **Germination speed:** the reference group (genotype 1, not primed) has a median of %s days (set so 20%% of it has not germinated by day %s).", fmt(om$median_time, 1), om$follow_up),
      sprintf("- **Seed rot (competing event, code 2):** its own hazard with no treatment effect, median %s days if nothing else happened. Outcome: %d germinated, %d rotted (%s%% of seeds), %d still waiting at day %s (censored).",
              fmt(om$competing$median_alone, 0), ev[["1"]], ev[["2"]], fmt(100 * ev[["2"]] / n, 0), ev[["0"]], om$follow_up),
      sprintf("- **Clusters:** trays share a log-normal frailty (SD %s on the log hazard), and so do the 10 seeds of one tray cell (SD %s).",
              fmt(om$frailty$block$sd_log_hazard, 2), fmt(om$frailty$plot$sd_log_hazard, 2)),
      "", "True values per group (conditional on an average tray), from `sim_truth.json` `outcome_model.groups`:", "",
      "| genotype | treatment | median, germination alone (days) | germinated by day 14 (CIF) | not germinated or rotted by day 14 | RMST to day 21 (days) |",
      "|---|---|---|---|---|---|", gtab, "",
      "### Choices made for this example, and how to change them", "",
      sprintf("- **%s.** The Biodes model is `treatment * treatment2`, so by default Simitall would plant a genotype x priming interaction the question never asked about. %s Option: `SHOWCASE_INTERACTION=model Rscript examples/showcase/make_showcase.R` plants it at the default size (`interaction_effect` in `simitall_from_biodes()`).",
              if (identical(ex$interaction_effect, 0)) "No interaction" else "Interaction planted",
              if (identical(ex$interaction_effect, 0)) "It is turned off here (`interaction_effect = 0`) so the story stays \"a delayed priming effect\"." else "Here it is planted at the model's default size."),
      sprintf("- **HR %s %s germination.** A hazard ratio below 1 means primed seeds germinate *more slowly* once the effect starts; above 1, faster. %s If priming should speed germination, the question should say HR 1.67, or run `SHOWCASE_PRIMING_HR=1.67 Rscript examples/showcase/make_showcase.R` (`params = list(suite_outcome = list(hazard_ratio = 1.67))`, which overrides the spec).",
              om$hazard_ratio, if (om$hazard_ratio < 1) "slows" else "speeds",
              if (grepl("^params", om$hazard_ratio_source)) "Here the HR comes from the option, not the question." else "0.6 is what the question typed, so that is what is planted."),
      "- **Survival fixes found while building this.** The survival goal's note in `sim_truth.json` still described the old log-normal model, and a stated hazard ratio raised a misleading \"no data on this scale\" note. Both were fixed in Simitall before these files were written.")
  }

  # README
  rel <- function(f) sprintf("[`%s`](%s)", f, f)
  data_files <- setdiff(list.files(d, recursive = TRUE), c(files, "README.md", "summary.png"))
  readme <- c(
    sprintf("# %s: simulated data", switch(nm, "01_gwas_sorghum" = "GWAS of plant height in a sorghum panel",
                                             "02_rnaseq_drought" = "RNA-seq: genotype x drought in Arabidopsis",
                                             "03_germination_survival" = "Seed priming and time to germination")),
    "", "**The question (typed into Biodes):**", "", paste0("> ", question), "",
    sprintf("Biodes turned it into %s `%s` design with the model `%s` ([Biodes folder](%s/%s)). Simitall simulated data for that exact layout, with seed 2026 and %s, and recorded everything it planted in [`sim_truth.json`](sim_truth.json).",
            if (grepl("^[aeiou]", spec$design$chosen)) "an" else "a", spec$design$chosen, spec$analysis$model, biodes_url, nm, ex$params_text),
    "", "## What was planted", "", planted, "",
    "## Summary figure", "", "![summary](summary.png)", "", fig_caption,
    "A quick look drawn by `make_showcase.R`; the full figure set for each example belongs to Biopaint.", "",
    "## Files", "",
    paste0("- From Biodes: ", paste(rel(files), collapse = ", ")),
    paste0("- From Simitall: ", paste(rel(data_files), collapse = ", ")),
    "- [`data_manifest.json`](data_manifest.json) tells Bioana which file holds which data and where its truth is.",
    "", "Assumptions Simitall made where the spec was silent (also in `sim_truth.json`):", "",
    paste0("- ", unlist(truth$assumptions)),
    "", "## Continue this example", "",
    sprintf("- **Before:** Biodes wrote the design: [BioDesignGPT/examples/showcase/%s](%s/%s)", nm, biodes_url, nm),
    sprintf("- **Next:** Bioana analyses these files with the model above and scores the result against the planted truth: [Bioana/examples/showcase/%s](%s/%s)", nm, bioana_url, nm),
    sprintf("- **Then:** Biopaint draws the figures: Biopaint/examples/showcase/%s (in progress)", nm),
    "", sprintf("Rebuild: `Rscript examples/showcase/make_showcase.R` from the repository root (simitall %s).", as.character(utils::packageVersion("simitall"))))
  writeLines(readme, file.path(d, "README.md"))
  message("wrote ", d)
}
