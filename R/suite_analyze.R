#' Analyze a Biodes study with the model its design calls for
#'
#' Reads a Biodes design folder after the outcome has been filled in (by the
#' experimenter, or by [simitall_from_biodes()]) and fits the spec's
#' `analysis.model`: `lm()` for fixed-effect designs, `lmer()` when the model
#' has random terms (alpha lattice, split plot). Subsamples (several rows per
#' experimental unit) are averaged per unit or nested as a random effect;
#' treating them as independent replicates is refused, because Biodes already
#' decided what the experimental unit is. For count assays with a `counts.tsv`
#' it also tests each gene with a linear model on log2-CPM using the same
#' design terms.
#'
#' When `sim_truth.json` is present (simulated data), estimates are compared
#' with the true values: `recovery.tsv` (per coefficient: truth, estimate, 95%
#' CI covered) and, for counts, true/false positive rates.
#'
#' @param biodes_dir Biodes design folder with `primary_outcome` filled in.
#' @param outdir Where to write results (default: the same folder).
#' @param subsamples `"aggregate"` (mean per unit, default) or `"nest"`
#'   (unit as a random effect). `"as_replicates"` is refused.
#' @param alpha Significance level for adjusted p-values and CI coverage.
#' @return Invisibly, a list with the fit, results and paths written.
#' @export
bioana_fit_biodes <- function(biodes_dir, outdir = biodes_dir,
                              subsamples = c("aggregate", "nest", "as_replicates"), alpha = 0.05) {
  subsamples <- match.arg(subsamples)
  h <- .suite_read_biodes(biodes_dir)
  p <- h$parsed
  spec <- h$spec
  meta <- h$meta
  warn <- character()
  truth_path <- file.path(biodes_dir, "sim_truth.json")
  truth <- if (file.exists(truth_path)) jsonlite::read_json(truth_path, simplifyVector = FALSE) else NULL
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  paths <- c()
  out <- list()

  y <- suppressWarnings(as.numeric(meta$primary_outcome))
  fit_outcome <- any(!is.na(y))
  summary <- list(model_from_spec = h$model, design = spec$design$chosen,
                  experimental_unit = spec$experimental_unit %||% "unit", handoff_version = spec$handoff_version)

  if (fit_outcome) {
    meta$y_ <- y
    unit <- .suite_unit_id(meta)
    # subsamples: rows of one unit that agree on every model variable (repeated
    # measures over time or period differ in a model variable, so they are not subsamples)
    cell <- do.call(paste, c(list(unit), lapply(p$vars, function(v) meta[[p$col[[v]]]]), sep = "\r"))
    has_subs <- anyDuplicated(cell[!is.na(unit)]) > 0
    f <- p$formula
    f[[2]] <- as.name("y_")
    if (has_subs && subsamples == "as_replicates") {
      stop(sprintf("These rows are subsamples of %d %s units. Biodes set the experimental unit to '%s', so subsamples are not independent replicates. Use subsamples = 'aggregate' or 'nest'.",
                   length(unique(unit)), spec$experimental_unit %||% "experimental", spec$experimental_unit %||% "unit"), call. = FALSE)
    }
    dropped <- meta$sample_id[is.na(y)]
    dat <- meta[!is.na(y), , drop = FALSE]
    dcell <- cell[!is.na(y)]
    if (has_subs && subsamples == "aggregate") {
      first <- dat[!duplicated(dcell), , drop = FALSE]
      first$y_ <- as.numeric(tapply(dat$y_, dcell, mean)[dcell[!duplicated(dcell)]])
      dat <- first
    } else if (has_subs && subsamples == "nest") {
      dat$unit_id_ <- dcell
      f <- stats::update(f, . ~ . + (1 | unit_id_))
    }
    frame <- .suite_frame(dat, p)
    has_random <- length(.suite_findbars(f)) > 0
    n_fixed <- ncol(stats::model.matrix(stats::reformulate(c("1", p$fixed_labels)), data = frame))
    max_sat <- .suite_rules("handoff")$satterthwaite_max_fixed_coefficients %||% Inf
    big <- has_random && n_fixed > max_sat
    if (big) warn <- c(warn, sprintf("%d fixed coefficients: fitted with lme4 without Satterthwaite df (limit %s in handoff.json); p-values use residual df %d",
                                     n_fixed, max_sat, nrow(frame) - n_fixed))
    fit <- if (has_random) {
      if (!big && requireNamespace("lmerTest", quietly = TRUE)) lmerTest::lmer(f, data = frame) else lme4::lmer(f, data = frame)
    } else {
      stats::lm(f, data = frame)
    }
    if (has_random && lme4::isSingular(fit)) warn <- c(warn, "singular fit: a random-effect variance was estimated as ~0")
    co <- if (has_random) summary(fit)$coefficients else summary(fit)$coefficients
    pcol <- grep("^Pr\\(", colnames(co), value = TRUE)
    # degrees of freedom per coefficient: residual df (lm), Satterthwaite (lmerTest), else n - rank
    dfv <- if ("df" %in% colnames(co)) co[, "df"] else if (!has_random) rep(stats::df.residual(fit), nrow(co)) else rep(max(1, nrow(frame) - nrow(co)), nrow(co))
    pv <- if (length(pcol)) co[, pcol] else 2 * stats::pt(-abs(co[, 3]), dfv)
    res <- data.frame(term = rownames(co), estimate = co[, 1], se = co[, 2], df = dfv, statistic = co[, 3], p = pv,
                      stringsAsFactors = FALSE, row.names = NULL)
    res$role <- vapply(res$term, function(t) .suite_coef_role(t, p), "")
    res$padj <- NA_real_
    test <- res$role == "treatment"
    res$padj[test] <- stats::p.adjust(res$p[test], "BH")
    rpath <- file.path(outdir, "results.tsv")
    utils::write.table(res, rpath, sep = "\t", row.names = FALSE, quote = FALSE)
    paths <- c(paths, results = rpath)
    aov <- tryCatch(as.data.frame(stats::anova(fit)), error = function(e) NULL)
    summary <- c(summary, list(
      method = if (has_random) class(fit)[1] else "lm", model_fitted = paste(deparse(f), collapse = ""),
      response = spec$primary_outcome$name %||% "primary_outcome",
      n_rows = nrow(meta), n_units_fitted = nrow(frame), subsamples = if (has_subs) subsamples else "none",
      dropped_samples = dropped,
      anova = if (is.null(aov)) NULL else lapply(split(cbind(term = rownames(aov), aov), seq_len(nrow(aov))), as.list)))
    out$fit <- fit
    out$results <- res
    if (!is.null(truth)) {
      rec <- .suite_recovery(res, truth, frame, p, alpha)
      if (!is.null(rec)) {
        rec_path <- file.path(outdir, "recovery.tsv")
        utils::write.table(rec$table, rec_path, sep = "\t", row.names = FALSE, quote = FALSE)
        paths <- c(paths, recovery = rec_path)
        summary$recovery <- rec$summary
        out$recovery <- rec
      }
    }
  }

  cpath <- file.path(biodes_dir, "counts.tsv")
  if (file.exists(cpath)) {
    de <- .suite_fit_counts(cpath, meta, p, alpha)
    dpath <- file.path(outdir, "de_results.tsv")
    utils::write.table(de$table, dpath, sep = "\t", row.names = FALSE, quote = FALSE)
    paths <- c(paths, de_results = dpath)
    summary$counts <- de$summary
    warn <- c(warn, de$warnings)
    if (!is.null(truth$counts)) {
      summary$counts$recovery <- .suite_count_recovery(de$table, truth$counts, alpha)
    }
    out$de <- de$table
  }
  if (!fit_outcome && !file.exists(cpath)) {
    stop(sprintf("Nothing to analyze in %s: the primary_outcome column is empty and there is no counts.tsv. Run simitall_from_biodes(\"%s\") first to simulate data, or fill primary_outcome with your measurements.",
                 biodes_dir, biodes_dir), call. = FALSE)
  }
  summary$warnings <- warn
  paths <- c(paths, summary = .suite_write_json(summary, file.path(outdir, "fit_summary.json")))
  out$summary <- summary
  out$paths <- paths
  if (!is.null(summary$recovery)) {
    message(sprintf("bioana: %s fitted; recovered %d/%d treatment coefficients within their %d%% CI (bias %.3g, RMSE %.3g).",
                    summary$method, summary$recovery$covered, summary$recovery$n_terms, round(100 * (1 - alpha)),
                    summary$recovery$bias, summary$recovery$rmse))
  }
  invisible(out)
}

# Which model variable does a coefficient name belong to?
.suite_coef_role <- function(term, p) {
  if (term == "(Intercept)") return("intercept")
  parts <- strsplit(term, ":", fixed = TRUE)[[1]]
  roles <- vapply(parts, function(x) {
    v <- if (startsWith(x, "I(")) intersect(all.vars(str2lang(x)), p$vars) else p$vars[startsWith(x, p$vars)]
    if (length(v)) p$role[[v[which.max(nchar(v))]]] else "other"
  }, "")
  if (all(roles == "treatment")) "treatment" else if (any(roles == "block")) "block" else roles[1]
}

# True coefficients under the same parametrisation: fit the fixed treatment
# part of the model to the noise-free treatment means.
.suite_recovery <- function(res, truth, frame, p, alpha) {
  fv <- names(truth$treatment_effects)
  nv <- names(truth$numeric_effects)
  tv <- c(fv, nv)
  if (!length(tv)) return(NULL)
  grid <- unique(frame[tv])
  mu <- rep(as.numeric(truth$intercept), nrow(grid))
  for (v in fv) mu <- mu + unlist(truth$treatment_effects[[v]])[as.character(grid[[v]])]
  for (v in nv) mu <- mu + as.numeric(truth$numeric_effects[[v]]) * grid[[v]]
  for (it in names(truth$interaction_effects)) mu <- mu + as.numeric(truth$interaction_effects[[it]])
  labels <- p$fixed_labels[vapply(p$fixed_labels, function(l) all(all.vars(str2lang(l)) %in% tv), logical(1))]
  if (!length(labels)) return(NULL)
  grid$true_mean_ <- mu
  tf <- stats::reformulate(labels, response = "true_mean_")
  true_coef <- stats::coef(stats::lm(tf, data = grid))
  keep <- intersect(names(true_coef)[names(true_coef) != "(Intercept)"], res$term)
  i <- match(keep, res$term)
  tab <- data.frame(term = keep, true = as.numeric(true_coef[keep]), estimate = res$estimate[i], se = res$se[i],
                    df = res$df[i])
  tab$true[abs(tab$true) < 1e-12] <- 0
  tab$error <- tab$estimate - tab$true
  tab$covered <- abs(tab$error) <= stats::qt(1 - alpha / 2, tab$df) * tab$se   # t interval with the coefficient's df
  list(table = tab, summary = list(
    n_terms = nrow(tab), covered = sum(tab$covered), coverage = round(mean(tab$covered), 3),
    bias = signif(mean(tab$error), 4), rmse = signif(sqrt(mean(tab$error^2)), 4),
    cor = if (nrow(tab) > 2 && stats::sd(tab$true) > 0) round(stats::cor(tab$true, tab$estimate), 4) else NA,
    nominal = 1 - alpha))
}

# Per-gene test on log-CPM with the fixed design terms (the core of limma-voom,
# without precision weights): one linear model per gene, t-tests on the
# treatment coefficients. Random terms such as (1 | rep:block) are not used
# for counts. Genes below min_cpm in most samples are not tested.
.suite_fit_counts <- function(cpath, meta, p, alpha, min_cpm = 1) {
  counts <- utils::read.delim(cpath, check.names = FALSE)
  ids <- intersect(meta$sample_id, names(counts))
  m <- as.matrix(counts[, ids, drop = FALSE])
  rownames(m) <- counts$gene_id
  frame <- .suite_frame(meta[match(ids, meta$sample_id), , drop = FALSE], p)
  rhs <- stats::reformulate(p$fixed_labels)
  X <- stats::model.matrix(rhs, data = frame)
  lib <- colSums(m)
  logcpm <- log2(sweep(m + 0.5, 2, lib + 1, "/") * 1e6)
  cpm <- sweep(m, 2, lib, "/") * 1e6
  keep <- rowSums(cpm >= min_cpm) >= ceiling(ncol(m) / 2)
  Y <- t(logcpm[keep, , drop = FALSE])
  qrx <- qr(X)
  beta <- qr.coef(qrx, Y)
  resid <- Y - X %*% beta
  df <- nrow(X) - qrx$rank
  s2 <- colSums(resid^2) / df
  xtxi <- chol2inv(qr.R(qrx))
  treat <- names(p$role)[p$role == "treatment"][1]
  rows <- grep(paste0("^", treat), colnames(X))
  tabs <- lapply(rows, function(j) {
    se <- sqrt(s2 * xtxi[j, j])
    tt <- beta[j, ] / se
    data.frame(gene_id = colnames(Y), contrast = colnames(X)[j], log2fc = beta[j, ], se = se,
               p = 2 * stats::pt(-abs(tt), df), stringsAsFactors = FALSE, row.names = NULL)
  })
  tab <- do.call(rbind, tabs)
  tab$padj <- stats::ave(tab$p, tab$contrast, FUN = function(x) stats::p.adjust(x, "BH"))
  warnings <- if (any(!keep)) sprintf("%d low-count genes not tested (CPM < %g in most samples)", sum(!keep), min_cpm) else character()
  list(table = tab, warnings = warnings, summary = list(
    method = "per-gene linear model on log2-CPM with t-tests (limma-voom core, no precision weights)",
    design = paste(deparse(rhs), collapse = ""), residual_df = df,
    random_terms_dropped = vapply(p$random, `[[`, "", "label"),
    n_genes_tested = sum(keep), n_samples = length(ids),
    significant = as.list(tapply(tab$padj < alpha, tab$contrast, sum, na.rm = TRUE))))
}

.suite_count_recovery <- function(tab, ct, alpha) {
  out <- list()
  for (cond in names(ct$de_genes)) {
    contrast <- paste0(ct$condition_term, cond)
    t <- tab[tab$contrast == contrast, , drop = FALSE]
    true_de <- names(ct$de_genes[[cond]])
    called <- t$gene_id[!is.na(t$padj) & t$padj < alpha]
    est <- t$log2fc[match(true_de, t$gene_id)]
    tru <- unlist(ct$de_genes[[cond]])
    out[[cond]] <- list(
      contrast = contrast, true_de = length(true_de), called = length(called),
      true_positives = length(intersect(called, true_de)),
      tpr = round(length(intersect(called, true_de)) / max(length(true_de), 1), 3),
      fdr = round(length(setdiff(called, true_de)) / max(length(called), 1), 3),
      log2fc_cor = if (sum(!is.na(est)) > 2) round(stats::cor(est, tru, use = "complete.obs"), 3) else NA)
  }
  out
}
