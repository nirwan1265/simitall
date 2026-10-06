# Single-cell RNA-seq analysis scored against simulation truth: clustering
# (does it recover the true cell types?) and condition differential
# expression (pseudobulk per donor versus per-cell tests, which treat cells
# as independent replicates and inflate false positives).

.simitall_scrna_counts <- function(object) {
  if (is.character(object)) object <- .simitall_read_scrna_object(object)
  counts <- object$counts
  if (is.null(rownames(counts))) rownames(counts) <- object$gene_metadata$gene_id
  if (is.null(colnames(counts))) colnames(counts) <- object$cell_metadata$cell_id
  list(counts = counts, cells = object$cell_metadata)
}

.simitall_scrna_lognorm <- function(counts) {
  lib <- Matrix::colSums(counts)
  log1p(t(t(as.matrix(counts)) / pmax(lib, 1)) * 1e4)
}

#' Cluster single cells and report the clusters
#'
#' Log-normalizes counts, keeps the most variable genes, runs PCA, and
#' clusters cells with k-means on the leading PCs. The number of clusters is
#' chosen by average silhouette width over a range of k.
#'
#' @param object simitall single-cell object or its prefix.
#' @param n_pcs Number of principal components used for clustering.
#' @param n_hvg Number of highly variable genes; `NULL` keeps 20% of genes
#'   (between 50 and 2,000).
#' @param k_range Candidate cluster numbers.
#' @param seed Random seed.
#' @return A list with `clusters` (data frame of cell_id and cluster), `pcs`,
#'   and `chosen_k`.
#' @export
cluster_scrnaseq <- function(object, n_pcs = 20L, n_hvg = NULL, k_range = 2:10, seed = 1L) {
  x <- .simitall_scrna_counts(object)
  logn <- .simitall_scrna_lognorm(x$counts)
  # Keep a minority of genes, as standard pipelines do (~2,000 of 20,000):
  # this enriches for cell-type markers over genome-wide variation.
  if (is.null(n_hvg)) n_hvg <- min(2000L, max(50L, round(0.2 * nrow(logn))))
  vars <- apply(logn, 1L, stats::var)
  keep <- order(vars, decreasing = TRUE)[seq_len(min(n_hvg, sum(vars > 0)))]
  scaled <- scale(t(logn[keep, , drop = FALSE]))
  scaled[!is.finite(scaled)] <- 0
  n_pcs <- max(2L, min(as.integer(n_pcs), ncol(scaled) - 1L, nrow(scaled) - 1L))
  pcs <- stats::prcomp(scaled, center = FALSE, rank. = n_pcs)$x
  set.seed(seed)
  sub <- if (nrow(pcs) > 1000L) sample(nrow(pcs), 1000L) else seq_len(nrow(pcs))
  d <- stats::dist(pcs[sub, , drop = FALSE])
  fits <- lapply(k_range[k_range < nrow(pcs)], function(k) stats::kmeans(pcs, centers = k, nstart = 10, iter.max = 50))
  dm <- as.matrix(d)
  widths <- vapply(fits, function(fit) .simitall_mean_silhouette(fit$cluster[sub], dm), numeric(1))
  best <- if (all(is.na(widths))) 1L else which.max(widths)
  list(clusters = data.frame(cell_id = rownames(pcs) %||% x$cells$cell_id, cluster = paste0("cluster_", fits[[best]]$cluster),
                             stringsAsFactors = FALSE),
       pcs = pcs, chosen_k = k_range[best], silhouette = stats::setNames(widths, k_range[seq_along(widths)]))
}

# Mean silhouette width from a full distance matrix (no extra packages).
.simitall_mean_silhouette <- function(labels, dm) {
  groups <- unique(labels)
  if (length(groups) < 2L) return(NA_real_)
  widths <- vapply(seq_along(labels), function(i) {
    same <- labels == labels[i]; same[i] <- FALSE
    if (!any(same)) return(0)
    a <- mean(dm[i, same])
    b <- min(vapply(setdiff(groups, labels[i]), function(g) mean(dm[i, labels == g]), numeric(1)))
    (b - a) / max(a, b)
  }, numeric(1))
  mean(widths)
}

.simitall_adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  comb <- function(n) n * (n - 1) / 2
  sum_ij <- sum(comb(tab)); sum_a <- sum(comb(rowSums(tab))); sum_b <- sum(comb(colSums(tab)))
  expected <- sum_a * sum_b / comb(sum(tab))
  max_index <- (sum_a + sum_b) / 2
  if (max_index == expected) return(1)
  (sum_ij - expected) / (max_index - expected)
}

.simitall_nmi <- function(a, b) {
  tab <- table(a, b) / length(a)
  pa <- rowSums(tab); pb <- colSums(tab)
  mi <- sum(tab[tab > 0] * log(tab[tab > 0] / outer(pa, pb)[tab > 0]))
  h <- function(p) -sum(p[p > 0] * log(p[p > 0]))
  if (h(pa) + h(pb) == 0) return(1)
  2 * mi / (h(pa) + h(pb))
}

#' Score clusters against the true cell types
#'
#' @param clusters Output of cluster_scrnaseq() (or its `clusters` table).
#' @param cell_metadata Cell metadata with `cell_id` and `cell_type`.
#' @return One-row data frame with the adjusted Rand index and NMI.
#' @export
score_scrnaseq_clusters <- function(clusters, cell_metadata) {
  if (is.list(clusters) && !is.data.frame(clusters)) clusters <- clusters$clusters
  truth <- cell_metadata$cell_type[match(clusters$cell_id, cell_metadata$cell_id)]
  data.frame(analysis = "clustering", true_cell_types = length(unique(truth)),
             clusters_found = length(unique(clusters$cluster)),
             adjusted_rand_index = .simitall_adjusted_rand(clusters$cluster, truth),
             nmi = .simitall_nmi(clusters$cluster, truth))
}

#' Test genes for condition differences in single-cell data
#'
#' Pseudobulk methods sum counts per donor and test donors as replicates
#' (edgeR, DESeq2, or a t-test on log-CPM). The per-cell Wilcoxon test treats
#' every cell as a replicate, which ignores donor variation and inflates false
#' positives; it is included as a contrast.
#'
#' @param object simitall single-cell object or its prefix.
#' @param method `"pseudobulk_deseq2"`, `"pseudobulk_edger"`,
#'   `"pseudobulk_ttest"`, or `"cell_wilcoxon"`.
#' @param fdr BH-adjusted threshold.
#' @return Data frame with gene_id, log2_fold_change, p_value, q_value, significant.
#' @export
test_scrnaseq_conditions <- function(object, method = c("pseudobulk_deseq2", "pseudobulk_edger", "pseudobulk_ttest", "cell_wilcoxon"),
                                     fdr = 0.05) {
  method <- match.arg(method)
  x <- .simitall_scrna_counts(object)
  cells <- x$cells[match(colnames(x$counts), x$cells$cell_id), , drop = FALSE]
  levels <- unique(cells$condition)
  if (length(levels) < 2L) stop("Condition testing needs two conditions")
  levels <- levels[1:2]
  if (identical(method, "cell_wilcoxon")) {
    logn <- .simitall_scrna_lognorm(x$counts)
    a <- cells$condition == levels[1L]; b <- cells$condition == levels[2L]
    stat <- t(apply(logn, 1L, function(v) {
      p <- if (stats::var(v[a | b]) == 0) 1 else suppressWarnings(stats::wilcox.test(v[b], v[a])$p.value)
      c(log2(mean(expm1(v[b])) + 1) - log2(mean(expm1(v[a])) + 1), p)
    }))
    result <- data.frame(gene_id = rownames(logn), log2_fold_change = stat[, 1L], p_value = stat[, 2L], stringsAsFactors = FALSE)
  } else {
    donors <- unique(cells$donor)
    pb <- vapply(donors, function(d) Matrix::rowSums(x$counts[, cells$donor == d, drop = FALSE]), numeric(nrow(x$counts)))
    rownames(pb) <- rownames(x$counts)
    group <- factor(cells$condition[match(donors, cells$donor)], levels = levels)
    keep <- !is.na(group); pb <- pb[, keep, drop = FALSE]; group <- droplevels(group[keep])
    if (min(table(group)) < 2L) stop("Pseudobulk testing needs at least two donors per condition")
    pb <- round(pb)
    result <- switch(method,
      pseudobulk_ttest = {
        logcpm <- log2(t(t(pb + 0.5) / (colSums(pb) + 1)) * 1e6)
        a <- group == levels(group)[1L]
        s <- t(apply(logcpm, 1L, function(v) {
          if (stats::sd(v[a]) + stats::sd(v[!a]) == 0) return(c(mean(v[!a]) - mean(v[a]), 1))
          c(mean(v[!a]) - mean(v[a]), stats::t.test(v[!a], v[a])$p.value)
        }))
        data.frame(log2_fold_change = s[, 1L], p_value = s[, 2L])
      },
      pseudobulk_edger = {
        if (!requireNamespace("edgeR", quietly = TRUE)) stop("method 'pseudobulk_edger' needs the edgeR package")
        y <- edgeR::calcNormFactors(edgeR::DGEList(pb, group = group))
        y <- edgeR::estimateDisp(y, stats::model.matrix(~group))
        et <- edgeR::exactTest(y)$table
        data.frame(log2_fold_change = et$logFC, p_value = et$PValue)
      },
      pseudobulk_deseq2 = {
        if (!requireNamespace("DESeq2", quietly = TRUE)) stop("method 'pseudobulk_deseq2' needs the DESeq2 package")
        dds <- DESeq2::DESeqDataSetFromMatrix(pb, data.frame(group = group), ~group)
        res <- as.data.frame(DESeq2::results(suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))))
        data.frame(log2_fold_change = res$log2FoldChange, p_value = ifelse(is.na(res$pvalue), 1, res$pvalue))
      })
    result <- cbind(data.frame(gene_id = rownames(pb), stringsAsFactors = FALSE), result)
  }
  result$q_value <- stats::p.adjust(result$p_value, "BH")
  result$significant <- result$q_value <= fdr
  result$method <- method
  result
}

#' Score condition DE calls against the simulated condition truth
#'
#' @param de Output of test_scrnaseq_conditions().
#' @param condition_truth Condition truth table or path (gene_id, effect_log2).
#' @param min_effect Smallest absolute true log2 effect counted as a real
#'   difference (tiny simulated effects are practically undetectable).
#' @return One-row data frame with recall, precision, F1, and false positives.
#' @export
score_scrnaseq_conditions <- function(de, condition_truth, min_effect = 0.25) {
  if (is.character(condition_truth)) condition_truth <- utils::read.delim(condition_truth, stringsAsFactors = FALSE)
  real <- unique(condition_truth$gene_id[abs(condition_truth$effect_log2) >= min_effect])
  tiny <- setdiff(unique(condition_truth$gene_id), real)
  called <- setdiff(de$gene_id[de$significant], tiny)
  tp <- length(intersect(called, real)); fp <- length(setdiff(called, real))
  recall <- if (length(real)) tp / length(real) else NA_real_
  precision <- if (length(called)) tp / length(called) else NA_real_
  data.frame(analysis = "condition DE", method = de$method[1L], true_de_genes = length(real),
             called = length(called), recovered = tp, false_positives = fp, recall = recall, precision = precision,
             f1 = if (is.na(recall) || is.na(precision) || recall + precision == 0) 0 else 2 * recall * precision / (recall + precision),
             note = paste0("genes with |true effect| < ", min_effect, " are excluded from both counts"))
}
