# m6A / MeRIP-seq simulation with known peaks and differentially methylated peaks.
#
# MeRIP-seq sequences two libraries per unit: the input (fragmented RNA, i.e.
# expression) and the m6A immunoprecipitate (IP). Counts are per window: m6A
# peaks and background windows on the same genes. For window w of gene g in
# unit j:
#   input  ~ NB(mean = depth_j * expr_gj)
#   IP     ~ NB(mean = depth_j * eff_j * expr_gj * (1 + (enrichment - 1) * m_wj))
# where m_wj is the share of transcripts methylated at that site (0 for
# background windows). Expression changes move IP and input together; m6A
# changes move IP only, so a method has to use the input to tell them apart.

#' Simulate m6A MeRIP-seq for a design, with known peaks and changes
#'
#' Writes IP and input counts per window (m6A peaks plus background windows on
#' the same genes) for each unit, with the truth: which windows are peaks, the
#' methylation level of each peak, which peaks change methylation between
#' conditions and by how much, and which genes change expression (these move IP
#' and input together).
#'
#' @param design Data frame with one row per unit: `sample`, `condition`,
#'   `batch`. The first condition level is the reference.
#' @param out_prefix Output prefix.
#' @param template Optional RNA-seq count template ([simitall_templates()],
#'   kind `counts`) whose gene means set the expression levels; default: a
#'   log-normal expression distribution.
#' @param template_group Template group (default: its `default_group`).
#' @param n_genes Genes carrying windows.
#' @param n_peaks m6A peaks (placed on random genes; a gene can have several).
#' @param n_background Background windows (no m6A) on the same genes.
#' @param window_bp Window width.
#' @param mean_input_count Mean input reads per window.
#' @param enrichment IP enrichment of a fully methylated site over background
#'   (IP/input ratio at m = 1).
#' @param m6a_shape1,m6a_shape2 Beta distribution of peak methylation levels.
#' @param dm_fraction Share of peaks whose methylation changes in each
#'   non-reference condition.
#' @param dm_delta Change in the methylated share at a changed peak (0-1).
#' @param expression_effect_fraction,expression_effect_sd Share of genes whose
#'   expression changes, and the SD of their log2 fold changes.
#' @param biological_sd Between-unit SD of peak methylation (logit scale).
#' @param expression_sd Between-unit SD of gene expression (log2 scale).
#' @param ip_efficiency_sd SD of the IP efficiency between units (log2 scale).
#' @param dispersion Negative binomial dispersion of both libraries.
#' @param batch_sd Batch SD of expression (log2 scale).
#' @param seed Random seed.
#' @return Invisibly, a list of paths (`windows`, `ip`, `input`, `peak_truth`,
#'   `gene_truth`) and `truth` (parameters and counts of planted changes).
#' @examples
#' d <- data.frame(sample = paste0("s", 1:6), condition = rep(c("ctrl", "ko"), each = 3), batch = "b1")
#' r <- simulate_merip(d, file.path(tempdir(), "m6a"), n_genes = 300, n_peaks = 200, n_background = 200)
#' head(read.delim(r$peak_truth))
#' @export
simulate_merip <- function(design, out_prefix, template = NULL, template_group = NULL,
                           n_genes = 2000L, n_peaks = 1500L, n_background = 1500L, window_bp = 200L,
                           mean_input_count = 50, enrichment = 20, m6a_shape1 = 2, m6a_shape2 = 3,
                           dm_fraction = 0.1, dm_delta = 0.25,
                           expression_effect_fraction = 0.1, expression_effect_sd = 1,
                           biological_sd = 0.3, expression_sd = 0.2, ip_efficiency_sd = 0.3,
                           dispersion = 0.05, batch_sd = 0.1, seed = 1) {
  design <- as.data.frame(design, stringsAsFactors = FALSE)
  if (is.null(design$batch)) design$batch <- "batch1"
  lev <- if (!is.null(attr(design, "levels"))) attr(design, "levels") else if (is.factor(design$condition)) levels(droplevels(design$condition)) else unique(as.character(design$condition))
  cond <- as.character(design$condition); n <- nrow(design)
  if (!dm_delta > 0 || dm_delta >= 1) stop("dm_delta must be between 0 and 1.", call. = FALSE)
  set.seed(seed)
  dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

  # ---- genes and their expression -----------------------------------------------------
  grp <- NULL
  if (!is.null(template)) {
    tp <- simitall_template(template)
    grp <- template_group %||% tp$info$default_group
    if (!grp %in% tp$samples$group) grp <- names(sort(table(tp$samples$group), decreasing = TRUE))[1]
    cols <- tp$samples$sample[tp$samples$group == grp]
    gm <- rowMeans(sweep(tp$counts[, cols, drop = FALSE], 2, colSums(tp$counts[, cols, drop = FALSE]), "/") * 1e6)
    gm <- gm[gm >= 1]                                         # expressed genes only: m6A is read on expressed transcripts
    pick <- sample(length(gm), min(n_genes, length(gm)))
    gene_id <- names(gm)[pick]; base_log2 <- log2(gm[pick])
  } else {
    gene_id <- sprintf("gene_%05d", seq_len(n_genes))
    base_log2 <- stats::rnorm(n_genes, 4, 1.5)                # log2 CPM-like levels
  }
  G <- length(gene_id)
  base_log2 <- base_log2 - mean(base_log2)                     # relative: mean_input_count sets the depth

  # ---- windows: peaks and background on the same genes ----------------------------------
  wg <- c(sample(G, n_peaks, replace = TRUE), sample(G, n_background, replace = TRUE))
  type <- rep(c("peak", "background"), c(n_peaks, n_background))
  o <- order(wg, type); wg <- wg[o]; type <- type[o]
  W <- length(wg)
  within <- stats::ave(seq_len(W), wg, FUN = seq_along)
  start <- (wg - 1) * 20000L + (within - 1) * 1000L + 1L
  win_id <- sprintf("w%06d", seq_len(W))
  windows <- data.frame(window_id = win_id, gene_id = gene_id[wg], chr = "chr1", start = start, end = start + window_bp - 1L,
                        type = type, stringsAsFactors = FALSE)
  is_peak <- type == "peak"
  m0 <- numeric(W); m0[is_peak] <- .clamp(stats::rbeta(sum(is_peak), m6a_shape1, m6a_shape2), 0.02, 0.98)

  # ---- planted changes ---------------------------------------------------------------------
  peak_idx <- which(is_peak)
  dm <- if (length(lev) > 1) sample(peak_idx, round(dm_fraction * length(peak_idx))) else integer()
  de <- if (length(lev) > 1) sample(G, round(expression_effect_fraction * G)) else integer()
  M_mean <- matrix(m0, W, n)                                   # expected methylated share per window and unit
  E_lfc <- matrix(0, G, n)                                     # expression log2 fold change per gene and unit
  peak_truth <- data.frame(window_id = win_id, gene_id = gene_id[wg], type = type, base_m6a = round(m0, 4),
                           log2_enrichment = round(log2(1 + (enrichment - 1) * m0), 4), dm = seq_len(W) %in% dm, stringsAsFactors = FALSE)
  gene_truth <- data.frame(gene_id = character(), condition = character(), log2fc = numeric(), stringsAsFactors = FALSE)
  for (l in lev[-1]) {
    j <- cond == l
    s <- ifelse(m0[dm] < 0.5, 1, -1)
    m1 <- .clamp(m0[dm] + s * dm_delta, 0.01, 0.99)
    M_mean[dm, j] <- m1
    dcol <- numeric(W); dcol[dm] <- m1 - m0[dm]
    peak_truth[[paste0("delta_", l)]] <- round(dcol, 4)
    # what a method sees: the change in log2(IP / input) at that window
    r <- numeric(W); r[dm] <- log2((1 + (enrichment - 1) * m1) / (1 + (enrichment - 1) * m0[dm]))
    peak_truth[[paste0("log2_ratio_change_", l)]] <- round(r, 4)
    lfc <- stats::rnorm(length(de), 0, expression_effect_sd)
    E_lfc[de, j] <- lfc
    gene_truth <- rbind(gene_truth, data.frame(gene_id = gene_id[de], condition = l, log2fc = round(lfc, 4), stringsAsFactors = FALSE))
  }

  # ---- unit-level variation and counts --------------------------------------------------------
  b <- unique(design$batch)
  batch_eff <- matrix(stats::rnorm(G * length(b), 0, batch_sd), G, length(b), dimnames = list(NULL, b))
  if (length(b)) batch_eff[, 1] <- 0
  expr_log2 <- base_log2 + E_lfc + batch_eff[, design$batch, drop = FALSE] + matrix(stats::rnorm(G * n, 0, expression_sd), G, n)
  M <- matrix(0, W, n)
  M[is_peak, ] <- .ilogit(.logit(.clamp(M_mean[is_peak, , drop = FALSE])) + matrix(stats::rnorm(sum(is_peak) * n, 0, biological_sd), sum(is_peak), n))
  depth <- exp(stats::rnorm(n, 0, 0.2))                         # library size per unit
  eff <- 2^stats::rnorm(n, 0, ip_efficiency_sd)                 # IP efficiency per unit
  mu_in <- mean_input_count * sweep(2^expr_log2[wg, , drop = FALSE], 2, depth, "*")
  mu_ip <- sweep(mu_in * (1 + (enrichment - 1) * M), 2, eff, "*") / (1 + (enrichment - 1) * mean(m0[is_peak]) * n_peaks / W)
  nb <- function(mu) matrix(stats::rnbinom(length(mu), mu = as.vector(mu), size = 1 / dispersion), nrow(mu), ncol(mu),
                            dimnames = list(win_id, design$sample))
  ip <- nb(mu_ip); input <- nb(mu_in)

  wr <- function(m, suffix) {
    p <- paste0(out_prefix, suffix); z <- gzfile(p, "w")
    utils::write.table(data.frame(window_id = rownames(m), m, check.names = FALSE), z, sep = "\t", row.names = FALSE, quote = FALSE)
    close(z); p
  }
  wp <- paste0(out_prefix, ".windows.tsv"); utils::write.table(windows, wp, sep = "\t", row.names = FALSE, quote = FALSE)
  ipp <- wr(ip, ".ip.tsv.gz"); inp <- wr(input, ".input.tsv.gz")
  pt <- paste0(out_prefix, ".peak_truth.tsv"); utils::write.table(peak_truth, pt, sep = "\t", row.names = FALSE, quote = FALSE)
  gt <- paste0(out_prefix, ".gene_truth.tsv"); utils::write.table(gene_truth, gt, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(list(windows = wp, ip = ipp, input = inp, peak_truth = pt, gene_truth = gt,
                 truth = list(engine = if (is.null(template)) "native" else "template_expression", template = template, template_group = grp,
                              n_genes = G, n_peaks = n_peaks, n_background = n_background, n_dm_peaks = length(dm),
                              n_de_genes = length(de), dm_delta = dm_delta, enrichment = enrichment,
                              mean_input_count = mean_input_count, biological_sd = biological_sd, reference = lev[1],
                              ip_efficiency = stats::setNames(as.list(round(eff, 3)), design$sample))))
}
