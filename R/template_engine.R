# Count simulation from real reference data (templates in inst/extdata/templates/).
#
# A template is a real raw-count matrix (features x samples) with a sample
# sheet. Its realism (count distribution, zeros, library sizes, feature-to-
# feature structure, biological variation between samples) carries into the
# simulation; the design effects are added with known values:
#
# * resample: each experimental unit gets its own real template sample from one
#   homogeneous group; effects are added by binomial thinning (Gerard 2020,
#   seqgendiff): counts of feature g in unit u are thinned with probability
#   2^(effect[g,u] - max_u effect[g,u]), so log2 fold changes between
#   conditions equal the programmed effects in expectation.
# * parametric: when there are more units than template samples, counts are
#   negative binomial with each feature's mean and dispersion estimated from
#   the template group, library sizes drawn from the template, and the effects
#   applied to the means.

.template_dir <- function() {
  p <- system.file("extdata", "templates", package = "simitall")
  if (!nzchar(p)) p <- file.path("inst", "extdata", "templates")
  p
}

#' Reference templates for count simulation
#'
#' Lists the real count matrices bundled as simulation templates (bulk RNA-seq,
#' small RNA-seq, 16S microbiome, ...), with organism, assays, source and
#' license.
#'
#' @return A data frame, one row per template.
#' @examples
#' simitall_templates()[, c("id", "organism", "assays")]
#' @export
simitall_templates <- function() {
  utils::read.delim(file.path(.template_dir(), "templates.tsv"), stringsAsFactors = FALSE,
                    colClasses = "character", quote = "")
}

#' Load one reference template
#'
#' @param id Template id from [simitall_templates()].
#' @return A list with `counts` (integer matrix, features x samples),
#'   `samples` (data frame with at least `sample` and `group`), `features`
#'   (annotation or `NULL`) and `info` (its row of `templates.tsv`).
#' @export
simitall_template <- function(id) {
  info <- simitall_templates()
  info <- info[info$id == id, , drop = FALSE]
  if (!nrow(info)) stop("Unknown template '", id, "'. See simitall_templates().", call. = FALSE)
  d <- file.path(.template_dir(), id)
  x <- utils::read.delim(file.path(d, "counts.tsv.gz"), check.names = FALSE, stringsAsFactors = FALSE)
  counts <- as.matrix(x[, -1, drop = FALSE]); rownames(counts) <- x[[1]]; storage.mode(counts) <- "integer"
  samples <- utils::read.delim(file.path(d, "samples.tsv"), stringsAsFactors = FALSE, colClasses = "character")
  fp <- file.path(d, "features.tsv.gz")
  features <- if (file.exists(fp)) utils::read.delim(fp, stringsAsFactors = FALSE) else NULL
  list(counts = counts, samples = samples, features = features, info = as.list(info))
}

#' Simulate counts for a design from a real reference template
#'
#' Gives every experimental unit a count profile built from real data (a
#' template from [simitall_templates()] or your own count matrix) and adds
#' treatment and batch effects of known size. Works for any count assay:
#' small RNA-seq, 16S/ITS microbiome, bulk RNA-seq, footprint or IP counts.
#'
#' @param template A template id, or a list with `counts` (features x samples
#'   integer matrix) and `samples` (data frame with `sample` and `group`).
#' @param design Data frame with one row per unit: `sample`, `condition`,
#'   `batch`. The first condition level (factor level order, else order of
#'   appearance) is the reference.
#' @param out_prefix Output prefix; writes `<prefix>.counts.tsv` and
#'   `<prefix>.condition_truth.tsv`.
#' @param condition_effect_fraction Share of features that respond to each
#'   non-reference condition.
#' @param condition_effect_sd SD of the true log2 fold changes.
#' @param batch_effect_sd SD of per-feature log2 batch effects.
#' @param template_group Which template group to draw from (a value of
#'   `samples$group`). Default: the template's `default_group`, else the
#'   largest group. One homogeneous group keeps the differences between units
#'   biological, not tissue or site differences.
#' @param mode `"auto"` (resample when there are at least as many template
#'   samples as units, else parametric), `"resample"` or `"parametric"`.
#' @param min_mean Features with a template mean below this are not chosen to
#'   respond (their fold changes could not be seen).
#' @param seed Random seed.
#' @return Invisibly, a list with `counts`, `condition_truth` (paths) and
#'   `truth` (template, group, mode, template samples used per unit).
#' @examples
#' d <- data.frame(sample = paste0("u", 1:6), condition = rep(c("ctrl", "heat"), each = 3), batch = "b1")
#' r <- simulate_counts_from_template("small_rnaseq_mouse_mirna_atlas", d, file.path(tempdir(), "srna"))
#' head(read.delim(r$condition_truth))
#' @export
simulate_counts_from_template <- function(template, design, out_prefix, condition_effect_fraction = 0.1,
                                          condition_effect_sd = 1, batch_effect_sd = 0.15,
                                          template_group = NULL, mode = c("auto", "resample", "parametric"),
                                          min_mean = 5, seed = 1) {
  mode <- match.arg(mode)
  tp <- if (is.character(template)) simitall_template(template) else template
  tid <- if (is.character(template)) template else (tp$info$id %||% "user_template")
  design <- as.data.frame(design, stringsAsFactors = FALSE)
  if (!all(c("sample", "condition") %in% names(design))) stop("design needs columns sample and condition.", call. = FALSE)
  if (is.null(design$batch)) design$batch <- "batch1"
  lev <- if (is.factor(design$condition)) levels(droplevels(design$condition)) else unique(as.character(design$condition))
  if (!is.null(attr(design, "levels"))) lev <- intersect(attr(design, "levels"), unique(as.character(design$condition)))
  cond <- as.character(design$condition); bat <- as.character(design$batch)
  set.seed(seed)

  # template group: one homogeneous set of real samples
  grp <- template_group %||% tp$info$default_group
  if (is.null(grp) || !nzchar(grp) || !grp %in% tp$samples$group) grp <- names(sort(table(tp$samples$group), decreasing = TRUE))[1]
  cols <- tp$samples$sample[tp$samples$group == grp]
  M <- tp$counts[, cols, drop = FALSE]
  M <- M[rowSums(M) > 0, , drop = FALSE]
  n <- nrow(design); G <- nrow(M)
  if (mode == "auto") mode <- if (ncol(M) >= n) "resample" else "parametric"
  if (mode == "resample" && ncol(M) < n) {
    stop(sprintf("resample needs one template sample per unit: group '%s' has %d samples for %d units (use mode = 'parametric').", grp, ncol(M), n), call. = FALSE)
  }

  # true effects (log2): condition effects on a share of the expressed features, batch effects on all
  lib <- colSums(M)
  q <- rowMeans(sweep(M, 2, lib / stats::median(lib), "/"))
  expressed <- which(q >= min_mean)
  if (!length(expressed)) expressed <- seq_len(G)
  E <- matrix(0, G, length(lev), dimnames = list(rownames(M), lev))
  truth <- data.frame(gene_id = character(), condition = character(), effect_log2 = numeric(), stringsAsFactors = FALSE)
  n_de <- min(length(expressed), round(G * condition_effect_fraction))
  if (n_de > 0 && length(lev) > 1) {
    sel <- sample(expressed, n_de)
    for (l in lev[-1]) {
      e <- stats::rnorm(n_de, 0, condition_effect_sd)
      E[sel, l] <- e
      truth <- rbind(truth, data.frame(gene_id = rownames(M)[sel], condition = l, effect_log2 = e, stringsAsFactors = FALSE))
    }
  }
  blev <- unique(bat)
  B <- matrix(if (length(blev) > 1) stats::rnorm(G * length(blev), 0, batch_effect_sd) else 0, G, length(blev), dimnames = list(NULL, blev))
  Tm <- E[, cond, drop = FALSE] + B[, bat, drop = FALSE]

  if (mode == "resample") {
    pick <- sample(colnames(M), n)
    base <- M[, pick, drop = FALSE]
    p <- 2^(Tm - apply(Tm, 1, max))
    Y <- matrix(stats::rbinom(length(base), as.vector(base), as.vector(p)), G, n)
  } else {
    pick <- rep(NA_character_, n)
    Mn <- sweep(M, 2, lib / stats::median(lib), "/")
    mu <- rowMeans(Mn); v <- apply(Mn, 1, stats::var)
    phi <- pmax((v - mu) / pmax(mu, 1e-8)^2, 1e-3)
    size_factor <- sample(lib, n, replace = TRUE) / stats::median(lib)
    Mu <- outer(mu, size_factor) * 2^Tm
    Y <- matrix(stats::rnbinom(length(Mu), mu = as.vector(Mu), size = rep(1 / phi, n)), G, n)
  }
  dimnames(Y) <- list(rownames(M), design$sample)
  dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)
  cp <- paste0(out_prefix, ".counts.tsv"); tp_ <- paste0(out_prefix, ".condition_truth.tsv")
  utils::write.table(data.frame(gene_id = rownames(Y), Y, check.names = FALSE), cp, sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(truth, tp_, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(list(counts = cp, condition_truth = tp_,
                 truth = list(template = tid, template_group = grp, mode = mode, n_template_samples = ncol(M),
                              n_features = G, reference = lev[1], template_sample_per_unit = stats::setNames(as.list(pick), design$sample),
                              source = tp$info$source %||% "user", license = tp$info$license %||% "user")))
}
