# DNA methylation simulation with planted differentially methylated regions (DMRs).
#
# Three engines share one truth format:
# * template: real bisulfite reads (methylation_counts templates, e.g. ENCODE
#   mouse WGBS). With enough template samples each unit is a real sample
#   ("resample"); otherwise each site's level, the between-sample variation
#   and the coverage are learned from the template ("parametric").
# * native: synthetic cytosines in CG / CHG / CHH contexts, grouped in regions
#   of one chromatin state (presets in inst/suite/methylation_presets.json).
# * array: beta values per probe learned from a beta template (450k/EPIC).
# A DMR is a run of consecutive sites whose methylation moves by `dmr_delta`
# in each non-reference condition: up when the region is mostly unmethylated,
# down otherwise.

.logit <- function(p) log(p / (1 - p))
.ilogit <- function(x) 1 / (1 + exp(-x))
.clamp <- function(p, lo = 0.005, hi = 0.995) pmin(pmax(p, lo), hi)

# Non-overlapping runs of `len` consecutive sites on one chromosome; returns site indices per DMR.
.meth_pick_runs <- function(chr, n_runs, len) {
  starts <- which(c(TRUE, chr[-1] == chr[-length(chr)]))      # any site can start a run if the run stays on its chromosome
  ok <- starts[starts + len - 1 <= length(chr)]
  ok <- ok[chr[ok] == chr[pmin(ok + len - 1, length(chr))]]
  runs <- list(); used <- logical(length(chr))
  for (s in sample(ok)) {
    if (length(runs) >= n_runs) break
    idx <- s:(s + len - 1)
    if (any(used[max(1, s - 1):min(length(chr), s + len)])) next
    used[idx] <- TRUE; runs[[length(runs) + 1]] <- idx
  }
  runs
}

#' Simulate DNA methylation for a design, with known DMRs
#'
#' Bisulfite sequencing (WGBS, RRBS, EM-seq) gives, per cytosine and unit, the
#' reads covering it and the methylated reads; methylation arrays give a beta
#' value per probe. A share of regions changes methylation between conditions
#' by a known amount, and the truth is written next to the data.
#'
#' @param design Data frame with one row per unit: `sample`, `condition`,
#'   `batch`. The first condition level is the reference.
#' @param out_prefix Output prefix.
#' @param engine `"auto"` (template of a `methylation_counts` kind when one is
#'   given, else native), `"template"`, `"native"` or `"array"`.
#' @param template Template id ([simitall_templates()]) for the template or
#'   array engines.
#' @param template_group Template group to learn from (default: its
#'   `default_group`).
#' @param kingdom `"plant"` or `"animal"`: the native preset (contexts and
#'   chromatin states).
#' @param contexts Contexts for the native engine (default: the preset's,
#'   CG/CHG/CHH for plants, CG for animals).
#' @param n_sites Cytosines for the native engine.
#' @param block_sites Sites per chromatin-state region (native engine).
#' @param mean_coverage Mean reads per cytosine. Default: the template's own
#'   depth, else 15.
#' @param n_dmrs Differentially methylated regions (array: positions).
#' @param dmr_sites Cytosines per DMR.
#' @param dmr_delta Methylation difference in a DMR (proportion, 0-1).
#' @param dmr_contexts Contexts a DMR changes (default: all simulated).
#' @param biological_sd Between-unit SD of regional methylation on the logit
#'   scale. Default: estimated from the template, else 0.3.
#' @param batch_sd Batch SD on the logit scale.
#' @param mode Template engine: `"auto"`, `"resample"` or `"parametric"`.
#' @param seed Random seed.
#' @return Invisibly, a list of paths (`sites`, `methylated`, `coverage` or
#'   `beta`, `dmr_truth`, `site_truth`) and `truth` (engine, template, DMR
#'   count, parameters used).
#' @examples
#' d <- data.frame(sample = paste0("p", 1:6), condition = rep(c("wt", "mut"), each = 3), batch = "b1")
#' r <- simulate_methylation(d, file.path(tempdir(), "meth"), engine = "native", kingdom = "plant", n_sites = 3000, n_dmrs = 10)
#' head(read.delim(r$dmr_truth))
#' @export
simulate_methylation <- function(design, out_prefix, engine = c("auto", "template", "native", "array"),
                                 template = NULL, template_group = NULL, kingdom = c("animal", "plant"),
                                 contexts = NULL, n_sites = 20000L, block_sites = 40L, mean_coverage = NULL,
                                 n_dmrs = 50L, dmr_sites = 20L, dmr_delta = 0.25, dmr_contexts = NULL,
                                 biological_sd = NULL, batch_sd = 0.1, mode = c("auto", "resample", "parametric"),
                                 seed = 1) {
  engine <- match.arg(engine); kingdom <- match.arg(kingdom); mode <- match.arg(mode)
  design <- as.data.frame(design, stringsAsFactors = FALSE)
  if (is.null(design$batch)) design$batch <- "batch1"
  lev <- if (!is.null(attr(design, "levels"))) attr(design, "levels") else if (is.factor(design$condition)) levels(droplevels(design$condition)) else unique(as.character(design$condition))
  cond <- as.character(design$condition); n <- nrow(design)
  tp <- if (!is.null(template)) simitall_template(template) else NULL
  if (engine == "auto") engine <- if (!is.null(tp) && identical(tp$info$kind, "beta")) "array" else if (!is.null(tp)) "template" else "native"
  if (engine %in% c("template", "array") && is.null(tp)) stop("engine '", engine, "' needs a template.", call. = FALSE)
  set.seed(seed)
  presets <- .suite_rules("methylation_presets")
  dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)
  grp <- NULL
  if (!is.null(tp)) {
    grp <- template_group %||% tp$info$default_group
    if (!grp %in% tp$samples$group) grp <- names(sort(table(tp$samples$group), decreasing = TRUE))[1]
  }

  # ---- array: beta per probe -------------------------------------------------------
  if (engine == "array") {
    B <- tp$beta[, tp$samples$sample[tp$samples$group == grp], drop = FALSE]
    mu <- .clamp(rowMeans(B), 0.01, 0.99); v <- pmax(apply(B, 1, stats::var), 1e-5)
    phi <- pmin(pmax(mu * (1 - mu) / v - 1, 2), 500)
    G <- nrow(B); mu_u <- matrix(mu, G, n)
    dmp <- sample(G, min(G, n_dmrs))
    truth <- data.frame(dmr_id = character(), feature_id = character(), condition = character(), delta = numeric(), base = numeric())
    for (l in lev[-1]) {
      s <- ifelse(mu[dmp] < 0.5, 1, -1)
      mu_u[dmp, cond == l] <- .clamp(mu[dmp] + s * dmr_delta, 0.01, 0.99)
      truth <- rbind(truth, data.frame(dmr_id = paste0("dmp_", seq_along(dmp)), feature_id = rownames(B)[dmp], condition = l,
                                       delta = .clamp(mu[dmp] + s * dmr_delta, 0.01, 0.99) - mu[dmp], base = mu[dmp]))
    }
    beta <- matrix(stats::rbeta(G * n, as.vector(mu_u * phi), as.vector((1 - mu_u) * phi)), G, n, dimnames = list(rownames(B), design$sample))
    bp <- paste0(out_prefix, ".beta.tsv.gz"); z <- gzfile(bp, "w")
    utils::write.table(data.frame(feature_id = rownames(beta), round(beta, 4), check.names = FALSE), z, sep = "\t", row.names = FALSE, quote = FALSE); close(z)
    tpth <- paste0(out_prefix, ".dmr_truth.tsv"); utils::write.table(truth, tpth, sep = "\t", row.names = FALSE, quote = FALSE)
    return(invisible(list(beta = bp, dmr_truth = tpth,
                          truth = list(engine = "array", template = tp$info$id, template_group = grp, n_probes = G, n_dmps = length(dmp),
                                       dmr_delta = dmr_delta, reference = lev[1], unit = "beta value"))))
  }

  # ---- sites and their baseline level ----------------------------------------------
  if (engine == "native") {
    ps <- presets[[kingdom]]
    ctx <- contexts %||% names(ps$context_share)
    share <- unlist(ps$context_share[ctx]); share <- share / sum(share)
    nb <- ceiling(n_sites / block_sites)
    state <- sample(names(ps$state_share), nb, replace = TRUE, prob = unlist(ps$state_share))
    block <- rep(seq_len(nb), each = block_sites)[seq_len(n_sites)]
    context <- sample(ctx, n_sites, replace = TRUE, prob = share)
    gap <- stats::rgeom(n_sites, 1 / unlist(presets$site_spacing_bp)[context]) + 2L
    sites <- data.frame(chr = "chr1", pos = cumsum(gap), strand = sample(c("+", "-"), n_sites, TRUE), context = context, region = block,
                        state = state[block], stringsAsFactors = FALSE)
    lvl <- mapply(function(s, c) ps$level[[s]][[c]], sites$state, sites$context)
    P0 <- .ilogit(.logit(.clamp(lvl)) + stats::rnorm(n_sites, 0, presets$site_sd_logit))
    bio <- biological_sd %||% 0.3
    cov_mu <- rep(mean_coverage %||% 15, n_sites) * stats::rgamma(n_sites, 4, 4)
  } else {
    cols <- tp$samples$sample[tp$samples$group == grp]
    N <- tp$counts[, cols, drop = FALSE]; M <- tp$methylated[, cols, drop = FALSE]
    f <- tp$features
    sites <- data.frame(chr = f$chr, pos = f$pos, strand = f$strand, context = f$context, stringsAsFactors = FALSE)
    sites$region <- paste(sites$chr, sites$pos %/% 1000)
    P0 <- .clamp((rowSums(M) + 0.5) / (rowSums(N) + 1))
    # between-sample variation: SD of logit levels at well-covered sites, minus the binomial part
    w <- rowSums(N >= 20) == ncol(N)
    bio <- biological_sd %||% if (sum(w) > 100 && ncol(N) > 1) {
      L <- .logit(.clamp((M[w, ] + 0.5) / (N[w, ] + 1)))
      tot <- stats::median(apply(L, 1, stats::var)); samp <- stats::median(rowMeans(1 / ((N[w, ] + 1) * .clamp(P0[w]) * (1 - .clamp(P0[w])))))
      sqrt(max(tot - samp, 0.01))
    } else 0.3
    cov_mu <- rowMeans(N)
    if (!is.null(mean_coverage)) cov_mu <- cov_mu * mean_coverage / mean(cov_mu)
    if (mode == "auto") mode <- if (ncol(N) >= n) "resample" else "parametric"
    if (mode == "resample" && ncol(N) < n) stop(sprintf("resample needs one template sample per unit: group '%s' has %d for %d units.", grp, ncol(N), n), call. = FALSE)
  }
  G <- nrow(sites)
  sites$feature_id <- paste0(sites$chr, ":", sites$pos, ":", sites$strand)

  # ---- DMRs ------------------------------------------------------------------------
  dctx <- dmr_contexts %||% unique(sites$context)
  runs <- .meth_pick_runs(sites$chr, n_dmrs, dmr_sites)
  dmr_of <- rep(NA_character_, G)
  Plev <- matrix(P0, G, length(lev), dimnames = list(NULL, lev))
  dtruth <- list()
  for (k in seq_along(runs)) {
    idx <- runs[[k]]; hit <- idx[sites$context[idx] %in% dctx]
    if (!length(hit)) next
    id <- sprintf("dmr_%03d", k); dmr_of[idx] <- id
    s <- if (mean(P0[hit]) < 0.5) 1 else -1
    for (l in lev[-1]) {
      Plev[hit, l] <- .clamp(P0[hit] + s * dmr_delta)
      dtruth[[length(dtruth) + 1]] <- data.frame(dmr_id = id, chr = sites$chr[idx[1]], start = min(sites$pos[idx]), end = max(sites$pos[idx]),
                                                 n_sites = length(hit), contexts = paste(sort(unique(sites$context[hit])), collapse = ","),
                                                 condition = l, delta = mean(Plev[hit, l] - P0[hit]), base = mean(P0[hit]))
    }
  }
  dtruth <- if (length(dtruth)) do.call(rbind, dtruth) else data.frame()

  # ---- per-unit methylated and covering reads --------------------------------------
  ci <- match(cond, lev)
  if (engine == "template" && mode == "resample") {
    pick <- sample(colnames(N), n)
    Nu <- N[, pick, drop = FALSE]; Mu <- M[, pick, drop = FALSE]
    if (!is.null(mean_coverage)) {                                   # thin both read types alike: level unchanged
      keep <- min(1, mean_coverage / mean(N)); um <- Nu - Mu
      Mu <- matrix(stats::rbinom(length(Mu), Mu, keep), G); um <- matrix(stats::rbinom(length(um), um, keep), G); Nu <- Mu + um
    }
    for (u in which(ci > 1)) {                                       # DMR sites: methylated reads redrawn at the shifted level
      d <- which(!is.na(dmr_of) & Plev[, ci[u]] != P0)
      p_hat <- (Mu[d, u] + 0.5) / (Nu[d, u] + 1)
      Mu[d, u] <- stats::rbinom(length(d), Nu[d, u], .clamp(p_hat + (Plev[d, ci[u]] - P0[d])))
    }
  } else {
    region <- as.integer(factor(sites$region))
    R <- matrix(stats::rnorm(max(region) * n, 0, bio), max(region), n)[region, , drop = FALSE]
    bl <- unique(design$batch); Bt <- matrix(stats::rnorm(max(region) * length(bl), 0, batch_sd), max(region), length(bl))[region, , drop = FALSE]
    Pu <- .ilogit(.logit(Plev[, ci, drop = FALSE]) + R + Bt[, match(design$batch, bl), drop = FALSE])
    Nu <- matrix(stats::rnbinom(G * n, mu = rep(cov_mu, n), size = presets$coverage_size), G, n)
    Mu <- matrix(stats::rbinom(G * n, Nu, Pu), G, n)
  }
  dimnames(Nu) <- dimnames(Mu) <- list(sites$feature_id, design$sample)
  storage.mode(Nu) <- storage.mode(Mu) <- "integer"

  wz <- function(x, f) { p <- paste0(out_prefix, f); z <- gzfile(p, "w"); utils::write.table(x, z, sep = "\t", row.names = FALSE, quote = FALSE); close(z); p }
  paths <- list(
    sites = wz(sites[, c("feature_id", "chr", "pos", "strand", "context")], ".sites.tsv.gz"),
    methylated = wz(data.frame(feature_id = sites$feature_id, Mu, check.names = FALSE), ".methylated.tsv.gz"),
    coverage = wz(data.frame(feature_id = sites$feature_id, Nu, check.names = FALSE), ".coverage.tsv.gz"),
    site_truth = wz(data.frame(feature_id = sites$feature_id, dmr_id = dmr_of, stats::setNames(as.data.frame(round(Plev, 4)), paste0("p_", lev)), check.names = FALSE), ".site_truth.tsv.gz"))
  paths$dmr_truth <- paste0(out_prefix, ".dmr_truth.tsv")
  utils::write.table(dtruth, paths$dmr_truth, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(c(paths, list(truth = list(
    engine = engine, mode = if (engine == "template") mode else NULL, template = tp$info$id, template_group = grp, kingdom = if (engine == "native") kingdom,
    contexts = unique(sites$context), n_sites = G, n_dmrs = length(unique(stats::na.omit(dmr_of))), dmr_sites = dmr_sites, dmr_delta = dmr_delta,
    dmr_contexts = dctx, biological_sd = round(bio, 3), mean_coverage = round(mean(Nu), 1), reference = lev[1]))))
}
