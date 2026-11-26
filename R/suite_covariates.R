# Covariates: measured numeric predictors (Bioana docs/contracts/regression.md; Biodes analysis.covariates).
# Generic terms covariate, covariate2, ... (role covariate), one value per unit in sample_metadata_template.tsv.
#
# Values: normal per unit with mean / sd from analysis.covariates (default 50 / 10), correlated covariate_cor
# between covariates. A covariate with transform log / log10 / log2 is log-normal with that mean and sd (so its
# log is normal); sqrt is normal floored at 0. Bioana fits the transformed value, so slopes are per unit of it.
# Effect: slope x (transformed value - its mean), added to the outcome's Gaussian scale. Slope per unit: the
# regression goal's size (scale per_unit), else the treatment effect per SD of the covariate. A covariate x
# covariate or treatment x covariate term gets a centred product (interaction_slope x the product).
# Planted problems for the diagnostics (suite_covariates): heteroscedastic noise, curvature, outliers.

.suite_covariate_terms <- function(p) names(p$role)[p$role == "covariate"]

.suite_simulate_covariates <- function(spec, p, meta, a, note) {
  terms <- .suite_covariate_terms(p)
  if (!length(terms)) return(list(meta = meta, info = list()))
  n <- nrow(meta); K <- length(terms)
  rho <- as.numeric(a$covariate_cor %||% 0.3)
  R <- matrix(rho, K, K); diag(R) <- 1
  if (K > 1 && min(eigen(R, only.values = TRUE)$values) <= 0) stop("covariate_cor gives a correlation matrix that is not positive definite.", call. = FALSE)
  Z <- matrix(stats::rnorm(n * K), n, K) %*% chol(R)
  covs <- spec$analysis$covariates %||% list()
  info <- list()
  for (k in seq_len(K)) {
    v <- terms[k]; col <- p$col[[v]]
    ce <- Filter(function(z) identical(z$term, v) || identical(z$column, col) || identical(z$name, col), covs)
    ce <- if (length(ce)) ce[[1]] else list()
    m <- as.numeric(ce$mean %||% a$covariate_mean %||% 50); s <- as.numeric(ce$sd %||% a$covariate_sd %||% 10)
    tf <- ce$transform %||% "none"
    raw <- if (tf %in% c("log", "log10", "log2")) {
      sl <- sqrt(log(1 + (s / m)^2)); exp(log(m) - sl^2 / 2 + sl * Z[, k])
    } else if (tf == "sqrt") pmax(m + s * Z[, k], 0) else m + s * Z[, k]
    meta[[col]] <- round(raw, 4)
    info[[v]] <- list(column = col, name = ce$name %||% col, unit = ce$unit, transform = tf, mean = m, sd = s)
  }
  if (K > 1) note(sprintf("covariates %s correlated %s (covariate_cor)", paste(terms, collapse = ", "), format(rho)))
  list(meta = meta, info = info, cor = if (K > 1) rho)
}

# Adds the covariate effects to y; returns y, the truth (slopes per unit of the fitted value) and planted problems.
.suite_covariate_effects <- function(y, frame, p, spec, a, eff, rsd, info, note, outcome_type = "continuous", logit_coef = NULL) {
  terms <- .suite_covariate_terms(p)
  if (!length(terms)) return(list(y = y, truth = NULL))
  ge <- spec$analysis$goal_effect$regression %||% list()
  gterms <- unlist(spec$analysis$goal_terms$regression) %||% terms
  x <- lapply(stats::setNames(terms, terms), function(v) frame[[v]])     # transformed values (as Bioana fits them)
  xc <- lapply(x, function(z) z - mean(z, na.rm = TRUE))
  binary <- identical(outcome_type, "binary") && !is.null(logit_coef)
  # slope on the Gaussian scale for a stated size: per unit (continuous) or log OR per unit (binary)
  stated <- function(v) {
    if (is.null(ge$size) || !v %in% gterms) return(NULL)
    sz <- as.numeric(ge$size)
    if (binary) { lor <- if (identical(ge$scale, "odds_ratio_per_unit")) log(sz) else sz; return(lor * rsd / logit_coef) }
    sz
  }
  # Biodes' a-priori f2 (pwr.f2.test) with no stated slope: equal standardized slopes explaining f2 x residual variance
  f2 <- .suite_apriori_effect(spec)$f2
  b_sd <- if (!is.null(f2) && is.null(ge$size)) { K <- length(terms); rho <- as.numeric(a$covariate_cor %||% 0.3)
    tot <- K + K * (K - 1) * (if (K > 1) rho else 0); rsd * sqrt(f2 / tot) } else eff
  if (!is.null(f2) && is.null(ge$size)) note(sprintf("covariate slopes from Biodes' a-priori f2 = %s (equal standardized slopes)", format(f2)))
  slope <- vapply(terms, function(v) { st <- stated(v); if (!is.null(st)) st else b_sd / stats::sd(x[[v]], na.rm = TRUE) }, 0)
  if (is.null(ge$size) && is.null(f2)) note(sprintf("covariate slopes: the treatment effect per SD of each covariate (%s)", paste(sprintf("%s %s", terms, signif(slope, 3)), collapse = ", ")))
  for (v in terms) y <- y + slope[[v]] * xc[[v]]
  # curvature: a quadratic term in the first covariate (its effect over +-1 SD bends by `curvature` x the slope)
  curv <- as.numeric(a$curvature %||% 0)
  if (curv != 0) { v <- terms[1]; sdx <- stats::sd(x[[v]], na.rm = TRUE); y <- y + curv * slope[[v]] * (xc[[v]]^2 - sdx^2) / sdx
    note(sprintf("planted curvature in %s (%s)", v, format(curv))) }
  # interactions with a covariate: centred products
  int_terms <- Filter(function(l) any(all.vars(str2lang(l)) %in% terms), p$fixed_labels[grepl(":", p$fixed_labels)])
  int_truth <- list()
  bi <- as.numeric(a$interaction_slope_ratio %||% 0.5)
  for (l in int_terms) {
    vs <- all.vars(str2lang(l))
    if (all(vs %in% terms)) {
      b <- bi * slope[[vs[1]]] / stats::sd(x[[vs[2]]], na.rm = TRUE)
      y <- y + b * xc[[vs[1]]] * xc[[vs[2]]]
      int_truth[[l]] <- list(coefficient = b, note = "per unit of each covariate (product term)")
    } else {
      cv <- intersect(vs, terms)[1]; tv <- setdiff(vs, terms)[1]
      lev <- levels(frame[[tv]]); sc <- (seq_along(lev) - 1) / max(1, length(lev) - 1)
      b <- bi * slope[[cv]]
      y <- y + b * sc[as.integer(frame[[tv]])] * xc[[cv]]
      int_truth[[l]] <- list(slope_change = stats::setNames(as.list(b * sc), lev), note = sprintf("change of the %s slope per level of %s vs the reference", cv, tv))
    }
  }
  # true coefficients of the model as Bioana fits it (uncentred covariates): main slopes absorb -b x mean(other)
  true <- as.list(slope)
  for (l in names(int_truth)) {
    vs <- all.vars(str2lang(l))
    if (all(vs %in% terms)) {
      b <- int_truth[[l]]$coefficient
      true[[vs[1]]] <- true[[vs[1]]] - b * mean(x[[vs[2]]], na.rm = TRUE)
      true[[vs[2]]] <- true[[vs[2]]] - b * mean(x[[vs[1]]], na.rm = TRUE)
      true[[l]] <- b
    }
  }
  out <- list(true_slopes = if (!binary) true, true_log_odds_ratios = if (binary) lapply(true, function(b) b * logit_coef / rsd),
              slope_scale = if (binary) "log odds per unit of the fitted (transformed) covariate" else "outcome per unit of the fitted (transformed) covariate",
              interactions = if (length(int_truth)) int_truth, covariates = info, curvature = if (curv != 0) curv)
  list(y = y, truth = out[!vapply(out, is.null, logical(1))], slope = slope, x = x)
}

# Planted problems on the residuals: heteroscedastic noise (SD grows with the first covariate), outliers.
.suite_covariate_noise <- function(eps, x, a, rsd, note) {
  het <- as.numeric(a$heteroscedasticity %||% 0)
  out <- list()
  if (het != 0 && length(x)) {
    z <- x[[1]]; zs <- (z - mean(z, na.rm = TRUE)) / stats::sd(z, na.rm = TRUE)
    eps <- eps * exp(het * zs) / sqrt(mean(exp(2 * het * zs), na.rm = TRUE))
    out$heteroscedasticity <- list(strength = het, over = names(x)[1], note = "residual SD x exp(strength x standardized covariate)")
    note(sprintf("planted heteroscedastic noise over %s (%s)", names(x)[1], format(het)))
  }
  k <- as.integer(a$outliers %||% 0)
  if (k > 0) {
    rows <- sort(sample(which(!is.na(eps)), min(k, sum(!is.na(eps)))))
    eps[rows] <- eps[rows] + sample(c(-1, 1), length(rows), replace = TRUE) * as.numeric(a$outlier_size %||% 5) * rsd
    out$planted_outliers <- as.list(rows)
    note(sprintf("planted %d outliers (%s residual SDs)", length(rows), format(as.numeric(a$outlier_size %||% 5))))
  }
  list(eps = eps, truth = out)
}
