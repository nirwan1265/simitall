# Categorical outcomes, traits for correlation, and raters for agreement
# (Bioana docs/contracts/categorical_correlation_agreement.md).
#
# Every model starts from the design's Gaussian scale: y_design (intercept + treatment effects) and y_mean
# (plus blocks, random terms, unit effects). `lp` = (y - mu) / rsd is that scale in residual SDs; the full treatment
# effect `eff` is eff / rsd of it.

# expected probabilities per group from a function of the group's lp
.suite_group_truth <- function(groups, y_design, mu, rsd, f) {
  g <- factor(groups)
  lp <- tapply((y_design - mu) / rsd, g, mean)
  lapply(stats::setNames(levels(g), levels(g)), function(k) f(lp[[k]]))
}

# Binary outcome on the logit scale: logit P(y = 1) = qlogis(prevalence) + coef * lp, with lp including unit and
# random effects (so paired designs carry a unit random effect on the logit scale). coef makes the full treatment
# effect an odds ratio `odds_ratio`; the ORs are conditional on the unit.
.suite_binary_logit <- function(a, y_mean, y_design, mu, rsd, eff, groups, cond_effects, random_sd, unit_sd = 0, block_effects = list()) {
  or <- as.numeric(a$odds_ratio %||% 3)
  coef <- if (eff != 0) log(or) / (eff / rsd) else 0
  a0 <- stats::qlogis(as.numeric(a$prevalence))
  lp <- (y_mean - mu) / rsd
  ok <- !is.na(y_mean)
  p <- stats::plogis(a0 + coef * lp)
  y <- ifelse(ok, stats::rbinom(length(p), 1, ifelse(ok, p, 0.5)), NA_integer_)
  # unit effects on the logit scale: random terms, the within-unit effect, and fixed blocks (each block one unit in a
  # paired design), as the SD of their logit-scale values
  re <- c(lapply(random_sd, function(s) round(abs(coef) * s / rsd, 4)), if (unit_sd > 0) list(unit = round(abs(coef) * unit_sd / rsd, 4)),
          lapply(block_effects, function(e) round(abs(coef) * stats::sd(unlist(e)) / rsd, 4)))
  list(type = "binary", y = as.integer(y), model = list(
    family = "logit", prevalence_reference = as.numeric(a$prevalence), intercept_logit = a0, odds_ratio = or,
    log_odds_per_outcome_sd = coef,
    conditional_odds_ratios = lapply(cond_effects, function(e) { e <- unlist(e); as.list(round(exp(coef * (e - e[1]) / rsd), 4)) }),
    random_effect_sd_logit = if (length(re)) re,
    group_probability = .suite_group_truth(groups[ok], y_design[ok], mu, rsd, function(l) round(stats::plogis(a0 + coef * l), 4)),
    note = "P(y = 1) = plogis(intercept + coef x (unit mean - intercept) / residual SD); odds ratios are conditional on the unit (and blocks); group probabilities at an average unit"))
}

# Nominal outcome (multinomial logit). Categories: primary_outcome.categories, else the labels of
# analysis.expected_ratio, else nominal_categories classes. Reference-group probabilities: expected_ratio (a
# segregation ratio, optionally distorted by segregation_distortion, the log-odds shift of the first class) or
# nominal_probabilities. Treatment: category j gets log-odds coef_j x lp vs the reference category, coef rising
# to log(nominal_odds_ratio) per full effect for the last category.
.suite_nominal <- function(a, spec, y_mean, y_design, mu, rsd, eff, groups, note) {
  er <- spec$analysis$expected_ratio
  cats <- as.character(unlist(spec$primary_outcome$categories %||% names(er)))
  if (!length(cats)) cats <- paste0("class_", seq_len(as.integer(a$nominal_categories %||% 3)))
  K <- length(cats)
  base <- if (length(er)) { w <- as.numeric(unlist(er[cats])); if (anyNA(w)) stop("analysis.expected_ratio labels must match primary_outcome.categories.", call. = FALSE); w }
          else { pr <- as.numeric(unlist(a$nominal_probabilities %||% rep(1, K))); if (length(pr) != K) rep(1, K) else pr }
  expected <- base / sum(base)
  d <- as.numeric(a$segregation_distortion %||% 0)
  truep <- expected * exp(d * (seq_len(K) == 1)); truep <- truep / sum(truep)
  if (d != 0) note(sprintf("nominal outcome: ratio distorted (first class log-odds %+g): %s instead of %s", d,
                           paste(round(truep, 3), collapse = ":"), paste(round(expected, 3), collapse = ":")))
  lor <- log(as.numeric(a$nominal_odds_ratio %||% 3))
  coef <- if (eff != 0 && K > 1) lor / (eff / rsd) * (seq_len(K) - 1) / (K - 1) else rep(0, K)
  alpha <- log(truep / truep[1])
  probs <- function(l) { z <- alpha + coef * l; e <- exp(z - max(z)); e / sum(e) }
  ok <- !is.na(y_mean)
  lp <- (y_mean - mu) / rsd
  y <- rep(NA_character_, length(y_mean))
  for (i in which(ok)) y[i] <- cats[sample.int(K, 1, prob = probs(lp[i]))]
  gp <- .suite_group_truth(groups[ok], y_design[ok], mu, rsd, function(l) as.list(stats::setNames(round(probs(l), 4), cats)))
  list(type = "nominal", y = y, model = list(
    family = "multinomial_logit", categories = as.list(cats), reference_category = cats[1],
    reference_probabilities = as.list(stats::setNames(round(truep, 4), cats)),
    expected_ratio = if (length(er)) er, true_ratio = if (length(er)) as.list(stats::setNames(round(truep, 4), cats)),
    segregation_distortion = if (length(er)) d, distorted = if (length(er)) d != 0,
    log_odds_vs_reference_per_outcome_sd = as.list(stats::setNames(round(coef, 4), cats)),
    group_probability = gp,
    note = "log P(class j) / P(reference class) = alpha_j + coef_j x lp; group probabilities at an average unit"))
}

# Secondary traits (spec secondary_outcomes [{name, unit, column, baseline_mean, residual_sd}]): each shares
# secondary_effect_share of the primary's design and random effects (in its own SDs) and has residuals correlated
# trait_residual_cor with the primary's and with each other (the within-design correlation).
.suite_secondary_traits <- function(spec, a, z0, y_mean, mu, rsd) {
  sec <- spec$secondary_outcomes %||% list()
  if (!length(sec)) return(NULL)
  K <- length(sec)
  rho <- as.numeric(a$trait_residual_cor %||% 0.5)
  share <- as.numeric(a$secondary_effect_share %||% 0.5)
  R <- matrix(rho, K + 1, K + 1); diag(R) <- 1
  if (min(eigen(R, only.values = TRUE)$values) <= 0) stop("trait_residual_cor gives a correlation matrix that is not positive definite.", call. = FALSE)
  Rs0 <- R[-1, 1, drop = FALSE]; Rss <- R[-1, -1, drop = FALSE] - Rs0 %*% t(Rs0)
  L <- chol(Rss)
  n <- length(z0); ok <- !is.na(z0)
  Z <- matrix(NA_real_, n, K)
  Z[ok, ] <- outer(z0[ok], drop(Rs0)) + matrix(stats::rnorm(sum(ok) * K), sum(ok), K) %*% L
  cols <- list(); info <- list()
  for (k in seq_len(K)) {
    s <- sec[[k]]; nm <- s$column %||% s$name %||% paste0("trait_", k)
    m_k <- as.numeric(s$baseline_mean %||% mu); sd_k <- as.numeric(s$residual_sd %||% rsd)
    cols[[nm]] <- round(m_k + share * (y_mean - mu) * sd_k / rsd + sd_k * Z[, k], 6)
    info[[nm]] <- list(name = s$name %||% nm, unit = s$unit, mean = m_k, residual_sd = sd_k, effect_share = share)
  }
  nms <- c("primary_outcome", names(cols))
  dimnames(R) <- list(nms, nms)
  list(columns = cols, truth = list(traits = info, within_design_correlation = stats::setNames(lapply(seq_len(nrow(R)), function(i) as.list(R[i, ])), nms),
    method = spec$analysis$correlation$method %||% "pearson",
    note = "traits = their mean + effect_share x the primary's design and random effects (in their SDs) + residuals with this correlation; raw correlations also carry the shared treatment effects"))
}

# expected Cohen's kappa of two raters with misclassification matrices M1, M2 (true class -> rated class) and
# true class shares pi; weights w (identity = unweighted)
.suite_expected_kappa <- function(pi, M1, M2, w = diag(length(pi))) {
  P <- t(M1) %*% diag(pi) %*% M2
  po <- sum(w * P); pe <- sum(w * outer(colSums(P), rowSums(P)))
  (po - pe) / (1 - pe)
}

# Raters (analysis.agreement {raters: [columns], scale, levels}): continuous = true value + rater bias + rater
# noise (truth: ICC2 from the variances); nominal / ordinal = the true class read through a misclassification
# matrix per rater (truth: expected pairwise kappa). The true value is the primary outcome.
.suite_ratings <- function(spec, a, truth_value, fam, note, rows = NULL) {
  ag <- spec$analysis$agreement
  if (is.null(ag)) return(NULL)
  raters <- as.character(unlist(ag$raters))
  # long format (Biodes 0.6.7: one row per unit and rater; rater_column, score_column, unit_column): the true value
  # is the unit's (its first row's simulated value); each row's score is that rater's reading of it
  long <- !length(raters) && !is.null(ag$rater_column) && !is.null(rows) && all(c(ag$rater_column, ag$unit_column %||% "") %in% names(rows))
  if (!length(raters) && !long) { note("analysis.agreement without raters or a rater_column in the template: not simulated"); return(NULL) }
  if (long) {
    unit <- as.character(rows[[ag$unit_column]]); rr <- as.character(rows[[ag$rater_column]])
    first <- match(unit, unit)
    truth_value <- truth_value[first]
    raters <- unique(rr)
  }
  scale <- ag$scale %||% if (fam$type %in% c("ordinal", "nominal", "binary")) if (fam$type == "nominal") "nominal" else "ordinal" else "continuous"
  ok <- !is.na(truth_value)
  n <- length(truth_value); cols <- list()
  if (scale == "continuous") {
    tv <- as.numeric(truth_value)
    sb <- as.numeric(a$rater_bias_sd %||% 0.3) * stats::sd(tv[ok]); se <- as.numeric(a$rater_noise_sd %||% 0.3) * stats::sd(tv[ok])
    bias <- stats::setNames(stats::rnorm(length(raters), 0, sb), raters)
    if (long) { score <- round(tv + unname(bias[rr]) + stats::rnorm(n, 0, se), 6); vt <- stats::var(tv[ok][!duplicated(unit[ok])]) }
    else { for (r in raters) cols[[r]] <- round(tv + bias[[r]] + stats::rnorm(n, 0, se), 6); vt <- stats::var(tv[ok]) }
    return(list(columns = cols, score = if (long) score, truth = list(scale = "continuous", raters = as.list(raters),
      format = if (long) "long" else "wide", true_value = if (long) "each unit's value before rating (not written)" else "primary_outcome",
      true_value_variance = vt, rater_bias_sd = sb, rater_noise_sd = se, rater_bias = as.list(round(bias, 4)),
      icc2 = vt / (vt + sb^2 + se^2), icc3 = vt / (vt + se^2),
      note = "rating = true value (primary_outcome) + rater bias + noise; ICC2 = two-way random, absolute agreement, single rater, from the planted variances")))
  }
  # categorical: true class from the primary outcome (or the continuous value cut into equal-share classes)
  lev <- as.character(unlist(ag$levels %||% spec$primary_outcome$categories))
  tc <- if (fam$type %in% c("ordinal", "nominal", "binary")) as.character(truth_value) else NULL
  if (is.null(tc)) {
    k <- max(2L, length(lev) %||% 0L); if (!length(lev)) lev <- as.character(seq_len(5)); k <- length(lev)
    cuts <- stats::quantile(as.numeric(truth_value[ok]), seq_len(k - 1) / k, names = FALSE)
    tc <- ifelse(ok, lev[findInterval(as.numeric(truth_value), cuts) + 1], NA)
  }
  if (!length(lev)) lev <- sort(unique(tc[ok]))
  if (fam$type == "ordinal" && all(grepl("^[0-9]+$", lev))) lev <- as.character(sort(as.integer(lev)))
  K <- length(lev); acc <- as.numeric(a$rater_accuracy %||% 0.8)
  M <- matrix(0, K, K, dimnames = list(lev, lev))
  for (i in seq_len(K)) {
    others <- if (scale == "ordinal") intersect(c(i - 1, i + 1), seq_len(K)) else setdiff(seq_len(K), i)
    M[i, i] <- acc; M[i, others] <- (1 - acc) / length(others)
  }
  read1 <- function(x) if (is.na(x)) NA_character_ else lev[sample.int(K, 1, prob = M[match(x, lev), ])]
  score <- NULL
  if (long) score <- vapply(tc, read1, "") else for (r in raters) cols[[r]] <- vapply(tc, read1, "")
  pi <- as.numeric(table(factor(if (long) tc[ok][!duplicated(unit[ok])] else tc[ok], lev))) / sum(if (long) !duplicated(unit[ok]) else ok)
  kap <- .suite_expected_kappa(pi, M, M)
  wq <- 1 - outer(seq_len(K), seq_len(K), function(i, j) (i - j)^2) / (K - 1)^2
  if (long && fam$type %in% c("ordinal", "count", "binary") && all(grepl("^-?[0-9]+$", lev))) score <- as.integer(score)
  list(columns = cols, score = score, truth = list(scale = scale, raters = as.list(raters), levels = as.list(lev),
    format = if (long) "long" else "wide", true_class = if (long) "each unit's class before rating (not written)" else "primary_outcome",
    rater_accuracy = acc, misclassification = stats::setNames(lapply(seq_len(K), function(i) as.list(M[i, ])), lev),
    true_class_shares = as.list(stats::setNames(round(pi, 4), lev)),
    expected_kappa = kap, expected_weighted_kappa_quadratic = if (scale == "ordinal") .suite_expected_kappa(pi, M, M, wq),
    note = sprintf("each rater reads the true class through the same misclassification matrix (correct with probability %s; errors %s); expected Cohen's kappa for any pair of raters",
                   format(acc), if (scale == "ordinal") "to an adjacent score" else "spread over the other classes")))
}
