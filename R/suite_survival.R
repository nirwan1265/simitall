# Survival outcomes (primary_outcome.type = "survival"; Bioana docs/contracts/survival.md).
# suite_outcome survival_model:
#   weibull_ph     - proportional hazards, Weibull baseline (shape weibull_shape): truth has the HR and the time ratio
#   exponential    - weibull_ph with shape 1
#   lognormal_aft  - log-normal accelerated failure time (time ratios; hazards not proportional)
#   delayed_effect - HR 1 until change_time, then the HR (non-PH; weighted log-rank fh_0_1, RMST)
#   crossing       - the HR until change_time, then 1/HR (curves cross)
#   early_effect   - the HR until change_time, then 1
# Biodes 0.6.6: analysis.survival.logrank_weight fh_0_1 (a stated delayed effect) -> delayed_effect, fh_1_0 (early)
# -> early_effect, unless params set survival_model; the follow-up length is read from primary_outcome.censoring
# ("administrative at 21 days ...") when follow_up is not given; analysis.survival.strata are recorded with their
# log-hazard shifts (fixed block effects).
# The treatment effect `eff` maps to the hazard ratio `hazard_ratio` (power.target_hazard_ratio, or the
# survival goal's size), so the last level of a treatment factor has that HR against its reference. Blocks and
# random terms become a shared log-normal frailty on the log hazard. Censoring: administrative at follow_up
# (default: the time that leaves censored_fraction without an event) plus random loss (loss_fraction).
# Competing events: primary_outcome.event_codes with codes 2+; each gets its own Weibull hazard, no treatment
# effect. Interval scoring (analysis.survival.life_table): times rounded up to the scoring grid.

.suite_num_list <- function(x) {
  if (is.null(x)) return(numeric())
  vapply(x, function(z) as.numeric(if (is.list(z)) z$value %||% z[[1]] else z), 0)
}

.suite_survival <- function(a, spec, y_obs, y_design, y_mean, mu, rsd, eff, groups, cond_effects, random_sd, note,
                            model_from_params = FALSE, block_effects = list(), median_default = FALSE, strata_data = NULL, cluster = NULL, hr_from_params = FALSE) {
  ss <- spec$analysis$survival %||% list()
  lw <- ss$logrank_weight
  model <- a$survival_model %||% "weibull_ph"
  if (!model_from_params && !is.null(lw)) {
    model <- switch(lw, fh_0_1 = "delayed_effect", fh_1_0 = "early_effect", model)
    note(sprintf("survival model %s: analysis.survival.logrank_weight %s (a stated %s effect)", model, lw, if (lw == "fh_0_1") "delayed" else "early"))
  }
  if (is.null(a$follow_up) && !is.null(spec$primary_outcome$censoring)) {
    fu <- suppressWarnings(as.numeric(sub(".*administrative at ([0-9.]+).*", "\\1", spec$primary_outcome$censoring)))
    if (!is.na(fu) && grepl("administrative at [0-9]", spec$primary_outcome$censoring)) a$follow_up <- fu
  }
  models <- c("weibull_ph", "exponential", "lognormal_aft", "delayed_effect", "crossing", "early_effect")
  if (!model %in% models) stop(sprintf("suite_outcome survival_model must be one of %s, not '%s'.", paste(models, collapse = ", "), model), call. = FALSE)
  ok <- !is.na(y_obs); n <- length(y_obs)
  ge <- spec$analysis$goal_effect$survival %||% list()
  med <- as.numeric(a$median_time); med_src <- if (median_default) "default" else "suite_outcome median_time"
  k <- if (model == "exponential") 1 else as.numeric(a$weibull_shape %||% 1.5)
  if (model == "lognormal_aft") k <- NA_real_
  # a stated follow-up with no stated median: the reference group's median is set so that the share without the
  # event at the end of follow-up is 1 - power.event_probability (or censored_fraction), not a fixed 30
  if (median_default && !is.null(a$follow_up)) {
    left <- if (!is.null(spec$power$event_probability)) 1 - as.numeric(spec$power$event_probability) else as.numeric(a$censored_fraction)
    left <- min(max(left, 0.01), 0.99)
    med <- if (model == "lognormal_aft") as.numeric(a$follow_up) / exp(a$time_sd * stats::qnorm(1 - left))
           else as.numeric(a$follow_up) * (log(2) / -log(left))^(1 / k)
    med_src <- sprintf("from follow-up %s: %s of the reference group without the event at its end (%s)", format(a$follow_up), format(left),
                       if (!is.null(spec$power$event_probability)) "1 - power.event_probability" else "censored_fraction")
    note(sprintf("survival median time %s (%s)", format(signif(med, 3)), med_src))
  }
  # hazard ratio of the full treatment effect (reference -> last level)
  hr <- as.numeric(a$hazard_ratio %||% 0.6); hr_src <- "suite_outcome hazard_ratio"
  if (!is.null(spec$power$target_hazard_ratio)) { hr <- as.numeric(spec$power$target_hazard_ratio); hr_src <- "power.target_hazard_ratio" }
  if (!is.null(ge$size) && identical(ge$scale, "hazard_ratio")) { hr <- as.numeric(ge$size); hr_src <- "analysis.goal_effect.survival (hazard_ratio)" }
  if (!is.null(ge$size) && identical(ge$scale, "time_ratio") && !is.na(k)) { hr <- as.numeric(ge$size)^(-k); hr_src <- "analysis.goal_effect.survival (time_ratio^-shape)" }
  # a hazard ratio set in params wins over the spec (e.g. to try the opposite direction of a stated HR)
  if (hr_from_params) { hr <- as.numeric(a$hazard_ratio); hr_src <- "params suite_outcome hazard_ratio (overrides the spec)"; note(sprintf("hazard ratio %s from params, overriding the spec", format(hr))) }
  coef <- if (eff != 0) log(hr) / (eff / rsd) else 0
  codes <- spec$primary_outcome$event_codes
  comp <- if (length(codes)) setdiff(suppressWarnings(as.integer(names(codes))), c(0L, 1L, NA)) else integer()

  # baseline (reference group, average cluster): Weibull with median `med`
  lam <- if (!is.na(k)) med / log(2)^(1 / k)
  H0 <- function(t) (t / lam)^k
  H0inv <- function(h) lam * h^(1 / k)
  t0 <- as.numeric(a$change_time %||% (med / 2))
  et <- coef * (y_design - mu) / rsd          # treatment part of the log hazard (time-varying in non-PH models)
  eo <- coef * (y_mean - y_design) / rsd      # blocks, random terms: frailty
  # strata that are not model terms (Biodes 0.6.6: stated blocks < 10 -> analysis.survival.strata, no + block):
  # each stratum its own baseline hazard, a log-hazard shift drawn with SD strata_log_hazard_sd
  # several subjects per unit (seeds in one tray cell; Biodes subsamples_per_unit, analysis.survival.cluster):
  # a shared log-normal frailty per unit
  cl_sd <- NULL
  if (!is.null(cluster)) {
    cl_sd <- as.numeric(a$cluster_log_hazard_sd %||% 0.3)
    eo <- eo + stats::rnorm(max(cluster), 0, cl_sd)[cluster]
  }
  strata_shift <- list()
  for (sc in setdiff(intersect(unlist(ss$strata), names(strata_data %||% list())), names(block_effects))) {
    lev <- unique(as.character(strata_data[[sc]]))
    sh <- stats::setNames(stats::rnorm(length(lev), 0, as.numeric(a$strata_log_hazard_sd %||% 0.3)), lev)
    eo <- eo + unname(sh[as.character(strata_data[[sc]])])
    strata_shift[[sc]] <- list(log_hazard_shift = as.list(round(sh, 4)), sd = as.numeric(a$strata_log_hazard_sd %||% 0.3),
                               note = "not in analysis.model: each stratum its own baseline hazard (fit with strata())")
  }
  draw1 <- function(et, eo) {
    E <- stats::rexp(length(et)) / exp(eo)
    switch(model,
      weibull_ph = , exponential = H0inv(E / exp(et)),
      delayed_effect = ifelse(E <= H0(t0), H0inv(E), H0inv(H0(t0) + (E - H0(t0)) / exp(et))),
      crossing = { a1 <- exp(et) * H0(t0); ifelse(E <= a1, H0inv(E / exp(et)), H0inv(H0(t0) + (E - a1) * exp(et))) },
      early_effect = { a1 <- exp(et) * H0(t0); ifelse(E <= a1, H0inv(E / exp(et)), H0inv(H0(t0) + (E - a1))) })
  }
  if (model == "lognormal_aft") {
    z <- (y_obs - mean(y_obs[ok])) / stats::sd(y_obs[ok])
    t1 <- med * exp(a$time_sd * z)
    sd_out <- stats::sd(y_obs[ok])
  } else t1 <- draw1(et, eo)
  # competing causes: Weibull (same shape, 1.5 for lognormal) with median competing_median_ratio x med
  kc <- if (is.na(k)) 1.5 else k
  lam_c <- med * as.numeric(a$competing_median_ratio %||% 2) / log(2)^(1 / kc)
  tc <- vapply(comp, function(cc) lam_c * stats::rexp(n)^(1 / kc), numeric(n))
  if (!length(comp)) tc <- matrix(numeric(), n, 0)
  tall <- cbind(t1, tc)
  tmin <- apply(tall, 1, min); cause <- c(1L, comp)[apply(tall, 1, which.min)]
  # censoring
  follow <- round(as.numeric(a$follow_up %||% stats::quantile(tmin[ok], 1 - a$censored_fraction, names = FALSE)), 3)   # times are written to 3 decimals
  lossf <- as.numeric(a$loss_fraction %||% 0)
  lost <- stats::runif(n) < lossf
  closs <- ifelse(lost, stats::runif(n, 0, follow), Inf)
  y <- pmin(tmin, follow, closs)
  ev <- ifelse(tmin <= pmin(follow, closs), cause, 0L)
  status <- ifelse(ev > 0, "event", ifelse(closs < follow & closs < tmin, "lost", "administrative"))
  w <- .suite_num_list(ss$life_table)[1]
  if (!is.na(w) && w > 0) {
    # checks every w: an event is seen at the next check; a check past the follow-up never happens, so an event
    # after the last check is censored there; a lost unit was last seen at the check before it was lost
    last <- floor(follow / w + 1e-9) * w
    seen <- ceiling(y / w - 1e-9) * w
    late <- ev > 0 & seen > last
    status[late] <- "administrative"
    ev[late] <- 0L
    y <- ifelse(ev > 0, seen, ifelse(status == "lost", floor(y / w + 1e-9) * w, last))
  }
  y[!ok] <- NA; ev[!ok] <- NA

  # truth per group (conditional on an average cluster): S(t), median, landmarks, RMST
  tau <- .suite_num_list(ss$rmst_tau)[1]; if (is.na(tau)) tau <- follow
  lms <- .suite_num_list(ss$landmarks); if (!length(lms)) lms <- signif(c(med, follow), 3)
  grid <- seq(0, max(tau, lms, follow, 4 * med) * 1.001, length.out = 8001)
  Hc <- if (length(comp)) length(comp) * (grid / lam_c)^kc else 0 * grid
  gid <- intersect(levels(factor(groups)), as.character(groups[ok])); groups <- as.character(groups); gt <- list()
  if (model == "lognormal_aft") eo <- rep(0, length(groups))
  gm <- tapply(y_design, groups, mean)
  for (g in gid) {
    e <- coef * (gm[[g]] - mu) / rsd
    H1 <- if (model == "lognormal_aft") {
      sdl <- a$time_sd * rsd / sd_out
      ml <- log(med) + a$time_sd * (gm[[g]] - mean(y_obs[ok])) / sd_out
      -stats::pnorm((log(pmax(grid, 1e-12)) - ml) / sdl, lower.tail = FALSE, log.p = TRUE)
    } else switch(model,
      weibull_ph = , exponential = exp(e) * H0(grid),
      delayed_effect = H0(pmin(grid, t0)) + exp(e) * pmax(H0(grid) - H0(t0), 0),
      crossing = exp(e) * H0(pmin(grid, t0)) + exp(-e) * pmax(H0(grid) - H0(t0), 0),
      early_effect = exp(e) * H0(pmin(grid, t0)) + pmax(H0(grid) - H0(t0), 0))
    at <- function(f, x) stats::approx(grid, f, xout = x, rule = 2)$y
    lk <- trimws(format(lms))
    rmst_of <- function(S) { i <- grid <= tau; sum(diff(grid[i]) * (S[i][-1] + S[i][-sum(i)]) / 2) }
    # measures for one frailty multiplier set fr (exp of the log-hazard shifts; 1 = an average cluster), averaged
    measures <- function(fr) {
      S1 <- rowMeans(exp(-outer(H1, fr))); S <- rowMeans(exp(-outer(H1, fr) - Hc))
      out <- list(median = if (any(S1 <= 0.5)) round(grid[which(S1 <= 0.5)[1]], 4) else NA,
                  event_free_survival_at = stats::setNames(as.list(round(at(S, lms), 4)), lk), event_free_rmst = round(rmst_of(S), 4))
      if (length(comp)) {
        h1 <- c(0, diff(H1))
        cif <- rowMeans(vapply(fr, function(f) { Sf <- exp(-f * H1 - Hc); cumsum(f * h1 * c(1, Sf[-length(Sf)])) }, numeric(length(grid))))
        out$cif_event_at <- stats::setNames(as.list(round(at(cif, lms), 4)), lk)
      }
      out
    }
    cond <- measures(1)
    one <- c(list(log_hazard = e), cond)
    # kept under the old names too (same values): free of any event, not the event of interest alone
    one$survival_at <- cond$event_free_survival_at; one$rmst <- cond$event_free_rmst
    one$measure_notes <- list(median = if (length(comp)) "median of the event of interest alone (net of competing events)" else "median time to the event",
                           event_free = if (length(comp)) "event_free_*: free of any event (the event of interest and the competing events); survival_at and rmst are the same values under their old names" else "event_free_*: free of the event; survival_at and rmst are the same values",
                           cif = if (length(comp)) "cif_event_at: cumulative incidence of the event of interest, with competing events")
    # marginal over the clusters as they were drawn in this data set (trays, cells, blocks): what the observed
    # group curves estimate; the conditional values above are for an average cluster (frailty 1)
    if (model != "lognormal_aft") {
      fr <- exp(eo[as.character(groups) == g & ok])
      if (length(fr) && stats::sd(log(fr)) > 0) one$marginal <- c(measures(fr), list(note = sprintf("averaged over the %d units' realised frailties (SD %.2f on the log hazard)", length(fr), stats::sd(log(fr)))))
    }
    gt[[g]] <- one
  }
  # hazard ratios / time ratios per level of each treatment factor (vs its reference)
  ph <- model %in% c("weibull_ph", "exponential")
  ratios <- lapply(cond_effects, function(e) { e <- unlist(e); d <- e - e[1]
    r <- list(hazard_ratio = if (model != "lognormal_aft") as.list(round(exp(coef * d / rsd), 4)))
    if (ph) r$time_ratio <- as.list(round(exp(-coef * d / rsd / k), 4))
    if (model == "lognormal_aft") r$time_ratio <- as.list(round(exp(a$time_sd * d / sd_out), 4))
    if (model == "delayed_effect") r$hazard_ratio_before_change <- as.list(rep(1, length(d)))
    if (model == "crossing") r$hazard_ratio_after_change <- as.list(round(exp(-coef * d / rsd), 4))
    if (model == "early_effect") r$hazard_ratio_after_change <- as.list(rep(1, length(d)))
    r })
  frailty <- if (model != "lognormal_aft" && length(random_sd)) lapply(random_sd, function(s) list(distribution = "lognormal",
    sd_log_hazard = round(abs(coef) * s / rsd, 4), variance = round((coef * s / rsd)^2, 4)))
  if (!is.null(cl_sd)) frailty[[spec$analysis$survival$cluster %||% "unit"]] <- list(distribution = "lognormal", sd_log_hazard = cl_sd, variance = cl_sd^2,
    level = "experimental unit (shared by its subsamples)", column = "unit_id")
  n_ok <- sum(ok)
  m <- list(family = model, proportional_hazards = ph, median_time = med, median_time_source = med_src,
            weibull_shape = if (!is.na(k)) k, hazard_ratio = if (model != "lognormal_aft") hr, hazard_ratio_source = if (model != "lognormal_aft") hr_src,
            log_hazard_per_outcome_sd = if (model != "lognormal_aft") coef,
            change_time = if (model %in% c("delayed_effect", "crossing", "early_effect")) t0,
            logrank_weight = lw, follow_up_source = if (!is.null(a$follow_up)) "suite_outcome follow_up or primary_outcome.censoring" else "quantile leaving censored_fraction",
            time_sd = if (model == "lognormal_aft") a$time_sd, scale = if (model == "lognormal_aft") list(sd_outcome = sd_out),
            follow_up = follow, censoring = list(
              administrative = list(follow_up = follow, n = sum(status[ok] == "administrative"), target_fraction = a$censored_fraction),
              loss = list(loss_fraction = lossf, n = sum(status[ok] == "lost"), times = "uniform on (0, follow_up)")),
            n_events = sum(ev[ok] > 0),
            event_codes = if (length(codes)) codes, competing = if (length(comp)) list(codes = as.list(comp), model = "Weibull cause-specific hazard, no treatment effect",
              median_alone = med * as.numeric(a$competing_median_ratio %||% 2), shape = kc, cause_specific_hazard_ratio = 1),
            frailty = frailty,
            strata = if (length(ss$strata)) c(lapply(stats::setNames(nm = intersect(unlist(ss$strata), names(block_effects))), function(b)
              list(log_hazard_shift = lapply(block_effects[[b]], function(e) round(coef * e / rsd, 4)), note = "fixed block effects: each stratum its own baseline hazard multiple")),
              strata_shift),
            life_table_width = if (!is.na(w)) w,
            landmarks = as.list(lms), rmst_tau = tau, groups = gt, ratios = ratios,
            note = paste0("primary_outcome = time to the first event or censoring; event 1 = the event, 0 = censored",
                          if (length(comp)) ", 2+ = competing events" else "",
                          ". Group truth (outcome_model.groups) is conditional on an average cluster (frailty 1); groups[].marginal averages over the realised cluster frailties, as the observed group curves do."))
  m <- m[!vapply(m, is.null, logical(1))]
  list(type = "survival", y = round(y, 3), event = as.integer(ev), model = m)
}
