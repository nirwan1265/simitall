# Survival outcomes (primary_outcome.type = "survival"; Bioana docs/contracts/survival.md).
# suite_outcome survival_model:
#   weibull_ph     - proportional hazards, Weibull baseline (shape weibull_shape): truth has the HR and the time ratio
#   exponential    - weibull_ph with shape 1
#   lognormal_aft  - log-normal accelerated failure time (time ratios; hazards not proportional)
#   delayed_effect - HR 1 until change_time, then the HR (non-PH; weighted log-rank fh_0_1, RMST)
#   crossing       - the HR until change_time, then 1/HR (curves cross)
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

.suite_survival <- function(a, spec, y_obs, y_design, y_mean, mu, rsd, eff, groups, cond_effects, random_sd, note) {
  model <- a$survival_model %||% "weibull_ph"
  models <- c("weibull_ph", "exponential", "lognormal_aft", "delayed_effect", "crossing")
  if (!model %in% models) stop(sprintf("suite_outcome survival_model must be one of %s, not '%s'.", paste(models, collapse = ", "), model), call. = FALSE)
  ok <- !is.na(y_obs); n <- length(y_obs)
  ss <- spec$analysis$survival %||% list()
  ge <- spec$analysis$goal_effect$survival %||% list()
  med <- as.numeric(a$median_time)
  k <- if (model == "exponential") 1 else as.numeric(a$weibull_shape %||% 1.5)
  if (model == "lognormal_aft") k <- NA_real_
  # hazard ratio of the full treatment effect (reference -> last level)
  hr <- as.numeric(a$hazard_ratio %||% 0.6); hr_src <- "suite_outcome hazard_ratio"
  if (!is.null(spec$power$target_hazard_ratio)) { hr <- as.numeric(spec$power$target_hazard_ratio); hr_src <- "power.target_hazard_ratio" }
  if (!is.null(ge$size) && identical(ge$scale, "hazard_ratio")) { hr <- as.numeric(ge$size); hr_src <- "analysis.goal_effect.survival (hazard_ratio)" }
  if (!is.null(ge$size) && identical(ge$scale, "time_ratio") && !is.na(k)) { hr <- as.numeric(ge$size)^(-k); hr_src <- "analysis.goal_effect.survival (time_ratio^-shape)" }
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
  draw1 <- function(et, eo) {
    E <- stats::rexp(length(et)) / exp(eo)
    switch(model,
      weibull_ph = , exponential = H0inv(E / exp(et)),
      delayed_effect = ifelse(E <= H0(t0), H0inv(E), H0inv(H0(t0) + (E - H0(t0)) / exp(et))),
      crossing = { a1 <- exp(et) * H0(t0); ifelse(E <= a1, H0inv(E / exp(et)), H0inv(H0(t0) + (E - a1) * exp(et))) })
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
  follow <- as.numeric(a$follow_up %||% stats::quantile(tmin[ok], 1 - a$censored_fraction, names = FALSE))
  lossf <- as.numeric(a$loss_fraction %||% 0)
  lost <- stats::runif(n) < lossf
  closs <- ifelse(lost, stats::runif(n, 0, follow), Inf)
  y <- pmin(tmin, follow, closs)
  ev <- ifelse(tmin <= pmin(follow, closs), cause, 0L)
  status <- ifelse(ev > 0, "event", ifelse(closs < follow & closs < tmin, "lost", "administrative"))
  w <- .suite_num_list(ss$life_table)[1]
  if (!is.na(w) && w > 0) y <- pmin(ceiling(y / w) * w, ceiling(follow / w) * w)
  y[!ok] <- NA; ev[!ok] <- NA

  # truth per group (conditional on an average cluster): S(t), median, landmarks, RMST
  tau <- .suite_num_list(ss$rmst_tau)[1]; if (is.na(tau)) tau <- follow
  lms <- .suite_num_list(ss$landmarks); if (!length(lms)) lms <- signif(c(med, follow), 3)
  grid <- seq(0, max(tau, lms, follow) * 1.001, length.out = 4001)
  Hc <- if (length(comp)) length(comp) * (grid / lam_c)^kc else 0 * grid
  gid <- unique(groups[ok]); gt <- list()
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
      crossing = exp(e) * H0(pmin(grid, t0)) + exp(-e) * pmax(H0(grid) - H0(t0), 0))
    S1 <- exp(-H1); S <- exp(-H1 - Hc)
    at <- function(f, x) stats::approx(grid, f, xout = x, rule = 2)$y
    rm <- { i <- grid <= tau; sum(diff(grid[i]) * (S[i][-1] + S[i][-sum(i)]) / 2) }
    one <- list(median = if (any(S1 <= 0.5)) round(grid[which(S1 <= 0.5)[1]], 4) else NA, log_hazard = e,
                survival_at = stats::setNames(as.list(round(at(S, lms), 4)), format(lms)), rmst = round(rm, 4))
    if (length(comp)) {
      h1 <- c(0, diff(H1)); cif <- cumsum(h1 * c(1, S[-length(S)]))
      one$cif_event_at <- stats::setNames(as.list(round(at(cif, lms), 4)), format(lms))
      one$median_note <- "median of the event of interest alone (net of competing events)"
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
    r })
  frailty <- if (model != "lognormal_aft" && length(random_sd)) lapply(random_sd, function(s) list(distribution = "lognormal",
    sd_log_hazard = round(abs(coef) * s / rsd, 4), variance = round((coef * s / rsd)^2, 4)))
  n_ok <- sum(ok)
  m <- list(family = model, proportional_hazards = ph, median_time = med,
            weibull_shape = if (!is.na(k)) k, hazard_ratio = if (model != "lognormal_aft") hr, hazard_ratio_source = if (model != "lognormal_aft") hr_src,
            log_hazard_per_outcome_sd = if (model != "lognormal_aft") coef,
            change_time = if (model %in% c("delayed_effect", "crossing")) t0,
            time_sd = if (model == "lognormal_aft") a$time_sd, scale = if (model == "lognormal_aft") list(sd_outcome = sd_out),
            follow_up = follow, censoring = list(
              administrative = list(follow_up = follow, n = sum(status[ok] == "administrative"), target_fraction = a$censored_fraction),
              loss = list(loss_fraction = lossf, n = sum(status[ok] == "lost"), times = "uniform on (0, follow_up)")),
            n_events = sum(ev[ok] > 0),
            event_codes = if (length(codes)) codes, competing = if (length(comp)) list(codes = as.list(comp), model = "Weibull cause-specific hazard, no treatment effect",
              median_alone = med * as.numeric(a$competing_median_ratio %||% 2), shape = kc, cause_specific_hazard_ratio = 1),
            frailty = frailty, life_table_width = if (!is.na(w)) w,
            landmarks = as.list(lms), rmst_tau = tau, groups = gt, ratios = ratios,
            note = paste0("primary_outcome = time to the first event or censoring; event 1 = the event, 0 = censored",
                          if (length(comp)) ", 2+ = competing events" else "",
                          ". Group truth is conditional on an average cluster (no frailty)."))
  m <- m[!vapply(m, is.null, logical(1))]
  list(type = "survival", y = round(y, 3), event = as.integer(ev), model = m)
}
