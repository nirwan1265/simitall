# Correlation among repeated measures of one unit (suite_outcome within_subject):
#   none          - only the shared unit effect (random intercept): compound symmetry, sphericity holds
#   ar1           - residuals of one unit follow AR(1) over the occasions (level order): corr = rho^|lag|
#   random_slope  - each unit gets its own slope over time (SD slope_sd_ratio x residual SD per full time range)
# ar1 and random_slope break sphericity (what Greenhouse-Geisser corrects).

# The repeated series: the occasion column (analysis.terms.time, else the first repeated factor) and, per row,
# the unit it belongs to (all layout columns but the occasion). NULL when no unit is measured more than once.
.suite_within_series <- function(meta, layout_cols, repeated, spec) {
  none <- function(why) structure(list(), reason = why)
  occ <- c(intersect(c(spec$analysis$terms$time, repeated), repeated), NA_character_)[1]
  if (is.na(occ) || !occ %in% names(meta)) return(none("the design has no repeated factor (analysis.terms.repeated_factors): units are measured once"))
  key_cols <- setdiff(intersect(layout_cols, names(meta)), occ)
  # a column unique on every row (sample_id, a row number) cannot hold a unit across occasions
  per_row <- key_cols[vapply(key_cols, function(k) !anyDuplicated(as.character(meta[[k]])), logical(1))]
  key_cols <- setdiff(key_cols, per_row)
  if (!length(key_cols)) return(none(sprintf("no layout column besides %s%s identifies a unit", occ, if (length(per_row)) sprintf(" and the per-row id %s", paste(per_row, collapse = ", ")) else "")))
  key <- do.call(paste, c(lapply(meta[key_cols], as.character), sep = "|"))
  if (!anyDuplicated(key)) return(none(sprintf("no unit is measured at more than one %s (units from %s)", occ, paste(key_cols, collapse = ", "))))
  f <- Filter(function(z) identical(z$name, occ), spec$factors %||% list())
  lev <- .suite_level_order(as.character(meta[[occ]]), if (length(f)) as.character(unlist(f[[1]]$levels)) else NULL)
  k <- match(as.character(meta[[occ]]), lev)
  # occasions on the time scale for slopes (real sampling times when the spec or labels give them)
  vals <- if (length(f)) suppressWarnings(as.numeric(unlist(f[[1]]$values))) else numeric()
  pos <- if (length(vals) == length(lev) && !anyNA(vals)) vals else suppressWarnings(as.numeric(sub("^[^0-9.-]*(-?[0-9.]+).*$", "\\1", lev)))
  if (anyNA(pos) || length(unique(pos)) < 2) pos <- seq_along(lev) - 1
  list(occasion = occ, levels = lev, index = k, position = pos[k], key = key, key_cols = key_cols, ignored_columns = as.list(per_row))
}

# Returns the residuals (n) and the slope part to add to the mean; truth describes what was done.
.suite_within_subject <- function(ws, n, rsd, args) {
  structure <- args$within_subject %||% "none"
  eps <- stats::rnorm(n, 0, rsd)
  slope_add <- numeric(n)
  if (!is.null(ws) && !length(ws)) {
    if (!identical(structure, "none"))
      warning(sprintf("within_subject = '%s' not applied: %s. Residuals are independent.", structure, attr(ws, "reason")), call. = FALSE)
    return(list(eps = eps, add = slope_add, truth = if (!identical(structure, "none")) list(structure = "independent", requested = structure, reason = attr(ws, "reason"))))
  }
  if (is.null(ws) || identical(structure, "none"))
    return(list(eps = eps, add = slope_add, truth = if (!is.null(ws)) list(structure = "random_intercept", over = ws$occasion,
      note = "units share one random effect over all occasions (compound symmetry; sphericity holds)")))
  units <- unique(ws$key)
  if (identical(structure, "ar1")) {
    rho <- as.numeric(args$within_rho %||% 0.7)
    for (u in units) {
      i <- which(ws$key == u); i <- i[order(ws$index[i])]
      lag <- c(0, diff(ws$index[i]))
      e <- numeric(length(i)); e[1] <- stats::rnorm(1)
      for (j in seq_along(i)[-1]) { r <- rho^lag[j]; e[j] <- r * e[j - 1] + sqrt(1 - r^2) * stats::rnorm(1) }
      eps[i] <- rsd * e
    }
    return(list(eps = eps, add = slope_add, truth = list(structure = "ar1", rho = rho, over = ws$occasion, lag = "level order (adjacent occasions lag 1)",
      units = length(units), unit_columns = as.list(ws$key_cols), residual_sd = rsd,
      note = sprintf("residuals of one unit: corr = %s^|lag| over %s; on top of the unit effect", format(rho), ws$occasion))))
  }
  if (identical(structure, "random_slope")) {
    ssd <- as.numeric(args$slope_sd_ratio %||% 1) * rsd
    span <- diff(range(ws$position)); if (span <= 0) span <- 1
    b <- stats::setNames(stats::rnorm(length(units), 0, ssd), units)
    slope_add <- unname(b[ws$key]) * (ws$position - min(ws$position)) / span
    return(list(eps = eps, add = slope_add, truth = list(structure = "random_slope", slope_sd = ssd, over = ws$occasion,
      slope_scale = "change over the full time range", units = length(units), unit_columns = as.list(ws$key_cols),
      note = sprintf("each unit's own slope over %s (SD %s over the full range) on top of the unit effect; variance grows with time", ws$occasion, format(signif(ssd, 3))))))
  }
  stop(sprintf("suite_outcome within_subject must be none, ar1 or random_slope, not '%s'.", structure), call. = FALSE)
}
