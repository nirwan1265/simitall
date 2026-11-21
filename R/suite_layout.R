# Build the field layout from Biodes' study details (spec) for a chosen design.
#
# How each design is laid out is a row of inst/suite/layouts.tsv; which spec
# factor plays which part (time, check, environment, period) is
# inst/suite/factor_roles.tsv. Nothing about a particular question is in this
# file: factor names, level counts, reps, block size and locations all come from
# the spec, and anything the spec does not give is either a declared default
# (reps, block size: listed in attr(book, "assumptions")) or a question back to
# the user (level counts are never invented).
#
# The layout's columns are Biodes' generic model terms (treatment, treatment2,
# block, rep, row, col, unit, subject, time, env, period, sequence, check), so
# write_design_package() lists them in analysis.terms.layout_columns.

#' Build a field layout from a Biodes study spec
#'
#' @param spec Biodes study spec (a list, e.g. `studyforge_ask(q)$spec`).
#' @param design Catalog design id (e.g. `"rcbd"`, `"split"`, `"met_fa"`).
#' @param seed Randomization seed.
#' @return A data frame (the field book) with attributes `design`, `seed`,
#'   `assumptions` and `factor_map`; or, when the spec lacks something the
#'   layout needs, a `simitall_layout_needs` object listing what to ask.
#' @export
simitall_layout_from_spec <- function(spec, design, seed = 42) {
  rules <- .suite_layout_rule(design)
  needs <- character()
  assumptions <- character()
  if (is.null(rules) || identical(rules$treatments, "none")) {
    return(.suite_needs(sprintf("Design '%s' has no layout rule%s.", design,
                                if (!is.null(rules) && nzchar(rules$note)) paste0(": ", rules$note) else " in inst/suite/layouts.tsv")))
  }
  f <- .suite_spec_factors(spec)
  treat <- f$treatment
  reps <- suppressWarnings(as.integer(spec$design$reps %||% NA))
  if (is.na(reps) || reps < 1) {
    reps <- as.integer(rules$default_reps)
    assumptions <- c(assumptions, sprintf("reps = %d (default for %s; the spec gives no design.reps)", reps, design))
  }
  need_levels <- function(x, what) {
    if (is.null(x)) {
      needs <<- c(needs, sprintf("%s: the spec has no such factor. Add it, e.g. answers = list(factors = \"<name>:<levels>, ...\").", what))
    } else if (is.na(x$n)) {
      needs <<- c(needs, sprintf("number of levels of '%s' (%s).", x$name, what))
    }
  }
  two <- rules$treatments %in% c("crossed", "split", "strip", "pairs", "coded")
  if (rules$treatments == "split") {
    main <- f$main %||% treat[[2]]
    other <- Filter(function(x) !identical(x$name, main$name), treat)[1][[1]]
    t1 <- if (rules$main_plot == "treatment2") other else main
    t2 <- if (rules$main_plot == "treatment2") main else other
  } else {
    t1 <- treat[1][[1]]
    t2 <- treat[2][[1]]
  }
  need_levels(t1, "the treatment factor")
  if (two) need_levels(t2, "the second treatment factor (treatment2)")
  if (rules$treatments == "entries_checks") need_levels(f$check, "the check entries")
  if (rules$within == "time") need_levels(f$time, "time points measured on each unit")
  n_env <- 1L
  if (rules$per_env == "yes") {
    n_env <- suppressWarnings(as.integer(spec$design$locations %||% f$env$n %||% length(spec$environments) %||% NA))
    if (is.na(n_env) || n_env < 1) needs <- c(needs, "number of environments (design.locations).")
  }
  if (length(needs)) return(.suite_needs(needs))

  set.seed(seed)
  lev <- function(x, coded = FALSE) {
    if (coded) return(format(seq(-1, 1, length.out = x$n), trim = TRUE))
    if (length(x$levels) == x$n) return(as.character(x$levels))
    paste0(x$name, "_", seq_len(x$n))
  }
  coded <- rules$treatments == "coded"
  L1 <- lev(t1, coded)
  L2 <- if (two) lev(t2, coded)
  cells <- switch(rules$treatments,
    pairs = {
      if (length(L1) < 2) return(.suite_needs("a diallel needs at least 2 parents."))
      cb <- utils::combn(L1, 2)
      data.frame(treatment = cb[1, ], treatment2 = cb[2, ], stringsAsFactors = FALSE)
    },
    entries_checks = data.frame(treatment = L1, stringsAsFactors = FALSE),
    one = data.frame(treatment = L1, stringsAsFactors = FALSE),
    expand.grid(treatment = L1, treatment2 = L2, stringsAsFactors = FALSE)
  )
  if (rules$within %in% c("crossover", "stepped") && rules$treatments != "one") {
    return(.suite_needs(sprintf("layout rule for '%s': crossover layouts take one treatment factor.", design)))
  }

  one_env <- function(env_label) {
    b <- .suite_arrange(cells, rules, reps, spec, design, checks = if (rules$treatments == "entries_checks") lev(f$check),
                        n_time = if (rules$within == "time") f$time$n, time_levels = if (rules$within == "time") lev(f$time))
    if (is.character(b)) return(b)
    assumptions <<- unique(c(assumptions, attr(b, "assumptions")))
    if (!is.null(env_label)) {
      b$env <- env_label
      for (k in intersect(c("block", "rep", "unit", "subject"), names(b))) b[[k]] <- paste0(env_label, "_", b[[k]])
    }
    b
  }
  envs <- if (rules$per_env == "yes") {
    if (!is.null(f$env) && length(f$env$levels) == n_env) as.character(f$env$levels) else paste0("env_", seq_len(n_env))
  }
  parts <- if (is.null(envs)) list(one_env(NULL)) else lapply(envs, one_env)
  bad <- Filter(is.character, parts)
  if (length(bad)) return(.suite_needs(bad[[1]]))
  book <- do.call(rbind, parts)
  book <- cbind(plot = seq_len(nrow(book)), book)
  rownames(book) <- NULL
  if (length(rules$note) && nzchar(rules$note)) assumptions <- c(assumptions, paste0(design, " layout: ", rules$note))
  attr(book, "design") <- design
  attr(book, "seed") <- seed
  attr(book, "assumptions") <- assumptions
  attr(book, "factor_map") <- Filter(Negate(is.null), list(treatment = t1$name, treatment2 = if (two) t2$name,
                                                            time = if (rules$within == "time") f$time$name,
                                                            check = if (rules$treatments == "entries_checks") f$check$name,
                                                            env = if (rules$per_env == "yes") f$env$name %||% "design.locations"))
  book
}

#' @export
print.simitall_layout_needs <- function(x, ...) {
  cat("The layout cannot be built yet; the study details lack:\n", paste0("  - ", x$needs, "\n"), sep = "")
  invisible(x)
}

.suite_needs <- function(x) structure(list(needs = x), class = "simitall_layout_needs")

.suite_layout_rule <- function(design) {
  p <- system.file("suite", "layouts.tsv", package = "simitall")
  if (!nzchar(p)) p <- file.path("inst", "suite", "layouts.tsv")
  tab <- utils::read.delim(p, colClasses = "character", na.strings = character(), quote = "")
  r <- tab[tab$design == design, , drop = FALSE]
  if (!nrow(r)) return(NULL)
  as.list(r[1, ])
}

# Spec factors by part: treatment factors in spec order, plus time, check, env,
# period (factor_roles.tsv) and the main-plot factor (handoff.json main_plot_applied_to).
.suite_spec_factors <- function(spec) {
  p <- system.file("suite", "factor_roles.tsv", package = "simitall")
  if (!nzchar(p)) p <- file.path("inst", "suite", "factor_roles.tsv")
  roles <- utils::read.delim(p, colClasses = "character", quote = "")
  main_applied <- unlist(.suite_rules("handoff")$main_plot_applied_to)
  out <- list(treatment = list())
  for (x in spec$factors %||% list()) {
    nm <- x$name %||% "factor"
    n <- suppressWarnings(as.integer(x$n_levels %||% if (length(x$levels)) length(x$levels) else NA))
    item <- list(name = nm, n = n, levels = unlist(x$levels))
    hit <- roles$role[vapply(roles$pattern, function(pt) grepl(pt, tolower(nm), perl = TRUE), logical(1))][1]
    if (!is.na(hit) && is.null(out[[hit]])) {
      out[[hit]] <- item
      next
    }
    if (is.null(out$main) && any(unlist(x$applied_to) %in% main_applied)) out$main <- item
    out$treatment[[length(out$treatment) + 1]] <- item
  }
  out
}

# Arrange treatment cells into reps/blocks, then add within-unit structure and ids.
.suite_arrange <- function(cells, rules, reps, spec, design, checks = NULL, n_time = NULL, time_levels = NULL) {
  assumptions <- character()
  shuffle <- function(d) d[sample(nrow(d)), , drop = FALSE]
  n <- nrow(cells)
  b <- switch(rules$blocking,
    none = {
      d <- shuffle(cells[rep(seq_len(n), reps), , drop = FALSE])
      d$rep <- as.character(stats::ave(seq_len(nrow(d)), do.call(paste, d[intersect(c("treatment", "treatment2"), names(d))]), FUN = seq_along))
      d
    },
    complete = {
      if (rules$treatments == "entries_checks") {
        # entries spread once over the blocks, every check in every block
        ent <- sample(cells$treatment)
        grp <- rep_len(seq_len(reps), length(ent))
        do.call(rbind, lapply(seq_len(reps), function(k) {
          d <- data.frame(treatment = c(ent[grp == k], checks), stringsAsFactors = FALSE)
          d$check <- ifelse(d$treatment %in% checks, "check", "entry")
          d <- shuffle(d)
          d$block <- as.character(k)
          d$rep <- as.character(k)
          d
        }))
      } else if (rules$treatments == "split") {
        do.call(rbind, lapply(seq_len(reps), function(k) {
          mains <- sample(unique(cells$treatment))
          do.call(rbind, lapply(seq_along(mains), function(m) {
            d <- shuffle(cells[cells$treatment == mains[m], , drop = FALSE])
            d$block <- as.character(k); d$rep <- as.character(k); d$main_plot <- paste0(k, "_", m)
            d
          }))
        }))
      } else if (rules$treatments == "strip") {
        do.call(rbind, lapply(seq_len(reps), function(k) {
          r1 <- sample(unique(cells$treatment)); r2 <- sample(unique(cells$treatment2))
          d <- cells
          d$row <- as.character(match(d$treatment, r1)); d$col <- as.character(match(d$treatment2, r2))
          d$block <- as.character(k); d$rep <- as.character(k)
          d[order(as.integer(d$row), as.integer(d$col)), , drop = FALSE]
        }))
      } else {
        do.call(rbind, lapply(seq_len(reps), function(k) {
          d <- shuffle(cells); d$block <- as.character(k); d$rep <- as.character(k); d
        }))
      }
    },
    incomplete = {
      k <- suppressWarnings(as.integer(spec$design$block_size %||% NA))
      if (is.na(k) || k < 2 || k >= n) {
        dv <- seq_len(n)[n %% seq_len(n) == 0 & seq_len(n) > 1 & seq_len(n) < n]
        k <- if (length(dv)) dv[which.min(abs(dv - sqrt(n)))] else max(2L, round(sqrt(n)))
        assumptions <- c(assumptions, sprintf("block size = %d (divisor of %d entries nearest its square root; the spec gives no design.block_size)", k, n))
      }
      do.call(rbind, lapply(seq_len(reps), function(r) {
        d <- shuffle(cells)
        d$rep <- as.character(r)
        d$block <- paste0(r, "_", ceiling(seq_len(nrow(d)) / k))
        d
      }))
    },
    latin = ,
    graeco = {
      if (rules$blocking == "graeco" && (n != length(unique(cells$treatment))^2 || length(unique(cells$treatment)) %% 2 == 0 ||
                                          length(unique(cells$treatment)) != length(unique(cells$treatment2)))) {
        return(sprintf("%s needs both factors with the same odd number of levels.", design))
      }
      t1 <- unique(cells$treatment)
      t2 <- if (rules$blocking == "graeco") unique(cells$treatment2)
      m <- length(t1)
      do.call(rbind, lapply(seq_len(reps), function(r) {
        g <- expand.grid(i = 0:(m - 1), j = 0:(m - 1))
        pr <- sample(m); pc <- sample(m); p1 <- sample(t1)
        d <- data.frame(row = paste0(r, "_", pr[g$i + 1]), col = paste0(r, "_", pc[g$j + 1]),
                        treatment = p1[(g$i + g$j) %% m + 1], stringsAsFactors = FALSE)
        if (!is.null(t2)) d$treatment2 <- sample(t2)[(g$i + 2 * g$j) %% m + 1]
        d$rep <- as.character(r)
        d
      }))
    },
    row_column = {
      plots <- if (rules$treatments == "entries_checks") n + length(checks) * reps else n
      if (rules$treatments == "entries_checks") {
        rules2 <- rules; rules2$blocking <- "complete"
        d <- .suite_arrange(cells, rules2, reps, spec, design, checks = checks)
        nc <- suppressWarnings(as.integer(spec$design$n_cols %||% NA))
        if (is.na(nc)) nc <- ceiling(sqrt(nrow(d)))
        d$row <- as.character(ceiling(seq_len(nrow(d)) / nc)); d$col <- as.character((seq_len(nrow(d)) - 1) %% nc + 1)
        d
      } else {
        nc <- suppressWarnings(as.integer(spec$design$n_cols %||% NA))
        if (is.na(nc)) {
          nc <- ceiling(sqrt(plots))
          assumptions <- c(assumptions, sprintf("field grid %d columns wide per rep (the spec gives no design.n_cols)", nc))
        }
        do.call(rbind, lapply(seq_len(reps), function(r) {
          d <- shuffle(cells)
          d$rep <- as.character(r)
          d$row <- paste0(r, "_", ceiling(seq_len(nrow(d)) / nc)); d$col <- paste0(r, "_", (seq_len(nrow(d)) - 1) %% nc + 1)
          d
        }))
      }
    }
  )
  if (is.character(b)) return(b)

  # the experimental unit of each row before within-unit structure
  b$unit <- paste0("u", seq_len(nrow(b)))
  if (rules$within == "time") {
    b <- b[rep(seq_len(nrow(b)), each = n_time), , drop = FALSE]
    b$time <- rep(time_levels, length.out = nrow(b))
  } else if (rules$within == "crossover") {
    tl <- unique(cells$treatment)
    m <- length(tl)
    if (m < 2) return("a crossover needs at least 2 treatments.")
    seqs <- lapply(0:(m - 1), function(s) tl[(s + 0:(m - 1)) %% m + 1])
    blocks <- if (rules$blocking == "complete") seq_len(reps) else 1L
    per_seq <- if (rules$blocking == "complete") 1L else reps
    b <- do.call(rbind, lapply(blocks, function(k) do.call(rbind, lapply(seq_len(m), function(s) do.call(rbind, lapply(seq_len(per_seq), function(j) {
      id <- paste0(if (rules$blocking == "complete") paste0(k, "_"), s, "_", j)
      data.frame(unit = id, sequence = paste(seqs[[s]], collapse = "-"), period = paste0("p", seq_len(m)),
                 treatment = seqs[[s]], block = as.character(k), rep = as.character(j), stringsAsFactors = FALSE)
    }))))))
  } else if (rules$within == "stepped") {
    tl <- unique(cells$treatment)
    if (length(tl) != 2) return("a stepped-wedge layout needs a two-level treatment (control, intervention).")
    steps <- reps
    b <- do.call(rbind, lapply(seq_len(steps), function(s) do.call(rbind, lapply(seq_len(steps + 1), function(p) {
      data.frame(unit = paste0("c", s), period = paste0("p", p), treatment = tl[if (p > s) 2 else 1], rep = "1", stringsAsFactors = FALSE)
    }))))
  }
  # subjects and repeated observations within the unit
  if (rules$subject == "nested") {
    sp <- as.integer(rules$subjects_per_unit)
    b <- b[rep(seq_len(nrow(b)), each = sp), , drop = FALSE]
    b$subject <- paste0(b$unit, "_s", stats::ave(seq_len(nrow(b)), paste(b$unit, b$time %||% ""), FUN = seq_along))
  } else if (rules$subject == "yes") {
    b$subject <- b$unit
  } else if (rules$subject %in% c("block", "block_batched")) {
    b$subject <- b$block
    if (rules$subject == "block_batched") b$block <- as.character(ceiling(as.integer(b$block) / 2))
  }
  obs <- as.integer(rules$obs_per_unit %||% "1")
  if (!is.na(obs) && obs > 1) b <- b[rep(seq_len(nrow(b)), each = obs), , drop = FALSE]
  keep <- intersect(c("rep", "block", "main_plot", "row", "col", "unit", "subject", "sequence", "period", "time",
                      "treatment", "treatment2", "check"), names(b))
  if (rules$unit == "no" && !rules$within %in% c("time", "crossover", "stepped") && rules$subject == "no") keep <- setdiff(keep, "unit")
  out <- b[, keep, drop = FALSE]
  rownames(out) <- NULL
  attr(out, "assumptions") <- assumptions
  out
}
