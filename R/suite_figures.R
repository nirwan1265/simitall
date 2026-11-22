#' Figures for a Biodes design folder (Biopaint, first version)
#'
#' Draws the standard figures from whatever files are in a design folder; the
#' design, factors and columns come from `study_spec.json` and the model, so
#' the same code serves every design:
#'
#' - `layout.png`: the field / bench / cage layout coloured by treatment
#'   (`field_layout.csv`);
#' - `outcome.png`: the outcome per experimental unit by treatment, with means
#'   and 95% intervals, coloured by the second factor and split by time or
#'   environment when the design has them (simulated or real
#'   `sample_metadata_template.tsv`);
#' - `recovery.png`: fitted effects with 95% intervals against the simulated
#'   truth (`recovery.tsv`);
#' - `power.png`: power against replication per model term (`power_sweep.tsv`);
#' - `volcano.png`: per-gene results for count assays (`de_results.tsv`).
#'
#' In an interactive session (RStudio, R console) the figures are also shown
#' in the plot window, all on one page; `separate = TRUE` gives one window per
#' figure. `plot(r)` on a [simitall_from_question()]
#' result and `plot(ps)` on a [simitall_power_sweep()] result show them again.
#'
#' Biopaint will move to its own package; this is the version inside Simitall.
#'
#' @param dir Design folder (or a `simitall_suite_run`).
#' @param outdir Where to write the PNGs. Default `<dir>/figures`.
#' @param which Figures to draw (default: all that the folder has data for).
#' @param show Show the figures in the plot window (default: in interactive
#'   sessions): all on one page, so they are visible in any plot window.
#' @param separate Show each figure on its own instead: a new plot window per
#'   figure (quartz, windows, X11), or one plot each in RStudio's Plots history.
#' @param save Write the PNG files.
#' @return Invisibly, the paths written (named by figure), with the ggplot
#'   objects in `attr(, "plots")`.
#' @export
biopaint_biodes <- function(dir, outdir = NULL, which = c("layout", "outcome", "recovery", "power", "volcano"),
                            show = interactive(), save = TRUE, separate = FALSE) {
  if (inherits(dir, "simitall_suite_run")) dir <- dir$dir
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Figures need ggplot2: install.packages('ggplot2')", call. = FALSE)
  outdir <- outdir %||% file.path(dir, "figures")
  if (isTRUE(save)) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  spec <- jsonlite::read_json(file.path(dir, "study_spec.json"), simplifyVector = FALSE)
  layout <- utils::read.csv(file.path(dir, "field_layout.csv"), colClasses = "character", check.names = FALSE)
  p <- tryCatch(.suite_parse_model(spec$analysis$model, spec, layout), error = function(e) NULL)
  cols <- .suite_figure_columns(p, layout, spec)
  title <- sprintf("%s: %s", spec$design$chosen %||% "design", spec$analysis$model %||% "")
  paths <- character()
  plots <- list()
  write_png <- save
  save <- function(plot, name, w = 7, h = 5) {
    plots[[name]] <<- plot
    if (isTRUE(write_png)) {
      f <- file.path(outdir, paste0(name, ".png"))
      ggplot2::ggsave(f, plot, width = w, height = h, dpi = 150, bg = "white")
      paths[name] <<- f
    }
  }
  if ("layout" %in% which) {
    g <- .suite_plot_layout(layout, cols, title)
    if (!is.null(g)) save(g$plot, "layout", g$w, g$h)
  }
  mpath <- file.path(dir, "sample_metadata_template.tsv")
  if ("outcome" %in% which && file.exists(mpath)) {
    meta <- utils::read.delim(mpath, check.names = FALSE, stringsAsFactors = FALSE)
    if ("primary_outcome" %in% names(meta) && any(!is.na(suppressWarnings(as.numeric(meta$primary_outcome))))) {
      save(.suite_plot_outcome(meta, cols, spec), "outcome")
    }
  }
  rpath <- file.path(dir, "recovery.tsv")
  if ("recovery" %in% which && file.exists(rpath)) save(.suite_plot_recovery(utils::read.delim(rpath, stringsAsFactors = FALSE)), "recovery")
  ppath <- file.path(dir, "power_sweep.tsv")
  if ("power" %in% which && file.exists(ppath)) {
    pj <- tryCatch(jsonlite::read_json(file.path(dir, "power_sweep.json")), error = function(e) list())
    save(.suite_plot_power(utils::read.delim(ppath, stringsAsFactors = FALSE), pj$target_power %||% 0.8, pj$alpha %||% 0.05), "power")
  }
  vpath <- file.path(dir, "de_results.tsv")
  if ("volcano" %in% which && file.exists(vpath)) {
    truth <- tryCatch(jsonlite::read_json(file.path(dir, "sim_truth.json"))$counts$de_genes, error = function(e) NULL)
    save(.suite_plot_volcano(utils::read.delim(vpath, stringsAsFactors = FALSE), truth), "volcano")
  }
  if (length(paths)) message("biopaint: wrote ", paste(basename(paths), collapse = ", "), " in ", outdir)
  if (isTRUE(show)) .suite_show_plots(plots, separate = separate)
  attr(paths, "plots") <- plots
  invisible(paths)
}

# Show figures. Default: all on one page (a grid), so every figure is visible in
# any plot window (RStudio Plots pane, quartz, windows, X11). separate = TRUE:
# one window per figure (a new device each), or one plot each in RStudio's history.
.suite_show_plots <- function(plots, separate = FALSE) {
  if (!length(plots)) return(invisible(plots))
  if (isTRUE(separate) || length(plots) == 1) {
    rstudio <- identical(names(grDevices::dev.cur()), "RStudioGD") || nzchar(Sys.getenv("RSTUDIO"))
    for (i in seq_along(plots)) {
      if (!rstudio && i > 1 && interactive()) grDevices::dev.new()
      print(plots[[i]])
    }
    return(invisible(plots))
  }
  n <- length(plots)
  nc <- if (n <= 2) n else if (n <= 4) 2 else 3
  nr <- ceiling(n / nc)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nr, nc)))
  for (i in seq_len(n)) {
    print(plots[[i]], vp = grid::viewport(layout.pos.row = (i - 1) %/% nc + 1, layout.pos.col = (i - 1) %% nc + 1))
  }
  grid::popViewport()
  message("biopaint: ", n, " figures on one page (", paste(names(plots), collapse = ", "),
          "). Enlarge the plot window, or plot(r, separate = TRUE) for one window per figure, or plot(r, which = \"outcome\") for one.")
  invisible(plots)
}

#' Show the figures of a suite run in the plot window
#'
#' @param x Result of [simitall_from_question()].
#' @param which Figures to show (default: all that the folder has data for).
#' @param separate One plot window per figure instead of one page.
#' @param ... Unused.
#' @export
plot.simitall_suite_run <- function(x, which = c("layout", "outcome", "recovery", "power", "volcano"), separate = FALSE, ...) {
  if (is.null(x$dir)) stop("This run has no design folder (it stopped before the layout was written).", call. = FALSE)
  invisible(biopaint_biodes(x$dir, which = which, show = TRUE, save = FALSE, separate = separate))
}

#' Show the power figure of a power sweep in the plot window
#'
#' @param x Result of [simitall_power_sweep()].
#' @param ... Unused.
#' @export
plot.simitall_power_sweep <- function(x, ...) {
  dir <- attr(x, "dir")
  if (is.null(dir)) stop("This power sweep has no design folder attached.", call. = FALSE)
  invisible(biopaint_biodes(dir, which = "power", show = TRUE, save = FALSE))
}

# Which layout column plays which part in the figures.
.suite_figure_columns <- function(p, layout, spec) {
  role_col <- function(var) if (!is.null(p) && var %in% names(p$col)) p$col[[var]] else if (var %in% names(layout)) var else NA
  treat_vars <- if (is.null(p)) intersect(c("treatment", "treatment2"), names(layout)) else names(p$role)[p$role == "treatment"]
  treat_vars <- setdiff(treat_vars, c("time", "check"))
  meaning <- spec$layout$factor_map %||% spec$analysis$terms %||% list()
  lab <- function(v) as.character(unlist(meaning[[v]])[1] %||% v)
  list(
    treatment = if (length(treat_vars)) role_col(treat_vars[1]) else NA,
    treatment_label = if (length(treat_vars)) lab(treat_vars[1]) else "treatment",
    treatment2 = if (length(treat_vars) > 1) role_col(treat_vars[2]) else NA,
    treatment2_label = if (length(treat_vars) > 1) lab(treat_vars[2]) else NA,
    facet = c(time = role_col("time"), env = role_col("env"))[!is.na(c(role_col("time"), role_col("env")))][1],
    block = c(role_col("block"), role_col("rep"))[!is.na(c(role_col("block"), role_col("rep")))][1],
    row = role_col("row"), col = role_col("col"), env = role_col("env")
  )
}

.suite_num_suffix <- function(x) suppressWarnings(as.numeric(sub("^.*_", "", as.character(x))))

.suite_plot_layout <- function(layout, cols, title) {
  if (is.na(cols$treatment)) return(NULL)
  d <- layout
  d$fill <- d[[cols$treatment]]
  # second factor written in the tile (short label: the part after the last "_", else the level)
  d$label <- if (!is.na(cols$treatment2)) { l <- as.character(d[[cols$treatment2]]); ifelse(grepl("_[^_]+$", l), sub("^.*_", "", l), l) } else ""
  many <- length(unique(d$fill)) > 24
  if (!is.na(cols$row) && !is.na(cols$col)) {
    d$x <- .suite_num_suffix(d[[cols$col]]); d$y <- .suite_num_suffix(d[[cols$row]])
    d$panel <- if ("rep" %in% names(d)) paste("rep", d$rep) else ""
  } else {
    blk <- if (!is.na(cols$block)) d[[cols$block]] else rep("all", nrow(d))
    d$y <- as.integer(factor(blk, levels = .suite_natural_sort(unique(blk))))
    d$x <- stats::ave(seq_len(nrow(d)), blk, FUN = seq_along)
    d$panel <- if (!is.na(cols$env)) d[[cols$env]] else ""
  }
  ny <- length(unique(d$y)); nx <- max(d$x, na.rm = TRUE)
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$fill)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    (if (!is.na(cols$treatment2) && nrow(d) <= 400) ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 2.4, colour = "grey15") else NULL) +
    ggplot2::scale_y_reverse(breaks = if (is.na(cols$row)) seq_len(ny) else ggplot2::waiver()) +
    ggplot2::labs(x = if (is.na(cols$col)) "position in block" else "column", y = if (is.na(cols$row)) "block" else "row",
                  fill = cols$treatment_label,
                  caption = if (!is.na(cols$treatment2)) sprintf("text in each tile: %s level", cols$treatment2_label) else NULL,
                  title = "Layout", subtitle = title) +
    ggplot2::theme_minimal(base_size = 10) + ggplot2::coord_equal()
  if (length(unique(d$panel)) > 1) g <- g + ggplot2::facet_wrap(~panel)
  if (many) g <- g + ggplot2::theme(legend.position = "none") + ggplot2::labs(caption = sprintf("%d %s levels (legend omitted)", length(unique(d$fill)), cols$treatment_label))
  list(plot = g, w = min(14, 3 + nx * 0.35 * max(1, length(unique(d$panel)))), h = min(12, 2.5 + ny * 0.45))
}

.suite_plot_outcome <- function(meta, cols, spec) {
  meta$y <- suppressWarnings(as.numeric(meta$primary_outcome))
  meta <- meta[!is.na(meta$y), , drop = FALSE]
  # one value per experimental unit: subsamples are averaged, as Bioana does
  if ("unit_id" %in% names(meta)) {
    keys <- c("unit_id", stats::na.omit(c(cols$treatment, cols$treatment2, cols$facet)))
    meta <- stats::aggregate(meta["y"], meta[keys], mean)
  }
  x <- cols$treatment
  meta$x <- factor(meta[[x]], levels = .suite_natural_sort(unique(meta[[x]])))
  many <- nlevels(meta$x) > 30
  if (many) meta$x <- stats::reorder(meta$x, meta$y, mean)
  aes_col <- if (!is.na(cols$treatment2)) ggplot2::aes(colour = .data[[cols$treatment2]]) else NULL
  if (many) {
    return(ggplot2::ggplot(meta, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::stat_summary(fun.data = ggplot2::mean_cl_normal, geom = "linerange", colour = "grey70", linewidth = 0.3) +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 0.8) +
      ggplot2::labs(x = sprintf("%s (%d levels, sorted by mean)", cols$treatment_label, nlevels(meta$x)), y = spec$primary_outcome$name %||% "outcome",
                    title = "Outcome per level", subtitle = "points: mean over units; lines: 95% interval") +
      ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(axis.text.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank()) +
      (if (!is.na(cols$facet)) ggplot2::facet_wrap(stats::as.formula(paste("~", cols$facet))) else NULL))
  }
  dodge <- ggplot2::position_dodge(width = 0.5)
  g <- ggplot2::ggplot(meta, ggplot2::aes(x = .data$x, y = .data$y, group = if (!is.na(cols$treatment2)) .data[[cols$treatment2]] else 1)) +
    aes_col +
    ggplot2::geom_point(position = if (!is.na(cols$treatment2)) ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5) else
      ggplot2::position_jitter(width = 0.15, height = 0), alpha = 0.45, size = 1.4) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_cl_normal, position = dodge, geom = "pointrange", size = 0.35) +
    ggplot2::labs(x = cols$treatment_label, y = spec$primary_outcome$name %||% "outcome",
                  colour = cols$treatment2_label, title = "Outcome per experimental unit",
                  subtitle = "points: units (subsamples averaged); bars: mean and 95% interval") +
    ggplot2::theme_minimal(base_size = 10)
  if (many) g <- g + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  else if (nlevels(meta$x) > 6) g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (!is.na(cols$facet)) g <- g + ggplot2::facet_wrap(stats::as.formula(paste("~", cols$facet)))
  g
}

.suite_plot_recovery <- function(r) {
  r$lo <- r$estimate - stats::qt(0.975, pmax(r$df, 1)) * r$se
  r$hi <- r$estimate + stats::qt(0.975, pmax(r$df, 1)) * r$se
  sub <- sprintf("%d of %d true effects inside their 95%% interval", sum(r$covered), nrow(r))
  if (nrow(r) > 40) {
    return(ggplot2::ggplot(r, ggplot2::aes(x = .data$true, y = .data$estimate, colour = .data$covered)) +
      ggplot2::geom_abline(linetype = 2, colour = "grey50") + ggplot2::geom_point(alpha = 0.7) +
      ggplot2::labs(x = "simulated (true) effect", y = "fitted effect", colour = "inside 95% CI", title = "Recovery of the simulated effects", subtitle = sub) +
      ggplot2::theme_minimal(base_size = 10))
  }
  r$term <- factor(r$term, levels = rev(r$term))
  ggplot2::ggplot(r, ggplot2::aes(y = .data$term)) +
    ggplot2::geom_pointrange(ggplot2::aes(x = .data$estimate, xmin = .data$lo, xmax = .data$hi, colour = .data$covered), size = 0.3) +
    ggplot2::geom_point(ggplot2::aes(x = .data$true), shape = 4, size = 2.5) +
    ggplot2::labs(x = "effect (x = simulated truth; point and bar = fitted, 95% CI)", y = NULL, colour = "inside 95% CI",
                  title = "Recovery of the simulated effects", subtitle = sub) +
    ggplot2::theme_minimal(base_size = 10)
}

.suite_plot_power <- function(ps, target = 0.8, alpha = 0.05) {
  ps$kind <- ifelse(is.na(ps$simulated_effect), "blocking term",
                    ifelse(ps$simulated_effect %in% c(TRUE, "TRUE"), "simulated effect", "no effect (false-positive rate)"))
  ps$label <- ifelse(ps$factor == ps$term, ps$term, paste0(ps$term, " (", ps$factor, ")"))
  ps <- ps[ps$kind != "blocking term", , drop = FALSE]
  ggplot2::ggplot(ps, ggplot2::aes(x = .data$reps, y = .data$power, colour = .data$label, linetype = .data$kind)) +
    ggplot2::geom_hline(yintercept = c(target, alpha), linetype = c(2, 3), colour = "grey50") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(0, .data$power - 2 * .data$mc_se), ymax = pmin(1, .data$power + 2 * .data$mc_se), fill = .data$label),
                         alpha = 0.12, colour = NA, show.legend = FALSE) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) + ggplot2::scale_x_continuous(breaks = unique(ps$reps)) +
    ggplot2::labs(x = "replicates (design.reps)", y = "power", colour = "model term", linetype = NULL,
                  title = "Power by replication", subtitle = sprintf("dashed: target %.0f%%; dotted: alpha %.2f; bands: +/- 2 Monte Carlo SE", 100 * target, alpha)) +
    ggplot2::theme_minimal(base_size = 10)
}

.suite_plot_volcano <- function(de, truth = NULL) {
  de$sig <- !is.na(de$padj) & de$padj < 0.05
  if (!is.null(truth)) {
    true_genes <- unique(unlist(lapply(truth, names)))
    de$status <- ifelse(de$gene_id %in% true_genes, "simulated DE gene", "null gene")
  } else de$status <- "gene"
  ggplot2::ggplot(de, ggplot2::aes(x = .data$log2fc, y = -log10(.data$p), colour = .data$sig, shape = .data$status)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.2) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), name = "FDR < 0.05") +
    ggplot2::labs(x = "log2 fold change", y = "-log10 p", shape = NULL, title = "Per-gene results",
                  subtitle = sprintf("%d of %d genes at FDR < 0.05", sum(de$sig), nrow(de))) +
    ggplot2::facet_wrap(~contrast) + ggplot2::theme_minimal(base_size = 10)
}
