# Breeding populations scored against Mendelian expectations. The simulator
# records which founder every chromosome segment came from (ancestry tracts);
# these functions measure each line's genome composition and heterozygosity
# from those tracts and compare the population mean with theory, so a breeder
# can see, for example, how much recurrent-parent genome a BC2 line carries.

#' Expected genome composition for a breeding scheme
#'
#' @param scheme `"F2"`, `"RIL"`, `"NIL"`, `"DH"`, or `"MAGIC"`.
#' @param backcrosses Backcross generations to parent 1 (NIL).
#' @param self_generations Selfing generations after the F1 (or after the
#'   last backcross for NIL, after the funnel for MAGIC).
#' @param n_founders Founders for MAGIC.
#' @param ril_mating `"SSD"` (selfing) or `"SIB"`; theory is given for SSD only.
#' @return List with expected `parent1_genome` and `heterozygosity`
#'   (`NA` where no simple closed form applies).
#' @export
breeding_expectations <- function(scheme, backcrosses = 0L, self_generations = 0L, n_founders = 2L, ril_mating = "SSD") {
  scheme <- toupper(scheme)
  switch(scheme,
    F2 = list(parent1_genome = 0.5, heterozygosity = 0.5),
    RIL = list(parent1_genome = 0.5, heterozygosity = if (toupper(ril_mating) == "SSD") 0.5^self_generations else NA_real_),
    NIL = list(parent1_genome = 1 - 0.5^(backcrosses + 1), heterozygosity = 0.5^(backcrosses + self_generations)),
    DH = list(parent1_genome = 0.5, heterozygosity = 0),
    MAGIC = list(parent1_genome = 1 / n_founders, heterozygosity = 0.5^self_generations),
    list(parent1_genome = NA_real_, heterozygosity = NA_real_))
}

.simitall_read_ancestry <- function(ancestry) {
  if (is.character(ancestry)) ancestry <- utils::read.delim(ancestry, stringsAsFactors = FALSE)
  need <- c("sample", "haplotype", "chromosome", "start", "end", "founder")
  if (any(!need %in% names(ancestry))) stop("Ancestry table needs columns: ", paste(need, collapse = ", "))
  ancestry
}

# Per-line genome fractions by founder and heterozygosity, from tracts.
.simitall_line_composition <- function(ancestry, founders) {
  lengths <- tapply(ancestry$end, ancestry$chromosome, max)
  total <- sum(lengths)
  rows <- lapply(split(ancestry, ancestry$sample), function(x) {
    founder_bp <- stats::setNames(numeric(length(founders)), founders)
    het_bp <- 0
    for (chr in unique(x$chromosome)) {
      h1 <- x[x$chromosome == chr & x$haplotype == 1L, , drop = FALSE]; h1 <- h1[order(h1$start), , drop = FALSE]
      h2 <- x[x$chromosome == chr & x$haplotype == 2L, , drop = FALSE]; h2 <- h2[order(h2$start), , drop = FALSE]
      cuts <- sort(unique(c(h1$start, h2$start, lengths[[chr]] + 1)))
      seg_len <- diff(cuts); seg_start <- cuts[-length(cuts)]
      f1 <- h1$founder[findInterval(seg_start, h1$start)]
      f2 <- h2$founder[findInterval(seg_start, h2$start)]
      het_bp <- het_bp + sum(seg_len[f1 != f2])
      for (f in founders) founder_bp[f] <- founder_bp[f] + (sum(seg_len[f1 == f]) + sum(seg_len[f2 == f])) / 2
    }
    data.frame(sample = x$sample[1L], t(founder_bp / total), heterozygosity = het_bp / total,
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Summarise a simulated breeding population against Mendelian theory
#'
#' @param ancestry Ancestry tract table or path (from simulate_breeding()).
#' @param scheme Breeding scheme (see breeding_expectations()).
#' @param parent1 Founder id of parent 1 (the recurrent parent for NIL).
#' @param backcrosses,self_generations,n_founders,ril_mating Scheme settings.
#' @return List with `lines` (per-line composition) and `metrics` (one row:
#'   observed vs expected parent-1 genome and heterozygosity, with a truth check).
#' @export
summarize_breeding_population <- function(ancestry, scheme, parent1 = "hap1", backcrosses = 0L, self_generations = 0L,
                                          n_founders = 2L, ril_mating = "SSD") {
  ancestry <- .simitall_read_ancestry(ancestry)
  founders <- sort(unique(ancestry$founder))
  if (!parent1 %in% founders) stop("parent1 '", parent1, "' is not among the founders in the ancestry table")
  founders <- c(parent1, setdiff(founders, parent1))
  lines <- .simitall_line_composition(ancestry, founders)
  expected <- breeding_expectations(scheme, backcrosses, self_generations, n_founders, ril_mating)
  check <- function(observed, target) {
    if (is.na(target)) return(c(se = stats::sd(observed) / sqrt(length(observed)), z = NA_real_, ok = NA))
    se <- stats::sd(observed) / sqrt(length(observed))
    diff <- mean(observed) - target
    z <- if (se > 0) diff / se else if (abs(diff) < 1e-9) 0 else Inf
    c(se = se, z = z, ok = abs(z) <= 3.5)
  }
  p1 <- check(lines[[parent1]], expected$parent1_genome)
  het <- check(lines$heterozygosity, expected$heterozygosity)
  agree <- c(p1[["ok"]], het[["ok"]]); agree <- agree[!is.na(agree)]
  metrics <- data.frame(
    analysis = "breeding population vs Mendelian theory", scheme = toupper(scheme), lines = nrow(lines),
    backcrosses = backcrosses, self_generations = self_generations,
    parent1_genome = mean(lines[[parent1]]), expected_parent1_genome = expected$parent1_genome, parent1_genome_z = p1[["z"]],
    heterozygosity = mean(lines$heterozygosity), expected_heterozygosity = expected$heterozygosity, heterozygosity_z = het[["z"]],
    founder_balance = if (length(founders) > 2L) max(abs(colMeans(lines[founders]) - 1 / length(founders))) else NA_real_,
    agrees_with_theory = if (length(agree)) all(as.logical(agree)) else NA,
    stringsAsFactors = FALSE)
  list(lines = lines, metrics = metrics)
}

#' Plot a breeding population: chromosome painting and genome composition
#'
#' @param ancestry Ancestry tract table or path.
#' @param summary Output of summarize_breeding_population().
#' @param figure_dir Output directory for PNG files.
#' @param n_lines Number of lines to paint.
#' @return Invisibly, the written paths.
#' @export
plot_breeding_population <- function(ancestry, summary, figure_dir, n_lines = 20L) {
  ancestry <- .simitall_read_ancestry(ancestry)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  founders <- setdiff(names(summary$lines), c("sample", "heterozygosity"))
  palette <- stats::setNames(c("#2F6B5A", "#D9984A", "#7A9CC6", "#B85252", "#8E6BB8", "#4FA3A5", "#C9A227", "#6B6B6B")[seq_along(founders)], founders)
  # Show lines spread across the whole range of parent-1 genome (not just the
  # top), so donor segments and variation between lines are visible.
  ordered <- summary$lines$sample[order(-summary$lines[[founders[1L]]])]
  shown <- ordered[unique(round(seq(1, length(ordered), length.out = min(n_lines, length(ordered)))))]
  chroms <- unique(ancestry$chromosome); lengths <- tapply(ancestry$end, ancestry$chromosome, max)[chroms]
  offsets <- c(0, cumsum(lengths + max(lengths) * 0.08))[seq_along(chroms)]
  painting <- file.path(figure_dir, "chromosome_painting.png")
  grDevices::png(painting, width = 1800, height = 1100, res = 170)
  graphics::par(mar = c(4, 7, 4.5, 1))
  graphics::plot(NA, xlim = c(0, max(offsets + lengths)), ylim = c(0.5, length(shown) + 0.5), axes = FALSE, xlab = "", ylab = "",
                 main = paste0("Which parent each chromosome segment came from (", length(shown), " lines across the range, simulation truth)"), cex.main = 0.9)
  for (i in seq_along(shown)) {
    x <- ancestry[ancestry$sample == shown[i], , drop = FALSE]
    for (h in 1:2) {
      y <- i + (h - 1.5) * 0.36
      xh <- x[x$haplotype == h, , drop = FALSE]
      off <- offsets[match(xh$chromosome, chroms)]
      graphics::rect(off + xh$start, y - 0.16, off + xh$end, y + 0.16, col = palette[xh$founder], border = NA)
    }
  }
  graphics::axis(2, at = seq_along(shown), labels = shown, las = 1, cex.axis = 0.6, tick = FALSE)
  graphics::axis(1, at = offsets + lengths / 2, labels = chroms, tick = FALSE, cex.axis = 0.8)
  graphics::legend("top", legend = paste0(founders, if (length(founders) == 2L) c(" (parent 1)", " (parent 2)") else ""),
                   fill = palette, border = NA, cex = 0.75, bty = "n", horiz = TRUE, inset = c(0, -0.045), xpd = TRUE)
  grDevices::dev.off()

  composition <- file.path(figure_dir, "genome_composition.png")
  m <- summary$metrics
  grDevices::png(composition, width = 1600, height = 800, res = 170)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  graphics::hist(summary$lines[[founders[1L]]], breaks = 20, col = "#2F6B5A", border = "white", xlim = c(0, 1),
                 main = paste0(founders[1L], " genome per line"), xlab = "Fraction of genome", cex.main = 0.9)
  if (!is.na(m$expected_parent1_genome)) graphics::abline(v = m$expected_parent1_genome, lwd = 2, lty = 2, col = "#B85252")
  graphics::hist(summary$lines$heterozygosity, breaks = 20, col = "#7A9CC6", border = "white", xlim = c(0, 1),
                 main = "Heterozygosity per line", xlab = "Fraction of genome heterozygous", cex.main = 0.9)
  if (!is.na(m$expected_heterozygosity)) graphics::abline(v = m$expected_heterozygosity, lwd = 2, lty = 2, col = "#B85252")
  graphics::mtext("Dashed line: Mendelian expectation", side = 1, line = -1.2, outer = TRUE, cex = 0.7)
  grDevices::dev.off()
  invisible(c(painting, composition))
}
