#' Read common biallelic marker frequencies from a VCF
#'
#' Reads `INFO/AF` from a plain or gzip-compressed VCF. The returned table is
#' intentionally limited to marker frequencies; it does not phase haplotypes or
#' claim that any marker is causal.
#'
#' @param vcf Path to a VCF or VCF.GZ file.
#' @param maf_min,maf_max Inclusive alternate-allele frequency bounds.
#' @param max_markers Maximum number of evenly spaced eligible markers to keep.
#'
#' @return A data frame with marker IDs, positions, and alternate-allele
#'   frequencies.
#' @export
read_human_marker_frequencies <- function(
    vcf,
    maf_min = 0.10,
    maf_max = 0.40,
    max_markers = 24L) {
  if (!is.character(vcf) || length(vcf) != 1L || !file.exists(vcf)) {
    stop("vcf must be one existing VCF or VCF.GZ path")
  }
  if (!is.finite(maf_min) || !is.finite(maf_max) || maf_min < 0 || maf_max > 1 || maf_min >= maf_max) {
    stop("maf_min and maf_max must satisfy 0 <= maf_min < maf_max <= 1")
  }
  connection <- if (grepl("\\.gz$", vcf, ignore.case = TRUE)) gzfile(vcf, "rt") else file(vcf, "rt")
  on.exit(close(connection), add = TRUE)
  lines <- readLines(connection, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]
  if (!length(lines)) stop("VCF contains no variant records")
  fields <- strsplit(lines, "\t", fixed = TRUE)
  keep <- vapply(fields, length, integer(1)) >= 8L
  fields <- fields[keep]
  markers <- data.frame(
    chromosome = vapply(fields, `[[`, character(1), 1L),
    position_bp = as.integer(vapply(fields, `[[`, character(1), 2L)),
    id = vapply(fields, `[[`, character(1), 3L),
    ref = vapply(fields, `[[`, character(1), 4L),
    alt = vapply(fields, `[[`, character(1), 5L),
    info = vapply(fields, `[[`, character(1), 8L),
    stringsAsFactors = FALSE
  )
  af_match <- regexec("(?:^|;)AF=([0-9.]+)(?:;|$)", markers$info, perl = TRUE)
  af_value <- regmatches(markers$info, af_match)
  markers$alternate_allele_frequency <- suppressWarnings(as.numeric(vapply(
    af_value, function(x) if (length(x) > 1L) x[2L] else NA_character_, character(1)
  )))
  markers <- markers[
    is.finite(markers$alternate_allele_frequency) &
      markers$alternate_allele_frequency >= maf_min &
      markers$alternate_allele_frequency <= maf_max &
      nchar(markers$ref) == 1L & nchar(markers$alt) == 1L &
      !grepl(",", markers$alt, fixed = TRUE),
    c("chromosome", "position_bp", "id", "ref", "alt", "alternate_allele_frequency"),
    drop = FALSE
  ]
  if (nrow(markers) < 2L) stop("VCF has fewer than two eligible common biallelic markers")
  max_markers <- as.integer(max_markers)
  if (!is.finite(max_markers) || max_markers < 2L) stop("max_markers must be at least two")
  if (nrow(markers) > max_markers) {
    index <- unique(round(seq(1L, nrow(markers), length.out = max_markers)))
    markers <- markers[index, , drop = FALSE]
  }
  markers$marker_id <- paste0("marker_", seq_len(nrow(markers)))
  rownames(markers) <- NULL
  markers
}

.simitall_human_gamete <- function(genotype) stats::rbinom(length(genotype), 1L, genotype / 2)
.simitall_human_child <- function(parent_a, parent_b) {
  .simitall_human_gamete(parent_a) + .simitall_human_gamete(parent_b)
}
.simitall_human_founders <- function(n, frequencies) {
  matrix(stats::rbinom(n * length(frequencies), 2L, rep(frequencies, each = n)), nrow = n, byrow = FALSE)
}

#' Simulate named human relationship templates from marker frequencies
#'
#' Creates unrelated, second-cousin, and first-cousin-descendant comparison
#' groups. Sex is recorded in the pedigree truth table. Genotypes are simulated
#' from independent common-marker frequencies, so this is a teaching simulation,
#' not a phased-haplotype or clinical model.
#'
#' @param markers Marker-frequency table returned by
#'   `read_human_marker_frequencies()`.
#' @param groups Canonical relationship templates: `unrelated`,
#'   `distantly_related`, and/or `first_cousin_descendant`.
#' @param n_per_group Number of descendant samples per group.
#' @param seed Random seed.
#' @param baseline_probability Baseline illustrative trait probability. Under
#'   `risk_model = "recessive"` this is the expected probability for a child of
#'   unrelated parents, and the model intercept is calibrated to hit it.
#' @param risk_model `"illustrative"` (default, unchanged behaviour) adds an
#'   explicit relatedness term to the liability. `"recessive"` has no
#'   relatedness term: risk comes only from children inheriting two copies of a
#'   risk allele at simulated recessive loci, so any difference between
#'   relationship groups emerges from inheritance through the pedigree.
#' @param n_risk_loci,risk_allele_frequency,risk_effect Recessive model only:
#'   number of simulated unlinked risk loci, their risk-allele frequency, and
#'   the log-odds increase per locus that is homozygous for the risk allele.
#'
#' @return A list containing `individuals`, `pedigree`, `marker_truth`, and
#'   `genotypes`. Under the recessive model, `individuals` also contains
#'   `n_homozygous_risk_loci`, `expected_probability` (the exact model value
#'   for the child's theoretical inbreeding coefficient), and a simulated
#'   `affected` status; the list also contains `model`.
#' @export
simulate_human_pedigree_groups <- function(
    markers,
    groups = c("unrelated", "distantly_related", "first_cousin_descendant"),
    n_per_group = 100L,
    seed = 1L,
    baseline_probability = 0.05,
    risk_model = c("illustrative", "recessive"),
    n_risk_loci = 10L,
    risk_allele_frequency = 0.05,
    risk_effect = 4) {
  risk_model <- match.arg(risk_model)
  required <- c("marker_id", "alternate_allele_frequency")
  if (!is.data.frame(markers) || !all(required %in% names(markers))) {
    stop("markers must contain marker_id and alternate_allele_frequency columns")
  }
  frequencies <- markers$alternate_allele_frequency
  if (length(frequencies) < 2L || any(!is.finite(frequencies) | frequencies <= 0 | frequencies >= 1)) {
    stop("marker frequencies must be finite values strictly between zero and one")
  }
  allowed <- c("unrelated", "distantly_related", "first_cousin_descendant")
  groups <- unique(as.character(groups))
  if (!length(groups) || any(!groups %in% allowed)) {
    stop("groups must contain unrelated, distantly_related, and/or first_cousin_descendant")
  }
  n_per_group <- as.integer(n_per_group)
  if (!is.finite(n_per_group) || n_per_group < 1L) stop("n_per_group must be at least one")
  if (!is.finite(baseline_probability) || baseline_probability <= 0 || baseline_probability >= 1) {
    stop("baseline_probability must be between zero and one")
  }
  set.seed(seed)
  L <- length(frequencies)
  recessive <- identical(risk_model, "recessive")
  if (recessive) {
    n_risk_loci <- as.integer(n_risk_loci)
    if (!is.finite(n_risk_loci) || n_risk_loci < 1L) stop("n_risk_loci must be at least one")
    if (!is.finite(risk_allele_frequency) || risk_allele_frequency <= 0 || risk_allele_frequency >= 0.5) {
      stop("risk_allele_frequency must be between zero and 0.5")
    }
    if (!is.finite(risk_effect) || risk_effect <= 0) stop("risk_effect must be positive")
    risk_frequencies <- rep(risk_allele_frequency, n_risk_loci)
    intercept <- .simitall_recessive_intercept(baseline_probability, risk_frequencies, risk_effect)
    # Risk loci are transmitted through the same pedigree as the markers.
    frequencies <- c(frequencies, risk_frequencies)
    effects <- rep(0, L)
  } else {
    effects <- c(0.35, stats::rnorm(L - 1L, 0.08, 0.03))
  }
  pedigree <- list()
  individual_rows <- list()
  genotype_rows <- list()
  pedigree_id <- 0L
  add_person <- function(id, sex, sire = NA_character_, dam = NA_character_, generation, group, role) {
    pedigree_id <<- pedigree_id + 1L
    pedigree[[pedigree_id]] <<- data.frame(
      id = id, sex = sex, sire = sire, dam = dam, generation = generation,
      group = group, role = role, stringsAsFactors = FALSE
    )
  }
  draw_founder <- function(id, sex, group, generation = 0L, role = "founder") {
    add_person(id, sex, generation = generation, group = group, role = role)
    .simitall_human_founders(1L, frequencies)[1L, ]
  }
  make_person <- function(id, sex, sire_id, sire, dam_id, dam, generation, group, role) {
    add_person(id, sex, sire = sire_id, dam = dam_id, generation = generation, group = group, role = role)
    .simitall_human_child(sire, dam)
  }
  sample_number <- 0L
  for (group in groups) {
    for (i in seq_len(n_per_group)) {
      key <- sprintf("%s_%03d", group, i)
      if (identical(group, "unrelated")) {
        father <- draw_founder(paste0(key, "_father"), "male", group)
        mother <- draw_founder(paste0(key, "_mother"), "female", group)
        child <- make_person(paste0(key, "_child"), sample(c("female", "male"), 1L),
          paste0(key, "_father"), father, paste0(key, "_mother"), mother, 1L, group, "target_descendant")
        F <- 0; relationship <- 0
      } else {
        gf <- draw_founder(paste0(key, "_grandfather"), "male", group)
        gm <- draw_founder(paste0(key, "_grandmother"), "female", group)
        sibling_a <- make_person(paste0(key, "_sibling_a"), "female", paste0(key, "_grandfather"), gf, paste0(key, "_grandmother"), gm, 1L, group, "shared_ancestor_branch")
        sibling_b <- make_person(paste0(key, "_sibling_b"), "male", paste0(key, "_grandfather"), gf, paste0(key, "_grandmother"), gm, 1L, group, "shared_ancestor_branch")
        partner_a <- draw_founder(paste0(key, "_partner_a"), "male", group)
        partner_b <- draw_founder(paste0(key, "_partner_b"), "female", group)
        cousin_a <- make_person(paste0(key, "_cousin_a"), "female", paste0(key, "_partner_a"), partner_a, paste0(key, "_sibling_a"), sibling_a, 2L, group, "cousin_branch")
        cousin_b <- make_person(paste0(key, "_cousin_b"), "male", paste0(key, "_sibling_b"), sibling_b, paste0(key, "_partner_b"), partner_b, 2L, group, "cousin_branch")
        if (identical(group, "first_cousin_descendant")) {
          child <- make_person(paste0(key, "_child"), sample(c("female", "male"), 1L), paste0(key, "_cousin_b"), cousin_b, paste0(key, "_cousin_a"), cousin_a, 3L, group, "target_descendant")
          F <- 1 / 16; relationship <- 1 / 8
        } else {
          external_a <- draw_founder(paste0(key, "_external_a"), "male", group)
          external_b <- draw_founder(paste0(key, "_external_b"), "female", group)
          second_a <- make_person(paste0(key, "_second_a"), "female", paste0(key, "_external_a"), external_a, paste0(key, "_cousin_a"), cousin_a, 3L, group, "second_cousin_branch")
          second_b <- make_person(paste0(key, "_second_b"), "male", paste0(key, "_cousin_b"), cousin_b, paste0(key, "_external_b"), external_b, 3L, group, "second_cousin_branch")
          child <- make_person(paste0(key, "_child"), sample(c("female", "male"), 1L), paste0(key, "_second_b"), second_b, paste0(key, "_second_a"), second_a, 4L, group, "target_descendant")
          F <- 1 / 64; relationship <- 1 / 32
        }
      }
      sample_number <- sample_number + 1L
      if (recessive) {
        # No relatedness term: risk depends only on inherited risk genotypes.
        markers_child <- child[seq_len(L)]
        n_hom <- sum(child[L + seq_len(n_risk_loci)] == 2L)
        genetic_score <- intercept + risk_effect * n_hom
        probability <- stats::plogis(genetic_score)
        individual_rows[[sample_number]] <- data.frame(
          sample = paste0(key, "_child"), group = group,
          theoretical_parental_relatedness = relationship, theoretical_inbreeding_F = F,
          marker_homozygosity = mean(markers_child == 0L | markers_child == 2L),
          synthetic_liability = genetic_score, synthetic_trait_probability = probability,
          n_homozygous_risk_loci = n_hom,
          expected_probability = .simitall_recessive_expected(F, intercept, risk_frequencies, risk_effect),
          affected = stats::rbinom(1L, 1L, probability),
          stringsAsFactors = FALSE
        )
      } else {
        # This is an explicit, illustrative relatedness component. It makes the
        # requested comparison visible without claiming a real human risk effect.
        genetic_score <- sum(child * effects) + 12 * F + stats::rnorm(1L, 0, 0.20)
        probability <- stats::plogis(stats::qlogis(baseline_probability) + genetic_score)
        individual_rows[[sample_number]] <- data.frame(
          sample = paste0(key, "_child"), group = group,
          theoretical_parental_relatedness = relationship, theoretical_inbreeding_F = F,
          marker_homozygosity = mean(child == 0L | child == 2L),
        synthetic_liability = genetic_score, synthetic_trait_probability = probability,
          stringsAsFactors = FALSE
        )
      }
      genotype_rows[[sample_number]] <- child
    }
  }
  individuals <- do.call(rbind, individual_rows)
  genotype <- do.call(rbind, genotype_rows)
  if (!recessive) {
    colnames(genotype) <- markers$marker_id
    return(list(
      individuals = individuals,
      pedigree = do.call(rbind, pedigree),
      marker_truth = transform(markers, synthetic_effect = effects),
      genotypes = genotype
    ))
  }
  risk_ids <- sprintf("risk_locus_%02d", seq_len(n_risk_loci))
  colnames(genotype) <- c(markers$marker_id, risk_ids)
  marker_truth <- rbind(
    data.frame(marker_id = markers$marker_id, alternate_allele_frequency = markers$alternate_allele_frequency,
               role = "neutral_marker", recessive_effect = 0, stringsAsFactors = FALSE),
    data.frame(marker_id = risk_ids, alternate_allele_frequency = risk_frequencies,
               role = "recessive_risk_locus", recessive_effect = risk_effect, stringsAsFactors = FALSE)
  )
  list(
    individuals = individuals,
    pedigree = do.call(rbind, pedigree),
    marker_truth = marker_truth,
    genotypes = genotype,
    model = list(risk_model = "recessive", intercept = intercept, n_risk_loci = n_risk_loci,
                 risk_allele_frequency = risk_allele_frequency, risk_effect = risk_effect,
                 baseline_probability = baseline_probability)
  )
}

# Exact distribution of the number of loci homozygous for the risk allele in a
# child with inbreeding coefficient F. Per locus, P(hom) = q^2 (1 - F) + q F.
.simitall_recessive_hom_distribution <- function(F, risk_frequencies) {
  p <- risk_frequencies^2 * (1 - F) + risk_frequencies * F
  dist <- 1
  for (pj in p) dist <- c(dist * (1 - pj), 0) + c(0, dist * pj)
  dist
}

.simitall_recessive_expected <- function(F, intercept, risk_frequencies, risk_effect) {
  dist <- .simitall_recessive_hom_distribution(F, risk_frequencies)
  sum(dist * stats::plogis(intercept + risk_effect * (seq_along(dist) - 1L)))
}

# Calibrate the intercept so a child of unrelated parents (F = 0) has exactly
# the requested baseline probability.
.simitall_recessive_intercept <- function(baseline_probability, risk_frequencies, risk_effect) {
  f <- function(a) .simitall_recessive_expected(0, a, risk_frequencies, risk_effect) - baseline_probability
  stats::uniroot(f, c(-30, 10), tol = 1e-10)$root
}

#' Simulate anonymous human marker frequencies
#'
#' Generates a marker-frequency table with the same columns as
#' [read_human_marker_frequencies()], for questions that name no region or
#' data set. The markers are synthetic and carry no biological meaning.
#'
#' @param n_markers Number of independent markers.
#' @param maf_min,maf_max Alternate-allele frequency range.
#' @param seed Random seed.
#' @return A data frame of synthetic markers.
#' @export
simulate_human_marker_frequencies <- function(n_markers = 24L, maf_min = 0.10, maf_max = 0.40, seed = 1L) {
  n_markers <- as.integer(n_markers)
  if (!is.finite(n_markers) || n_markers < 2L) stop("n_markers must be at least two")
  if (!is.finite(maf_min) || !is.finite(maf_max) || maf_min <= 0 || maf_max >= 1 || maf_min >= maf_max) {
    stop("maf_min and maf_max must satisfy 0 < maf_min < maf_max < 1")
  }
  set.seed(seed)
  data.frame(
    chromosome = "synthetic", position_bp = seq_len(n_markers) * 10000L,
    id = sprintf("synthetic_marker_%02d", seq_len(n_markers)), ref = "A", alt = "G",
    alternate_allele_frequency = round(stats::runif(n_markers, maf_min, maf_max), 4),
    marker_id = paste0("marker_", seq_len(n_markers)),
    stringsAsFactors = FALSE
  )
}

#' Summarise disease probability by relationship group
#'
#' For output of `simulate_human_pedigree_groups(risk_model = "recessive")`.
#' Compares the simulated mean probability in each group with the exact model
#' expectation. The truth check is an exact binomial test on how many children
#' are homozygous for a risk allele at any locus; `agrees_with_theory` is
#' `TRUE` when its p-value is at least 0.001.
#'
#' @param run Output of [simulate_human_pedigree_groups()] with the recessive model.
#' @return A data frame with one row per relationship group.
#' @export
summarize_family_risk <- function(run) {
  ind <- run$individuals
  if (is.null(run$model) || !"expected_probability" %in% names(ind)) {
    stop("summarize_family_risk() needs a run simulated with risk_model = 'recessive'")
  }
  risk_frequencies <- rep(run$model$risk_allele_frequency, run$model$n_risk_loci)
  groups <- unique(ind$group)
  rows <- lapply(groups, function(g) {
    x <- ind[ind$group == g, , drop = FALSE]
    F <- x$theoretical_inbreeding_F[1L]
    # Truth check: the number of children homozygous for a risk allele at any
    # locus is exactly Binomial(n, p) under the model, so an exact test stays
    # calibrated even when cases are rare (normal-theory checks do not).
    p_any <- 1 - .simitall_recessive_hom_distribution(F, risk_frequencies)[1L]
    observed <- sum(x$n_homozygous_risk_loci > 0L)
    check_p <- stats::binom.test(observed, nrow(x), p_any)$p.value
    data.frame(
      group = g, children = nrow(x), inbreeding_F = F,
      expected_probability = x$expected_probability[1L],
      simulated_probability = mean(x$synthetic_trait_probability),
      simulated_se = stats::sd(x$synthetic_trait_probability) / sqrt(nrow(x)),
      expected_risk_homozygous_children = p_any * nrow(x),
      observed_risk_homozygous_children = observed,
      truth_check_p = check_p,
      agrees_with_theory = check_p >= 0.001,
      affected = sum(x$affected), affected_rate = mean(x$affected),
      mean_marker_homozygosity = mean(x$marker_homozygosity),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  ref <- out$expected_probability[out$group == "unrelated"]
  out$relative_risk_vs_unrelated <- if (length(ref)) out$expected_probability / ref else NA_real_
  out
}

#' Plot disease probability by relationship group
#'
#' Bars show the exact model expectation; points with 95% intervals show the
#' simulated families, so agreement between them is visible at a glance.
#'
#' @param summary Output of [summarize_family_risk()].
#' @param path PNG path.
#' @param condition Label for the plotted condition.
#' @return Invisibly, `path`.
#' @export
plot_family_risk <- function(summary, path, condition = "condition") {
  pretty <- c(unrelated = "Unrelated parents", distantly_related = "Second cousins",
              first_cousin_descendant = "First cousins")
  labels <- unname(pretty[summary$group]); labels[is.na(labels)] <- summary$group[is.na(labels)]
  pct <- 100 * summary$expected_probability
  sim <- 100 * summary$simulated_probability
  half <- 100 * 1.96 * summary$simulated_se
  top <- max(pct, sim + half) * 1.3
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = 1600, height = 1000, res = 170)
  on.exit(grDevices::dev.off(), add = TRUE)
  mids <- graphics::barplot(pct, names.arg = labels, ylim = c(0, top),
                            col = c("#3D7A8A", "#D9984A", "#B85252")[seq_along(pct)], border = NA,
                            ylab = "Probability (%)", main = paste0("Simulated ", condition, " probability by parents' relationship"))
  graphics::arrows(mids, pmax(0, sim - half), mids, sim + half, angle = 90, code = 3, length = 0.05)
  graphics::points(mids, sim, pch = 21, bg = "white", cex = 1.3)
  # Labels sit above both the bar and its interval so they never overlap.
  graphics::text(mids, pmax(pct, sim + half), sprintf("%.2f%%  (x%.2f)", pct, summary$relative_risk_vs_unrelated),
                 pos = 3, cex = 0.95, font = 2)
  graphics::mtext("Bars: exact model expectation. Points: simulated families with 95% intervals. Synthetic model, not clinical risk.",
                  side = 1, line = 3, cex = 0.75)
  invisible(path)
}

#' Write human-pedigree comparison figures
#'
#' @param individuals Output `individuals` table from
#'   `simulate_human_pedigree_groups()`.
#' @param out_dir Output directory for PNG files.
#' @param prefix File prefix.
#'
#' @return Invisibly returns generated figure paths.
#' @export
plot_human_pedigree_diagnostics <- function(individuals, out_dir, prefix = "human_pedigree") {
  needed <- c("group", "synthetic_trait_probability", "marker_homozygosity")
  if (!is.data.frame(individuals) || !all(needed %in% names(individuals))) {
    stop("individuals must contain group, synthetic_trait_probability, and marker_homozygosity")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  groups <- unique(individuals$group)
  pretty <- c(unrelated = "Unrelated", distantly_related = "Distantly related", first_cousin_descendant = "First-cousin descendants")
  labels <- unname(pretty[groups]); labels[is.na(labels)] <- groups[is.na(labels)]
  probability_path <- file.path(out_dir, paste0(prefix, "_probability_comparison.png"))
  distribution_path <- file.path(out_dir, paste0(prefix, "_trait_distribution.png"))
  homozygosity_path <- file.path(out_dir, paste0(prefix, "_homozygosity.png"))
  summary <- stats::aggregate(individuals$synthetic_trait_probability, list(group = individuals$group), mean)
  spread <- stats::aggregate(individuals$synthetic_trait_probability, list(group = individuals$group), stats::sd)
  n <- table(individuals$group)
  summary$se <- spread$x / sqrt(as.numeric(n[summary$group]))
  summary_labels <- unname(pretty[summary$group])
  summary_labels[is.na(summary_labels)] <- summary$group[is.na(summary_labels)]
  png(probability_path, width = 1500, height = 950, res = 170)
  mids <- barplot(summary$x, names.arg = summary_labels, ylim = c(0, max(0.05, summary$x + 2 * summary$se)), col = c("#3D7A8A", "#D9984A", "#B85252"), ylab = "Mean illustrative synthetic probability", main = "Synthetic probability by relationship group")
  arrows(mids, pmax(0, summary$x - 1.96 * summary$se), mids, summary$x + 1.96 * summary$se, angle = 90, code = 3, length = 0.06)
  mtext("Error bars: 95% normal-approximation intervals; not clinical risk estimates", side = 1, line = 3, cex = 0.8)
  dev.off()
  png(distribution_path, width = 1500, height = 950, res = 170)
  boxplot(synthetic_trait_probability ~ group, data = individuals, names = labels, col = c("#A7D3E0", "#F2D19B", "#E9B2A8"), ylab = "Illustrative synthetic probability", main = "Synthetic trait-probability distribution")
  dev.off()
  png(homozygosity_path, width = 1500, height = 950, res = 170)
  boxplot(marker_homozygosity ~ group, data = individuals, names = labels, col = c("#A7D3E0", "#F2D19B", "#E9B2A8"), ylab = "Fraction of markers homozygous", main = "Marker homozygosity by relationship group")
  dev.off()
  invisible(c(probability = probability_path, distribution = distribution_path, homozygosity = homozygosity_path))
}
