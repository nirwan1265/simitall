# Genetics layer of the Biodes bridge: genotypes for the lines in the layout
# and their genetic values under a phenotype model.
#
# The lines are the levels of the model term that holds the genotype
# (inst/suite/genetics.json genotype_pattern). Genotypes come from a diversity
# panel (suite_genotype_panel) or a cross simulated with simulate_breeding()
# (suite_genotype_cross, schemes in inst/suite/mating_schemes.tsv). Genetic
# values follow a row of inst/suite/phenotype_models.tsv and are scaled to the
# requested plot-level heritability against the design's residual SD.

# Model term that holds the lines, or NA.
.suite_genotype_term <- function(spec, p) {
  rules <- .suite_rules("genetics")
  pat <- rules$genotype_pattern
  aterms <- spec$analysis$terms %||% list()
  treat <- names(p$role)[p$role == "treatment"]
  for (v in treat) {
    meaning <- aterms[[v]]
    cands <- c(if (is.character(meaning)) meaning, v, p$col[[v]])
    if (any(grepl(pat, .suite_norm(cands)))) return(v)
  }
  # treatment = first spec factor (model_terms.json from_spec_factor)
  fac <- vapply(spec$factors %||% list(), function(x) x$name %||% "", "")
  if (length(fac) && grepl(pat, .suite_norm(fac[1])) && "treatment" %in% treat) return("treatment")
  NA_character_
}

.suite_is_inbred <- function(spec, setting) {
  if (!identical(tolower(as.character(setting)), "auto")) return(list(value = isTRUE(as.logical(setting)), assumed = FALSE))
  rules <- .suite_rules("genetics")
  sys <- .suite_norm(spec$system)
  if (length(sys) && any(sys %in% unlist(rules$outbred_systems))) return(list(value = FALSE, assumed = FALSE))
  list(value = TRUE, assumed = !(length(sys) && any(sys %in% unlist(rules$inbred_systems))))
}

# Number from a spec value: 6, "6", "F6" (filial generation), or a list of names (its length).
.suite_count <- function(x, filial = FALSE) {
  if (is.null(x) || !length(x)) return(NULL)
  if (length(unlist(x)) > 1) return(length(unlist(x)))
  s <- as.character(unlist(x))
  if (grepl("^[Ff][0-9]+$", s)) return(as.integer(sub("^[Ff]", "", s)) - if (filial) 1L else 0L)
  n <- suppressWarnings(as.integer(s))
  if (is.na(n)) NULL else n
}

# Diversity panel with Hardy-Weinberg (outbred) or homozygous (inbred) lines
# drawn from Balding-Nichols subpopulations. Markers are independent (no LD).
.suite_panel_dosage <- function(lines, n_markers, n_chromosomes, chromosome_length_bp, n_subpops, fst, maf_min, inbred) {
  n <- length(lines)
  subpop <- sample(rep_len(seq_len(n_subpops), n))
  anc <- stats::runif(n_markers, 0.05, 0.5)
  a <- anc * (1 - fst) / fst; b <- (1 - anc) * (1 - fst) / fst
  sf <- matrix(vapply(seq_len(n_subpops), function(k) stats::rbeta(n_markers, a, b), numeric(n_markers)), nrow = n_markers)
  pf <- sf[, subpop, drop = FALSE]
  dos <- if (inbred) 2L * (matrix(stats::runif(n_markers * n), n_markers) < pf)
         else (matrix(stats::runif(n_markers * n), n_markers) < pf) + (matrix(stats::runif(n_markers * n), n_markers) < pf)
  storage.mode(dos) <- "integer"
  freq <- rowMeans(dos) / 2
  keep <- pmin(freq, 1 - freq) >= maf_min
  if (sum(keep) < 20L) stop("Too few common markers in the simulated panel; raise n_markers or the number of lines.", call. = FALSE)
  dos <- dos[keep, , drop = FALSE]
  chr <- sort(sample(seq_len(n_chromosomes), nrow(dos), replace = TRUE))
  pos <- unlist(lapply(seq_len(n_chromosomes), function(k) sort(sample.int(chromosome_length_bp, sum(chr == k)))))
  colnames(dos) <- lines
  variants <- data.frame(id = paste0("chr", chr, "_", pos), seqname = paste0("chr", chr), pos = pos, stringsAsFactors = FALSE)
  list(dosage = dos, variants = variants, structure = data.frame(line = lines, subpopulation = paste0("pop", subpop)))
}

.suite_write_vcf <- function(dosage, variants, path, source = "simitall") {
  gt <- c("0/0", "0/1", "1/1")[dosage + 1L]
  dim(gt) <- dim(dosage)
  contigs <- unique(variants$seqname)
  header <- c("##fileformat=VCFv4.2", paste0("##source=", source),
              paste0("##contig=<ID=", contigs, ">"),
              "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
              paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", colnames(dosage)), collapse = "\t"))
  body <- paste(variants$seqname, variants$pos, variants$id, "A", "G", ".", "PASS", ".", "GT",
                apply(gt, 1L, paste, collapse = "\t"), sep = "\t")
  writeLines(c(header, body), path)
  path
}

# Genotypes for the lines. Returns list(dosage, variants, truth_files, info).
.suite_simulate_genotypes <- function(route, lines, spec, params, outdir, seed, note) {
  gdir <- file.path(outdir, "genetic_truth")
  dir.create(gdir, showWarnings = FALSE, recursive = TRUE)
  if (route$simulator == "suite_genotype_panel") {
    a <- .suite_args("suite_genotype_panel", spec, params)$args
    inb <- .suite_is_inbred(spec, .suite_args("suite_genetics", spec, params)$args$inbred %||% "auto")
    if (inb$assumed) note(sprintf("lines treated as inbred: system '%s' is not listed in genetics.json", paste(spec$system %||% "?", collapse = ",")))
    set.seed(seed + 101L)
    g <- .suite_panel_dosage(lines, a$n_markers, a$n_chromosomes, a$chromosome_length_bp, a$n_subpops, a$fst, a$maf_min, inb$value)
    sp <- file.path(gdir, "population_structure.tsv")
    utils::write.table(g$structure, sp, sep = "\t", row.names = FALSE, quote = FALSE)
    return(list(dosage = g$dosage, variants = g$variants, files = c(population_structure = sp),
                info = list(source = "diversity_panel", inbred = inb$value, n_subpops = a$n_subpops, fst = a$fst)))
  }
  # cross: founders -> simulate_breeding()
  a <- .suite_args("suite_genotype_cross", spec, params)$args
  schemes <- .suite_tsv("mating_schemes")
  mating <- .suite_norm(spec$population$mating %||% character())
  row <- schemes[schemes$mating %in% mating, , drop = FALSE]
  if (!nrow(row)) row <- schemes[schemes$mating == "default_mapping", , drop = FALSE]
  row <- row[1, ]
  if (!length(mating)) note(sprintf("cross scheme %s (population.mating not given; mating_schemes.tsv default_mapping)", row$scheme))
  n_founders <- .suite_count(a$n_founders) %||% as.integer(row$n_founders)
  if (row$scheme %in% c("F2", "RIL", "DH", "NIL")) n_founders <- 2L
  if (row$scheme == "MAGIC" && !n_founders %in% c(4L, 8L, 16L)) stop("MAGIC needs 4, 8 or 16 founders.", call. = FALSE)
  self_g <- .suite_count(a$self_generations, filial = TRUE) %||% as.integer(row$self_generations)
  bc_g <- .suite_count(a$backcross_generations) %||% as.integer(row$backcross_generations)
  tmp <- tempfile("simitall_cross_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  panel <- file.path(tmp, "founders.fa")
  suppressMessages(utils::capture.output(generate_random_haplotype_panel(
    out_fa = panel, n_haplotypes = max(2L, n_founders), length = as.integer(a$chromosome_length_bp),
    n_chromosomes = as.integer(a$n_chromosomes), snp_rate = a$snp_rate, indel_rate = 0, seed = seed + 102L)))
  map <- file.path(tmp, "map.tsv")
  utils::write.table(data.frame(chromosome = rep(paste0("chr", seq_len(a$n_chromosomes)), each = 2),
                                pos_bp = rep(c(1L, as.integer(a$chromosome_length_bp)), a$n_chromosomes),
                                cM = rep(c(0, a$cM_per_chromosome), a$n_chromosomes)),
                     map, sep = "\t", row.names = FALSE, quote = FALSE)
  prefix <- file.path(tmp, "cross")
  founder_ids <- paste0("hap", seq_len(n_founders))
  bargs <- list(haplotype_fa = panel, out_prefix = prefix, scheme = row$scheme, n_offspring = length(lines),
                self_generations = self_g, backcross_generations = bc_g, ril_mating = row$ril_mating,
                recomb_map_in = map, vcf_out = paste0(prefix, ".vcf"), seed = seed + 103L)
  if (row$scheme %in% c("NAM", "MAGIC")) bargs$founders <- paste(founder_ids, collapse = ",") else bargs$parents <- "hap1,hap2"
  suppressMessages(utils::capture.output(do.call(simulate_breeding, bargs)))
  gt <- .simitall_read_genotypes(paste0(prefix, ".vcf"))
  dos <- gt$genotype
  dos[is.na(dos)] <- 0
  if (ncol(dos) != length(lines)) stop(sprintf("simulate_breeding() returned %d lines for %d genotype levels.", ncol(dos), length(lines)), call. = FALSE)
  freq <- rowMeans(dos) / 2
  keep <- freq > 0 & freq < 1
  if (sum(keep) < 20L) stop("Too few polymorphic markers in the cross; raise chromosome_length_bp or snp_rate.", call. = FALSE)
  old <- colnames(dos)
  dos <- dos[keep, , drop = FALSE]; storage.mode(dos) <- "integer"
  colnames(dos) <- lines
  variants <- gt$variants[keep, c("id", "seqname", "pos")]
  variants$id <- paste0(variants$seqname, "_", variants$pos)
  anc <- utils::read.delim(paste0(prefix, ".ancestry.tsv"), stringsAsFactors = FALSE)
  if ("sample" %in% names(anc)) anc$line <- lines[match(anc$sample, old)]
  ap <- file.path(gdir, "founder_ancestry.tsv")
  utils::write.table(anc, ap, sep = "\t", row.names = FALSE, quote = FALSE)
  meta <- utils::read.delim(paste0(prefix, ".meta.tsv"), stringsAsFactors = FALSE)
  if ("sample" %in% names(meta)) meta$line <- lines[match(meta$sample, old)]
  mp <- file.path(gdir, "cross_lines.tsv")
  utils::write.table(meta, mp, sep = "\t", row.names = FALSE, quote = FALSE)
  list(dosage = dos, variants = variants, files = c(founder_ancestry = ap, cross_lines = mp),
       info = list(source = "cross", scheme = row$scheme, mating = row$mating, n_founders = n_founders,
                   self_generations = self_g, backcross_generations = bc_g,
                   expected = tryCatch(breeding_expectations(row$scheme, if (row$scheme == "NIL") bc_g else 0L,
                                                             if (row$scheme %in% c("RIL", "NIL", "MAGIC", "NAM")) self_g else 0L,
                                                             n_founders, row$ril_mating), error = function(e) NULL)))
}

# QTN effects under a phenotype model. Genetic values are scaled so that the
# main genetic variance is (1 - gxe_share) * Vg and the GxE variance gxe_share * Vg,
# with Vg = h2 / (1 - h2) * residual_sd^2 (plot-level heritability).
.suite_genetic_values <- function(dosage, variants, spec, params, residual_sd, env_levels, seed, note) {
  ga <- .suite_args("suite_genetics", spec, params)$args
  models <- .suite_tsv("phenotype_models")
  mid <- ga$phenotype_model %||% "additive"
  m <- models[models$model == mid, , drop = FALSE]
  if (!nrow(m)) stop(sprintf("Unknown phenotype_model '%s'; rows of phenotype_models.tsv: %s", mid, paste(models$model, collapse = ", ")), call. = FALSE)
  pick <- function(a, type) ga[[a]] %||% .suite_parse_value(m[[a]], type)
  n_add <- pick("n_add_qtn", "integer"); n_dom <- pick("n_dom_qtn", "integer"); dd <- pick("degree_of_dom", "numeric")
  n_epi <- pick("n_epi_pairs", "integer"); gxe_share <- pick("gxe_share", "numeric")
  dist <- pick("effect_distribution", "character")
  h2 <- as.numeric(ga$heritability)
  if (!(h2 > 0 && h2 < 1)) stop("heritability must be between 0 and 1 (exclusive).", call. = FALSE)
  if (gxe_share > 0 && length(env_levels) < 2) {
    note("phenotype model asks for GxE but the design has no second treatment factor; GxE share set to 0")
    gxe_share <- 0
  }
  if (n_dom > 0 && !any(dosage == 1L)) note("dominance QTN have no effect: the lines carry no heterozygous markers")
  set.seed(seed + 201L)
  freq <- rowMeans(dosage) / 2
  maf <- pmin(freq, 1 - freq)
  elig <- which(maf >= (ga$qtn_maf_min %||% 0.1))
  need <- n_add + n_dom + 2L * n_epi
  if (length(elig) < max(1L, n_add)) {
    # backcross / NIL lines: the donor allele is rare everywhere; take the commonest polymorphic markers instead
    poly <- which(maf > 0)
    if (length(poly) < max(1L, n_add)) stop(sprintf("Only %d polymorphic markers for %d QTN.", length(poly), n_add), call. = FALSE)
    elig <- poly[order(-maf[poly])][seq_len(min(length(poly), max(need, 5L * n_add)))]
    note(sprintf("only %d markers have MAF >= %s (backcross or NIL lines?); QTN drawn among the %d commonest polymorphic markers (MAF %.3f-%.3f)",
                 sum(maf >= (ga$qtn_maf_min %||% 0.1)), ga$qtn_maf_min %||% 0.1, length(elig), min(maf[elig]), max(maf[elig])))
  }
  draw <- function(k) {
    if (k < 1) return(numeric())
    e <- switch(dist, geometric = 0.9^(seq_len(k) - 1), equal = rep(1, k), normal = abs(stats::rnorm(k)),
                stop("effect_distribution must be geometric, normal or equal", call. = FALSE))
    e * sample(c(-1, 1), k, replace = TRUE)
  }
  qtn <- sample(elig, min(length(elig), need))
  take <- function(k) { if (k < 1) return(integer()); out <- qtn[seq_len(min(k, length(qtn)))]; qtn <<- qtn[-seq_len(min(k, length(qtn)))]; out }
  add_i <- take(n_add); dom_i <- take(n_dom); epi_i <- take(2L * n_epi)
  X <- dosage - 1  # -1/0/1 coding
  rows <- list()
  g_add <- numeric(ncol(dosage)); a_eff <- draw(length(add_i))
  if (length(add_i)) g_add <- as.numeric(crossprod(X[add_i, , drop = FALSE], a_eff))
  g_dom <- numeric(ncol(dosage)); ad <- d_eff <- numeric()
  if (length(dom_i)) {
    ad <- draw(length(dom_i)); d_eff <- dd * abs(ad)
    g_dom <- as.numeric(crossprod(X[dom_i, , drop = FALSE], ad)) + as.numeric(crossprod((dosage[dom_i, , drop = FALSE] == 1L) * 1, d_eff))
  }
  g_epi <- numeric(ncol(dosage))
  epi_pairs <- if (length(epi_i) >= 2) matrix(epi_i[seq_len(2 * (length(epi_i) %/% 2))], ncol = 2, byrow = TRUE) else matrix(integer(), ncol = 2)
  e_eff <- draw(nrow(epi_pairs))
  for (k in seq_len(nrow(epi_pairs))) g_epi <- g_epi + e_eff[k] * X[epi_pairs[k, 1], ] * X[epi_pairs[k, 2], ]
  g <- g_add + g_dom + g_epi
  vg <- h2 / (1 - h2) * residual_sd^2
  if (stats::var(g) <= 0) stop("The QTN produced no genetic variance among the lines.", call. = FALSE)
  s_main <- sqrt((1 - gxe_share) * vg / stats::var(g))
  g <- (g - mean(g)) * s_main
  # GxE: the additive QTN effects change linearly across the other factor's levels
  gxe <- NULL
  if (gxe_share > 0) {
    z <- seq(-1, 1, length.out = length(env_levels))
    b <- draw(length(add_i))
    gx <- as.numeric(crossprod(X[add_i, , drop = FALSE], b))
    gx <- gx - mean(gx)
    vgx <- stats::var(as.numeric(outer(gx, z)))
    s_gx <- if (vgx > 0) sqrt(gxe_share * vg / vgx) else 0
    gxe <- outer(gx * s_gx, z)
    dimnames(gxe) <- list(colnames(dosage), env_levels)
  }
  qt <- function(idx, type, eff, dom = NA_real_, partner = NA_character_) if (length(idx)) data.frame(
    marker_id = variants$id[idx], seqname = variants$seqname[idx], pos = variants$pos[idx], type = type,
    effect = eff, dominance = dom, maf = maf[idx], partner = partner, stringsAsFactors = FALSE)
  qtn_tab <- rbind(
    qt(add_i, "additive", a_eff * s_main),
    if (length(dom_i)) qt(dom_i, "dominance", ad * s_main, d_eff * s_main),
    if (nrow(epi_pairs)) qt(epi_pairs[, 1], "epistatic", e_eff * s_main, partner = variants$id[epi_pairs[, 2]]),
    if (!is.null(gxe)) qt(add_i, "gxe_slope", b * s_gx)
  )
  names(g) <- colnames(dosage)
  list(g = g, gxe = gxe, qtn = qtn_tab,
       info = list(phenotype_model = mid, heritability = h2, genetic_variance = vg,
                   realized_main_variance = stats::var(g),
                   realized_gxe_variance = if (is.null(gxe)) 0 else stats::var(as.numeric(gxe)),
                   n_add_qtn = length(add_i), n_dom_qtn = length(dom_i), degree_of_dom = dd,
                   n_epi_pairs = nrow(epi_pairs), gxe_share = gxe_share, effect_distribution = dist,
                   scale = "additive coding -1/0/1; heterozygote adds degree_of_dom x |a| for dominance QTN"))
}
