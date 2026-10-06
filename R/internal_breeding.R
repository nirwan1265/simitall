main_11_simulate_breeding <- function(args = commandArgs(trailingOnly = TRUE)) {
  args <- args

  usage <- function() {
    cat(paste0(
      "Usage: simulate_breeding --haplotype_fa <path> --out_prefix <path> [options]\n\n",
      "Multi-chromosome FASTA headers use >founder|chromosome. Legacy\n",
      "single-record founder panels are treated as chromosome chr1.\n\n",
      "Core options:\n",
      "  --parents <csv>                 two founder IDs or indices\n",
      "  --founders <csv>                founder IDs/indices for NAM/MAGIC\n",
      "  --n_founders <int>              founders sampled when omitted (default 4)\n",
      "  --n_offspring <int>             final population size (default 100)\n",
      "  --scheme <F2|MAGIC|NAM|RIL|NIL|DH>\n",
      "  --sequence <tokens>             e.g. F1,SELF:3,SIB:2,BC:P1:2,DH\n",
      "  --self_generations <int>        RIL/NIL/MAGIC selfing cycles (default 6)\n",
      "  --backcross_generations <int>   NIL backcross cycles (default 3)\n",
      "  --ril_mating <SSD|SIB>          default SSD\n\n",
      "Recombination:\n",
      "  --recomb_map_in <path>          TSV with chromosome,pos_bp,cM\n",
      "  --recomb_rate_mean <float>      cM/Mb mean (default 1.0)\n",
      "  --recomb_rate_sd <float>        cM/Mb SD (default 0.3)\n",
      "  --recomb_hotspots <int>         hotspots per chromosome (default 3)\n",
      "  --recomb_hotspot_mult <float>   hotspot multiplier (default 5.0)\n",
      "  --interference_shape <float>    gamma interference shape (default 1)\n\n",
      "Truth outputs:\n",
      "  --vcf_out <path>                chromosome-aware VCF\n",
      "  --breakpoints_out <path>        founder-ancestry breakpoints TSV\n",
      "  --ancestry_out <path>           founder-ancestry tracts TSV\n",
      "  --recomb_map_out <path>         realized maps TSV\n",
      "  --graph_out <path>              Mermaid crossing graph\n"
    ))
    stop("Missing required breeding arguments", call. = FALSE)
  }

  get_arg <- function(flag, default = NULL) {
    if (!(flag %in% args)) return(default)
    idx <- match(flag, args)
    if (idx == length(args)) return(default)
    args[idx + 1L]
  }
  has_flag <- function(flag) flag %in% args

  hap_fa <- get_arg("--haplotype_fa")
  out_prefix <- get_arg("--out_prefix")
  if (is.null(hap_fa) || is.null(out_prefix)) usage()

  parents_arg <- get_arg("--parents", NA_character_)
  founders_arg <- get_arg("--founders", NA_character_)
  n_founders <- as.integer(get_arg("--n_founders", 4))
  n_offspring <- as.integer(get_arg("--n_offspring", 100))
  sequence_spec <- get_arg("--sequence", "F1,SELF:1")
  scheme <- toupper(get_arg("--scheme", NA_character_))
  self_generations <- as.integer(
    get_arg("--self_generations", get_arg("--generations", 6))
  )
  backcross_generations <- as.integer(get_arg("--backcross_generations", 3))
  ril_mating <- toupper(get_arg("--ril_mating", "SSD"))

  vcf_out <- get_arg("--vcf_out", NA_character_)
  graph_out <- get_arg("--graph_out", NA_character_)
  graph_format <- tolower(get_arg("--graph_format", "mmd"))
  breakpoints_out <- get_arg(
    "--breakpoints_out", paste0(out_prefix, ".breakpoints.tsv")
  )
  ancestry_out <- get_arg("--ancestry_out", paste0(out_prefix, ".ancestry.tsv"))
  recomb_map_out <- get_arg(
    "--recomb_map_out", paste0(out_prefix, ".recombination_map.tsv")
  )
  sv_truth_out <- get_arg("--sv_truth_out", paste0(out_prefix, ".sv_truth.tsv"))

  recomb_map_in <- get_arg("--recomb_map_in", NA_character_)
  recomb_rate_mean <- as.numeric(get_arg("--recomb_rate_mean", 1.0))
  recomb_rate_sd <- as.numeric(get_arg("--recomb_rate_sd", 0.3))
  recomb_hotspots <- as.integer(get_arg("--recomb_hotspots", 3))
  recomb_hotspot_mult <- as.numeric(get_arg("--recomb_hotspot_mult", 5.0))
  interference_shape <- as.numeric(get_arg("--interference_shape", 1.0))

  fix_locus <- get_arg("--fix_locus", NA_character_)
  fix_allele <- get_arg("--fix_allele", NA_character_)
  background_selection <- has_flag("--background_selection")
  selection_pool <- as.integer(get_arg("--selection_pool", 50))
  marker_step <- as.integer(get_arg("--marker_step", 1000))
  introgression_target_len <- as.numeric(
    get_arg("--introgression_target_len", NA_character_)
  )
  genotype_error <- as.numeric(get_arg("--genotype_error", 0))
  missing_rate <- as.numeric(get_arg("--missing_rate", 0))
  ascertainment <- tolower(get_arg("--ascertainment", "founders"))
  sv_rate <- as.numeric(get_arg("--sv_rate", 0))
  sv_maxlen <- as.integer(get_arg("--sv_maxlen", 1000))
  selection_loci_arg <- get_arg("--selection_loci", NA_character_)
  selection_model <- tolower(get_arg("--selection_model", "add"))
  selection_strength <- as.numeric(get_arg("--selection_strength", 0))
  distortion_rate <- as.numeric(get_arg("--distortion_rate", 0))
  seed <- as.integer(get_arg("--seed", 1))

  if (n_offspring < 1L) stop("n_offspring must be positive")
  if (self_generations < 0L || backcross_generations < 0L) {
    stop("generation counts cannot be negative")
  }
  if (!ril_mating %in% c("SSD", "SIB")) {
    stop("ril_mating must be SSD or SIB")
  }
  set.seed(seed)

  parse_header <- function(header) {
    token <- strsplit(header, "[[:space:]]+")[[1L]][1L]
    if (grepl("::", token, fixed = TRUE)) {
      fields <- strsplit(token, "::", fixed = TRUE)[[1L]]
      return(c(founder = fields[1L], chromosome = paste(fields[-1L], collapse = "::")))
    }
    if (grepl("|", token, fixed = TRUE)) {
      fields <- strsplit(token, "|", fixed = TRUE)[[1L]]
      return(c(founder = fields[1L], chromosome = paste(fields[-1L], collapse = "|")))
    }
    founder_match <- regmatches(header, regexec("founder=([^ ;]+)", header))[[1L]]
    chrom_match <- regmatches(
      header, regexec("(?:chromosome|chrom|chr)=([^ ;]+)", header)
    )[[1L]]
    if (length(founder_match) > 1L && length(chrom_match) > 1L) {
      return(c(founder = founder_match[2L], chromosome = chrom_match[2L]))
    }
    c(founder = token, chromosome = "chr1")
  }

  read_panel <- function(path) {
    if (!file.exists(path)) stop("Haplotype FASTA not found: ", path)
    lines <- readLines(path, warn = FALSE)
    header_idx <- grep("^>", lines)
    if (!length(header_idx)) stop("No FASTA records found in: ", path)
    headers <- substring(lines[header_idx], 2L)
    parsed <- t(vapply(headers, parse_header, character(2L)))
    seqs <- vapply(seq_along(header_idx), function(i) {
      start <- header_idx[i] + 1L
      end <- if (i < length(header_idx)) header_idx[i + 1L] - 1L else length(lines)
      if (start > end) return("")
      toupper(paste(lines[start:end], collapse = ""))
    }, character(1L))
    records <- data.frame(
      founder = parsed[, "founder"], chromosome = parsed[, "chromosome"],
      sequence = seqs, stringsAsFactors = FALSE
    )
    if (any(!nzchar(records$founder)) || any(!nzchar(records$chromosome))) {
      stop("FASTA founder and chromosome labels cannot be empty")
    }
    if (anyDuplicated(records[c("founder", "chromosome")])) {
      stop("Each founder/chromosome pair must occur exactly once")
    }
    founders <- unique(records$founder)
    if (length(founders) < 2L) stop("Need at least two founders")
    chromosomes <- records$chromosome[records$founder == founders[1L]]
    if (!length(chromosomes)) stop("No chromosomes found")

    genomes <- stats::setNames(vector("list", length(founders)), founders)
    chromosome_lengths <- stats::setNames(integer(length(chromosomes)), chromosomes)
    for (founder in founders) {
      block <- records[records$founder == founder, , drop = FALSE]
      if (!setequal(block$chromosome, chromosomes)) {
        stop("Every founder must contain the same chromosome set; problem: ", founder)
      }
      block <- block[match(chromosomes, block$chromosome), , drop = FALSE]
      genomes[[founder]] <- stats::setNames(block$sequence, chromosomes)
      lengths <- nchar(block$sequence)
      if (founder == founders[1L]) {
        chromosome_lengths[] <- lengths
      } else if (any(lengths != chromosome_lengths)) {
        stop("Aligned sequence lengths differ for founder: ", founder)
      }
    }
    if (any(chromosome_lengths < 2L)) stop("Chromosomes must contain at least 2 bp")
    list(
      founder_ids = founders, chromosomes = chromosomes,
      chromosome_lengths = chromosome_lengths, genomes = genomes
    )
  }

  panel <- read_panel(hap_fa)
  founder_ids_all <- panel$founder_ids
  chromosomes <- panel$chromosomes
  chromosome_lengths <- panel$chromosome_lengths
  genomes <- panel$genomes

  resolve_ids <- function(value, count = NULL) {
    if (!is.na(value) && nzchar(value)) {
      fields <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
      idx <- vapply(fields, function(field) {
        if (grepl("^[0-9]+$", field)) as.integer(field) else match(field, founder_ids_all)
      }, integer(1L))
      if (anyNA(idx) || any(idx < 1L) || any(idx > length(founder_ids_all))) {
        stop("Founder IDs or indices were not found: ", value)
      }
      return(idx)
    }
    if (is.null(count)) return(integer())
    if (count > length(founder_ids_all)) stop("n_founders exceeds available founders")
    sample(seq_along(founder_ids_all), count, replace = FALSE)
  }

  multi_parent_scheme <- !is.na(scheme) && scheme %in% c("NAM", "MAGIC")
  founder_idx <- if (multi_parent_scheme || !is.na(founders_arg)) {
    resolve_ids(founders_arg, n_founders)
  } else {
    integer()
  }
  if (multi_parent_scheme && length(founder_idx) < 2L) {
    stop(scheme, " requires at least two founders")
  }
  parent_idx <- if (!is.na(parents_arg)) {
    resolved <- resolve_ids(parents_arg)
    if (length(resolved) != 2L) stop("parents must contain exactly two values")
    resolved
  } else if (multi_parent_scheme) {
    founder_idx[1:2]
  } else {
    sample(seq_along(founder_ids_all), 2L, replace = FALSE)
  }
  if (!length(founder_idx)) founder_idx <- parent_idx

  normalize_map <- function(map, chromosome, chromosome_length) {
    map <- map[order(map$pos_bp), c("pos_bp", "cM"), drop = FALSE]
    map$pos_bp <- as.integer(map$pos_bp)
    map$cM <- as.numeric(map$cM)
    map <- map[is.finite(map$pos_bp) & is.finite(map$cM), , drop = FALSE]
    map$pos_bp <- pmax(1L, pmin(chromosome_length, map$pos_bp))
    map <- map[!duplicated(map$pos_bp), , drop = FALSE]
    if (!nrow(map)) stop("Empty recombination map for chromosome ", chromosome)
    map$cM <- map$cM - map$cM[1L]
    if (map$pos_bp[1L] > 1L) {
      map <- rbind(data.frame(pos_bp = 1L, cM = 0), map)
    }
    if (tail(map$pos_bp, 1L) < chromosome_length) {
      final_cm <- tail(map$cM, 1L)
      map <- rbind(map, data.frame(pos_bp = chromosome_length, cM = final_cm))
    }
    if (is.unsorted(map$cM, strictly = FALSE)) {
      stop("cM positions must be nondecreasing for chromosome ", chromosome)
    }
    map$chromosome <- chromosome
    map[, c("chromosome", "pos_bp", "cM")]
  }

  random_map <- function(chromosome, chromosome_length) {
    n_intervals <- max(2L, min(1000L, ceiling(chromosome_length / 1000)))
    positions <- unique(as.integer(round(seq(
      1L, chromosome_length, length.out = n_intervals
    ))))
    rates <- pmax(
      0, stats::rnorm(length(positions) - 1L, recomb_rate_mean, recomb_rate_sd)
    )
    if (recomb_hotspots > 0L && length(rates)) {
      hot <- sample(seq_along(rates), min(recomb_hotspots, length(rates)))
      rates[hot] <- rates[hot] * recomb_hotspot_mult
    }
    data.frame(
      chromosome = chromosome,
      pos_bp = positions,
      cM = c(0, cumsum(rates * diff(positions) / 1e6)),
      stringsAsFactors = FALSE
    )
  }

  if (!is.na(recomb_map_in)) {
    if (!file.exists(recomb_map_in)) stop("Recombination map not found: ", recomb_map_in)
    raw_map <- utils::read.delim(
      recomb_map_in, stringsAsFactors = FALSE, check.names = FALSE
    )
    chrom_col <- intersect(c("chromosome", "chrom", "chr", "CHROM"), names(raw_map))
    if (!all(c("pos_bp", "cM") %in% names(raw_map))) {
      stop("Recombination map requires pos_bp and cM columns")
    }
    if (!length(chrom_col)) {
      if (length(chromosomes) > 1L) {
        stop("A multi-chromosome recombination map requires a chromosome column")
      }
      raw_map$chromosome <- chromosomes[1L]
    } else {
      names(raw_map)[match(chrom_col[1L], names(raw_map))] <- "chromosome"
    }
    missing_maps <- setdiff(chromosomes, unique(raw_map$chromosome))
    if (length(missing_maps)) {
      stop("Recombination map is missing: ", paste(missing_maps, collapse = ", "))
    }
    recomb_maps <- stats::setNames(lapply(chromosomes, function(chromosome) {
      normalize_map(
        raw_map[raw_map$chromosome == chromosome, , drop = FALSE], chromosome,
        chromosome_lengths[[chromosome]]
      )
    }), chromosomes)
  } else {
    recomb_maps <- stats::setNames(lapply(chromosomes, function(chromosome) {
      random_map(chromosome, chromosome_lengths[[chromosome]])
    }), chromosomes)
  }
  recomb_map_table <- do.call(rbind, recomb_maps)
  rownames(recomb_map_table) <- NULL
  dir.create(dirname(recomb_map_out), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    recomb_map_table, recomb_map_out, sep = "\t", quote = FALSE,
    row.names = FALSE
  )

  founder_ancestry <- function(founder_id) {
    stats::setNames(lapply(chromosomes, function(chromosome) {
      data.frame(
        start = 1L, end = chromosome_lengths[[chromosome]],
        founder = founder_id, stringsAsFactors = FALSE
      )
    }), chromosomes)
  }

  founder_individual <- function(index, family = "") {
    founder_id <- founder_ids_all[index]
    genome <- genomes[[founder_id]]
    ancestry <- founder_ancestry(founder_id)
    list(
      h1 = genome, h2 = genome, a1 = ancestry, a2 = ancestry,
      generation = "P", scheme = "founder", family = family
    )
  }

  merge_tracts <- function(tracts) {
    if (!nrow(tracts)) return(tracts)
    tracts <- tracts[order(tracts$start), , drop = FALSE]
    out <- tracts[1L, , drop = FALSE]
    if (nrow(tracts) > 1L) {
      for (i in 2:nrow(tracts)) {
        last <- nrow(out)
        if (out$founder[last] == tracts$founder[i] &&
            out$end[last] + 1L >= tracts$start[i]) {
          out$end[last] <- max(out$end[last], tracts$end[i])
        } else {
          out <- rbind(out, tracts[i, , drop = FALSE])
        }
      }
    }
    rownames(out) <- NULL
    out
  }

  slice_tracts <- function(tracts, start, end) {
    keep <- tracts$end >= start & tracts$start <= end
    out <- tracts[keep, , drop = FALSE]
    if (!nrow(out)) return(out)
    out$start <- pmax(out$start, start)
    out$end <- pmin(out$end, end)
    out
  }

  replace_ancestry <- function(tracts, start, end, founder) {
    before <- if (start > 1L) slice_tracts(tracts, 1L, start - 1L) else tracts[0, ]
    middle <- data.frame(start = start, end = end, founder = founder)
    after <- if (end < max(tracts$end)) {
      slice_tracts(tracts, end + 1L, max(tracts$end))
    } else tracts[0, ]
    merge_tracts(rbind(before, middle, after))
  }

  draw_crossovers <- function(map, chromosome_length) {
    total_cm <- tail(map$cM, 1L)
    if (!is.finite(total_cm) || total_cm <= 0) return(integer())
    expected <- total_cm / 100
    if (interference_shape <= 1) {
      n_xo <- stats::rpois(1L, expected)
      if (!n_xo) return(integer())
      xo_cm <- sort(stats::runif(n_xo, 0, total_cm))
    } else {
      mean_distance <- total_cm / (expected + 1)
      scale <- mean_distance / interference_shape
      xo_cm <- numeric()
      position <- 0
      while (position < total_cm) {
        position <- position + stats::rgamma(
          1L, shape = interference_shape, scale = scale
        )
        if (position < total_cm) xo_cm <- c(xo_cm, position)
      }
      if (!length(xo_cm)) return(integer())
    }
    positions <- as.integer(round(stats::approx(
      x = map$cM, y = map$pos_bp, xout = xo_cm,
      method = "linear", ties = "ordered", rule = 2
    )$y))
    sort(unique(pmax(1L, pmin(chromosome_length - 1L, positions))))
  }

  recombine_chromosome <- function(seq1, seq2, anc1, anc2, crossovers) {
    chromosome_length <- nchar(seq1)
    starts <- c(1L, crossovers + 1L)
    ends <- c(crossovers, chromosome_length)
    use_first <- stats::runif(1L) < 0.5
    pieces <- character(length(starts))
    ancestry <- anc1[0, ]
    for (i in seq_along(starts)) {
      source_seq <- if (use_first) seq1 else seq2
      source_anc <- if (use_first) anc1 else anc2
      pieces[i] <- substr(source_seq, starts[i], ends[i])
      ancestry <- rbind(
        ancestry, slice_tracts(source_anc, starts[i], ends[i])
      )
      use_first <- !use_first
    }
    list(sequence = paste0(pieces, collapse = ""), ancestry = merge_tracts(ancestry))
  }

  make_gamete <- function(parent) {
    sequences <- stats::setNames(character(length(chromosomes)), chromosomes)
    ancestry <- stats::setNames(vector("list", length(chromosomes)), chromosomes)
    crossovers <- stats::setNames(vector("list", length(chromosomes)), chromosomes)
    for (chromosome in chromosomes) {
      xo <- draw_crossovers(
        recomb_maps[[chromosome]], chromosome_lengths[[chromosome]]
      )
      recombinant <- recombine_chromosome(
        parent$h1[[chromosome]], parent$h2[[chromosome]],
        parent$a1[[chromosome]], parent$a2[[chromosome]], xo
      )
      sequences[[chromosome]] <- recombinant$sequence
      ancestry[[chromosome]] <- recombinant$ancestry
      crossovers[[chromosome]] <- xo
    }
    list(sequence = sequences, ancestry = ancestry, crossovers = crossovers)
  }

  parse_locus <- function(value) {
    if (is.na(value) || !nzchar(value)) return(NULL)
    normalized <- gsub("-", ":", value, fixed = TRUE)
    fields <- strsplit(normalized, ":", fixed = TRUE)[[1L]]
    if (length(fields) == 2L) {
      chromosome <- chromosomes[1L]
      coordinates <- suppressWarnings(as.integer(fields))
    } else if (length(fields) == 3L) {
      chromosome <- fields[1L]
      coordinates <- suppressWarnings(as.integer(fields[2:3]))
    } else {
      stop("Locus must use start:end or chromosome:start:end")
    }
    if (!chromosome %in% chromosomes || anyNA(coordinates) ||
        coordinates[1L] > coordinates[2L]) {
      stop("Invalid locus: ", value)
    }
    coordinates[1L] <- max(1L, coordinates[1L])
    coordinates[2L] <- min(chromosome_lengths[[chromosome]], coordinates[2L])
    list(chromosome = chromosome, start = coordinates[1L], end = coordinates[2L])
  }

  parse_positions <- function(value) {
    if (is.na(value) || !nzchar(value)) {
      return(data.frame(chromosome = character(), position = integer()))
    }
    fields <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
    rows <- lapply(fields, function(field) {
      parts <- strsplit(field, ":", fixed = TRUE)[[1L]]
      if (length(parts) == 1L) {
        chromosome <- chromosomes[1L]
        position <- as.integer(parts[1L])
      } else {
        chromosome <- parts[1L]
        position <- as.integer(parts[2L])
      }
      if (!chromosome %in% chromosomes || is.na(position) || position < 1L ||
          position > chromosome_lengths[[chromosome]]) {
        stop("Invalid selected locus: ", field)
      }
      data.frame(chromosome = chromosome, position = position)
    })
    do.call(rbind, rows)
  }

  fixed_locus <- parse_locus(fix_locus)
  selected_loci <- parse_positions(selection_loci_arg)
  P1 <- founder_individual(parent_idx[1L])
  P2 <- founder_individual(parent_idx[2L])
  fixed_founder <- NULL
  if (!is.na(fix_allele)) {
    if (tolower(fix_allele) == "donor") {
      fixed_founder <- founder_ids_all[parent_idx[2L]]
    } else if (tolower(fix_allele) == "recipient") {
      fixed_founder <- founder_ids_all[parent_idx[1L]]
    } else if (fix_allele %in% founder_ids_all) {
      fixed_founder <- fix_allele
    } else {
      stop("fix_allele must be donor, recipient, or a founder ID")
    }
  }

  apply_fixed_locus <- function(gamete1, gamete2) {
    if (is.null(fixed_locus) || is.null(fixed_founder)) {
      return(list(g1 = gamete1, g2 = gamete2))
    }
    chromosome <- fixed_locus$chromosome
    start <- fixed_locus$start
    end <- fixed_locus$end
    donor_sequence <- genomes[[fixed_founder]][[chromosome]]
    replacement <- substr(donor_sequence, start, end)
    replace_sequence <- function(sequence) {
      paste0(
        substr(sequence, 1L, start - 1L), replacement,
        substr(sequence, end + 1L, nchar(sequence))
      )
    }
    gamete1$sequence[[chromosome]] <- replace_sequence(
      gamete1$sequence[[chromosome]]
    )
    gamete2$sequence[[chromosome]] <- replace_sequence(
      gamete2$sequence[[chromosome]]
    )
    gamete1$ancestry[[chromosome]] <- replace_ancestry(
      gamete1$ancestry[[chromosome]], start, end, fixed_founder
    )
    gamete2$ancestry[[chromosome]] <- replace_ancestry(
      gamete2$ancestry[[chromosome]], start, end, fixed_founder
    )
    list(g1 = gamete1, g2 = gamete2)
  }

  donor_allele_count <- function(gamete1, gamete2) {
    if (!nrow(selected_loci)) return(0L)
    donor <- genomes[[founder_ids_all[parent_idx[2L]]]]
    counts <- vapply(seq_len(nrow(selected_loci)), function(i) {
      chromosome <- selected_loci$chromosome[i]
      position <- selected_loci$position[i]
      allele <- substr(donor[[chromosome]], position, position)
      as.integer(substr(gamete1$sequence[[chromosome]], position, position) == allele) +
        as.integer(substr(gamete2$sequence[[chromosome]], position, position) == allele)
    }, integer(1L))
    if (selection_model == "dom") return(sum(counts > 0L))
    if (selection_model == "rec") return(sum(counts == 2L))
    sum(counts)
  }

  apply_distortion <- function(gamete1, gamete2) {
    if (!nrow(selected_loci) || distortion_rate <= 0) {
      return(list(g1 = gamete1, g2 = gamete2))
    }
    donor_id <- founder_ids_all[parent_idx[2L]]
    donor <- genomes[[donor_id]]
    for (i in seq_len(nrow(selected_loci))) {
      if (stats::runif(1L) >= distortion_rate) next
      chromosome <- selected_loci$chromosome[i]
      position <- selected_loci$position[i]
      allele <- substr(donor[[chromosome]], position, position)
      replace_base <- function(sequence) {
        paste0(
          substr(sequence, 1L, position - 1L), allele,
          substr(sequence, position + 1L, nchar(sequence))
        )
      }
      gamete1$sequence[[chromosome]] <- replace_base(
        gamete1$sequence[[chromosome]]
      )
      gamete2$sequence[[chromosome]] <- replace_base(
        gamete2$sequence[[chromosome]]
      )
      gamete1$ancestry[[chromosome]] <- replace_ancestry(
        gamete1$ancestry[[chromosome]], position, position, donor_id
      )
      gamete2$ancestry[[chromosome]] <- replace_ancestry(
        gamete2$ancestry[[chromosome]], position, position, donor_id
      )
    }
    list(g1 = gamete1, g2 = gamete2)
  }

  make_child <- function(parent1, parent2, generation, child_scheme, family = "") {
    attempts <- 0L
    repeat {
      attempts <- attempts + 1L
      g1 <- make_gamete(parent1)
      g2 <- make_gamete(parent2)
      distorted <- apply_distortion(g1, g2)
      fixed <- apply_fixed_locus(distorted$g1, distorted$g2)
      risk <- donor_allele_count(fixed$g1, fixed$g2)
      fitness <- if (selection_strength > 0) exp(-selection_strength * risk) else 1
      if (stats::runif(1L) <= fitness || attempts >= 10000L) break
    }
    list(
      h1 = fixed$g1$sequence, h2 = fixed$g2$sequence,
      a1 = fixed$g1$ancestry, a2 = fixed$g2$ancestry,
      generation = generation, scheme = child_scheme, family = family
    )
  }

  make_children <- function(
      n, parent1, parent2, generation, child_scheme, family = "") {
    lapply(seq_len(n), function(i) {
      make_child(parent1, parent2, generation, child_scheme, family)
    })
  }

  marker_positions <- stats::setNames(lapply(chromosomes, function(chromosome) {
    seq.int(1L, chromosome_lengths[[chromosome]], by = max(1L, marker_step))
  }), chromosomes)
  donor_id <- founder_ids_all[parent_idx[2L]]

  donor_fraction <- function(individual, exclude_locus = NULL) {
    matches <- 0L
    total <- 0L
    for (chromosome in chromosomes) {
      positions <- marker_positions[[chromosome]]
      if (!is.null(exclude_locus) && chromosome == exclude_locus$chromosome) {
        positions <- positions[
          positions < exclude_locus$start | positions > exclude_locus$end
        ]
      }
      if (!length(positions)) next
      donor_sequence <- genomes[[donor_id]][[chromosome]]
      for (position in positions) {
        donor_base <- substr(donor_sequence, position, position)
        matches <- matches +
          as.integer(substr(individual$h1[[chromosome]], position, position) == donor_base) +
          as.integer(substr(individual$h2[[chromosome]], position, position) == donor_base)
        total <- total + 2L
      }
    }
    if (!total) return(0)
    matches / total
  }

  introgression_length <- function(individual, locus) {
    if (is.null(locus)) return(NA_real_)
    chromosome <- locus$chromosome
    lengths <- vapply(list(individual$a1, individual$a2), function(ancestry) {
      tracts <- ancestry[[chromosome]]
      hit <- tracts$start <= locus$start & tracts$end >= locus$end &
        tracts$founder == fixed_founder
      if (!any(hit)) return(NA_real_)
      tract <- tracts[which(hit)[1L], , drop = FALSE]
      tract$end - tract$start + 1L
    }, numeric(1L))
    mean(lengths, na.rm = TRUE)
  }

  parse_token <- function(token) {
    parts <- strsplit(token, ":", fixed = TRUE)[[1L]]
    list(
      type = toupper(parts[1L]),
      arg1 = if (length(parts) >= 2L) parts[2L] else NA_character_,
      arg2 = if (length(parts) >= 3L) parts[3L] else NA_character_
    )
  }
  expand_tokens <- function(tokens) {
    expanded <- character()
    for (token in tokens) {
      parsed <- parse_token(token)
      if (parsed$type %in% c("SELF", "SIB") && !is.na(parsed$arg1)) {
        expanded <- c(expanded, rep(parsed$type, as.integer(parsed$arg1)))
      } else if (parsed$type == "BC" && !is.na(parsed$arg2)) {
        expanded <- c(
          expanded, rep(paste("BC", parsed$arg1, sep = ":"),
                        as.integer(parsed$arg2))
        )
      } else {
        expanded <- c(expanded, token)
      }
    }
    expanded
  }

  sequence_tokens <- trimws(strsplit(sequence_spec, ",", fixed = TRUE)[[1L]])
  current_population <- NULL

  if (!is.na(scheme) && scheme == "MAGIC") {
    # Each final line comes from its own funnel (independent intercrosses), as
    # in real MAGIC populations. A single shared funnel would make every line a
    # selfed descendant of one plant, so founder contributions would not
    # average 1/n_founders across lines.
    run_funnel <- function() {
      parents <- lapply(founder_idx, founder_individual)
      round <- 1L
      while (length(parents) > 1L) {
        next_round <- list()
        pair_index <- 1L
        i <- 1L
        while (i <= length(parents)) {
          if (i == length(parents)) {
            next_round[[length(next_round) + 1L]] <- parents[[i]]
          } else {
            next_round[[length(next_round) + 1L]] <- make_child(
              parents[[i]], parents[[i + 1L]], paste0("MAGIC_R", round),
              "MAGIC", paste0("pair", pair_index)
            )
            pair_index <- pair_index + 1L
          }
          i <- i + 2L
        }
        parents <- next_round
        round <- round + 1L
      }
      parents[[1L]]
    }
    current_population <- lapply(seq_len(n_offspring), function(i) run_funnel())
    sequence_tokens <- rep("SELF", self_generations)
  } else if (!is.na(scheme) && scheme == "NAM") {
    common <- founder_individual(founder_idx[1L])
    family_parents <- lapply(seq_along(founder_idx[-1L]), function(i) {
      make_child(
        common, founder_individual(founder_idx[i + 1L]), "NAM_F1", "NAM",
        paste0("NAM_F", i)
      )
    })
    family_sizes <- rep(n_offspring %/% length(family_parents), length(family_parents))
    family_sizes[seq_len(n_offspring %% length(family_parents))] <-
      family_sizes[seq_len(n_offspring %% length(family_parents))] + 1L
    current_population <- unlist(lapply(seq_along(family_parents), function(i) {
      make_children(
        family_sizes[i], family_parents[[i]], family_parents[[i]],
        "NAM_F2", "NAM", paste0("NAM_F", i)
      )
    }), recursive = FALSE)
    sequence_tokens <- character()
  } else if (!is.na(scheme) && scheme == "F2") {
    sequence_tokens <- c("F1", "SELF")
  } else if (!is.na(scheme) && scheme == "RIL") {
    sequence_tokens <- c(
      "F1", rep(if (ril_mating == "SIB") "SIB" else "SELF", self_generations)
    )
  } else if (!is.na(scheme) && scheme == "NIL") {
    sequence_tokens <- c(
      "F1", rep("BC:P1", backcross_generations), rep("SELF", self_generations)
    )
  } else if (!is.na(scheme) && scheme == "DH") {
    sequence_tokens <- c("F1", "DH")
  }

  sequence_tokens <- expand_tokens(sequence_tokens)
  filial_generation <- 0L
  backcross_generation <- 0L
  sib_generation <- 0L

  for (token in sequence_tokens) {
    parsed <- parse_token(token)
    type <- parsed$type
    if (type == "F1") {
      count <- if (!is.na(parsed$arg1)) as.integer(parsed$arg1) else n_offspring
      filial_generation <- 1L
      current_population <- make_children(count, P1, P2, "F1", "F1")
    } else if (type == "SELF") {
      if (is.null(current_population)) stop("SELF requires an existing population")
      filial_generation <- filial_generation + 1L
      generation <- if (filial_generation > 0L) {
        paste0("F", filial_generation)
      } else "SELF"
      current_population <- lapply(seq_len(n_offspring), function(i) {
        parent <- current_population[[((i - 1L) %% length(current_population)) + 1L]]
        make_child(parent, parent, generation, "SELF", parent$family)
      })
    } else if (type == "SIB") {
      if (is.null(current_population)) stop("SIB requires an existing population")
      sib_generation <- sib_generation + 1L
      generation <- paste0("SIB", sib_generation)
      current_population <- lapply(seq_len(n_offspring), function(i) {
        chosen <- sample(seq_along(current_population), 2L, replace = TRUE)
        family <- current_population[[chosen[1L]]]$family
        make_child(
          current_population[[chosen[1L]]], current_population[[chosen[2L]]],
          generation, "SIB", family
        )
      })
    } else if (type == "BC") {
      if (is.null(current_population)) stop("BC requires an existing population")
      backcross_generation <- backcross_generation + 1L
      parent_tag <- if (!is.na(parsed$arg1)) toupper(parsed$arg1) else "P1"
      recurrent <- if (parent_tag == "P2") P2 else P1
      generation <- paste0("BC", backcross_generation)
      current_population <- lapply(seq_len(n_offspring), function(i) {
        source <- current_population[[((i - 1L) %% length(current_population)) + 1L]]
        if (!background_selection) {
          return(make_child(source, recurrent, generation, paste0("BC:", parent_tag)))
        }
        candidates <- lapply(seq_len(max(1L, selection_pool)), function(j) {
          make_child(source, recurrent, generation, paste0("BC:", parent_tag))
        })
        scores <- vapply(candidates, function(candidate) {
          score <- donor_fraction(candidate, fixed_locus)
          if (!is.na(introgression_target_len) && !is.null(fixed_locus)) {
            observed <- introgression_length(candidate, fixed_locus)
            if (is.finite(observed)) {
              score <- score + abs(observed - introgression_target_len) /
                max(1, introgression_target_len)
            }
          }
          score
        }, numeric(1L))
        candidates[[which.min(scores)]]
      })
    } else if (type == "DH") {
      if (is.null(current_population)) stop("DH requires an existing population")
      current_population <- lapply(seq_len(n_offspring), function(i) {
        parent <- current_population[[((i - 1L) %% length(current_population)) + 1L]]
        gamete <- make_gamete(parent)
        fixed <- apply_fixed_locus(gamete, gamete)
        list(
          h1 = fixed$g1$sequence, h2 = fixed$g1$sequence,
          a1 = fixed$g1$ancestry, a2 = fixed$g1$ancestry,
          generation = "DH", scheme = "DH", family = parent$family
        )
      })
    } else {
      stop("Unknown breeding sequence token: ", type)
    }
  }

  if (is.null(current_population) || !length(current_population)) {
    stop("Breeding simulation produced no offspring")
  }

  sample_ids <- paste0("sample", seq_along(current_population))
  metadata <- data.frame(
    sample = sample_ids,
    generation = vapply(current_population, `[[`, character(1L), "generation"),
    scheme = vapply(current_population, `[[`, character(1L), "scheme"),
    family = vapply(current_population, `[[`, character(1L), "family"),
    stringsAsFactors = FALSE
  )

  make_sv_truth <- function() {
    n_events <- as.integer(sv_rate * sum(chromosome_lengths))
    empty <- data.frame(
      sv_id = character(), chromosome = character(), position_bp = integer(),
      end_bp = integer(), type = character(), length = integer(),
      alt_founders = character(), inserted_sequence = character(),
      stringsAsFactors = FALSE
    )
    if (n_events <= 0L) return(empty)
    bases <- c("A", "C", "G", "T")
    rows <- vector("list", n_events)
    for (i in seq_len(n_events)) {
      chromosome <- sample(
        chromosomes, 1L, prob = as.numeric(chromosome_lengths)
      )
      chromosome_length <- chromosome_lengths[[chromosome]]
      type <- sample(c("DEL", "INS", "DUP"), 1L)
      event_length <- sample(
        seq_len(max(1L, min(sv_maxlen, chromosome_length - 1L))), 1L
      )
      max_start <- if (type %in% c("DEL", "DUP")) {
        chromosome_length - event_length + 1L
      } else chromosome_length
      position <- sample(seq_len(max_start), 1L)
      alt_count <- sample(seq_len(max(1L, length(founder_ids_all) - 1L)), 1L)
      alt_founders <- sample(founder_ids_all, alt_count, replace = FALSE)
      inserted <- if (type == "INS") {
        paste(sample(bases, event_length, replace = TRUE), collapse = "")
      } else ""
      rows[[i]] <- data.frame(
        sv_id = paste0("sv", i), chromosome = chromosome,
        position_bp = position,
        end_bp = if (type == "INS") position else position + event_length - 1L,
        type = type, length = event_length,
        alt_founders = paste(alt_founders, collapse = ","),
        inserted_sequence = inserted, stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }

  sv_truth <- make_sv_truth()
  dir.create(dirname(sv_truth_out), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    sv_truth, sv_truth_out, sep = "\t", quote = FALSE, row.names = FALSE
  )

  ancestry_founder_at <- function(tracts, position) {
    hit <- which(tracts$start <= position & tracts$end >= position)
    if (!length(hit)) return(NA_character_)
    tracts$founder[hit[1L]]
  }

  carries_sv <- function(ancestry, event) {
    founder <- ancestry_founder_at(
      ancestry[[event$chromosome]], event$position_bp
    )
    founder %in% strsplit(event$alt_founders, ",", fixed = TRUE)[[1L]]
  }

  apply_sv_events <- function(sequence, ancestry, chromosome) {
    events <- sv_truth[sv_truth$chromosome == chromosome, , drop = FALSE]
    if (!nrow(events)) return(sequence)
    events <- events[order(events$position_bp, decreasing = TRUE), , drop = FALSE]
    for (i in seq_len(nrow(events))) {
      event <- events[i, , drop = FALSE]
      if (!carries_sv(ancestry, event)) next
      position <- event$position_bp
      event_length <- event$length
      if (event$type == "DEL") {
        sequence <- paste0(
          substr(sequence, 1L, position - 1L),
          substr(sequence, position + event_length, nchar(sequence))
        )
      } else if (event$type == "INS") {
        sequence <- paste0(
          substr(sequence, 1L, position), event$inserted_sequence,
          substr(sequence, position + 1L, nchar(sequence))
        )
      } else if (event$type == "DUP") {
        duplicated <- substr(sequence, position, position + event_length - 1L)
        sequence <- paste0(
          substr(sequence, 1L, position + event_length - 1L), duplicated,
          substr(sequence, position + event_length, nchar(sequence))
        )
      }
    }
    sequence
  }

  write_fasta <- function(path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    connection <- file(path, "w")
    on.exit(close(connection))
    multi_chromosome <- length(chromosomes) > 1L
    for (i in seq_along(current_population)) {
      for (haplotype in c("h1", "h2")) {
        hap_label <- if (haplotype == "h1") "hap1" else "hap2"
        for (chromosome in chromosomes) {
          identifier <- paste0(sample_ids[i], "_", hap_label)
          if (multi_chromosome) identifier <- paste(identifier, chromosome, sep = "|")
          writeLines(paste0(">", identifier), connection)
          sequence <- apply_sv_events(
            current_population[[i]][[haplotype]][[chromosome]],
            current_population[[i]][[paste0("a", if (haplotype == "h1") 1 else 2)]],
            chromosome
          )
          starts <- seq.int(1L, nchar(sequence), by = 80L)
          writeLines(vapply(starts, function(start) {
            substr(sequence, start, min(nchar(sequence), start + 79L))
          }, character(1L)), connection)
        }
      }
    }
  }

  out_fasta <- paste0(out_prefix, ".fa")
  out_metadata <- paste0(out_prefix, ".meta.tsv")
  write_fasta(out_fasta)
  dir.create(dirname(out_metadata), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    metadata, out_metadata, sep = "\t", quote = FALSE, row.names = FALSE
  )

  ancestry_rows <- list()
  breakpoint_rows <- list()
  row_index <- 1L
  breakpoint_index <- 1L
  for (i in seq_along(current_population)) {
    for (haplotype_number in 1:2) {
      ancestry_name <- paste0("a", haplotype_number)
      for (chromosome in chromosomes) {
        tracts <- current_population[[i]][[ancestry_name]][[chromosome]]
        ancestry_rows[[row_index]] <- data.frame(
          sample = sample_ids[i], haplotype = haplotype_number,
          chromosome = chromosome, start = tracts$start, end = tracts$end,
          founder = tracts$founder, stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L
        if (nrow(tracts) > 1L) {
          breakpoint_rows[[breakpoint_index]] <- data.frame(
            sample = sample_ids[i], haplotype = haplotype_number,
            chromosome = chromosome, position_bp = tracts$end[-nrow(tracts)],
            left_founder = tracts$founder[-nrow(tracts)],
            right_founder = tracts$founder[-1L], stringsAsFactors = FALSE
          )
          breakpoint_index <- breakpoint_index + 1L
        }
      }
    }
  }
  ancestry_table <- do.call(rbind, ancestry_rows)
  if (length(breakpoint_rows)) {
    breakpoint_table <- do.call(rbind, breakpoint_rows)
  } else {
    breakpoint_table <- data.frame(
      sample = character(), haplotype = integer(), chromosome = character(),
      position_bp = integer(), left_founder = character(),
      right_founder = character(), stringsAsFactors = FALSE
    )
  }
  dir.create(dirname(ancestry_out), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(breakpoints_out), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    ancestry_table, ancestry_out, sep = "\t", quote = FALSE, row.names = FALSE
  )
  utils::write.table(
    breakpoint_table, breakpoints_out, sep = "\t", quote = FALSE,
    row.names = FALSE
  )

  if (!is.na(graph_out)) {
    dir.create(dirname(graph_out), recursive = TRUE, showWarnings = FALSE)
    lines <- c(
      "flowchart LR",
      paste0("P1[\"", founder_ids_all[parent_idx[1L]], "\"] --> CROSS[\"",
             ifelse(is.na(scheme), sequence_spec, scheme), "\"]"),
      paste0("P2[\"", founder_ids_all[parent_idx[2L]], "\"] --> CROSS"),
      paste0("CROSS --> FINAL[\"", length(current_population), " final lines; ",
             length(chromosomes), " chromosome(s)\"]")
    )
    writeLines(lines, graph_out)
    if (graph_format %in% c("svg", "png") && nzchar(Sys.which("mmdc"))) {
      rendered <- sub("\\.mmd$", paste0(".", graph_format), graph_out)
      system2("mmdc", c("-i", graph_out, "-o", rendered),
              stdout = FALSE, stderr = FALSE)
    }
  }

  if (!is.na(vcf_out)) {
    dir.create(dirname(vcf_out), recursive = TRUE, showWarnings = FALSE)
    relevant_founders <- if (!is.na(scheme) && scheme %in% c("NAM", "MAGIC")) {
      founder_idx
    } else parent_idx
    variants <- list()
    variant_index <- 1L
    for (chromosome in chromosomes) {
      reference <- genomes[[founder_ids_all[parent_idx[1L]]]][[chromosome]]
      for (position in seq_len(chromosome_lengths[[chromosome]])) {
        final_alleles <- unique(unlist(lapply(current_population, function(individual) {
          c(
            substr(individual$h1[[chromosome]], position, position),
            substr(individual$h2[[chromosome]], position, position)
          )
        })))
        reference_allele <- substr(reference, position, position)
        alternate_alleles <- setdiff(final_alleles, reference_allele)
        if (!length(alternate_alleles)) next
        if (ascertainment == "founders") {
          founder_alleles <- unique(vapply(relevant_founders, function(index) {
            substr(genomes[[founder_ids_all[index]]][[chromosome]], position, position)
          }, character(1L)))
          if (length(founder_alleles) < 2L) next
        }
        variants[[variant_index]] <- data.frame(
          chromosome = chromosome, position = position,
          ref = reference_allele, alt = paste(alternate_alleles, collapse = ","),
          type = "SNP", stringsAsFactors = FALSE
        )
        variant_index <- variant_index + 1L
      }
    }
    variant_table <- if (length(variants)) do.call(rbind, variants) else NULL

    writeLines("##fileformat=VCFv4.2", vcf_out)
    for (chromosome in chromosomes) {
      .simitall_append_lines(
        paste0("##contig=<ID=", chromosome, ",length=",
               chromosome_lengths[[chromosome]], ">"), vcf_out
      )
    }
    .simitall_append_lines(
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>", vcf_out
    )
    .simitall_append_lines(
      "##INFO=<ID=TYPE,Number=1,Type=String,Description=Variant type>", vcf_out
    )
    .simitall_append_lines(
      "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=Structural variant type>",
      vcf_out
    )
    .simitall_append_lines(
      "##INFO=<ID=END,Number=1,Type=Integer,Description=End coordinate>",
      vcf_out
    )
    .simitall_append_lines(
      "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=SV length>",
      vcf_out
    )
    for (i in seq_along(sample_ids)) {
      family <- if (nzchar(metadata$family[i])) metadata$family[i] else "NA"
      .simitall_append_lines(
        paste0("##SAMPLE=<ID=", sample_ids[i], ",FAMILY=", family, ">"),
        vcf_out
      )
    }
    header <- c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", sample_ids
    )
    .simitall_append_lines(paste(header, collapse = "\t"), vcf_out)

    if (!is.null(variant_table)) {
      for (i in seq_len(nrow(variant_table))) {
        chromosome <- variant_table$chromosome[i]
        position <- variant_table$position[i]
        allele_levels <- c(
          variant_table$ref[i],
          strsplit(variant_table$alt[i], ",", fixed = TRUE)[[1L]]
        )
        genotypes <- vapply(current_population, function(individual) {
          alleles <- c(
            substr(individual$h1[[chromosome]], position, position),
            substr(individual$h2[[chromosome]], position, position)
          )
          genotype <- paste(match(alleles, allele_levels) - 1L, collapse = "/")
          if (missing_rate > 0 && stats::runif(1L) < missing_rate) {
            genotype <- "./."
          } else if (genotype_error > 0 && stats::runif(1L) < genotype_error) {
            genotype <- paste(
              sample(seq_along(allele_levels) - 1L, 2L, replace = TRUE),
              collapse = "/"
            )
          }
          genotype
        }, character(1L))
        row <- c(
          chromosome, position, paste0("var", i), variant_table$ref[i],
          variant_table$alt[i], ".", "PASS",
          paste0("TYPE=", variant_table$type[i]), "GT", genotypes
        )
        .simitall_append_lines(paste(row, collapse = "\t"), vcf_out)
      }
    }

    if (nrow(sv_truth)) {
      for (i in seq_len(nrow(sv_truth))) {
        event <- sv_truth[i, , drop = FALSE]
        reference <- substr(
          genomes[[founder_ids_all[parent_idx[1L]]]][[event$chromosome]],
          event$position_bp, event$position_bp
        )
        genotypes <- vapply(current_population, function(individual) {
          alleles <- c(
            as.integer(carries_sv(individual$a1, event)),
            as.integer(carries_sv(individual$a2, event))
          )
          genotype <- paste(alleles, collapse = "/")
          if (missing_rate > 0 && stats::runif(1L) < missing_rate) {
            genotype <- "./."
          } else if (genotype_error > 0 && stats::runif(1L) < genotype_error) {
            genotype <- paste(sample(0:1, 2L, replace = TRUE), collapse = "/")
          }
          genotype
        }, character(1L))
        signed_length <- if (event$type == "DEL") -event$length else event$length
        info <- paste0(
          "TYPE=SV;SVTYPE=", event$type, ";END=", event$end_bp,
          ";SVLEN=", signed_length
        )
        row <- c(
          event$chromosome, event$position_bp, event$sv_id, reference,
          paste0("<", event$type, ">"), ".", "PASS", info, "GT", genotypes
        )
        .simitall_append_lines(paste(row, collapse = "\t"), vcf_out)
      }
    }
  }

  cat("Wrote\n")
  cat("  FASTA:       ", out_fasta, "\n", sep = "")
  cat("  metadata:    ", out_metadata, "\n", sep = "")
  cat("  ancestry:    ", ancestry_out, "\n", sep = "")
  cat("  breakpoints: ", breakpoints_out, "\n", sep = "")
  cat("  recomb map:  ", recomb_map_out, "\n", sep = "")
  cat("  SV truth:    ", sv_truth_out, "\n", sep = "")

  invisible(list(
    fasta = out_fasta, metadata = out_metadata, vcf = vcf_out,
    ancestry = ancestry_out, breakpoints = breakpoints_out,
    recombination_map = recomb_map_out, sv_truth = sv_truth_out
  ))
}
