main_11_simulate_breeding <- function(args = commandArgs(trailingOnly = TRUE)) {

# Simulate breeding populations from haplotype FASTA panels
# Supports F1, selfing, backcross, and arbitrary sequences.

args <- args
usage <- function() {
  cat("Usage: 11_simulate_breeding.R --haplotype_fa <path> --out_prefix <path> [options]
")
  cat("
Required:
")
  cat("  --haplotype_fa <path>        multi-FASTA panel (same length)
")
  cat("  --out_prefix <path>          output prefix (.fa, .meta.tsv)
")
  cat("  --vcf_out <path>             output VCF of bred lines
")
  cat("  --graph_out <path>           output Mermaid crossing graph
")
  cat("  --graph_format <mmd|svg|png>  default mmd (requires mmdc for svg/png)
")
  cat("
Parents / population:
")
  cat("  --parents <csv>              two haplotype IDs or indices (default random)
")
  cat("  --n_offspring <int>           default 100
")
  cat("  --founders <csv>             founder hap IDs/indices for MAGIC/NAM
")
  cat("  --n_founders <int>            number of founders (if founders not provided)
")
  cat("
Breeding sequence:
")
  cat("  --sequence <list>            e.g., F1,SELF:3,SIB:2,BC:P1:2,DH
")
  cat("  --scheme <F2|MAGIC|NAM|RIL|NIL|DH> optional high-level scheme
")
  cat("  --self_generations <int>      for RIL/NIL (default 6)
")
  cat("  --backcross_generations <int> for NIL (default 3)
")
  cat("  --ril_mating <SSD|SIB>        default SSD
")
  cat("  --fix_locus <start:end>        force locus from donor/recipient
")
  cat("  --fix_allele <donor|recipient|hapID>
")
  cat("  --background_selection        select NILs by minimal donor background
")
  cat("  --selection_pool <int>        candidates per generation (default 50)
")
  cat("  --marker_step <int>           marker spacing for selection (default 1000)
")
  cat("  --introgression_target_len <int> target donor tract length (bp)
")
  cat("  --genotype_error <float>      per-genotype error rate (default 0)
")
  cat("  --missing_rate <float>        per-genotype missingness (default 0)
  --qc_out <path>               write QC report (TSV)
  --qc_ld_bins <int>            LD bins (default 20)
  --qc_ld_maxdist <int>         LD max distance bp (default 100000)
")
  cat("  --ascertainment <founders|all> default founders
")
  cat("  --sv_rate <float>             structural variant rate (default 0)
")
  cat("  --sv_maxlen <int>             max SV length (default 1000)
")
  cat("
Recombination map:
")
  cat("  --recomb_map_in <path>        recomb TSV (pos_bp,cM)
")
  cat("  --recomb_rate_mean <float>    cM/Mb mean (default 1.0)
")
  cat("  --recomb_rate_sd <float>      cM/Mb SD (default 0.3)
")
  cat("  --recomb_hotspots <int>       number of hotspots (default 3)
")
  cat("  --recomb_hotspot_mult <float> hotspot multiplier (default 5.0)
")
  cat("  --interference_shape <float>  gamma shape for crossover interference (default 1.0)
  --selection_loci <csv>        selected loci positions (bp)
  --selection_model <add|dom|rec> default add
  --selection_strength <float>  fitness penalty per risk allele (default 0)
  --distortion_rate <float>     segregation distortion rate (default 0)
")
  cat("
Misc:
")
  cat("  --seed <int>                  default 1
")
  quit(status = 1)
}


get_arg <- function(flag, default = NULL) {
  if (!(flag %in% args)) return(default)
  idx <- match(flag, args)
  if (idx == length(args)) return(default)
  args[idx + 1]
}

has_flag <- function(flag) flag %in% args

parse_founders <- function(ids, founders_arg, n_founders) {
  if (!is.na(founders_arg)) {
    parts <- strsplit(founders_arg, ",")[[1]]
    idx <- sapply(parts, function(p) {
      if (grepl("^[0-9]+$", p)) as.integer(p) else match(p, ids)
    })
    if (any(is.na(idx))) stop("Founder IDs not found")
    return(idx)
  }
  if (n_founders > length(ids)) stop("n_founders > available haplotypes")
  sample(seq_along(ids), n_founders, replace = FALSE)
}


hap_fa <- get_arg("--haplotype_fa")
out_prefix <- get_arg("--out_prefix")
if (is.null(hap_fa) || is.null(out_prefix)) usage()

parents_arg <- get_arg("--parents", NA)
n_offspring <- as.integer(get_arg("--n_offspring", 100))
sequence <- get_arg("--sequence", "F1,SELF:1")
scheme <- get_arg("--scheme", NA)
founders_arg <- get_arg("--founders", NA)
n_founders <- as.integer(get_arg("--n_founders", 4))
vcf_out <- get_arg("--vcf_out", NA)
graph_out <- get_arg("--graph_out", NA)
graph_format <- get_arg("--graph_format", "mmd")

recomb_map_in <- get_arg("--recomb_map_in", NA)
recomb_rate_mean <- as.numeric(get_arg("--recomb_rate_mean", 1.0))
recomb_rate_sd <- as.numeric(get_arg("--recomb_rate_sd", 0.3))
recomb_hotspots <- as.integer(get_arg("--recomb_hotspots", 3))
recomb_hotspot_mult <- as.numeric(get_arg("--recomb_hotspot_mult", 5.0))

interference_shape <- as.numeric(get_arg("--interference_shape", 1.0))
self_generations <- as.integer(
  get_arg("--self_generations", get_arg("--generations", 6))
)
backcross_generations <- as.integer(get_arg("--backcross_generations", 3))
ril_mating <- toupper(get_arg("--ril_mating", "SSD"))
fix_locus <- get_arg("--fix_locus", NA)
fix_allele <- get_arg("--fix_allele", NA)
background_selection <- has_flag("--background_selection")
selection_pool <- as.integer(get_arg("--selection_pool", 50))
marker_step <- as.integer(get_arg("--marker_step", 1000))
introgression_target_len <- as.numeric(
  get_arg("--introgression_target_len", NA)
)
genotype_error <- as.numeric(get_arg("--genotype_error", 0))
missing_rate <- as.numeric(get_arg("--missing_rate", 0))
ascertainment <- get_arg("--ascertainment", "founders")
sv_rate <- as.numeric(get_arg("--sv_rate", 0))
sv_maxlen <- as.integer(get_arg("--sv_maxlen", 1000))
selection_loci <- get_arg("--selection_loci", NA)
selection_model <- tolower(get_arg("--selection_model", "add"))
selection_strength <- as.numeric(get_arg("--selection_strength", 0))
distortion_rate <- as.numeric(get_arg("--distortion_rate", 0))

seed <- as.integer(get_arg("--seed", 1))
set.seed(seed)

read_fasta_multi <- function(path) {
  if (!file.exists(path)) {
    stop("Haplotype FASTA not found: ", path)
  }
  lines <- readLines(path)
  header_idx <- grep("^>", lines)
  if (!length(header_idx)) {
    stop("No FASTA records found in: ", path)
  }

  ids <- sub("^>", "", lines[header_idx])
  ids <- sub("\\s.*$", "", ids)
  seqs <- vapply(seq_along(header_idx), function(i) {
    start <- header_idx[i] + 1L
    end <- if (i < length(header_idx)) header_idx[i + 1L] - 1L else length(lines)
    if (start > end) return("")
    toupper(paste(lines[start:end], collapse = ""))
  }, character(1))

  list(ids = ids, seqs = seqs)
}

read_recomb_map <- function(path, genome_length) {
  map <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("pos_bp", "cM")
  if (!all(required %in% names(map))) {
    stop("Recombination map must contain columns: pos_bp and cM")
  }
  map <- map[order(map$pos_bp), required]
  map <- map[!duplicated(map$pos_bp), , drop = FALSE]
  map$pos_bp <- pmax(1, pmin(genome_length, as.integer(map$pos_bp)))
  map$cM <- as.numeric(map$cM)
  if (map$pos_bp[1L] > 1L) {
    map <- rbind(data.frame(pos_bp = 1L, cM = 0), map)
  }
  if (tail(map$pos_bp, 1L) < genome_length) {
    map <- rbind(
      map,
      data.frame(
        pos_bp = genome_length,
        cM = tail(map$cM, 1L)
      )
    )
  }
  map
}

make_recomb_map <- function(
    genome_length,
    rate_mean,
    rate_sd,
    hotspots,
    hotspot_multiplier) {
  n_intervals <- max(2L, min(1000L, ceiling(genome_length / 1000)))
  positions <- unique(as.integer(round(seq(1, genome_length, length.out = n_intervals))))
  rates <- pmax(0, stats::rnorm(length(positions) - 1L, rate_mean, rate_sd))

  if (hotspots > 0L && length(rates)) {
    hot <- sample(seq_along(rates), min(hotspots, length(rates)))
    rates[hot] <- rates[hot] * hotspot_multiplier
  }

  delta_mb <- diff(positions) / 1e6
  data.frame(
    pos_bp = positions,
    cM = c(0, cumsum(rates * delta_mb))
  )
}

parse_locus <- function(value) {
  if (is.na(value) || !nzchar(value)) {
    return(NULL)
  }
  parts <- as.integer(strsplit(value, "[:-]")[[1L]])
  if (length(parts) != 2L || anyNA(parts) || parts[1L] > parts[2L]) {
    stop("--fix_locus must use start:end coordinates")
  }
  c(start = parts[1L], end = parts[2L])
}

apply_fix_locus <- function(h1, h2, donor, locus) {
  if (is.null(donor) || is.null(locus)) {
    return(list(h1 = h1, h2 = h2))
  }

  start <- max(1L, locus[["start"]])
  end <- min(nchar(donor), locus[["end"]])
  segment <- substr(donor, start, end)
  replace_segment <- function(sequence) {
    paste0(
      substr(sequence, 1L, start - 1L),
      segment,
      substr(sequence, end + 1L, nchar(sequence))
    )
  }

  list(h1 = replace_segment(h1), h2 = replace_segment(h2))
}

make_gamete <- function(h1, h2, recomb_map, interference_shape = 1) {
  total_cM <- tail(recomb_map$cM, 1)
  expected_xo <- total_cM / 100
  if (!is.finite(total_cM) || total_cM <= 0) {
    return(if (runif(1) < 0.5) h1 else h2)
  }
  if (interference_shape <= 1) {
    n_xo <- rpois(1, expected_xo)
    if (n_xo == 0) return(if (runif(1) < 0.5) h1 else h2)
    xo_cM <- sort(runif(n_xo, 0, total_cM))
  } else {
    # gamma interference model
    mean_dist <- ifelse(expected_xo > 0, total_cM / (expected_xo + 1), total_cM)
    scale <- mean_dist / interference_shape
    xo_cM <- c()
    pos <- 0
    while (pos < total_cM) {
      step <- rgamma(1, shape = interference_shape, scale = scale)
      pos <- pos + step
      if (pos < total_cM) xo_cM <- c(xo_cM, pos)
    }
    if (length(xo_cM) == 0) return(if (runif(1) < 0.5) h1 else h2)
  }

  xo_pos <- sapply(xo_cM, function(x) {
    idx <- max(which(recomb_map$cM <= x))
    recomb_map$pos_bp[idx]
  })

  # build recombinant haplotype
  segments <- c(1, xo_pos, glen + 1)
  use_h1 <- runif(1) < 0.5
  out <- character()
  for (i in seq_len(length(segments) - 1)) {
    start <- segments[i]
    end <- segments[i+1] - 1
    if (end < start) next
    if (use_h1) out <- c(out, substr(h1, start, end)) else out <- c(out, substr(h2, start, end))
    use_h1 <- !use_h1
  }
  paste(out, collapse = "")
}

# selection helpers
parse_positions <- function(x) {
  if (is.na(x)) return(integer())
  as.integer(strsplit(x, ",")[[1]])
}

selection_loci_pos <- parse_positions(selection_loci)

selection_score <- function(h1, h2, loci, model) {
  if (length(loci) == 0) return(0)
  risk <- 0
  for (pos in loci) {
    a1 <- substr(h1, pos, pos); a2 <- substr(h2, pos, pos)
    # risk allele defined as donor (P2) allele at that pos
    if (a1 == substr(P2$h1, pos, pos)) risk <- risk + 1
    if (a2 == substr(P2$h1, pos, pos)) risk <- risk + 1
  }
  if (model == "dom") return(as.integer(risk > 0))
  if (model == "rec") return(as.integer(risk >= 2))
  risk
}

# Load haplotypes
panel <- read_fasta_multi(hap_fa)
ids <- panel$ids
seqs <- panel$seqs
if (length(seqs) < 2) stop("Need at least 2 haplotypes")

lens <- nchar(seqs)
if (length(unique(lens)) != 1) stop("All haplotypes must have same length")

glen <- lens[1]

recomb_map <- if (!is.na(recomb_map_in) && file.exists(recomb_map_in)) {
  read_recomb_map(recomb_map_in, glen)
} else {
  make_recomb_map(glen, recomb_rate_mean, recomb_rate_sd, recomb_hotspots, recomb_hotspot_mult)
}

# resolve parents
founder_idx <- parse_founders(ids, founders_arg, n_founders)
founder_ids <- ids[founder_idx]
founder_seqs <- seqs[founder_idx]

if (is.na(parents_arg)) {
  parent_idx <- sample(seq_along(ids), 2, replace = FALSE)
} else {
  parts <- strsplit(parents_arg, ",")[[1]]
  if (length(parts) != 2) stop("--parents must have two values")
  # allow indices or IDs
  parent_idx <- sapply(parts, function(p) {
    if (grepl("^[0-9]+$", p)) as.integer(p) else match(p, ids)
  })
  if (any(is.na(parent_idx))) stop("parent IDs not found")
}

P1 <- list(id = ids[parent_idx[1]], h1 = seqs[parent_idx[1]], h2 = seqs[parent_idx[1]])
P2 <- list(id = ids[parent_idx[2]], h1 = seqs[parent_idx[2]], h2 = seqs[parent_idx[2]])

# fixed locus donor
loc <- parse_locus(fix_locus)
fixed_donor <- NULL
if (!is.na(fix_allele)) {
  if (fix_allele == "donor") fixed_donor <- P2$h1
  else if (fix_allele == "recipient") fixed_donor <- P1$h1
  else {
    idx <- match(fix_allele, ids)
    if (is.na(idx)) stop("fix_allele hapID not found")
    fixed_donor <- seqs[idx]
  }
}

# selection helpers
marker_positions <- seq(1, glen, by = max(1, marker_step))
donor_seq <- if (!is.null(fixed_donor)) fixed_donor else P2$h1

donor_fraction <- function(h1, h2, donor, positions, locus = NULL) {
  idx <- positions
  if (!is.null(locus)) {
    idx <- positions[positions < locus["start"] | positions > locus["end"]]
  }
  if (length(idx) == 0) return(0)
  m1 <- sum(sapply(idx, function(i) substr(h1, i, i) == substr(donor, i, i)))
  m2 <- sum(sapply(idx, function(i) substr(h2, i, i) == substr(donor, i, i)))
  (m1 + m2) / (2 * length(idx))
}

donor_tract_len <- function(h, donor, locus) {
  if (is.null(locus)) return(NA_integer_)
  s <- locus["start"]; e <- locus["end"]
  # expand left
  left <- s
  while (left > 1 && substr(h, left-1, left-1) == substr(donor, left-1, left-1)) left <- left - 1
  right <- e
  while (right < nchar(h) && substr(h, right+1, right+1) == substr(donor, right+1, right+1)) right <- right + 1
  right - left + 1
}

# parse sequence tokens
parse_token <- function(tok) {
  parts <- strsplit(tok, ":")[[1]]
  list(type = parts[1], arg1 = ifelse(length(parts) >= 2, parts[2], NA), arg2 = ifelse(length(parts) >= 3, parts[3], NA))
}

sequence_tokens <- strsplit(sequence, ",")[[1]]
sequence_tokens <- trimws(sequence_tokens)

expand_generations <- function(tokens) {
  expanded <- character()
  for (token in tokens) {
    parsed <- parse_token(token)
    if (parsed$type %in% c("SELF", "SIB") && !is.na(parsed$arg1)) {
      expanded <- c(expanded, rep(parsed$type, as.integer(parsed$arg1)))
    } else if (parsed$type == "BC" && !is.na(parsed$arg2)) {
      expanded <- c(
        expanded,
        rep(paste("BC", parsed$arg1, sep = ":"), as.integer(parsed$arg2))
      )
    } else {
      expanded <- c(expanded, token)
    }
  }
  expanded
}

# population: list of individuals (h1,h2)
pop <- list()
meta <- data.frame(sample = character(), generation = character(), scheme = character(), family = character(), stringsAsFactors = FALSE)

make_offspring <- function(n, parentA, parentB, gen_label, scheme, family = "") {
  out <- vector("list", n)
  i <- 1
  while (i <= n) {
    g1 <- make_gamete(parentA$h1, parentA$h2, recomb_map, interference_shape)
    g2 <- make_gamete(parentB$h1, parentB$h2, recomb_map, interference_shape)
    offs <- apply_fix_locus(g1, g2, fixed_donor, loc)

    # segregation distortion: bias toward donor allele at selected loci
    if (!is.na(distortion_rate) && distortion_rate > 0 && length(selection_loci_pos) > 0) {
      for (pos in selection_loci_pos) {
        if (runif(1) < distortion_rate) {
          donor_a <- substr(P2$h1, pos, pos)
          offs$h1 <- paste0(substr(offs$h1, 1, pos-1), donor_a, substr(offs$h1, pos+1, nchar(offs$h1)))
          offs$h2 <- paste0(substr(offs$h2, 1, pos-1), donor_a, substr(offs$h2, pos+1, nchar(offs$h2)))
        }
      }
    }

    # selection: accept/reject offspring based on risk allele count
    if (!is.na(selection_strength) && selection_strength > 0 && length(selection_loci_pos) > 0) {
      risk <- selection_score(offs$h1, offs$h2, selection_loci_pos, selection_model)
      fitness <- exp(-selection_strength * risk)
      if (runif(1) > fitness) next
    }

    out[[i]] <- list(h1 = offs$h1, h2 = offs$h2)
    meta <<- rbind(meta, data.frame(sample = paste0(gen_label, "_", i), generation = gen_label, scheme = scheme, family = family, stringsAsFactors = FALSE))
    i <- i + 1
  }
  out
}

# Prebuilt schemes
if (!is.na(scheme)) {
  if (scheme == "MAGIC") {
    parents <- lapply(seq_along(founder_seqs), function(i) list(id = founder_ids[i], h1 = founder_seqs[i], h2 = founder_seqs[i]))
    while (length(parents) > 1) {
      a <- parents[[1]]; b <- parents[[2]]
      f1 <- make_offspring(1, a, b, "MAGIC_F1", "MAGIC")[[1]]
      parents <- c(list(list(h1 = f1$h1, h2 = f1$h2)), parents[-c(1,2)])
    }
    current_pop <- list(parents[[1]])
    sequence_tokens <- paste0("SELF:", self_generations)
  } else if (scheme == "NAM") {
    common <- list(id = founder_ids[1], h1 = founder_seqs[1], h2 = founder_seqs[1])
    fams <- list()
    for (i in 2:length(founder_seqs)) {
      other <- list(id = founder_ids[i], h1 = founder_seqs[i], h2 = founder_seqs[i])
      f1 <- make_offspring(1, common, other, paste0("NAM_F1_", i), "NAM")[[1]]
      fams[[i-1]] <- f1
    }
    per_fam <- max(1, floor(n_offspring / length(fams)))
    current_pop <- list()
    for (i in seq_along(fams)) {
      parent <- fams[[i]]
      for (k in seq_len(per_fam)) {
        g1 <- make_gamete(parent$h1, parent$h2, recomb_map, interference_shape)
        g2 <- make_gamete(parent$h1, parent$h2, recomb_map, interference_shape)
        current_pop[[length(current_pop)+1]] <- list(h1 = g1, h2 = g2)
        meta <- rbind(meta, data.frame(sample = paste0("NAM", i, "_", k), generation = "NAM", scheme = "NAM", family = paste0("F", i), stringsAsFactors = FALSE))
      }
    }
    sequence_tokens <- character()
  } else if (scheme == "F2") {
    sequence_tokens <- c("F1", "SELF")
  } else if (scheme == "RIL") {
    if (toupper(ril_mating) == "SIB") {
      sequence_tokens <- c("F1", paste0("SIB:", self_generations))
    } else {
      sequence_tokens <- c("F1", paste0("SELF:", self_generations))
    }
  } else if (scheme == "NIL") {
    sequence_tokens <- c("F1", paste0("BC:P1:", backcross_generations), paste0("SELF:", self_generations))
  } else if (scheme == "DH") {
    sequence_tokens <- c("F1", "DH")
  }
}

# Run sequence
sequence_tokens <- expand_generations(sequence_tokens)
if (!exists("current_pop", inherits = FALSE)) {
  current_pop <- NULL
}
for (tok in sequence_tokens) {
  t <- parse_token(tok)
  typ <- t$type
  n_gen <- ifelse(!is.na(t$arg2), as.integer(t$arg2), ifelse(!is.na(t$arg1) && typ %in% c("F1", "SELF", "SIB", "DH"), as.integer(t$arg1), n_offspring))

  if (typ == "F1") {
    n <- ifelse(!is.na(t$arg1), as.integer(t$arg1), n_offspring)
    current_pop <- make_offspring(n, P1, P2, "F1", "F1")
  } else if (typ == "SELF") {
    if (is.null(current_pop)) stop("SELF requires existing population")
    n <- ifelse(!is.na(t$arg1), as.integer(t$arg1), n_offspring)
    # selfing: sample individuals and self
    new_pop <- vector("list", n)
    for (i in seq_len(n)) {
      parent <- current_pop[[sample(seq_along(current_pop), 1)]]
      g1 <- make_gamete(parent$h1, parent$h2, recomb_map, interference_shape)
      g2 <- make_gamete(parent$h1, parent$h2, recomb_map, interference_shape)
      offs <- apply_fix_locus(g1, g2, fixed_donor, loc)
      new_pop[[i]] <- list(h1 = offs$h1, h2 = offs$h2)
      meta <- rbind(meta, data.frame(sample = paste0("SELF", "_", i), generation = "SELF", scheme = "SELF", family = "", stringsAsFactors = FALSE))
    }
    current_pop <- new_pop
  } else if (typ == "SIB") {
    if (is.null(current_pop)) stop("SIB requires existing population")
    n <- ifelse(!is.na(t$arg1), as.integer(t$arg1), n_offspring)
    new_pop <- vector("list", n)
    for (i in seq_len(n)) {
      pidx <- sample(seq_along(current_pop), 2, replace = TRUE)
      p1 <- current_pop[[pidx[1]]]
      p2 <- current_pop[[pidx[2]]]
      g1 <- make_gamete(p1$h1, p1$h2, recomb_map, interference_shape)
      g2 <- make_gamete(p2$h1, p2$h2, recomb_map, interference_shape)
      offs <- apply_fix_locus(g1, g2, fixed_donor, loc)
      new_pop[[i]] <- list(h1 = offs$h1, h2 = offs$h2)
      meta <- rbind(meta, data.frame(sample = paste0("SIB", "_", i), generation = "SIB", scheme = "SIB", family = "", stringsAsFactors = FALSE))
    }
    current_pop <- new_pop
  } else if (typ == "DH") {
    if (is.null(current_pop)) stop("DH requires existing population")
    n <- ifelse(!is.na(t$arg1), as.integer(t$arg1), n_offspring)
    new_pop <- vector("list", n)
    for (i in seq_len(n)) {
      parent <- current_pop[[sample(seq_along(current_pop), 1)]]
      g1 <- make_gamete(parent$h1, parent$h2, recomb_map, interference_shape)
      offs <- apply_fix_locus(g1, g1, fixed_donor, loc)
      new_pop[[i]] <- list(h1 = offs$h1, h2 = offs$h2)
      meta <- rbind(meta, data.frame(sample = paste0("DH", "_", i), generation = "DH", scheme = "DH", family = "", stringsAsFactors = FALSE))
    }
    current_pop <- new_pop
  } else if (typ == "BC") {
    if (is.null(current_pop)) stop("BC requires existing population")
    parent_tag <- ifelse(!is.na(t$arg1), t$arg1, "P1")
    back_parent <- if (parent_tag == "P2") P2 else P1
    n <- ifelse(!is.na(t$arg2), as.integer(t$arg2), n_offspring)
    new_pop <- vector("list", n)
    for (i in seq_len(n)) {
      if (background_selection) {
        pool <- vector("list", selection_pool)
        scores <- numeric(selection_pool)
        for (k in seq_len(selection_pool)) {
          parent <- current_pop[[sample(seq_along(current_pop), 1)]]
          g1 <- make_gamete(parent$h1, parent$h2, recomb_map, interference_shape)
          g2 <- make_gamete(back_parent$h1, back_parent$h2, recomb_map, interference_shape)
          offs <- apply_fix_locus(g1, g2, fixed_donor, loc)
          pool[[k]] <- list(h1 = offs$h1, h2 = offs$h2)
          scores[k] <- donor_fraction(offs$h1, offs$h2, donor_seq, marker_positions, loc)
          if (!is.na(introgression_target_len)) {
            tlen <- mean(c(donor_tract_len(g1, donor_seq, loc), donor_tract_len(g2, donor_seq, loc)), na.rm = TRUE)
            scores[k] <- scores[k] + abs(tlen - introgression_target_len) / introgression_target_len
          }
        }
        best <- which.min(scores)
        new_pop[[i]] <- pool[[best]]
      } else {
        parent <- current_pop[[sample(seq_along(current_pop), 1)]]
        g1 <- make_gamete(parent$h1, parent$h2, recomb_map, interference_shape)
        g2 <- make_gamete(back_parent$h1, back_parent$h2, recomb_map, interference_shape)
        offs <- apply_fix_locus(g1, g2, fixed_donor, loc)
        new_pop[[i]] <- list(h1 = offs$h1, h2 = offs$h2)
      }
      meta <- rbind(meta, data.frame(sample = paste0("BC", "_", i), generation = "BC", scheme = paste0("BC:", parent_tag), family = "", stringsAsFactors = FALSE))
    }
    current_pop <- new_pop
  } else {
    stop("Unknown token in --sequence: ", typ)
  }
}

apply_sv_tracts <- function(hap_seqs, sv_rate, sv_maxlen) {
  if (is.na(sv_rate) || sv_rate <= 0) return(hap_seqs)
  L <- nchar(hap_seqs[1])
  n_sv <- as.integer(sv_rate * L)
  if (n_sv <= 0) return(hap_seqs)
  for (k in seq_len(n_sv)) {
    pos <- sample(2:(L-1), 1)
    len <- sample(1:sv_maxlen, 1)
    typ <- sample(c("DEL","DUP","INS"), 1)
    for (i in seq_along(hap_seqs)) {
      h <- hap_seqs[i]
      if (typ == "DEL") {
        h <- paste0(substr(h, 1, pos-1), substr(h, pos+len, nchar(h)))
      } else if (typ == "DUP") {
        seg <- substr(h, pos, pos+len-1)
        h <- paste0(substr(h, 1, pos-1), seg, seg, substr(h, pos+len, nchar(h)))
      } else {
        ins <- paste(sample(c("A","C","G","T"), len, replace=TRUE), collapse="")
        h <- paste0(substr(h, 1, pos-1), ins, substr(h, pos, nchar(h)))
      }
      hap_seqs[i] <- h
    }
  }
  hap_seqs
}

# Write output FASTA (diploid as two haplotypes per sample)
write_fasta <- function(path, ids, seqs, width = 80) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, "w")
  on.exit(close(con))
  for (i in seq_along(ids)) {
    writeLines(paste0(">", ids[i]), con)
    s <- seqs[i]
    for (j in seq(1, nchar(s), by = width)) {
      writeLines(substr(s, j, min(nchar(s), j + width - 1)), con)
    }
  }
}

out_fa <- paste0(out_prefix, ".fa")
out_meta <- paste0(out_prefix, ".meta.tsv")

if (is.null(current_pop) || !length(current_pop)) {
  stop("Breeding simulation produced no offspring.")
}

# Keep metadata aligned one-to-one with the final population.
if (nrow(meta) >= length(current_pop)) {
  meta <- utils::tail(meta, length(current_pop))
} else {
  meta <- data.frame(
    sample = paste0("sample", seq_along(current_pop)),
    generation = if (is.na(scheme)) "custom" else scheme,
    scheme = if (is.na(scheme)) "custom" else scheme,
    family = "",
    stringsAsFactors = FALSE
  )
}
meta$sample <- paste0("sample", seq_along(current_pop))

hap_ids <- character()
hap_seqs <- character()
for (i in seq_along(current_pop)) {
  hap_ids <- c(hap_ids, paste0("sample", i, "_hap1"), paste0("sample", i, "_hap2"))
  hap_seqs <- c(hap_seqs, current_pop[[i]]$h1, current_pop[[i]]$h2)
}

hap_seqs <- apply_sv_tracts(hap_seqs, sv_rate, sv_maxlen)
write_fasta(out_fa, hap_ids, hap_seqs)
# Mermaid crossing graph (expanded)
if (!is.na(graph_out)) {
  dir.create(dirname(graph_out), showWarnings = FALSE, recursive = TRUE)
  con <- file(graph_out, "w")
  writeLines("flowchart LR", con)
  # parents
  writeLines(paste0("P1[\"", P1$id, "\"]"), con)
  writeLines(paste0("P2[\"", P2$id, "\"]"), con)
  if (!is.na(scheme) && scheme == "MAGIC") {
    writeLines("P1 --> MAGIC[\"MAGIC\"]", con)
    writeLines("P2 --> MAGIC", con)
  } else if (!is.na(scheme) && scheme == "NAM") {
    writeLines("P1 --> NAM[\"NAM\"]", con)
  } else {
    writeLines("P1 --> F1[\"F1\"]", con)
    writeLines("P2 --> F1", con)
  }
  # expanded generations from meta
  if (nrow(meta) > 0) {
    for (i in seq_len(nrow(meta))) {
      node <- meta$sample[i]
      gen <- meta$generation[i]
      fam <- meta$family[i]
      label <- ifelse(is.na(fam) || fam == "", gen, paste0(gen, ":", fam))
      writeLines(paste0(gen, "[", label, "]"), con)
      writeLines(paste0(gen, " --> ", node, "[", node, "]"), con)
    }
  }
  close(con)
  # optional render via mermaid-cli (mmdc)
  if (graph_format %in% c("svg", "png")) {
    out_path <- sub("\\.mmd$", paste0(".", graph_format), graph_out)
    cmd <- paste("mmdc -i", shQuote(graph_out), "-o", shQuote(out_path))
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  }
}

write.table(meta, out_meta, sep = "\t", row.names = FALSE, quote = FALSE)

# Write VCF (SNPs + optional SVs)
if (!is.na(vcf_out)) {
  ref <- hap_seqs[1]
  glen <- nchar(ref)
  bases <- strsplit(ref, "")[[1]]
  variants <- data.frame(pos = integer(), ref = character(), alt = character(), type = character(), stringsAsFactors = FALSE)
  for (pos in seq_len(glen)) {
    alleles <- unique(sapply(hap_seqs, function(h) substr(h, pos, pos)))
    if (length(alleles) > 1) {
      refb <- bases[pos]
      altb <- setdiff(alleles, refb)
      if (length(altb) > 0) {
        variants <- rbind(variants, data.frame(pos = pos, ref = refb, alt = altb[1], type = "SNP", stringsAsFactors = FALSE))
      }
    }
  }
  # ascertainment: keep only variants polymorphic between founders
  if (tolower(ascertainment) == "founders" && nrow(variants) > 0L) {
    p1 <- P1$h1; p2 <- P2$h1
    keep <- vapply(
      variants$pos,
      function(pos) substr(p1, pos, pos) != substr(p2, pos, pos),
      logical(1L)
    )
    variants <- variants[keep, , drop = FALSE]
  }

  # add SVs as symbolic alleles
  if (!is.na(sv_rate) && sv_rate > 0) {
    n_sv <- as.integer(sv_rate * glen)
    if (n_sv > 0) {
      sv_pos <- sample(seq_len(glen), n_sv, replace = FALSE)
      for (pos in sv_pos) {
        sv_type <- ifelse(runif(1) < 0.5, "DEL", "DUP")
        variants <- rbind(variants, data.frame(pos = pos, ref = "N", alt = paste0("<", sv_type, ">"), type = "SV", stringsAsFactors = FALSE))
      }
    }
  }

  writeLines("##fileformat=VCFv4.2", vcf_out)
  .simitall_append_lines(
    paste0("##contig=<ID=chr1,length=", glen, ">"),
    vcf_out
  )
  .simitall_append_lines(
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    vcf_out
  )

  samples <- unique(meta$sample)
  # per-sample labels
  for (i in seq_along(samples)) {
    fam <- ifelse(is.na(meta$family[i]) || meta$family[i] == "", "NA", meta$family[i])
    .simitall_append_lines(
      paste0("##SAMPLE=<ID=", samples[i], ",FAMILY=", fam, ">"),
      vcf_out
    )
  }

  header <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", samples)
  .simitall_append_lines(paste(header, collapse = "\t"), vcf_out)

  for (i in seq_len(nrow(variants))) {
    pos <- variants$pos[i]
    gt <- character(length(samples))
    for (s in seq_along(samples)) {
      h1 <- hap_seqs[(s-1)*2 + 1]
      h2 <- hap_seqs[(s-1)*2 + 2]
      a1 <- substr(h1, pos, pos)
      a2 <- substr(h2, pos, pos)
      if (a1 == variants$ref[i] && a2 == variants$ref[i]) gt[s] <- "0/0"
      else if (a1 == variants$alt[i] && a2 == variants$alt[i]) gt[s] <- "1/1"
      else gt[s] <- "0/1"

      # missingness / error
      if (!is.na(missing_rate) && runif(1) < missing_rate) {
        gt[s] <- "./."
      } else if (!is.na(genotype_error) && runif(1) < genotype_error) {
        gt[s] <- sample(c("0/0","0/1","1/1"), 1)
      }
    }
    row <- c("chr1", variants$pos[i], paste0("var", i), variants$ref[i], variants$alt[i], ".", "PASS", variants$type[i], "GT", gt)
    .simitall_append_lines(paste(row, collapse = "\t"), vcf_out)
  }
}

cat("Wrote\n")
cat("  FASTA: ", out_fa, "\n", sep = "")
cat("  META:  ", out_meta, "\n", sep = "")
invisible(list(fasta = out_fa, metadata = out_meta, vcf = vcf_out))
}
