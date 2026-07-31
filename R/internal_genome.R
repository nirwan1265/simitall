main_01_make_tandem_repeats <- function(args = commandArgs(trailingOnly = TRUE)) {

# Insert tandem repeats and/or motif repeats into a genome FASTA.
# Modes: tandem | motif | both
# Supports random genome generation, GC-preserve, spacing distributions.
# Outputs: TSV, BED, GFF3, summary JSON, header labels, name mapping.

args <- args
usage <- function() {
  cat("Usage: 01_make_tandem_repeats.R --out_fa <path> [options]\n")
  cat("\nRequired:\n")
  cat("  --out_fa <path>\n")
  cat("\nInput genome:\n")
  cat("  --in_fa <path>                (required unless --random_genome)\n")
  cat("  --contig <name>                optional contig name\n")
  cat("  --random_genome                generate random genome\n")
  cat("  --random_length <int>          default 100000\n")
  cat("  --random_gc <float>            default 0.5\n")
  cat("\nRepeat modes:\n")
  cat("  --mode <tandem|motif|both>      default tandem\n")
  cat("\nTandem params:\n")
  cat("  --n_events <int>               default 10\n")
  cat("  --seg_len <int>                default 1000\n")
  cat("  --copies <int>                 default 5\n")
  cat("\nMotif params:\n")
  cat("  --motif <string>               default ATTA\n")
  cat("  --motif_repeat <int>           default 1000\n")
  cat("  --motif_events <int>           default 5\n")
  cat("  --motif_mode <insert|replace>  default insert\n")
  cat("  --motif_gc_target <float>      generate random motif with target GC\n")
  cat("  --motif_len <int>              default 12 (if using motif_gc_target)\n")
  cat("\nSpacing / GC:\n")
  cat("  --min_spacing <int>            default 0\n")
  cat("  --spacing_distribution <uniform|poisson|fixed>  default uniform\n")
  cat("  --spacing_mean <float>         default 10000 (poisson)\n")
  cat("  --spacing_fixed <int>          default 1000 (fixed interval)\n")
  cat("  --gc_preserve                  attempt GC-preserving motif insert\n")
  cat("\nOutputs:\n")
  cat("  --coords_tsv <path>            write insertion coordinates TSV\n")
  cat("  --coords_bed <path>            write BED coordinates\n")
  cat("  --coords_gff3 <path>           write GFF3 coordinates\n")
  cat("  --summary_json <path>          write JSON summary of parameters + counts\n")
  cat("  --label_blocks                 add block summary to FASTA header\n")
  cat("  --name_map_tsv <path>          write mapping of original name to output name\n")
  cat("\nMisc:\n")
  cat("  --seed <int>                   default 1\n")
  quit(status = 1)
}

get_arg <- function(flag, default = NULL) {
  if (!(flag %in% args)) return(default)
  idx <- match(flag, args)
  if (idx == length(args)) return(default)
  args[idx + 1]
}

has_flag <- function(flag) {
  flag %in% args
}

out_fa <- get_arg("--out_fa")
if (is.null(out_fa)) usage()

in_fa <- get_arg("--in_fa")
contig <- get_arg("--contig", NA)
random_genome <- has_flag("--random_genome")
random_length <- as.integer(get_arg("--random_length", 100000))
random_gc <- as.numeric(get_arg("--random_gc", 0.5))

ploidy <- as.integer(get_arg("--ploidy", 2))
ploidy_mode <- get_arg("--ploidy_mode", "identical")
ploidy_snp_rate <- as.numeric(get_arg("--ploidy_snp_rate", 0.001))
ploidy_indel_rate <- as.numeric(get_arg("--ploidy_indel_rate", 0.0001))
ploidy_indel_maxlen <- as.integer(get_arg("--ploidy_indel_maxlen", 3))
ploidy_ref_copy <- as.integer(get_arg("--ploidy_ref_copy", 1))

mode <- get_arg("--mode", "tandem")

n_events <- as.integer(get_arg("--n_events", 10))
seg_len <- as.integer(get_arg("--seg_len", 1000))
copies <- as.integer(get_arg("--copies", 5))

motif <- get_arg("--motif", "ATTA")
motif_repeat <- as.integer(get_arg("--motif_repeat", 1000))
motif_events <- as.integer(get_arg("--motif_events", 5))
motif_mode <- get_arg("--motif_mode", "insert")

min_spacing <- as.integer(get_arg("--min_spacing", 0))
spacing_distribution <- get_arg("--spacing_distribution", "uniform")
spacing_mean <- as.numeric(get_arg("--spacing_mean", 10000))
spacing_fixed <- as.integer(get_arg("--spacing_fixed", 1000))

gc_preserve <- has_flag("--gc_preserve")
seed <- as.integer(get_arg("--seed", 1))

coords_tsv <- get_arg("--coords_tsv", NA)
coords_bed <- get_arg("--coords_bed", NA)
coords_gff3 <- get_arg("--coords_gff3", NA)
coords_tsv_per_copy <- get_arg("--coords_tsv_per_copy", NA)
coords_bed_per_copy <- get_arg("--coords_bed_per_copy", NA)
coords_gff3_per_copy <- get_arg("--coords_gff3_per_copy", NA)
summary_json <- get_arg("--summary_json", NA)
label_blocks <- has_flag("--label_blocks")
name_map_tsv <- get_arg("--name_map_tsv", NA)

motif_gc_target <- get_arg("--motif_gc_target")
motif_len <- as.integer(get_arg("--motif_len", 12))

set.seed(seed)

# ploidy checks
if (is.na(ploidy) || ploidy < 1) stop("--ploidy must be >= 1")
if (!(ploidy_mode %in% c("identical", "diverged"))) stop("--ploidy_mode must be identical or diverged")
if (ploidy_snp_rate < 0 || ploidy_snp_rate > 1) stop("--ploidy_snp_rate must be 0-1")
if (ploidy_indel_rate < 0 || ploidy_indel_rate > 1) stop("--ploidy_indel_rate must be 0-1")
if (is.na(ploidy_indel_maxlen) || ploidy_indel_maxlen < 1) ploidy_indel_maxlen <- 1

read_fasta_one <- function(path, contig = NA) {
  if (!file.exists(path)) stop("FASTA not found: ", path)
  lines <- readLines(path)
  hdr_idx <- grep("^>", lines)
  if (length(hdr_idx) == 0) stop("No FASTA header found")
  if (is.na(contig)) {
    name <- sub("^>", "", strsplit(lines[hdr_idx[1]], " ")[[1]][1])
    seq <- paste(lines[(hdr_idx[1] + 1):(if (length(hdr_idx) > 1) hdr_idx[2] - 1 else length(lines))], collapse = "")
    return(list(name = name, seq = toupper(seq)))
  }
  for (i in seq_along(hdr_idx)) {
    name <- sub("^>", "", strsplit(lines[hdr_idx[i]], " ")[[1]][1])
    if (name == contig) {
      start <- hdr_idx[i] + 1
      end <- if (i < length(hdr_idx)) hdr_idx[i + 1] - 1 else length(lines)
      seq <- paste(lines[start:end], collapse = "")
      return(list(name = name, seq = toupper(seq)))
    }
  }
  stop("Contig not found: ", contig)
}

write_fasta <- function(path, name, seq, width = 80) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, "w")
  on.exit(close(con))
  writeLines(paste0(">", name), con)
  for (i in seq(1, nchar(seq), by = width)) {
    writeLines(substr(seq, i, min(nchar(seq), i + width - 1)), con)
  }
}

write_fasta_multi <- function(path, names, seqs, width = 80) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, "w")
  on.exit(close(con))
  for (i in seq_along(names)) {
    writeLines(paste0(">", names[i]), con)
    for (j in seq(1, nchar(seqs[i]), by = width)) {
      writeLines(substr(seqs[i], j, min(nchar(seqs[i]), j + width - 1)), con)
    }
  }
}

gc_fraction <- function(s) {
  s <- toupper(s)
  if (nchar(s) == 0) return(0)
  gc <- sum(strsplit(s, "")[[1]] %in% c("G", "C"))
  gc / nchar(s)
}

make_random_genome <- function(length, gc_target) {
  if (length <= 0) stop("--random_length must be > 0")
  if (gc_target < 0 || gc_target > 1) stop("--random_gc must be between 0 and 1")
  bases <- sapply(runif(length), function(x) {
    if (x < gc_target) sample(c("G", "C"), 1) else sample(c("A", "T"), 1)
  })
  paste(bases, collapse = "")
}

mutate_sequence <- function(seq, snp_rate, indel_rate, indel_maxlen) {
  bases <- c("A", "C", "G", "T")
  chars <- strsplit(seq, "")[[1]]
  out <- character()
  i <- 1
  n <- length(chars)
  while (i <= n) {
    base <- chars[i]
    if (runif(1) < snp_rate) {
      base <- sample(bases[bases != base], 1)
    }
    if (runif(1) < indel_rate) {
      if (runif(1) < 0.5) {
        ins_len <- sample(1:indel_maxlen, 1)
        ins <- sample(bases, ins_len, replace = TRUE)
        out <- c(out, base, ins)
        i <- i + 1
        next
      } else {
        del_len <- sample(1:indel_maxlen, 1)
        i <- i + del_len
        next
      }
    }
    out <- c(out, base)
    i <- i + 1
  }
  paste(out, collapse = "")
}

make_gc_matched_motif <- function(length, gc_target) {
  if (length <= 0) stop("--motif_len must be > 0")
  if (gc_target < 0 || gc_target > 1) stop("--motif_gc_target must be between 0 and 1")
  bases <- sapply(runif(length), function(x) {
    if (x < gc_target) sample(c("G", "C"), 1) else sample(c("A", "T"), 1)
  })
  paste(bases, collapse = "")
}

pick_positions_uniform <- function(len, n, min_spacing) {
  if (n <= 0) return(integer(0))
  if (min_spacing <= 0) return(sample(0:len, n, replace = TRUE))
  pos <- integer(0)
  attempts <- 0
  while (length(pos) < n && attempts < max(1000, n * 500)) {
    attempts <- attempts + 1
    p <- sample(0:len, 1)
    if (all(abs(p - pos) >= min_spacing)) pos <- c(pos, p)
  }
  if (length(pos) < n) stop("Could not place positions with min_spacing; reduce spacing or n")
  pos
}

pick_positions_poisson <- function(len, n, mean_gap, min_spacing) {
  if (n <= 0) return(integer(0))
  if (mean_gap <= 0) return(pick_positions_uniform(len, n, min_spacing))
  pos <- integer(0)
  p <- sample(0:len, 1)
  pos <- c(pos, p)
  attempts <- 0
  while (length(pos) < n && attempts < max(1000, n * 500)) {
    attempts <- attempts + 1
    gap <- rexp(1, rate = 1 / mean_gap)
    p <- as.integer(p + gap)
    if (p > len) p <- p %% max(1, len)
    if (all(abs(p - pos) >= min_spacing)) pos <- c(pos, p)
  }
  if (length(pos) < n) stop("Could not place positions with poisson spacing; reduce spacing or n")
  pos
}

pick_positions_fixed <- function(len, n, fixed_gap, min_spacing) {
  if (n <= 0) return(integer(0))
  if (fixed_gap <= 0) return(pick_positions_uniform(len, n, min_spacing))
  start <- sample(0:len, 1)
  pos <- start + (0:(n - 1)) * fixed_gap
  pos <- pos %% max(1, len)
  if (min_spacing > 0) {
    filtered <- integer(0)
    for (p in pos) {
      if (all(abs(p - filtered) >= min_spacing)) filtered <- c(filtered, p)
    }
    if (length(filtered) < n) stop("Fixed spacing failed with min_spacing; reduce min_spacing or n")
    return(filtered)
  }
  pos
}

remove_gc_matched_segment <- function(seq, target_gc, seg_len, tol = 0.05, max_tries = 2000) {
  if (seg_len <= 0 || seg_len > nchar(seq)) return(seq)
  for (i in 1:max_tries) {
    start <- sample(1:(nchar(seq) - seg_len + 1), 1)
    seg <- substr(seq, start, start + seg_len - 1)
    if (abs(gc_fraction(seg) - target_gc) <= tol) {
      return(paste0(substr(seq, 1, start - 1), substr(seq, start + seg_len, nchar(seq))))
    }
  }
  start <- sample(1:(nchar(seq) - seg_len + 1), 1)
  paste0(substr(seq, 1, start - 1), substr(seq, start + seg_len, nchar(seq)))
}

if (random_genome) {
  name <- "RandomGenome"
  seq <- make_random_genome(random_length, random_gc)
} else {
  if (is.null(in_fa)) stop("--in_fa is required unless --random_genome is set")
  res <- read_fasta_one(in_fa, ifelse(is.na(contig), NA, contig))
  name <- res$name
  seq <- res$seq
}

if (mode %in% c("tandem", "both")) {
  if (copies < 2) stop("--copies must be >= 2")
  if (seg_len <= 0 || seg_len > nchar(seq)) stop("--seg_len must be between 1 and genome length")
}

if (mode %in% c("motif", "both")) {
  if (!is.null(motif_gc_target)) {
    motif <- make_gc_matched_motif(motif_len, as.numeric(motif_gc_target))
  }
  if (motif_repeat < 1 || motif_events < 1) stop("--motif_repeat and --motif_events must be >= 1")
}

get_positions <- function(len, n) {
  if (spacing_distribution == "poisson") {
    return(pick_positions_poisson(len, n, spacing_mean, min_spacing))
  } else if (spacing_distribution == "fixed") {
    return(pick_positions_fixed(len, n, spacing_fixed, min_spacing))
  }
  pick_positions_uniform(len, n, min_spacing)
}

coords <- data.frame(
  mode = character(),
  index = integer(),
  insert_pos = integer(),
  start = integer(),
  end = integer(),
  block_len = integer(),
  type = character(),
  stringsAsFactors = FALSE
)

s <- seq
if (mode %in% c("motif", "both")) {
  motif_block <- paste(rep(toupper(motif), motif_repeat), collapse = "")
  positions <- sort(get_positions(nchar(s), motif_events), decreasing = TRUE)
  idx <- 1
  for (p in positions) {
    start <- p + 1
    end <- start + nchar(motif_block) - 1
    if (motif_mode == "replace") {
      s <- paste0(substr(s, 1, p), motif_block, substr(s, p + nchar(motif_block) + 1, nchar(s)))
    } else {
      s <- paste0(substr(s, 1, p), motif_block, substr(s, p + 1, nchar(s)))
      if (gc_preserve) {
        s <- remove_gc_matched_segment(s, gc_fraction(motif_block), nchar(motif_block))
      }
    }
    coords <- rbind(coords, data.frame(mode = "motif", index = idx, insert_pos = p, start = start, end = end,
                                       block_len = nchar(motif_block), type = "motif", stringsAsFactors = FALSE))
    idx <- idx + 1
  }
}

if (mode %in% c("tandem", "both")) {
  positions <- sort(get_positions(nchar(s), n_events), decreasing = TRUE)
  idx <- 1
  for (p in positions) {
    start <- p + 1
    seg_start <- sample(1:(nchar(s) - seg_len + 1), 1)
    seg <- substr(s, seg_start, seg_start + seg_len - 1)
    block <- paste(rep(seg, copies), collapse = "")
    end <- start + nchar(block) - 1
    s <- paste0(substr(s, 1, p), block, substr(s, p + 1, nchar(s)))
    coords <- rbind(coords, data.frame(mode = "tandem", index = idx, insert_pos = p, start = start, end = end,
                                       block_len = nchar(block), type = "tandem", stringsAsFactors = FALSE))
    idx <- idx + 1
  }
}

out_name <- paste0(
  name,
  "|mode", mode,
  "|tandem_events", n_events, "_len", seg_len, "_copies", copies,
  "|motif_", motif, "_rep", motif_repeat, "_events", motif_events,
  "|motifmode", motif_mode, "_gc_preserve", as.integer(gc_preserve),
  "|spacing_", spacing_distribution, "_mean", spacing_mean, "_fixed", spacing_fixed, "_min", min_spacing
)
if (random_genome) {
  out_name <- paste0(out_name, "|random_len", random_length, "_gc", random_gc)
}
if (!is.null(motif_gc_target)) {
  out_name <- paste0(out_name, "|motif_gc_target", motif_gc_target, "_len", motif_len)
}

if (!is.na(ploidy)) {
  out_name <- paste0(out_name, "|ploidy", ploidy, "_", ploidy_mode)
}
if (label_blocks) {
  out_name <- paste0(out_name, "|blocks_motif", sum(coords$type == "motif"), "_tandem", sum(coords$type == "tandem"))
}

# build ploidy copies
contig_names <- c(out_name)
contig_seqs <- c(s)
if (ploidy > 1) {
  contig_names <- paste0(out_name, "|copy", seq_len(ploidy))
  contig_seqs <- character(length(contig_names))
  contig_seqs[1] <- s
  if (ploidy_mode == "identical") {
    for (i in 2:ploidy) contig_seqs[i] <- s
  } else {
    for (i in 2:ploidy) contig_seqs[i] <- mutate_sequence(s, ploidy_snp_rate, ploidy_indel_rate, ploidy_indel_maxlen)
  }
}

if (length(contig_names) == 1) {
  write_fasta(out_fa, contig_names[1], contig_seqs[1])
} else {
  write_fasta_multi(out_fa, contig_names, contig_seqs)
}


coord_chrom <- if (ploidy <= 1) out_name else paste0(out_name, "|copy", ploidy_ref_copy)
coords$contig <- coord_chrom

if (!is.na(coords_tsv)) {
  dir.create(dirname(coords_tsv), showWarnings = FALSE, recursive = TRUE)
  write.table(coords, coords_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
}

if (!is.na(coords_bed)) {
  dir.create(dirname(coords_bed), showWarnings = FALSE, recursive = TRUE)
  bed <- data.frame(
    chrom = coord_chrom,
    chromStart = coords$start - 1,
    chromEnd = coords$end,
    name = paste0(coords$type, "_", coords$index),
    score = 0,
    strand = ".",
    stringsAsFactors = FALSE
  )
  write.table(bed, coords_bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# per-copy coordinate outputs
if (!is.na(coords_tsv_per_copy) && ploidy > 1) {
  for (i in seq_len(ploidy)) {
    outp <- paste0(coords_tsv_per_copy, ".copy", i)
    dir.create(dirname(outp), showWarnings = FALSE, recursive = TRUE)
    cc <- coords
    cc$contig <- paste0(out_name, "|copy", i)
    write.table(cc, outp, sep = "	", row.names = FALSE, quote = FALSE)
  }
}
if (!is.na(coords_bed_per_copy) && ploidy > 1) {
  for (i in seq_len(ploidy)) {
    outp <- paste0(coords_bed_per_copy, ".copy", i)
    dir.create(dirname(outp), showWarnings = FALSE, recursive = TRUE)
    bed <- data.frame(
      chrom = paste0(out_name, "|copy", i),
      chromStart = coords$start - 1,
      chromEnd = coords$end,
      name = paste0(coords$type, "_", coords$index),
      score = 0,
      strand = "."
    )
    write.table(bed, outp, sep = "	", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}
if (!is.na(coords_gff3_per_copy) && ploidy > 1) {
  for (i in seq_len(ploidy)) {
    outp <- paste0(coords_gff3_per_copy, ".copy", i)
    dir.create(dirname(outp), showWarnings = FALSE, recursive = TRUE)
    writeLines("##gff-version 3", outp)
    gff <- data.frame(
      seqid = paste0(out_name, "|copy", i),
      source = "simulator",
      type = coords$type,
      start = coords$start,
      end = coords$end,
      score = ".",
      strand = ".",
      phase = ".",
      attributes = paste0("ID=", coords$type, "_", coords$index)
    )
    write.table(gff, outp, sep = "	", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  }
}

if (!is.na(coords_gff3)) {
  dir.create(dirname(coords_gff3), showWarnings = FALSE, recursive = TRUE)
  writeLines("##gff-version 3", coords_gff3)
  gff <- data.frame(
    seqid = coord_chrom,
    source = "simulator",
    type = coords$type,
    start = coords$start,
    end = coords$end,
    score = ".",
    strand = ".",
    phase = ".",
    attributes = paste0("ID=", coords$type, "_", coords$index, ";Length=", coords$block_len),
    stringsAsFactors = FALSE
  )
  write.table(
    gff,
    coords_gff3,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    append = TRUE
  )
}

if (!is.na(summary_json)) {
  dir.create(dirname(summary_json), showWarnings = FALSE, recursive = TRUE)
  summary <- list(
    params = list(
      in_fa = in_fa,
      out_fa = out_fa,
      mode = mode,
      n_events = n_events,
      seg_len = seg_len,
      copies = copies,
      motif = motif,
      motif_repeat = motif_repeat,
      motif_events = motif_events,
      motif_mode = motif_mode,
      motif_gc_target = motif_gc_target,
      motif_len = motif_len,
      min_spacing = min_spacing,
      spacing_distribution = spacing_distribution,
      spacing_mean = spacing_mean,
      spacing_fixed = spacing_fixed,
      gc_preserve = gc_preserve,
      seed = seed,
      random_genome = random_genome,
      random_length = random_length,
      random_gc = random_gc,
      ploidy = ploidy,
      ploidy_mode = ploidy_mode,
      ploidy_snp_rate = ploidy_snp_rate,
      ploidy_indel_rate = ploidy_indel_rate,
      ploidy_indel_maxlen = ploidy_indel_maxlen,
      ploidy_ref_copy = ploidy_ref_copy
    ),
    counts = list(
      motif_blocks = sum(coords$type == "motif"),
      tandem_blocks = sum(coords$type == "tandem"),
      total_blocks = nrow(coords),
      ploidy = ploidy,
      copy_lengths = sapply(contig_seqs, nchar)
    )
  )
  jsonlite::write_json(
    summary,
    path = summary_json,
    auto_unbox = TRUE,
    pretty = TRUE,
    na = "null"
  )
}

if (!is.na(name_map_tsv)) {
  dir.create(dirname(name_map_tsv), showWarnings = FALSE, recursive = TRUE)
  map <- data.frame(original_name = name, output_name = contig_names, stringsAsFactors = FALSE)
  write.table(map, name_map_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
}
}
