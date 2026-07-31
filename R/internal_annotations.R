main_01b_generate_annotations <- function(args = commandArgs(trailingOnly = TRUE)) {

# Generate gene models, CDS, operons, TSS, promoters, rRNA/tRNA clusters, terminators.
# Optional splice-like exon/intron features for demonstration.
# Optional riboswitches, CRISPR arrays, origin/terminus GC-skew signal.
# Optional plasmid contig simulation + plasmid GFF3.
# Optional riboswitch/CRISPR FASTA outputs.
# Optional plasmid marker genes (quick random or curated FASTA).

args <- args
usage <- function() {
  cat("Usage: 01b_generate_annotations.R --genome_fa <path> --out_gff3 <path> [options]\n")
  cat("\nRequired:\n")
  cat("  --genome_fa <path>            genome FASTA (single-contig preferred)\n")
  cat("  --out_gff3 <path>             output GFF3\n")
  cat("\nInput annotation (optional):\n")
  cat("  --in_gff3 <path>              if provided, use this GFF3 and extract sequences\n")
  cat("\nRandom model options (used if no --in_gff3):\n")
  cat("  --n_genes <int>               default 400\n")
  cat("  --gene_len_mean <int>         default 900\n")
  cat("  --gene_len_sd <int>           default 200\n")
  cat("  --intergenic_min <int>        default 20\n")
  cat("  --intergenic_max <int>        default 200\n")
  cat("  --operon_prob <float>         probability next gene joins current operon (default 0.6)\n")
  cat("  --tss_per_operon <int>        default 1\n")
  cat("  --promoter_len <int>          default 60\n")
  cat("  --strand_bias <float>         probability of '+' strand (default 0.5)\n")
  cat("\nAdditional features:\n")
  cat("  --rrna_clusters <int>         number of rRNA clusters (default 1)\n")
  cat("  --trna_per_cluster <int>      tRNAs per cluster (default 10)\n")
  cat("  --terminator_prob <float>     probability of terminator after operon (default 0.8)\n")
  cat("  --riboswitch_count <int>      number of riboswitch features (default 10)\n")
  cat("  --crispr_count <int>          number of CRISPR arrays (default 1)\n")
  cat("  --crispr_len <int>            length of CRISPR array (default 1200)\n")
  cat("  --gc_skew_signal              add origin/terminus GC-skew features\n")
  cat("\nPlasmid simulation:\n")
  cat("  --plasmid_count <int>         number of plasmids (default 0)\n")
  cat("  --plasmid_length <int>        length of each plasmid (default 20000)\n")
  cat("  --plasmid_gc <float>          plasmid GC fraction (default 0.5)\n")
  cat("  --plasmid_fa_out <path>       write plasmid FASTA (multi-contig)\n")
  cat("  --plasmid_gff3_out <path>     write plasmid GFF3\n")
  cat("\nPlasmid markers:\n")
  cat("  --marker_mode <quick|curated|both> default both\n")
  cat("  --markers_fa <path>           curated marker FASTA\n")
  cat("  --markers_per_plasmid <int>   default 2\n")
  cat("  --marker_insert_mode <replace|insert> default replace\n")
  cat("\nRegulatory elements (optional):\n")
  cat("  --regulatory_panel <path>     TSV panel of enhancers/silencers\n")
  cat("  --regulatory_default_human    use the bundled human regulatory panel\n")
  cat("  --regulatory_fa <path>        optional FASTA of element sequences\n")
  cat("  --regulatory_fa_out <path>    write enhancer/silencer FASTA\n")
  cat("  --regulatory_tsv_out <path>   write regulatory element table\n")
  cat("  --enhancer_gc <float>         default 0.55 (synthetic elements)\n")
  cat("  --silencer_gc <float>         default 0.45 (synthetic elements)\n")
  cat("\nSplice-like features (optional):\n")
  cat("  --add_splicing                add exon/intron features for a subset of genes\n")
  cat("  --splice_prob <float>         probability a gene gets exons (default 0.1)\n")
  cat("  --min_exons <int>             default 2\n")
  cat("  --max_exons <int>             default 4\n")
  cat("\nOutputs:\n")
  cat("  --genes_fa <path>             output gene FASTA\n")
  cat("  --cds_fa <path>               output CDS FASTA (same as genes here)\n")
  cat("  --promoters_fa <path>         output promoter FASTA\n")
  cat("  --tss_tsv <path>              output TSS table\n")
  cat("  --operons_tsv <path>          output operon membership\n")
  cat("  --rrna_fa <path>              output rRNA FASTA\n")
  cat("  --trna_fa <path>              output tRNA FASTA\n")
  cat("  --riboswitch_fa <path>        output riboswitch FASTA\n")
  cat("  --crispr_fa <path>            output CRISPR FASTA\n")
  cat("  --regulatory_fa_out <path>    output enhancer/silencer FASTA\n")
  cat("\nMisc:\n")
  cat("  --seed <int>                  default 1\n")
  quit(status = 1)
}

get_arg <- function(flag, default = NULL) {
  if (!(flag %in% args)) return(default)
  idx <- match(flag, args)
  if (idx == length(args)) return(default)
  args[idx + 1]
}

has_flag <- function(flag) flag %in% args

# Required
fasta_path <- get_arg("--genome_fa")
out_gff3 <- get_arg("--out_gff3")
if (is.null(fasta_path) || is.null(out_gff3)) usage()

in_gff3 <- get_arg("--in_gff3")

# Random model options
n_genes <- as.integer(get_arg("--n_genes", 400))
len_mean <- as.integer(get_arg("--gene_len_mean", 900))
len_sd <- as.integer(get_arg("--gene_len_sd", 200))
ig_min <- as.integer(get_arg("--intergenic_min", 20))
ig_max <- as.integer(get_arg("--intergenic_max", 200))
operon_prob <- as.numeric(get_arg("--operon_prob", 0.6))
tss_per_operon <- as.integer(get_arg("--tss_per_operon", 1))
prom_len <- as.integer(get_arg("--promoter_len", 60))
strand_bias <- as.numeric(get_arg("--strand_bias", 0.5))

# Additional features
rrna_clusters <- as.integer(get_arg("--rrna_clusters", 1))
trna_per_cluster <- as.integer(get_arg("--trna_per_cluster", 10))
terminator_prob <- as.numeric(get_arg("--terminator_prob", 0.8))
riboswitch_count <- as.integer(get_arg("--riboswitch_count", 10))
crispr_count <- as.integer(get_arg("--crispr_count", 1))
crispr_len <- as.integer(get_arg("--crispr_len", 1200))
gc_skew_signal <- has_flag("--gc_skew_signal")

# Plasmids
plasmid_count <- as.integer(get_arg("--plasmid_count", 0))
plasmid_length <- as.integer(get_arg("--plasmid_length", 20000))
plasmid_gc <- as.numeric(get_arg("--plasmid_gc", 0.5))
plasmid_fa_out <- get_arg("--plasmid_fa_out", NA)
plasmid_gff3_out <- get_arg("--plasmid_gff3_out", NA)

regulatory_panel <- get_arg("--regulatory_panel", NA)
regulatory_default_human <- has_flag("--regulatory_default_human")
if (regulatory_default_human) {
  regulatory_panel <- .simitall_extdata(
    "regulatory",
    "regulatory_panel_human.tsv"
  )
}
regulatory_fa <- get_arg("--regulatory_fa", NA)
regulatory_fa_out <- get_arg("--regulatory_fa_out", NA)
enhancer_gc <- as.numeric(get_arg("--enhancer_gc", 0.55))
silencer_gc <- as.numeric(get_arg("--silencer_gc", 0.45))

# Plasmid markers
marker_mode <- get_arg("--marker_mode", "both")
markers_fa <- get_arg("--markers_fa", NA)
markers_per_plasmid <- as.integer(get_arg("--markers_per_plasmid", 2))
marker_insert_mode <- get_arg("--marker_insert_mode", "replace")

# Splicing-like features
add_splicing <- has_flag("--add_splicing")
splice_prob <- as.numeric(get_arg("--splice_prob", 0.1))
min_exons <- as.integer(get_arg("--min_exons", 2))
max_exons <- as.integer(get_arg("--max_exons", 4))

# Outputs
genes_fa <- get_arg("--genes_fa", NA)
cds_fa <- get_arg("--cds_fa", NA)
prom_fa <- get_arg("--promoters_fa", NA)
tss_tsv <- get_arg("--tss_tsv", NA)
operons_tsv <- get_arg("--operons_tsv", NA)
rrna_fa <- get_arg("--rrna_fa", NA)
trna_fa <- get_arg("--trna_fa", NA)
riboswitch_fa <- get_arg("--riboswitch_fa", NA)
crispr_fa <- get_arg("--crispr_fa", NA)

seed <- as.integer(get_arg("--seed", 1))
set.seed(seed)

read_fasta_one <- function(path) {
  lines <- readLines(path)
  hdr_idx <- grep("^>", lines)
  if (length(hdr_idx) == 0) stop("No FASTA header found")
  name <- sub("^>", "", strsplit(lines[hdr_idx[1]], " ")[[1]][1])
  seq <- paste(lines[(hdr_idx[1] + 1):(if (length(hdr_idx) > 1) hdr_idx[2] - 1 else length(lines))], collapse = "")
  list(name = name, seq = toupper(seq))
}

read_fasta_multi <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    return(list(ids = character(), seqs = character()))
  }
  lines <- readLines(path)
  idx <- grep("^>", lines)
  if (length(idx) == 0) return(list(ids = character(), seqs = character()))
  ids <- sub("^>", "", lines[idx])
  seqs <- sapply(seq_along(idx), function(i) {
    start <- idx[i] + 1
    end <- if (i < length(idx)) idx[i + 1] - 1 else length(lines)
    paste(lines[start:end], collapse = "")
  })
  list(ids = ids, seqs = toupper(seqs))
}

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

revcomp <- function(seq) {
  comp <- chartr("ACGTacgt", "TGCAtgca", seq)
  paste(rev(strsplit(comp, "")[[1]]), collapse = "")
}

parse_gff3 <- function(path) {
  lines <- readLines(path)
  lines <- lines[!grepl("^#", lines)]
  df <- read.delim(text = lines, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  colnames(df) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  df
}

make_random_genome <- function(length, gc_target) {
  if (length <= 0) stop("plasmid length must be > 0")
  if (gc_target < 0 || gc_target > 1) stop("plasmid GC must be between 0 and 1")
  bases <- sapply(runif(length), function(x) {
    if (x < gc_target) sample(c("G", "C"), 1) else sample(c("A", "T"), 1)
  })
  paste(bases, collapse = "")
}

# Quick marker panel (name -> length)
quick_markers <- list(
  blaTEM1 = 861,
  nptII = 795,
  cat = 660,
  tetA = 1200,
  aadA = 789,
  aacC1 = 528,
  repA = 900,
  parA = 900,
  parB = 750,
  oriV = 400
)

# Build marker pool
get_marker_pool <- function() {
  ids <- character()
  seqs <- character()
  if (marker_mode %in% c("curated", "both")) {
    cur <- read_fasta_multi(markers_fa)
    if (length(cur$ids) > 0) {
      ids <- c(ids, cur$ids)
      seqs <- c(seqs, cur$seqs)
    }
  }
  if (marker_mode %in% c("quick", "both") || length(ids) == 0) {
    for (k in names(quick_markers)) {
      len <- quick_markers[[k]]
      seq <- make_random_genome(len, 0.5)
      ids <- c(ids, k)
      seqs <- c(seqs, seq)
    }
  }
  list(ids = ids, seqs = seqs)
}

# Insert markers into plasmid sequence
insert_markers <- function(seq, marker_ids, marker_seqs) {
  glen <- nchar(seq)
  feats <- data.frame(id = character(), start = integer(), end = integer(), stringsAsFactors = FALSE)
  if (length(marker_ids) == 0) return(list(seq = seq, feats = feats))

  chosen <- sample(seq_along(marker_ids), min(markers_per_plasmid, length(marker_ids)), replace = FALSE)
  for (i in chosen) {
    mseq <- marker_seqs[i]
    mlen <- nchar(mseq)
    if (mlen >= glen) next
    pos <- sample(1:(glen - mlen + 1), 1)
    if (marker_insert_mode == "replace") {
      seq <- paste0(substr(seq, 1, pos - 1), mseq, substr(seq, pos + mlen, glen))
    } else {
      seq <- paste0(substr(seq, 1, pos - 1), mseq, substr(seq, pos, glen))
      glen <- nchar(seq)
    }
    feats <- rbind(feats, data.frame(id = marker_ids[i], start = pos, end = pos + mlen - 1, stringsAsFactors = FALSE))
  }
  list(seq = seq, feats = feats)
}

if (is.null(in_gff3)) {
  ref <- read_fasta_one(fasta_path)
  seqname <- ref$name
  genome <- ref$seq
  glen <- nchar(genome)

  gene_lens <- pmax(150, as.integer(rnorm(n_genes, len_mean, len_sd)))
  gene_lens[gene_lens %% 3 != 0] <- gene_lens[gene_lens %% 3 != 0] + (3 - gene_lens[gene_lens %% 3 != 0])

  starts <- integer(n_genes)
  ends <- integer(n_genes)
  strands <- character(n_genes)
  pos <- 1
  for (i in 1:n_genes) {
    gap <- sample(ig_min:ig_max, 1)
    pos <- pos + gap
    if (pos + gene_lens[i] > glen) break
    starts[i] <- pos
    ends[i] <- pos + gene_lens[i] - 1
    strands[i] <- ifelse(runif(1) < strand_bias, "+", "-")
    pos <- ends[i]
  }
  valid <- which(starts > 0)
  starts <- starts[valid]
  ends <- ends[valid]
  strands <- strands[valid]
  n_genes <- length(starts)

  operon_id <- integer(n_genes)
  current_op <- 1
  operon_id[1] <- current_op
  for (i in 2:n_genes) {
    if (strands[i] == strands[i-1] && runif(1) < operon_prob) {
      operon_id[i] <- current_op
    } else {
      current_op <- current_op + 1
      operon_id[i] <- current_op
    }
  }

  tss_rows <- data.frame(operon = integer(), tss = integer(), strand = character(), stringsAsFactors = FALSE)
  prom_rows <- data.frame(id = character(), start = integer(), end = integer(), strand = character(), stringsAsFactors = FALSE)

  for (op in unique(operon_id)) {
    first_idx <- which(operon_id == op)[1]
    if (strands[first_idx] == "+") {
      tss_pos <- max(1, starts[first_idx] - sample(20:50, 1))
      prom_start <- max(1, tss_pos - prom_len + 1)
      prom_end <- tss_pos
    } else {
      tss_pos <- min(glen, ends[first_idx] + sample(20:50, 1))
      prom_start <- tss_pos
      prom_end <- min(glen, tss_pos + prom_len - 1)
    }
    tss_rows <- rbind(tss_rows, data.frame(operon = op, tss = tss_pos, strand = strands[first_idx]))
    prom_rows <- rbind(prom_rows, data.frame(id = paste0("promoter_", op), start = prom_start, end = prom_end, strand = strands[first_idx]))
  }

  gff <- data.frame(
    seqid = character(), source = character(), type = character(), start = integer(), end = integer(),
    score = character(), strand = character(), phase = character(), attributes = character(),
    stringsAsFactors = FALSE
  )

  for (i in 1:n_genes) {
    gid <- paste0("gene_", i)
    gff <- rbind(gff, data.frame(seqid = seqname, source = "simulator", type = "gene",
                                 start = starts[i], end = ends[i], score = ".", strand = strands[i], phase = ".",
                                 attributes = paste0("ID=", gid, ";Operon=", operon_id[i])))
    gff <- rbind(gff, data.frame(seqid = seqname, source = "simulator", type = "CDS",
                                 start = starts[i], end = ends[i], score = ".", strand = strands[i], phase = "0",
                                 attributes = paste0("ID=cds_", i, ";Parent=", gid)))


  # Regulatory elements (enhancers/silencers)
  reg_table <- data.frame(
    gene_symbol = character(), target_gene_index = integer(), element_type = character(),
    offset_bp = integer(), length_bp = integer(), start = integer(), end = integer(),
    sequence_id = character(), stringsAsFactors = FALSE
  )
  reg_ids <- character()
  reg_seqs <- character()
  if (!is.na(regulatory_panel) && file.exists(regulatory_panel)) {
    reg <- read.delim(regulatory_panel, header = TRUE, sep = "	", stringsAsFactors = FALSE, check.names = FALSE)
    reg_fa <- list(ids = character(), seqs = character())
    if (!is.na(regulatory_fa) && file.exists(regulatory_fa)) {
      reg_fa <- read_fasta_multi(regulatory_fa)
    }

    if (nrow(reg) > 0) {
      for (i in seq_len(nrow(reg))) {
        tgi <- as.integer(reg$target_gene_index[i])
        if (is.na(tgi) || tgi < 1 || tgi > n_genes) next
        etype <- tolower(reg$element_type[i])
        if (!(etype %in% c("enhancer", "silencer"))) next
        offset <- as.integer(reg$offset_bp[i])
        if (is.na(offset)) offset <- 0
        elen <- as.integer(reg$length_bp[i])
        if (is.na(elen) || elen < 1) next

        if (offset < 0) {
          estart <- starts[tgi] + offset
        } else {
          estart <- ends[tgi] + offset
        }
        eend <- estart + elen - 1
        if (estart < 1 || eend > glen) next

        seq_id <- if (!is.null(reg$sequence_id)) reg$sequence_id[i] else NA
        seq_val <- NA
        if (!is.na(seq_id) && length(reg_fa$ids) > 0) {
          idx <- which(reg_fa$ids == seq_id)[1]
          if (!is.na(idx)) seq_val <- reg_fa$seqs[idx]
        }
        if (is.na(seq_val) || seq_val == "") {
          gc_t <- if (!is.null(reg$gc_target) && !is.na(reg$gc_target[i])) reg$gc_target[i] else if (etype == "enhancer") enhancer_gc else silencer_gc
          seq_val <- make_random_genome(elen, gc_t)
        }

        rid <- paste0(etype, "_", i)
        gsym <- if (!is.null(reg$gene_symbol)) reg$gene_symbol[i] else paste0("gene_", tgi)

        reg_ids <- c(reg_ids, paste0(rid, "|", gsym))
        reg_seqs <- c(reg_seqs, seq_val)

        reg_table <- rbind(reg_table, data.frame(
          gene_symbol = gsym,
          target_gene_index = tgi,
          element_type = etype,
          offset_bp = offset,
          length_bp = elen,
          start = estart,
          end = eend,
          sequence_id = ifelse(is.na(seq_id), "", seq_id),
          stringsAsFactors = FALSE
        ))

        gff <- rbind(gff, data.frame(
          seqid = seqname, source = "simulator", type = etype,
          start = estart, end = eend, score = ".", strand = strands[tgi], phase = ".",
          attributes = paste0("ID=", rid, ";Parent=gene_", tgi, ";Name=", gsym)
        ))
      }

      if (!is.na(regulatory_fa_out) && length(reg_ids) > 0) {
        write_fasta(regulatory_fa_out, reg_ids, reg_seqs)
      }
    }
  }

  }

  # rRNA + tRNA clusters
  rrna_rows <- data.frame(id = character(), start = integer(), end = integer(), strand = character(), stringsAsFactors = FALSE)
  trna_rows <- data.frame(id = character(), start = integer(), end = integer(), strand = character(), stringsAsFactors = FALSE)

  if (rrna_clusters > 0) {
    cluster_start <- max(1, floor(glen / 2) - 5000)
    for (c in 1:rrna_clusters) {
      rstart <- cluster_start + (c - 1) * 5000
      rlen <- 1500
      rend <- min(glen, rstart + rlen - 1)
      rrna_rows <- rbind(rrna_rows, data.frame(id = paste0("rRNA_", c), start = rstart, end = rend, strand = "+"))
      for (t in 1:trna_per_cluster) {
        tstart <- rend + 50 + (t - 1) * 100
        tend <- min(glen, tstart + 80 - 1)
        trna_rows <- rbind(trna_rows, data.frame(id = paste0("tRNA_", c, "_", t), start = tstart, end = tend, strand = "+"))
      }
    }
  }

  for (j in seq_len(nrow(prom_rows))) {
    gff <- rbind(gff, data.frame(seqid = seqname, source = "simulator", type = "promoter",
                                 start = prom_rows$start[j], end = prom_rows$end[j], score = ".",
                                 strand = prom_rows$strand[j], phase = ".",
                                 attributes = paste0("ID=", prom_rows$id[j])))
  }

  for (j in seq_len(nrow(tss_rows))) {
    gff <- rbind(gff, data.frame(seqid = seqname, source = "simulator", type = "TSS",
                                 start = tss_rows$tss[j], end = tss_rows$tss[j], score = ".",
                                 strand = tss_rows$strand[j], phase = ".",
                                 attributes = paste0("ID=TSS_", j, ";Operon=", tss_rows$operon[j])))
  }

  # terminators after operons
  term_id <- 1
  for (op in unique(operon_id)) {
    if (runif(1) < terminator_prob) {
      last_idx <- tail(which(operon_id == op), 1)
      if (strands[last_idx] == "+") {
        tstart <- min(glen, ends[last_idx] + 20)
        tend <- min(glen, tstart + 30)
      } else {
        tend <- max(1, starts[last_idx] - 20)
        tstart <- max(1, tend - 30)
      }
      gff <- rbind(gff, data.frame(seqid = seqname, source = "simulator", type = "terminator",
                                   start = tstart, end = tend, score = ".", strand = strands[last_idx], phase = ".",
                                   attributes = paste0("ID=term_", term_id, ";Operon=", op)))
      term_id <- term_id + 1
    }
  }

  # riboswitches
  ribo_rows <- data.frame(id = character(), start = integer(), end = integer(), strand = character(), stringsAsFactors = FALSE)
  if (riboswitch_count > 0) {
    for (i in 1:riboswitch_count) {
      pos <- sample(1:glen, 1)
      rstart <- max(1, pos - 40)
      rend <- min(glen, pos + 40)
      ribo_rows <- rbind(ribo_rows, data.frame(id = paste0("riboswitch_", i), start = rstart, end = rend, strand = "."))
      gff <- rbind(gff, data.frame(seqid = seqname, source = "simulator", type = "riboswitch",
                                   start = rstart, end = rend, score = ".", strand = ".", phase = ".",
                                   attributes = paste0("ID=riboswitch_", i)))
    }
  }

  # CRISPR arrays
  crispr_rows <- data.frame(id = character(), start = integer(), end = integer(), strand = character(), stringsAsFactors = FALSE)
  if (crispr_count > 0) {
    for (i in 1:crispr_count) {
      pos <- sample(1:(glen - crispr_len), 1)
      cstart <- pos
      cend <- pos + crispr_len - 1
      crispr_rows <- rbind(crispr_rows, data.frame(id = paste0("CRISPR_", i), start = cstart, end = cend, strand = "+"))
      gff <- rbind(gff, data.frame(seqid = seqname, source = "simulator", type = "CRISPR",
                                   start = cstart, end = cend, score = ".", strand = "+", phase = ".",
                                   attributes = paste0("ID=CRISPR_", i)))
    }
  }

  # origin/terminus GC-skew signal
  if (gc_skew_signal) {
    ori <- max(1, floor(glen * 0.05))
    ter <- max(1, floor(glen * 0.55))
    gff <- rbind(gff, data.frame(seqid = seqname, source = "simulator", type = "origin_of_replication",
                                 start = ori, end = ori + 100, score = ".", strand = ".", phase = ".",
                                 attributes = "ID=oriC"))
    gff <- rbind(gff, data.frame(seqid = seqname, source = "simulator", type = "replication_terminus",
                                 start = ter, end = ter + 100, score = ".", strand = ".", phase = ".",
                                 attributes = "ID=ter"))
  }

  # rRNA/tRNA features
  if (nrow(rrna_rows) > 0) {
    for (j in seq_len(nrow(rrna_rows))) {
      gff <- rbind(gff, data.frame(seqid = seqname, source = "simulator", type = "rRNA",
                                   start = rrna_rows$start[j], end = rrna_rows$end[j], score = ".",
                                   strand = rrna_rows$strand[j], phase = ".",
                                   attributes = paste0("ID=", rrna_rows$id[j])))
    }
  }
  if (nrow(trna_rows) > 0) {
    for (j in seq_len(nrow(trna_rows))) {
      gff <- rbind(gff, data.frame(seqid = seqname, source = "simulator", type = "tRNA",
                                   start = trna_rows$start[j], end = trna_rows$end[j], score = ".",
                                   strand = trna_rows$strand[j], phase = ".",
                                   attributes = paste0("ID=", trna_rows$id[j])))
    }
  }

  # optional splice-like exon/intron features
  if (add_splicing) {
    for (i in 1:n_genes) {
      if (runif(1) < splice_prob) {
        n_exons <- sample(min_exons:max_exons, 1)
        gstart <- starts[i]
        gend <- ends[i]
        exon_len <- floor((gend - gstart + 1) / n_exons)
        for (e in 1:n_exons) {
          estart <- gstart + (e - 1) * exon_len
          eend <- if (e == n_exons) gend else estart + exon_len - 1
          gff <- rbind(gff, data.frame(seqid = seqname, source = "simulator", type = "exon",
                                       start = estart, end = eend, score = ".", strand = strands[i], phase = ".",
                                       attributes = paste0("Parent=gene_", i)))
          if (e < n_exons) {
            istart <- eend + 1
            iend <- min(gend, istart + sample(20:100, 1))
            gff <- rbind(gff, data.frame(seqid = seqname, source = "simulator", type = "intron",
                                         start = istart, end = iend, score = ".", strand = strands[i], phase = ".",
                                         attributes = paste0("Parent=gene_", i)))
          }
        }
      }
    }
  }

  dir.create(dirname(out_gff3), showWarnings = FALSE, recursive = TRUE)
  writeLines("##gff-version 3", out_gff3)
  write.table(gff, out_gff3, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)

  if (!is.na(genes_fa) || !is.na(cds_fa)) {
    ids <- paste0("gene_", seq_len(n_genes))
    seqs <- mapply(function(s, e, st) {
      sub <- substr(genome, s, e)
      if (st == "-") revcomp(sub) else sub
    }, starts, ends, strands, SIMPLIFY = TRUE)
    if (!is.na(genes_fa)) write_fasta(genes_fa, ids, seqs)
    if (!is.na(cds_fa)) write_fasta(cds_fa, paste0("cds_", seq_len(n_genes)), seqs)
  }

  if (!is.na(prom_fa)) {
    pseqs <- mapply(function(s, e, st) {
      sub <- substr(genome, s, e)
      if (st == "-") revcomp(sub) else sub
    }, prom_rows$start, prom_rows$end, prom_rows$strand, SIMPLIFY = TRUE)
    write_fasta(prom_fa, prom_rows$id, pseqs)
  }

  if (!is.na(rrna_fa) && nrow(rrna_rows) > 0) {
    rseqs <- mapply(function(s, e, st) {
      sub <- substr(genome, s, e)
      if (st == "-") revcomp(sub) else sub
    }, rrna_rows$start, rrna_rows$end, rrna_rows$strand, SIMPLIFY = TRUE)
    write_fasta(rrna_fa, rrna_rows$id, rseqs)
  }

  if (!is.na(trna_fa) && nrow(trna_rows) > 0) {
    tseqs <- mapply(function(s, e, st) {
      sub <- substr(genome, s, e)
      if (st == "-") revcomp(sub) else sub
    }, trna_rows$start, trna_rows$end, trna_rows$strand, SIMPLIFY = TRUE)
    write_fasta(trna_fa, trna_rows$id, tseqs)
  }

  if (!is.na(riboswitch_fa) && nrow(ribo_rows) > 0) {
    rseqs <- mapply(function(s, e, st) {
      sub <- substr(genome, s, e)
      if (st == "-") revcomp(sub) else sub
    }, ribo_rows$start, ribo_rows$end, ribo_rows$strand, SIMPLIFY = TRUE)
    write_fasta(riboswitch_fa, ribo_rows$id, rseqs)
  }

  if (!is.na(crispr_fa) && nrow(crispr_rows) > 0) {
    cseqs <- mapply(function(s, e, st) {
      sub <- substr(genome, s, e)
      if (st == "-") revcomp(sub) else sub
    }, crispr_rows$start, crispr_rows$end, crispr_rows$strand, SIMPLIFY = TRUE)
    write_fasta(crispr_fa, crispr_rows$id, cseqs)
  }

  if (!is.na(tss_tsv)) {
    dir.create(dirname(tss_tsv), showWarnings = FALSE, recursive = TRUE)
    write.table(tss_rows, tss_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  }

  if (!is.na(operons_tsv)) {
    dir.create(dirname(operons_tsv), showWarnings = FALSE, recursive = TRUE)
    op_df <- data.frame(gene_id = paste0("gene_", seq_len(n_genes)), operon = operon_id, stringsAsFactors = FALSE)
    write.table(op_df, operons_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  }

  # Plasmid simulation + GFF3 + markers
  if (plasmid_count > 0 && !is.na(plasmid_fa_out)) {
    plasmid_ids <- paste0("plasmid_", seq_len(plasmid_count))
    plasmid_seqs <- replicate(plasmid_count, make_random_genome(plasmid_length, plasmid_gc))

    marker_pool <- get_marker_pool()
    plasmid_gff <- data.frame(seqid = character(), source = character(), type = character(),
                              start = integer(), end = integer(), score = character(), strand = character(),
                              phase = character(), attributes = character(), stringsAsFactors = FALSE)

    for (i in seq_len(plasmid_count)) {
      inserted <- insert_markers(plasmid_seqs[i], marker_pool$ids, marker_pool$seqs)
      plasmid_seqs[i] <- inserted$seq
      if (nrow(inserted$feats) > 0) {
        for (j in seq_len(nrow(inserted$feats))) {
          plasmid_gff <- rbind(plasmid_gff, data.frame(
            seqid = plasmid_ids[i], source = "simulator", type = "marker",
            start = inserted$feats$start[j], end = inserted$feats$end[j],
            score = ".", strand = ".", phase = ".",
            attributes = paste0("ID=", inserted$feats$id[j])
          ))
        }
      }
    }

    write_fasta(plasmid_fa_out, plasmid_ids, plasmid_seqs)

    if (!is.na(plasmid_gff3_out)) {
      dir.create(dirname(plasmid_gff3_out), showWarnings = FALSE, recursive = TRUE)
      writeLines("##gff-version 3", plasmid_gff3_out)
      for (i in seq_len(plasmid_count)) {
        gline <- data.frame(seqid = plasmid_ids[i], source = "simulator", type = "plasmid",
                            start = 1, end = nchar(plasmid_seqs[i]), score = ".", strand = ".", phase = ".",
                            attributes = paste0("ID=", plasmid_ids[i]))
        write.table(gline, plasmid_gff3_out, sep = "\t", row.names = FALSE, col.names = FALSE,
                    quote = FALSE, append = TRUE)
      }
      if (nrow(plasmid_gff) > 0) {
        write.table(plasmid_gff, plasmid_gff3_out, sep = "\t", row.names = FALSE, col.names = FALSE,
                    quote = FALSE, append = TRUE)
      }
    }
  }

} else {
  ref <- read_fasta_one(fasta_path)
  seqname <- ref$name
  genome <- ref$seq
  gff <- parse_gff3(in_gff3)

  dir.create(dirname(out_gff3), showWarnings = FALSE, recursive = TRUE)
  writeLines("##gff-version 3", out_gff3)
  write.table(gff, out_gff3, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)

  genes <- gff[gff$type == "gene", ]
  cds <- gff[gff$type == "CDS", ]
  promoters <- gff[gff$type == "promoter", ]
  tss <- gff[gff$type == "TSS", ]

  if (!is.na(genes_fa) && nrow(genes) > 0) {
    ids <- paste0("gene_", seq_len(nrow(genes)))
    seqs <- mapply(function(s, e, st) {
      sub <- substr(genome, s, e)
      if (st == "-") revcomp(sub) else sub
    }, genes$start, genes$end, genes$strand, SIMPLIFY = TRUE)
    write_fasta(genes_fa, ids, seqs)
  }

  if (!is.na(cds_fa) && nrow(cds) > 0) {
    ids <- paste0("cds_", seq_len(nrow(cds)))
    seqs <- mapply(function(s, e, st) {
      sub <- substr(genome, s, e)
      if (st == "-") revcomp(sub) else sub
    }, cds$start, cds$end, cds$strand, SIMPLIFY = TRUE)
    write_fasta(cds_fa, ids, seqs)
  }

  if (!is.na(prom_fa) && nrow(promoters) > 0) {
    ids <- paste0("promoter_", seq_len(nrow(promoters)))
    seqs <- mapply(function(s, e, st) {
      sub <- substr(genome, s, e)
      if (st == "-") revcomp(sub) else sub
    }, promoters$start, promoters$end, promoters$strand, SIMPLIFY = TRUE)
    write_fasta(prom_fa, ids, seqs)
  }

  if (!is.na(tss_tsv) && nrow(tss) > 0) {
    dir.create(dirname(tss_tsv), showWarnings = FALSE, recursive = TRUE)
    out <- data.frame(tss = tss$start, strand = tss$strand, stringsAsFactors = FALSE)
    write.table(out, tss_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}
}
