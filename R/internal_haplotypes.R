main_11a_generate_random_haplotype_panel <- function(args = commandArgs(trailingOnly = TRUE)) {

# Generate a random haplotype panel (multi-FASTA)
# Each haplotype is a mutated copy of a base sequence

args <- args
usage <- function() {
  cat("Usage: 11a_generate_random_haplotype_panel.R --out_fa <path> [options]\n")
  cat("\nOptions:\n")
  cat("  --n_haps <int>            default 8\n")
  cat("  --length <int>            default 100000\n")
  cat("  --gc <float>              default 0.5\n")
  cat("  --snp_rate <float>        default 0.001\n")
  cat("  --indel_rate <float>      default 0.0001\n")
  cat("  --indel_maxlen <int>      default 3\n")
  cat("  --seed <int>              default 1\n")
  quit(status = 1)
}

get_arg <- function(flag, default = NULL) {
  if (!(flag %in% args)) return(default)
  idx <- match(flag, args)
  if (idx == length(args)) return(default)
  args[idx + 1]
}

out_fa <- get_arg("--out_fa")
if (is.null(out_fa)) usage()

n_haps <- as.integer(get_arg("--n_haps", 8))
len <- as.integer(get_arg("--length", 100000))
gc <- as.numeric(get_arg("--gc", 0.5))
snp_rate <- as.numeric(get_arg("--snp_rate", 0.001))
indel_rate <- as.numeric(get_arg("--indel_rate", 0.0001))
indel_maxlen <- as.integer(get_arg("--indel_maxlen", 3))
seed <- as.integer(get_arg("--seed", 1))
set.seed(seed)

make_random_genome <- function(length, gc_target) {
  bases <- c("A", "C", "G", "T")
  probs <- c((1 - gc_target)/2, gc_target/2, gc_target/2, (1 - gc_target)/2)
  paste(sample(bases, length, replace = TRUE, prob = probs), collapse = "")
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

base <- make_random_genome(len, gc)
ids <- paste0("hap", seq_len(n_haps))
seqs <- character(n_haps)
seqs[1] <- base
if (n_haps > 1) {
  for (i in 2:n_haps) {
    seqs[i] <- mutate_sequence(base, snp_rate, indel_rate, indel_maxlen)
  }
}

write_fasta(out_fa, ids, seqs)
cat("Wrote panel:", out_fa, "\n")
}
