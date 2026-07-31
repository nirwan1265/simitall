main_10_simulate_gwas_cohort <- function(args = commandArgs(trailingOnly = TRUE)) {

# Simulate a GWAS cohort from a reference genome.
# Outputs: VCF, genotype matrix TSV, phenotype TSV.

args <- args
usage <- function() {
  cat("Usage: 10_simulate_gwas_cohort.R --genome_fa <path> --out_prefix <path> [options]\n")
  cat("\nRequired:\n")
  cat("  --genome_fa <path>           reference FASTA (single-contig preferred)\n")
  cat("  --out_prefix <path>          output prefix (writes .vcf, .geno.tsv, .pheno.tsv)\n")
  cat("\nCohort:\n")
  cat("  --n_samples <int>            default 100\n")
  cat("  --ploidy <int>               default 2 (1 for haploid)\n")
  cat("\nVariants:\n")
  cat("  --snp_rate <float>           per-bp SNP rate (default 0.001)\n")
  cat("  --indel_rate <float>         per-bp indel rate (default 0.0001)\n")
  cat("  --indel_maxlen <int>         max indel length (default 3)\n")
  cat("  --af_beta1 <float>           allele-freq Beta alpha (default 0.8)\n")
  cat("  --af_beta2 <float>           allele-freq Beta beta (default 0.8)\n")
  cat("\nPhenotype:\n")
  cat("  --phenotype <none|binary|quantitative> default quantitative\n")
  cat("  --n_causal <int>             number of causal variants (default 20)\n")
  cat("  --effect_sd <float>          effect size SD (default 0.5)\n")
  cat("  --case_frac <float>          target case fraction for binary (default 0.5)\n")
  cat("\nOutputs:\n")
  cat("  --geno_tsv_out <path>         override genotype TSV path\n")
  cat("  --pheno_tsv_out <path>        override phenotype TSV path\n")
  cat("\nMisc:\n")
  cat("  --seed <int>                 default 1\n")
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
out_prefix <- get_arg("--out_prefix")
if (is.null(fasta_path) || is.null(out_prefix)) usage()

# Params
n_samples <- as.integer(get_arg("--n_samples", 100))
ploidy <- as.integer(get_arg("--ploidy", 2))

snp_rate <- as.numeric(get_arg("--snp_rate", 0.001))
indel_rate <- as.numeric(get_arg("--indel_rate", 0.0001))
indel_maxlen <- as.integer(get_arg("--indel_maxlen", 3))

af_beta1 <- as.numeric(get_arg("--af_beta1", 0.8))
af_beta2 <- as.numeric(get_arg("--af_beta2", 0.8))

ld_block_size <- as.integer(get_arg("--ld_block_size", 0))
ld_haplotypes <- as.integer(get_arg("--ld_haplotypes", 6))

recomb_map_out <- get_arg("--recomb_map_out", NA)
recomb_map_in <- get_arg("--recomb_map_in", NA)
recomb_rate_mean <- as.numeric(get_arg("--recomb_rate_mean", 1.0))
recomb_rate_sd <- as.numeric(get_arg("--recomb_rate_sd", 0.3))
recomb_hotspots <- as.integer(get_arg("--recomb_hotspots", 3))
recomb_hotspot_mult <- as.numeric(get_arg("--recomb_hotspot_mult", 5.0))

n_pops <- as.integer(get_arg("--n_pops", 1))
pop_sizes_raw <- get_arg("--pop_sizes", NA)
fst <- as.numeric(get_arg("--fst", 0.0))
pop_effect_shift <- as.numeric(get_arg("--pop_effect_shift", 0.0))

phenotype <- get_arg("--phenotype", "quantitative")
n_causal <- as.integer(get_arg("--n_causal", 20))
effect_sd <- as.numeric(get_arg("--effect_sd", 0.5))
case_frac <- as.numeric(get_arg("--case_frac", 0.5))

seed <- as.integer(get_arg("--seed", 1))
set.seed(seed)

if (ploidy < 1) stop("--ploidy must be >= 1")
if (!(phenotype %in% c("none", "binary", "quantitative"))) stop("--phenotype must be none|binary|quantitative")
if (snp_rate < 0 || indel_rate < 0) stop("variant rates must be >= 0")
if (indel_maxlen < 1) indel_maxlen <- 1
if (is.na(ld_block_size) || ld_block_size < 0) ld_block_size <- 0
if (is.na(ld_haplotypes) || ld_haplotypes < 2) ld_haplotypes <- 2
if (is.na(n_pops) || n_pops < 1) n_pops <- 1
if (fst < 0) fst <- 0
if (recomb_rate_mean < 0) recomb_rate_mean <- 0
if (recomb_rate_sd < 0) recomb_rate_sd <- 0
if (is.na(recomb_hotspots) || recomb_hotspots < 0) recomb_hotspots <- 0
if (recomb_hotspot_mult < 1) recomb_hotspot_mult <- 1

# Output paths
vcf_out <- paste0(out_prefix, ".vcf")
geno_out <- get_arg("--geno_tsv_out", paste0(out_prefix, ".geno.tsv"))
pheno_out <- get_arg("--pheno_tsv_out", paste0(out_prefix, ".pheno.tsv"))

# FASTA reader (single contig)
read_fasta_one <- function(path) {
  lines <- readLines(path)
  hdr_idx <- grep("^>", lines)
  if (length(hdr_idx) == 0) stop("No FASTA header found")
  name <- sub("^>", "", strsplit(lines[hdr_idx[1]], " ")[[1]][1])
  seq <- paste(lines[(hdr_idx[1] + 1):(if (length(hdr_idx) > 1) hdr_idx[2] - 1 else length(lines))], collapse = "")
  list(name = name, seq = toupper(seq))
}

ref <- read_fasta_one(fasta_path)
seqname <- ref$name
seq <- ref$seq
glen <- nchar(seq)

# Build synthetic recombination map (TSV: pos_bp, cM)
make_recomb_map <- function(glen, rate_mean, rate_sd, n_hotspots, hotspot_mult) {
  base_rate <- rate_mean / 1e6
  step <- 50000
  positions <- seq(1, glen, by = step)
  if (tail(positions, 1) != glen) positions <- c(positions, glen)

  rates <- pmax(0, rnorm(length(positions), mean = base_rate, sd = rate_sd / 1e6))

  if (n_hotspots > 0) {
    hs_pos <- sample(positions, min(n_hotspots, length(positions)), replace = FALSE)
    for (i in seq_along(positions)) {
      if (positions[i] %in% hs_pos) rates[i] <- rates[i] * hotspot_mult
    }
  }

  cM <- numeric(length(positions))
  for (i in 2:length(positions)) {
    dist <- positions[i] - positions[i-1]
    cM[i] <- cM[i-1] + rates[i] * dist
  }

  data.frame(pos_bp = positions, cM = cM, stringsAsFactors = FALSE)
}

read_recomb_map <- function(path, glen) {
  m <- read.delim(path, header = TRUE, sep = "	", stringsAsFactors = FALSE)
  if (!("pos_bp" %in% names(m)) || !("cM" %in% names(m))) stop("recomb map must have pos_bp and cM columns")
  m <- m[order(m$pos_bp), ]
  if (m$pos_bp[1] != 1) {
    m <- rbind(data.frame(pos_bp = 1, cM = 0), m)
  }
  if (tail(m$pos_bp, 1) < glen) {
    m <- rbind(m, data.frame(pos_bp = glen, cM = tail(m$cM, 1)))
  }
  m
}

# Recombination map
if (!is.na(recomb_map_in) && file.exists(recomb_map_in)) {
  recomb_map <- read_recomb_map(recomb_map_in, glen)
} else {
  recomb_map <- make_recomb_map(
    glen,
    recomb_rate_mean,
    recomb_rate_sd,
    recomb_hotspots,
    recomb_hotspot_mult
  )
  if (!is.na(recomb_map_out)) {
    dir.create(dirname(recomb_map_out), showWarnings = FALSE, recursive = TRUE)
    write.table(
      recomb_map,
      recomb_map_out,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
}

parse_pop_sizes <- function(n_samples, n_pops, pop_sizes_raw) {
  if (is.na(pop_sizes_raw)) {
    base <- floor(n_samples / n_pops)
    sizes <- rep(base, n_pops)
    sizes[1] <- sizes[1] + (n_samples - sum(sizes))
    return(sizes)
  }
  parts <- as.numeric(strsplit(pop_sizes_raw, ",")[[1]])
  if (length(parts) != n_pops) stop("--pop_sizes must have n_pops values")
  if (all(parts <= 1)) {
    sizes <- round(parts / sum(parts) * n_samples)
    sizes[1] <- sizes[1] + (n_samples - sum(sizes))
  } else {
    sizes <- as.integer(parts)
    if (sum(sizes) != n_samples) stop("--pop_sizes counts must sum to n_samples")
  }
  sizes
}

# Determine variant counts
n_snps <- as.integer(floor(snp_rate * glen))
n_indels <- as.integer(floor(indel_rate * glen))

# Helper: choose non-overlapping indel positions
occupied <- rep(FALSE, glen)

pick_indel_pos <- function(len) {
  tries <- 0
  while (tries < 1000) {
    pos <- sample.int(glen - len + 1, 1)
    if (!any(occupied[pos:(pos+len-1)])) {
      occupied[pos:(pos+len-1)] <<- TRUE
      return(pos)
    }
    tries <- tries + 1
  }
  NA
}

# SNP positions (avoid occupied)
if (n_snps > 0) {
  free_pos <- which(!occupied)
  if (length(free_pos) < n_snps) n_snps <- length(free_pos)
  snp_pos <- sort(sample(free_pos, n_snps, replace = FALSE))
  occupied[snp_pos] <- TRUE
} else {
  snp_pos <- integer()
}

# Indels
indel_pos <- integer()
indel_len <- integer()
indel_type <- character() # INS or DEL
if (n_indels > 0) {
  for (i in seq_len(n_indels)) {
    len <- sample(1:indel_maxlen, 1)
    pos <- pick_indel_pos(len)
    if (is.na(pos)) next
    indel_pos <- c(indel_pos, pos)
    indel_len <- c(indel_len, len)
    indel_type <- c(indel_type, ifelse(runif(1) < 0.5, "INS", "DEL"))
  }
}

# Build variant table
variants <- data.frame(
  pos = c(snp_pos, indel_pos),
  type = c(rep("SNP", length(snp_pos)), indel_type),
  len = c(rep(1, length(snp_pos)), indel_len),
  stringsAsFactors = FALSE
)
if (nrow(variants) == 0) stop("No variants generated; increase rates or genome length")
variants <- variants[order(variants$pos), ]

# Generate REF/ALT
bases <- c("A", "C", "G", "T")
ref_allele <- character(nrow(variants))
alt_allele <- character(nrow(variants))

for (i in seq_len(nrow(variants))) {
  pos <- variants$pos[i]
  if (variants$type[i] == "SNP") {
    refb <- substr(seq, pos, pos)
    altb <- sample(bases[bases != refb], 1)
    ref_allele[i] <- refb
    alt_allele[i] <- altb
  } else {
    if (variants$type[i] == "DEL") {
      refb <- substr(seq, pos, pos + variants$len[i])
      altb <- substr(seq, pos, pos)
      ref_allele[i] <- refb
      alt_allele[i] <- altb
    } else {
      refb <- substr(seq, pos, pos)
      ins <- paste(sample(bases, variants$len[i], replace = TRUE), collapse = "")
      ref_allele[i] <- refb
      alt_allele[i] <- paste0(refb, ins)
    }
  }
}

variants$ref <- ref_allele
variants$alt <- alt_allele
variants$id <- paste0("var", seq_len(nrow(variants)))

# Allele frequencies and genotypes
p <- rbeta(nrow(variants), af_beta1, af_beta2)

# Population-specific allele frequencies (Fst drift)
# Beta approximation: p_k ~ Beta(p*(1-fst)/fst, (1-p)*(1-fst)/fst)
if (fst > 0 && n_pops > 1) {
  alpha <- p * (1 - fst) / fst
  beta <- (1 - p) * (1 - fst) / fst
  p_pop <- sapply(seq_len(n_pops), function(k) rbeta(length(p), alpha, beta))
} else {
  p_pop <- matrix(p, nrow = length(p), ncol = n_pops)
}

# Genotype matrix: rows = variants, cols = samples
samples <- paste0("sample_", seq_len(n_samples))
pop_sizes <- parse_pop_sizes(n_samples, n_pops, pop_sizes_raw)
pop_ids <- rep(paste0("pop", seq_len(n_pops)), times = pop_sizes)
if (length(pop_ids) != n_samples) stop("population sizes do not match n_samples")
G <- matrix(0, nrow = nrow(variants), ncol = n_samples)
pop_index <- match(pop_ids, paste0("pop", seq_len(n_pops)))
if (ld_block_size > 0) {
  block_id <- floor((variants$pos - 1) / ld_block_size) + 1
  for (b in sort(unique(block_id))) {
    idx <- which(block_id == b)
    # haplotypes: rows=variants in block, cols=haplotypes
    H <- matrix(0, nrow = length(idx), ncol = ld_haplotypes)
    for (j in seq_along(idx)) {
      H[j, ] <- rbinom(ld_haplotypes, 1, p_pop[idx[j], pop_idx])
    }
    # assign haplotypes per sample and sum per ploidy
    for (s in seq_len(n_samples)) {
      pop_idx <- pop_index[s]
      # re-simulate haplotypes per population to reflect drift
      H <- matrix(0, nrow = length(idx), ncol = ld_haplotypes)
      for (j in seq_along(idx)) {
        H[j, ] <- rbinom(ld_haplotypes, 1, p_pop[idx[j], pop_idx])
      }
      h_idx <- sample.int(ld_haplotypes, ploidy, replace = TRUE)
      G[idx, s] <- rowSums(H[, h_idx, drop = FALSE])
    }
  }
} else {
  for (i in seq_len(nrow(variants))) {
    G[i, ] <- sapply(seq_len(n_samples), function(s) rbinom(1, ploidy, p_pop[i, pop_index[s]]))
  }
}

# Phenotype
pheno <- data.frame(sample = samples, pop = pop_ids, trait = NA, stringsAsFactors = FALSE)
if (phenotype != "none") {
  n_causal <- min(n_causal, nrow(variants))
  causal_idx <- sample(seq_len(nrow(variants)), n_causal, replace = FALSE)
  effects <- rnorm(n_causal, mean = 0, sd = effect_sd)
  genetic_score <- as.numeric(t(effects) %*% G[causal_idx, , drop = FALSE])

  if (phenotype == "quantitative") {
    pop_shift <- (pop_index - 1) * pop_effect_shift
      pheno$trait <- genetic_score + pop_shift + rnorm(n_samples, 0, 1)
  } else if (phenotype == "binary") {
    # logistic model, adjust intercept to get target case fraction
    f <- function(intercept) {
      pop_shift <- (pop_index - 1) * pop_effect_shift
      pcase <- 1 / (1 + exp(-(intercept + genetic_score + pop_shift)))
      mean(pcase) - case_frac
    }
    lo <- -10; hi <- 10
    for (k in 1:30) {
      mid <- (lo + hi)/2
      if (f(mid) > 0) hi <- mid else lo <- mid
    }
    intercept <- (lo + hi)/2
    pop_shift <- (pop_index - 1) * pop_effect_shift
      pcase <- 1 / (1 + exp(-(intercept + genetic_score + pop_shift)))
    pheno$trait <- rbinom(n_samples, 1, pcase)
  }
}

# Write VCF
writeLines("##fileformat=VCFv4.2", vcf_out)
.simitall_append_lines(
  paste0("##contig=<ID=", seqname, ",length=", glen, ">"),
  vcf_out
)
.simitall_append_lines(
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
  vcf_out
)
header <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", samples)
# add per-sample population labels as additional header lines
for (i in seq_along(samples)) {
  .simitall_append_lines(
    paste0("##SAMPLE=<ID=", samples[i], ",POP=", pop_ids[i], ">"),
    vcf_out
  )
}
.simitall_append_lines(paste(header, collapse = "\t"), vcf_out)

for (i in seq_len(nrow(variants))) {
  gt <- character(n_samples)
  if (ploidy == 1) {
    gt <- ifelse(G[i, ] == 0, "0", "1")
  } else {
    gt <- ifelse(G[i, ] == 0, "0/0",
                 ifelse(G[i, ] == 1, "0/1", "1/1"))
  }
  row <- c(seqname, variants$pos[i], variants$id[i], variants$ref[i], variants$alt[i], ".", "PASS", ".", "GT", gt)
  .simitall_append_lines(paste(row, collapse = "\t"), vcf_out)
}

# Write genotype TSV (rows=variants, cols=samples)
geno_df <- data.frame(id = variants$id, pos = variants$pos, ref = variants$ref, alt = variants$alt, G, stringsAsFactors = FALSE)
colnames(geno_df)[5:ncol(geno_df)] <- samples
write.table(geno_df, geno_out, sep = "\t", row.names = FALSE, quote = FALSE)

# Write phenotype TSV
if (phenotype != "none") {
  write.table(pheno, pheno_out, sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("Wrote:\n")
cat("  VCF:   ", vcf_out, "\n", sep = "")
cat("  Geno:  ", geno_out, "\n", sep = "")
if (phenotype != "none") cat("  Pheno: ", pheno_out, "\n", sep = "")
}
