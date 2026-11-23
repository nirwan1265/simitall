# Build the simulation templates in inst/extdata/templates/ from data-raw/templates_incoming/.
# Each template: counts.tsv.gz (features x samples, raw integer counts), samples.tsv (design),
# optional features.tsv.gz (annotation). Registered in inst/extdata/templates/templates.tsv.
# Run from the package root: Rscript data-raw/templates/make_templates.R
inc <- "data-raw/templates_incoming"
out <- "inst/extdata/templates"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rd <- function(f, ...) utils::read.csv(file.path(inc, f), check.names = FALSE, stringsAsFactors = FALSE, ...)
put <- function(id, counts, samples, features = NULL) {
  d <- file.path(out, id); dir.create(d, showWarnings = FALSE)
  stopifnot(identical(colnames(counts), samples$sample), all(counts == round(counts)), all(counts >= 0))
  storage.mode(counts) <- "integer"
  z <- gzfile(file.path(d, "counts.tsv.gz"), "w")
  utils::write.table(data.frame(feature_id = rownames(counts), counts, check.names = FALSE), z, sep = "\t", row.names = FALSE, quote = FALSE); close(z)
  utils::write.table(samples, file.path(d, "samples.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  if (!is.null(features)) {
    z <- gzfile(file.path(d, "features.tsv.gz"), "w")
    utils::write.table(features, z, sep = "\t", row.names = FALSE, quote = FALSE); close(z)
  }
  cat(id, ":", nrow(counts), "features x", ncol(counts), "samples\n")
}
keep_detected <- function(m, min_samples) m[rowSums(m > 0) >= min_samples, , drop = FALSE]

# human bulk RNA-seq: airway (4 airway smooth-muscle cell lines x dexamethasone)
a <- as.matrix(rd("airway_counts.csv", row.names = 1)); s <- rd("airway_samples.csv")
a <- keep_detected(a, 2)
put("rnaseq_human_airway", a,
    data.frame(sample = s[[1]], group = s$dex, cell_line = s$cell, run = s$Run))

# human / mouse miRNA: miRNA Tissue Atlas 2025, counts recovered exactly from isomiR RPMM
for (sp in c(hsa = "human", mmu = "mouse")) {
  k <- names(which(c(hsa = "human", mmu = "mouse") == sp))
  m <- as.matrix(rd(paste0(k, "_mirna_counts_3tissues.csv"), row.names = 1)); s <- rd(paste0(k, "_mirna_samples.csv"))
  put(paste0("small_rnaseq_", sp, "_mirna_atlas"), m,
      data.frame(sample = s$sample, group = s$group, donor = s$donor, library_size = s$library_size))
}

# plant 16S (soybean, Mendeley doi:10.17632/7y4gbwk44k.1, CC BY 4.0): soil / rhizosphere / root x 5 time points x 9 reps, rarefied to 14,510 reads
o <- rd("plant_16s_otu_counts.csv", row.names = 1); g <- rd("group.csv")
o <- t(as.matrix(o)); g <- g[match(colnames(o), g[[1]]), ]
o <- keep_detected(o, 5)
put("microbiome_plant_root_16s", o,
    data.frame(sample = colnames(o), group = g$zone, time = g$times, week = g$weeks, plot_group = g$group))

# human 16S: HMP V35 (V3-V5). Gut: 60 stool samples; body sites: the first 60 V35 samples. Taxonomy as features.
tx_all <- rd("hmp16s_taxonomy.csv"); names(tx_all)[1] <- "feature_id"
hmp <- function(id, counts_file, samples_file, min_samples) {
  h <- as.matrix(rd(counts_file, row.names = 1)); s <- rd(samples_file)
  h <- keep_detected(h, min_samples)
  s <- s[match(colnames(h), s[[1]]), ]
  tx <- tx_all[match(rownames(h), tx_all$feature_id), ]
  put(id, h,
      data.frame(sample = colnames(h), group = gsub(" ", "_", tolower(s$HMP_BODY_SUBSITE)), body_site = s$HMP_BODY_SITE, sex = s$SEX),
      tx[, c("feature_id", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS")])
}
hmp("microbiome_human_gut_hmp_16s", "hmp16s_stool_counts.csv", "hmp16s_stool_samples.csv", 3)
hmp("microbiome_human_body_sites_hmp_16s", "hmp16s_counts.csv", "hmp16s_samples.csv", 3)

# plant miRNA: PmiREN 3.1 (Zenodo), raw counts (miRNA_annotation/<sp>-mature.count); tissue of each library from
# miRNA_info/<sp>_mature.lib.exp.tsv (library names) and _mature.tissue.exp.tsv (labels, same column order).
# Libraries are matched by name; when the count file numbers its columns, by expression values, accepted only when
# the match is one-to-one with r >= 0.9. Species without a reliable match are left out (rice, potato, soybean).
pm <- file.path(inc, "pmiren")
norm_tissue <- function(x) {
  x <- tolower(trimws(gsub("\\.[0-9]+$", "", x)))
  x <- sub("^leaves$", "leaf", x); x <- sub("^seedlings$", "seedling", x); x <- sub("^seeds$", "seed", x); x <- sub("^roots$", "root", x)
  gsub("[^a-z0-9]+", "_", x)
}
zc <- function(a) { a <- sweep(a, 2, colMeans(a)); sweep(a, 2, pmax(sqrt(colSums(a^2)), 1e-9), "/") }
pmiren <- function(sp, organism) {
  cnt <- utils::read.delim(file.path(pm, paste0(sp, "-mature.count")), check.names = FALSE, stringsAsFactors = FALSE)
  tis <- utils::read.delim(file.path(pm, paste0(sp, "_mature.tissue.exp.tsv")), check.names = FALSE, stringsAsFactors = FALSE)
  lib <- strsplit(readLines(file.path(pm, paste0(sp, "_lib_names.tsv")), n = 1), "\t")[[1]][-1]
  labels <- names(tis)[-1]
  cols <- names(cnt)[-(1:2)]
  key <- ifelse(grepl("__", cols), sub(".*__", "", cols), cols)
  lab <- stats::setNames(labels[match(key, lib)], cols)
  how <- "library names"
  if (mean(!is.na(lab)) < 0.5) {
    k1 <- tolower(cnt$MIR); k2 <- tolower(tis$MIR)
    com <- intersect(k1[!duplicated(k1)], k2[!duplicated(k2)])
    C <- log1p(as.matrix(cnt[match(com, k1), cols])); T <- log1p(as.matrix(tis[match(com, k2), -1]))
    R <- crossprod(zc(C), zc(T)); best <- apply(R, 1, which.max); bv <- apply(R, 1, max)
    if (anyDuplicated(best) || min(bv) < 0.9) { message(sp, ": no reliable library -> tissue match; skipped"); return(invisible(NULL)) }
    lab <- stats::setNames(labels[best], cols); how <- sprintf("expression values (one-to-one, min r %.2f)", min(bv))
  }
  g <- norm_tissue(lab)
  ok <- !is.na(lab) & !grepl("^(srr|drr|err|gsm)[0-9]", g)
  m <- rowsum(as.matrix(cnt[, cols[ok]]), cnt$MIR)               # precursor rows of one miRNA name summed
  g <- g[ok]
  big <- names(which(table(g) >= 3)); keep <- g %in% big
  m <- keep_detected(m[, keep, drop = FALSE], 3); g <- g[keep]
  colnames(m) <- make.unique(colnames(m))
  dg <- if (sum(g == "leaf") >= 6) "leaf" else names(sort(table(g), decreasing = TRUE))[1]
  id <- paste0("small_rnaseq_", tolower(sub("_.*", "", sp)), "_", tolower(sub(".*_", "", sp)), "_pmiren")
  put(id, m, data.frame(sample = colnames(m), group = g, library = colnames(m)))
  data.frame(id = id, organism = organism, default_group = dg, how = how, n = ncol(m), groups = paste(names(sort(table(g), decreasing = TRUE))[1:min(5, length(unique(g)))], collapse = ", "))
}
pm_rows <- do.call(rbind, list(
  pmiren("Arabidopsis_thaliana", "Arabidopsis thaliana"), pmiren("Zea_mays", "Zea mays"), pmiren("Sorghum_bicolor", "Sorghum bicolor"),
  pmiren("Hordeum_vulgare", "Hordeum vulgare"), pmiren("Solanum_lycopersicum", "Solanum lycopersicum"),
  pmiren("Oryza_sativa", "Oryza sativa"), pmiren("Solanum_tuberosum", "Solanum tuberosum"), pmiren("Glycine_max", "Glycine max")))
print(pm_rows)

# mouse WGBS: ENCODE "methylation state at CpG" (bedMethyl), chr19:10-15 Mb, CpGs covered in all samples.
# methylated reads = round(coverage x percent / 100). Heart: ENCFF483COZ, ENCFF308OBC, ENCFF774SHU; liver: ENCFF157DQJ, ENCFF567QQL, ENCFF911IVZ.
wg <- file.path(inc, "encode_wgbs")
acc <- c(ENCFF483COZ = "heart", ENCFF308OBC = "heart", ENCFF774SHU = "heart", ENCFF157DQJ = "liver", ENCFF567QQL = "liver", ENCFF911IVZ = "liver")
w <- lapply(names(acc), function(a) {
  x <- utils::read.delim(file.path(wg, paste0(a, ".chr19_10-15Mb.tsv")), header = FALSE, col.names = c("chr", "pos", "strand", "cov", "pct"))
  x$key <- paste(x$chr, x$pos, x$strand); x
})
k <- Reduce(intersect, lapply(w, `[[`, "key"))
cov <- sapply(w, function(x) x$cov[match(k, x$key)]); pct <- sapply(w, function(x) x$pct[match(k, x$key)])
met <- round(cov * pct / 100); colnames(cov) <- colnames(met) <- names(acc)
sites <- do.call(rbind, strsplit(k, " ")); o <- order(as.integer(sites[, 2]))
sites <- data.frame(feature_id = paste0(sites[o, 1], ":", sites[o, 2], ":", sites[o, 3]), chr = sites[o, 1], pos = as.integer(sites[o, 2]), strand = sites[o, 3], context = "CG")
cov <- cov[o, ]; met <- met[o, ]; rownames(cov) <- rownames(met) <- sites$feature_id
d <- file.path(out, "methylation_mouse_encode_wgbs"); dir.create(d, showWarnings = FALSE)
wz <- function(m, f) { z <- gzfile(file.path(d, f), "w"); utils::write.table(data.frame(feature_id = rownames(m), m, check.names = FALSE), z, sep = "\t", row.names = FALSE, quote = FALSE); close(z) }
storage.mode(cov) <- storage.mode(met) <- "integer"
wz(cov, "counts.tsv.gz"); wz(met, "methylated.tsv.gz")
z <- gzfile(file.path(d, "features.tsv.gz"), "w"); utils::write.table(sites, z, sep = "\t", row.names = FALSE, quote = FALSE); close(z)
utils::write.table(data.frame(sample = names(acc), group = unname(acc), accession = names(acc)), file.path(d, "samples.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("methylation_mouse_encode_wgbs :", nrow(cov), "CpGs x", ncol(cov), "samples\n")

# human methylation array: TCGA TGCT HumanMethylation450 beta values (UCSC Xena), every 80th probe, 20 tumours.
b <- utils::read.delim(file.path(inc, "tcga_tgct_450k_slice.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
bm <- as.matrix(b[, -1]); rownames(bm) <- b[[1]]; bm <- bm[stats::complete.cases(bm), ]
d <- file.path(out, "methylation_array_human_tcga_tgct"); dir.create(d, showWarnings = FALSE)
z <- gzfile(file.path(d, "beta.tsv.gz"), "w"); utils::write.table(data.frame(feature_id = rownames(bm), round(bm, 4), check.names = FALSE), z, sep = "\t", row.names = FALSE, quote = FALSE); close(z)
utils::write.table(data.frame(sample = colnames(bm), group = ifelse(grepl("-05$", colnames(bm)), "additional_primary", "primary_tumor")), file.path(d, "samples.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("methylation_array_human_tcga_tgct :", nrow(bm), "probes x", ncol(bm), "samples (beta)\n")
