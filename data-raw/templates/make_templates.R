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

# plant 16S: soil / rhizosphere / root x 5 time points x 9 reps, rarefied to 14,510 reads
o <- rd("plant_16s_otu_counts.csv", row.names = 1); g <- rd("group.csv")
o <- t(as.matrix(o)); g <- g[match(colnames(o), g[[1]]), ]
o <- keep_detected(o, 5)
put("microbiome_plant_root_16s", o,
    data.frame(sample = colnames(o), group = g$zone, time = g$times, week = g$weeks, plot_group = g$group))

# human 16S: HMP (V3-V5), body sites; taxonomy as features
h <- as.matrix(rd("hmp16s_counts.csv", row.names = 1)); s <- rd("hmp16s_samples.csv"); tx <- rd("hmp16s_taxonomy.csv")
h <- keep_detected(h, 3)
s <- s[match(colnames(h), s[[1]]), ]
tx <- tx[match(rownames(h), tx[[1]]), ]; names(tx)[1] <- "feature_id"
put("microbiome_human_hmp_16s", h,
    data.frame(sample = colnames(h), group = gsub(" ", "_", tolower(s$HMP_BODY_SUBSITE)), body_site = s$HMP_BODY_SITE, sex = s$SEX),
    tx[, c("feature_id", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS")])
