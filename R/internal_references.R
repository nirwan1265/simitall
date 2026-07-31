main_01d_fetch_ref_files <- function(args = commandArgs(trailingOnly = TRUE)) {

# Download reference FASTA files using NCBI accessions

args <- args
outdir <- ifelse(length(args) >= 1, args[1], "inst/extdata/ref_files")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

fetch_ncbi <- function(acc, out) {
  base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  url <- paste0(base, "?db=nuccore&id=", acc, "&rettype=fasta&retmode=text")
  fa <- tryCatch(readLines(url), error = function(e) NULL)
  if (is.null(fa)) stop("Failed to fetch: ", acc)
  writeLines(fa, out)
}

fetch_ncbi("NC_000913.3", file.path(outdir, "ecoli_k12_mg1655_chr.fa"))
fetch_ncbi("NC_000001.11", file.path(outdir, "human_grch38_chr1.fa"))
fetch_ncbi("NC_050096.1", file.path(outdir, "maize_b73_refseq_chr1.fa"))
fetch_ncbi("NC_003070.9", file.path(outdir, "arabidopsis_thaliana_chr1.fa"))

cat("Wrote reference FASTA files to ", outdir, "\n", sep = "")
}
