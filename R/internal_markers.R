main_01c_fetch_marker_panel <- function(args = commandArgs(trailingOnly = TRUE)) {

# Fetch curated marker panel FASTA from NCBI using TSV list

args <- args
ts_path <- ifelse(length(args) >= 1, args[1], "data/markers_curated.tsv")
out_fa <- ifelse(length(args) >= 2, args[2], "data/markers_curated.fa")

if (!file.exists(ts_path)) stop("TSV not found: ", ts_path)

dir.create(dirname(out_fa), showWarnings = FALSE, recursive = TRUE)

fetch_one <- function(name, acc, start, end, strand) {
  base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  url <- paste0(base, "?db=nuccore&id=", acc, "&rettype=fasta&retmode=text")
  if (!is.na(start) && !is.na(end) && start != "" && end != "") {
    url <- paste0(url, "&seq_start=", start, "&seq_stop=", end)
  }
  fa <- tryCatch(paste(readLines(url), collapse = "\n"), error = function(e) "")
  if (fa == "") return(NA)
  hdr <- paste0(">", name, "|", acc,
               ifelse(!is.na(start) && start != "" && !is.na(end) && end != "", paste0("|", start, "-", end), ""),
               ifelse(!is.na(strand) && strand != "", paste0("|strand=", strand), ""))
  seq <- paste(grep("^[^>]", strsplit(fa, "\n")[[1]], value = TRUE), collapse = "")
  if (strand == "-") {
    comp <- chartr("ACGTacgt", "TGCAtgca", seq)
    seq <- paste(rev(strsplit(comp, "")[[1]]), collapse = "")
  }
  paste(hdr, seq, sep = "\n")
}

x <- read.delim(ts_path, sep = "\t", stringsAsFactors = FALSE)
fa_entries <- c()
for (i in seq_len(nrow(x))) {
  fa <- fetch_one(x$name[i], x$accession[i], x$start[i], x$end[i], x$strand[i])
  if (!is.na(fa)) fa_entries <- c(fa_entries, fa, "")
}
writeLines(fa_entries, out_fa)
cat("Wrote curated markers: ", out_fa, "\n", sep = "")
}
