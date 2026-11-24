# usage: Rscript read.R <lib or ""> <out.tsv>
a <- commandArgs(TRUE); if (nzchar(a[1])) .libPaths(c(a[1], .libPaths()))
suppressMessages(library(studyforge))
sp <- "/home/claude/BioDesignGPT/training/paper_questions/splits"
qs <- unlist(lapply(c("test","dev","dev2","test_draw4","test_draw5","test_heldout","test_newkinds","train"), function(f)
  lapply(readLines(file.path(sp, paste0(f, ".jsonl"))), jsonlite::fromJSON, simplifyVector = FALSE)), recursive = FALSE)
qs <- qs[!duplicated(vapply(qs, `[[`, "", "id"))]
fmt <- function(s) paste(vapply(s$factors %||% list(), function(f) sprintf("%s(%s)", f$name %||% "?", if (is.null(f$n_levels)) "?" else f$n_levels), ""), collapse = " | ")
`%||%` <- function(x, y) if (is.null(x)) y else x
t0 <- Sys.time()
out <- do.call(rbind, lapply(qs, function(q) {
  s <- tryCatch(suppressMessages(suppressWarnings(studyforge::extract_spec_local(q$prompt))), error = function(e) NULL)
  data.frame(id = q$id, factors = if (is.null(s)) "ERROR" else fmt(s),
             open = if (is.null(s)) NA else sum(vapply(s$factors %||% list(), function(f) is.null(f$n_levels), TRUE)), stringsAsFactors = FALSE)
}))
write.table(out, a[2], sep = "\t", row.names = FALSE, quote = FALSE)
cat(nrow(out), "prompts,", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
