# Biodes repeated-measures reading (.sf_repeated_on_unit) on every reviewed spec. usage: Rscript scan_repeated_factors.R out.tsv
`%||%` <- function(x, y) if (is.null(x)) y else x
suppressMessages(library(studyforge)); sf <- asNamespace("studyforge")
sp <- "/home/claude/BioDesignGPT/training/paper_questions/splits"
# test_draw5 is Biodes' locked test split: never read. Other test splits: report counts only, no ids or prompts.
qs <- unlist(lapply(c("test","dev","dev2","test_draw4","test_heldout","test_newkinds","train"), function(f)
  lapply(readLines(file.path(sp, paste0(f, ".jsonl"))), function(l) { q <- jsonlite::fromJSON(l, simplifyVector = FALSE); q$split <- f; q })), recursive = FALSE)
qs <- qs[!duplicated(vapply(qs, `[[`, "", "id"))]
out <- do.call(rbind, lapply(qs, function(q) {
  s <- q$expected_spec; txt <- sf$.sf_words_to_digits(tolower(q$prompt))
  s2 <- tryCatch(sf$.sf_repeated_on_unit(s, txt), error = function(e) { s$repeated_factors <- "ERROR"; s })
  occ <- intersect(vapply(s$factors %||% list(), function(f) f$name %||% "", ""), c("time", "stage", "age"))
  data.frame(id = q$id, split = q$split, occasion = paste(occ, collapse = ","), unit = s$experimental_unit %||% "",
             rf = if (is.null(s2$repeated_factors)) "(absent)" else paste0("[", paste(unlist(s2$repeated_factors), collapse = ","), "]"),
             factors = paste(vapply(s$factors %||% list(), function(f) sprintf("%s(%s)", f$name %||% "?", f$n_levels %||% "?"), ""), collapse = " | "),
             prompt = gsub("[\t\n]", " ", q$prompt), stringsAsFactors = FALSE)
}))
write.table(out, commandArgs(TRUE)[1] %||% "repeated_scan.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
print(table(out$rf)); cat("with an occasion factor:", sum(nzchar(out$occasion)), "\n"); print(table(out$rf[nzchar(out$occasion)]))
cat("set but no occasion factor:", sum(out$rf != "(absent)" & !nzchar(out$occasion)), "\n")
