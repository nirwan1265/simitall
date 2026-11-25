`%||%` <- function(x, y) if (is.null(x)) y else x
sp <- "/home/claude/BioDesignGPT/training/paper_questions/splits"
qs <- unlist(lapply(c("test","dev","dev2","test_draw4","test_heldout","test_newkinds","train")   # test_draw5 is locked: never read, function(f)
  lapply(readLines(file.path(sp, paste0(f, ".jsonl"))), jsonlite::fromJSON, simplifyVector = FALSE)), recursive = FALSE)
qs <- qs[!duplicated(vapply(qs, `[[`, "", "id"))]
gold <- setNames(lapply(qs, function(q) { f <- q$expected_spec$factors; if (!length(f)) return(setNames(integer(), character()))
  setNames(vapply(f, function(z) if (is.null(z$n_levels)) NA_integer_ else as.integer(z$n_levels), 1L), vapply(f, function(z) z$name %||% "?", "")) }), vapply(qs, `[[`, "", "id"))
`%||%` <- function(x, y) if (is.null(x)) y else x
parse <- function(s) { if (!nzchar(s)) return(setNames(integer(), character())); p <- strsplit(s, " | ", fixed = TRUE)[[1]]
  setNames(suppressWarnings(as.integer(sub(".*\\((.*)\\)$", "\\1", p))), sub("\\(.*$", "", p)) }
score <- function(file) {
  x <- read.delim(file, colClasses = "character", quote = ""); r <- lapply(seq_len(nrow(x)), function(i) {
    g <- gold[[x$id[i]]]; g <- g[!is.na(g)]; p <- parse(x$factors[i])
    if (!length(g)) return(c(gold = 0, right = 0, wrong = 0, open = 0, missing = 0, all_right = NA))
    pv <- p[names(g)]
    c(gold = length(g), right = sum(pv == g, na.rm = TRUE), wrong = sum(!is.na(pv) & pv != g), open = sum(names(g) %in% names(p) & is.na(pv)),
      missing = sum(!names(g) %in% names(p)), all_right = all(!is.na(pv) & pv == g))
  })
  m <- do.call(rbind, r); s <- colSums(m[, 1:5]); c(s, specs_all_counts_right = sum(m[, 6], na.rm = TRUE), specs_with_gold = sum(!is.na(m[, 6])))
}
o <- score("old.tsv"); n <- score("new.tsv"); print(rbind(old = o, new = n))
