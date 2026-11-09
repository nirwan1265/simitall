#!/usr/bin/env Rscript
# Compare the organism/technique lexicon with the current local agent on the
# prompt -> StudySpec examples.
#   Rscript inst/agent/prompt_eval/score_lexicon.R [dev|test]
# Needs results/all_by_example.tsv from score_examples.R for the agent's picks.
suppressMessages(devtools::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
split <- if (length(args)) args[1] else NULL
ex <- .simitall_eval_examples()
if (!is.null(split)) ex <- Filter(function(e) identical(e$split, split), ex)
agent <- utils::read.delim("inst/agent/prompt_eval/results/all_by_example.tsv", colClasses = "character")

rows <- lapply(ex, function(e) {
  s <- e$expected_spec
  o <- .simitall_detect_organism(e$prompt)
  t <- .simitall_detect_techniques(e$prompt)
  a <- agent[agent$id == e$id, ]
  old_org <- .simitall_agent_species(e$prompt)
  data.frame(
    id = e$id, split = e$split, batch = e$batch, workflow = e$workflow,
    want_org = s$organism %||% NA, lex_org = o$organism, old_org = if (identical(old_org, "unspecified")) NA else old_org,
    org_should_ask = "organism" %in% unlist(e$must_ask),
    want_primary = s$primary_analysis %||% NA, lex_primary = t$primary_analysis,
    old_primary = if (nrow(a)) a$got_primary[1] else NA,
    want_unsupported = identical(s$execution_status, "unsupported"),
    lex_unsupported = identical(t$execution_status, "unsupported"),
    lex_flags = paste(.simitall_lexicon_checks(e$prompt, o$organism), collapse = "; "),
    stringsAsFactors = FALSE)
})
r <- do.call(rbind, rows)
pct <- function(x) round(100 * mean(x), 1)
has_org <- !is.na(r$want_org)
no_org <- is.na(r$want_org)
has_pri <- !is.na(r$want_primary)
cat(sprintf("Examples: %d (%s)\n\n", nrow(r), split %||% "all"))
cat("ORGANISM                                     current agent   lexicon\n")
cat(sprintf("  named organism found correctly (n=%3d)   %6.1f%%       %6.1f%%\n", sum(has_org),
            pct(r$old_org[has_org] %in% r$want_org[has_org] & !is.na(r$old_org[has_org]) & r$old_org[has_org] == r$want_org[has_org]),
            pct(!is.na(r$lex_org[has_org]) & r$lex_org[has_org] == r$want_org[has_org])))
cat(sprintf("  no organism invented when none named (n=%3d) %6.1f%%     %6.1f%%\n", sum(no_org),
            pct(is.na(r$old_org[no_org])), pct(is.na(r$lex_org[no_org]))))
cat("\nANALYSIS (primary)\n")
cat(sprintf("  correct primary analysis (n=%3d)         %6.1f%%       %6.1f%%\n", sum(has_pri),
            pct(!is.na(r$old_primary[has_pri]) & r$old_primary[has_pri] == r$want_primary[has_pri]),
            pct(!is.na(r$lex_primary[has_pri]) & r$lex_primary[has_pri] == r$want_primary[has_pri])))
cat(sprintf("  wrong analysis chosen (should be none or other) %4d          %4d\n",
            sum(!is.na(r$old_primary[has_pri]) & r$old_primary[has_pri] != r$want_primary[has_pri]),
            sum(!is.na(r$lex_primary[has_pri]) & r$lex_primary[has_pri] != r$want_primary[has_pri])))
cat(sprintf("  unsupported technique recognised (n=%3d)          lexicon %5.1f%%\n", sum(r$want_unsupported),
            pct(r$lex_unsupported[r$want_unsupported])))
cat(sprintf("  supported request wrongly called unsupported            %4d\n", sum(r$lex_unsupported & !r$want_unsupported)))
cat("\nBiology flags raised:\n")
f <- r[nzchar(r$lex_flags), c("id", "lex_flags")]
for (i in seq_len(nrow(f))) cat("  ", f$id[i], ":", f$lex_flags[i], "\n")
dir.create("inst/agent/prompt_eval/results", showWarnings = FALSE)
out <- file.path("inst/agent/prompt_eval/results", paste0("lexicon_", split %||% "all", ".tsv"))
utils::write.table(r, out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nWrote", out, "\n")
