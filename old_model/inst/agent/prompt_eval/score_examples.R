#!/usr/bin/env Rscript
# Score the local agent against the prompt -> StudySpec examples.
#   Rscript inst/agent/prompt_eval/score_examples.R            # all rows
#   Rscript inst/agent/prompt_eval/score_examples.R test       # frozen test split only
# Writes inst/agent/prompt_eval/results/<split>_by_example.tsv and prints the summary.
suppressMessages(devtools::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
split <- if (length(args)) args[1] else NULL
res <- .simitall_evaluate_examples(split = split)
dir.create("inst/agent/prompt_eval/results", showWarnings = FALSE, recursive = TRUE)
f <- file.path("inst/agent/prompt_eval/results", paste0(split %||% "all", "_by_example.tsv"))
utils::write.table(res$by_example, f, sep = "\t", row.names = FALSE, quote = FALSE)
str(res$summary)
cat("\nPer workflow (reading recall %, primary analysis correct %):\n")
b <- res$by_example
agg <- do.call(rbind, lapply(split(b, b$workflow), function(d) data.frame(
  workflow = d$workflow[1], n = nrow(d),
  reading = if (sum(d$reading_n)) round(100 * sum(d$reading_hit) / sum(d$reading_n)) else NA,
  primary = round(100 * mean(d$primary_ok, na.rm = TRUE)),
  invented = sum(d$invented > 0), defaulted = sum(d$defaulted_not_asked > 0))))
print(agg[order(agg$reading), ], row.names = FALSE)
cat("\nWrote", f, "\n")
