#!/usr/bin/env Rscript
# Score the lexicon reader alone and lexicon + trained model on the
# prompt -> StudySpec examples, with the same scorer as score_examples.R.
#   python inst/python/prompt_model.py --model inst/models/simitall-tagger \
#       --examples inst/agent/prompt_eval/study_spec_examples_draft.jsonl --out preds.json
#   Rscript inst/agent/prompt_eval/score_model.R preds.json [dev|test]
suppressMessages(devtools::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: score_model.R preds.json [dev|test]")
preds <- jsonlite::read_json(args[1])
split <- if (length(args) > 1) args[2] else NULL
readers <- list(
  lexicon = .simitall_lexicon_read_prompt,
  hybrid = function(q) .simitall_model_merge(.simitall_lexicon_read_prompt(q), preds[[q]], q)
)
res <- lapply(readers, function(f) .simitall_evaluate_examples(split = split, view = f))
tab <- do.call(cbind, lapply(res, function(r) unlist(r$summary)))
cat(sprintf("Split: %s   predictions: %s\n\n", split %||% "all", args[1]))
print(tab)
dir.create("inst/agent/prompt_eval/results", showWarnings = FALSE)
for (n in names(res)) utils::write.table(res[[n]]$by_example, file.path("inst/agent/prompt_eval/results",
  paste0("model_", n, "_", split %||% "all", ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
b <- merge(res$lexicon$by_example[, c("id", "reading_hit", "primary_ok")], res$hybrid$by_example[, c("id", "reading_hit", "primary_ok", "invented")],
           by = "id", suffixes = c("_lex", "_hyb"))
cat("\nexamples the model improved:", sum(b$reading_hit_hyb > b$reading_hit_lex | (b$primary_ok_hyb %in% TRUE & !b$primary_ok_lex %in% TRUE)),
    "  made worse:", sum(b$reading_hit_hyb < b$reading_hit_lex | (b$primary_ok_lex %in% TRUE & !b$primary_ok_hyb %in% TRUE)), "\n")
