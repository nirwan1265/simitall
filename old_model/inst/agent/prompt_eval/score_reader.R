#!/usr/bin/env Rscript
# Score three readers on the prompt -> StudySpec examples with the same scorer:
#   agent    current local agent (simitall_ask, provider = "local")
#   lexicon  the four lexicon tables only (.simitall_lexicon_read_prompt)
#   merged   lexicon first, agent values for fields the lexicon left empty
#   Rscript inst/agent/prompt_eval/score_reader.R [dev|test]
# Agent views are cached in results/agent_views.rds (the agent is slow).
suppressMessages(devtools::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
split <- if (length(args)) args[1] else NULL
dir.create("inst/agent/prompt_eval/results", showWarnings = FALSE)
cache_file <- "inst/agent/prompt_eval/results/agent_views.rds"
cache <- if (file.exists(cache_file)) readRDS(cache_file) else list()
agent_view <- function(q) {
  if (is.null(cache[[q]])) {
    a <- tryCatch(suppressMessages(suppressWarnings(simitall_ask(q, provider = "local"))), error = function(e) list())
    cache[[q]] <<- .simitall_eval_view(a)
  }
  cache[[q]]
}
merge_views <- function(lex, ag) {
  out <- lex
  for (k in c("goal", "organism", "primary_analysis", "region")) if (is.null(out[[k]]) && !is.null(ag[[k]])) out[[k]] <- ag[[k]]
  for (k in c("design", "parameters")) for (f in names(ag[[k]])) if (is.null(out[[k]][[f]])) out[[k]][[f]] <- ag[[k]][[f]]
  if (!length(out$secondary_analyses) && length(ag$secondary_analyses) && identical(out$primary_analysis, ag$primary_analysis)) {
    out$secondary_analyses <- ag$secondary_analyses
  }
  out$defaulted <- ag$defaulted
  out
}
readers <- list(
  agent = agent_view,
  lexicon = .simitall_lexicon_read_prompt,
  merged = function(q) merge_views(.simitall_lexicon_read_prompt(q), agent_view(q))
)
res <- lapply(readers, function(f) .simitall_evaluate_examples(split = split, view = f))
saveRDS(cache, cache_file)
tab <- do.call(cbind, lapply(res, function(r) unlist(r$summary)))
cat(sprintf("Split: %s\n\n", split %||% "all"))
print(tab)
for (n in names(res)) utils::write.table(res[[n]]$by_example, file.path("inst/agent/prompt_eval/results",
  paste0("reader_", n, "_", split %||% "all", ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
