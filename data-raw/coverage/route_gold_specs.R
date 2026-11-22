# Simulation coverage of Biodes' paper questions (see docs/SIMULATION_COVERAGE.md).
# Set the two repo paths; questions = test.jsonl + the first 236 of dev.jsonl (1000).
BIODES <- Sys.getenv("BIODES_DIR", "~/Documents/Github/BioDesignGPT")
SIMITALL <- Sys.getenv("SIMITALL_DIR", "~/Documents/Github/simitall")
.split <- function(f) lapply(readLines(file.path(BIODES, "training/paper_questions/splits", f)), jsonlite::fromJSON, simplifyVector = FALSE)
QUESTIONS <- function() utils::head(c(.split("test.jsonl"), .split("dev.jsonl")), 1000)
suppressMessages(devtools::load_all(SIMITALL, quiet = TRUE))
qs <- QUESTIONS()
pat <- simitall:::.suite_rules("genetics")$genotype_pattern
res <- lapply(qs, function(q) {
  e <- q$expected_spec
  rt <- simitall:::.suite_match_routes(e)
  assays <- simitall:::.suite_norm(unlist(e$assays))
  fac <- vapply(e$factors %||% list(), function(f) f$name %||% "", "")
  has_geno <- any(grepl(pat, simitall:::.suite_norm(fac)))
  gen <- rt[rt$layer == "genetics", , drop = FALSE]
  ar <- rt[rt$layer == "assay", , drop = FALSE]
  # every requested assay -> its route status (assays that only switch on genetics count as genetics)
  tab <- simitall:::.suite_tsv("simulators")
  st <- vapply(assays, function(a) {
    hit <- tab[tab$field == "assays" & vapply(strsplit(tab$values, "|", fixed = TRUE), function(v) a %in% v, TRUE), , drop = FALSE]
    if (!nrow(hit)) "unrouted" else paste0(hit$layer[1], ":", hit$status[1])
  }, "")
  data.frame(id = q$id, species = q$species %||% "", technique = q$technique %||% "",
             assays = paste(assays, collapse = ","), assay_status = paste(paste0(assays, "=", st), collapse = ","),
             n_assays = length(assays), n_ready = sum(grepl("ready", st)),
             genetics = nrow(gen) > 0, genotype_factor = has_geno,
             mating = paste(unlist(e$population$mating), collapse = ","), study_type = e$study_type %||% "",
             outcome_type = e$primary_outcome$type %||% "", system = e$system %||% "",
             factors = paste(fac, collapse = ","), gaps = paste(unlist(q$gaps), collapse = ","),
             stringsAsFactors = FALSE)
})
g <- do.call(rbind, res)
g$verdict <- with(g, ifelse(n_assays == 0, "outcome_only",
                    ifelse(n_ready == n_assays, ifelse(genetics & !genotype_factor, "ready_but_needs_genotype_factor", "all_ready"),
                    ifelse(n_ready > 0, "partial", "none_ready"))))
utils::write.table(g, "route_gold_specs.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
print(table(g$verdict))
a <- unlist(strsplit(g$assay_status, ",")); a <- a[nzchar(a)]
print(sort(table(a), decreasing = TRUE))
