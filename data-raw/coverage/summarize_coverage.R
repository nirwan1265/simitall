# Verdict per question from run_gold_specs.R output (see docs/SIMULATION_COVERAGE.md).
# usage: QUESTION_SET=second Rscript summarize_coverage.R out.tsv run_a.tsv [run_b.tsv ...]
SIMITALL <- Sys.getenv("SIMITALL_DIR", "~/Documents/Github/simitall")
args <- commandArgs(TRUE); out <- args[1]
suppressMessages(devtools::load_all(SIMITALL, quiet = TRUE))
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character", na.strings = character())
x <- do.call(rbind, lapply(args[-1], rd))
tab <- simitall:::.suite_tsv("simulators")
tab <- tab[tab$field == "assays", , drop = FALSE]
status_of <- function(a) {
  hit <- tab[vapply(strsplit(tab$values, "|", fixed = TRUE), function(v) a %in% v, TRUE), , drop = FALSE]
  if (!nrow(hit)) "unrouted" else if (hit$status[1] %in% c("ready", "treatment")) "ready" else hit$status[1]
}
# what routing sees: coarse assays plus Biodes assay_kinds (a kind can replace its coarse label)
needs <- lapply(seq_len(nrow(x)), function(i) {
  a <- strsplit(x$b_assays[i], " | ", fixed = TRUE)[[1]]
  k <- if (!is.null(x$assay_kinds)) strsplit(x$assay_kinds[i], ",", fixed = TRUE)[[1]] else character()
  simitall:::.suite_assay_values(list(assays = as.list(a[nzchar(a)]), assay_kinds = as.list(k[nzchar(k)])))
})
# what the run itself reported as not simulated (e.g. genetics with no genotype factor) overrides the route table
skipped <- lapply(strsplit(x$not_simulated, " | ", fixed = TRUE), function(z) {
  z <- z[nzchar(z)]; st <- sub(".*:", "", z)
  stats::setNames(rep(st, lengths(strsplit(sub(":[^:]*$", "", z), ","))), unlist(strsplit(sub(":[^:]*$", "", z), ",")))
})
st_list <- Map(function(a, sk) { s <- vapply(a, status_of, ""); hit <- intersect(a, names(sk)); s[hit] <- sk[hit]; s }, needs, skipped)
x$needs <- vapply(needs, paste, "", collapse = ",")
x$needs_status <- vapply(seq_along(needs), function(i) if (length(needs[[i]])) paste0(needs[[i]], "=", st_list[[i]], collapse = ",") else "", "")
n_ready <- vapply(st_list, function(s) sum(s == "ready"), 0)
n_all <- lengths(needs)
x$verdict <- ifelse(x$stage == "simitall", "not_run_simitall_error",
             ifelse(x$stage != "ok", "not_run_design",
             ifelse(n_all == 0, "simulated_outcome_question",
             ifelse(n_ready == n_all, "simulated_all",
             ifelse(n_ready > 0, "simulated_part", "no_simulator_for_data_type")))))
x$blocker <- ifelse(x$verdict %in% c("not_run_design", "not_run_simitall_error"), x$error,
             ifelse(x$verdict %in% c("simulated_part", "no_simulator_for_data_type"),
                    vapply(st_list, function(s) paste0("not simulated: ", paste(paste0(names(s)[s != "ready"], " (", s[s != "ready"], ")"), collapse = ", ")), ""), ""))
cols <- c("id", "species", "technique", "prompt", "needs", "needs_status", "reps_assumed", "solve_for", "defaults_filled", "genetic_unit", "outcome_type", "assay_kinds", "verdict", "b_design", "b_factors", "blocker", "simulated", "not_simulated")
utils::write.table(x[, cols], out, sep = "\t", row.names = FALSE, quote = FALSE)
print(table(x$verdict))
cat("with layout defaults:", sum(nzchar(x$defaults_filled) & x$verdict %in% c("simulated_all", "simulated_outcome_question")), "of the fully simulated\n")
cat("fully simulated:", sum(x$verdict %in% c("simulated_all", "simulated_outcome_question")), "/", nrow(x), "\n")
