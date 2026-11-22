# Simulation coverage of Biodes' paper questions (see docs/SIMULATION_COVERAGE.md).
# Set the two repo paths; questions = test.jsonl + the first 236 of dev.jsonl (1000).
BIODES <- Sys.getenv("BIODES_DIR", "~/Documents/Github/BioDesignGPT")
SIMITALL <- Sys.getenv("SIMITALL_DIR", "~/Documents/Github/simitall")
.split <- function(f) lapply(readLines(file.path(BIODES, "training/paper_questions/splits", f)), jsonlite::fromJSON, simplifyVector = FALSE)
QUESTIONS <- function() utils::head(c(.split("test.jsonl"), .split("dev.jsonl")), 1000)
# usage: Rscript <this file> <start> <end> <out.tsv>
args <- commandArgs(TRUE); i0 <- as.integer(args[1]); i1 <- as.integer(args[2]); out <- args[3]
suppressMessages({devtools::load_all(BIODES, quiet = TRUE); devtools::load_all(SIMITALL, quiet = TRUE)})
qs <- QUESTIONS()
small <- list(simulate_rnaseq_experiment = list(n_genes = 200), simulate_rnaseq_from_gwas = list(n_genes = 200, n_cis_eqtl = 20, n_trans_eqtl = 5),
              simulate_scrnaseq_experiment = list(n_genes = 150, cells_per_replicate = 20),
              simulate_chipseq = list(n_peaks = 30, n_reads = 3000), simulate_atacseq = list(n_peaks = 30, n_reads = 3000),
              suite_peak_genome = list(chromosome_length = 60000L, n_genes = 40L),
              suite_genotype_panel = list(n_markers = 600), suite_genotype_cross = list(chromosome_length_bp = 3000L, n_chromosomes = 2L, snp_rate = 0.05))
clean <- function(x) gsub("[\t\r\n]+", " ", substr(paste(x, collapse = " | "), 1, 400))
rows <- list()
for (i in i0:i1) {
  q <- qs[[i]]
  d <- file.path(tempdir(), q$id)
  r <- list(id = q$id, species = q$species %||% "", technique = q$technique %||% "", prompt = clean(q$prompt))
  t0 <- Sys.time()
  run <- tryCatch({
    utils::capture.output(x <- suppressMessages(simitall_from_question(q$prompt, outdir = d, simulate = FALSE, use_suggestions = FALSE,
                                                                       confirm = FALSE, figures = FALSE, analyze = FALSE)))
    x
  }, error = function(e) e)
  if (inherits(run, "error")) {
    r$stage <- "biodes_or_layout"; r$error <- clean(conditionMessage(run))
  } else {
    s <- run$answer$spec
    r$b_assays <- clean(unlist(s$assays)); r$b_outcome_type <- clean(s$primary_outcome$type %||% "")
    r$b_mating <- clean(unlist(s$population$mating)); r$b_study_type <- clean(s$study_type %||% "")
    r$b_design <- clean(run$design %||% s$design$chosen %||% ""); r$b_factors <- clean(vapply(s$factors %||% list(), function(f) paste0(f$name %||% "?", "(", f$n_levels %||% "?", ")"), ""))
    r$layout_issues <- clean(run$layout_issues); r$open_n <- length(run$open_questions)
    if (!file.exists(file.path(d, "study_spec.json")) || !file.exists(file.path(d, "field_layout.csv"))) {
      r$stage <- "no_design_folder"; r$error <- clean(c(run$open_questions))
    } else {
      rt <- tryCatch(simitall:::.suite_match_routes(jsonlite::read_json(file.path(d, "study_spec.json"))), error = function(e) NULL)
      if (!is.null(rt)) { r$routes <- clean(paste0(rt$route, ":", rt$status)) }
      sim <- tryCatch({ utils::capture.output(p <- suppressMessages(simitall_from_biodes(d, params = small))); p }, error = function(e) e)
      if (inherits(sim, "error")) { r$stage <- "simitall"; r$error <- clean(conditionMessage(sim)) }
      else {
        r$stage <- "ok"
        tr <- jsonlite::read_json(file.path(d, "sim_truth.json"))
        r$simulated <- clean(c(tr$outcome_type, if (!is.null(tr$genetics)) paste0("genetics:", tr$genetics$genotype_source$source), names(tr$assays)))
        r$not_simulated <- clean(vapply(tr$not_simulated %||% list(), function(z) paste0(z$requested, ":", z$status), ""))
        r$outcome_simulated <- isTRUE(tr$outcome_simulated)
      }
    }
  }
  r$seconds <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  unlink(d, recursive = TRUE)
  rows[[length(rows) + 1]] <- r
}
cols <- c("id","species","technique","prompt","stage","error","b_assays","b_outcome_type","b_mating","b_study_type","b_design","b_factors","layout_issues","open_n","routes","simulated","not_simulated","outcome_simulated","seconds")
df <- do.call(rbind, lapply(rows, function(r) as.data.frame(setNames(lapply(cols, function(k) if (is.null(r[[k]])) "" else as.character(r[[k]])), cols))))
utils::write.table(df, out, sep = "\t", row.names = FALSE, quote = FALSE)
