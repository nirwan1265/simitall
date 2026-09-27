#!/usr/bin/env Rscript

# Minimal GenomeAdmixR smoke test for simitall.
# This simulates a small one-chromosome, three-founder ancestry scenario. It is
# a simulation of known ancestry truth, not local-ancestry inference on real data.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1L]] else "results/ancestry_tract_smoke"

# Allow the repository's ignored development library without affecting a normal
# installed-package workflow. Resolve relative to this script when run by
# Rscript so testthat can execute it from its own working directory.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else NA_character_
repo_root <- if (!is.na(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", "..", ".."), mustWork = FALSE)
} else {
  getwd()
}
local_lib <- file.path(repo_root, ".Rlib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

if (!requireNamespace("GenomeAdmixR", quietly = TRUE)) {
  stop(
    "GenomeAdmixR is required for this example. Install it with ",
    "install.packages('GenomeAdmixR') or use the simitall ancestry profile once available.",
    call. = FALSE
  )
}

set.seed(20260926)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

markers <- seq(0, 1, length.out = 21L)
simulation <- GenomeAdmixR::simulate_admixture(
  module = GenomeAdmixR::ancestry_module(
    number_of_founders = 3L,
    morgan = 1,
    markers = markers,
    track_junctions = TRUE
  ),
  pop_size = 20L,
  total_runtime = 5L,
  num_threads = 1L,
  verbose = FALSE
)

initial <- as.data.frame(simulation$initial_frequency)
final <- as.data.frame(simulation$final_frequency)
if (!all(c("location", "ancestor", "frequency") %in% names(final))) {
  stop("GenomeAdmixR returned an unexpected final_frequency schema.", call. = FALSE)
}

by_ancestor <- split(final$frequency, final$ancestor)
summary <- data.frame(
  ancestor = as.integer(names(by_ancestor)),
  mean_frequency = vapply(by_ancestor, mean, numeric(1)),
  min_frequency = vapply(by_ancestor, min, numeric(1)),
  max_frequency = vapply(by_ancestor, max, numeric(1)),
  stringsAsFactors = FALSE
)

utils::write.csv(initial, file.path(out_dir, "initial_ancestry_frequency.csv"), row.names = FALSE)
utils::write.csv(final, file.path(out_dir, "final_ancestry_frequency.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "ancestry_frequency_summary.csv"), row.names = FALSE)

png(file.path(out_dir, "final_ancestry_frequency.png"), width = 1200, height = 800, res = 150)
ancestors <- sort(unique(final$ancestor))
colors <- grDevices::hcl.colors(length(ancestors), "Dark 3")
plot(
  NA,
  xlim = range(final$location), ylim = c(0, 1),
  xlab = "Genetic position (Morgans)", ylab = "Founder ancestry frequency",
  main = "GenomeAdmixR smoke test: ancestry after five generations"
)
for (i in seq_along(ancestors)) {
  rows <- final$ancestor == ancestors[[i]]
  lines(final$location[rows], final$frequency[rows], type = "l", lwd = 2, col = colors[[i]])
}
legend("topright", legend = paste("Founder", ancestors), col = colors, lwd = 2, bty = "n")
dev.off()

metadata <- list(
  package = "GenomeAdmixR",
  version = as.character(utils::packageVersion("GenomeAdmixR")),
  seed = 20260926L,
  model = "one chromosome, three founders, random-mating ancestry simulation",
  markers = length(markers),
  population_size = 20L,
  generations = 5L,
  outputs = c(
    "initial_ancestry_frequency.csv",
    "final_ancestry_frequency.csv",
    "ancestry_frequency_summary.csv",
    "final_ancestry_frequency.png"
  ),
  limitation = "This smoke test simulates ancestry truth only. It does not phase real data, infer local ancestry, or run an ancestry-aware GWAS."
)
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(metadata, file.path(out_dir, "metadata.json"), pretty = TRUE, auto_unbox = TRUE)
}

cat("GenomeAdmixR ancestry smoke test completed. Outputs: ", normalizePath(out_dir), "\n", sep = "")
