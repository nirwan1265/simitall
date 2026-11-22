# Build a Biodes-style design folder (study_spec.json, field_layout.csv,
# sample_metadata_template.tsv) the way Biodes write_design_package() does,
# with the layout from simitall_layout_from_spec(). Used by the simulator
# route tests so each route can be tried on a real layout.
make_design_folder <- function(design, model, factors, reps = 3, extra = list(), seed = 42) {
  spec <- c(list(system = "field", organism = "Zea mays", factors = factors,
                 design = list(reps = reps, chosen = design, randomization_seed = seed),
                 primary_outcome = list(name = "yield", type = "continuous"),
                 power = list(residual_sd = 10, minimum_detectable_difference = 10)), extra)
  book <- simitall_layout_from_spec(spec, design, seed = seed)
  if (inherits(book, "simitall_layout_needs")) stop(paste(book$needs, collapse = "; "))
  terms <- list(layout_columns = intersect(c("treatment", "treatment2", "block", "rep", "row", "col", "unit", "subject", "time", "env"), names(book)))
  terms$treatment <- factors[[1]]$name
  if (length(factors) > 1) terms$treatment2 <- factors[[2]]$name
  spec$analysis <- list(model = model, terms = terms)
  spec$handoff_version <- "0.4.0"
  dir <- file.path(tempfile("design_"), design)
  dir.create(dir, recursive = TRUE)
  jsonlite::write_json(spec, file.path(dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(as.data.frame(book), file.path(dir, "field_layout.csv"), row.names = FALSE)
  meta <- as.data.frame(book)
  meta$sample_id <- sprintf("SF_%s_%04d", design, seq_len(nrow(meta)))
  meta$date_measured <- ""; meta$primary_outcome <- ""; meta$notes <- ""
  utils::write.table(meta, file.path(dir, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  dir
}

design_factor <- function(name, n, levels = NULL) {
  c(list(name = name, n_levels = n), if (!is.null(levels)) list(levels = as.list(levels)))
}
