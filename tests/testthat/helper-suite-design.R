# Build a design folder the way Biodes does: studyforge::layout_from_spec()
# then studyforge::write_design_package(). Tests that use it skip when Biodes
# (R package 'studyforge') is not installed.
make_design_folder <- function(design, model = NULL, factors, reps = 3, extra = list(), seed = 42) {
  testthat::skip_if_not_installed("studyforge")
  spec <- utils::modifyList(list(system = "field", organism = "Zea mays", factors = factors,
                 design = list(reps = reps, chosen = design, randomization_seed = seed),
                 primary_outcome = list(name = "yield", type = "continuous"),
                 power = list(residual_sd = 10, minimum_detectable_difference = 10)), extra)
  book <- studyforge::layout_from_spec(spec, seed = seed)
  dir <- file.path(tempfile("design_"), design)
  suppressMessages(studyforge::write_design_package(book, spec, dir))
  dir
}

design_factor <- function(name, n, levels = NULL, applied_to = NULL) {
  c(list(name = name, n_levels = n), if (!is.null(levels)) list(levels = as.list(levels)),
    if (!is.null(applied_to)) list(applied_to = applied_to))
}

# A design folder from a full spec (Biodes layout + package writer).
biodes_folder <- function(spec, name = "design", seed = 42) {
  testthat::skip_if_not_installed("studyforge")
  book <- studyforge::layout_from_spec(spec, seed = seed)
  dir <- file.path(tempfile("design_"), name)
  suppressMessages(studyforge::write_design_package(book, spec, dir))
  dir
}
