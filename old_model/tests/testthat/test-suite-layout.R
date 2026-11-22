# Layouts built from Biodes' study details (inst/suite/layouts.tsv), for every
# design in Biodes' catalog (fixtures/biodes_terms/catalog_models.tsv).

layout_spec <- function(levels1 = 5, levels2 = 3) {
  list(factors = list(list(name = "genotype", n_levels = levels1),
                      list(name = "water", n_levels = levels2, applied_to = "main_plot"),
                      list(name = "time", n_levels = 3),
                      list(name = "check", n_levels = 2)),
       design = list(reps = 3, locations = 2))
}

write_layout_package <- function(book, spec, model, dir) {
  dir.create(dir, recursive = TRUE)
  utils::write.csv(book, file.path(dir, "field_layout.csv"), row.names = FALSE)
  meta <- cbind(book, sample_id = sprintf("S%04d", seq_len(nrow(book))), primary_outcome = "")
  utils::write.table(meta, file.path(dir, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  spec$handoff_version <- "0.4.0"
  spec$design$chosen <- attr(book, "design")
  spec$design$randomization_seed <- attr(book, "seed")
  spec$power <- list(minimum_detectable_difference = 4, residual_sd = 1)
  spec$analysis <- list(model = model, terms = list(layout_columns = intersect(
    c("treatment", "treatment2", "block", "rep", "row", "col", "unit", "subject", "time", "env", "period", "sequence", "check"), names(book))))
  jsonlite::write_json(spec, file.path(dir, "study_spec.json"), auto_unbox = TRUE)
  dir
}

test_that("every catalog design with a layout rule gives a layout holding its model's terms, which simulates and fits", {
  skip_if_not_installed("lme4")
  cat_models <- read.delim(test_path("fixtures", "biodes_terms", "catalog_models.tsv"), stringsAsFactors = FALSE)
  root <- tempfile("layouts_")
  no_rule <- failed <- character()
  for (i in seq_len(nrow(cat_models))) {
    id <- cat_models$id[i]
    spec <- if (id == "graeco_latin_square") layout_spec(5, 5) else if (id %in% c("stepped_wedge_cluster", "randomized_matched_pairs")) layout_spec(2, 3) else layout_spec()
    book <- simitall_layout_from_spec(spec, id, seed = 1)
    if (inherits(book, "simitall_layout_needs")) {
      no_rule <- c(no_rule, id)
      next
    }
    vars <- setdiff(all.vars(stats::as.formula(cat_models$model[i])), "y")
    expect_true(all(vars %in% names(book)), info = paste(id, "lacks", paste(setdiff(vars, names(book)), collapse = ", ")))
    d <- write_layout_package(book, spec, cat_models$model[i], file.path(root, id))
    res <- tryCatch({
      suppressMessages(suppressWarnings(simitall_from_biodes(d)))
      suppressMessages(suppressWarnings(bioana_fit_biodes(d)))
    }, error = function(e) e)
    if (inherits(res, "error")) failed <- c(failed, paste0(id, ": ", conditionMessage(res)))
  }
  expect_equal(failed, character())
  expect_setequal(no_rule, c("mapping_population", "split_split_plot"))   # declared in layouts.tsv, with the reason
})

test_that("layouts follow the study's factors and ask instead of inventing level counts", {
  s <- list(factors = list(list(name = "genotype", n_levels = 6), list(name = "water")), design = list(reps = 4))
  need <- simitall_layout_from_spec(s, "factorial")
  expect_s3_class(need, "simitall_layout_needs")
  expect_match(need$needs, "number of levels of 'water'")
  expect_s3_class(simitall_layout_from_spec(list(factors = list()), "rcbd"), "simitall_layout_needs")

  b <- simitall_layout_from_spec(s, "rcbd")
  expect_equal(nrow(b), 24)
  expect_setequal(unique(b$treatment), paste0("genotype_", 1:6))
  expect_equal(sort(unique(b$block)), as.character(1:4))
  expect_true(all(table(b$block, b$treatment) == 1))

  sp <- list(factors = list(list(name = "rootstock", n_levels = 4), list(name = "irrigation", n_levels = 2, applied_to = "main_plot")),
             design = list(reps = 3))
  b <- simitall_layout_from_spec(sp, "split")
  expect_setequal(unique(b$treatment), paste0("irrigation_", 1:2))      # main plots
  expect_setequal(unique(b$treatment2), paste0("rootstock_", 1:4))
  expect_true(all(tapply(b$treatment, b$main_plot_id, function(x) length(unique(x))) == 1))

  aug <- simitall_layout_from_spec(list(factors = list(list(name = "genotype", n_levels = 20), list(name = "check", n_levels = 3)),
                                        design = list(reps = 4)), "augmented")
  expect_equal(nrow(aug), 20 + 3 * 4)
  expect_true(all(table(aug$treatment[aug$check == "check"]) == 4))

  met <- simitall_layout_from_spec(list(factors = list(list(name = "genotype", n_levels = 10)), design = list(reps = 2, locations = 3)), "met_fa")
  expect_equal(nrow(met), 10 * 2 * 3)
  expect_equal(length(unique(met$env)), 3)

  noreps <- simitall_layout_from_spec(list(factors = list(list(name = "genotype", n_levels = 8))), "alpha")
  expect_true(any(grepl("reps = 2", attr(noreps, "assumptions"))))
  expect_true(any(grepl("block size = 2", attr(noreps, "assumptions"))))
})
