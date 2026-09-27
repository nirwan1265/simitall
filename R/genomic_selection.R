.simitall_marker_matrix <- function(genotype_file, training_samples = NULL) {
  genotype_data <- .simitall_read_genotypes(genotype_file)
  dosage <- t(genotype_data$genotype)
  storage.mode(dosage) <- "numeric"
  if (is.null(training_samples)) training_samples <- rownames(dosage)
  training_index <- match(training_samples, rownames(dosage))
  if (anyNA(training_index)) {
    stop(
      "Training samples not found in genotype data: ",
      paste(training_samples[is.na(training_index)], collapse = ", ")
    )
  }
  training <- dosage[training_index, , drop = FALSE]
  means <- colMeans(training, na.rm = TRUE)
  usable <- is.finite(means) & vapply(seq_len(ncol(training)), function(j) {
    value <- training[, j]
    value <- value[is.finite(value)]
    length(value) >= 2L && stats::var(value) > 0
  }, logical(1L))
  if (!any(usable)) stop("No polymorphic markers remain in the training set")
  dosage <- dosage[, usable, drop = FALSE]
  means <- means[usable]
  for (j in seq_len(ncol(dosage))) {
    missing <- !is.finite(dosage[, j])
    if (any(missing)) dosage[missing, j] <- means[j]
  }
  list(
    dosage = dosage,
    centered = sweep(dosage, 2L, means, "-"),
    means = means,
    marker_ids = colnames(dosage),
    variants = genotype_data$variants[
      match(colnames(dosage), genotype_data$variants$id),
      ,
      drop = FALSE
    ]
  )
}

.simitall_bglr_model <- function(model) {
  switch(
    model,
    bayes_ridge = "BRR",
    bayes_a = "BayesA",
    bayes_b = "BayesB",
    bayes_c = "BayesC",
    bayesian_lasso = "BL",
    stop("Unknown BGLR model: ", model)
  )
}

#' Calculate a realized additive genomic relationship matrix
#'
#' Read VCF or genotype-TSV dosages, convert diploid genotypes to the
#' `-1/0/1` coding expected by [rrBLUP::A.mat()], and calculate a sample-by-
#' sample additive relationship matrix.
#'
#' @param genotype_file VCF, VCF.GZ, or genotype TSV path.
#' @param min_maf Optional minimum minor-allele frequency.
#' @param max_missing Optional maximum marker missing-data proportion.
#' @param impute_method Imputation method accepted by [rrBLUP::A.mat()].
#' @param out_file Optional `.rds` or TSV output path.
#'
#' @return A genomic relationship matrix, invisibly when `out_file` is used.
#' @examples
#' \dontrun{
#' K <- calculate_genomic_kinship("results/gwas/demo.geno.tsv")
#' }
#' @export
calculate_genomic_kinship <- function(
    genotype_file,
    min_maf = NULL,
    max_missing = NULL,
    impute_method = c("mean", "EM"),
    out_file = NULL) {
  .simitall_require_package("rrBLUP", "genomic relationship estimation")
  impute_method <- match.arg(impute_method)
  genotype <- .simitall_read_genotypes(genotype_file)$genotype
  marker_matrix <- t(genotype) - 1
  relationship <- rrBLUP::A.mat(
    marker_matrix,
    min.MAF = min_maf,
    max.missing = max_missing,
    impute.method = impute_method
  )
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    if (grepl("\\.rds$", out_file, ignore.case = TRUE)) {
      saveRDS(relationship, out_file)
    } else {
      output <- data.frame(
        sample = rownames(relationship),
        relationship,
        check.names = FALSE
      )
      utils::write.table(
        output, out_file, sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
    return(invisible(relationship))
  }
  relationship
}

#' Fit a genomic prediction model
#'
#' Fit GBLUP, RR-BLUP, Bayesian marker regression, or random forest models to
#' phenotyped training samples while retaining predictions for every genotyped
#' candidate. GBLUP and RR-BLUP use `rrBLUP`; Bayesian models use `BGLR`; and
#' random forest uses `ranger`.
#'
#' @param genotype_file VCF, VCF.GZ, or genotype TSV path.
#' @param phenotype Data frame or phenotype TSV containing sample IDs and a
#'   numeric trait.
#' @param trait Numeric phenotype column.
#' @param sample_id Sample-ID column in `phenotype`.
#' @param model Prediction model: `"gblup"`, `"rrblup"`, `"bayes_ridge"`,
#'   `"bayes_a"`, `"bayes_b"`, `"bayes_c"`, `"bayesian_lasso"`, or
#'   `"random_forest"`.
#' @param training_samples Optional sample IDs used for model training. By
#'   default, every genotyped sample with a finite phenotype is used.
#' @param n_iter,burn_in,thin MCMC controls passed to [BGLR::BGLR()].
#' @param num_trees Number of trees passed to [ranger::ranger()].
#' @param seed Random-number seed.
#' @param verbose Show backend progress where supported.
#'
#' @return A `simitall_genomic_model` object containing fitted backend objects,
#'   marker preprocessing information, and candidate predictions.
#' @examples
#' \dontrun{
#' fit <- fit_genomic_model(
#'   genotype_file = "results/gwas/demo.geno.tsv",
#'   phenotype = "results/gwas/demo.pheno.tsv",
#'   model = "gblup"
#' )
#' predict_genomic_values(fit)
#' }
#' @export
fit_genomic_model <- function(
    genotype_file,
    phenotype,
    trait = "trait",
    sample_id = "sample",
    model = c(
      "gblup", "rrblup", "bayes_ridge", "bayes_a", "bayes_b",
      "bayes_c", "bayesian_lasso", "random_forest"
    ),
    training_samples = NULL,
    n_iter = 12000L,
    burn_in = 2000L,
    thin = 5L,
    num_trees = 500L,
    seed = 1L,
    verbose = FALSE) {
  model <- match.arg(model)
  genotype_data <- .simitall_read_genotypes(genotype_file)
  sample_ids <- colnames(genotype_data$genotype)
  phenotype <- .simitall_prepare_phenotypes(
    phenotype, sample_ids, sample_id, trait
  )
  available <- sample_ids[is.finite(phenotype[[trait]])]
  if (is.null(training_samples)) training_samples <- available
  training_samples <- unique(as.character(training_samples))
  unavailable <- setdiff(training_samples, available)
  if (length(unavailable)) {
    stop(
      "Training samples need finite phenotypes and genotypes: ",
      paste(unavailable, collapse = ", ")
    )
  }
  if (length(training_samples) < 10L) {
    stop("At least ten training samples are required")
  }
  marker_data <- .simitall_marker_matrix(genotype_file, training_samples)
  training_index <- match(training_samples, sample_ids)
  y <- phenotype[[trait]][training_index]
  set.seed(as.integer(seed))

  fit <- NULL
  marker_effects <- NULL
  predicted_trait <- rep(NA_real_, length(sample_ids))
  gebv <- rep(NA_real_, length(sample_ids))
  names(predicted_trait) <- names(gebv) <- sample_ids
  intercept <- mean(y)
  backend <- NULL

  if (model == "gblup") {
    .simitall_require_package("rrBLUP", "GBLUP genomic prediction")
    relationship <- rrBLUP::A.mat(marker_data$dosage - 1)
    training_data <- data.frame(
      sample = training_samples,
      phenotype = y,
      stringsAsFactors = FALSE
    )
    fit <- rrBLUP::kin.blup(
      data = training_data,
      geno = "sample",
      pheno = "phenotype",
      K = relationship,
      PEV = TRUE
    )
    gebv[names(fit$g)] <- as.numeric(fit$g)
    predicted_trait <- intercept + gebv
    backend <- "rrBLUP::kin.blup"
  } else if (model == "rrblup") {
    .simitall_require_package("rrBLUP", "RR-BLUP genomic prediction")
    fit <- rrBLUP::mixed.solve(
      y = y,
      Z = marker_data$centered[training_index, , drop = FALSE]
    )
    marker_effects <- as.numeric(fit$u)
    names(marker_effects) <- marker_data$marker_ids
    intercept <- as.numeric(fit$beta[1L])
    gebv[] <- as.numeric(marker_data$centered %*% marker_effects)
    predicted_trait <- intercept + gebv
    backend <- "rrBLUP::mixed.solve"
  } else if (startsWith(model, "bayes") || model == "bayesian_lasso") {
    .simitall_require_package("BGLR", "Bayesian genomic prediction")
    if (n_iter <= burn_in || burn_in < 0L || thin < 1L) {
      stop("BGLR requires n_iter > burn_in >= 0 and thin >= 1")
    }
    y_all <- rep(NA_real_, length(sample_ids))
    y_all[training_index] <- y
    save_prefix <- tempfile("simitall-bglr-")
    on.exit(
      unlink(list.files(
        dirname(save_prefix),
        pattern = paste0("^", basename(save_prefix)),
        full.names = TRUE
      )),
      add = TRUE
    )
    fit <- BGLR::BGLR(
      y = y_all,
      ETA = list(markers = list(
        X = marker_data$centered,
        model = .simitall_bglr_model(model)
      )),
      nIter = as.integer(n_iter),
      burnIn = as.integer(burn_in),
      thin = as.integer(thin),
      verbose = isTRUE(verbose),
      saveAt = save_prefix
    )
    marker_effects <- as.numeric(fit$ETA[[1L]]$b)
    names(marker_effects) <- marker_data$marker_ids
    intercept <- as.numeric(fit$mu)
    gebv[] <- as.numeric(marker_data$centered %*% marker_effects)
    predicted_trait[] <- as.numeric(fit$yHat)
    backend <- paste0("BGLR::BGLR(", .simitall_bglr_model(model), ")")
  } else if (model == "random_forest") {
    .simitall_require_package("ranger", "random-forest genomic prediction")
    fit <- ranger::ranger(
      x = as.data.frame(
        marker_data$centered[training_index, , drop = FALSE],
        check.names = FALSE
      ),
      y = y,
      num.trees = as.integer(num_trees),
      seed = as.integer(seed),
      num.threads = 1L,
      write.forest = TRUE,
      verbose = isTRUE(verbose)
    )
    predicted_trait[] <- predict(
      fit,
      data = as.data.frame(marker_data$centered, check.names = FALSE)
    )$predictions
    intercept <- mean(y)
    gebv <- predicted_trait - intercept
    backend <- "ranger::ranger"
  }

  predictions <- data.frame(
    sample = sample_ids,
    observed = phenotype[[trait]],
    gebv = as.numeric(gebv[sample_ids]),
    predicted_trait = as.numeric(predicted_trait[sample_ids]),
    training_sample = sample_ids %in% training_samples,
    stringsAsFactors = FALSE
  )
  object <- list(
    model = model,
    backend = backend,
    trait = trait,
    sample_id = sample_id,
    sample_ids = sample_ids,
    training_samples = training_samples,
    marker_ids = marker_data$marker_ids,
    marker_means = marker_data$means,
    marker_effects = marker_effects,
    intercept = intercept,
    fit = fit,
    predictions = predictions,
    seed = as.integer(seed)
  )
  class(object) <- "simitall_genomic_model"
  object
}

#' Predict genomic values from a fitted simitall model
#'
#' Return stored predictions for all candidates used during model fitting or
#' predict a new genotype file for marker-effect and random-forest models.
#' GBLUP candidates must be included when [fit_genomic_model()] is called so
#' their relationships to the training population can be estimated jointly.
#'
#' @param object A model returned by [fit_genomic_model()].
#' @param genotype_file Optional new VCF or genotype TSV. If omitted, returns
#'   predictions stored in `object`.
#'
#' @return A data frame with sample IDs, GEBVs, predicted trait values, and
#'   training-set membership where known.
#' @examples
#' \dontrun{
#' predictions <- predict_genomic_values(fit)
#' }
#' @export
predict_genomic_values <- function(object, genotype_file = NULL) {
  if (!inherits(object, "simitall_genomic_model")) {
    stop("object must be returned by fit_genomic_model()")
  }
  if (is.null(genotype_file)) return(object$predictions)
  if (object$model == "gblup") {
    stop(
      "GBLUP candidates must be present in the genotype file used for model ",
      "fitting; call predict_genomic_values() without genotype_file"
    )
  }
  genotype_data <- .simitall_read_genotypes(genotype_file)
  missing_markers <- setdiff(object$marker_ids, rownames(genotype_data$genotype))
  if (length(missing_markers)) {
    stop(
      "New genotype data are missing ", length(missing_markers),
      " model markers"
    )
  }
  dosage <- t(genotype_data$genotype[object$marker_ids, , drop = FALSE])
  for (j in seq_len(ncol(dosage))) {
    missing <- !is.finite(dosage[, j])
    if (any(missing)) dosage[missing, j] <- object$marker_means[j]
  }
  centered <- sweep(dosage, 2L, object$marker_means, "-")
  if (object$model == "random_forest") {
    predicted <- predict(
      object$fit,
      data = as.data.frame(centered, check.names = FALSE)
    )$predictions
    gebv <- predicted - object$intercept
  } else {
    gebv <- as.numeric(centered %*% object$marker_effects)
    predicted <- object$intercept + gebv
  }
  data.frame(
    sample = rownames(centered),
    observed = NA_real_,
    gebv = as.numeric(gebv),
    predicted_trait = as.numeric(predicted),
    training_sample = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Cross-validate genomic prediction
#'
#' Perform repeated random or group-wise cross-validation and calculate
#' predictive correlation, rank correlation, RMSE, MAE, and calibration slope.
#' Group-wise splitting keeps all members of a family, generation, or
#' subpopulation in the same fold and is preferred when random splitting would
#' leak close relatives between training and testing sets.
#'
#' @inheritParams fit_genomic_model
#' @param out_prefix Optional prefix for prediction, metric, and JSON outputs.
#' @param folds Number of cross-validation folds.
#' @param repeats Number of independently randomized fold assignments.
#' @param split `"random"` for individual folds or `"group"` to hold out whole
#'   groups.
#' @param group Optional phenotype column used when `split = "group"`.
#' @param true_value Optional phenotype column containing simulated true
#'   breeding values.
#' @param ... Additional arguments passed to [fit_genomic_model()].
#'
#' @return Invisibly returns fold predictions, fold metrics, summary metrics,
#'   and output paths when requested.
#' @examples
#' \dontrun{
#' cv <- cross_validate_genomic_prediction(
#'   "results/gwas/demo.geno.tsv",
#'   "results/gwas/demo.pheno.tsv",
#'   model = "gblup",
#'   folds = 5,
#'   true_value = "true_breeding_value"
#' )
#' }
#' @export
cross_validate_genomic_prediction <- function(
    genotype_file,
    phenotype,
    trait = "trait",
    sample_id = "sample",
    model = "gblup",
    out_prefix = NULL,
    folds = 5L,
    repeats = 1L,
    split = c("random", "group"),
    group = NULL,
    true_value = NULL,
    seed = 1L,
    ...) {
  split <- match.arg(split)
  genotype_samples <- colnames(.simitall_read_genotypes(genotype_file)$genotype)
  required <- c(if (split == "group") group, true_value)
  phenotype_table <- .simitall_prepare_phenotypes(
    phenotype,
    genotype_samples,
    sample_id,
    trait,
    required_columns = required
  )
  eligible <- genotype_samples[is.finite(phenotype_table[[trait]])]
  folds <- as.integer(folds)
  repeats <- as.integer(repeats)
  if (folds < 2L || folds > length(eligible)) {
    stop("folds must be between two and the number of phenotyped samples")
  }
  if (repeats < 1L) stop("repeats must be at least one")
  if (split == "group" && (is.null(group) || !nzchar(group))) {
    stop("group must name a phenotype column for group-wise splitting")
  }

  predictions <- list()
  metrics <- list()
  row_index <- 1L
  set.seed(as.integer(seed))
  for (repeat_index in seq_len(repeats)) {
    if (split == "random") {
      assignment <- sample(rep(seq_len(folds), length.out = length(eligible)))
      names(assignment) <- eligible
    } else {
      group_values <- phenotype_table[[group]][match(
        eligible, phenotype_table[[sample_id]]
      )]
      if (anyNA(group_values)) stop("Group values cannot be missing")
      groups <- sample(unique(as.character(group_values)))
      group_fold <- rep(seq_len(folds), length.out = length(groups))
      names(group_fold) <- groups
      assignment <- group_fold[as.character(group_values)]
      names(assignment) <- eligible
    }
    for (fold_index in seq_len(folds)) {
      test_samples <- names(assignment)[assignment == fold_index]
      train_samples <- setdiff(eligible, test_samples)
      fit <- fit_genomic_model(
        genotype_file = genotype_file,
        phenotype = phenotype_table,
        trait = trait,
        sample_id = sample_id,
        model = model,
        training_samples = train_samples,
        seed = seed + repeat_index * 1000L + fold_index,
        ...
      )
      predicted <- predict_genomic_values(fit)
      predicted <- predicted[
        match(test_samples, predicted$sample),
        ,
        drop = FALSE
      ]
      observed <- phenotype_table[[trait]][match(
        test_samples, phenotype_table[[sample_id]]
      )]
      true_values <- if (is.null(true_value)) {
        rep(NA_real_, length(test_samples))
      } else {
        as.numeric(phenotype_table[[true_value]][match(
          test_samples, phenotype_table[[sample_id]]
        )])
      }
      fold_predictions <- data.frame(
        repeat_id = repeat_index,
        fold = fold_index,
        sample = test_samples,
        observed = observed,
        true_value = true_values,
        gebv = predicted$gebv,
        predicted_trait = predicted$predicted_trait,
        stringsAsFactors = FALSE
      )
      predictions[[row_index]] <- fold_predictions
      correlation <- suppressWarnings(stats::cor(
        fold_predictions$observed,
        fold_predictions$predicted_trait,
        use = "complete.obs"
      ))
      rank_correlation <- suppressWarnings(stats::cor(
        fold_predictions$observed,
        fold_predictions$predicted_trait,
        method = "spearman",
        use = "complete.obs"
      ))
      truth_accuracy <- if (all(is.na(true_values))) {
        NA_real_
      } else {
        suppressWarnings(stats::cor(
          true_values, fold_predictions$gebv, use = "complete.obs"
        ))
      }
      residual <- fold_predictions$observed -
        fold_predictions$predicted_trait
      calibration <- tryCatch(
        stats::coef(stats::lm(
          observed ~ predicted_trait,
          data = fold_predictions
        ))[[2L]],
        error = function(error) NA_real_
      )
      metrics[[row_index]] <- data.frame(
        repeat_id = repeat_index,
        fold = fold_index,
        training_samples = length(train_samples),
        testing_samples = length(test_samples),
        predictive_correlation = correlation,
        rank_correlation = rank_correlation,
        truth_accuracy = truth_accuracy,
        rmse = sqrt(mean(residual^2, na.rm = TRUE)),
        mae = mean(abs(residual), na.rm = TRUE),
        calibration_slope = calibration,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  predictions <- do.call(rbind, predictions)
  metrics <- do.call(rbind, metrics)
  summary <- data.frame(
    model = model,
    split = split,
    folds = folds,
    repeats = repeats,
    predictive_correlation = mean(
      metrics$predictive_correlation, na.rm = TRUE
    ),
    rank_correlation = mean(metrics$rank_correlation, na.rm = TRUE),
    truth_accuracy = if (all(is.na(metrics$truth_accuracy))) {
      NA_real_
    } else {
      mean(metrics$truth_accuracy, na.rm = TRUE)
    },
    rmse = mean(metrics$rmse, na.rm = TRUE),
    mae = mean(metrics$mae, na.rm = TRUE),
    calibration_slope = mean(metrics$calibration_slope, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  paths <- list()
  if (!is.null(out_prefix)) {
    paths <- list(
      predictions = paste0(out_prefix, ".cv_predictions.tsv"),
      metrics = paste0(out_prefix, ".cv_metrics.tsv"),
      summary = paste0(out_prefix, ".cv_summary.json")
    )
    dir.create(dirname(paths$predictions), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      predictions, paths$predictions, sep = "\t", quote = FALSE,
      row.names = FALSE
    )
    utils::write.table(
      metrics, paths$metrics, sep = "\t", quote = FALSE, row.names = FALSE
    )
    .simitall_write_json(as.list(summary[1L, ]), paths$summary)
  }
  invisible(list(
    predictions = predictions,
    metrics = metrics,
    summary = summary,
    paths = paths
  ))
}

#' Select parents from genomic estimated breeding values
#'
#' Rank candidates by GEBV and optionally penalize genomic coancestry through a
#' greedy diversity-aware selection algorithm. Family representation can be
#' constrained to prevent one family from supplying all selected parents.
#'
#' @param predictions Prediction table, TSV, or model returned by
#'   [fit_genomic_model()]. The table must contain `sample` and `gebv`.
#' @param n_parents Number of parents to select. When `NULL`, calculated from
#'   `proportion`.
#' @param proportion Fraction selected when `n_parents` is `NULL`.
#' @param direction Select larger (`"max"`) or smaller (`"min"`) GEBVs.
#' @param kinship Optional named genomic relationship matrix.
#' @param diversity_penalty Non-negative penalty applied to mean relationship
#'   with parents already selected.
#' @param family_data Optional data frame or TSV containing sample and family
#'   columns.
#' @param sample_id,family_col Column names in `family_data`.
#' @param max_per_family Maximum selected parents from any family.
#'
#' @return A selected-parent table ordered by selection rank.
#' @examples
#' \dontrun{
#' selected <- select_parents(
#'   predict_genomic_values(fit),
#'   n_parents = 10,
#'   kinship = K,
#'   diversity_penalty = 0.25
#' )
#' }
#' @export
select_parents <- function(
    predictions,
    n_parents = NULL,
    proportion = 0.1,
    direction = c("max", "min"),
    kinship = NULL,
    diversity_penalty = 0,
    family_data = NULL,
    sample_id = "sample",
    family_col = "family",
    max_per_family = Inf) {
  direction <- match.arg(direction)
  if (inherits(predictions, "simitall_genomic_model")) {
    predictions <- predict_genomic_values(predictions)
  } else {
    predictions <- .simitall_read_analysis_table(predictions, "Prediction")
  }
  if (any(!c("sample", "gebv") %in% names(predictions))) {
    stop("predictions must contain sample and gebv columns")
  }
  predictions <- predictions[is.finite(predictions$gebv), , drop = FALSE]
  if (!nrow(predictions)) stop("No candidates have finite GEBVs")
  if (is.null(n_parents)) {
    if (!is.numeric(proportion) || length(proportion) != 1L ||
        proportion <= 0 || proportion > 1) {
      stop("proportion must be one number in (0, 1]")
    }
    n_parents <- max(1L, ceiling(nrow(predictions) * proportion))
  }
  n_parents <- as.integer(n_parents)
  if (n_parents < 1L || n_parents > nrow(predictions)) {
    stop("n_parents must be between one and the number of candidates")
  }
  if (!is.finite(diversity_penalty) || diversity_penalty < 0) {
    stop("diversity_penalty must be non-negative")
  }
  if (!is.null(family_data)) {
    family_data <- .simitall_read_analysis_table(family_data, "Family metadata")
    if (any(!c(sample_id, family_col) %in% names(family_data))) {
      stop("family_data does not contain the requested columns")
    }
    predictions$family <- as.character(family_data[[family_col]][match(
      predictions$sample, family_data[[sample_id]]
    )])
  } else {
    predictions$family <- NA_character_
  }
  if (!is.null(kinship)) {
    kinship <- as.matrix(kinship)
    if (is.null(rownames(kinship)) || is.null(colnames(kinship))) {
      stop("kinship must have sample IDs as row and column names")
    }
    missing <- setdiff(predictions$sample, rownames(kinship))
    if (length(missing)) stop("kinship is missing candidate samples")
  }

  scale_value <- stats::sd(predictions$gebv)
  base_score <- if (!is.finite(scale_value) || scale_value == 0) {
    rep(0, nrow(predictions))
  } else {
    as.numeric(scale(predictions$gebv))
  }
  if (direction == "min") base_score <- -base_score
  selected_rows <- integer()
  selection_scores <- numeric()
  coancestry_values <- numeric()
  family_counts <- integer()

  while (length(selected_rows) < n_parents) {
    available <- setdiff(seq_len(nrow(predictions)), selected_rows)
    if (is.finite(max_per_family) && any(!is.na(predictions$family[available]))) {
      allowed <- vapply(available, function(index) {
        family <- predictions$family[index]
        if (is.na(family)) return(TRUE)
        current <- if (family %in% names(family_counts)) family_counts[[family]] else 0L
        current < max_per_family
      }, logical(1L))
      available <- available[allowed]
    }
    if (!length(available)) {
      stop("Family constraints prevent selection of n_parents candidates")
    }
    coancestry <- rep(0, length(available))
    if (!is.null(kinship) && length(selected_rows)) {
      selected_ids <- predictions$sample[selected_rows]
      coancestry <- vapply(available, function(index) {
        mean(kinship[predictions$sample[index], selected_ids], na.rm = TRUE)
      }, numeric(1L))
    }
    score <- base_score[available] - diversity_penalty * coancestry
    chosen_position <- which.max(score)
    chosen <- available[chosen_position]
    selected_rows <- c(selected_rows, chosen)
    selection_scores <- c(selection_scores, score[chosen_position])
    coancestry_values <- c(coancestry_values, coancestry[chosen_position])
    family <- predictions$family[chosen]
    if (!is.na(family)) {
      if (!family %in% names(family_counts)) family_counts[family] <- 0L
      family_counts[family] <- family_counts[family] + 1L
    }
  }
  selected <- predictions[selected_rows, , drop = FALSE]
  selected$selection_rank <- seq_len(nrow(selected))
  selected$selection_score <- selection_scores
  selected$mean_kinship_to_previous <- coancestry_values
  rownames(selected) <- NULL
  selected
}

#' Design crosses among selected parents
#'
#' Construct unique parent pairs using random, assortative, disassortative, or
#' minimum-kinship mating. This function creates a crossing plan; offspring are
#' generated separately by the breeding simulator.
#'
#' @param selected_parents Table returned by [select_parents()] or compatible
#'   data frame/TSV containing `sample` and `gebv`.
#' @param n_crosses Number of unique crosses. Defaults to all parent pairs.
#' @param mating Pairing strategy.
#' @param kinship Named relationship matrix required for
#'   `"minimum_kinship"`.
#' @param seed Random-number seed used for random mating.
#' @param out_file Optional crossing-plan TSV path.
#'
#' @return A crossing-plan data frame.
#' @examples
#' \dontrun{
#' crosses <- design_parent_crosses(
#'   selected,
#'   n_crosses = 20,
#'   mating = "minimum_kinship",
#'   kinship = K
#' )
#' }
#' @export
design_parent_crosses <- function(
    selected_parents,
    n_crosses = NULL,
    mating = c("random", "minimum_kinship", "assortative", "disassortative"),
    kinship = NULL,
    seed = 1L,
    out_file = NULL) {
  mating <- match.arg(mating)
  parents <- .simitall_read_analysis_table(selected_parents, "Selected parent")
  if (any(!c("sample", "gebv") %in% names(parents))) {
    stop("selected_parents must contain sample and gebv columns")
  }
  if (nrow(parents) < 2L) stop("At least two selected parents are required")
  pair_index <- utils::combn(seq_len(nrow(parents)), 2L)
  crosses <- data.frame(
    parent1 = parents$sample[pair_index[1L, ]],
    parent2 = parents$sample[pair_index[2L, ]],
    parent1_gebv = parents$gebv[pair_index[1L, ]],
    parent2_gebv = parents$gebv[pair_index[2L, ]],
    stringsAsFactors = FALSE
  )
  crosses$parent_mean_gebv <- rowMeans(crosses[, c(
    "parent1_gebv", "parent2_gebv"
  )])
  crosses$gebv_difference <- abs(
    crosses$parent1_gebv - crosses$parent2_gebv
  )
  crosses$kinship <- NA_real_
  if (!is.null(kinship)) {
    kinship <- as.matrix(kinship)
    missing <- setdiff(parents$sample, rownames(kinship))
    if (length(missing)) stop("kinship is missing selected parents")
    crosses$kinship <- vapply(seq_len(nrow(crosses)), function(i) {
      kinship[crosses$parent1[i], crosses$parent2[i]]
    }, numeric(1L))
  }
  if (mating == "minimum_kinship" && all(is.na(crosses$kinship))) {
    stop("minimum_kinship mating requires a named kinship matrix")
  }
  ordering <- switch(
    mating,
    random = {
      set.seed(as.integer(seed))
      sample(seq_len(nrow(crosses)))
    },
    minimum_kinship = order(
      crosses$kinship, -crosses$parent_mean_gebv, na.last = TRUE
    ),
    assortative = order(
      crosses$gebv_difference, -crosses$parent_mean_gebv
    ),
    disassortative = order(
      -crosses$gebv_difference, -crosses$parent_mean_gebv
    )
  )
  crosses <- crosses[ordering, , drop = FALSE]
  if (is.null(n_crosses)) n_crosses <- nrow(crosses)
  n_crosses <- as.integer(n_crosses)
  if (n_crosses < 1L || n_crosses > nrow(crosses)) {
    stop("n_crosses must be between one and the number of unique parent pairs")
  }
  crosses <- crosses[seq_len(n_crosses), , drop = FALSE]
  crosses$cross_id <- paste0("cross", seq_len(nrow(crosses)))
  crosses$mating_strategy <- mating
  crosses <- crosses[, c(
    "cross_id", "parent1", "parent2", "parent1_gebv", "parent2_gebv",
    "parent_mean_gebv", "gebv_difference", "kinship", "mating_strategy"
  )]
  rownames(crosses) <- NULL
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      crosses, out_file, sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  crosses
}

#' Run one genomic-selection round
#'
#' Fit a genomic model, predict candidates, select parents, and produce a
#' crossing plan. This function completes one decision round and records all
#' artifacts. It does not yet generate recombinant offspring automatically;
#' the crossing plan can be passed to a breeding workflow in the next step.
#'
#' @inheritParams fit_genomic_model
#' @param out_prefix Prefix for model, prediction, selected-parent, crossing,
#'   kinship, and summary outputs.
#' @param candidate_samples Optional IDs eligible for selection. The default is
#'   every genotyped sample.
#' @param n_parents Number of parents selected.
#' @param n_crosses Number of unique crosses designed.
#' @param mating Crossing strategy passed to [design_parent_crosses()].
#' @param diversity_penalty Kinship penalty passed to [select_parents()].
#' @param family_data Optional sample-family metadata.
#' @param family_col Family column in `family_data`.
#' @param max_per_family Maximum selected parents per family.
#' @param ... Additional arguments passed to [fit_genomic_model()].
#'
#' @return Invisibly returns the fitted model, prediction table, selected
#'   parents, crossing plan, kinship matrix, summary, and output paths.
#' @examples
#' \dontrun{
#' selection <- run_genomic_selection(
#'   "results/gwas/demo.geno.tsv",
#'   "results/gwas/demo.pheno.tsv",
#'   out_prefix = "results/selection/cycle1",
#'   model = "gblup",
#'   n_parents = 10,
#'   n_crosses = 20,
#'   diversity_penalty = 0.25
#' )
#' }
#' @export
run_genomic_selection <- function(
    genotype_file,
    phenotype,
    out_prefix,
    trait = "trait",
    sample_id = "sample",
    model = "gblup",
    candidate_samples = NULL,
    n_parents = 10L,
    n_crosses = 20L,
    mating = "minimum_kinship",
    diversity_penalty = 0,
    family_data = NULL,
    family_col = "family",
    max_per_family = Inf,
    seed = 1L,
    ...) {
  n_parents <- as.integer(n_parents)
  n_crosses <- as.integer(n_crosses)
  if (length(n_parents) != 1L || is.na(n_parents) || n_parents < 2L) {
    stop("n_parents must be at least two so a crossing plan can be made")
  }
  if (length(n_crosses) != 1L || is.na(n_crosses) || n_crosses < 1L) {
    stop("n_crosses must be at least one")
  }
  fit <- fit_genomic_model(
    genotype_file = genotype_file,
    phenotype = phenotype,
    trait = trait,
    sample_id = sample_id,
    model = model,
    seed = seed,
    ...
  )
  predictions <- predict_genomic_values(fit)
  if (!is.null(candidate_samples)) {
    missing <- setdiff(candidate_samples, predictions$sample)
    if (length(missing)) stop("candidate_samples are missing from predictions")
    predictions <- predictions[
      match(candidate_samples, predictions$sample),
      ,
      drop = FALSE
    ]
  }
  kinship <- calculate_genomic_kinship(genotype_file)
  selected <- select_parents(
    predictions = predictions,
    n_parents = n_parents,
    kinship = kinship,
    diversity_penalty = diversity_penalty,
    family_data = family_data,
    sample_id = sample_id,
    family_col = family_col,
    max_per_family = max_per_family
  )
  max_crosses <- choose(nrow(selected), 2L)
  n_crosses <- min(as.integer(n_crosses), max_crosses)
  crosses <- design_parent_crosses(
    selected,
    n_crosses = n_crosses,
    mating = mating,
    kinship = kinship,
    seed = seed
  )
  paths <- list(
    model = paste0(out_prefix, ".model.rds"),
    predictions = paste0(out_prefix, ".predictions.tsv"),
    selected = paste0(out_prefix, ".selected_parents.tsv"),
    crosses = paste0(out_prefix, ".crosses.tsv"),
    kinship = paste0(out_prefix, ".kinship.rds"),
    summary = paste0(out_prefix, ".summary.json")
  )
  dir.create(dirname(paths$model), recursive = TRUE, showWarnings = FALSE)
  saveRDS(fit, paths$model)
  saveRDS(kinship, paths$kinship)
  utils::write.table(
    predictions, paths$predictions, sep = "\t", quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    selected, paths$selected, sep = "\t", quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    crosses, paths$crosses, sep = "\t", quote = FALSE, row.names = FALSE
  )
  summary <- list(
    model = model,
    trait = trait,
    training_samples = length(fit$training_samples),
    candidate_samples = nrow(predictions),
    selected_parents = nrow(selected),
    crosses = nrow(crosses),
    mating = mating,
    diversity_penalty = diversity_penalty,
    mean_candidate_gebv = mean(predictions$gebv),
    mean_selected_gebv = mean(selected$gebv),
    selection_differential = mean(selected$gebv) - mean(predictions$gebv),
    seed = as.integer(seed)
  )
  .simitall_write_json(summary, paths$summary)
  invisible(list(
    model = fit,
    predictions = predictions,
    selected = selected,
    crosses = crosses,
    kinship = kinship,
    summary = summary,
    paths = paths
  ))
}
