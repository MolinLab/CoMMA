elastic_net_lambda_ratios <- function(elastic_config) {
  minimum <- as.numeric(elastic_config$lambda_ratio_min)
  points <- as.integer(elastic_config$lambda_ratio_points)
  if (!is.finite(minimum) || minimum <= 0 || minimum >= 1 || points < 3L) {
    stop("Invalid Elastic Net lambda-ratio configuration.", call. = FALSE)
  }
  exp(seq(log(1), log(minimum), length.out = points))
}

elastic_net_lambda_max <- function(z, y, alpha) {
  if (!is.matrix(z) || nrow(z) != length(y) || anyNA(z) || anyNA(y)) {
    stop("Elastic Net lambda calculation requires aligned complete inputs.", call. = FALSE)
  }
  effective_alpha <- max(as.numeric(alpha), 0.001)
  centered_y <- y - mean(y)
  gradient <- abs(as.numeric(crossprod(z, centered_y))) / length(y)
  maximum <- max(gradient) / effective_alpha
  if (!is.finite(maximum) || maximum <= 0) {
    stop("Elastic Net lambda maximum is not positive and finite.", call. = FALSE)
  }
  maximum
}

fit_elastic_net_path <- function(z, y, alpha, lambda, max_iterations = 100000L) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for the Elastic Net baseline.", call. = FALSE)
  }
  glmnet::glmnet(
    x = z,
    y = y,
    family = "gaussian",
    alpha = alpha,
    lambda = lambda,
    intercept = TRUE,
    standardize = FALSE,
    control = list(
      thresh = 1e-8,
      maxit = as.integer(max_iterations)
    )
  )
}

tune_elastic_net_inner <- function(
    dataset,
    phenotype,
    inner_assignment,
    elastic_config,
    relative_loss_floor = 0.01) {
  alphas <- sort(unique(as.numeric(unlist(elastic_config$development_alpha))))
  if (!length(alphas) || any(alphas < 0 | alphas > 1)) {
    stop("Elastic Net alpha candidates must lie in [0, 1].", call. = FALSE)
  }
  lambda_ratios <- elastic_net_lambda_ratios(elastic_config)
  candidates <- expand.grid(
    alpha = alphas,
    lambda_ratio = lambda_ratios,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  inner_levels <- sort(unique(inner_assignment$inner_fold))
  fold_rmse <- matrix(
    NA_real_,
    nrow = nrow(candidates),
    ncol = length(inner_levels)
  )

  for (inner_index in seq_along(inner_levels)) {
    validation_ids <- intersect(
      inner_assignment$sample_id[inner_assignment$inner_fold == inner_levels[[inner_index]]],
      phenotype$sample_id
    )
    training_ids <- intersect(
      inner_assignment$sample_id[inner_assignment$inner_fold != inner_levels[[inner_index]]],
      phenotype$sample_id
    )
    if (length(training_ids) < 3L || length(validation_ids) < 1L) {
      stop("An inner split has too few phenotype observations.", call. = FALSE)
    }
    training_rows <- match(training_ids, rownames(dataset$genotype))
    validation_rows <- match(validation_ids, rownames(dataset$genotype))
    phenotype_training_rows <- match(training_ids, phenotype$sample_id)
    phenotype_validation_rows <- match(validation_ids, phenotype$sample_id)
    preprocessor <- fit_genotype_preprocessor(
      dataset$genotype[training_rows, , drop = FALSE]
    )
    z_training <- apply_genotype_preprocessor(
      dataset$genotype[training_rows, , drop = FALSE],
      preprocessor
    )
    z_validation <- apply_genotype_preprocessor(
      dataset$genotype[validation_rows, , drop = FALSE],
      preprocessor
    )
    y_training <- phenotype$value[phenotype_training_rows]
    y_validation <- phenotype$value[phenotype_validation_rows]

    for (alpha in alphas) {
      lambda_maximum <- elastic_net_lambda_max(z_training, y_training, alpha)
      lambda <- lambda_maximum * lambda_ratios
      fit <- fit_elastic_net_path(
        z_training,
        y_training,
        alpha,
        lambda,
        max_iterations = elastic_config$max_iterations
      )
      predicted <- as.matrix(stats::predict(fit, newx = z_validation, s = lambda))
      if (!identical(dim(predicted), c(length(validation_ids), length(lambda)))) {
        stop("Elastic Net prediction path has unexpected dimensions.", call. = FALSE)
      }
      alpha_rows <- which(candidates$alpha == alpha)
      fold_rmse[alpha_rows, inner_index] <- sqrt(colMeans(
        (predicted - matrix(y_validation, nrow = length(y_validation), ncol = length(lambda)))^2
      ))
    }
  }

  candidates$mean_rmse <- rowMeans(fold_rmse)
  candidates$se_rmse <- apply(fold_rmse, 1L, stats::sd) / sqrt(ncol(fold_rmse))
  best_index <- order(candidates$mean_rmse, -candidates$lambda_ratio, -candidates$alpha)[[1L]]
  tolerance <- max(
    candidates$se_rmse[[best_index]],
    as.numeric(relative_loss_floor) * candidates$mean_rmse[[best_index]]
  )
  threshold <- candidates$mean_rmse[[best_index]] + tolerance
  eligible <- which(candidates$mean_rmse <= threshold)
  selected_index <- eligible[order(
    -candidates$lambda_ratio[eligible],
    -candidates$alpha[eligible],
    candidates$mean_rmse[eligible]
  )[[1L]]]
  candidates$best <- seq_len(nrow(candidates)) == best_index
  candidates$within_one_se <- candidates$mean_rmse <= threshold
  candidates$selected <- seq_len(nrow(candidates)) == selected_index
  candidates$selection_threshold <- threshold
  candidates$inner_folds <- length(inner_levels)

  list(
    selected = candidates[selected_index, , drop = FALSE],
    archive = candidates,
    provenance = list(
      preprocessing_fitted_inside_each_inner_training_split = TRUE,
      lambda_scale = "fold-specific lambda_max times a common ratio grid",
      alpha_zero_lambda_max_surrogate = 0.001,
      selection_rule = "one-standard-error with 1% relative-loss floor; prefer stronger normalized penalty then higher alpha"
    )
  )
}

evaluate_elastic_net_outer_folds <- function(
    dataset,
    folds,
    inner_folds,
    phenotype,
    task_id,
    tuning_config) {
  validate_genomic_dataset(dataset)
  if (anyDuplicated(phenotype$sample_id)) {
    stop("Phenotype sample IDs must be unique within a task.", call. = FALSE)
  }
  phenotype <- phenotype[is.finite(phenotype$value), , drop = FALSE]
  combinations <- unique(folds[, c("scheme", "repeat_id", "fold"), drop = FALSE])
  predictions <- list()
  metrics <- list()
  selected_rows <- list()
  archive_rows <- list()
  coefficient_rows <- list()

  for (combination_index in seq_len(nrow(combinations))) {
    scheme <- combinations$scheme[[combination_index]]
    repeat_id <- combinations$repeat_id[[combination_index]]
    outer_fold <- combinations$fold[[combination_index]]
    current <- folds[
      folds$scheme == scheme & folds$repeat_id == repeat_id,
      ,
      drop = FALSE
    ]
    split <- fold_split_ids(current, repeat_id, outer_fold)
    train_ids <- intersect(split$train, phenotype$sample_id)
    test_ids <- intersect(split$test, phenotype$sample_id)
    current_inner <- inner_folds[
      inner_folds$outer_scheme == scheme &
        inner_folds$outer_repeat_id == repeat_id &
        inner_folds$outer_fold == outer_fold,
      ,
      drop = FALSE
    ]
    tuning <- tune_elastic_net_inner(
      dataset,
      phenotype[phenotype$sample_id %in% train_ids, , drop = FALSE],
      current_inner,
      tuning_config$elastic_net,
      relative_loss_floor = tuning_config$near_optimal$relative_loss_floor
    )

    train_rows <- match(train_ids, rownames(dataset$genotype))
    test_rows <- match(test_ids, rownames(dataset$genotype))
    phenotype_train_rows <- match(train_ids, phenotype$sample_id)
    phenotype_test_rows <- match(test_ids, phenotype$sample_id)
    preprocessor <- fit_genotype_preprocessor(dataset$genotype[train_rows, , drop = FALSE])
    z_train <- apply_genotype_preprocessor(
      dataset$genotype[train_rows, , drop = FALSE], preprocessor
    )
    z_test <- apply_genotype_preprocessor(
      dataset$genotype[test_rows, , drop = FALSE], preprocessor
    )
    y_train <- phenotype$value[phenotype_train_rows]
    y_test <- phenotype$value[phenotype_test_rows]
    true_test <- if ("true_additive_value" %in% names(phenotype)) {
      phenotype$true_additive_value[phenotype_test_rows]
    } else {
      NULL
    }
    selected_alpha <- tuning$selected$alpha[[1L]]
    selected_ratio <- tuning$selected$lambda_ratio[[1L]]
    lambda_maximum <- elastic_net_lambda_max(z_train, y_train, selected_alpha)
    selected_lambda <- lambda_maximum * selected_ratio
    fit <- fit_elastic_net_path(
      z_train,
      y_train,
      selected_alpha,
      selected_lambda,
      max_iterations = tuning_config$elastic_net$max_iterations
    )
    predicted <- as.numeric(stats::predict(fit, newx = z_test, s = selected_lambda))
    coefficients <- as.matrix(stats::coef(fit, s = selected_lambda))
    marker_effects <- coefficients[-1L, 1L]
    nonzero <- marker_effects != 0
    current_metrics <- regression_metrics(y_test, predicted, true_test)
    material <- dataset$samples$breeding_material_type[
      match(test_ids, dataset$samples$sample_id)
    ]

    predictions[[combination_index]] <- data.frame(
      task_id = task_id,
      scheme = scheme,
      repeat_id = repeat_id,
      fold = outer_fold,
      sample_id = test_ids,
      breeding_material_type = material,
      observed = y_test,
      predicted = predicted,
      true_genetic_value = if (is.null(true_test)) NA_real_ else true_test,
      stringsAsFactors = FALSE
    )
    metrics[[combination_index]] <- data.frame(
      task_id = task_id,
      scheme = scheme,
      repeat_id = repeat_id,
      fold = outer_fold,
      n_train = length(train_ids),
      n_test = length(test_ids),
      n_markers = ncol(z_train),
      alpha = selected_alpha,
      lambda_ratio = selected_ratio,
      lambda = selected_lambda,
      n_nonzero = sum(nonzero),
      breeding_material_type = paste(sort(unique(material)), collapse = ";"),
      current_metrics,
      stringsAsFactors = FALSE
    )
    selected_rows[[combination_index]] <- data.frame(
      task_id = task_id,
      scheme = scheme,
      repeat_id = repeat_id,
      fold = outer_fold,
      alpha = selected_alpha,
      lambda_ratio = selected_ratio,
      lambda = selected_lambda,
      inner_mean_rmse = tuning$selected$mean_rmse[[1L]],
      inner_se_rmse = tuning$selected$se_rmse[[1L]],
      selection_threshold = tuning$selected$selection_threshold[[1L]],
      n_nonzero = sum(nonzero),
      stringsAsFactors = FALSE
    )
    archive <- tuning$archive
    archive$task_id <- task_id
    archive$scheme <- scheme
    archive$repeat_id <- repeat_id
    archive$fold <- outer_fold
    archive_rows[[combination_index]] <- archive[, c(
      "task_id", "scheme", "repeat_id", "fold", "alpha", "lambda_ratio",
      "mean_rmse", "se_rmse", "best", "within_one_se", "selected",
      "selection_threshold", "inner_folds"
    )]
    n_selected_markers <- sum(nonzero)
    coefficient_rows[[combination_index]] <- data.frame(
      task_id = rep(task_id, n_selected_markers),
      scheme = rep(scheme, n_selected_markers),
      repeat_id = rep(repeat_id, n_selected_markers),
      fold = rep(outer_fold, n_selected_markers),
      marker_id = preprocessor$marker_ids[nonzero],
      standardized_effect = marker_effects[nonzero],
      alpha = rep(selected_alpha, n_selected_markers),
      lambda = rep(selected_lambda, n_selected_markers),
      stringsAsFactors = FALSE
    )
  }

  list(
    predictions = do.call(rbind, predictions),
    metrics = do.call(rbind, metrics),
    selected_hyperparameters = do.call(rbind, selected_rows),
    tuning_archive = do.call(rbind, archive_rows),
    nonzero_coefficients = do.call(rbind, coefficient_rows),
    provenance = list(
      outer_test_used_for_preprocessing_or_tuning = FALSE,
      inner_preprocessing_fitted_on_inner_training_only = TRUE,
      final_preprocessing_fitted_on_outer_training_only = TRUE,
      package = paste0("glmnet ", as.character(utils::packageVersion("glmnet")))
    )
  )
}

summarize_elastic_net_development <- function(predictions, metrics, selected_tasks) {
  groups <- unique(predictions[, c("task_id", "scheme"), drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (index in seq_len(nrow(groups))) {
    current_predictions <- predictions[
      predictions$task_id == groups$task_id[[index]] &
        predictions$scheme == groups$scheme[[index]],
      ,
      drop = FALSE
    ]
    current_metrics <- metrics[
      metrics$task_id == groups$task_id[[index]] &
        metrics$scheme == groups$scheme[[index]],
      ,
      drop = FALSE
    ]
    pooled <- regression_metrics(
      current_predictions$observed,
      current_predictions$predicted,
      current_predictions$true_genetic_value
    )
    task_row <- selected_tasks[selected_tasks$task_id == groups$task_id[[index]], , drop = FALSE]
    rows[[index]] <- data.frame(
      task_id = groups$task_id[[index]],
      architecture = task_row$architecture[[1L]],
      target_h2 = task_row$target_h2[[1L]],
      scheme = groups$scheme[[index]],
      n_oof = nrow(current_predictions),
      breeding_material_type = paste(
        sort(unique(current_predictions$breeding_material_type)), collapse = ";"
      ),
      mean_fold_rmse = mean(current_metrics$rmse),
      sd_fold_rmse = stats::sd(current_metrics$rmse),
      mean_alpha = mean(current_metrics$alpha),
      mean_lambda_ratio = mean(current_metrics$lambda_ratio),
      mean_nonzero = mean(current_metrics$n_nonzero),
      pooled,
      stringsAsFactors = FALSE
    )
  }
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

run_elastic_net_development_baseline <- function(
    dataset,
    outer_folds,
    inner_folds,
    simulation_panel,
    tuning_config) {
  selected_tasks <- select_development_simulation_tasks(simulation_panel$registry)
  task_results <- vector("list", nrow(selected_tasks))
  for (task_index in seq_len(nrow(selected_tasks))) {
    task_id <- selected_tasks$task_id[[task_index]]
    phenotype <- simulation_panel$phenotypes[
      simulation_panel$phenotypes$task_id == task_id,
      ,
      drop = FALSE
    ]
    task_results[[task_index]] <- evaluate_elastic_net_outer_folds(
      dataset,
      outer_folds$folds,
      inner_folds,
      phenotype,
      task_id,
      tuning_config
    )
  }
  predictions <- do.call(rbind, lapply(task_results, `[[`, "predictions"))
  metrics <- do.call(rbind, lapply(task_results, `[[`, "metrics"))
  list(
    selected_tasks = selected_tasks,
    predictions = predictions,
    metrics = metrics,
    selected_hyperparameters = do.call(rbind, lapply(task_results, `[[`, "selected_hyperparameters")),
    tuning_archive = do.call(rbind, lapply(task_results, `[[`, "tuning_archive")),
    nonzero_coefficients = do.call(rbind, lapply(task_results, `[[`, "nonzero_coefficients")),
    summary = summarize_elastic_net_development(predictions, metrics, selected_tasks),
    provenance = list(
      profile = "one calibration task per architecture at h2=0.5",
      locked_benchmark_analyzed = FALSE,
      stored_inner_folds_reused = TRUE
    )
  )
}

write_elastic_net_development_outputs <- function(
    baseline,
    processed_path = file.path("data", "processed", "tassel_tutorial", "elastic_net_development.rds"),
    table_directory = file.path("results", "tables")) {
  dir.create(dirname(processed_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
  output_paths <- c(
    baseline = processed_path,
    metrics = file.path(table_directory, "tassel_tutorial_elastic_net_development_metrics.tsv"),
    summary = file.path(table_directory, "tassel_tutorial_elastic_net_development_summary.tsv"),
    selected = file.path(table_directory, "tassel_tutorial_elastic_net_selected_hyperparameters.tsv"),
    tuning = file.path(table_directory, "tassel_tutorial_elastic_net_tuning_archive.tsv")
  )
  saveRDS(baseline, output_paths[["baseline"]], version = 3)
  tables <- list(
    metrics = baseline$metrics,
    summary = baseline$summary,
    selected = baseline$selected_hyperparameters,
    tuning = baseline$tuning_archive
  )
  for (table_name in names(tables)) {
    utils::write.table(
      tables[[table_name]],
      output_paths[[table_name]],
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  unname(output_paths)
}
