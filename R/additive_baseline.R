fit_genotype_preprocessor <- function(genotype_train, max_markers = NULL) {
  if (!is.matrix(genotype_train) || !is.numeric(genotype_train)) {
    stop("genotype_train must be a numeric matrix.", call. = FALSE)
  }
  if (is.null(colnames(genotype_train))) {
    stop("genotype_train must have marker IDs.", call. = FALSE)
  }
  marker_means <- colMeans(genotype_train, na.rm = TRUE)
  callable <- is.finite(marker_means)
  if (!any(callable)) stop("No marker is callable in the training split.", call. = FALSE)
  filled <- genotype_train[, callable, drop = FALSE]
  marker_means <- marker_means[callable]
  missing <- which(is.na(filled), arr.ind = TRUE)
  if (nrow(missing)) filled[missing] <- marker_means[missing[, "col"]]
  marker_sd <- apply(filled, 2L, stats::sd)
  variable <- is.finite(marker_sd) & marker_sd > 0
  if (!any(variable)) stop("No variable marker remains in the training split.", call. = FALSE)

  variable_ids <- colnames(filled)[variable]
  selected_ids <- variable_ids
  if (!is.null(max_markers)) {
    max_markers <- as.integer(max_markers)
    if (length(max_markers) != 1L || is.na(max_markers) || max_markers < 1L) {
      stop("max_markers must be NULL or one positive integer.", call. = FALSE)
    }
    if (length(variable_ids) > max_markers) {
      variable_sd <- marker_sd[variable]
      selected_order <- order(-variable_sd, variable_ids)
      selected_ids <- variable_ids[selected_order[seq_len(max_markers)]]
    }
  }
  selected <- match(selected_ids, colnames(filled))

  list(
    marker_ids = selected_ids,
    center = marker_means[selected],
    scale = marker_sd[selected],
    n_input_markers = ncol(genotype_train),
    n_callable_markers = sum(callable),
    n_variable_markers = sum(variable),
    n_selected_markers = length(selected_ids),
    provenance = list(
      fitted_on_training_only = TRUE,
      imputation = "training marker mean",
      scaling = "training sample standard deviation",
      marker_filter = if (is.null(max_markers)) {
        "all training-variable markers"
      } else {
        paste0(
          "up to ", max_markers,
          " markers selected by training-only sample variance; ties by marker ID"
        )
      },
      phenotype_used_for_marker_filter = FALSE
    )
  )
}

apply_genotype_preprocessor <- function(genotype, preprocessor) {
  if (!is.matrix(genotype) || !is.numeric(genotype)) {
    stop("genotype must be a numeric matrix.", call. = FALSE)
  }
  marker_rows <- match(preprocessor$marker_ids, colnames(genotype))
  if (anyNA(marker_rows)) stop("Prediction genotype is missing fitted marker IDs.", call. = FALSE)
  aligned <- genotype[, marker_rows, drop = FALSE]
  colnames(aligned) <- preprocessor$marker_ids
  missing <- which(is.na(aligned), arr.ind = TRUE)
  if (nrow(missing)) aligned[missing] <- preprocessor$center[missing[, "col"]]
  standardized <- sweep(
    sweep(aligned, 2L, preprocessor$center, FUN = "-"),
    2L,
    preprocessor$scale,
    FUN = "/"
  )
  if (any(!is.finite(standardized))) {
    stop("Non-finite standardized genotype value.", call. = FALSE)
  }
  standardized
}

profile_gblup_reml_delta <- function(y, relationship) {
  if (!is.numeric(y) || anyNA(y) || any(!is.finite(y))) {
    stop("y must contain finite, non-missing values.", call. = FALSE)
  }
  n <- length(y)
  if (!is.matrix(relationship) || !identical(dim(relationship), c(n, n))) {
    stop("relationship must be a square matrix aligned to y.", call. = FALSE)
  }
  decomposition <- eigen((relationship + t(relationship)) / 2, symmetric = TRUE)
  tolerance <- max(1, max(abs(decomposition$values))) * 1e-8
  if (any(decomposition$values < -tolerance)) {
    stop("relationship is not positive semidefinite.", call. = FALSE)
  }
  eigenvalues <- pmax(decomposition$values, 0)
  rotated_y <- as.numeric(crossprod(decomposition$vectors, y))
  rotated_x <- as.numeric(crossprod(decomposition$vectors, rep(1, n)))
  residual_df <- n - 1L

  negative_reml <- function(log_delta) {
    delta <- exp(log_delta)
    inverse_variance <- 1 / (eigenvalues + delta)
    xt_vinv_x <- sum(rotated_x^2 * inverse_variance)
    if (!is.finite(xt_vinv_x) || xt_vinv_x <= 0) return(Inf)
    intercept <- sum(rotated_x * rotated_y * inverse_variance) / xt_vinv_x
    residual <- rotated_y - rotated_x * intercept
    quadratic <- sum(residual^2 * inverse_variance)
    if (!is.finite(quadratic) || quadratic <= 0) return(Inf)
    log_likelihood <- -0.5 * (
      residual_df * (log(2 * pi) + 1 + log(quadratic / residual_df)) +
        sum(log(eigenvalues + delta)) + log(xt_vinv_x)
    )
    -log_likelihood
  }
  optimum <- stats::optimize(
    negative_reml,
    interval = log(c(1e-5, 1e5)),
    tol = 1e-8
  )
  list(
    delta = exp(optimum$minimum),
    reml_log_likelihood = -optimum$objective,
    boundary = exp(optimum$minimum) < 1.01e-5 || exp(optimum$minimum) > 0.99e5
  )
}

fit_gblup <- function(z_train, y_train, delta = NULL) {
  if (!is.matrix(z_train) || !is.numeric(z_train) || nrow(z_train) != length(y_train)) {
    stop("z_train and y_train are not aligned.", call. = FALSE)
  }
  if (anyNA(z_train) || anyNA(y_train) || ncol(z_train) < 1L) {
    stop("GBLUP requires complete standardized training inputs.", call. = FALSE)
  }
  p <- ncol(z_train)
  relationship <- tcrossprod(z_train) / p
  if (is.null(delta)) {
    reml <- profile_gblup_reml_delta(y_train, relationship)
    delta <- reml$delta
    shrinkage_source <- "REML"
  } else {
    delta <- as.numeric(delta)
    if (length(delta) != 1L || !is.finite(delta) || delta <= 0) {
      stop("delta must be one positive number.", call. = FALSE)
    }
    reml <- NULL
    shrinkage_source <- "fixed"
  }
  covariance_shape <- relationship + diag(delta, nrow(relationship))
  inverse_one <- solve(covariance_shape, rep(1, nrow(relationship)))
  intercept <- sum(inverse_one * y_train) / sum(inverse_one)
  alpha <- solve(covariance_shape, y_train - intercept)
  marker_effects <- as.numeric(crossprod(z_train, alpha)) / p
  names(marker_effects) <- colnames(z_train)
  list(
    intercept = intercept,
    alpha = alpha,
    marker_effects = marker_effects,
    z_train = z_train,
    delta = delta,
    lambda_marker = p * delta,
    n_markers = p,
    shrinkage_source = shrinkage_source,
    reml = reml
  )
}

predict_gblup <- function(model, z_test) {
  if (!is.matrix(z_test) || ncol(z_test) != model$n_markers || anyNA(z_test)) {
    stop("z_test is incomplete or does not match the fitted markers.", call. = FALSE)
  }
  relationship_test_train <- tcrossprod(z_test, model$z_train) / model$n_markers
  dual_prediction <- as.numeric(model$intercept + relationship_test_train %*% model$alpha)
  marker_prediction <- as.numeric(model$intercept + z_test %*% model$marker_effects)
  if (!isTRUE(all.equal(dual_prediction, marker_prediction, tolerance = 1e-8))) {
    stop("GBLUP and marker RR-BLUP predictions are not numerically equivalent.", call. = FALSE)
  }
  dual_prediction
}

regression_metrics <- function(observed, predicted, true_genetic_value = NULL) {
  if (length(observed) != length(predicted) || !length(observed) ||
      anyNA(observed) || anyNA(predicted)) {
    stop("Observed and predicted values must be complete and aligned.", call. = FALSE)
  }
  residual <- observed - predicted
    sst <- sum((observed - mean(observed))^2)
    safe_correlation <- function(x, y, method = "pearson") {
      complete <- is.finite(x) & is.finite(y)
      x <- x[complete]
      y <- y[complete]
      if (length(x) < 2L) return(NA_real_)
      x_sd <- stats::sd(x)
      y_sd <- stats::sd(y)
      if (!is.finite(x_sd) || !is.finite(y_sd) || x_sd == 0 || y_sd == 0) {
        return(NA_real_)
      }
      suppressWarnings(stats::cor(x, y, method = method))
    }
  output <- data.frame(
    rmse = sqrt(mean(residual^2)),
    mae = mean(abs(residual)),
    r2 = if (sst > 0) 1 - sum(residual^2) / sst else NA_real_,
    pearson_r = safe_correlation(observed, predicted),
    spearman_r = safe_correlation(observed, predicted, method = "spearman"),
    stringsAsFactors = FALSE
  )
  output$true_genetic_r <- if (is.null(true_genetic_value)) {
    NA_real_
  } else {
    safe_correlation(true_genetic_value, predicted)
  }
  output
}

evaluate_gblup_outer_folds <- function(
    dataset,
    folds,
    phenotype,
    task_id,
    fixed_h2 = NULL) {
  validate_genomic_dataset(dataset)
  required_phenotype <- c("sample_id", "value")
  if (!is.data.frame(phenotype) || !all(required_phenotype %in% names(phenotype))) {
    stop("phenotype must contain sample_id and value.", call. = FALSE)
  }
  if (anyDuplicated(phenotype$sample_id)) {
    stop("Phenotype sample IDs must be unique within a task.", call. = FALSE)
  }
  phenotype <- phenotype[is.finite(phenotype$value), , drop = FALSE]
  combinations <- unique(folds[, c("scheme", "repeat_id", "fold"), drop = FALSE])
  predictions <- list()
  metrics <- list()
  row_index <- 0L
  modes <- "reml"
  if (!is.null(fixed_h2)) modes <- c(modes, "known_h2")

  for (combination_index in seq_len(nrow(combinations))) {
    current <- folds[
      folds$scheme == combinations$scheme[[combination_index]] &
        folds$repeat_id == combinations$repeat_id[[combination_index]],
      ,
      drop = FALSE
    ]
    split <- fold_split_ids(current, combinations$repeat_id[[combination_index]], combinations$fold[[combination_index]])
    train_ids <- intersect(split$train, phenotype$sample_id)
    test_ids <- intersect(split$test, phenotype$sample_id)
    if (length(train_ids) < 3L || length(test_ids) < 1L) {
      stop("A task has too few train/test observations in an outer fold.", call. = FALSE)
    }
    train_rows <- match(train_ids, rownames(dataset$genotype))
    test_rows <- match(test_ids, rownames(dataset$genotype))
    phenotype_train_rows <- match(train_ids, phenotype$sample_id)
    phenotype_test_rows <- match(test_ids, phenotype$sample_id)
    preprocessor <- fit_genotype_preprocessor(dataset$genotype[train_rows, , drop = FALSE])
    z_train <- apply_genotype_preprocessor(dataset$genotype[train_rows, , drop = FALSE], preprocessor)
    z_test <- apply_genotype_preprocessor(dataset$genotype[test_rows, , drop = FALSE], preprocessor)
    y_train <- phenotype$value[phenotype_train_rows]
    y_test <- phenotype$value[phenotype_test_rows]
    true_test <- if ("true_additive_value" %in% names(phenotype)) {
      phenotype$true_additive_value[phenotype_test_rows]
    } else {
      NULL
    }

    for (mode in modes) {
      delta <- if (mode == "known_h2") (1 - fixed_h2) / fixed_h2 else NULL
      model <- fit_gblup(z_train, y_train, delta = delta)
      predicted <- predict_gblup(model, z_test)
      row_index <- row_index + 1L
      material <- unique(dataset$samples$breeding_material_type[
        match(test_ids, dataset$samples$sample_id)
      ])
      predictions[[row_index]] <- data.frame(
        task_id = task_id,
        scheme = combinations$scheme[[combination_index]],
        repeat_id = combinations$repeat_id[[combination_index]],
        fold = combinations$fold[[combination_index]],
        shrinkage = mode,
        sample_id = test_ids,
        breeding_material_type = dataset$samples$breeding_material_type[
          match(test_ids, dataset$samples$sample_id)
        ],
        observed = y_test,
        predicted = predicted,
        true_genetic_value = if (is.null(true_test)) NA_real_ else true_test,
        stringsAsFactors = FALSE
      )
      current_metrics <- regression_metrics(y_test, predicted, true_test)
      metrics[[row_index]] <- data.frame(
        task_id = task_id,
        scheme = combinations$scheme[[combination_index]],
        repeat_id = combinations$repeat_id[[combination_index]],
        fold = combinations$fold[[combination_index]],
        shrinkage = mode,
        n_train = length(train_ids),
        n_test = length(test_ids),
        n_markers = model$n_markers,
        delta = model$delta,
        lambda_marker = model$lambda_marker,
        reml_boundary = isTRUE(model$reml$boundary),
        breeding_material_type = paste(sort(material), collapse = ";"),
        current_metrics,
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    predictions = do.call(rbind, predictions),
    metrics = do.call(rbind, metrics),
    provenance = list(
      outer_test_used_for_preprocessing = FALSE,
      preprocessing = "marker filtering, mean imputation, and scaling fitted inside each outer-training split",
      fixed_h2_rule = "delta=(1-h2)/h2; lambda_marker=p*delta",
      equivalence_checked = "every prediction checked between dual GBLUP and marker RR-BLUP"
    )
  )
}

select_development_simulation_tasks <- function(registry) {
  eligible <- registry[
    registry$meta_partition == "calibration" & abs(registry$target_h2 - 0.5) < 1e-12,
    ,
    drop = FALSE
  ]
  eligible <- eligible[order(eligible$architecture, eligible$replicate_id), , drop = FALSE]
  eligible[!duplicated(eligible$architecture), , drop = FALSE]
}

summarize_gblup_development <- function(predictions, metrics, selected_tasks) {
  groups <- unique(predictions[, c("task_id", "scheme", "shrinkage"), drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (index in seq_len(nrow(groups))) {
    current_predictions <- predictions[
      predictions$task_id == groups$task_id[[index]] &
        predictions$scheme == groups$scheme[[index]] &
        predictions$shrinkage == groups$shrinkage[[index]],
      ,
      drop = FALSE
    ]
    current_metrics <- metrics[
      metrics$task_id == groups$task_id[[index]] &
        metrics$scheme == groups$scheme[[index]] &
        metrics$shrinkage == groups$shrinkage[[index]],
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
      shrinkage = groups$shrinkage[[index]],
      n_oof = nrow(current_predictions),
      breeding_material_type = paste(
        sort(unique(current_predictions$breeding_material_type)),
        collapse = ";"
      ),
      mean_fold_rmse = mean(current_metrics$rmse),
      sd_fold_rmse = stats::sd(current_metrics$rmse),
      mean_delta = mean(current_metrics$delta),
      reml_boundary_folds = sum(current_metrics$reml_boundary),
      pooled,
      stringsAsFactors = FALSE
    )
  }
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

run_gblup_development_baseline <- function(dataset, outer_folds, simulation_panel) {
  selected_tasks <- select_development_simulation_tasks(simulation_panel$registry)
  if (!nrow(selected_tasks)) stop("No calibration h2=0.5 simulation task is available.", call. = FALSE)
  task_results <- vector("list", nrow(selected_tasks))
  for (task_index in seq_len(nrow(selected_tasks))) {
    task_id <- selected_tasks$task_id[[task_index]]
    phenotype <- simulation_panel$phenotypes[
      simulation_panel$phenotypes$task_id == task_id,
      ,
      drop = FALSE
    ]
    task_results[[task_index]] <- evaluate_gblup_outer_folds(
      dataset,
      outer_folds$folds,
      phenotype,
      task_id,
      fixed_h2 = selected_tasks$target_h2[[task_index]]
    )
  }
  predictions <- do.call(rbind, lapply(task_results, `[[`, "predictions"))
  metrics <- do.call(rbind, lapply(task_results, `[[`, "metrics"))
  list(
    selected_tasks = selected_tasks,
    predictions = predictions,
    metrics = metrics,
    summary = summarize_gblup_development(predictions, metrics, selected_tasks),
    provenance = list(
      profile = "one calibration task per architecture at h2=0.5",
      locked_benchmark_analyzed = FALSE,
      learners = c("GBLUP", "equivalent marker RR-BLUP"),
      shrinkage = c("REML", "known simulation h2")
    )
  )
}

write_gblup_development_outputs <- function(
    baseline,
    processed_path = file.path("data", "processed", "tassel_tutorial", "gblup_development.rds"),
    table_directory = file.path("results", "tables")) {
  dir.create(dirname(processed_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
  output_paths <- c(
    baseline = processed_path,
    selected_tasks = file.path(table_directory, "tassel_tutorial_gblup_selected_tasks.tsv"),
    metrics = file.path(table_directory, "tassel_tutorial_gblup_development_metrics.tsv"),
    summary = file.path(table_directory, "tassel_tutorial_gblup_development_summary.tsv")
  )
  saveRDS(baseline, output_paths[["baseline"]], version = 3)
  utils::write.table(
    baseline$selected_tasks,
    output_paths[["selected_tasks"]],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    baseline$metrics,
    output_paths[["metrics"]],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    baseline$summary,
    output_paths[["summary"]],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  unname(output_paths)
}
