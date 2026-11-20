xgboost_candidate_design <- function(xgb_config, development_config, seed) {
  n <- as.integer(development_config$development_candidates %||% 9L)
  if (length(n) != 1L || is.na(n) || n < 3L) {
    stop("XGBoost development_candidates must be at least three.", call. = FALSE)
  }
  bounds <- function(name) {
    value <- xgb_config[[name]]
    c(lower = as.numeric(value$lower), upper = as.numeric(value$upper))
  }
  eta <- bounds("eta")
  depth <- bounds("max_depth")
  child <- bounds("min_child_weight")
  subsample <- bounds("subsample")
  colsample <- bounds("colsample_bytree")
  lambda <- bounds("lambda")
  alpha <- bounds("alpha")
  gamma <- bounds("gamma")
  all_bounds <- c(eta, depth, child, subsample, colsample, lambda, alpha, gamma)
  if (any(!is.finite(all_bounds)) || eta[[1L]] <= 0 || child[[1L]] <= 0 ||
      colsample[[1L]] <= 0 || lambda[[1L]] <= 0 || alpha[[1L]] <= 0 ||
      subsample[[1L]] <= 0 || subsample[[2L]] > 1 || colsample[[2L]] > 1) {
    stop("XGBoost tuning bounds are invalid.", call. = FALSE)
  }
  design <- with_seed_preserved(as.integer(seed), {
    stratified <- function() sample((seq_len(n) - stats::runif(n)) / n)
    log_value <- function(bound) exp(
      log(bound[[1L]]) + stratified() * (log(bound[[2L]]) - log(bound[[1L]]))
    )
    linear_value <- function(bound) {
      bound[[1L]] + stratified() * (bound[[2L]] - bound[[1L]])
    }
    data.frame(
      candidate_id = sprintf("xgb_%03d", seq_len(n)),
      eta = log_value(eta),
      max_depth = pmin(
        as.integer(depth[[2L]]),
        pmax(as.integer(depth[[1L]]), floor(linear_value(depth) + 0.5))
      ),
      min_child_weight = log_value(child),
      subsample = linear_value(subsample),
      colsample_bytree = log_value(colsample),
      lambda = log_value(lambda),
      alpha = log_value(alpha),
      gamma = linear_value(gamma),
      stringsAsFactors = FALSE
    )
  })
  design$max_depth <- as.integer(design$max_depth)
  design
}

xgboost_parameters <- function(candidate, development_config, seed) {
  list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
    booster = "gbtree",
    tree_method = as.character(development_config$tree_method %||% "hist"),
    eta = as.numeric(candidate$eta[[1L]]),
    max_depth = as.integer(candidate$max_depth[[1L]]),
    min_child_weight = as.numeric(candidate$min_child_weight[[1L]]),
    subsample = as.numeric(candidate$subsample[[1L]]),
    colsample_bytree = as.numeric(candidate$colsample_bytree[[1L]]),
    lambda = as.numeric(candidate$lambda[[1L]]),
    alpha = as.numeric(candidate$alpha[[1L]]),
    gamma = as.numeric(candidate$gamma[[1L]]),
    nthread = as.integer(development_config$internal_threads %||% 1L),
    seed = as.integer(seed),
    verbosity = 0L
  )
}

fit_xgboost_model <- function(x, y, candidate, nrounds, development_config, seed) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required for the XGBoost baseline.", call. = FALSE)
  }
  if (!is.matrix(x) || nrow(x) != length(y) || anyNA(x) || anyNA(y)) {
    stop("XGBoost inputs must be aligned and complete.", call. = FALSE)
  }
  dtrain <- xgboost::xgb.DMatrix(data = x, label = as.numeric(y))
  timing <- system.time({
    fit <- xgboost::xgb.train(
      params = xgboost_parameters(candidate, development_config, seed),
      data = dtrain,
      nrounds = as.integer(nrounds),
      verbose = 0L
    )
  })
  list(fit = fit, timing = timing)
}

predict_xgboost_model <- function(model, x) {
  prediction <- as.numeric(stats::predict(
    model,
    newdata = xgboost::xgb.DMatrix(data = x)
  ))
  if (length(prediction) != nrow(x) || any(!is.finite(prediction))) {
    stop("XGBoost produced invalid predictions.", call. = FALSE)
  }
  prediction
}

xgboost_shap_contributions <- function(model, x, predictions, tolerance = 1e-5) {
  contributions <- stats::predict(
    model,
    newdata = xgboost::xgb.DMatrix(data = x),
    predcontrib = TRUE,
    approxcontrib = FALSE
  )
  contributions <- as.matrix(contributions)
  if (!identical(dim(contributions), c(nrow(x), ncol(x) + 1L))) {
    stop("XGBoost TreeSHAP returned unexpected dimensions.", call. = FALSE)
  }
  marker_contributions <- contributions[, seq_len(ncol(x)), drop = FALSE]
  colnames(marker_contributions) <- colnames(x)
  error <- max(abs(rowSums(contributions) - predictions))
  if (!is.finite(error) || error > tolerance) {
    stop("XGBoost TreeSHAP additivity check failed.", call. = FALSE)
  }
  list(marker_contributions = marker_contributions, max_abs_additivity_error = error)
}

evaluate_xgboost_candidate_inner <- function(
    dataset,
    phenotype,
    inner_assignment,
    candidate,
    nrounds,
    development_config,
    seed) {
  inner_levels <- sort(unique(inner_assignment$inner_fold))
  rows <- vector("list", length(inner_levels))
  for (inner_index in seq_along(inner_levels)) {
    validation_ids <- intersect(
      inner_assignment$sample_id[
        inner_assignment$inner_fold == inner_levels[[inner_index]]
      ],
      phenotype$sample_id
    )
    training_ids <- intersect(
      inner_assignment$sample_id[
        inner_assignment$inner_fold != inner_levels[[inner_index]]
      ],
      phenotype$sample_id
    )
    if (length(training_ids) < 3L || length(validation_ids) < 1L) {
      stop("An XGBoost inner split has too few phenotype observations.", call. = FALSE)
    }
    training_rows <- match(training_ids, rownames(dataset$genotype))
    validation_rows <- match(validation_ids, rownames(dataset$genotype))
    y_training <- phenotype$value[match(training_ids, phenotype$sample_id)]
    y_validation <- phenotype$value[match(validation_ids, phenotype$sample_id)]
    preprocessor <- fit_genotype_preprocessor(
      dataset$genotype[training_rows, , drop = FALSE],
      max_markers = development_config$max_markers %||% NULL
    )
    x_training <- apply_genotype_preprocessor(
      dataset$genotype[training_rows, , drop = FALSE], preprocessor
    )
    x_validation <- apply_genotype_preprocessor(
      dataset$genotype[validation_rows, , drop = FALSE], preprocessor
    )
    fitted <- fit_xgboost_model(
      x_training,
      y_training,
      candidate,
      nrounds,
      development_config,
      seed = as.integer(seed + inner_index)
    )
    predicted <- predict_xgboost_model(fitted$fit, x_validation)
    rows[[inner_index]] <- data.frame(
      inner_fold = inner_levels[[inner_index]],
      rmse = sqrt(mean((y_validation - predicted)^2)),
      elapsed_seconds = fitted$timing[["elapsed"]],
      cpu_seconds = fitted$timing[["user.self"]] + fitted$timing[["sys.self"]],
      model_bytes = as.numeric(utils::object.size(fitted$fit)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

tune_xgboost_inner <- function(
    dataset,
    phenotype,
    inner_assignment,
    xgb_config,
    development_config,
    relative_loss_floor = 0.01,
    seed = 1L) {
  candidates <- xgboost_candidate_design(xgb_config, development_config, seed)
  fidelities <- sort(unique(as.integer(unlist(development_config$fidelities))))
  reduction <- as.integer(development_config$reduction_factor %||% 3L)
  finalists <- as.integer(development_config$finalists %||% 3L)
  if (!length(fidelities) || anyNA(fidelities) || any(fidelities < 1L) ||
      reduction < 2L || finalists < 1L || finalists > nrow(candidates)) {
    stop("XGBoost multi-fidelity configuration is invalid.", call. = FALSE)
  }
  active <- seq_len(nrow(candidates))
  archive_rows <- list()
  archive_index <- 0L
  final_stage_rows <- NULL

  for (stage in seq_along(fidelities)) {
    stage_rows <- vector("list", length(active))
    for (active_index in seq_along(active)) {
      candidate_index <- active[[active_index]]
      fold_values <- evaluate_xgboost_candidate_inner(
        dataset,
        phenotype,
        inner_assignment,
        candidates[candidate_index, , drop = FALSE],
        nrounds = fidelities[[stage]],
        development_config = development_config,
        seed = as.integer(seed + stage * 100000L + candidate_index * 100L)
      )
      candidate <- candidates[candidate_index, , drop = FALSE]
      candidate$stage <- stage
      candidate$nrounds <- fidelities[[stage]]
      candidate$mean_rmse <- mean(fold_values$rmse)
      candidate$se_rmse <- stats::sd(fold_values$rmse) / sqrt(nrow(fold_values))
      candidate$total_elapsed_seconds <- sum(fold_values$elapsed_seconds)
      candidate$total_cpu_seconds <- sum(fold_values$cpu_seconds)
      candidate$max_model_bytes <- max(fold_values$model_bytes)
      candidate$inner_folds <- nrow(fold_values)
      stage_rows[[active_index]] <- candidate
    }
    stage_table <- do.call(rbind, stage_rows)
    keep_n <- if (stage == length(fidelities)) {
      nrow(stage_table)
    } else {
      max(finalists, ceiling(nrow(stage_table) / reduction))
    }
    retained_order <- order(
      stage_table$mean_rmse,
      stage_table$total_elapsed_seconds,
      stage_table$candidate_id
    )
    retained <- retained_order[seq_len(min(keep_n, length(retained_order)))]
    stage_table$retained_for_next_stage <- seq_len(nrow(stage_table)) %in% retained
    for (row_index in seq_len(nrow(stage_table))) {
      archive_index <- archive_index + 1L
      archive_rows[[archive_index]] <- stage_table[row_index, , drop = FALSE]
    }
    active <- match(stage_table$candidate_id[retained], candidates$candidate_id)
    final_stage_rows <- stage_table
  }

  best_index <- order(
    final_stage_rows$mean_rmse,
    final_stage_rows$total_elapsed_seconds
  )[[1L]]
  tolerance <- max(
    final_stage_rows$se_rmse[[best_index]],
    as.numeric(relative_loss_floor) * final_stage_rows$mean_rmse[[best_index]]
  )
  threshold <- final_stage_rows$mean_rmse[[best_index]] + tolerance
  eligible <- which(final_stage_rows$mean_rmse <= threshold)
  selected_index <- eligible[order(
    final_stage_rows$max_depth[eligible],
    -final_stage_rows$min_child_weight[eligible],
    final_stage_rows$colsample_bytree[eligible],
    final_stage_rows$mean_rmse[eligible]
  )[[1L]]]
  archive <- do.call(rbind, archive_rows)
  archive$best_final <- archive$stage == length(fidelities) &
    archive$candidate_id == final_stage_rows$candidate_id[[best_index]]
  archive$within_one_se_final <- archive$stage == length(fidelities) &
    archive$mean_rmse <= threshold
  archive$selected <- archive$stage == length(fidelities) &
    archive$candidate_id == final_stage_rows$candidate_id[[selected_index]]
  archive$selection_threshold <- threshold
  selected <- final_stage_rows[selected_index, , drop = FALSE]
  selected$selection_threshold <- threshold
  list(
    selected = selected,
    archive = archive,
    provenance = list(
      preprocessing_fitted_inside_each_inner_training_split = TRUE,
      strategy = "successive halving with a one-standard-error final selection",
      candidate_design = "seeded stratified space-filling broad design",
      internal_threads = as.integer(development_config$internal_threads %||% 1L)
    )
  )
}

evaluate_xgboost_outer_folds <- function(
    dataset,
    folds,
    inner_folds,
    phenotype,
    task_id,
    xgb_config,
    development_config,
    relative_loss_floor = 0.01,
    seed = 1L,
    compute_shap = TRUE) {
  validate_genomic_dataset(dataset)
  phenotype <- phenotype[is.finite(phenotype$value), , drop = FALSE]
  if (anyDuplicated(phenotype$sample_id)) {
    stop("Phenotype sample IDs must be unique within an XGBoost task.", call. = FALSE)
  }
  combinations <- unique(folds[, c("scheme", "repeat_id", "fold"), drop = FALSE])
  predictions <- metrics <- selected_rows <- archive_rows <- shap_rows <-
    vector("list", nrow(combinations))

  for (combination_index in seq_len(nrow(combinations))) {
    scheme <- combinations$scheme[[combination_index]]
    repeat_id <- combinations$repeat_id[[combination_index]]
    outer_fold <- combinations$fold[[combination_index]]
    current <- folds[
      folds$scheme == scheme & folds$repeat_id == repeat_id, , drop = FALSE
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
    fold_seed <- as.integer(seed + combination_index * 10000L)
    tuning <- tune_xgboost_inner(
      dataset,
      phenotype[phenotype$sample_id %in% train_ids, , drop = FALSE],
      current_inner,
      xgb_config,
      development_config,
      relative_loss_floor,
      seed = fold_seed
    )
    train_rows <- match(train_ids, rownames(dataset$genotype))
    test_rows <- match(test_ids, rownames(dataset$genotype))
    y_train <- phenotype$value[match(train_ids, phenotype$sample_id)]
    y_test <- phenotype$value[match(test_ids, phenotype$sample_id)]
    true_test <- if ("true_additive_value" %in% names(phenotype)) {
      phenotype$true_additive_value[match(test_ids, phenotype$sample_id)]
    } else NULL
    preprocessor <- fit_genotype_preprocessor(
      dataset$genotype[train_rows, , drop = FALSE],
      max_markers = development_config$max_markers %||% NULL
    )
    x_train <- apply_genotype_preprocessor(
      dataset$genotype[train_rows, , drop = FALSE], preprocessor
    )
    x_test <- apply_genotype_preprocessor(
      dataset$genotype[test_rows, , drop = FALSE], preprocessor
    )
    fitted <- fit_xgboost_model(
      x_train,
      y_train,
      tuning$selected,
      nrounds = tuning$selected$nrounds[[1L]],
      development_config = development_config,
      seed = fold_seed + 9000L
    )
    predicted <- predict_xgboost_model(fitted$fit, x_test)
    shap <- NULL
    if (isTRUE(compute_shap)) {
      shap <- xgboost_shap_contributions(fitted$fit, x_test, predicted)
      fitted_mean_abs_shap <- colMeans(abs(shap$marker_contributions))
      fitted_mean_signed_shap <- colMeans(shap$marker_contributions)
      all_marker_ids <- colnames(dataset$genotype)
      mean_abs_shap <- stats::setNames(rep(0, length(all_marker_ids)), all_marker_ids)
      mean_signed_shap <- stats::setNames(rep(0, length(all_marker_ids)), all_marker_ids)
      mean_abs_shap[names(fitted_mean_abs_shap)] <- fitted_mean_abs_shap
      mean_signed_shap[names(fitted_mean_signed_shap)] <- fitted_mean_signed_shap
      shap_order <- order(-mean_abs_shap, names(mean_abs_shap))
      shap_rank <- integer(length(mean_abs_shap))
      shap_rank[shap_order] <- seq_along(shap_order)
    }
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
      n_markers = ncol(x_train),
      n_input_markers = ncol(dataset$genotype),
      nrounds = tuning$selected$nrounds[[1L]],
      eta = tuning$selected$eta[[1L]],
      max_depth = tuning$selected$max_depth[[1L]],
      min_child_weight = tuning$selected$min_child_weight[[1L]],
      subsample = tuning$selected$subsample[[1L]],
      colsample_bytree = tuning$selected$colsample_bytree[[1L]],
      lambda = tuning$selected$lambda[[1L]],
      alpha = tuning$selected$alpha[[1L]],
      gamma = tuning$selected$gamma[[1L]],
      elapsed_seconds = fitted$timing[["elapsed"]],
      cpu_seconds = fitted$timing[["user.self"]] + fitted$timing[["sys.self"]],
      model_bytes = as.numeric(utils::object.size(fitted$fit)),
      shap_additivity_max_abs_error = if (is.null(shap)) {
        NA_real_
      } else {
        shap$max_abs_additivity_error
      },
      breeding_material_type = paste(sort(unique(material)), collapse = ";"),
      current_metrics,
      stringsAsFactors = FALSE
    )
    selected <- tuning$selected
    selected$task_id <- task_id
    selected$scheme <- scheme
    selected$repeat_id <- repeat_id
    selected$fold <- outer_fold
    selected_rows[[combination_index]] <- selected
    archive <- tuning$archive
    archive$task_id <- task_id
    archive$scheme <- scheme
    archive$repeat_id <- repeat_id
    archive$fold <- outer_fold
    archive_rows[[combination_index]] <- archive
    shap_rows[[combination_index]] <- if (isTRUE(compute_shap)) {
      data.frame(
        task_id = task_id,
        scheme = scheme,
        repeat_id = repeat_id,
        fold = outer_fold,
        marker_id = names(mean_abs_shap),
        mean_abs_shap = as.numeric(mean_abs_shap),
        mean_signed_shap = as.numeric(mean_signed_shap),
        shap_rank = as.integer(shap_rank),
        n_held_out = length(test_ids),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        task_id = character(), scheme = character(), repeat_id = integer(),
        fold = integer(), marker_id = character(), mean_abs_shap = numeric(),
        mean_signed_shap = numeric(), shap_rank = integer(), n_held_out = integer(),
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    predictions = do.call(rbind, predictions),
    metrics = do.call(rbind, metrics),
    selected_hyperparameters = do.call(rbind, selected_rows),
    tuning_archive = do.call(rbind, archive_rows),
    oof_shap = do.call(rbind, shap_rows),
    provenance = list(
      outer_test_used_for_preprocessing_or_tuning = FALSE,
      TreeSHAP_computed_on_held_out_outer_samples = isTRUE(compute_shap),
      TreeSHAP_candidate_selection_for_interactions_performed = FALSE,
      marker_cap_fitted_on_training_only = development_config$max_markers %||% NULL,
      internal_threads = as.integer(development_config$internal_threads %||% 1L),
      package = paste0("xgboost ", as.character(utils::packageVersion("xgboost")))
    )
  )
}

summarize_xgboost_development <- function(predictions, metrics, selected_tasks) {
  groups <- unique(predictions[, c("task_id", "scheme"), drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (index in seq_len(nrow(groups))) {
    task_id <- groups$task_id[[index]]
    scheme <- groups$scheme[[index]]
    current_predictions <- predictions[
      predictions$task_id == task_id & predictions$scheme == scheme,
      , drop = FALSE
    ]
    current_metrics <- metrics[
      metrics$task_id == task_id & metrics$scheme == scheme,
      , drop = FALSE
    ]
    task_row <- selected_tasks[selected_tasks$task_id == task_id, , drop = FALSE]
    pooled <- regression_metrics(
      current_predictions$observed,
      current_predictions$predicted,
      current_predictions$true_genetic_value
    )
    rows[[index]] <- data.frame(
      task_id = task_id,
      architecture = task_row$architecture[[1L]],
      target_h2 = task_row$target_h2[[1L]],
      scheme = scheme,
      n_oof = nrow(current_predictions),
      breeding_material_type = paste(
        sort(unique(current_predictions$breeding_material_type)), collapse = ";"
      ),
      mean_fold_rmse = mean(current_metrics$rmse),
      sd_fold_rmse = stats::sd(current_metrics$rmse),
      total_final_fit_seconds = sum(current_metrics$elapsed_seconds),
      max_shap_additivity_error = max(current_metrics$shap_additivity_max_abs_error),
      pooled,
      stringsAsFactors = FALSE
    )
  }
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

summarize_xgboost_shap <- function(oof_shap, selected_tasks, causal_effects, top_k = 50L) {
  oof_shap$top_k <- oof_shap$shap_rank <= as.integer(top_k)
  groups <- unique(oof_shap[, c("task_id", "scheme", "marker_id"), drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (index in seq_len(nrow(groups))) {
    current <- oof_shap[
      oof_shap$task_id == groups$task_id[[index]] &
        oof_shap$scheme == groups$scheme[[index]] &
        oof_shap$marker_id == groups$marker_id[[index]],
      , drop = FALSE
    ]
    rows[[index]] <- data.frame(
      task_id = groups$task_id[[index]],
      scheme = groups$scheme[[index]],
      marker_id = groups$marker_id[[index]],
      weighted_mean_abs_shap = stats::weighted.mean(
        current$mean_abs_shap, current$n_held_out
      ),
      weighted_mean_signed_shap = stats::weighted.mean(
        current$mean_signed_shap, current$n_held_out
      ),
      mean_shap_rank = mean(current$shap_rank),
      top_k_frequency = mean(current$top_k),
      sign_positive_fraction = mean(current$mean_signed_shap > 0),
      n_outer_folds = nrow(current),
      stringsAsFactors = FALSE
    )
  }
  shap_summary <- do.call(rbind, rows)
  recovery_rows <- list()
  recovery_index <- 0L
  for (task_id in unique(oof_shap$task_id)) {
    task_row <- selected_tasks[selected_tasks$task_id == task_id, , drop = FALSE]
    causal_ids <- causal_effects$marker_id[
      causal_effects$base_simulation_id == task_row$base_simulation_id[[1L]]
    ]
    for (scheme in unique(oof_shap$scheme[oof_shap$task_id == task_id])) {
      current <- oof_shap[
        oof_shap$task_id == task_id & oof_shap$scheme == scheme,
        , drop = FALSE
      ]
      folds <- unique(current[, c("repeat_id", "fold"), drop = FALSE])
      for (fold_index in seq_len(nrow(folds))) {
        fold_values <- current[
          current$repeat_id == folds$repeat_id[[fold_index]] &
            current$fold == folds$fold[[fold_index]],
          , drop = FALSE
        ]
        causal_ranks <- fold_values$shap_rank[match(causal_ids, fold_values$marker_id)]
        for (k in c(10L, 50L, 100L)) {
          recovery_index <- recovery_index + 1L
          recovery_rows[[recovery_index]] <- data.frame(
            task_id = task_id,
            architecture = task_row$architecture[[1L]],
            scheme = scheme,
            repeat_id = folds$repeat_id[[fold_index]],
            fold = folds$fold[[fold_index]],
            top_k = k,
            exact_causal_recall = mean(causal_ranks <= k),
            median_causal_rank = stats::median(causal_ranks),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  list(shap_summary = shap_summary, causal_recovery = do.call(rbind, recovery_rows))
}

run_xgboost_development_baseline <- function(
    dataset,
    outer_folds,
    inner_folds,
    simulation_panel,
    xgb_config,
    development_config,
    relative_loss_floor = 0.01,
    seed = 1L) {
  selected_tasks <- select_development_simulation_tasks(simulation_panel$registry)
  task_results <- vector("list", nrow(selected_tasks))
  for (task_index in seq_len(nrow(selected_tasks))) {
    task_id <- selected_tasks$task_id[[task_index]]
    phenotype <- simulation_panel$phenotypes[
      simulation_panel$phenotypes$task_id == task_id, , drop = FALSE
    ]
    task_results[[task_index]] <- evaluate_xgboost_outer_folds(
      dataset,
      outer_folds$folds,
      inner_folds,
      phenotype,
      task_id,
      xgb_config,
      development_config,
      relative_loss_floor,
      seed = as.integer(seed + task_index * 1000000L)
    )
  }
  predictions <- do.call(rbind, lapply(task_results, `[[`, "predictions"))
  metrics <- do.call(rbind, lapply(task_results, `[[`, "metrics"))
  oof_shap <- do.call(rbind, lapply(task_results, `[[`, "oof_shap"))
  shap_output <- summarize_xgboost_shap(
    oof_shap,
    selected_tasks,
    simulation_panel$causal_effects,
    top_k = development_config$shap_top_k %||% 50L
  )
  list(
    selected_tasks = selected_tasks,
    predictions = predictions,
    metrics = metrics,
    selected_hyperparameters = do.call(
      rbind, lapply(task_results, `[[`, "selected_hyperparameters")
    ),
    tuning_archive = do.call(rbind, lapply(task_results, `[[`, "tuning_archive")),
    oof_shap = oof_shap,
    shap_summary = shap_output$shap_summary,
    shap_causal_recovery = shap_output$causal_recovery,
    summary = summarize_xgboost_development(predictions, metrics, selected_tasks),
    provenance = list(
      profile = "one calibration task per architecture at h2=0.5",
      locked_benchmark_analyzed = FALSE,
      stored_inner_folds_reused = TRUE,
      TreeSHAP_is_out_of_fold = TRUE,
      interaction_analysis_performed = FALSE,
      internal_threads = as.integer(development_config$internal_threads %||% 1L)
    )
  )
}

write_xgboost_development_figure <- function(
    baseline,
    path = file.path("results", "figures", "tassel_tutorial_xgboost_development.pdf")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(path, width = 10, height = 7, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (metric in c("true_genetic_r", "rmse")) {
    values <- xtabs(
      stats::as.formula(paste(metric, "~ architecture + scheme")),
      data = baseline$summary
    )
    graphics::barplot(
      values, beside = TRUE,
      col = c("#4472C4", "#70AD47", "#ED7D31"),
      ylab = metric,
      main = paste("XGBoost", metric, "by architecture and CV scheme"),
      legend.text = rownames(values),
      args.legend = list(x = "topright", bty = "n")
    )
  }
  graphics::boxplot(
    eta ~ scheme,
    data = baseline$selected_hyperparameters,
    col = c("#4472C4", "#70AD47"),
    ylab = "Selected learning rate", xlab = "Outer CV scheme",
    main = "Nested-CV XGBoost learning-rate choices"
  )
  for (task_id in unique(baseline$shap_summary$task_id)) {
    for (scheme in unique(baseline$shap_summary$scheme[
      baseline$shap_summary$task_id == task_id
    ])) {
      current <- baseline$shap_summary[
        baseline$shap_summary$task_id == task_id &
          baseline$shap_summary$scheme == scheme,
        , drop = FALSE
      ]
      current <- current[order(-current$weighted_mean_abs_shap, current$marker_id),
                           , drop = FALSE]
      current <- utils::head(current, 20L)
      graphics::par(mar = c(8, 4.5, 3, 1))
      graphics::barplot(
        current$weighted_mean_abs_shap,
        names.arg = current$marker_id,
        las = 2, cex.names = 0.65, col = "#ED7D31",
        ylab = "Weighted mean absolute held-out SHAP",
        main = paste(task_id, scheme, "top 20 OOF TreeSHAP markers")
      )
    }
  }
  invisible(path)
}

write_xgboost_development_outputs <- function(
    baseline,
    processed_path = file.path(
      "data", "processed", "tassel_tutorial", "xgboost_development.rds"
    ),
    table_directory = file.path("results", "tables"),
    figure_directory = file.path("results", "figures")) {
  dir.create(dirname(processed_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    rds = processed_path,
    metrics = file.path(table_directory, "tassel_tutorial_xgboost_development_metrics.tsv"),
    summary = file.path(table_directory, "tassel_tutorial_xgboost_development_summary.tsv"),
    selected = file.path(table_directory, "tassel_tutorial_xgboost_selected_hyperparameters.tsv"),
    tuning = file.path(table_directory, "tassel_tutorial_xgboost_tuning_archive.tsv"),
    shap = file.path(table_directory, "tassel_tutorial_xgboost_shap_summary.tsv"),
    recovery = file.path(table_directory, "tassel_tutorial_xgboost_shap_causal_recovery.tsv")
  )
  saveRDS(baseline, paths[["rds"]], version = 3)
  tables <- list(
    metrics = baseline$metrics,
    summary = baseline$summary,
    selected = baseline$selected_hyperparameters,
    tuning = baseline$tuning_archive,
    shap = baseline$shap_summary,
    recovery = baseline$shap_causal_recovery
  )
  for (name in names(tables)) {
    utils::write.table(
      tables[[name]], paths[[name]], sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  figure <- write_xgboost_development_figure(
    baseline,
    file.path(figure_directory, "tassel_tutorial_xgboost_development.pdf")
  )
  unname(c(paths, figure))
}
