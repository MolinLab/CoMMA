random_forest_candidate_design <- function(rf_config, seed) {
  n <- as.integer(rf_config$development_candidates %||% 12L)
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("Random Forest development_candidates must be at least two.", call. = FALSE)
  }
  mtry_lower <- as.numeric(rf_config$mtry_ratio$lower)
  mtry_upper <- as.numeric(rf_config$mtry_ratio$upper)
  node_lower <- as.integer(rf_config$min_node_size$lower)
  node_upper <- as.integer(rf_config$min_node_size$upper)
  sample_lower <- as.numeric(rf_config$sample_fraction$lower)
  sample_upper <- as.numeric(rf_config$sample_fraction$upper)
  depths <- as.integer(unlist(rf_config$max_depth$values))
  if (!is.finite(mtry_lower) || !is.finite(mtry_upper) || mtry_lower <= 0 ||
      mtry_upper > 1 || mtry_lower >= mtry_upper || node_lower < 1L ||
      node_upper < node_lower || sample_lower <= 0 || sample_upper > 1 ||
      sample_lower > sample_upper || !length(depths) || anyNA(depths) ||
      any(depths < 0L)) {
    stop("Random Forest tuning bounds are invalid.", call. = FALSE)
  }
  design <- with_seed_preserved(as.integer(seed), {
    stratified <- function() sample((seq_len(n) - stats::runif(n)) / n)
    mtry_unit <- stratified()
    sample_unit <- stratified()
    node_unit <- stratified()
    depth_values <- sample(rep(depths, length.out = n))
    data.frame(
      candidate_id = sprintf("rf_%03d", seq_len(n)),
      mtry_ratio = exp(
        log(mtry_lower) + mtry_unit * (log(mtry_upper) - log(mtry_lower))
      ),
      min_node_size = pmin(
        node_upper,
        pmax(node_lower, floor(node_lower + node_unit * (node_upper - node_lower + 1L)))
      ),
      sample_fraction = sample_lower + sample_unit * (sample_upper - sample_lower),
      max_depth = depth_values,
      stringsAsFactors = FALSE
    )
  })
  design$min_node_size <- as.integer(design$min_node_size)
  design$max_depth <- as.integer(design$max_depth)
  design
}

fit_random_forest_model <- function(
    x,
    y,
    candidate,
    num_trees,
    seed,
    importance = "none") {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required for the Random Forest baseline.", call. = FALSE)
  }
  if (!is.matrix(x) || nrow(x) != length(y) || anyNA(x) || anyNA(y)) {
    stop("Random Forest inputs must be aligned and complete.", call. = FALSE)
  }
  num_trees <- as.integer(num_trees)
  if (num_trees < 1L) stop("Random Forest num_trees must be positive.", call. = FALSE)
  mtry <- max(1L, min(ncol(x), as.integer(round(
    as.numeric(candidate$mtry_ratio[[1L]]) * ncol(x)
  ))))
  arguments <- list(
    x = as.data.frame(x, check.names = FALSE),
    y = as.numeric(y),
    num.trees = num_trees,
    mtry = mtry,
    min.node.size = as.integer(candidate$min_node_size[[1L]]),
    sample.fraction = as.numeric(candidate$sample_fraction[[1L]]),
    replace = FALSE,
    splitrule = "variance",
    importance = importance,
    write.forest = TRUE,
    num.threads = 1L,
    seed = as.integer(seed),
    verbose = FALSE
  )
  depth <- as.integer(candidate$max_depth[[1L]])
  if (depth > 0L) arguments$max.depth <- depth
  fit <- do.call(ranger::ranger, arguments)
  attr(fit, "realized_mtry") <- mtry
  fit
}

predict_random_forest_model <- function(model, x) {
  prediction <- stats::predict(
    model,
    data = as.data.frame(x, check.names = FALSE),
    num.threads = 1L
  )$predictions
  prediction <- as.numeric(prediction)
  if (length(prediction) != nrow(x) || any(!is.finite(prediction))) {
    stop("Random Forest produced invalid predictions.", call. = FALSE)
  }
  prediction
}

tune_random_forest_inner <- function(
    dataset,
    phenotype,
    inner_assignment,
    rf_config,
    relative_loss_floor = 0.01,
    seed = 1L) {
  candidates <- random_forest_candidate_design(rf_config, seed)
  inner_levels <- sort(unique(inner_assignment$inner_fold))
  fold_rmse <- matrix(NA_real_, nrow(candidates), length(inner_levels))
  fold_elapsed <- matrix(NA_real_, nrow(candidates), length(inner_levels))
  fold_cpu <- matrix(NA_real_, nrow(candidates), length(inner_levels))
  fold_model_bytes <- matrix(NA_real_, nrow(candidates), length(inner_levels))

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
      stop("A Random Forest inner split has too few phenotype observations.",
           call. = FALSE)
    }
    training_rows <- match(training_ids, rownames(dataset$genotype))
    validation_rows <- match(validation_ids, rownames(dataset$genotype))
    y_training <- phenotype$value[match(training_ids, phenotype$sample_id)]
    y_validation <- phenotype$value[match(validation_ids, phenotype$sample_id)]
    preprocessor <- fit_genotype_preprocessor(
      dataset$genotype[training_rows, , drop = FALSE],
      max_markers = rf_config$max_markers %||% NULL
    )
    x_training <- apply_genotype_preprocessor(
      dataset$genotype[training_rows, , drop = FALSE], preprocessor
    )
    x_validation <- apply_genotype_preprocessor(
      dataset$genotype[validation_rows, , drop = FALSE], preprocessor
    )

    for (candidate_index in seq_len(nrow(candidates))) {
      timing <- system.time({
        fit <- fit_random_forest_model(
          x_training,
          y_training,
          candidates[candidate_index, , drop = FALSE],
          num_trees = rf_config$tuning_num_trees %||% rf_config$num_trees,
          seed = as.integer(seed + 1000L * inner_index + candidate_index),
          importance = "none"
        )
        predicted <- predict_random_forest_model(fit, x_validation)
      })
      fold_rmse[candidate_index, inner_index] <- sqrt(mean(
        (y_validation - predicted)^2
      ))
      fold_elapsed[candidate_index, inner_index] <- timing[["elapsed"]]
      fold_cpu[candidate_index, inner_index] <- timing[["user.self"]] +
        timing[["sys.self"]]
      fold_model_bytes[candidate_index, inner_index] <- as.numeric(utils::object.size(fit))
    }
  }

  candidates$mean_rmse <- rowMeans(fold_rmse)
  candidates$se_rmse <- apply(fold_rmse, 1L, stats::sd) / sqrt(ncol(fold_rmse))
  candidates$total_elapsed_seconds <- rowSums(fold_elapsed)
  candidates$total_cpu_seconds <- rowSums(fold_cpu)
  candidates$max_model_bytes <- apply(fold_model_bytes, 1L, max)
  candidates$tuning_num_trees <- as.integer(
    rf_config$tuning_num_trees %||% rf_config$num_trees
  )
  best_index <- order(candidates$mean_rmse, candidates$total_elapsed_seconds)[[1L]]
  tolerance <- max(
    candidates$se_rmse[[best_index]],
    as.numeric(relative_loss_floor) * candidates$mean_rmse[[best_index]]
  )
  threshold <- candidates$mean_rmse[[best_index]] + tolerance
  eligible <- which(candidates$mean_rmse <= threshold)
  depth_complexity <- ifelse(candidates$max_depth == 0L, Inf, candidates$max_depth)
  selected_index <- eligible[order(
    candidates$mtry_ratio[eligible],
    -candidates$min_node_size[eligible],
    candidates$sample_fraction[eligible],
    depth_complexity[eligible],
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
      candidate_design = "seeded stratified space-filling random design",
      selection_rule = paste(
        "one-standard-error with 1% relative-loss floor; prefer lower mtry,",
        "larger nodes, lower sample fraction, then shallower finite depth"
      ),
      internal_threads = 1L
    )
  )
}

evaluate_random_forest_outer_folds <- function(
    dataset,
    folds,
    inner_folds,
    phenotype,
    task_id,
    tuning_config,
    seed = 1L,
    compute_importance = TRUE) {
  validate_genomic_dataset(dataset)
  if (anyDuplicated(phenotype$sample_id)) {
    stop("Phenotype sample IDs must be unique within a task.", call. = FALSE)
  }
  phenotype <- phenotype[is.finite(phenotype$value), , drop = FALSE]
  combinations <- unique(folds[, c("scheme", "repeat_id", "fold"), drop = FALSE])
  predictions <- metrics <- selected_rows <- archive_rows <- importance_rows <-
    vector("list", nrow(combinations))
  rf_config <- tuning_config$random_forest

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
    fold_seed <- as.integer(seed + combination_index * 10000L)
    tuning <- tune_random_forest_inner(
      dataset,
      phenotype[phenotype$sample_id %in% train_ids, , drop = FALSE],
      current_inner,
      rf_config,
      relative_loss_floor = tuning_config$near_optimal$relative_loss_floor,
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
      max_markers = rf_config$max_markers %||% NULL
    )
    x_train <- apply_genotype_preprocessor(
      dataset$genotype[train_rows, , drop = FALSE], preprocessor
    )
    x_test <- apply_genotype_preprocessor(
      dataset$genotype[test_rows, , drop = FALSE], preprocessor
    )
    final_timing <- system.time({
      fit <- fit_random_forest_model(
        x_train,
        y_train,
        tuning$selected,
        num_trees = rf_config$num_trees,
        seed = fold_seed + 9000L,
        importance = if (isTRUE(compute_importance)) "permutation" else "none"
      )
      predicted <- predict_random_forest_model(fit, x_test)
    })
    current_metrics <- regression_metrics(y_test, predicted, true_test)
    material <- dataset$samples$breeding_material_type[
      match(test_ids, dataset$samples$sample_id)
    ]
    realized_mtry <- attr(fit, "realized_mtry")
    importance <- fit$variable.importance
    if (isTRUE(compute_importance)) {
      if (is.null(importance) || !identical(names(importance), colnames(x_train))) {
        stop("Random Forest permutation importance is missing or misaligned.",
             call. = FALSE)
      }
      importance_rank <- rank(-importance, ties.method = "min")
    }

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
      num_trees = as.integer(rf_config$num_trees),
      mtry_ratio = tuning$selected$mtry_ratio[[1L]],
      mtry = realized_mtry,
      min_node_size = tuning$selected$min_node_size[[1L]],
      sample_fraction = tuning$selected$sample_fraction[[1L]],
      max_depth = tuning$selected$max_depth[[1L]],
      elapsed_seconds = final_timing[["elapsed"]],
      cpu_seconds = final_timing[["user.self"]] + final_timing[["sys.self"]],
      model_bytes = as.numeric(utils::object.size(fit)),
      breeding_material_type = paste(sort(unique(material)), collapse = ";"),
      current_metrics,
      stringsAsFactors = FALSE
    )
    selected_rows[[combination_index]] <- data.frame(
      task_id = task_id,
      scheme = scheme,
      repeat_id = repeat_id,
      fold = outer_fold,
      candidate_id = tuning$selected$candidate_id[[1L]],
      mtry_ratio = tuning$selected$mtry_ratio[[1L]],
      mtry = realized_mtry,
      min_node_size = tuning$selected$min_node_size[[1L]],
      sample_fraction = tuning$selected$sample_fraction[[1L]],
      max_depth = tuning$selected$max_depth[[1L]],
      inner_mean_rmse = tuning$selected$mean_rmse[[1L]],
      inner_se_rmse = tuning$selected$se_rmse[[1L]],
      selection_threshold = tuning$selected$selection_threshold[[1L]],
      stringsAsFactors = FALSE
    )
    archive <- tuning$archive
    archive$task_id <- task_id
    archive$scheme <- scheme
    archive$repeat_id <- repeat_id
    archive$fold <- outer_fold
    archive_rows[[combination_index]] <- archive[, c(
      "task_id", "scheme", "repeat_id", "fold", "candidate_id",
      "mtry_ratio", "min_node_size", "sample_fraction", "max_depth",
      "mean_rmse", "se_rmse", "best", "within_one_se", "selected",
      "selection_threshold", "inner_folds", "tuning_num_trees",
      "total_elapsed_seconds", "total_cpu_seconds", "max_model_bytes"
    )]
    importance_rows[[combination_index]] <- if (isTRUE(compute_importance)) {
      data.frame(
        task_id = task_id,
        scheme = scheme,
        repeat_id = repeat_id,
        fold = outer_fold,
        marker_id = names(importance),
        permutation_importance = as.numeric(importance),
        importance_rank = as.integer(importance_rank),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        task_id = character(), scheme = character(), repeat_id = integer(),
        fold = integer(), marker_id = character(),
        permutation_importance = numeric(), importance_rank = integer(),
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    predictions = do.call(rbind, predictions),
    metrics = do.call(rbind, metrics),
    selected_hyperparameters = do.call(rbind, selected_rows),
    tuning_archive = do.call(rbind, archive_rows),
    permutation_importance = do.call(rbind, importance_rows),
    provenance = list(
      outer_test_used_for_preprocessing_tuning_or_importance = FALSE,
      inner_preprocessing_fitted_on_inner_training_only = TRUE,
      final_preprocessing_fitted_on_outer_training_only = TRUE,
      importance_fitted_on_outer_training_only = isTRUE(compute_importance),
      marker_cap_fitted_on_training_only = rf_config$max_markers %||% NULL,
      internal_threads = 1L,
      package = paste0("ranger ", as.character(utils::packageVersion("ranger")))
    )
  )
}

summarize_random_forest_development <- function(predictions, metrics, selected_tasks) {
  groups <- unique(predictions[, c("task_id", "scheme"), drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (index in seq_len(nrow(groups))) {
    task_id <- groups$task_id[[index]]
    scheme <- groups$scheme[[index]]
    current_predictions <- predictions[
      predictions$task_id == task_id & predictions$scheme == scheme,
      ,
      drop = FALSE
    ]
    current_metrics <- metrics[
      metrics$task_id == task_id & metrics$scheme == scheme,
      ,
      drop = FALSE
    ]
    pooled <- regression_metrics(
      current_predictions$observed,
      current_predictions$predicted,
      current_predictions$true_genetic_value
    )
    task_row <- selected_tasks[selected_tasks$task_id == task_id, , drop = FALSE]
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
      mean_mtry_ratio = mean(current_metrics$mtry_ratio),
      mean_min_node_size = mean(current_metrics$min_node_size),
      mean_sample_fraction = mean(current_metrics$sample_fraction),
      total_final_fit_seconds = sum(current_metrics$elapsed_seconds),
      pooled,
      stringsAsFactors = FALSE
    )
  }
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

summarize_random_forest_importance <- function(
    importance,
    selected_tasks,
    causal_effects,
    top_k = 50L) {
  top_k <- as.integer(top_k)
  importance$top_k <- importance$importance_rank <= top_k
  groups <- unique(importance[, c("task_id", "scheme", "marker_id"), drop = FALSE])
  summary_rows <- vector("list", nrow(groups))
  for (index in seq_len(nrow(groups))) {
    current <- importance[
      importance$task_id == groups$task_id[[index]] &
        importance$scheme == groups$scheme[[index]] &
        importance$marker_id == groups$marker_id[[index]],
      ,
      drop = FALSE
    ]
    summary_rows[[index]] <- data.frame(
      task_id = groups$task_id[[index]],
      scheme = groups$scheme[[index]],
      marker_id = groups$marker_id[[index]],
      mean_permutation_importance = mean(current$permutation_importance),
      sd_permutation_importance = stats::sd(current$permutation_importance),
      mean_importance_rank = mean(current$importance_rank),
      top_k_frequency = mean(current$top_k),
      positive_importance_fraction = mean(current$permutation_importance > 0),
      n_outer_folds = nrow(current),
      stringsAsFactors = FALSE
    )
  }
  importance_summary <- do.call(rbind, summary_rows)
  recovery_rows <- list()
  recovery_index <- 0L
  for (task_id in unique(importance$task_id)) {
    task_row <- selected_tasks[selected_tasks$task_id == task_id, , drop = FALSE]
    causal_ids <- causal_effects$marker_id[
      causal_effects$base_simulation_id == task_row$base_simulation_id[[1L]]
    ]
    for (scheme in unique(importance$scheme[importance$task_id == task_id])) {
      current <- importance[
        importance$task_id == task_id & importance$scheme == scheme,
        ,
        drop = FALSE
      ]
      folds <- unique(current[, c("repeat_id", "fold"), drop = FALSE])
      for (fold_index in seq_len(nrow(folds))) {
        fold_values <- current[
          current$repeat_id == folds$repeat_id[[fold_index]] &
            current$fold == folds$fold[[fold_index]],
          ,
          drop = FALSE
        ]
        causal_ranks <- fold_values$importance_rank[
          match(causal_ids, fold_values$marker_id)
        ]
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
  list(
    importance_summary = importance_summary,
    causal_recovery = do.call(rbind, recovery_rows)
  )
}

run_random_forest_development_baseline <- function(
    dataset,
    outer_folds,
    inner_folds,
    simulation_panel,
    tuning_config,
    seed = 1L) {
  selected_tasks <- select_development_simulation_tasks(simulation_panel$registry)
  task_results <- vector("list", nrow(selected_tasks))
  for (task_index in seq_len(nrow(selected_tasks))) {
    task_id <- selected_tasks$task_id[[task_index]]
    phenotype <- simulation_panel$phenotypes[
      simulation_panel$phenotypes$task_id == task_id,
      ,
      drop = FALSE
    ]
    task_results[[task_index]] <- evaluate_random_forest_outer_folds(
      dataset,
      outer_folds$folds,
      inner_folds,
      phenotype,
      task_id,
      tuning_config,
      seed = as.integer(seed + task_index * 1000000L)
    )
  }
  predictions <- do.call(rbind, lapply(task_results, `[[`, "predictions"))
  metrics <- do.call(rbind, lapply(task_results, `[[`, "metrics"))
  importance <- do.call(rbind, lapply(task_results, `[[`, "permutation_importance"))
  importance_output <- summarize_random_forest_importance(
    importance,
    selected_tasks,
    simulation_panel$causal_effects,
    top_k = tuning_config$random_forest$importance_top_k %||% 50L
  )
  list(
    selected_tasks = selected_tasks,
    predictions = predictions,
    metrics = metrics,
    selected_hyperparameters = do.call(
      rbind, lapply(task_results, `[[`, "selected_hyperparameters")
    ),
    tuning_archive = do.call(rbind, lapply(task_results, `[[`, "tuning_archive")),
    permutation_importance = importance,
    importance_summary = importance_output$importance_summary,
    importance_causal_recovery = importance_output$causal_recovery,
    summary = summarize_random_forest_development(
      predictions, metrics, selected_tasks
    ),
    provenance = list(
      profile = "one calibration task per architecture at h2=0.5",
      locked_benchmark_analyzed = FALSE,
      stored_inner_folds_reused = TRUE,
      importance_is_predictive_not_gwas = TRUE,
      internal_threads = 1L
    )
  )
}

write_random_forest_development_figure <- function(
    baseline,
    path = file.path("results", "figures", "tassel_tutorial_random_forest_development.pdf")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  summary <- baseline$summary
  grDevices::pdf(path, width = 10, height = 7, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (metric in c("true_genetic_r", "rmse")) {
    matrix_values <- xtabs(
      stats::as.formula(paste(metric, "~ architecture + scheme")),
      data = summary
    )
    graphics::barplot(
      matrix_values,
      beside = TRUE,
      col = c("#4472C4", "#70AD47", "#ED7D31"),
      ylab = metric,
      main = paste("Random Forest", metric, "by architecture and CV scheme"),
      legend.text = rownames(matrix_values),
      args.legend = list(x = "topright", bty = "n")
    )
  }
  graphics::boxplot(
    mtry_ratio ~ scheme,
    data = baseline$selected_hyperparameters,
    col = c("#4472C4", "#70AD47"),
    ylab = "Selected mtry / marker count",
    xlab = "Outer CV scheme",
    main = "Nested-CV Random Forest mtry choices"
  )
  for (task_id in unique(baseline$importance_summary$task_id)) {
    for (scheme in unique(baseline$importance_summary$scheme[
      baseline$importance_summary$task_id == task_id
    ])) {
      current <- baseline$importance_summary[
        baseline$importance_summary$task_id == task_id &
          baseline$importance_summary$scheme == scheme,
        ,
        drop = FALSE
      ]
      current <- current[order(
        -current$mean_permutation_importance, current$marker_id
      ), , drop = FALSE]
      current <- utils::head(current, 20L)
      graphics::par(mar = c(8, 4.5, 3, 1))
      graphics::barplot(
        current$mean_permutation_importance,
        names.arg = current$marker_id,
        las = 2,
        cex.names = 0.65,
        col = "#4472C4",
        ylab = "Mean permutation importance",
        main = paste(task_id, scheme, "top 20 markers")
      )
    }
  }
  invisible(path)
}

write_random_forest_development_outputs <- function(
    baseline,
    processed_path = file.path(
      "data", "processed", "tassel_tutorial", "random_forest_development.rds"
    ),
    table_directory = file.path("results", "tables"),
    figure_directory = file.path("results", "figures")) {
  dir.create(dirname(processed_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    rds = processed_path,
    metrics = file.path(
      table_directory, "tassel_tutorial_random_forest_development_metrics.tsv"
    ),
    summary = file.path(
      table_directory, "tassel_tutorial_random_forest_development_summary.tsv"
    ),
    selected = file.path(
      table_directory, "tassel_tutorial_random_forest_selected_hyperparameters.tsv"
    ),
    tuning = file.path(
      table_directory, "tassel_tutorial_random_forest_tuning_archive.tsv"
    ),
    importance = file.path(
      table_directory, "tassel_tutorial_random_forest_importance_summary.tsv"
    ),
    recovery = file.path(
      table_directory, "tassel_tutorial_random_forest_importance_causal_recovery.tsv"
    )
  )
  saveRDS(baseline, paths[["rds"]], version = 3)
  tables <- list(
    metrics = baseline$metrics,
    summary = baseline$summary,
    selected = baseline$selected_hyperparameters,
    tuning = baseline$tuning_archive,
    importance = baseline$importance_summary,
    recovery = baseline$importance_causal_recovery
  )
  for (name in names(tables)) {
    utils::write.table(
      tables[[name]], paths[[name]], sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  figure <- write_random_forest_development_figure(
    baseline,
    file.path(figure_directory, "tassel_tutorial_random_forest_development.pdf")
  )
  unname(c(paths, figure))
}
