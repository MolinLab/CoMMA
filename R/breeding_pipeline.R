# Breeder-facing orchestration over the existing model and native imputation code.
# Evaluation, resampling and tuning scripts are deliberately not dependencies.
breeding_pipeline_sources <- function() unique(c(gbs_parallel_sources(), "R/genotype_processing.R", "R/gbs_preprocessing.R",
  "R/additive_baseline.R", "R/elastic_net_baseline.R", "R/random_forest_baseline.R", "R/xgboost_baseline.R",
  "R/breeder_prediction_models.R", "R/population_structure.R", "R/gwas.R", "R/breeding_pipeline.R"))

breeding_timed_step <- function(record, step, code, method = NA_character_) {
  if (is.null(record)) return(force(code))
  start <- proc.time()[["elapsed"]]
  status <- "incomplete"
  on.exit(record(step, method, proc.time()[["elapsed"]] - start, status))
  value <- force(code)
  status <- "complete"
  value
}

breeding_metadata <- function(metadata, ids) {
  if (!is.data.frame(metadata) || !"sample_id" %in% names(metadata) || anyNA(metadata$sample_id) || anyDuplicated(metadata$sample_id)) stop("Unique sample metadata IDs are required.")
  i <- match(ids, metadata$sample_id); if (anyNA(i)) stop("Genotype IDs are missing from sample metadata.")
  out <- metadata[i, , drop = FALSE]; rownames(out) <- NULL
  if (!"source_group" %in% names(out)) out$source_group <- "unspecified"
  if (anyNA(out$source_group) || any(!nzchar(out$source_group))) stop("Missing source groups; use an explicit 'unspecified' label.")
  out
}

breeding_check_config <- function(config) {
  if (config$schema_version != 1L || !is.finite(config$sample_missing_max) || config$sample_missing_max <= 0 || config$sample_missing_max > 1 ||
      config$imputation$ld < 2L || config$imputation$k < 1L || config$imputation$ld > config$model$max_ld_sites ||
      config$imputation$maximum_residual_fraction < 0 || config$imputation$maximum_residual_fraction > 1) stop("Invalid breeding pipeline settings.")
  for (name in c("predictor_qc", "analysis_qc")) {
    s <- config[[name]]
    if (any(!is.finite(unlist(s))) || s$maf_min < 0 || s$maf_min >= .5 || s$marker_missing_max <= 0 || s$marker_missing_max > 1 ||
        s$marker_heterozygosity_max <= 0 || s$marker_heterozygosity_max > 1 || s$min_mac < 0) stop("Invalid marker QC settings.")
  }
  if (config$analysis_qc$maf_min < config$predictor_qc$maf_min || config$analysis_qc$marker_missing_max > config$predictor_qc$marker_missing_max ||
      config$analysis_qc$marker_heterozygosity_max > config$predictor_qc$marker_heterozygosity_max || config$analysis_qc$min_mac < config$predictor_qc$min_mac)
    stop("Analysis QC must be nested within predictor QC.")
  gbs_parallel_budget(config$resources, 0, .25, 1L); invisible(TRUE)
}

breeding_complete_calls <- function(imputer, raw, analysis_ids, config, checkpoint) {
  x <- raw[, analysis_ids, drop = FALSE]; missing <- which(is.na(x), arr.ind = TRUE)
  allocation <- NULL
  if (nrow(missing)) {
    targets <- data.frame(sample_id = rownames(x)[missing[, 1L]], marker_id = colnames(x)[missing[, 2L]])
    grid <- data.frame(ld = config$imputation$ld, k = config$imputation$k, id = "fixed")
    pred <- predict_gbs_targets_parallel(imputer, raw, targets, grid, config$resources, checkpoint)
    pos <- cbind(match(pred$targets$sample_id, rownames(x)), match(pred$targets$marker_id, colnames(x)))
    if (anyNA(pos)) stop("Imputation target identities changed.")
    x[pos] <- pred$prediction[, 1L]; allocation <- pred$allocation
  }
  residual <- is.na(x); fraction <- rowMeans(residual)
  if (any(fraction > config$imputation$maximum_residual_fraction)) stop("Residual missingness exceeds the configured limit for: ", paste(head(rownames(x)[fraction > config$imputation$maximum_residual_fraction], 5), collapse = ", "))
  means <- 2 * imputer$marker_statistics$alt_frequency[match(analysis_ids, imputer$marker_ids)]
  for (j in which(colSums(residual) > 0L)) x[residual[, j], j] <- means[j]
  if (any(!is.finite(x))) stop("Analysis dosages remain incomplete.")
  list(genotype = x, audit = data.frame(sample_id = rownames(x), missing_analysis_calls = rowSums(is.na(raw[, analysis_ids, drop = FALSE])),
    residual_mean_fills = rowSums(residual), residual_fraction = fraction), allocation = allocation)
}

breeding_fit_preprocessor <- function(reference, markers, config, checkpoint) {
  breeding_check_config(config); gbs_validate_calls(reference)
  gbs_parallel_budget(config$resources, as.double(length(reference)) * 8 * 12 / 1024^3 + .5, .25, 1L)
  reference <- reference[sort(rownames(reference)), sort(colnames(reference)), drop = FALSE]
  markers <- gbs_preprocessing_map(markers, colnames(reference))
  initial <- gbs_reference_qc(reference, config$predictor_qc); pool <- initial$marker_id[initial$retained]
  if (length(pool) < 3L) stop("Too few initial predictor SNPs.")
  missing <- rowMeans(is.na(reference[, pool, drop = FALSE])); retained <- missing <= config$sample_missing_max
  audit <- data.frame(sample_id = rownames(reference), natural_missingness = missing, retained = retained,
    reason = ifelse(retained, "retained", "missingness_above_ceiling"), pool = "initial_training_predictor_QC")
  if (sum(retained) < 3L) stop("Too few training lines after the recorded line exclusion.")
  reference <- reference[retained, , drop = FALSE]
  predictor <- gbs_reference_qc(reference, config$predictor_qc); analysis <- gbs_reference_qc(reference, config$analysis_qc)
  p <- predictor$marker_id[predictor$retained]; a <- analysis$marker_id[analysis$retained]
  if (length(p) < 3L || !length(a) || any(!a %in% p)) stop("Insufficient nested predictor/analysis pools.")
  raw <- reference[, p, drop = FALSE]
  imputer <- fit_gbs_tassel_parallel(raw, markers, config$model, config$resources, a, file.path(checkpoint, "LD"))
  training <- breeding_complete_calls(imputer, raw, a, config, file.path(checkpoint, "training_calls"))
  center <- colMeans(training$genotype); scale <- apply(training$genotype, 2L, stats::sd)
  if (any(!is.finite(scale) | scale <= 1e-8)) stop("An analysis SNP is constant after preprocessing.")
  training$z <- sweep(sweep(training$genotype, 2L, center, "-"), 2L, scale, "/")
  structure(list(imputer = imputer, predictor_ids = p, analysis_ids = a, training_ids = rownames(raw), center = center, scale = scale,
    training = training, markers = markers[match(a, markers$marker_id), ], initial_QC = initial, predictor_QC = predictor, analysis_QC = analysis,
    line_audit = audit, config = config, provenance = list(phenotype_used = FALSE, fitted_on_training_only = TRUE,
      line_filter_applied_once_before_marker_refit = TRUE, fixed_imputation = config$imputation)), class = "breeding_pipeline_preprocessor")
}

breeding_apply_preprocessor <- function(query, markers, preprocessor, checkpoint, resources = NULL) {
  if (!inherits(preprocessor, "breeding_pipeline_preprocessor")) stop("Wrong saved preprocessing object.")
  gbs_validate_calls(query); gbs_preprocessing_map(markers, preprocessor$predictor_ids, preprocessor$imputer$markers)
  raw <- query[, preprocessor$predictor_ids, drop = FALSE]
  shared <- intersect(rownames(raw), preprocessor$training_ids)
  if (length(shared) && !identical(unname(raw[shared, , drop = FALSE]), unname(preprocessor$imputer$reference[shared, , drop = FALSE]))) stop("A training sample ID was reused with different raw calls.")
  config <- preprocessor$config; if (!is.null(resources)) config$resources <- resources
  breeding_check_config(config)
  missing <- rowMeans(is.na(raw)); keep <- missing <= config$sample_missing_max
  audit <- data.frame(sample_id = rownames(raw), natural_missingness = missing, retained = keep,
    reason = ifelse(keep, "retained", "missingness_above_ceiling"), pool = "frozen_training_predictor_pool")
  if (!any(keep)) return(list(genotype = matrix(numeric(), 0, length(preprocessor$analysis_ids), dimnames = list(character(), preprocessor$analysis_ids)),
    z = matrix(numeric(), 0, length(preprocessor$analysis_ids), dimnames = list(character(), preprocessor$analysis_ids)), line_audit = audit, audit = NULL))
  result <- breeding_complete_calls(preprocessor$imputer, raw[keep, , drop = FALSE], preprocessor$analysis_ids, config, checkpoint)
  result$z <- sweep(sweep(result$genotype, 2L, preprocessor$center, "-"), 2L, preprocessor$scale, "/")
  result$line_audit <- audit; result
}

breeding_save_preprocessor <- function(reference, markers, config, path, checkpoint, shared = NULL,
                                      reference_hash = NULL, marker_hash = NULL) {
  cached <- NULL
  if (!is.null(shared)) {
    key <- digest::digest(list(reference_hash %||% digest::digest(reference, algo = "sha256"),
      marker_hash %||% digest::digest(markers, algo = "sha256"), config), algo = "sha256")
    cached <- shared$preprocessors[[key]]
  }
  if (is.null(cached)) {
    prep <- breeding_fit_preprocessor(reference, markers, config, checkpoint)
  } else {
    source <- gbs_project_path(cached$path)
    if (!identical(digest::digest(source, file = TRUE, algo = "sha256"), cached$sha256)) {
      stop("Shared genotype preprocessing changed: ", cached$path)
    }
    prep <- readRDS(source)
    if (!inherits(prep, "breeding_pipeline_preprocessor") || !identical(prep$config, config) ||
        !identical(prep$line_audit$sample_id, sort(rownames(reference)))) {
      stop("Shared genotype preprocessing does not match this training panel.")
    }
    cat("Reusing genotype preprocessing from ", cached$path, "\n", sep = "")
  }
  gbs_parallel_atomic(prep, path)
  if (!is.null(shared)) {
    if (is.null(cached)) {
      shared$preprocessors[[key]] <- list(path = path, sha256 = digest::digest(path, file = TRUE, algo = "sha256"))
    }
    shared$usage[[path]] <- data.frame(preprocessor = path, source = if (is.null(cached)) path else cached$path,
      reused = !is.null(cached), lines_before_qc = nrow(reference), lines_after_qc = length(prep$training_ids))
  }
  prep
}

breeding_prepare_panel <- function(reference, markers, output, config = "config/breeding_pipeline.yml", shared = NULL) {
  if (is.character(config)) config <- yaml::read_yaml(gbs_project_path(config))
  gbs_project_path(output, must_exist = FALSE); if (dir.exists(output)) stop("Preparation output exists; choose a new directory.")
  dir.create(output, recursive = TRUE)
  prep <- breeding_save_preprocessor(reference, markers, config, file.path(output, "preprocessor.rds"),
                                    file.path(output, "checkpoints"), shared)
  data.table::fwrite(prep$line_audit, file.path(output, "training_line_QC.tsv"), sep = "\t")
  data.table::fwrite(prep$training$audit, file.path(output, "imputation_support.tsv"), sep = "\t")
  prep
}

breeding_learner_status <- function(learners) {
  supported <- c("gblup_reml", "elastic_net", "random_forest", "xgboost", "gaussian_rkhs")
  if (!length(learners) || anyDuplicated(learners) || any(!learners %in% supported))
    stop("This entry currently supports GBLUP, Elastic Net, Random Forest, XGBoost and RKHS. BayesC still requires the production convergence wrapper; BayesB/NN remain optional pending validation. No learner is substituted.")
  packages <- c(elastic_net = "glmnet", random_forest = "ranger", xgboost = "xgboost")
  needed <- unname(packages[intersect(learners, names(packages))])
  if (length(needed) && any(!vapply(needed, requireNamespace, logical(1L), quietly = TRUE))) stop("Restore the project-local model packages before fitting: ", paste(needed, collapse = ", "))
  invisible(TRUE)
}

breeding_fit_learner <- function(learner, z, y, config, seed) {
  estimated <- as.double(length(z)) * 8 * 6 / 1024^3 + as.double(nrow(z))^2 * 8 * 12 / 1024^3 + .25
  gbs_parallel_budget(config$resources, estimated, .25, 1L)
  settings <- config$learners[[learner]]
  fit <- switch(learner,
    gblup_reml = fit_gblup(z, unname(y)),
    gaussian_rkhs = fit_gaussian_rkhs(z, unname(y), settings$gamma, settings$lambda),
    elastic_net = {
      lambda <- elastic_net_lambda_max(z, unname(y), settings$alpha) * settings$lambda_ratio
      list(fit = fit_elastic_net_path(z, unname(y), settings$alpha, lambda), lambda = lambda)
    },
    random_forest = fit_random_forest_model(z, unname(y), as.data.frame(settings[setdiff(names(settings), "num_trees")]), settings$num_trees, seed),
    xgboost = {
      booster <- fit_xgboost_model(z, unname(y), as.data.frame(settings[setdiff(names(settings), "nrounds")]), settings$nrounds,
        list(internal_threads = 1L, tree_method = "hist"), seed)$fit
      list(serialized_booster = xgboost::xgb.save.raw(booster))
    }, stop("Unknown learner."))
  list(learner = learner, fit = fit, marker_ids = colnames(z), training_ids = rownames(z), seed = as.integer(seed))
}

breeding_predict_learner <- function(model, z) {
  if (!identical(colnames(z), model$marker_ids) || any(!is.finite(z))) stop("Saved learner marker schema differs.")
  predicted <- switch(model$learner,
    gblup_reml = predict_gblup(model$fit, z), gaussian_rkhs = predict_gaussian_rkhs(model$fit, z),
    elastic_net = as.numeric(stats::predict(model$fit$fit, newx = z, s = model$fit$lambda)),
    random_forest = predict_random_forest_model(model$fit, z),
    xgboost = predict_xgboost_model(xgboost::xgb.load.raw(model$fit$serialized_booster), z), stop("Unknown saved learner."))
  if (length(predicted) != nrow(z) || any(!is.finite(predicted))) stop("Invalid saved-model predictions.")
  as.numeric(predicted)
}

breeding_model_parameters <- function(model) {
  stopifnot(inherits(model, "breeding_trait_model"))
  configured_names <- list(gblup_reml = character(), elastic_net = c("alpha", "lambda_ratio"),
    gaussian_rkhs = c("gamma", "lambda"), random_forest = c("mtry_ratio", "min_node_size", "sample_fraction", "max_depth", "num_trees"),
    xgboost = c("eta", "max_depth", "min_child_weight", "subsample", "colsample_bytree", "lambda", "alpha", "gamma", "nrounds"))
  do.call(rbind, lapply(model$models, function(learner) {
    id <- learner$learner
    fit <- learner$fit
    settings <- model$request$config$learners[[id]][configured_names[[id]]]
    derived <- switch(id, gblup_reml = list(delta = fit$delta, lambda_marker = fit$lambda_marker),
      elastic_net = list(lambda = fit$lambda), gaussian_rkhs = list(distance_scale = fit$distance_scale),
      random_forest = list(mtry = attr(fit, "realized_mtry")), xgboost = list())
    fixed <- switch(id, random_forest = list(replace = FALSE, splitrule = "variance"),
      xgboost = list(tree_method = "hist", objective = "reg:squarederror"), list())
    values <- c(settings, derived, fixed)
    if (any(lengths(values) != 1L)) stop("Missing fitted or configured settings for ", id, ".")
    source <- c(rep("configured", length(settings)), rep("training_derived", length(derived)), rep("pipeline_fixed", length(fixed)))
    if (identical(id, "gblup_reml")) source[names(values) == "delta"] <- "estimated_REML"
    data.frame(learner = id, parameter = names(values), value = vapply(values, as.character, character(1L)),
      source = source, training_lines = length(learner$training_ids), markers = length(learner$marker_ids), seed = learner$seed)
  }))
}

breeding_fit_models <- function(reference, markers, phenotype, sample_metadata, output,
                               learners = c("gblup_reml", "elastic_net", "random_forest", "xgboost", "gaussian_rkhs"),
                               config = "config/breeding_pipeline.yml", seed = 2026103062L, resume = FALSE,
                               shared = NULL, timing = NULL) {
  breeding_learner_status(learners)
  if (is.character(config)) config <- yaml::read_yaml(gbs_project_path(config))
  breeding_check_config(config); gbs_validate_calls(reference)
  metadata <- breeding_metadata(sample_metadata, rownames(reference))
  y <- gbs_training_response(phenotype, rownames(reference))
  paths <- breeding_pipeline_sources()
  request <- list(reference_sha256 = digest::digest(reference, algo = "sha256"), markers_sha256 = digest::digest(markers, algo = "sha256"),
    training_phenotype_sha256 = digest::digest(y, algo = "sha256"), metadata_sha256 = digest::digest(metadata, algo = "sha256"),
    learners = learners, config = config, seed = as.integer(seed), R_version = as.character(getRversion()),
    package_versions = setNames(vapply(c("yaml", "digest", "data.table", "glmnet", "ranger", "xgboost"), function(p)
      tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_), character(1L)),
      c("yaml", "digest", "data.table", "glmnet", "ranger", "xgboost")), sources = data.frame(path = paths,
      sha256 = vapply(paths, digest::digest, character(1L), file = TRUE, algo = "sha256")))
  gbs_project_path(output, must_exist = FALSE)
  if (dir.exists(output)) {
    if (!resume || !file.exists(file.path(output, "request.rds")) || !identical(request, readRDS(file.path(output, "request.rds")))) stop("Output exists or resume inputs changed.")
  } else { dir.create(output, recursive = TRUE); gbs_parallel_atomic(request, file.path(output, "request.rds")) }
  if (file.exists(file.path(output, "completion.rds"))) stop("Completed fit is protected; use its saved model.")
  lock <- file.path(output, ".active_lock"); if (!dir.create(lock, showWarnings = FALSE)) stop("Another or interrupted fitting controller owns this output.")
  saveRDS(list(pid = Sys.getpid(), started = Sys.time()), file.path(lock, "owner.rds")); on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  prep_path <- file.path(output, "preprocessor.rds")
  prep <- breeding_timed_step(timing, "training_preprocessing", {
    if (file.exists(prep_path)) {
      if (digest::digest(prep_path, file = TRUE, algo = "sha256") != readLines(file.path(output, "preprocessor_sha256.txt"), warn = FALSE)) stop("Saved preprocessing cache changed.")
      readRDS(prep_path)
    } else {
      prepared <- breeding_save_preprocessor(reference, markers, config, prep_path, file.path(output, "checkpoints/preprocessing"),
        shared, request$reference_sha256, request$markers_sha256)
      writeLines(digest::digest(prep_path, file = TRUE, algo = "sha256"), file.path(output, "preprocessor_sha256.txt"))
      prepared
    }
  })
  y <- gbs_training_response(phenotype, prep$training_ids); metadata <- breeding_metadata(sample_metadata, prep$training_ids)
  baselines <- list(global_mean = mean(y), cohort_means = tapply(as.numeric(y), metadata$source_group, mean))
  models <- list()
  for (i in seq_along(learners)) {
    path <- file.path(output, paste0("learner_", learners[i], ".rds")); checksum <- paste0(path, ".sha256")
    if (file.exists(path)) {
      if (!file.exists(checksum) || digest::digest(path, file = TRUE, algo = "sha256") != readLines(checksum, warn = FALSE)) stop("Learner checkpoint changed.")
      model <- readRDS(path)
    } else {
      model <- breeding_timed_step(timing, "fit_learner",
        breeding_fit_learner(learners[i], prep$training$z, y, config, as.integer(seed + i)), learners[i])
      gbs_parallel_atomic(model, path); writeLines(digest::digest(path, file = TRUE, algo = "sha256"), checksum)
    }
    models[[learners[i]]] <- model
  }
  object <- structure(list(preprocessing = prep, models = models, baselines = baselines, metadata = metadata,
    training_ids = prep$training_ids, request = request, provenance = list(training_only = TRUE,
      tuning_performed = FALSE, parameter_source = "caller_supplied_fixed_settings", outer_test_phenotypes_used = FALSE)), class = "breeding_trait_model")
  path <- file.path(output, "breeding_model.rds"); gbs_parallel_atomic(object, path)
  data.table::fwrite(breeding_model_parameters(object), file.path(output, "hyperparameters.tsv"), sep = "\t")
  data.table::fwrite(prep$line_audit, file.path(output, "training_line_QC.tsv"), sep = "\t")
  data.table::fwrite(prep$predictor_QC, file.path(output, "predictor_QC.tsv"), sep = "\t")
  data.table::fwrite(prep$analysis_QC, file.path(output, "analysis_QC.tsv"), sep = "\t")
  data.table::fwrite(prep$training$audit, file.path(output, "training_imputation_support.tsv"), sep = "\t")
  gbs_parallel_atomic(list(status = "breeding_models_fitted", model_sha256 = digest::digest(path, file = TRUE, algo = "sha256"),
    learners = learners, rows = length(prep$training_ids), markers = length(prep$analysis_ids)), file.path(output, "completion.rds"))
  object
}

breeding_predict_models <- function(model, query, markers, sample_metadata, output, resources = NULL, timing = NULL) {
  if (is.character(model)) model <- readRDS(gbs_project_path(model))
  if (!inherits(model, "breeding_trait_model")) stop("A saved breeding trait model is required.")
  # Register package-specific prediction methods in a fresh R session as well.
  breeding_learner_status(names(model$models))
  gbs_project_path(output, must_exist = FALSE); if (dir.exists(output)) stop("Prediction output exists; choose a new directory.")
  dir.create(output, recursive = TRUE)
  metadata <- breeding_metadata(sample_metadata, rownames(query))
  if (is.null(resources)) resources <- model$preprocessing$config$resources
  prepared <- breeding_timed_step(timing, "query_preprocessing",
    breeding_apply_preprocessor(query, markers, model$preprocessing, file.path(output, "checkpoints"), resources))
  rows <- lapply(c("training_mean", "cohort_mean", names(model$models)), function(learner) {
    prediction <- rep(NA_real_, nrow(query)); ids <- rownames(prepared$z)
    if (length(ids)) {
      values <- if (learner == "training_mean") rep(model$baselines$global_mean, length(ids)) else if (learner == "cohort_mean") {
        x <- unname(model$baselines$cohort_means[metadata$source_group[match(ids, metadata$sample_id)]])
        x[is.na(x)] <- model$baselines$global_mean; x
      } else breeding_timed_step(timing, "predict_learner",
        breeding_predict_learner(model$models[[learner]], prepared$z), learner)
      prediction[match(ids, rownames(query))] <- values
    }
    data.frame(sample_id = rownames(query), learner = learner, prediction = prediction,
      eligible = prepared$line_audit$retained, reason = prepared$line_audit$reason,
      cohort_baseline_source = if (learner == "cohort_mean") ifelse(metadata$source_group %in% names(model$baselines$cohort_means),
        "training_cohort_mean", "training_global_mean_for_unseen_cohort") else "not_applicable")
  })
  predictions <- do.call(rbind, rows); rownames(predictions) <- NULL
  data.table::fwrite(predictions, file.path(output, "predictions.tsv"), sep = "\t")
  data.table::fwrite(prepared$line_audit, file.path(output, "query_line_QC.tsv"), sep = "\t")
  if (!is.null(prepared$audit)) data.table::fwrite(prepared$audit, file.path(output, "query_imputation_support.tsv"), sep = "\t")
  receipt <- list(status = "saved_breeding_models_applied", phenotype_values_read = FALSE, preprocessing_refitted = FALSE,
    query_rows = nrow(query), eligible_rows = sum(prepared$line_audit$retained), assessment_ids = rownames(query),
    training_ids = model$preprocessing$line_audit$sample_id)
  gbs_parallel_atomic(receipt, file.path(output, "prediction_receipt.rds"))
  list(predictions = predictions, prepared = prepared, receipt = receipt)
}

breeding_score_predictions <- function(predictions, phenotype, assessment_ids, selection_fraction = .10, direction = "higher") {
  if (!is.data.frame(predictions) || !all(c("sample_id", "learner", "prediction", "eligible") %in% names(predictions)) ||
      anyDuplicated(predictions[c("sample_id", "learner")]) || anyNA(predictions[c("sample_id", "learner", "eligible")]) ||
      anyDuplicated(assessment_ids) || !length(assessment_ids) || !direction %in% c("higher", "lower") ||
      !is.finite(selection_fraction) || selection_fraction <= 0 || selection_fraction > 1) stop("Invalid explicit assessment or selection policy.")
  observed <- gbs_training_response(phenotype, assessment_ids)
  eligibility <- NULL; rows <- list()
  for (learner in unique(predictions$learner)) {
    x <- predictions[predictions$learner == learner, , drop = FALSE]
    if (!setequal(x$sample_id, assessment_ids)) stop("Every learner must account for every registered assessment ID.")
    x <- x[match(assessment_ids, x$sample_id), , drop = FALSE]
    if (is.null(eligibility)) eligibility <- x$eligible else if (!identical(eligibility, x$eligible)) stop("Learners cannot be compared on different eligible assessment sets.")
    if (any(!is.finite(x$prediction[x$eligible])) || any(!is.na(x$prediction[!x$eligible]))) stop("Prediction eligibility and values disagree.")
    if (any(x$eligible)) {
      y <- observed[x$eligible]; predicted <- x$prediction[x$eligible]
      metric <- regression_metrics(as.numeric(y), predicted)
      n_top <- max(1L, as.integer(ceiling(length(y) * selection_fraction)))
      sign <- if (direction == "higher") -1 else 1
      true_top <- names(y)[head(order(sign * y, names(y)), n_top)]
      predicted_top <- names(y)[head(order(sign * predicted, names(y)), n_top)]
      overlap <- length(intersect(true_top, predicted_top)) / n_top
    } else {
      metric <- data.frame(rmse = NA_real_, mae = NA_real_, r2 = NA_real_, pearson_r = NA_real_, spearman_r = NA_real_, true_genetic_r = NA_real_)
      n_top <- 0L; overlap <- NA_real_
    }
    rows[[learner]] <- data.frame(learner = learner, n_registered = length(assessment_ids), n_scored = sum(x$eligible),
      n_QC_excluded = sum(!x$eligible), status = if (any(x$eligible)) "scored" else "no_predictions_after_QC", metric,
      selection_fraction = selection_fraction, selection_direction = direction, selection_n = n_top, top_fraction_overlap = overlap)
  }
  result <- do.call(rbind, rows); rownames(result) <- NULL; result
}

breeding_gwas_by_marker <- function(scans) {
  keys <- c("trait", "marker_id", "chromosome", "position_bp")
  metrics <- c("p_value", "bh_fdr", "significant_bonferroni", "significant_fdr", "effect", "standard_error", "bonferroni_threshold")
  if (anyNA(scans[c(keys, "method")]) || anyDuplicated(scans[c("trait", "marker_id", "method")]) ||
      anyDuplicated(unique(scans[keys])[c("trait", "marker_id")])) stop("GWAS results disagree on marker identity or coordinates.")
  wide <- data.table::dcast(data.table::as.data.table(scans), trait + marker_id + chromosome + position_bp ~ method,
                          value.var = metrics)
  columns <- unlist(lapply(unique(scans$method), function(method) paste(metrics, method, sep = "_")), use.names = FALSE)
  data.table::setcolorder(wide, c(keys, columns))
  wide
}

breeding_run_gwas <- function(model, phenotype, sample_metadata, trait, output, n_pcs = 2L, include_source_group = TRUE,
                             models = c("MLM", "FarmCPU", "BLINK"), report_context = FALSE) {
  if (is.character(model)) model <- readRDS(gbs_project_path(model))
  prep <- if (inherits(model, "breeding_pipeline_preprocessor")) model else if (inherits(model, "breeding_trait_model")) model$preprocessing else
    stop("A fitted genotype-preparation object or breeding trait model is required.")
  ids <- prep$training_ids; metadata <- breeding_metadata(sample_metadata, ids)
  if (!"breeding_material_type" %in% names(metadata)) metadata$breeding_material_type <- "unknown"
  if (!"breeding_material_type_source" %in% names(metadata)) metadata$breeding_material_type_source <- "not_supplied_not_inferred"
  map <- prep$markers; map$ref <- substr(map$alleles, 1L, 1L); map$alt <- substr(map$alleles, 3L, 3L)
  response <- gbs_training_response(phenotype, ids)
  dataset <- new_genomic_dataset(prep$training$genotype, metadata, map,
    data.frame(sample_id = ids, trait = trait, value = unname(response)), provenance = list(fractional_dosage_allowed = TRUE))
  structure <- compute_population_structure(dataset, n_pcs = max(3L, as.integer(n_pcs)))
  selection <- list(selected_by_trait = data.frame(trait = trait, selected_n_pcs = as.integer(n_pcs)))
  inputs <- prepare_gapit_numeric_inputs(dataset, structure, selection, trait)
  if (include_source_group && length(unique(metadata$source_group)) > 1L) {
    covariates <- stats::model.matrix(~ factor(metadata$source_group))[, -1L, drop = FALSE]
    colnames(covariates) <- paste0("source_group_", seq_len(ncol(covariates)))
    inputs$CV <- if (is.null(inputs$CV)) data.frame(Taxa = ids, covariates) else cbind(inputs$CV, covariates)
  }
  chromosomes <- as.character(inputs$GM$Chromosome)
  numeric_chr <- suppressWarnings(as.numeric(chromosomes))
  if (anyNA(numeric_chr)) {
    wheat <- paste0(rep(1:7, each = 3), rep(c("A", "B", "D"), 7))
    if (any(!chromosomes %in% wheat)) stop("Provide numeric chromosomes or explicit wheat 1A-7D labels for GAPIT.")
    inputs$GM$Chromosome <- match(chromosomes, wheat)
  } else inputs$GM$Chromosome <- numeric_chr
  validate_gapit_numeric_inputs(inputs)
  gbs_project_path(output, must_exist = FALSE); if (dir.exists(output)) stop("GWAS output exists; preserve it.")
  dir.create(output, recursive = TRUE)
  result <- run_gapit_models(inputs, models = models, output_directory = output)
  result$results$engine_chromosome <- result$results$chromosome
  result$results$chromosome <- map$chromosome[match(result$results$marker_id, map$marker_id)]
  result$provenance$fixed_source_group_covariates <- include_source_group
  result$provenance$descriptive_association_not_GP_feature_selection <- TRUE
  saveRDS(result, file.path(output, "GWAS_results.rds"))
  data.table::fwrite(result$results, file.path(output, "GWAS_results.tsv"), sep = "\t")
  data.table::fwrite(breeding_gwas_by_marker(result$results), file.path(output, "GWAS_by_marker.tsv"), sep = "\t", na = "NA")
  data.table::fwrite(result$summary, file.path(output, "GWAS_summary.tsv"), sep = "\t")
  if (report_context) result$report_context <- list(
    scores = structure$scores, variance = structure$variance, covariates = inputs$CV)
  result
}
