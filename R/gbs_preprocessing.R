# SPDX-License-Identifier: LGPL-2.1-only
# Training-reference preprocessing around the native TASSEL-rule implementation.
# Imputation predictor and analysis marker pools are distinct. No phenotype input.

gbs_validate_preprocessing_config <- function(config) {
  needed <- c("maf_min", "marker_missing_max", "marker_heterozygosity_max", "min_mac")
  for (name in c("predictor_qc", "analysis_qc")) {
    x <- config[[name]]
    if (!is.list(x) || !all(needed %in% names(x)) || any(!is.finite(unlist(x[needed]))) ||
        x$maf_min < 0 || x$maf_min >= .5 || x$marker_missing_max <= 0 || x$marker_missing_max > 1 ||
        x$marker_heterozygosity_max <= 0 || x$marker_heterozygosity_max > 1 || x$min_mac < 0) stop("Invalid predictor/analysis QC settings.")
  }
  p <- config$predictor_qc; a <- config$analysis_qc
  if (a$maf_min < p$maf_min || a$marker_missing_max > p$marker_missing_max ||
      a$marker_heterozygosity_max > p$marker_heterozygosity_max || a$min_mac < p$min_mac) stop("Analysis QC must be at least as strict as predictor QC.")
  if (!identical(config$analysis_marker_policy, "observed_training_reference_only") ||
      !identical(config$reference_sample_failure_action, "error_without_dropping_lines") ||
      !config$selection$mode %in% c("fixed", "automatic") ||
      !config$residual$policy %in% c("reference_mean", "reference_mode", "error", "keep_na") ||
      !is.finite(config$residual$maximum_fraction_per_sample) || config$residual$maximum_fraction_per_sample < 0 ||
      config$residual$maximum_fraction_per_sample > 1 ||
      !is.finite(config$reference_sample_missing_max_after_predictor_QC) ||
      config$reference_sample_missing_max_after_predictor_QC <= 0 || config$reference_sample_missing_max_after_predictor_QC > 1 ||
      !is.finite(config$scaling$minimum_training_sd) || config$scaling$minimum_training_sd <= 0 ||
      config$resources$global_cpu_budget != 1L || !is.finite(config$resources$maximum_memory_gib) ||
      config$resources$maximum_memory_gib <= 0) stop("Invalid preprocessing contract or resource settings.")
  invisible(TRUE)
}

gbs_preprocessing_map <- function(markers, ids, expected = NULL) {
  if (!is.data.frame(markers) || !all(c("marker_id", "chromosome", "position_bp") %in% names(markers)) ||
      anyNA(markers$marker_id) || anyDuplicated(markers$marker_id)) stop("Invalid marker identity metadata.")
  if (!"alleles" %in% names(markers) && all(c("ref", "alt") %in% names(markers))) markers$alleles <- paste(markers$ref, markers$alt, sep = "/")
  if (!"alleles" %in% names(markers)) stop("Numeric calls require explicit allele orientation.")
  index <- match(ids, markers$marker_id)
  if (anyNA(index)) stop("Genotype marker IDs lack metadata.")
  map <- markers[index, , drop = FALSE]; rownames(map) <- NULL
  gbs_tassel_byte_codes(map$alleles)
  if (!is.null(expected)) {
    for (field in c("marker_id", "alleles", "chromosome", "position_bp")) {
      if (anyNA(map[[field]]) || anyNA(expected[[field]]) || any(map[[field]] != expected[[field]])) stop("Query marker metadata/orientation differ: ", field)
    }
    if ("strand" %in% names(map) && "strand" %in% names(expected)) {
      known <- map$strand %in% c("+", "-") & expected$strand %in% c("+", "-")
      if (any(map$strand[known] != expected$strand[known])) stop("Query marker strands differ.")
    }
  }
  map
}

gbs_impute_analysis_targets <- function(imputer, query, analysis_ids, high_ld_sites, knn_taxa,
                                         residual, raw_reference_means, raw_reference_modes) {
  gbs_validate_calls(query)
  if (anyNA(match(imputer$marker_ids, colnames(query))) || any(!analysis_ids %in% imputer$marker_ids)) stop("Missing fitted predictor/analysis marker IDs.")
  pool <- query[, imputer$marker_ids, drop = FALSE]
  shared <- intersect(rownames(imputer$reference), rownames(pool))
  if (length(shared)) {
    old <- unname(imputer$reference[shared, , drop = FALSE]); new <- unname(pool[shared, , drop = FALSE])
    if (!identical(is.na(old), is.na(new)) || any(old[!is.na(old)] != new[!is.na(new)])) stop("A training ID was reused with different raw calls.")
  }
  before <- pool[, analysis_ids, drop = FALSE]
  missing <- which(is.na(before), arr.ind = TRUE)
  after <- before
  counts <- c(ld_knn = 0L, monomorphic = 0L, abstained_before_residual = 0L)
  if (nrow(missing)) {
    targets <- data.frame(sample_id = rownames(before)[missing[, 1]], marker_id = analysis_ids[missing[, 2]])
    cache <- build_gbs_tassel_rank_cache(imputer, pool, high_ld_sites, targets, knn_taxa)
    cache$provenance$scope <- "complete_analysis_target_imputation_auxiliary_predictor_calls_not_reconstructed"
    predicted <- predict_gbs_tassel_rank_cache(cache, knn_taxa)
    after <- predicted[, analysis_ids, drop = FALSE]
    dispositions <- table(factor(cache$targets$disposition, levels = c("ld_knn", "monomorphic", "abstained")))
    counts <- c(ld_knn = unname(dispositions[[1]]), monomorphic = unname(dispositions[[2]]), abstained_before_residual = unname(dispositions[[3]]))
  }
  remaining <- is.na(after)
  fraction <- rowMeans(remaining)
  if (any(fraction > residual$maximum_fraction_per_sample)) {
    bad <- rownames(after)[fraction > residual$maximum_fraction_per_sample]
    stop("Residual missing-call fraction exceeds the prespecified limit for: ", paste(head(bad, 5), collapse = ", "))
  }
  if (residual$policy == "error" && any(remaining)) stop("LD imputation left unsupported analysis calls; residual policy is error.")
  replaced <- remaining & residual$policy %in% c("reference_mean", "reference_mode")
  values <- if (residual$policy == "reference_mean") raw_reference_means else raw_reference_modes
  if (any(replaced)) for (j in seq_along(analysis_ids)) after[replaced[, j], j] <- values[[analysis_ids[[j]]]]
  gbs_validate_imputation_output(after, before)
  list(genotype = after, counts = counts,
    sample_audit = data.frame(sample_id = rownames(after), originally_missing_analysis = rowSums(is.na(before)),
      unsupported_after_LD = rowSums(remaining), residual_fraction = fraction,
      residual_replacements = rowSums(replaced), final_missing = rowSums(is.na(after))),
    marker_audit = data.frame(marker_id = analysis_ids, originally_missing_analysis = colSums(is.na(before)),
      unsupported_after_LD = colSums(remaining), residual_replacements = colSums(replaced), final_missing = colSums(is.na(after))),
    provenance = list(auxiliary_markers_imputed_for_output = FALSE, residual_policy = residual$policy,
      raw_observed_reference_fallback_values = TRUE, observed_analysis_calls_unchanged = TRUE,
      complete = !anyNA(after)))
}

fit_gbs_preprocessor <- function(reference, markers, config, calibration_output = NULL) {
  gbs_validate_preprocessing_config(config)
  estimated <- as.double(length(reference)) * 8 * 12 / 1024^3 + .15
  if (estimated > config$resources$maximum_memory_gib) stop("Preprocessing fit exceeds memory guard before copying input.")
  gbs_validate_calls(reference)
  if (nrow(reference) < 3L) stop("At least three training lines are required.")
  reference <- reference[sort(rownames(reference)), sort(colnames(reference)), drop = FALSE]
  markers <- gbs_preprocessing_map(markers, colnames(reference))
  calibration <- NULL; calibration_split <- NULL
  if (config$selection$mode == "automatic") {
    automatic <- config$automatic
    if (!is.finite(automatic$query_fraction) || automatic$query_fraction <= 0 || automatic$query_fraction >= .5 ||
        !identical(automatic$no_eligible_candidate_action, "error")) stop("Invalid internal calibration design.")
    hash <- digest::digest(paste(rownames(reference), collapse = "|"), algo = "xxhash32", serialize = FALSE)
    seed <- as.integer((as.double(automatic$seed) + strtoi(substr(hash, 1, 7), 16)) %% .Machine$integer.max)
    n_query <- max(automatic$minimum_query_lines, floor(nrow(reference) * automatic$query_fraction))
    if (nrow(reference) - n_query < automatic$minimum_reference_lines) stop("Training partition is too small for the internal calibration split.")
    query_ids <- with_seed_preserved(seed, sample(rownames(reference), n_query))
    reference_ids <- setdiff(rownames(reference), query_ids)
    internal_reference <- reference[reference_ids, , drop = FALSE]
    internal_qc <- gbs_reference_qc(internal_reference, config$predictor_qc)
    internal_ids <- internal_qc$marker_id[internal_qc$retained]
    if (length(internal_ids) < 3L) stop("Internal training reference leaves too few predictor markers.")
    internal_analysis_qc <- gbs_reference_qc(internal_reference, config$analysis_qc)
    internal_analysis_ids <- internal_analysis_qc$marker_id[internal_analysis_qc$retained]
    if (!length(internal_analysis_ids) || any(!internal_analysis_ids %in% internal_ids)) stop("Internal reference leaves no valid analysis targets.")
    search <- yaml::read_yaml(gbs_project_path(automatic$search_config_file))
    search$resources$maximum_memory_gib <- config$resources$maximum_memory_gib
    search$resources$global_cpu_budget <- 1L
    search$masking$seed <- as.integer((as.double(seed) + 1) %% .Machine$integer.max)
    settings <- config$model; settings$max_ld_sites <- search$search$maximum_LD_sites
    settings$maximum_memory_gib <- config$resources$maximum_memory_gib
    internal_model <- fit_gbs_tassel_native(internal_reference[, internal_ids, drop = FALSE], markers, settings)
    auxiliary <- setdiff(internal_ids, internal_analysis_ids)
    exclusion_estimate <- as.double(length(query_ids)) * length(auxiliary) * 160 / 1024^3 +
      internal_model$provenance$estimated_fit_memory_gib + as.double(length(reference)) * 8 * 6 / 1024^3 + .10
    if (exclusion_estimate > config$resources$maximum_memory_gib) stop("Analysis-target calibration exclusions exceed the memory guard.")
    excluded <- if (length(auxiliary)) expand.grid(sample_id = query_ids, marker_id = auxiliary, stringsAsFactors = FALSE) else NULL
    calibration <- tune_gbs_tassel_neighbours(internal_model, reference[query_ids, internal_ids, drop = FALSE],
      search, excluded_cells = excluded, output_dir = calibration_output)
    if (is.null(calibration$selection$selected)) stop("No internal calibration candidate met coverage; no fixed-setting fallback was applied.")
    high_ld <- calibration$selection$selected$ld; k <- calibration$selection$selected$k
    calibration_split <- list(reference_ids = reference_ids, query_ids = query_ids, seed = seed,
      predictor_ids = internal_ids, analysis_ids = internal_analysis_ids,
      tuning_masks_analysis_targets_only = TRUE, auxiliary_calls_remain_predictors = TRUE,
      QC_fitted_on_internal_reference_only = TRUE)
  } else {
    high_ld <- config$selection$fixed_high_ld_sites; k <- config$selection$fixed_knn_taxa
  }
  if (length(high_ld) != 1L || !is.finite(high_ld) || high_ld != floor(high_ld) || high_ld < 2 || high_ld > 200 ||
      length(k) != 1L || !is.finite(k) || k != floor(k) || k < 2 || k > 200) stop("Invalid selected LD/k settings.")
  predictor_qc <- gbs_reference_qc(reference, config$predictor_qc)
  analysis_qc <- gbs_reference_qc(reference, config$analysis_qc)
  pool_ids <- predictor_qc$marker_id[predictor_qc$retained]
  analysis_ids <- analysis_qc$marker_id[analysis_qc$retained]
  if (length(pool_ids) < 3L || !length(analysis_ids) || any(!analysis_ids %in% pool_ids)) stop("Invalid/empty fitted predictor and analysis pools.")
  sample_missing <- rowMeans(is.na(reference[, pool_ids, drop = FALSE]))
  if (any(sample_missing > config$reference_sample_missing_max_after_predictor_QC)) stop("A training line fails reference sample QC; no lines were silently dropped.")
  settings <- config$model; settings$max_ld_sites <- high_ld; settings$maximum_memory_gib <- config$resources$maximum_memory_gib
  imputer <- fit_gbs_tassel_native(reference[, pool_ids, drop = FALSE], markers, settings)
  means <- colMeans(reference[, analysis_ids, drop = FALSE], na.rm = TRUE)
  modes <- vapply(analysis_ids, function(id) {
    counts <- tabulate(reference[, id] + 1L, nbins = 3L)
    winners <- which(counts == max(counts)); j <- match(id, imputer$marker_ids)
    winners[which.min(imputer$byte_codes[j, winners])] - 1L
  }, numeric(1L))
  training <- gbs_impute_analysis_targets(imputer, reference[, pool_ids, drop = FALSE], analysis_ids,
    high_ld, k, config$residual, means, modes)
  center <- colMeans(training$genotype, na.rm = TRUE)
  scale <- apply(training$genotype, 2L, stats::sd, na.rm = TRUE)
  variable <- is.finite(center) & is.finite(scale) & scale > config$scaling$minimum_training_sd
  if (!all(variable)) stop("An analysis marker is constant/undefined after training preprocessing; revise training QC explicitly.")
  training$z <- sweep(sweep(training$genotype, 2L, center, "-"), 2L, scale, "/")
  ledger <- predictor_qc
  ledger$passes_analysis_QC <- analysis_qc$retained[match(ledger$marker_id, analysis_qc$marker_id)]
  ledger$role <- ifelse(ledger$passes_analysis_QC, "analysis_and_imputation", ifelse(ledger$retained, "imputation_predictor_only", "excluded"))
  structure(list(imputer = imputer, predictor_ids = pool_ids, analysis_ids = analysis_ids,
    markers = markers[match(analysis_ids, markers$marker_id), , drop = FALSE],
    raw_reference_means = means, raw_reference_modes = modes, center = center, scale = scale,
    selected = list(high_ld_sites = high_ld, knn_taxa = k), training = training,
    training_ids = rownames(reference), marker_ledger = ledger, config = config,
    calibration = calibration, calibration_split = calibration_split,
    provenance = list(method = paste0("native_TASSEL_LD_kNN_residual_", config$residual$policy),
      fitted_on_training_reference_only = TRUE, analysis_QC_on_observed_training_calls = TRUE,
      imputation_pool_MAF = config$predictor_qc$maf_min, analysis_MAF = config$analysis_qc$maf_min,
      phenotype_used = FALSE, query_used_for_fit = FALSE,
      hyperparameter_source = if (config$selection$mode == "automatic") "internal_training_reference_mask_calibration" else "predeclared_fixed",
      reference_sha256 = digest::digest(reference, algo = "sha256"), config_sha256 = digest::digest(config, algo = "sha256"),
      memory_guard_gib = config$resources$maximum_memory_gib, CPU_budget = 1L,
      fractional_dosage_allowed = config$residual$policy == "reference_mean")), class = "gbs_preprocessor")
}

apply_gbs_preprocessor <- function(query, markers, preprocessor, require_complete = TRUE) {
  if (!inherits(preprocessor, "gbs_preprocessor")) stop("Invalid fitted GBS preprocessor.")
  estimate <- preprocessor$imputer$provenance$estimated_fit_memory_gib + as.double(length(query)) * 8 * 6 / 1024^3
  if (estimate > preprocessor$config$resources$maximum_memory_gib) stop("Reference plus query exceeds the preprocessing memory guard.")
  gbs_validate_calls(query)
  gbs_preprocessing_map(markers, preprocessor$predictor_ids, preprocessor$imputer$markers)
  result <- gbs_impute_analysis_targets(preprocessor$imputer, query, preprocessor$analysis_ids,
    preprocessor$selected$high_ld_sites, preprocessor$selected$knn_taxa, preprocessor$config$residual,
    preprocessor$raw_reference_means, preprocessor$raw_reference_modes)
  result$z <- sweep(sweep(result$genotype, 2L, preprocessor$center, "-"), 2L, preprocessor$scale, "/")
  if (require_complete && any(!is.finite(result$z))) stop("Learner input is incomplete under the configured residual policy.")
  result
}

prepare_gbs_fold <- function(genotype, markers, split, config, calibration_output = NULL) {
  for (name in c("train", "test")) {
    ids <- split[[name]]
    if (!is.character(ids) || !length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
        anyNA(match(ids, rownames(genotype)))) stop("Invalid explicit training/test split IDs.")
  }
  if (length(intersect(split$train, split$test))) stop("Training/test split IDs overlap.")
  fitted <- fit_gbs_preprocessor(genotype[split$train, , drop = FALSE], markers, config, calibration_output)
  train <- fitted$training
  for (name in c("genotype", "z")) train[[name]] <- train[[name]][split$train, , drop = FALSE]
  train$sample_audit <- train$sample_audit[match(split$train, train$sample_audit$sample_id), , drop = FALSE]
  if (any(!is.finite(train$z))) stop("Training learner input is incomplete.")
  test <- apply_gbs_preprocessor(genotype[split$test, , drop = FALSE], markers, fitted)
  list(id = split$id %||% "explicit_split", split = split, preprocessing = fitted, train = train, test = test,
    provenance = list(train_test_disjoint = TRUE, phenotype_used = FALSE, same_marker_order = identical(colnames(train$z), colnames(test$z))))
}

gbs_training_response <- function(phenotype, training_ids) {
  if (!is.data.frame(phenotype) || !all(c("sample_id", "value") %in% names(phenotype)) ||
      anyNA(phenotype$sample_id) || anyDuplicated(phenotype$sample_id)) stop("Invalid phenotype IDs.")
  index <- match(training_ids, phenotype$sample_id)
  if (anyNA(index) || any(!is.finite(phenotype$value[index]))) stop("Training response is missing/nonfinite.")
  setNames(phenotype$value[index], training_ids)
}
