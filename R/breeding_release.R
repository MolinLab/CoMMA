# Final fits and explanations are explicitly selected and run separately.
breeding_release_choices <- function(selected, available, label, allow_empty = TRUE) {
  if (!is.character(selected) || anyNA(selected) || anyDuplicated(selected) ||
      any(!selected %in% available) || (!allow_empty && !length(selected))) {
    stop("Choose ", label, " from: ", paste(available, collapse = ", "), call. = FALSE)
  }
  selected
}

breeding_assessment_split <- function(metadata, fraction, seed) {
  if (!is.numeric(fraction) || length(fraction) != 1L || !is.finite(fraction) || fraction <= 0 || fraction >= .5) {
    stop("assessment_fraction must be greater than zero and less than 0.5.")
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
    stop("The assessment seed must be a non-negative integer.")
  }
  x <- metadata[order(metadata$sample_id, method = "radix"), , drop = FALSE]
  assert_unique_nonmissing(x$sample_id, "Assessment panel IDs")
  genetic <- x$genetic_stratum
  groups <- if (is.null(x$fingerprint_group)) x$sample_id else as.character(x$fingerprint_group)
  eligible <- if (is.null(x$assessment_eligible)) rep(TRUE, nrow(x)) else as.logical(x$assessment_eligible)
  if (anyNA(eligible) || anyNA(groups) || any(!nzchar(groups)) || (!is.null(genetic) && (anyNA(genetic) || any(!nzchar(genetic))))) {
    stop("Assessment metadata must give valid strata, fingerprint groups and eligibility for every line.")
  }
  stratum <- if (is.null(genetic)) x$source_group else paste(x$source_group, genetic, sep = " / ")
  members <- split(seq_len(nrow(x)), groups)
  group_strata <- vapply(members, function(i) stratum[i[1L]], character(1L))
  can_assess <- vapply(members, function(i) all(eligible[i]), logical(1L))
  levels <- sort(unique(group_strata), method = "radix")
  group_count <- table(factor(group_strata, levels = levels))
  eligible_count <- table(factor(group_strata[can_assess], levels = levels))
  capacity <- pmin(eligible_count, group_count - 1L)
  target_lines <- max(1L, round(fraction * nrow(x)))
  eligible_lines <- sum(lengths(members)[can_assess])
  if (target_lines > eligible_lines || nrow(x) - target_lines < 3L) stop("Too few eligible training or assessment lines for the requested split.")
  target_groups <- max(1L, round(target_lines * sum(can_assess) / eligible_lines))
  if (target_groups > sum(capacity)) stop("The requested holdout cannot keep a training group in every stratum. Use fewer strata or a smaller holdout.")
  expected <- target_groups * as.numeric(eligible_count) / sum(eligible_count)
  take <- pmin(floor(expected), as.integer(capacity))
  while (sum(take) < target_groups) {
    deficit <- ifelse(take < capacity, expected - take, -Inf)
    take[which.max(deficit)] <- take[which.max(deficit)] + 1L
  }
  selected <- with_seed_preserved(seed, unlist(lapply(seq_along(levels), function(i) {
    available <- names(members)[can_assess & group_strata == levels[i]]
    available[sample.int(length(available), size = take[i])]
  }), use.names = FALSE))
  is_test <- groups %in% selected
  if (!any(is_test) || sum(!is_test) < 3L || any(is_test & !eligible)) stop("Invalid generated holdout partition.")
  data.frame(sample_id = x$sample_id, source_group = x$source_group, stratum = stratum,
             allocation_stratum = unname(group_strata[groups]), fingerprint_group = groups,
             assessment_eligible = eligible, role = ifelse(is_test, "assessment", "training"))
}

breeding_run_assessment <- function(reference, markers, observed, metadata, registry, learners, config, seed, output, trait, unit, plots,
                                    shared = NULL, timing = NULL) {
  train_ids <- registry$sample_id[registry$role == "training"]
  query_ids <- registry$sample_id[registry$role == "assessment"]
  dir.create(output)
  registry_file <- file.path(output, "assessment_split.tsv")
  data.table::fwrite(registry, registry_file, sep = "\t")
  gp <- file.path(output, "GP")
  fitted <- breeding_fit_models(reference[train_ids, , drop = FALSE], markers,
    observed[match(train_ids, observed$sample_id), , drop = FALSE], metadata, gp, learners, config, seed, shared = shared, timing = timing)
  prediction_dir <- file.path(output, "predictions")
  prediction <- breeding_predict_models(fitted, reference[query_ids, , drop = FALSE], markers, metadata,
                                        prediction_dir, config$resources, timing = timing)
  response <- observed[match(query_ids, observed$sample_id), , drop = FALSE]
  metrics <- breeding_timed_step(timing, "score_predictions", breeding_score_predictions(prediction$predictions, response, query_ids))
  metrics_file <- file.path(output, "assessment_metrics.tsv")
  data.table::fwrite(metrics, metrics_file, sep = "\t")
  receipt <- list(assessment_ids = query_ids, training_ids = train_ids, split_file = registry_file, model_dir = gp,
    predictions_dir = prediction_dir, metrics_file = metrics_file, seed = seed, phenotype_values_used_to_assign_split = FALSE)
  if (isTRUE(plots)) receipt$plot_report <- breeding_timed_step(timing, "prediction_plots",
    breeding_plot_report(file.path(output, "plots"), predictions = prediction, observed = response, trait = trait, unit = unit))
  gbs_parallel_atomic(receipt, file.path(output, "assessment_completion.rds"))
  receipt
}

breeding_run_release_gwas <- function(model, phenotype, sample_metadata, trait, output, n_pcs = 2L,
                                      include_source_group = TRUE, models = c("MLM", "FarmCPU", "BLINK"),
                                      gwas_candidates = NULL, plots = TRUE, unit = "trait units", timing = NULL) {
  if (isTRUE(plots)) breeding_plot_requirements()
  breeding_check_gwas_candidates(gwas_candidates, models)
  if (is.character(model)) model <- readRDS(gbs_project_path(model))
  prep <- if (inherits(model, "breeding_pipeline_preprocessor")) model else if (inherits(model, "breeding_trait_model"))
    model$preprocessing else stop("A fitted genotype-preparation object or breeding trait model is required.")
  result <- breeding_run_gwas(prep, phenotype, sample_metadata, trait, output, n_pcs, include_source_group, models, report_context = TRUE)
  if (!is.null(timing)) for (i in seq_len(nrow(result$summary))) {
    timing("gwas_model", result$summary$method[i], result$summary$elapsed_seconds[i], "complete")
  }
  if (!is.null(gwas_candidates)) {
    raw <- prep$imputer$reference
    if (!is.matrix(raw) || !identical(rownames(raw), prep$training_ids) || !all(prep$analysis_ids %in% colnames(raw))) {
      stop("The shortlist needs the original genotype calls stored by the fitted preprocessor.")
    }
    ranking <- breeding_gwas_rank_regions(result$results, prep$markers)
    result$candidate_shortlist <- breeding_gwas_clump_regions(ranking, raw, gwas_candidates)
    result$candidate_shortlist_file <- file.path(output, "GWAS_candidate_shortlist.tsv")
    temp <- tempfile("shortlist_", tmpdir = output)
    on.exit(unlink(temp), add = TRUE)
    data.table::fwrite(result$candidate_shortlist, temp, sep = "\t", na = "NA")
    if (!file.rename(temp, result$candidate_shortlist_file)) stop("Could not save the candidate shortlist.")
  }
  context <- result$report_context
  ids <- prep$training_ids
  evidence <- breeding_gwas_plot_evidence(prep, gbs_training_response(phenotype, ids), breeding_metadata(sample_metadata, ids),
    context, context$covariates, result$results, result$candidate_shortlist)
  gbs_parallel_atomic(evidence, file.path(output, "report_data.rds"))
  result$report_context <- NULL
  if (isTRUE(plots)) result$plot_report <- breeding_timed_step(timing, "gwas_plots",
    breeding_plot_report(file.path(output, "plots"), gwas_output = output, unit = unit))
  result
}

breeding_check_trait_methods <- function(gp, gwas, gwas_candidates, assessment_fraction, plots) {
  breeding_release_choices(gp, c("gblup_reml", "elastic_net", "random_forest",
    "xgboost", "gaussian_rkhs"), "GP learners")
  breeding_release_choices(gwas, c("MLM", "FarmCPU", "BLINK"), "GWAS methods")
  if (!length(c(gp, gwas))) stop("Select at least one GP learner or GWAS method.")
  breeding_check_gwas_candidates(gwas_candidates, gwas)
  if (length(gp)) breeding_learner_status(gp)
  if (length(gwas) && !isTRUE(gapit_backend_status()$compatible)) stop("Restore the pinned GAPIT library for GWAS.")
  if (!is.null(assessment_fraction) && !length(gp)) stop("Automatic assessment requires at least one GP learner.")
  if ((length(gwas) || !is.null(assessment_fraction)) && isTRUE(plots)) breeding_plot_requirements()
}

breeding_trait_panel <- function(sample_ids, phenotypes, trait, sample_metadata, minimum_lines, assessment_fraction, seed,
                                 sort_ids = FALSE) {
  rows <- phenotypes[which(phenotypes$trait == trait), c("sample_id", "value"), drop = FALSE]
  assert_unique_nonmissing(as.character(rows$sample_id), paste("Sample IDs for", trait))
  if (any(is.infinite(rows$value))) stop("Selected trait contains infinite phenotype values: ", trait)
  observed <- rows[!is.na(rows$value), , drop = FALSE]
  ids <- intersect(sample_ids, as.character(observed$sample_id))
  if (length(ids) < minimum_lines) stop(trait, " requires at least ", minimum_lines, " matched lines with finite phenotypes.")
  if (sort_ids) ids <- sort(ids, method = "radix")
  observed <- observed[match(ids, observed$sample_id), , drop = FALSE]
  observed$sample_id <- ids
  metadata <- breeding_metadata(sample_metadata %||% data.frame(sample_id = ids), ids)
  registry <- if (!is.null(assessment_fraction)) breeding_assessment_split(metadata, assessment_fraction, seed)
  list(observed = observed, ids = ids, metadata = metadata, registry = registry, n_observed = sum(!is.na(rows$value)))
}

breeding_check_gwas_design <- function(config, n_lines) {
  n_pcs <- config$GWAS$n_pcs
  use_groups <- config$GWAS$include_source_group
  if (!is.numeric(n_pcs) || length(n_pcs) != 1L || !is.finite(n_pcs) ||
      n_pcs != floor(n_pcs) || n_pcs < 0 || n_pcs >= n_lines ||
      !is.logical(use_groups) || length(use_groups) != 1L || is.na(use_groups)) {
    stop("Provide a valid GWAS PC count and include_source_group setting.")
  }
}

breeding_run_trait <- function(reference, markers = NULL, phenotypes, trait, output,
                               gp = character(), gwas = character(), sample_metadata = NULL,
                               config = "config/breeding_pipeline.yml", seed = 2026103062L,
                               gwas_candidates = NULL, input_memory_gib = 2, plots = TRUE, unit = "trait units",
                               assessment_fraction = NULL, shared = NULL) {
  breeding_check_trait_methods(gp, gwas, gwas_candidates, assessment_fraction, plots)
  if (!is.character(trait) || length(trait) != 1L || is.na(trait) || !nzchar(trait)) stop("Supply one trait name.")
  phenotypes <- breeding_input_phenotypes(phenotypes, trait)
  if (!is.data.frame(phenotypes) || !all(c("sample_id", "trait", "value") %in% names(phenotypes)) ||
      !is.numeric(phenotypes$value)) stop("Supply one trait and phenotype columns sample_id, trait and numeric value.")
  imported <- NULL
  if (is.character(reference) && !is.matrix(reference)) {
    imported <- breeding_read_genotypes(reference, markers, maximum_memory_gib = input_memory_gib)
    reference <- imported$genotype
    markers <- imported$markers
    imported$genotype <- imported$markers <- NULL
  } else markers <- breeding_input_markers(markers)
  sample_metadata <- breeding_input_metadata(sample_metadata)
  if (is.null(shared)) gbs_validate_calls(reference)
  panel <- shared$panels[[trait]] %||% breeding_trait_panel(rownames(reference), phenotypes, trait, sample_metadata,
    if (length(gwas)) 5L else 3L, assessment_fraction, seed)
  ids <- panel$ids; observed <- panel$observed; metadata <- panel$metadata; registry <- panel$registry
  if (is.character(config)) config <- yaml::read_yaml(gbs_project_path(config))
  breeding_check_config(config)
  if (length(gwas)) breeding_check_gwas_design(config, length(ids))
  output_path <- gbs_project_path(output, must_exist = FALSE)
  if (!dir.create(output_path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Choose a new trait output folder.")
  }
  selection <- list(trait = trait, gp = gp, gwas = gwas, matched_phenotyped_lines = length(ids),
    genotype_lines_without_selected_trait = nrow(reference) - length(ids),
    observed_trait_lines_without_genotype = panel$n_observed - length(ids),
    cross_validation_performed = FALSE, assessment_fraction = assessment_fraction, seed = seed)
  if (!is.null(gwas_candidates)) selection$gwas_candidates <- gwas_candidates
  if (!is.null(imported)) {
    selection$genotype_import <- c(imported$provenance, list(excluded_markers = nrow(imported$exclusions)))
    if (nrow(imported$exclusions)) {
      path <- file.path(output, "genotype_import_exclusions.tsv")
      data.table::fwrite(imported$exclusions, path, sep = "\t")
      selection$genotype_import$exclusions_file <- path
    }
  }
  gbs_parallel_atomic(selection, file.path(output, "trait_selection.rds"))
  reference <- reference[ids, , drop = FALSE]
  recorder <- function(partition) {
    if (is.null(shared$timing_file)) return(NULL)
    breeding_timing_recorder(shared$timing_file, trait, partition)
  }
  assessment_timing <- recorder("assessment")
  assessment <- if (!is.null(registry)) breeding_timed_step(assessment_timing, "assessment_total",
    breeding_run_assessment(reference, markers, observed, metadata, registry, gp, config, seed,
      file.path(output, "assessment"), trait, unit, plots, shared, timing = assessment_timing))
  gp_file <- gwas_file <- NULL
  if (length(gp)) {
    gp_output <- file.path(output, "GP")
    final_timing <- recorder("final")
    prepared <- breeding_timed_step(final_timing, "final_gp_total",
      breeding_fit_models(reference, markers, observed, metadata, gp_output,
        learners = gp, config = config, seed = seed, shared = shared, timing = final_timing))
    gp_file <- file.path(gp_output, "breeding_model.rds")
  }
  if (length(gwas)) {
    gwas_timing <- recorder("GWAS")
    if (!length(gp)) prepared <- breeding_timed_step(gwas_timing, "training_preprocessing",
      breeding_prepare_panel(reference, markers, file.path(output, "preprocessing"), config, shared))
    gwas_output <- file.path(output, "GWAS")
    association <- breeding_timed_step(gwas_timing, "gwas_total", breeding_run_release_gwas(prepared, observed, metadata, trait, gwas_output,
      n_pcs = config$GWAS$n_pcs, include_source_group = config$GWAS$include_source_group, models = gwas, gwas_candidates = gwas_candidates,
      plots = plots, unit = unit, timing = gwas_timing))
    gwas_file <- file.path(gwas_output, "GWAS_results.tsv")
  }
  receipt <- c(selection, list(status = "complete", model_file = gp_file, gwas_file = gwas_file, assessment = assessment))
  if (!is.null(gwas_candidates)) receipt$candidate_shortlist_file <- association$candidate_shortlist_file
  if (length(gwas) && isTRUE(plots)) receipt$plot_report <- association$plot_report
  gbs_parallel_atomic(receipt, file.path(output, "trait_completion.rds"))
  receipt
}

breeding_timing_recorder <- function(path, trait = NA_character_, partition = "batch") {
  force(path); force(trait); force(partition)
  function(step, method, elapsed, status) {
    row <- data.frame(trait = trait, partition = partition, step = step, method = method,
      elapsed_seconds = unname(elapsed), status = status)
    data.table::fwrite(row, path, sep = "\t", append = file.exists(path), na = "NA")
  }
}

breeding_run_traits <- function(reference, markers = NULL, phenotypes, output, traits = NULL,
                                gp = character(), gwas = character(), sample_metadata = NULL,
                                config = "config/breeding_pipeline.yml", seed = 2026103062L,
                                gwas_candidates = NULL, input_memory_gib = 2, plots = TRUE, unit = "trait units",
                                assessment_fraction = NULL) {
  breeding_check_trait_methods(gp, gwas, gwas_candidates, assessment_fraction, plots)
  table <- breeding_input_trait_table(phenotypes, traits)
  traits <- table$traits
  if (!is.character(unit) || anyNA(unit) || any(!nzchar(unit)) ||
      (length(unit) != 1L && (anyDuplicated(names(unit)) || !setequal(names(unit), traits)))) {
    stop("Supply one unit or a named unit for each selected trait.")
  }
  units <- if (length(unit) == 1L) rep(unit, length(traits)) else unit[traits]
  if (is.character(config)) config <- yaml::read_yaml(gbs_project_path(config))
  breeding_check_config(config)
  output_path <- gbs_project_path(output, must_exist = FALSE)
  if (!dir.create(output_path, recursive = TRUE, showWarnings = FALSE)) stop("Choose a new multi-trait output folder.")
  timing_file <- file.path(output, "timing.tsv")
  batch_timing <- breeding_timing_recorder(timing_file)
  batch_start <- proc.time()[["elapsed"]]
  batch_status <- "incomplete"
  on.exit(batch_timing("batch_total", NA_character_, proc.time()[["elapsed"]] - batch_start, batch_status))
  imported <- NULL
  if (is.character(reference)) {
    imported <- breeding_timed_step(batch_timing, "genotype_import",
      breeding_read_genotypes(reference, markers, maximum_memory_gib = input_memory_gib))
    reference <- imported$genotype
    markers <- imported$markers
    imported$genotype <- imported$markers <- NULL
  } else markers <- breeding_input_markers(markers)
  gbs_validate_calls(reference)
  sample_metadata <- breeding_input_metadata(sample_metadata)
  shared <- new.env(parent = emptyenv())
  shared$preprocessors <- shared$usage <- list()
  shared$timing_file <- timing_file
  shared$panels <- breeding_timed_step(batch_timing, "panel_and_split_setup", setNames(lapply(traits, function(trait) {
    panel <- breeding_trait_panel(rownames(reference), table$phenotypes, trait, sample_metadata, if (length(gwas)) 5L else 3L,
                                  assessment_fraction, seed, sort_ids = TRUE)
    if (length(gwas)) breeding_check_gwas_design(config, length(panel$ids))
    panel
  }), traits))
  folders <- file.path(output, sprintf("%03d_%s", seq_along(traits), substr(gsub("[^A-Za-z0-9_-]+", "_", traits), 1L, 64L)))
  manifest <- data.frame(trait = traits, output = folders, unit = unname(units),
    matched_lines = vapply(shared$panels, function(x) length(x$ids), integer(1L)),
    assessment_lines = vapply(shared$panels, function(x) sum(x$registry$role == "assessment"), integer(1L)), row.names = NULL)
  manifest_file <- file.path(output, "trait_manifest.tsv")
  data.table::fwrite(manifest, manifest_file, sep = "\t")
  if (!is.null(imported) && nrow(imported$exclusions)) {
    data.table::fwrite(imported$exclusions, file.path(output, "genotype_import_exclusions.tsv"), sep = "\t")
  }
  gbs_parallel_atomic(list(genotype_import = imported$provenance, input_lines = nrow(reference), input_markers = ncol(reference),
    config = config, seed = seed, assessment_fraction = assessment_fraction), file.path(output, "multi_trait_request.rds"))
  receipts <- setNames(lapply(seq_along(traits), function(i) {
    cat(sprintf("Trait %d/%d: %s\n", i, length(traits), traits[i]))
    breeding_timed_step(breeding_timing_recorder(timing_file, traits[i], "trait"), "trait_total",
      breeding_run_trait(reference, markers, table$phenotypes, traits[i], folders[i], gp, gwas, sample_metadata,
        config, seed, gwas_candidates, plots = plots, unit = units[[i]], assessment_fraction = assessment_fraction, shared = shared))
  }), traits)
  usage_file <- file.path(output, "preprocessing_usage.tsv")
  data.table::fwrite(do.call(rbind, unname(shared$usage)), usage_file, sep = "\t")
  result <- list(status = "complete", traits = receipts, manifest_file = manifest_file, preprocessing_usage_file = usage_file,
    timing_file = timing_file,
    n_preprocessing_fits = length(shared$preprocessors), n_preprocessing_reuses = length(shared$usage) - length(shared$preprocessors))
  if (length(gp)) {
    parameters <- Map(breeding_trait_parameters, folders, receipts)
    result$hyperparameters_file <- file.path(output, "hyperparameters.tsv")
    data.table::fwrite(do.call(rbind, parameters), result$hyperparameters_file, sep = "\t")
  }
  gbs_parallel_atomic(result, file.path(output, "multi_trait_completion.rds"))
  batch_status <- "complete"
  result
}

breeding_load_model <- function(path) {
  if (dir.exists(gbs_project_path(path))) path <- file.path(path, "breeding_model.rds")
  file <- gbs_project_path(path, must_exist = FALSE)
  receipt_file <- gbs_project_path(file.path(dirname(path), "completion.rds"), must_exist = FALSE)
  if (!file.exists(file) || !file.exists(receipt_file)) stop("A completed GP fit folder is required.")
  receipt <- readRDS(receipt_file)
  if (!identical(receipt$status, "breeding_models_fitted") ||
      !identical(receipt$model_sha256, digest::digest(file = file, algo = "sha256"))) {
    stop("Saved model differs from its completion checksum.")
  }
  model <- readRDS(file)
  if (!inherits(model, "breeding_trait_model") || !identical(names(model$models), receipt$learners) ||
      !identical(model$training_ids, model$preprocessing$training_ids)) {
    stop("Saved model and completion metadata differ.")
  }
  model
}

breeding_export_imputed <- function(preprocessor, output) {
  destination <- gbs_project_path(output, must_exist = FALSE)
  if (file.exists(destination)) stop("Imputed export exists. Choose a new directory.")
  source_hash <- NULL
  if (is.character(preprocessor)) {
    path <- gbs_project_path(preprocessor)
    source_hash <- digest::digest(path, file = TRUE, algo = "sha256")
    checksum <- file.path(dirname(path), "preprocessor_sha256.txt")
    if (file.exists(checksum) && !identical(source_hash, readLines(checksum, warn = FALSE))) {
      stop("Saved preprocessing checksum does not match.")
    }
    preprocessor <- readRDS(path)
  }
  if (!inherits(preprocessor, "breeding_pipeline_preprocessor")) stop("Supply a fitted pipeline preprocessor or its RDS file.")
  calls <- preprocessor$training$genotype
  markers <- preprocessor$markers
  lines <- preprocessor$line_audit
  support <- preprocessor$training$audit
  if (!is.matrix(calls) || !identical(rownames(calls), preprocessor$training_ids) ||
      !identical(colnames(calls), preprocessor$analysis_ids) || !identical(colnames(calls), markers$marker_id) ||
      !identical(rownames(calls), lines$sample_id[lines$retained]) || !identical(rownames(calls), support$sample_id)) {
    stop("The saved imputed matrix, sample IDs and marker map disagree.")
  }
  if (!dir.create(destination, recursive = TRUE, showWarnings = FALSE)) stop("Could not create imputed export.")
  files <- file.path(output, c("imputed_genotypes.rds", "markers.tsv", "training_line_QC.tsv", "imputation_support.tsv"))
  gbs_parallel_atomic(calls, files[1L])
  data.table::fwrite(markers, files[2L], sep = "\t")
  data.table::fwrite(lines, files[3L], sep = "\t")
  data.table::fwrite(support, files[4L], sep = "\t")
  receipt <- list(matrix_file = files[1L], rows = nrow(calls), markers = ncol(calls),
    scale = "unscaled_second_allele_dosage", panel = "fitted_training_lines", preprocessor_sha256 = source_hash,
    file_sha256 = vapply(files, digest::digest, character(1L), file = TRUE, algo = "sha256"))
  gbs_parallel_atomic(receipt, file.path(output, "export_completion.rds"))
  invisible(receipt)
}

breeding_trait_parameters <- function(trait_output, receipt) {
  if (!length(receipt$model_file)) return(NULL)
  directories <- c(final = file.path(trait_output, "GP"))
  if (!is.null(receipt$assessment)) directories <- c(directories, assessment = file.path(trait_output, "assessment", "GP"))
  do.call(rbind, lapply(names(directories), function(stage) {
    model_dir <- directories[[stage]]
    file <- gbs_project_path(file.path(model_dir, "hyperparameters.tsv"), must_exist = FALSE)
    parameters <- if (file.exists(file)) utils::read.delim(file, colClasses = c(value = "character")) else {
      breeding_model_parameters(breeding_load_model(model_dir))
    }
    data.frame(trait = receipt$trait, fit = stage, parameters)
  }))
}

breeding_release_rf_explainer <- function() {
  # Use the verified matrix path for coalition predictions. Keep the saved
  # model and the functions bound to the completed studies intact.
  backend <- environment(breeding_coordinate_shapley)
  scope <- new.env(parent = backend)
  scope$breeding_predict_learner <- function(model, z) {
    if (!identical(colnames(z), model$marker_ids) || any(!is.finite(z))) {
      stop("Saved forest marker schema differs.")
    }
    value <- as.numeric(stats::predict(model$fit, data = z, num.threads = 1L)$predictions)
    if (length(value) != nrow(z) || any(!is.finite(value))) stop("Invalid saved forest predictions.")
    value
  }
  for (name in c("breeding_group_path", "breeding_group_shapley", "breeding_coordinate_shapley")) {
    fun <- get(name, envir = backend, inherits = TRUE)
    environment(fun) <- scope
    scope[[name]] <- fun
  }
  scope$breeding_coordinate_shapley
}

breeding_release_explanation_precision <- function(result) {
  metrics <- breeding_coordinate_precision(result, top = 10L)
  minimum_overlap <- min(metrics$top, 8L)
  minimum_rank <- 0.90
  maximum_ratio <- 0.10
  passed_overlap <- metrics$top_overlap >= minimum_overlap
  passed_rank <- metrics$constant_region_importance ||
    (is.finite(metrics$rank_spearman) && metrics$rank_spearman >= minimum_rank)
  passed_mcse <- metrics$top_max_local_MCSE_to_mean_abs <= maximum_ratio
  status <- if (result$paired_permutations == 0L) "exact" else
    if (passed_overlap && passed_rank && passed_mcse) "passed" else "review"
  data.frame(status = status, paired_permutations = result$paired_permutations, metrics,
    minimum_overlap = minimum_overlap, minimum_rank_spearman = minimum_rank,
    maximum_mcse_ratio = maximum_ratio, passed_overlap = passed_overlap,
    passed_rank = passed_rank, passed_mcse = passed_mcse)
}

breeding_explain_release_models <- function(model, query, markers = NULL, learners, output,
                                             query_rows = 24L, background_rows = 8L, pairs = 16L,
                                             region_width_bp = 10000000, seed = 1L, resources = NULL,
                                             input_memory_gib = 2, plots = TRUE, trait = NULL, unit = "trait units",
                                             grouping = c("coordinate", "ld"), ld_options = list()) {
  grouping <- match.arg(grouping)
  ld <- if (grouping == "ld") breeding_ld_options(ld_options)
  if (isTRUE(plots)) breeding_plot_requirements()
  fit_directory <- if (dir.exists(gbs_project_path(model))) model else dirname(model)
  fit_receipt <- gbs_project_path(file.path(dirname(fit_directory), "trait_completion.rds"), must_exist = FALSE)
  if (is.null(trait) && file.exists(fit_receipt)) trait <- readRDS(fit_receipt)$trait
  model <- breeding_load_model(model)
  learners <- breeding_release_choices(learners, names(model$models), "saved learners", allow_empty = FALSE)
  breeding_learner_status(learners)
  imported <- breeding_read_query(query, markers, model, input_memory_gib)
  query <- imported$genotype
  markers <- imported$markers
  gbs_validate_calls(query)
  valid_count <- function(x, minimum) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= minimum &&
      x <= .Machine$integer.max && x == floor(x)
  }
  if (!valid_count(query_rows, 1) || !valid_count(background_rows, 1) ||
      !valid_count(pairs, 2) || !valid_count(seed, 0)) stop("Invalid explanation sample size, permutation count or seed.")
  prep <- model$preprocessing
  if (is.null(resources)) resources <- prep$config$resources
  partition_started <- proc.time()[["elapsed"]]
  partition <- if (grouping == "ld") breeding_ld_shap_partition(prep, ld, seed) else
    breeding_coordinate_groups(prep$analysis_ids, prep$markers, prep$markers, region_width_bp)
  partition_elapsed <- proc.time()[["elapsed"]] - partition_started
  if (!dir.create(gbs_project_path(output, must_exist = FALSE), recursive = TRUE, showWarnings = FALSE)) {
    stop("Choose a new explanation output folder.")
  }
  query_ids <- with_seed_preserved(seed, {
    ids <- sort(rownames(query))
    if (length(ids) > query_rows) sort(sample(ids, query_rows)) else ids
  })
  training <- prep$training$z
  background_ids <- with_seed_preserved(seed, {
    ids <- sort(rownames(training))
    if (length(ids) > background_rows) sort(sample(ids, background_rows)) else ids
  })
  transformed <- breeding_apply_preprocessor(query[query_ids, , drop = FALSE], markers, prep,
                                             file.path(output, "checkpoints/preprocessing"), resources)
  data.table::fwrite(transformed$line_audit, file.path(output, "query_line_QC.tsv"), sep = "\t")
  if (!nrow(transformed$z)) stop("None of the sampled query lines passed saved-model QC.")
  background <- training[background_ids, , drop = FALSE]
  reports <- lapply(learners, function(learner) {
    directory <- file.path(output, learner)
    dir.create(directory)
    matrix_rf <- learner == "random_forest"
    explain <- if (matrix_rf) breeding_release_rf_explainer() else breeding_coordinate_shapley
    started <- proc.time()[["elapsed"]]
    result <- explain(model$models[[learner]], training, transformed$z,
      background, partition$groups, pairs, seed, file.path(directory, "pair_checkpoints"))
    elapsed <- proc.time()[["elapsed"]] - started
    gbs_parallel_atomic(result, file.path(directory, "explanations.rds"))
    list(timing = data.frame(learner = learner, elapsed_seconds = elapsed),
      importance = data.frame(learner = learner, result$summary),
      precision = data.frame(learner = learner, method = result$method,
        prediction_backend = if (matrix_rf) "ranger_matrix" else "saved_learner",
        breeding_release_explanation_precision(result)))
  })
  importance <- do.call(rbind, lapply(reports, `[[`, "importance"))
  precision <- do.call(rbind, lapply(reports, `[[`, "precision"))
  precision_file <- file.path(output, "explanation_precision.tsv")
  data.table::fwrite(do.call(rbind, lapply(reports, `[[`, "timing")), file.path(output, "explanation_timing.tsv"), sep = "\t")
  data.table::fwrite(importance, file.path(output, "regional_importance.tsv"), sep = "\t")
  data.table::fwrite(precision, precision_file, sep = "\t")
  gbs_parallel_atomic(prep$markers, file.path(output, "report_markers.rds"))
  gbs_parallel_atomic(partition$coverage, file.path(output, "report_regions.rds"))
  if (grouping == "ld") {
    data.table::fwrite(partition$membership, file.path(output, "ld_group_members.tsv"), sep = "\t")
    data.table::fwrite(partition$distances, file.path(output, "ld_distances.tsv"), sep = "\t")
    data.table::fwrite(partition$bins, file.path(output, "ld_decay.tsv"), sep = "\t", na = "NA")
  }
  receipt <- list(learners = learners, sampled_query_ids = query_ids, explained_query_ids = rownames(transformed$z),
    background_ids = background_ids, pairs_for_nonlinear_learners = pairs,
    region_width_bp = if (grouping == "coordinate") region_width_bp else NA_real_,
    seed = seed, phenotype_values_read = FALSE, models_refitted = FALSE,
    precision_file = precision_file, precision = precision, trait = trait, grouping = grouping,
    grouping_elapsed_seconds = partition_elapsed, ld_options = ld,
    ld_distance_reference = if (grouping == "ld") "https://www.cog-genomics.org/plink/2.0/ld")
  gbs_parallel_atomic(receipt, file.path(output, "explanation_completion.rds"))
  if (isTRUE(plots)) receipt$plot_report <- breeding_plot_report(file.path(output, "plots"), explanation_output = output, unit = unit)
  receipt
}
