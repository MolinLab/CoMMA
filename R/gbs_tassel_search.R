# SPDX-License-Identifier: LGPL-2.1-only
# Calibration acceleration for the TASSEL 5.2.96 biallelic native adaptation.
# Algorithm attribution and license: THIRD_PARTY_NOTICES.md, licenses/LGPL-2.1.txt.
# Original reference implementation remains in gbs_tassel_native.R unchanged.

build_gbs_tassel_rank_cache <- function(model, query, high_ld_sites, target_cells,
                                         maximum_neighbours = 200L) {
  if (!inherits(model, "gbs_tassel_native")) stop("Invalid fitted reference model.")
  valid_integer <- function(x, lo, hi) length(x) == 1L && !is.na(x) && is.finite(x) &&
    x == floor(x) && x >= lo && x <= hi
  if (!valid_integer(high_ld_sites, 2L, model$max_ld_sites) ||
      !valid_integer(maximum_neighbours, 2L, 200L)) stop("Invalid LD/cache neighbour limits.")
  if (!is.data.frame(target_cells) || !all(c("sample_id", "marker_id") %in% names(target_cells)) ||
      !nrow(target_cells) || anyNA(target_cells[c("sample_id", "marker_id")]) ||
      anyDuplicated(target_cells[c("sample_id", "marker_id")])) stop("Invalid calibration target IDs.")
  estimate <- model$provenance$estimated_fit_memory_gib +
    (as.double(length(query)) * 8 * 5 + as.double(nrow(target_cells)) * maximum_neighbours * 24) / 1024^3 + 0.10
  if (!is.finite(estimate) || estimate > model$settings$maximum_memory_gib) stop("Rank cache exceeds its memory guard.")
  gbs_validate_calls(query)
  columns <- match(model$marker_ids, colnames(query))
  if (anyNA(columns)) stop("Query lacks fitted marker IDs.")
  shared <- intersect(rownames(model$reference), rownames(query))
  if (length(shared)) {
    old <- unname(model$reference[shared, , drop = FALSE])
    new <- unname(query[shared, columns, drop = FALSE])
    if (!identical(is.na(old), is.na(new)) || any(old[!is.na(old)] != new[!is.na(new)])) stop("Overlapping IDs require identical masked calls.")
  }
  target <- target_cells[c("sample_id", "marker_id")]
  target <- target[order(target$marker_id, target$sample_id), , drop = FALSE]
  rownames(target) <- NULL
  r <- match(target$sample_id, rownames(query)); c <- match(target$marker_id, colnames(query))
  j_target <- match(target$marker_id, model$marker_ids)
  if (anyNA(r) || anyNA(c) || anyNA(j_target) || any(!is.na(query[cbind(r, c)]))) stop("Targets must be missing query cells with fitted IDs.")
  query_aligned <- query[, columns, drop = FALSE]
  donors <- vector("list", nrow(target))
  target$eligible_donors <- 0L
  target$fixed_hardcall <- NA_real_
  target$disposition <- "abstained"
  started <- proc.time()
  distance_calculations <- 0L
  for (j in unique(j_target)) {
    positions <- which(j_target == j)
    frequency <- model$marker_statistics$alt_frequency[[j]]
    if (is.finite(frequency) && frequency %in% c(0, 1)) {
      target$fixed_hardcall[positions] <- 2 * frequency
      target$disposition[positions] <- "monomorphic"
      next
    }
    predictors <- head(model$neighbours[[j]], high_ld_sites)
    predictors <- predictors[order(model$markers$chromosome[predictors],
      model$markers$position_bp[predictors], model$marker_ids[predictors])]
    if (!length(predictors)) next
    donor_values <- model$reference[, predictors, drop = FALSE]
    donor_coverage <- rowMeans(!is.na(donor_values))
    for (position in positions) {
      qvalues <- query_aligned[r[[position]], predictors]
      expected <- donor_coverage * mean(!is.na(qvalues)) * length(predictors)
      eligible <- which(expected > 10 & !is.na(model$reference[, j]))
      if (!length(eligible)) next
      distances <- gbs_tassel_distance(qvalues, donor_values[eligible, , drop = FALSE])$distance
      distance_calculations <- distance_calculations + 1L
      good <- is.finite(distances)
      eligible <- eligible[good]; distances <- distances[good]
      target$eligible_donors[[position]] <- length(eligible)
      if (!length(eligible)) next
      calls <- model$reference[eligible, j]
      rank <- head(order(distances, model$byte_codes[j, calls + 1L], rownames(model$reference)[eligible]), maximum_neighbours)
      donors[[position]] <- list(calls = calls[rank], distances = distances[rank], byte_codes = model$byte_codes[j, ])
      target$disposition[[position]] <- "ld_knn"
    }
  }
  structure(list(query = query, targets = target, donors = donors, high_ld_sites = high_ld_sites,
    maximum_neighbours = maximum_neighbours,
    provenance = list(reference_sha256 = model$provenance$fitted_reference_sha256,
      query_sha256 = digest::digest(query, algo = "sha256"),
      target_ID_sha256 = digest::digest(target[c("sample_id", "marker_id")], algo = "sha256"),
      estimated_memory_gib = estimate, maximum_memory_gib = model$settings$maximum_memory_gib,
      distance_calculations = distance_calculations,
      elapsed_seconds = unname((proc.time() - started)[["elapsed"]]),
      CPU_seconds = unname((proc.time() - started)[["user.self"]] + (proc.time() - started)[["sys.self"]]),
      truth_values_received = FALSE, phenotype_used = FALSE,
      scope = "masked_calibration_targets_only_not_final_full_imputation")), class = "gbs_tassel_rank_cache")
}

predict_gbs_tassel_rank_cache <- function(cache, knn_taxa) {
  if (!inherits(cache, "gbs_tassel_rank_cache") || length(knn_taxa) != 1L || !is.finite(knn_taxa) ||
      knn_taxa != floor(knn_taxa) || knn_taxa < 2 || knn_taxa > cache$maximum_neighbours) stop("Invalid cached prediction setting.")
  if (digest::digest(cache$query, algo = "sha256") != cache$provenance$query_sha256) stop("Cached masked query changed.")
  result <- cache$query
  for (i in seq_len(nrow(cache$targets))) {
    value <- cache$targets$fixed_hardcall[[i]]
    if (!is.null(cache$donors[[i]])) {
      donor <- cache$donors[[i]]
      selected <- seq_len(min(knn_taxa, length(donor$calls)))
      value <- gbs_tassel_vote(donor$calls[selected], donor$distances[selected], cache$high_ld_sites, donor$byte_codes)
    }
    result[cache$targets$sample_id[[i]], cache$targets$marker_id[[i]]] <- value
  }
  gbs_validate_imputation_output(result, cache$query)
  attr(result, "imputation_diagnostics") <- list(scope = cache$provenance$scope,
    counts = table(factor(cache$targets$disposition, levels = c("ld_knn", "monomorphic", "abstained"))),
    requested_k = knn_taxa, available_donor_counts = cache$targets$eligible_donors,
    no_distance_recomputation = TRUE, truth_values_received = FALSE)
  result
}

gbs_search_mask_plans <- function(query, statistics, settings, excluded_cells = NULL) {
  gbs_validate_calls(query)
  repeats <- settings$repeats
  if (length(repeats) != 1L || !is.finite(repeats) || repeats < 2 || repeats > 5 || repeats != floor(repeats) ||
      !is.finite(settings$fraction_per_phase) || settings$fraction_per_phase <= 0 ||
      settings$fraction_per_phase * 2 * repeats >= 1 || settings$minimum_remaining_calls < 2 ||
      length(settings$seed) != 1L || !is.finite(settings$seed) || settings$seed < 0 || settings$seed > .Machine$integer.max) stop("Invalid repeated-mask settings.")
  estimate <- as.double(length(query)) * 8 * (repeats + 5) / 1024^3 + 0.10
  if (estimate > settings$maximum_memory_gib) stop("Search masking exceeds memory guard.")
  statistics <- statistics[match(colnames(query), statistics$marker_id), , drop = FALSE]
  if (anyNA(statistics$marker_id) || anyDuplicated(statistics$marker_id)) stop("Reference statistics are not ID-aligned.")
  exclusions <- character()
  if (!is.null(excluded_cells)) {
    if (!is.data.frame(excluded_cells) || !all(c("sample_id", "marker_id") %in% names(excluded_cells)) ||
        anyNA(excluded_cells[c("sample_id", "marker_id")])) stop("Invalid excluded target IDs.")
    exclusions <- unique(paste(excluded_cells$sample_id, excluded_cells$marker_id, sep = "\037"))
  }
  pieces <- list()
  with_seed_preserved(settings$seed, {
    for (id in sort(colnames(query))) {
      j <- match(id, colnames(query))
      frequency <- statistics$alt_frequency[[j]]
      if (!is.finite(frequency) || frequency <= 0 || frequency >= 1) next
      observed <- which(!is.na(query[, j]))
      available <- observed[!paste(rownames(query)[observed], id, sep = "\037") %in% exclusions]
      available <- available[order(rownames(query)[available])]
      k <- min(max(1L, floor(length(observed) * settings$fraction_per_phase)),
        floor(length(available) / (2 * repeats)),
        floor((length(observed) - settings$minimum_remaining_calls) / (2 * repeats)))
      if (k < 1L) next
      hash <- digest::digest(paste(settings$seed, id, sep = "|"), algo = "xxhash32", serialize = FALSE)
      set.seed(as.integer(strtoi(substr(hash, 1, 7), 16)))
      chosen <- available[sample.int(length(available), 2 * repeats * k)]
      values <- query[chosen, j]
      pieces[[length(pieces) + 1L]] <- data.frame(sample_id = rownames(query)[chosen], marker_id = id,
        truth = values, phase = rep(c("tuning", "verification"), repeats, each = k),
        repeat_id = rep(seq_len(repeats), each = 2 * k), reference_maf = statistics$maf[[j]],
        reference_mac = statistics$mac[[j]],
        nonmajor_genotype = if (frequency <= 0.5) values > 0 else values < 2)
    }
  })
  if (!length(pieces)) stop("No markers can support fresh protected masks.")
  ledger <- do.call(rbind, pieces)
  keys <- with(ledger, paste(sample_id, marker_id, sep = "\037"))
  if (anyDuplicated(keys) || any(keys %in% exclusions)) stop("Mask phase/repeat/exclusion collision.")
  if (any(table(ledger$phase, ledger$repeat_id) > settings$maximum_masks_per_phase)) stop("Mask cap exceeded; choose a smaller panel explicitly.")
  verification <- ledger[ledger$phase == "verification", , drop = FALSE]
  plans <- lapply(seq_len(repeats), function(repeat_id) {
    truth <- ledger[ledger$repeat_id == repeat_id, , drop = FALSE]
    hidden <- rbind(truth[truth$phase == "tuning", , drop = FALSE], verification)
    masked <- query
    masked[cbind(match(hidden$sample_id, rownames(query)), match(hidden$marker_id, colnames(query)))] <- NA_real_
    list(masked = masked, truth = truth, repeat_id = repeat_id)
  })
  list(plans = plans, full_ledger = ledger, provenance = list(seed = settings$seed,
    all_verification_cells_hidden_in_every_context = TRUE, excluded_target_count = length(exclusions),
    excluded_cells_remain_ordinary_observed_predictors_when_present = TRUE,
    estimated_memory_gib = estimate, natural_missing_truth_known = FALSE))
}

gbs_search_grid <- function(ld, k) {
  if (!length(ld) || !length(k) || anyNA(c(ld, k)) || any(!is.finite(c(ld, k))) ||
      any(c(ld, k) != floor(c(ld, k))) || any(c(ld, k) < 2) || any(c(ld, k) > 200) || anyDuplicated(ld) || anyDuplicated(k)) stop("Invalid bounded search axes.")
  grid <- expand.grid(ld = as.integer(ld), k = as.integer(k))
  grid$id <- sprintf("LD%03d_k%03d", grid$ld, grid$k)
  grid[order(grid$ld, grid$k), , drop = FALSE]
}

gbs_search_select <- function(tuning, grid, settings) {
  required <- c("candidate_id", "repeat_id", "n_masked", "n_imputed", "coverage", "hardcall_error_called")
  if (!all(required %in% names(tuning)) || !setequal(unique(tuning$candidate_id), grid$id) ||
      anyDuplicated(tuning[c("candidate_id", "repeat_id")]) || nrow(tuning) != nrow(grid) * settings$repeats) stop("Incomplete/duplicate tuning grid.")
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    x <- tuning[tuning$candidate_id == grid$id[[i]], ]
    if (!setequal(x$repeat_id, seq_len(settings$repeats))) stop("Tuning repeat identities differ.")
    eligible <- all(is.finite(x$coverage)) && all(x$coverage >= settings$minimum_coverage_each_repeat) &&
      all(is.finite(x$hardcall_error_called)) && sum(x$n_imputed) > 0
    data.frame(candidate_id = grid$id[[i]], ld = grid$ld[[i]], k = grid$k[[i]], eligible = eligible,
      pooled_error = if (eligible) sum(x$n_imputed * x$hardcall_error_called) / sum(x$n_imputed) else NA_real_,
      worst_repeat_error = if (eligible) max(x$hardcall_error_called) else NA_real_)
  })
  table <- do.call(rbind, rows)
  valid <- table[table$eligible, , drop = FALSE]
  selected <- if (nrow(valid)) valid[order(valid$pooled_error, valid$ld, valid$k), , drop = FALSE][1, ] else NULL
  table$near_best <- if (is.null(selected)) FALSE else table$eligible &
    table$pooled_error <= selected$pooled_error + settings$near_best_absolute_error_tolerance
  table$near_best[is.na(table$near_best)] <- FALSE
  list(selected = selected, table = table, verification_used_for_selection = FALSE,
    objective = "pooled_tuning_hardcall_error_after_per_repeat_coverage_gate",
    near_best_band_is_practical_tolerance_not_confidence_interval = TRUE)
}

gbs_search_refinement <- function(coarse, selection, limits) {
  if (is.null(selection$selected)) return(coarse[FALSE, , drop = FALSE])
  neighbours <- function(value, axis, maximum) {
    axis <- sort(unique(axis)); position <- match(value, axis)
    lower <- if (position == 1L) max(2L, floor(value / 2)) else floor((axis[[position - 1L]] + value) / 2)
    upper <- if (position == length(axis)) min(maximum, value * 2L) else ceiling((value + axis[[position + 1L]]) / 2)
    sort(unique(c(lower, value, upper)))
  }
  more <- gbs_search_grid(neighbours(selection$selected$ld, coarse$ld, limits$maximum_LD_sites),
    neighbours(selection$selected$k, coarse$k, limits$maximum_neighbours))
  more[!more$id %in% coarse$id, , drop = FALSE]
}

tune_gbs_tassel_neighbours <- function(model, query, config, excluded_cells = NULL,
                                        output_dir = NULL, progress = NULL) {
  if (!inherits(model, "gbs_tassel_native")) stop("Invalid fitted reference model.")
  gbs_validate_calls(query)
  if (length(intersect(rownames(model$reference), rownames(query)))) stop("Tuning requires disjoint reference and query IDs.")
  if (anyNA(match(model$marker_ids, colnames(query)))) stop("Query lacks fitted markers.")
  query <- query[, model$marker_ids, drop = FALSE]
  search <- config$search
  if (!identical(search$objective, "pooled_tuning_hardcall_error") ||
      !identical(search$tie_break, "smaller_LD_then_smaller_k")) stop("Unsupported selection objective or tie policy.")
  if (search$maximum_LD_sites > model$max_ld_sites || search$maximum_neighbours > 200L ||
      search$maximum_candidates > 64L || config$resources$global_cpu_budget != 1L ||
      config$resources$maximum_memory_gib > model$settings$maximum_memory_gib ||
      search$minimum_coverage_each_repeat <= 0 || search$minimum_coverage_each_repeat > 1 ||
      search$near_best_absolute_error_tolerance < 0) stop("Invalid search limits or allocation.")
  if (!is.null(progress) && !is.function(progress)) stop("Progress must be a function.")
  if (!is.null(output_dir)) {
    output_dir <- gbs_project_path(output_dir, must_exist = FALSE)
    if (dir.exists(output_dir)) stop("Search output exists; choose a new directory.")
    dir.create(output_dir, recursive = TRUE)
    dir.create(file.path(output_dir, "rank_cache")); dir.create(file.path(output_dir, "predictions"))
  }
  mask_settings <- config$masking
  mask_settings$maximum_memory_gib <- config$resources$maximum_memory_gib
  masks <- gbs_search_mask_plans(query, model$marker_statistics, mask_settings, excluded_cells)
  selection_settings <- search; selection_settings$repeats <- config$masking$repeats
  coarse <- gbs_search_grid(search$high_LD_sites, search$nearest_neighbor_taxa)
  if (nrow(coarse) > search$maximum_candidates) stop("Coarse grid exceeds candidate budget.")
  cache_store <- new.env(parent = emptyenv()); prediction_store <- new.env(parent = emptyenv())
  cache_bytes <- 0; prediction_bytes <- 0
  timing <- list(); rows <- list(); started <- proc.time()
  check_time <- function() {
    if ((proc.time() - started)[["elapsed"]] > config$resources$maximum_elapsed_seconds) stop("Search elapsed budget reached; incomplete grid cannot select settings.")
  }
  prediction_key <- function(id, repeat_id) paste0(id, "_repeat", repeat_id)
  save_prediction <- function(key, prediction) {
    if (is.null(output_dir)) {
      prediction_bytes <<- prediction_bytes + as.numeric(object.size(prediction))
      if ((cache_bytes + prediction_bytes) / 1024^3 + model$provenance$estimated_fit_memory_gib > config$resources$maximum_memory_gib) stop("In-memory search exceeds allocation; use disk checkpoints.")
      assign(key, prediction, envir = prediction_store)
    } else saveRDS(prediction, file.path(output_dir, "predictions", paste0(key, ".rds")))
  }
  get_prediction <- function(key) if (is.null(output_dir)) get(key, envir = prediction_store) else
    readRDS(file.path(output_dir, "predictions", paste0(key, ".rds")))
  evaluate <- function(grid, stage) {
    if (!nrow(grid)) return(invisible(NULL))
    for (ld in sort(unique(grid$ld))) for (repeat_id in seq_len(config$masking$repeats)) {
      check_time()
      key <- paste0("LD", ld, "_repeat", repeat_id)
      plan <- masks$plans[[repeat_id]]
      cache_path <- if (!is.null(output_dir)) file.path(output_dir, "rank_cache", paste0(key, ".rds")) else NULL
      available <- if (is.null(output_dir)) exists(key, envir = cache_store, inherits = FALSE) else file.exists(cache_path)
      if (available) cache <- if (is.null(output_dir)) get(key, envir = cache_store) else readRDS(cache_path) else {
        cache <- build_gbs_tassel_rank_cache(model, plan$masked, ld, plan$truth[c("sample_id", "marker_id")], search$maximum_neighbours)
        timing[[length(timing) + 1L]] <<- data.frame(ld = ld, repeat_id = repeat_id,
          rank_elapsed_seconds = cache$provenance$elapsed_seconds, rank_CPU_seconds = cache$provenance$CPU_seconds,
          distance_calculations = cache$provenance$distance_calculations,
          target_cells = nrow(cache$targets), cache_bytes = as.numeric(object.size(cache)))
        if (is.null(output_dir)) {
          cache_bytes <<- cache_bytes + as.numeric(object.size(cache))
          if ((cache_bytes + prediction_bytes) / 1024^3 + model$provenance$estimated_fit_memory_gib > config$resources$maximum_memory_gib) stop("Rank-cache memory budget reached; use disk checkpoints.")
          assign(key, cache, envir = cache_store)
        } else saveRDS(cache, cache_path)
      }
      for (i in which(grid$ld == ld)) {
        check_time()
        vote_started <- proc.time()
        prediction <- predict_gbs_tassel_rank_cache(cache, grid$k[[i]])
        vote_seconds <- unname((proc.time() - vote_started)[["elapsed"]])
        save_prediction(prediction_key(grid$id[[i]], repeat_id), prediction)
        # Only tuning truth is scored until coarse and refinement selection freeze.
        truth <- plan$truth[plan$truth$phase == "tuning", , drop = FALSE]
        score <- gbs_score_imputation(prediction, truth)$overall
        row <- data.frame(candidate_id = grid$id[[i]], ld = ld, k = grid$k[[i]], repeat_id = repeat_id,
          stage = stage, score, vote_elapsed_seconds = vote_seconds)
        rows[[length(rows) + 1L]] <<- row
        if (!is.null(progress)) progress(row)
      }
    }
  }
  evaluate(coarse, "coarse")
  coarse_selection <- gbs_search_select(do.call(rbind, rows), coarse, selection_settings)
  refinement <- if (isTRUE(search$refinement_enabled)) gbs_search_refinement(coarse, coarse_selection, search) else coarse[FALSE, ]
  combined <- rbind(coarse, refinement)
  if (nrow(combined) > search$maximum_candidates) stop("Refinement exceeds candidate budget.")
  if (!is.null(output_dir)) saveRDS(list(coarse_selection = coarse_selection, refinement_grid = refinement,
    verification_used = FALSE), file.path(output_dir, "refinement_decision.rds"))
  evaluate(refinement, "refinement")
  tuning <- do.call(rbind, rows)
  selection <- gbs_search_select(tuning, combined, selection_settings)
  if (!is.null(output_dir)) saveRDS(selection, file.path(output_dir, "frozen_selection.rds"))
  verification <- list()
  if (!is.null(selection$selected)) {
    for (repeat_id in seq_len(config$masking$repeats)) {
      truth <- masks$plans[[repeat_id]]$truth
      truth <- truth[truth$phase == "verification", , drop = FALSE]
      prediction <- get_prediction(prediction_key(selection$selected$candidate_id, repeat_id))
      verification[[repeat_id]] <- gbs_score_imputation(prediction, truth)
    }
  }
  result <- list(status = if (is.null(selection$selected)) "no_candidate_met_coverage_gate" else "coarse_refined_selection_and_verification_complete",
    selection = selection, coarse_selection = coarse_selection, coarse_grid = coarse, refinement_grid = refinement,
    tuning = tuning, mask_plan = masks, verification = verification, rank_timing = do.call(rbind, timing),
    provenance = list(reference_sha256 = model$provenance$fitted_reference_sha256,
      query_sha256 = digest::digest(query, algo = "sha256"), config_sha256 = digest::digest(config, algo = "sha256"),
      phenotype_used = FALSE, rare_recovery_used_as_selection_gate = FALSE,
      clustering_used = FALSE, target_only_calibration = TRUE,
      elapsed_seconds = unname((proc.time() - started)[["elapsed"]]),
      selected_at_search_boundary = if (is.null(selection$selected)) NA else
        selection$selected$ld %in% c(2L, search$maximum_LD_sites) || selection$selected$k %in% c(2L, search$maximum_neighbours)))
  if (!is.null(output_dir)) {
    saveRDS(result, file.path(output_dir, "search_result.rds"))
    utils::write.table(tuning, file.path(output_dir, "tuning.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    utils::write.table(selection$table, file.path(output_dir, "selection.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    utils::write.table(masks$full_ledger, file.path(output_dir, "mask_ledger.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    utils::write.table(result$rank_timing, file.path(output_dir, "rank_timing.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  }
  result
}
