# Backend-independent genotype masking, calibration and verification.
# This module never launches an external executable or reads phenotypes.

gbs_validate_calls <- function(x) {
  if (!is.matrix(x) || !is.numeric(x) || !length(x) ||
      is.null(rownames(x)) || is.null(colnames(x))) {
    stop("Genotypes must be a nonempty numeric matrix with sample and marker IDs.")
  }
  for (ids in dimnames(x)) {
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
        any(grepl("\037", ids, fixed = TRUE))) stop("Invalid genotype IDs.")
  }
  observed <- x[!is.na(x)]
  if (any(!is.finite(observed)) || any(!observed %in% 0:2)) {
    stop("Masking requires original observed 0/1/2 calls, with missing calls as NA.")
  }
  invisible(TRUE)
}

gbs_make_mask_plan <- function(genotype, fraction_per_phase = 0.05,
                               seed = 2026103001L,
                               minimum_remaining_calls = 5L,
                               maximum_masks_per_phase = 50000L,
                               maximum_memory_gib = 4,
                               reference_statistics = NULL) {
  # Bound allocations before validation or making matrix copies.
  estimated_gib <- (as.double(length(genotype)) * 8 * 5) / 1024^3 + 0.10
  if (!is.finite(maximum_memory_gib) || maximum_memory_gib <= 0 ||
      estimated_gib > maximum_memory_gib) stop("Masking memory budget exceeded.")
  gbs_validate_calls(genotype)
  if (length(fraction_per_phase) != 1L || !is.finite(fraction_per_phase) ||
      fraction_per_phase <= 0 || fraction_per_phase >= 0.5 ||
      length(seed) != 1L || !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
      minimum_remaining_calls < 2L || maximum_masks_per_phase < 1L) {
    stop("Invalid masking fraction, seed, or call limits.")
  }
  x <- genotype[order(rownames(genotype)), order(colnames(genotype)), drop = FALSE]
  count <- colSums(!is.na(x))
  frequency <- colMeans(x, na.rm = TRUE) / 2
  reference_count <- count
  if (!is.null(reference_statistics)) {
    needed <- c("marker_id", "alt_frequency", "n_observed")
    if (!is.data.frame(reference_statistics) ||
        !all(needed %in% names(reference_statistics)) || anyDuplicated(reference_statistics$marker_id)) {
      stop("Invalid independent reference statistics.")
    }
    align <- match(colnames(x), reference_statistics$marker_id)
    if (anyNA(align)) stop("Mask markers lack independent reference statistics.")
    frequency <- reference_statistics$alt_frequency[align]
    reference_count <- reference_statistics$n_observed[align]
    if (any(!is.finite(frequency)) || any(frequency < 0 | frequency > 1) ||
        any(!is.finite(reference_count)) || any(reference_count < 1)) stop("Invalid reference frequencies/counts.")
  }
  eligible <- which(count >= minimum_remaining_calls + 2L &
    is.finite(frequency) & frequency > 0 & frequency < 1)
  if (!length(eligible)) stop("No polymorphic marker has enough observed calls to mask.")
  truth <- with_seed_preserved(seed, {
    if (length(eligible) > maximum_masks_per_phase) {
      eligible <- sort(sample(eligible, maximum_masks_per_phase))
    }
    per_marker_cap <- max(1L, floor(maximum_masks_per_phase / length(eligible)))
    pieces <- vector("list", length(eligible))
    for (j in seq_along(eligible)) {
      marker <- eligible[[j]]
      marker_hash <- digest::digest(paste(seed, colnames(x)[[marker]], sep = "|"),
        algo = "xxhash32", serialize = FALSE)
      marker_seed <- as.integer((as.double(seed) + strtoi(substr(marker_hash, 1L, 7L), 16L)) %% .Machine$integer.max)
      set.seed(marker_seed)
      observed <- which(!is.na(x[, marker]))
      k <- min(max(1L, floor(length(observed) * fraction_per_phase)),
        floor((length(observed) - minimum_remaining_calls) / 2), per_marker_cap)
      chosen <- sample(observed, 2L * k)
      value <- x[chosen, marker]
      pieces[[j]] <- data.frame(sample_id = rownames(x)[chosen],
        marker_id = colnames(x)[[marker]], truth = value,
        phase = rep(c("tuning", "verification"), each = k),
        reference_maf = min(frequency[[marker]], 1 - frequency[[marker]]),
        reference_mac = 2 * reference_count[[marker]] * min(frequency[[marker]], 1 - frequency[[marker]]),
        nonmajor_genotype = if (frequency[[marker]] <= 0.5) value > 0 else value < 2,
        stringsAsFactors = FALSE)
    }
    do.call(rbind, pieces)
  })
  keys <- paste(truth$sample_id, truth$marker_id, sep = "\037")
  if (anyDuplicated(keys)) stop("Tuning and verification masks overlap.")
  masked <- genotype
  masked[cbind(match(truth$sample_id, rownames(masked)),
    match(truth$marker_id, colnames(masked)))] <- NA_real_
  structure(list(masked = masked, truth = truth,
    provenance = list(seed = as.integer(seed), fraction_per_phase = fraction_per_phase,
      minimum_remaining_calls = minimum_remaining_calls,
      marker_ID_seeded_masks = TRUE,
      mask_type = "observed_sample_by_marker_calls_stratified_by_marker",
      reference_MAF_MAC_for_evaluation_only = TRUE,
      natural_missing_genotypes_not_known = TRUE,
      frequency_source = if (is.null(reference_statistics)) "masked_cohort_observed_calls_for_evaluation_only" else
        "independent_reference_lines",
      estimated_memory_gib = estimated_gib,
      input_sha256 = digest::digest(genotype, algo = "sha256"),
      mask_sha256 = digest::digest(truth, algo = "sha256"))),
    class = "gbs_imputation_mask_plan")
}

gbs_validate_imputation_output <- function(prediction, masked) {
  if (!is.matrix(prediction) || !is.numeric(prediction) ||
      !identical(dimnames(prediction), dimnames(masked)) ||
      !identical(dim(prediction), dim(masked))) stop("Imputation output IDs or dimensions changed.")
  values <- prediction[!is.na(prediction)]
  if (any(!is.finite(values)) || any(values < 0 | values > 2)) {
    stop("Imputed dosages must be finite values in [0,2], or NA abstentions.")
  }
  known <- !is.na(masked)
  if (anyNA(prediction[known]) || any(prediction[known] != masked[known])) {
    stop("The imputer changed an observed unmasked genotype call.")
  }
  invisible(TRUE)
}

gbs_score_imputation <- function(prediction, truth) {
  if (!is.matrix(prediction) || !is.numeric(prediction) ||
      is.null(rownames(prediction)) || is.null(colnames(prediction)) ||
      anyDuplicated(rownames(prediction)) || anyDuplicated(colnames(prediction)) ||
      any(!is.finite(prediction[!is.na(prediction)]))) stop("Invalid prediction matrix.")
  required <- c("sample_id", "marker_id", "truth", "reference_maf", "nonmajor_genotype")
  if (!is.data.frame(truth) || !nrow(truth) || !all(required %in% names(truth)) ||
      anyNA(truth[, required]) || any(!truth$truth %in% 0:2) ||
      anyDuplicated(paste(truth$sample_id, truth$marker_id, sep = "\037"))) stop("Invalid mask truth ledger.")
  rows <- match(truth$sample_id, rownames(prediction))
  cols <- match(truth$marker_id, colnames(prediction))
  if (anyNA(rows) || anyNA(cols)) stop("Mask IDs absent from prediction.")
  predicted <- prediction[cbind(rows, cols)]
  valid <- is.finite(predicted)
  if (any(predicted[valid] < 0 | predicted[valid] > 2)) stop("Invalid predicted dosage.")
  hard_call <- floor(predicted + 0.5)
  correct <- valid & hard_call == truth$truth
  correct[is.na(correct)] <- FALSE
  summarize <- function(idx) {
    called <- idx[valid[idx]]
    nonmajor <- idx[truth$nonmajor_genotype[idx]]
    class_recall <- vapply(sort(unique(truth$truth[idx])), function(g)
      mean(correct[idx[truth$truth[idx] == g]]), numeric(1L))
    data.frame(n_masked = length(idx), n_imputed = length(called),
      coverage = mean(valid[idx]),
      dosage_rmse_called = if (length(called)) sqrt(mean((predicted[called] - truth$truth[called])^2)) else NA_real_,
      dosage_mae_called = if (length(called)) mean(abs(predicted[called] - truth$truth[called])) else NA_real_,
      hardcall_error_called = if (length(called)) mean(!correct[called]) else NA_real_,
      genotype_class_macro_recall = mean(class_recall),
      n_nonmajor = length(nonmajor),
      nonmajor_genotype_recall = if (length(nonmajor)) mean(correct[nonmajor]) else NA_real_,
      low_nonmajor_mask_count = length(nonmajor) < 20L,
      low_mask_count = length(idx) < 30L)
  }
  overall <- summarize(seq_len(nrow(truth)))
  bins <- cut(truth$reference_maf, c(-Inf, 0.01, 0.05, Inf), right = FALSE,
    labels = c("below_0.01", "0.01_to_below_0.05", "0.05_or_more"))
  strata <- do.call(rbind, lapply(levels(bins), function(bin) {
    idx <- which(bins == bin)
    if (!length(idx)) return(NULL)
    data.frame(maf_bin = bin, summarize(idx))
  }))
  by_mac <- if ("reference_mac" %in% names(truth)) {
    mac_bins <- cut(truth$reference_mac, c(-Inf, 5, 20, Inf), right = FALSE,
      labels = c("below_5", "5_to_below_20", "20_or_more"))
    do.call(rbind, lapply(levels(mac_bins), function(bin) {
      idx <- which(mac_bins == bin)
      if (!length(idx)) return(NULL)
      data.frame(mac_bin = bin, summarize(idx))
    }))
  } else NULL
  by_genotype <- do.call(rbind, lapply(sort(unique(truth$truth)), function(g)
    data.frame(true_genotype = g, summarize(which(truth$truth == g)))))
  marker_groups <- split(seq_len(nrow(truth)), truth$marker_id)
  per_marker <- do.call(rbind, lapply(names(marker_groups), function(marker) {
    idx <- marker_groups[[marker]]
    idx <- idx[valid[idx]]
    r <- if (length(idx) >= 3L && stats::sd(predicted[idx]) > 0 &&
      stats::sd(truth$truth[idx]) > 0) stats::cor(predicted[idx], truth$truth[idx]) else NA_real_
    data.frame(marker_id = marker, n_called_masks = length(idx),
      dosage_r = r, dosage_r_squared = r^2)
  }))
  list(overall = overall, by_maf = strata, by_mac = by_mac,
    by_genotype = by_genotype, by_marker = per_marker)
}

gbs_impute_training_mean <- function(masked, settings = list()) {
  gbs_validate_calls(masked)
  n <- colSums(!is.na(masked))
  frequency <- colMeans(masked, na.rm = TRUE) / 2
  maf <- pmin(frequency, 1 - frequency)
  mac <- 2 * n * maf
  heterozygosity <- colSums(masked == 1, na.rm = TRUE) / n
  keep <- is.finite(maf) & n > 0 &
    maf > (settings$maf_min %||% 0) & mac >= (settings$min_mac %||% 0) &
    colMeans(is.na(masked)) < (settings$marker_missing_max %||% 1) &
    heterozygosity < (settings$marker_heterozygosity_max %||% 1.01)
  output <- masked
  for (j in which(keep)) output[is.na(output[, j]), j] <- 2 * frequency[[j]]
  attr(output, "retained_marker_ids") <- colnames(masked)[keep]
  output
}

gbs_build_imputation_candidates <- function(config, profile_id,
    methods = c("mean", "native_ldknn")) {
  profiles <- config$qc_profiles
  if (length(profile_id) != 1L || profile_id == "shared" ||
      !profile_id %in% names(profiles) ||
      any(!methods %in% c("mean", "TASSEL_LinkImpute", "native_ldknn")) || anyDuplicated(methods)) {
    stop("Choose one named QC profile and supported imputation methods.")
  }
  shared <- profiles$shared
  settings <- list(maf_min = profiles[[profile_id]]$maf_min,
    min_mac = shared$min_mac, marker_missing_max = shared$marker_missing_max,
    marker_heterozygosity_max = shared$marker_heterozygosity_max)
  if (!is.finite(settings$maf_min) || settings$maf_min < 0 || settings$maf_min >= 0.5) {
    stop("Invalid configured MAF cutoff.")
  }
  candidates <- list()
  if ("mean" %in% methods) candidates[[length(candidates) + 1L]] <- list(
    id = paste0(profile_id, "_mean"), method = "mean", settings = settings)
  for (method in intersect(c("native_ldknn", "TASSEL_LinkImpute"), methods)) {
    search <- if (method == "native_ldknn") config$native_ldknn else config$linkimpute
    grid <- expand.grid(high_LD_sites = unlist(search$high_LD_sites),
      nearest_neighbor_taxa = unlist(search$nearest_neighbor_taxa))
    if (any(!is.finite(as.matrix(grid))) ||
        any(grid$high_LD_sites < 2 | grid$high_LD_sites > 2000) ||
        any(grid$nearest_neighbor_taxa < 2 | grid$nearest_neighbor_taxa > 200) ||
        anyDuplicated(grid)) stop("Invalid LD/kNN search grid.")
    for (i in seq_len(nrow(grid))) {
      candidates[[length(candidates) + 1L]] <- list(id = sprintf("%s_%s_L%04d_K%03d",
        profile_id, method, grid$high_LD_sites[[i]], grid$nearest_neighbor_taxa[[i]]),
        method = method, settings = c(settings, list(
          high_LD_sites = grid$high_LD_sites[[i]],
          nearest_neighbor_taxa = grid$nearest_neighbor_taxa[[i]],
          maximum_physical_distance_bp = search$maximum_physical_distance_bp)))
    }
  }
  if (!length(candidates) || length(candidates) > config$resources$maximum_candidates_per_QC_profile) {
    stop("Configured search is empty or exceeds its candidate budget.")
  }
  candidates
}

gbs_calibrate_imputer <- function(plan, candidates, backends,
                                 minimum_coverage = 0.98,
                                 maximum_candidates = 32L,
                                 progress = NULL) {
  if (!inherits(plan, "gbs_imputation_mask_plan") ||
      !is.list(candidates) || !length(candidates) || length(candidates) > maximum_candidates ||
      !is.finite(minimum_coverage) || minimum_coverage <= 0 || minimum_coverage > 1) {
    stop("Invalid mask plan, candidate list or calibration limits.")
  }
  if (!is.null(progress) && !is.function(progress)) stop("Progress callback must be a function or NULL.")
  ids <- vapply(candidates, function(x) x$id, character(1L))
  methods <- vapply(candidates, function(x) x$method, character(1L))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Candidate IDs must be unique.")
  if (any(!methods %in% names(backends)) ||
      any(!vapply(backends[unique(methods)], is.function, logical(1L)))) {
    stop("A requested imputation backend is unavailable; no candidate was run.")
  }
  # Search LD/kNN within one QC profile. Compare MAF profiles separately.
  qc_keys <- c("maf_min", "min_mac", "marker_missing_max", "marker_heterozygosity_max")
  profile <- lapply(candidates, function(x) (x$settings %||% list())[qc_keys])
  if (any(!vapply(profile, identical, logical(1L), profile[[1L]]))) {
    stop("Calibrate each QC/MAF profile separately; filtering changes must remain explicit.")
  }
  tuning_truth <- plan$truth[plan$truth$phase == "tuning", , drop = FALSE]
  verification_truth <- plan$truth[plan$truth$phase == "verification", , drop = FALSE]
  if (!nrow(tuning_truth) || !nrow(verification_truth)) stop("Both mask phases are required.")
  best <- NULL
  rows <- vector("list", length(candidates))
  diagnostics <- vector("list", length(candidates))
  for (i in order(ids)) {
    candidate <- candidates[[i]]
    start <- proc.time()
    # The backend receives neither truth ledger nor unmasked input nor phenotype.
    prediction <- backends[[candidate$method]](plan$masked, candidate$settings %||% list())
    elapsed <- proc.time() - start
    gbs_validate_imputation_output(prediction, plan$masked)
    score <- gbs_score_imputation(prediction, tuning_truth)
    rows[[i]] <- data.frame(candidate_id = candidate$id, method = candidate$method,
      score$overall, elapsed_seconds = unname(elapsed[["elapsed"]]),
      CPU_seconds = unname(elapsed[["user.self"]] + elapsed[["sys.self"]]))
    diagnostics[[i]] <- list(candidate_id = candidate$id, by_maf = score$by_maf,
      by_mac = score$by_mac, by_genotype = score$by_genotype, by_marker = score$by_marker,
      backend_diagnostics = attr(prediction, "imputation_diagnostics"))
    eligible <- score$overall$coverage >= minimum_coverage &&
      is.finite(score$overall$dosage_rmse_called)
    if (eligible && (is.null(best) || score$overall$dosage_rmse_called < best$rmse)) {
      best <- list(id = candidate$id, method = candidate$method,
        settings = candidate$settings, rmse = score$overall$dosage_rmse_called,
        prediction = prediction)
    }
    if (!is.null(progress)) progress(candidate$id, rows[[i]])
  }
  tuning <- do.call(rbind, rows)
  if (is.null(best)) return(list(status = "no_candidate_met_coverage_gate", tuning = tuning,
    tuning_diagnostics = diagnostics, selected = NULL, verification = NULL))
  # Reveal verification scores only after the selected setting is fixed.
  verification <- gbs_score_imputation(best$prediction, verification_truth)
  verification_passed <- verification$overall$coverage >= minimum_coverage
  list(status = if (verification_passed) "selected_and_assessed_on_verification_masks" else
      "selected_but_verification_coverage_failed", tuning = tuning,
    tuning_diagnostics = diagnostics,
    selected = best[c("id", "method", "settings")], verification = verification,
    provenance = list(selection = "minimum_tuning_dosage_RMSE_after_coverage_gate",
      tie_break = "candidate_ID_lexical_order", verification_used_for_selection = FALSE,
      verification_coverage_passed = verification_passed,
      natural_missing_accuracy_established = FALSE,
      mask_sha256 = plan$provenance$mask_sha256))
}
