# Native R LD-aware nearest-neighbour imputation, version 1.
# This is an explicit implementation, not a bit-for-bit TASSEL reproduction.

gbs_reference_marker_statistics <- function(genotype) {
  gbs_validate_calls(genotype)
  n <- colSums(!is.na(genotype))
  allele_count <- colSums(genotype, na.rm = TRUE)
  frequency <- allele_count / (2 * n)
  maf <- pmin(frequency, 1 - frequency)
  data.frame(marker_id = colnames(genotype), n_observed = n,
    missingness = colMeans(is.na(genotype)), alt_frequency = frequency,
    maf = maf, mac = 2 * n * maf,
    heterozygosity = colSums(genotype == 1L, na.rm = TRUE) / n,
    stringsAsFactors = FALSE)
}

gbs_reference_qc <- function(genotype, settings = list()) {
  stats <- gbs_reference_marker_statistics(genotype)
  stats$fail_missingness <- stats$missingness >= (settings$marker_missing_max %||% 0.40)
  stats$fail_maf <- !is.finite(stats$maf) | stats$maf <= (settings$maf_min %||% 0.01)
  stats$fail_mac <- !is.finite(stats$mac) | stats$mac < (settings$min_mac %||% 0)
  stats$fail_heterozygosity <- !is.finite(stats$heterozygosity) |
    stats$heterozygosity >= (settings$marker_heterozygosity_max %||% 0.05)
  stats$retained <- !stats$fail_missingness & !stats$fail_maf &
    !stats$fail_mac & !stats$fail_heterozygosity
  stats
}

gbs_native_memory_estimate <- function(n_reference, n_markers, block_size, max_ld_sites) {
  (as.double(n_reference) * n_markers * 8 * 5 +
    as.double(block_size) * n_markers * 8 * 4 +
    as.double(n_markers) * max_ld_sites * 16) / 1024^3 + 0.10
}

fit_gbs_native_ldknn <- function(reference, markers, settings = list()) {
  if (!is.matrix(reference) || !is.numeric(reference) || !length(reference)) {
    stop("Reference must be a nonempty numeric genotype matrix.")
  }
  n <- nrow(reference)
  p <- ncol(reference)
  block_size <- as.integer(settings$block_size %||% 32L)
  max_ld <- as.integer(settings$max_ld_sites %||% 120L)
  memory <- settings$maximum_memory_gib %||% 2
  estimated <- gbs_native_memory_estimate(n, p, block_size, max_ld)
  if (!is.finite(estimated) || !is.finite(memory) || estimated > memory || memory <= 0) {
    stop("Native imputation fit exceeds its memory budget before allocation.")
  }
  gbs_validate_calls(reference)
  if (n < 3L || block_size < 1L || max_ld < 2L) stop("Insufficient reference lines or invalid LD limits.")
  required <- c("marker_id", "chromosome", "position_bp")
  if (!is.data.frame(markers) || !all(required %in% names(markers)) ||
      anyNA(markers$marker_id) || anyDuplicated(markers$marker_id)) stop("Invalid marker metadata.")
  marker_index <- match(colnames(reference), markers$marker_id)
  if (anyNA(marker_index)) stop("Reference markers are missing metadata.")
  markers <- markers[marker_index, , drop = FALSE]
  rownames(markers) <- NULL
  max_distance <- settings$max_distance_bp %||% 10000000
  min_pair <- as.integer(settings$min_pair_calls %||% 20L)
  min_r2 <- settings$min_ld_r2 %||% 0
  if (!is.finite(max_distance) || (max_distance < 0 && max_distance != -1) ||
      min_pair < 3L || !is.finite(min_r2) || min_r2 < 0 || min_r2 > 1) stop("Invalid LD settings.")
  if (max_distance >= 0 && (anyNA(markers$chromosome) ||
      any(!nzchar(markers$chromosome)) || anyNA(markers$position_bp) ||
      any(!is.finite(markers$position_bp)) || any(markers$position_bp <= 0))) {
    stop("A physical LD window requires explicit valid chromosome/position metadata.")
  }
  qc <- gbs_reference_qc(reference, settings)
  keep <- which(qc$retained)
  if (length(keep) < 3L) stop("Fewer than three reference-qualified markers remain.")
  z <- reference[, keep, drop = FALSE]
  map <- markers[keep, , drop = FALSE]
  p <- ncol(z)
  present <- !is.na(z)
  neighbors <- weights <- vector("list", p)
  names(neighbors) <- names(weights) <- colnames(z)
  index_order <- if (max_distance < 0) order(map$marker_id) else
    order(map$chromosome, map$position_bp, map$marker_id)
  started <- proc.time()
  for (start in seq.int(1L, p, by = block_size)) {
    targets <- index_order[seq.int(start, min(p, start + block_size - 1L))]
    allowed <- lapply(targets, function(j) {
      if (max_distance < 0) setdiff(seq_len(p), j) else
        which(map$chromosome == map$chromosome[[j]] &
          abs(map$position_bp - map$position_bp[[j]]) <= max_distance & seq_len(p) != j)
    })
    pool <- sort(unique(unlist(allowed, use.names = FALSE)))
    if (!length(pool)) next
    pair_n <- crossprod(present[, targets, drop = FALSE], present[, pool, drop = FALSE])
    correlations <- suppressWarnings(stats::cor(z[, targets, drop = FALSE],
      z[, pool, drop = FALSE], use = "pairwise.complete.obs"))
    for (k in seq_along(targets)) {
      j <- targets[[k]]
      local <- match(allowed[[k]], pool)
      r2 <- correlations[k, local]^2
      valid <- is.finite(r2) & pair_n[k, local] >= min_pair & r2 > min_r2
      candidates <- allowed[[k]][valid]
      r2 <- r2[valid]
      if (!length(candidates)) next
      rank <- order(-r2, map$marker_id[candidates])
      rank <- head(rank, max_ld)
      neighbors[[j]] <- candidates[rank]
      weights[[j]] <- r2[rank]
    }
  }
  elapsed <- proc.time() - started
  structure(list(reference = z, marker_ids = colnames(z), markers = map,
    marker_statistics = qc, settings = settings, neighbors = neighbors,
    ld_weights = weights, max_ld_sites = max_ld,
    training_means = colMeans(z, na.rm = TRUE),
    training_modes = vapply(seq_len(ncol(z)), function(j) {
      counts <- tabulate(z[, j] + 1L, nbins = 3L)
      which.max(counts) - 1L
    }, numeric(1L)),
    provenance = list(algorithm = "native_R_LD_weighted_kNN_v1",
      source_software = "project_native_not_TASSEL",
      reference_ids = rownames(z), reference_sha256 = digest::digest(reference, algo = "sha256"),
      metadata_sha256 = digest::digest(markers, algo = "sha256"),
      settings_sha256 = digest::digest(settings, algo = "sha256"),
      phenotype_used = FALSE, query_used_for_fit = FALSE,
      estimated_fit_memory_gib = estimated,
      LD_fit_elapsed_seconds = unname(elapsed[["elapsed"]]),
      LD_fit_CPU_seconds = unname(elapsed[["user.self"]] + elapsed[["sys.self"]]),
      R_version = as.character(getRversion()))), class = "gbs_native_ldknn")
}

predict_gbs_native_ldknn <- function(model, query, high_ld_sites = 60L,
    knn_taxa = 20L, output = c("dosage", "hardcall"),
    minimum_shared_sites = 3L, minimum_shared_fraction = 0.25,
    minimum_neighbors = 2L, fallback = c("reference", "abstain"),
    keep_cell_diagnostics = FALSE) {
  output <- match.arg(output)
  fallback <- match.arg(fallback)
  if (!inherits(model, "gbs_native_ldknn")) stop("Invalid native imputation model.")
  if (as.double(length(query)) * 8 * 4 / 1024^3 +
      model$provenance$estimated_fit_memory_gib >
      (model$settings$maximum_memory_gib %||% 2)) stop("Reference plus query exceeds memory budget.")
  gbs_validate_calls(query)
  if (anyNA(match(model$marker_ids, colnames(query)))) stop("Query is missing fitted marker IDs.")
  if (!is.finite(high_ld_sites) || high_ld_sites < 2L || high_ld_sites > model$max_ld_sites ||
      !is.finite(knn_taxa) || knn_taxa < 1L || minimum_shared_sites < 1L ||
      minimum_shared_fraction <= 0 || minimum_shared_fraction > 1 || minimum_neighbors < 1L) {
    stop("Invalid prediction/neighbor settings.")
  }
  if (as.double(length(query)) * 8 * 4 / 1024^3 + 0.10 >
      (model$settings$maximum_memory_gib %||% 2)) stop("Query output exceeds memory budget.")
  result <- query
  counts <- c(ld_knn = 0L, mean_fallback = 0L, mode_fallback = 0L, abstained = 0L)
  cells <- list()
  reference <- model$reference
  for (j in seq_along(model$marker_ids)) {
    qcol <- match(model$marker_ids[[j]], colnames(query))
    missing_rows <- which(is.na(query[, qcol]))
    if (!length(missing_rows)) next
    predictors <- head(model$neighbors[[j]], high_ld_sites)
    ld_weights <- head(model$ld_weights[[j]], high_ld_sites)
    donors <- which(!is.na(reference[, j]))
    required <- max(minimum_shared_sites, ceiling(length(predictors) * minimum_shared_fraction))
    for (r in missing_rows) {
      value <- NA_real_
      reason <- "abstained"
      used_neighbors <- 0L
      if (length(predictors) >= required && length(donors) >= minimum_neighbors) {
        query_values <- query[r, match(model$marker_ids[predictors], colnames(query))]
        donor_values <- reference[donors, predictors, drop = FALSE]
        shared <- !is.na(donor_values) &
          matrix(!is.na(query_values), nrow(donor_values), length(query_values), byrow = TRUE)
        distance <- abs(sweep(donor_values, 2L, query_values, "-")) / 2
        distance[!shared] <- 0
        numerator <- rowSums(sweep(distance, 2L, ld_weights, "*"))
        denominator <- as.numeric(shared %*% ld_weights)
        eligible <- rowSums(shared) >= required & denominator > 0
        eligible <- eligible & rownames(reference)[donors] != rownames(query)[[r]]
        positions <- which(eligible)
        if (length(positions) >= minimum_neighbors) {
          d <- numerator[positions] / denominator[positions]
          rank <- order(d, rownames(reference)[donors[positions]])
          rank <- head(rank, knn_taxa)
          selected <- positions[rank]
          d <- d[rank]
          exact <- d < 1e-12
          if (any(exact)) { selected <- selected[exact]; d <- d[exact] }
          if (length(selected) >= minimum_neighbors) {
            w <- if (all(d < 1e-12)) rep(1, length(d)) else 1 / pmax(d, 1e-8)
            g <- reference[donors[selected], j]
            if (output == "dosage") value <- sum(w * g) / sum(w) else {
              votes <- vapply(0:2, function(call) sum(w[g == call]), numeric(1L))
              winners <- which(abs(votes - max(votes)) < 1e-10 * max(1, max(votes)))
              if (length(winners) == 1L) value <- winners - 1L
            }
            if (is.finite(value)) { reason <- "ld_knn"; used_neighbors <- length(selected) }
          }
        }
      }
      if (!is.finite(value) && fallback == "reference") {
        value <- if (output == "dosage") model$training_means[[j]] else model$training_modes[[j]]
        reason <- if (output == "dosage") "mean_fallback" else "mode_fallback"
      }
      result[r, qcol] <- value
      counts[[reason]] <- counts[[reason]] + 1L
      if (keep_cell_diagnostics) cells[[length(cells) + 1L]] <- data.frame(
        sample_id = rownames(query)[[r]], marker_id = model$marker_ids[[j]],
        method = reason, n_predictor_sites = length(predictors),
        n_neighbors_used = used_neighbors, dosage = value)
    }
  }
  gbs_validate_imputation_output(result, query)
  attr(result, "retained_marker_ids") <- model$marker_ids
  attr(result, "imputation_diagnostics") <- list(counts = counts,
    cells = if (length(cells)) do.call(rbind, cells) else NULL,
    output = output, fallback = fallback, reference_sha256 = model$provenance$reference_sha256,
    query_used_for_fit = FALSE, independent_queries = TRUE)
  result
}

calibrate_gbs_native_imputation <- function(reference, query, markers, config,
    profile_id = "low_frequency_candidate", native_settings = list(), progress = NULL) {
  gbs_validate_calls(reference)
  gbs_validate_calls(query)
  if (length(intersect(rownames(reference), rownames(query)))) {
    stop("Imputation calibration reference and query IDs must be disjoint.")
  }
  candidates <- gbs_build_imputation_candidates(config, profile_id, c("mean", "native_ldknn"))
  search <- config$native_ldknn
  fit_settings <- utils::modifyList(list(max_ld_sites = max(unlist(search$high_LD_sites)),
    max_distance_bp = search$maximum_physical_distance_bp,
    min_pair_calls = search$min_pair_calls, min_ld_r2 = search$min_ld_r2,
    block_size = search$block_size, maximum_memory_gib = config$resources$maximum_memory_gib),
    native_settings)
  fit_settings <- utils::modifyList(fit_settings, candidates[[1L]]$settings)
  model <- fit_gbs_native_ldknn(reference, markers, fit_settings)
  missing <- setdiff(model$marker_ids, colnames(query))
  if (length(missing)) stop("Query lacks reference-qualified calibration markers.")
  query <- query[, model$marker_ids, drop = FALSE]
  plan <- gbs_make_mask_plan(query,
    fraction_per_phase = config$masking$fraction_per_phase,
    seed = config$masking$seed,
    minimum_remaining_calls = config$masking$minimum_remaining_calls_per_marker,
    maximum_masks_per_phase = config$masking$maximum_masks_per_phase,
    maximum_memory_gib = config$resources$maximum_memory_gib,
    reference_statistics = model$marker_statistics)
  backends <- list(mean = function(masked, settings)
      predict_gbs_reference_mean(model$reference, masked, settings),
    native_ldknn = function(masked, settings) predict_gbs_native_ldknn(model, masked,
      high_ld_sites = settings$high_LD_sites, knn_taxa = settings$nearest_neighbor_taxa,
      output = search$output, minimum_shared_sites = search$minimum_shared_sites,
      minimum_shared_fraction = search$minimum_shared_fraction,
      minimum_neighbors = search$minimum_neighbors, fallback = search$fallback))
  calibration <- gbs_calibrate_imputer(plan, candidates, backends,
    minimum_coverage = config$selection$minimum_masked_call_coverage,
    maximum_candidates = config$resources$maximum_candidates_per_QC_profile,
    progress = progress)
  list(model = model, mask_plan = plan, calibration = calibration,
    shared_LD_fit_elapsed_seconds = model$provenance$LD_fit_elapsed_seconds,
    provenance = list(reference_query_disjoint = TRUE, phenotype_used = FALSE,
      reference_sha256 = digest::digest(reference, algo = "sha256"),
      query_sha256 = digest::digest(query, algo = "sha256"),
      profile_id = profile_id, config_sha256 = digest::digest(config, algo = "sha256")))
}

predict_gbs_reference_mean <- function(reference, query, settings = list()) {
  gbs_validate_calls(reference)
  gbs_validate_calls(query)
  qc <- gbs_reference_qc(reference, settings)
  ids <- qc$marker_id[qc$retained]
  if (anyNA(match(ids, colnames(query)))) stop("Query lacks reference-qualified markers.")
  result <- query
  means <- colMeans(reference[, ids, drop = FALSE], na.rm = TRUE)
  for (id in ids) result[is.na(result[, id]), id] <- means[[id]]
  attr(result, "retained_marker_ids") <- ids
  result
}
