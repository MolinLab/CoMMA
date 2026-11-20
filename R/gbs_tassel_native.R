# SPDX-License-Identifier: LGPL-2.1-only
# Native biallelic adaptation of TASSEL 5.2.96 LDKNNiImputationPlugin.
# Algorithm: Daniel Money; TASSEL conversion: Ed Buckler and contributors.
# R implementation and reference/query extension: this project, 2026.
# See THIRD_PARTY_NOTICES.md and docs/tassel_native_algorithm_audit.md.
# No Java is launched. Numerical equivalence to TASSEL is not yet established.

gbs_tassel_float <- function(x) {
  readBin(writeBin(as.numeric(x), raw(), size = 4L), "numeric", n = length(x), size = 4L)
}

gbs_tassel_byte_codes <- function(alleles) {
  codes <- c(A = 0L, C = 1L, G = 2L, T = 3L)
  pairs <- strsplit(toupper(alleles), "/", fixed = TRUE)
  if (any(lengths(pairs) != 2L) || any(vapply(pairs, function(z)
      any(!z %in% names(codes)) || z[[1L]] == z[[2L]], logical(1L)))) {
    stop("TASSEL native mode requires explicit distinct biallelic A/C/G/T metadata.")
  }
  result <- t(vapply(pairs, function(z) {
    a <- unname(codes[z[[1L]]]); b <- unname(codes[z[[2L]]])
    c(17L * a, 16L * min(a, b) + max(a, b), 17L * b)
  }, numeric(3L)))
  colnames(result) <- c("0", "1", "2")
  result
}

gbs_tassel_r2_counts <- function(a, b, c, d, min_count = 20L, min_minor = 2L) {
  n <- a + b + c + d
  p <- (a + b) / n
  q <- (a + c) / n
  value <- ((a / n) * (d / n) - (b / n) * (c / n))^2 /
    (p * (1 - p) * q * (1 - q))
  value[n < min_count | c + d < min_minor | b + d < min_minor |
    !is.finite(value)] <- NA_real_
  # TASSEL LDResult stores r2 as float, which affects predictor ranking ties.
  matrix(gbs_tassel_float(value), nrow = nrow(a), ncol = ncol(a))
}

gbs_tassel_distance <- function(query, donors) {
  if (!is.matrix(donors) || ncol(donors) != length(query)) stop("Invalid distance dimensions.")
  shared <- !is.na(donors) & matrix(!is.na(query), nrow(donors), ncol(donors), byrow = TRUE)
  count <- rowSums(shared)
  differences <- abs(sweep(donors, 2L, query, "-"))
  differences[!shared] <- 0
  distance <- rowSums(differences) / (2 * count)
  distance[count < 10L] <- NA_real_
  list(distance = distance, shared_calls = count)
}

gbs_tassel_vote <- function(genotypes, distances, requested_ld_sites, byte_codes) {
  if (!length(genotypes)) return(NA_real_)
  if (length(genotypes) != length(distances) || anyNA(genotypes) ||
      any(!genotypes %in% 0:2) || any(!is.finite(distances)) ||
      any(distances < 0 | distances > 1)) stop("Invalid TASSEL neighbour vote.")
  weights <- 1 / (1 + requested_ld_sites * distances)
  votes <- vapply(0:2, function(g) sum(weights[genotypes == g]), numeric(1L))
  winners <- which(votes == max(votes))
  as.numeric(winners[which.min(byte_codes[winners])] - 1L)
}

fit_gbs_tassel_native <- function(reference, markers, settings = list()) {
  n <- nrow(reference); p <- ncol(reference)
  maximum_memory <- settings$maximum_memory_gib %||% 2
  max_ld <- settings$max_ld_sites %||% 60L
  block_size <- settings$block_size %||% 16L
  max_distance <- settings$max_distance_bp %||% 10000000
  chromosome_policy <- settings$chromosome_policy %||% "tassel_position_only"
  integer_setting <- function(x, low, high) length(x) == 1L && !is.na(x) &&
    is.finite(x) && x == floor(x) && x >= low && x <= high
  if (!integer_setting(max_ld, 2L, 2000L) || !integer_setting(block_size, 1L, 1024L) ||
      length(max_distance) != 1L || !is.finite(max_distance) ||
      (max_distance < 0 && max_distance != -1) ||
      length(chromosome_policy) != 1L ||
      !chromosome_policy %in% c("tassel_position_only", "same_chromosome")) stop("Invalid TASSEL fit settings.")
  estimated <- as.double(length(reference)) * 8 * 12 / 1024^3 +
    as.double(block_size) * p * 8 * 10 / 1024^3 + as.double(p) * max_ld * 16 / 1024^3 + 0.10
  if (length(maximum_memory) != 1L || !is.finite(maximum_memory) || maximum_memory <= 0 ||
      !is.finite(estimated) || estimated > maximum_memory) stop("TASSEL native fit exceeds memory guard.")
  gbs_validate_calls(reference)
  if (!is.data.frame(markers) || !all(c("marker_id", "alleles", "chromosome", "position_bp") %in% names(markers)) ||
      anyNA(markers$marker_id) || anyDuplicated(markers$marker_id)) stop("Invalid TASSEL marker metadata.")
  indices <- match(colnames(reference), markers$marker_id)
  if (anyNA(indices)) stop("Missing reference marker metadata.")
  map <- markers[indices, , drop = FALSE]
  rownames(map) <- NULL
  if (max_distance >= 0 && (anyNA(map$position_bp) || any(!is.finite(map$position_bp)) ||
      any(map$position_bp <= 0) || anyNA(map$chromosome) || any(!nzchar(map$chromosome)))) stop("Invalid physical map.")
  byte_codes <- gbs_tassel_byte_codes(map$alleles)
  rownames(byte_codes) <- map$marker_id
  statistics <- gbs_reference_marker_statistics(reference)
  frequency <- statistics$alt_frequency
  # Equal-frequency alleles follow nucleotide order: the smaller code is major.
  alt_is_minor <- frequency < 0.5 | (frequency == 0.5 & byte_codes[, 1L] < byte_codes[, 3L])
  alt_is_minor[!is.finite(frequency)] <- TRUE
  major <- reference < 2L
  minor <- reference > 0L
  flip <- which(!alt_is_minor)
  if (length(flip)) {
    major[, flip] <- reference[, flip, drop = FALSE] > 0L
    minor[, flip] <- reference[, flip, drop = FALSE] < 2L
  }
  major[is.na(major)] <- FALSE
  minor[is.na(minor)] <- FALSE
  neighbours <- r2_values <- vector("list", p)
  cutoff_tie <- logical(p)
  names(neighbours) <- names(r2_values) <- map$marker_id
  started <- proc.time()
  for (start in seq.int(1L, p, by = block_size)) {
    targets <- seq.int(start, min(p, start + block_size - 1L))
    allowed <- lapply(targets, function(j) {
      ok <- seq_len(p) != j
      if (max_distance >= 0) ok <- ok & abs(map$position_bp - map$position_bp[[j]]) <= max_distance
      if (chromosome_policy == "same_chromosome") ok <- ok & map$chromosome == map$chromosome[[j]]
      which(ok)
    })
    pool <- sort(unique(unlist(allowed, use.names = FALSE)))
    if (!length(pool)) next
    a <- crossprod(major[, targets, drop = FALSE], major[, pool, drop = FALSE])
    b <- crossprod(minor[, targets, drop = FALSE], major[, pool, drop = FALSE])
    c <- crossprod(major[, targets, drop = FALSE], minor[, pool, drop = FALSE])
    d <- crossprod(minor[, targets, drop = FALSE], minor[, pool, drop = FALSE])
    r2 <- gbs_tassel_r2_counts(a, b, c, d)
    for (k in seq_along(targets)) {
      j <- targets[[k]]
      values <- r2[k, match(allowed[[k]], pool)]
      valid <- is.finite(values)
      candidates <- allowed[[k]][valid]
      values <- values[valid]
      rank <- order(-values, map$marker_id[candidates])
      if (length(rank) > max_ld) cutoff_tie[[j]] <- values[rank[[max_ld]]] == values[rank[[max_ld + 1L]]]
      rank <- head(rank, max_ld)
      neighbours[[j]] <- candidates[rank]
      r2_values[[j]] <- values[rank]
    }
  }
  structure(list(reference = reference, markers = map, marker_ids = map$marker_id,
    marker_statistics = statistics, neighbours = neighbours, r2 = r2_values,
    byte_codes = byte_codes, LD_cutoff_tie = cutoff_tie,
    max_ld_sites = as.integer(max_ld), settings = list(max_distance_bp = max_distance,
      chromosome_policy = chromosome_policy, maximum_memory_gib = maximum_memory,
      block_size = block_size),
    provenance = list(algorithm = "native_R_TASSEL_5.2.96_biallelic_v1",
      numerical_TASSEL_parity_verified = FALSE, input_QC_applied_internally = FALSE,
      LD_boundary_tie_policy = "marker_ID_lexical_order_not_Guava_queue",
      fitted_reference_sha256 = digest::digest(reference, algo = "sha256"),
      marker_metadata_sha256 = digest::digest(map, algo = "sha256"),
      phenotype_used = FALSE, query_used_for_fit = FALSE,
      estimated_fit_memory_gib = estimated,
      elapsed_fit_seconds = unname((proc.time() - started)[["elapsed"]]),
      R_version = as.character(getRversion()))), class = "gbs_tassel_native")
}

predict_gbs_tassel_native <- function(model, query, high_ld_sites = 60L, knn_taxa = 20L,
                                       keep_cell_diagnostics = FALSE) {
  if (!inherits(model, "gbs_tassel_native")) stop("Invalid TASSEL native model.")
  estimated <- model$provenance$estimated_fit_memory_gib + as.double(length(query)) * 8 * 5 / 1024^3
  if (!is.finite(estimated) || estimated > model$settings$maximum_memory_gib) stop("TASSEL native prediction exceeds memory guard.")
  gbs_validate_calls(query)
  valid_integer <- function(x, low, high) length(x) == 1L && !is.na(x) && is.finite(x) &&
    x == floor(x) && x >= low && x <= high
  if (!valid_integer(high_ld_sites, 2L, model$max_ld_sites) || !valid_integer(knn_taxa, 2L, 200L)) stop("Invalid LD or k setting.")
  if (anyNA(match(model$marker_ids, colnames(query)))) stop("Query is missing fitted marker IDs.")
  query_aligned <- query[, model$marker_ids, drop = FALSE]
  shared_ids <- intersect(rownames(model$reference), rownames(query))
  if (length(shared_ids)) {
    old <- unname(model$reference[shared_ids, , drop = FALSE])
    new <- unname(query_aligned[shared_ids, , drop = FALSE])
    if (!identical(is.na(old), is.na(new)) || any(old[!is.na(old)] != new[!is.na(new)])) {
      stop("Overlapping reference/query IDs must have identical raw masked calls; otherwise use disjoint IDs.")
    }
  }
  result <- query
  counts <- c(ld_knn = 0L, monomorphic = 0L, abstained = 0L)
  details <- list()
  for (j in seq_along(model$marker_ids)) {
    qcol <- match(model$marker_ids[[j]], colnames(query))
    missing_rows <- which(is.na(query[, qcol]))
    if (!length(missing_rows)) next
    predictors <- head(model$neighbours[[j]], high_ld_sites)
    cutoff_tie <- if (length(model$r2[[j]]) > high_ld_sites)
      model$r2[[j]][[high_ld_sites]] == model$r2[[j]][[high_ld_sites + 1L]] else
      high_ld_sites == model$max_ld_sites && model$LD_cutoff_tie[[j]]
    # TASSEL genotype-table filters expose positions in genomic order.
    predictors <- predictors[order(model$markers$chromosome[predictors],
      model$markers$position_bp[predictors], model$marker_ids[predictors])]
    frequency <- model$marker_statistics$alt_frequency[[j]]
    monomorphic <- is.finite(frequency) && frequency %in% c(0, 1)
    donor_values <- model$reference[, predictors, drop = FALSE]
    donor_coverage <- if (length(predictors)) rowMeans(!is.na(donor_values)) else rep(0, nrow(donor_values))
    for (r in missing_rows) {
      value <- NA_real_; reason <- "abstained"; n_donors <- 0L
      if (monomorphic) {
        value <- 2 * frequency; reason <- "monomorphic"
      } else if (length(predictors)) {
        qvalues <- query_aligned[r, predictors]
        expected_overlap <- donor_coverage * mean(!is.na(qvalues)) * length(predictors)
        eligible <- which(expected_overlap > 10 & !is.na(model$reference[, j]) &
          rownames(model$reference) != rownames(query)[[r]])
        if (length(eligible)) {
          distances <- gbs_tassel_distance(qvalues, donor_values[eligible, , drop = FALSE])$distance
          good <- is.finite(distances)
          eligible <- eligible[good]; distances <- distances[good]
          if (length(eligible)) {
            calls <- model$reference[eligible, j]
            rank <- head(order(distances, model$byte_codes[j, calls + 1L],
              rownames(model$reference)[eligible]), knn_taxa)
            value <- gbs_tassel_vote(calls[rank], distances[rank], high_ld_sites, model$byte_codes[j, ])
            n_donors <- length(rank); reason <- "ld_knn"
          }
        }
      }
      result[r, qcol] <- value
      counts[[reason]] <- counts[[reason]] + 1L
      if (keep_cell_diagnostics) details[[length(details) + 1L]] <- data.frame(
        sample_id = rownames(query)[[r]], marker_id = model$marker_ids[[j]],
        method = reason, n_predictors = length(predictors), n_donors = n_donors,
        LD_cutoff_tie = cutoff_tie, hardcall = value)
    }
  }
  gbs_validate_imputation_output(result, query)
  attr(result, "imputation_diagnostics") <- list(counts = counts,
    cells = if (length(details)) do.call(rbind, details) else NULL,
    hardcalls_only = TRUE, reference_fallback_used = FALSE, independent_queries = TRUE,
    query_used_for_fit = FALSE, numerical_TASSEL_parity_verified = FALSE)
  result
}

impute_gbs_tassel_joint <- function(masked, markers, high_ld_sites = 60L, knn_taxa = 20L,
                                    settings = list()) {
  settings$max_ld_sites <- high_ld_sites
  model <- fit_gbs_tassel_native(masked, markers, settings)
  list(model = model, prediction = predict_gbs_tassel_native(model, masked,
    high_ld_sites, knn_taxa, keep_cell_diagnostics = TRUE))
}

calibrate_gbs_tassel_native <- function(model, query, calibration, progress = NULL) {
  if (!inherits(model, "gbs_tassel_native")) stop("Invalid fitted reference model.")
  gbs_validate_calls(query)
  if (length(intersect(rownames(model$reference), rownames(query)))) stop("Calibration requires disjoint reference/query IDs.")
  if (anyNA(match(model$marker_ids, colnames(query)))) stop("Calibration query lacks fitted markers.")
  grid <- expand.grid(ld = calibration$high_LD_sites, k = calibration$nearest_neighbor_taxa)
  if (!nrow(grid) || nrow(grid) > 32L || any(!is.finite(as.matrix(grid))) ||
      any(grid$ld != floor(grid$ld) | grid$ld < 2 | grid$ld > model$max_ld_sites) ||
      any(grid$k != floor(grid$k) | grid$k < 2 | grid$k > 200) || anyDuplicated(grid)) stop("Invalid bounded calibration grid.")
  candidates <- lapply(seq_len(nrow(grid)), function(i) list(
    id = sprintf("tassel_native_LD%04d_k%03d", grid$ld[[i]], grid$k[[i]]),
    method = "tassel_native", profile_id = "fixed_reference_QC",
    settings = list(high_LD_sites = grid$ld[[i]], nearest_neighbor_taxa = grid$k[[i]])))
  query <- query[, model$marker_ids, drop = FALSE]
  stats <- model$marker_statistics
  maskable <- is.finite(stats$alt_frequency) & stats$alt_frequency > 0 & stats$alt_frequency < 1 & stats$n_observed > 0
  if (!any(maskable)) stop("No reference-polymorphic markers are available for calibration.")
  plan <- gbs_make_mask_plan(query[, maskable, drop = FALSE],
    fraction_per_phase = calibration$fraction_per_phase,
    seed = calibration$seed, minimum_remaining_calls = 5L,
    maximum_masks_per_phase = calibration$maximum_masks_per_phase,
    maximum_memory_gib = model$settings$maximum_memory_gib,
    reference_statistics = stats[maskable, , drop = FALSE])
  # Keep monomorphic/unobserved reference sites in backend inputs, without
  # inventing reference frequencies or masking undefined truth strata.
  masked <- query
  masked[cbind(match(plan$truth$sample_id, rownames(masked)),
    match(plan$truth$marker_id, colnames(masked)))] <- NA_real_
  plan$masked <- masked
  plan$provenance$input_sha256 <- digest::digest(query, algo = "sha256")
  plan$provenance$markers_excluded_from_masks <- stats$marker_id[!maskable]
  result <- gbs_calibrate_imputer(plan, candidates,
    list(tassel_native = function(masked, settings) predict_gbs_tassel_native(model, masked,
      high_ld_sites = settings$high_LD_sites, knn_taxa = settings$nearest_neighbor_taxa)),
    minimum_coverage = calibration$minimum_coverage, maximum_candidates = 32L, progress = progress)
  baseline <- plan$masked
  for (j in seq_along(model$marker_ids)) {
    baseline[is.na(baseline[, j]), j] <- mean(model$reference[, j], na.rm = TRUE)
  }
  list(calibration = result, mask_plan = plan,
    mean_baseline_tuning = gbs_score_imputation(baseline, plan$truth[plan$truth$phase == "tuning", ]),
    mean_baseline_verification = gbs_score_imputation(baseline, plan$truth[plan$truth$phase == "verification", ]),
    reference_query_disjoint = TRUE, phenotype_used = FALSE)
}
