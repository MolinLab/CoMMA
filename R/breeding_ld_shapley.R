# Estimate LD from original calls on the saved model's discovery lines.
breeding_ld_options <- function(options) {
  defaults <- list(scope = "chromosome", decay_r2 = .20, group_r2 = .20,
    max_distance_bp = 50000000, fallback_distance_bp = 250000,
    min_joint_calls = 20L, pairs_per_bin = 128L, min_pairs_per_bin = 30L)
  if (!is.list(options) || anyDuplicated(names(options)) ||
      length(options) != length(names(options)) || any(!names(options) %in% names(defaults))) {
    stop("Use named LD settings from the explanation configuration.")
  }
  x <- utils::modifyList(defaults, options, keep.null = TRUE)
  scalar <- function(value) is.numeric(value) && length(value) == 1L && is.finite(value)
  thresholds <- x[c("decay_r2", "group_r2")]
  counts <- x[c("max_distance_bp", "fallback_distance_bp", "min_joint_calls", "pairs_per_bin", "min_pairs_per_bin")]
  if (!identical(x$scope, "chromosome") && !identical(x$scope, "panel")) stop("Set LD scope to chromosome or panel.")
  if (!all(vapply(thresholds, function(v) scalar(v) && v > 0 && v < 1, logical(1L))) ||
      !all(vapply(counts, function(v) scalar(v) && v >= 1 && v == floor(v) && v <= .Machine$integer.max, logical(1L))) ||
      x$min_joint_calls < 3 || x$min_pairs_per_bin > x$pairs_per_bin || x$fallback_distance_bp > x$max_distance_bp) {
    stop("Invalid LD thresholds, distance bounds or sample counts.")
  }
  x
}

breeding_ld_pair <- function(a, b, minimum) {
  called <- is.finite(a) & is.finite(b)
  n <- sum(called)
  value <- if (n >= minimum) suppressWarnings(stats::cor(a[called], b[called])^2) else NA_real_
  c(r2 = value, joint_calls = n)
}

breeding_ld_bin <- function(raw, positions, indices, low, high, options) {
  # Count eligible pairs without allocating a pairwise matrix. A sampled rank
  # identifies a unique pair in the concatenated, sorted neighbor ranges.
  lower <- findInterval(positions + low, positions)
  upper <- findInterval(positions + high, positions)
  cumulative <- c(0, cumsum(as.double(upper - lower)))
  total <- tail(cumulative, 1L)
  n <- min(total, options$pairs_per_bin)
  row <- data.frame(bin_from_bp = low, bin_to_bp = high, possible_pairs = total,
    sampled_pairs = n, valid_pairs = 0, mean_distance_bp = NA_real_, mean_r2 = NA_real_)
  if (!n) return(row)
  sampled <- if (n == total) seq_len(n) else sample.int(total, n)
  first <- findInterval(sampled - 1, cumulative)
  second <- lower[first] + sampled - cumulative[first]
  r2 <- vapply(seq_along(first), function(i) {
    breeding_ld_pair(raw[, indices[first[i]]], raw[, indices[second[i]]], options$min_joint_calls)[[1L]]
  }, numeric(1L))
  keep <- is.finite(r2)
  row$valid_pairs <- sum(keep)
  if (any(keep)) {
    row$mean_distance_bp <- mean((positions[second] - positions[first])[keep])
    row$mean_r2 <- mean(r2[keep])
  }
  row
}

breeding_ld_crossing <- function(bins, options) {
  supported <- which(bins$valid_pairs >= options$min_pairs_per_bin)
  bins$smoothed_r2 <- NA_real_
  result <- list(distance = NA_real_, status = "insufficient_supported_bins")
  if (length(supported) >= 3L) {
    distance <- bins$mean_distance_bp[supported]
    smooth <- -stats::isoreg(distance, -bins$mean_r2[supported])$yf
    bins$smoothed_r2[supported] <- smooth
    crossing <- which(smooth <= options$decay_r2)
    if (!length(crossing)) {
      result <- list(distance = tail(distance, 1L), status = "above_threshold_at_last_supported_bin")
    } else if (crossing[1L] == 1L) {
      result <- list(distance = distance[1L], status = "below_threshold_at_first_supported_bin")
    } else {
      i <- crossing[1L]
      fraction <- (smooth[i - 1L] - options$decay_r2) / (smooth[i - 1L] - smooth[i])
      result <- list(distance = distance[i - 1L] + fraction * (distance[i] - distance[i - 1L]), status = "crossed_threshold")
    }
  }
  c(result, list(bins = bins, supported_bins = length(supported)))
}

breeding_ld_distances <- function(raw, map, ordered, options, seed) {
  grid <- 250000 * 2^(0:16)
  edges <- c(0, grid[grid < options$max_distance_bp], options$max_distance_bp)
  bins <- with_seed_preserved(seed, lapply(ordered, function(j) {
    do.call(rbind, lapply(seq_len(length(edges) - 1L), function(k) {
      breeding_ld_bin(raw, map$position_bp[j], j, edges[k], edges[k + 1L], options)
    }))
  }))
  pooled <- do.call(rbind, lapply(seq_len(length(edges) - 1L), function(k) {
    rows <- do.call(rbind, lapply(bins, function(x) x[k, ]))
    usable <- rows$valid_pairs > 0L
    weighted <- function(name) if (any(usable)) stats::weighted.mean(rows[[name]][usable], rows$valid_pairs[usable]) else NA_real_
    data.frame(bin_from_bp = edges[k], bin_to_bp = edges[k + 1L],
      possible_pairs = sum(rows$possible_pairs), sampled_pairs = sum(rows$sampled_pairs), valid_pairs = sum(rows$valid_pairs),
      mean_distance_bp = weighted("mean_distance_bp"), mean_r2 = weighted("mean_r2"))
  }))
  chromosomes <- lapply(bins, breeding_ld_crossing, options = options)
  panel <- breeding_ld_crossing(pooled, options)
  result <- lapply(names(chromosomes), function(chr) {
    local <- chromosomes[[chr]]
    use_local <- options$scope == "chromosome" && is.finite(local$distance)
    chosen <- if (use_local) local else panel
    fallback <- !is.finite(chosen$distance)
    data.frame(chromosome = chr, distance_bp = if (fallback) options$fallback_distance_bp else max(1, round(chosen$distance)),
      distance_source = if (fallback) "fallback" else if (use_local) "chromosome" else "panel",
      status = chosen$status, chromosome_status = local$status,
      chromosome_supported_bins = local$supported_bins, panel_supported_bins = panel$supported_bins)
  })
  by_chr <- Map(function(name, x) data.frame(scope = "chromosome", chromosome = name, x$bins), names(chromosomes), chromosomes)
  list(distances = do.call(rbind, result), bins = do.call(rbind, c(by_chr,
    list(data.frame(scope = "panel", chromosome = NA_character_, panel$bins)))))
}

breeding_ld_shap_partition <- function(prep, options, seed) {
  reference <- prep$imputer$reference
  markers <- prep$markers
  ids <- prep$analysis_ids
  if (!is.matrix(reference) || !identical(rownames(reference), prep$training_ids) ||
      !identical(colnames(prep$training$z), ids) || !all(ids %in% colnames(reference)) ||
      !identical(as.character(markers$marker_id), ids)) stop("LD grouping needs the saved discovery calls and marker map.")
  raw <- reference[, ids, drop = FALSE]
  map <- breeding_coordinate_map(markers)
  ordered <- split(seq_len(nrow(map)), map$chromosome)
  ordered <- lapply(ordered, function(j) j[order(map$position_bp[j], map$marker_id[j], method = "radix")])
  estimated <- breeding_ld_distances(raw, map, ordered, options, seed)
  assigned <- character(length(ids))
  anchor <- character(length(ids))
  distance <- r2 <- joint <- numeric(length(ids))
  for (chr in names(ordered)) {
    j <- ordered[[chr]]
    p <- map$position_bp[j]
    radius <- estimated$distances$distance_bp[match(chr, estimated$distances$chromosome)]
    ends <- findInterval(p + radius, p)
    for (i in seq_along(j)) {
      lead <- j[i]
      if (nzchar(assigned[lead])) next
      nearby <- if (ends[i] > i) j[seq.int(i + 1L, ends[i])] else integer()
      nearby <- nearby[!nzchar(assigned[nearby])]
      linkage <- vapply(nearby, function(k) breeding_ld_pair(raw[, lead], raw[, k], options$min_joint_calls),
                        c(r2 = 0, joint_calls = 0))
      keep <- which(is.finite(linkage["r2", ]) & linkage["r2", ] >= options$group_r2)
      members <- c(lead, nearby[keep])
      assigned[members] <- paste0("LD:", ids[lead])
      anchor[members] <- ids[lead]
      distance[members] <- radius
      r2[members] <- c(1, linkage["r2", keep])
      joint[members] <- c(sum(is.finite(raw[, lead])), linkage["joint_calls", keep])
    }
  }
  if (any(!nzchar(assigned))) stop("LD groups do not cover the fitted marker set.")
  groups <- split(seq_along(ids), assigned)
  coverage <- do.call(rbind, lapply(names(groups), function(name) {
    j <- groups[[name]]
    data.frame(region_id = name, chromosome = map$chromosome[j[1L]], region_start = min(map$position_bp[j]),
      region_end = max(map$position_bp[j]), anchor_marker_id = anchor[j[1L]], source_marker_count = length(j),
      fitted_marker_count = length(j), active_in_fitted_model = TRUE)
  }))
  membership <- data.frame(group_id = assigned, marker_id = ids, chromosome = map$chromosome, position_bp = map$position_bp,
    anchor_marker_id = anchor, search_distance_bp = distance, anchor_r2 = r2, joint_calls = joint)
  c(list(groups = groups, coverage = coverage, membership = membership), estimated)
}
