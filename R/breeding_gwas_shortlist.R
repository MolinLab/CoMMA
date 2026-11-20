breeding_gwas_shortlist_methods <- c("MLM", "FarmCPU", "BLINK")

breeding_check_gwas_candidates <- function(max_candidates, models) {
  if (is.null(max_candidates)) return(invisible(NULL))
  if (!is.numeric(max_candidates) || length(max_candidates) != 1L || !is.finite(max_candidates) ||
      max_candidates < 1 || max_candidates != floor(max_candidates)) {
    stop("Set gwas_candidates to a positive whole number, such as 10.", call. = FALSE)
  }
  if (anyDuplicated(toupper(models)) || !setequal(toupper(models), toupper(breeding_gwas_shortlist_methods))) {
    stop("A candidate shortlist requires MLM, FarmCPU and BLINK in the same GWAS run.", call. = FALSE)
  }
}

breeding_gwas_rank_regions <- function(scans, markers) {
  map <- breeding_coordinate_map(markers)
  methods <- breeding_gwas_shortlist_methods
  breeding_check_gwas_candidates(1L, unique(scans$method))
  stopifnot(!anyDuplicated(scans[c("method", "marker_id")]), is.numeric(scans$p_value),
            all(is.na(scans$p_value) | (scans$p_value >= 0 & scans$p_value <= 1)))
  p <- vapply(methods, function(method) {
    x <- scans[scans$method == method, , drop = FALSE]
    if (!setequal(x$marker_id, map$marker_id)) stop("GWAS methods must share the same mapped SNPs.")
    x$p_value[match(map$marker_id, x$marker_id)]
  }, numeric(nrow(map)))
  dim(p) <- c(nrow(map), length(methods))
  colnames(p) <- methods
  available_p <- p
  available_p[is.na(available_p)] <- Inf
  regions <- split(seq_len(nrow(map)), map$region_id)
  region_p <- t(vapply(regions, function(j) apply(available_p[j, , drop = FALSE], 2L, min), numeric(length(methods))))
  colnames(region_p) <- methods
  rank_methods <- function(values) {
    for (method in methods) {
      available <- is.finite(values[, method])
      values[available, method] <- rank(values[available, method], ties.method = "average")
    }
    values
  }
  region_ranks <- rank_methods(region_p)
  region_median <- apply(region_ranks, 1L, stats::median)
  eligible <- which(is.finite(region_median))
  ordered <- eligible[order(region_median[eligible], region_ranks[eligible, "MLM"], names(regions)[eligible], method = "radix")]

  marker_ranks <- rank_methods(available_p)
  marker_median <- apply(marker_ranks, 1L, stats::median)
  # A window can be supported by two scans whose tested markers do not overlap.
  fallback <- ifelse(is.finite(marker_median), 0, apply(marker_ranks, 1L, min))
  leads <- vapply(regions[ordered], function(j) {
    j[order(marker_median[j], fallback[j], marker_ranks[j, "MLM"], map$marker_id[j], method = "radix")[1L]]
  }, integer(1L))
  ranks <- region_ranks[ordered, , drop = FALSE]
  ranks[!is.finite(ranks)] <- NA_real_
  colnames(ranks) <- paste0("region_rank_", methods)
  lead_p <- p[leads, , drop = FALSE]
  colnames(lead_p) <- paste0("lead_p_", methods)
  data.frame(region_id = map$region_id[leads], chromosome = map$chromosome[leads],
             region_start = map$region_start[leads], region_end = map$region_end[leads],
             lead_marker_id = map$marker_id[leads], lead_position_bp = map$position_bp[leads],
             median_region_rank = region_median[ordered], ranks, lead_p, row.names = NULL)
}

breeding_gwas_clump_regions <- function(ranking, raw, max_candidates) {
  breeding_check_gwas_candidates(max_candidates, breeding_gwas_shortlist_methods)
  stopifnot(!is.null(max_candidates), is.matrix(raw), is.numeric(raw), !is.null(rownames(raw)),
            !anyDuplicated(rownames(raw)), !anyDuplicated(colnames(raw)),
            !anyDuplicated(ranking$region_id), all(ranking$lead_marker_id %in% colnames(raw)))
  kept <- integer()
  unresolved <- integer(nrow(ranking))
  for (i in seq_len(nrow(ranking))) {
    nearby <- kept[ranking$chromosome[kept] == ranking$chromosome[i] &
                     abs(ranking$lead_position_bp[kept] - ranking$lead_position_bp[i]) <= 1e7]
    if (length(nearby)) {
      a <- raw[, ranking$lead_marker_id[i]]
      r2 <- vapply(nearby, function(j) {
        b <- raw[, ranking$lead_marker_id[j]]
        called <- is.finite(a) & is.finite(b)
        if (sum(called) < 20L) return(NA_real_)
        suppressWarnings(stats::cor(a[called], b[called])^2)
      }, numeric(1L))
      if (any(r2 >= .2, na.rm = TRUE)) next
      unresolved[i] <- sum(is.na(r2))
    }
    kept <- c(kept, i)
    if (length(kept) == max_candidates) break
  }
  data.frame(priority = seq_along(kept), ranking[kept, , drop = FALSE],
             unresolved_nearby_ld_pairs = unresolved[kept], row.names = NULL)
}
