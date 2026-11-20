mean_impute_genotype_for_structure <- function(genotype) {
  if (!is.matrix(genotype) || !is.numeric(genotype)) {
    stop("genotype must be a numeric matrix.", call. = FALSE)
  }
  marker_means <- colMeans(genotype, na.rm = TRUE)
  if (any(!is.finite(marker_means))) {
    stop("Structure calculation cannot impute a marker with no callable genotypes.", call. = FALSE)
  }
  imputed <- genotype
  missing <- which(is.na(imputed), arr.ind = TRUE)
  if (nrow(missing)) {
    imputed[missing] <- marker_means[missing[, "col"]]
  }
  list(genotype = imputed, marker_means = marker_means)
}

prepare_structure_matrix <- function(dataset) {
  validate_genomic_dataset(dataset)
  imputed <- mean_impute_genotype_for_structure(dataset$genotype)
  allele_frequency <- imputed$marker_means / 2
  centered <- sweep(imputed$genotype, 2L, 2 * allele_frequency, FUN = "-")
  marker_scale <- sqrt(2 * allele_frequency * (1 - allele_frequency))
  usable <- is.finite(marker_scale) & marker_scale > 0
  if (!all(usable)) {
    stop("Post-QC structure input still contains zero-variance markers.", call. = FALSE)
  }
  standardized <- sweep(centered, 2L, marker_scale, FUN = "/")
  list(
    genotype = imputed$genotype,
    marker_means = imputed$marker_means,
    allele_frequency = allele_frequency,
    centered = centered,
    standardized = standardized
  )
}

ld_prune_structure_markers <- function(
    dataset,
    standardized,
    window_bp = 10000000,
    r2_max = 0.20) {
  validate_genomic_dataset(dataset)
  require_marker_map(dataset, "Physical-window LD pruning")
  if (!is.matrix(standardized) ||
      !identical(colnames(standardized), colnames(dataset$genotype))) {
    stop("standardized must align exactly to dataset markers.", call. = FALSE)
  }
  window_bp <- as.numeric(window_bp)
  r2_max <- as.numeric(r2_max)
  if (length(window_bp) != 1L || !is.finite(window_bp) || window_bp <= 0) {
    stop("window_bp must be one positive number.", call. = FALSE)
  }
  if (length(r2_max) != 1L || !is.finite(r2_max) || r2_max <= 0 || r2_max >= 1) {
    stop("r2_max must be strictly between zero and one.", call. = FALSE)
  }

  markers <- dataset$markers
  chromosome <- as.character(markers$chromosome)
  chromosome_numeric <- suppressWarnings(as.numeric(chromosome))
  order_index <- order(
    is.na(chromosome_numeric), chromosome_numeric, chromosome,
    markers$position_bp, markers$marker_id
  )
  kept <- integer()
  audit <- data.frame(
    marker_id = markers$marker_id,
    chromosome = chromosome,
    position_bp = markers$position_bp,
    kept_for_pca = FALSE,
    max_local_r2 = NA_real_,
    blocking_marker_id = NA_character_,
    stringsAsFactors = FALSE
  )

  for (marker_index in order_index) {
    local_kept <- kept[
      chromosome[kept] == chromosome[[marker_index]] &
        markers$position_bp[kept] >= markers$position_bp[[marker_index]] - window_bp
    ]
    if (!length(local_kept)) {
      kept <- c(kept, marker_index)
      audit$kept_for_pca[[marker_index]] <- TRUE
      next
    }

    current <- standardized[, marker_index]
    comparison <- standardized[, local_kept, drop = FALSE]
    denominator <- sqrt(sum(current^2) * colSums(comparison^2))
    correlations <- as.numeric(crossprod(current, comparison)) / denominator
    correlations[!is.finite(correlations)] <- 0
    r2 <- correlations^2
    maximum_index <- which.max(r2)
    audit$max_local_r2[[marker_index]] <- r2[[maximum_index]]
    if (r2[[maximum_index]] <= r2_max) {
      kept <- c(kept, marker_index)
      audit$kept_for_pca[[marker_index]] <- TRUE
    } else {
      audit$blocking_marker_id[[marker_index]] <- markers$marker_id[local_kept[[maximum_index]]]
    }
  }

  kept <- sort(kept)
  list(marker_ids = markers$marker_id[kept], audit = audit)
}

compute_population_structure <- function(dataset, n_pcs = 10L, pca_config = list()) {
  validate_genomic_dataset(dataset)
  prepared <- prepare_structure_matrix(dataset)
  pruning <- ld_prune_structure_markers(
    dataset,
    prepared$standardized,
    window_bp = pca_config$ld_window_bp %||% 10000000,
    r2_max = pca_config$ld_r2_max %||% 0.20
  )
  pca_matrix <- prepared$standardized[, pruning$marker_ids, drop = FALSE]
  n_pcs <- min(as.integer(n_pcs), nrow(pca_matrix) - 1L, ncol(pca_matrix))
  if (n_pcs < 3L) stop("At least three principal components are required.", call. = FALSE)

  pca <- stats::prcomp(
    pca_matrix,
    center = FALSE,
    scale. = FALSE,
    rank. = n_pcs,
    retx = TRUE
  )
  score_matrix <- pca$x[, seq_len(n_pcs), drop = FALSE]
  loading_matrix <- pca$rotation[, seq_len(n_pcs), drop = FALSE]
  scores <- data.frame(
    sample_id = rownames(score_matrix),
    breeding_material_type = dataset$samples$breeding_material_type[
      match(rownames(score_matrix), dataset$samples$sample_id)
    ],
    breeding_material_type_source = dataset$samples$breeding_material_type_source[
      match(rownames(score_matrix), dataset$samples$sample_id)
    ],
    score_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  loadings <- data.frame(
    marker_id = rownames(loading_matrix),
    loading_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  component_variance <- pca$sdev^2
  variance <- data.frame(
    component = paste0("PC", seq_len(n_pcs)),
    variance = component_variance[seq_len(n_pcs)],
    proportion = component_variance[seq_len(n_pcs)] / sum(component_variance),
    cumulative_proportion = cumsum(component_variance[seq_len(n_pcs)] / sum(component_variance)),
    stringsAsFactors = FALSE
  )

  grm_denominator <- sum(2 * prepared$allele_frequency * (1 - prepared$allele_frequency))
  if (!is.finite(grm_denominator) || grm_denominator <= 0) {
    stop("Invalid genomic relationship matrix denominator.", call. = FALSE)
  }
  grm <- tcrossprod(prepared$centered) / grm_denominator
  dimnames(grm) <- list(rownames(dataset$genotype), rownames(dataset$genotype))
  grm <- (grm + t(grm)) / 2

  squared_distance <- outer(diag(grm), diag(grm), "+") - 2 * grm
  squared_distance[squared_distance < 0 & squared_distance > -1e-10] <- 0
  if (any(squared_distance < 0)) {
    stop("GRM-derived squared genetic distance contains negative values.", call. = FALSE)
  }
  genetic_distance <- sqrt(squared_distance)
  dimnames(genetic_distance) <- dimnames(grm)

  off_diagonal <- grm[lower.tri(grm)]
  probabilities <- c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)
  relatedness_summary <- data.frame(
    statistic = c("mean", paste0("q", sprintf("%02d", probabilities * 100))),
    value = c(mean(off_diagonal), unname(stats::quantile(off_diagonal, probabilities))),
    stringsAsFactors = FALSE
  )

  list(
    scores = scores,
    loadings = loadings,
    variance = variance,
    grm = grm,
    genetic_distance = genetic_distance,
    allele_frequency = setNames(prepared$allele_frequency, colnames(dataset$genotype)),
    imputation_means = setNames(prepared$marker_means, colnames(dataset$genotype)),
    pca_marker_audit = pruning$audit,
    relatedness_summary = relatedness_summary,
    provenance = list(
      purpose = "genotype-only population structure and fold construction",
      imputation = "global marker mean for structure only",
      scaling = "centered dosage divided by sqrt(2p(1-p)) for PCA",
      pca_ld_pruning = paste0(
        "deterministic physical-window pruning: ",
        pca_config$ld_window_bp %||% 10000000, " bp; r2 <= ",
        pca_config$ld_r2_max %||% 0.20
      ),
      pca_markers_input = ncol(dataset$genotype),
      pca_markers_retained = length(pruning$marker_ids),
      grm = "VanRaden centered relationship matrix",
      predictive_preprocessing_reuse = FALSE
    )
  )
}

make_genomic_null_blocks <- function(markers, block_bp) {
  block_bp <- as.numeric(block_bp)
  if (length(block_bp) != 1L || !is.finite(block_bp) || block_bp <= 0) {
    stop("block_bp must be one positive number.", call. = FALSE)
  }
  paste(markers$chromosome, floor((markers$position_bp - 1) / block_bp), sep = ":")
}

screen_gwas_pc_candidates <- function(dataset, structure, selection_config, seed = 20260930L) {
  validate_genomic_dataset(dataset)
  required <- c("max_pcs", "null_permutations", "null_quantile", "null_block_bp")
  if (!all(required %in% names(selection_config))) {
    stop("GWAS PC selection configuration is incomplete.", call. = FALSE)
  }
  marker_ids <- structure$pca_marker_audit$marker_id[
    structure$pca_marker_audit$kept_for_pca
  ]
  marker_rows <- match(marker_ids, dataset$markers$marker_id)
  if (anyNA(marker_rows)) stop("PCA marker audit is not aligned to the dataset.", call. = FALSE)

  prepared <- prepare_structure_matrix(dataset)
  pca_matrix <- prepared$standardized[, marker_ids, drop = FALSE]
  max_pcs <- min(
    as.integer(selection_config$max_pcs),
    nrow(pca_matrix) - 1L,
    ncol(pca_matrix),
    nrow(structure$variance)
  )
  permutations <- as.integer(selection_config$null_permutations)
  null_quantile <- as.numeric(selection_config$null_quantile)
  if (max_pcs < 1L || permutations < 20L) {
    stop("GWAS PC selection needs at least one PC and 20 null permutations.", call. = FALSE)
  }
  if (!is.finite(null_quantile) || null_quantile <= 0.5 || null_quantile >= 1) {
    stop("null_quantile must be strictly between 0.5 and 1.", call. = FALSE)
  }

  blocks <- make_genomic_null_blocks(
    dataset$markers[marker_rows, , drop = FALSE],
    selection_config$null_block_bp
  )
  block_columns <- split(seq_along(marker_ids), blocks)
  null_eigenvalues <- matrix(NA_real_, nrow = permutations, ncol = max_pcs)
  null_eigenvalues <- with_seed_preserved(seed, {
    for (permutation_index in seq_len(permutations)) {
      null_matrix <- pca_matrix
      for (columns in block_columns) {
        row_permutation <- sample.int(nrow(pca_matrix))
        null_matrix[, columns] <- pca_matrix[row_permutation, columns, drop = FALSE]
      }
      singular_values <- La.svd(null_matrix, nu = 0L, nv = 0L)$d
      null_eigenvalues[permutation_index, ] <-
        singular_values[seq_len(max_pcs)]^2 / (nrow(pca_matrix) - 1L)
    }
    null_eigenvalues
  })

  observed <- structure$variance$variance[seq_len(max_pcs)]
  threshold <- apply(
    null_eigenvalues,
    2L,
    stats::quantile,
    probs = null_quantile,
    names = FALSE,
    type = 8
  )
  exceeds <- observed > threshold
  selected_n <- if (all(exceeds)) {
    max_pcs
  } else {
    which(!exceeds)[[1L]] - 1L
  }
  comparison <- data.frame(
    component = paste0("PC", seq_len(max_pcs)),
    observed_eigenvalue = observed,
    null_mean = colMeans(null_eigenvalues),
    null_sd = apply(null_eigenvalues, 2L, stats::sd),
    null_quantile = threshold,
    exceeds_null = exceeds,
    supported_contiguous = seq_len(max_pcs) <= selected_n,
    stringsAsFactors = FALSE
  )

  list(
    supported_n_pcs = as.integer(selected_n),
    selection_saturated = selected_n == max_pcs,
    comparison = comparison,
    provenance = list(
      phenotype_values_used = FALSE,
      method = "contiguous block-permutation parallel analysis after LD pruning",
      rule = paste0(
        "retain consecutive leading eigenvalues above the ",
        100 * null_quantile, "th percentile of the genotype-only null"
      ),
      null_blocks = paste0(selection_config$null_block_bp, " bp chromosome blocks"),
      null_permutations = permutations,
      seed = as.integer(seed),
      final_observed_gwas_inflation_used_for_selection = FALSE
    )
  )
}

profile_null_mixed_model <- function(y, fixed_effects, grm) {
  if (!is.numeric(y) || anyNA(y) || any(!is.finite(y))) {
    stop("y must contain finite, non-missing numeric values.", call. = FALSE)
  }
  fixed_effects <- as.matrix(fixed_effects)
  if (!is.numeric(fixed_effects) || nrow(fixed_effects) != length(y) || anyNA(fixed_effects)) {
    stop("fixed_effects must be a complete numeric matrix aligned to y.", call. = FALSE)
  }
  if (!is.matrix(grm) || !is.numeric(grm) ||
      !identical(dim(grm), c(length(y), length(y))) || anyNA(grm)) {
    stop("grm must be a complete square numeric matrix aligned to y.", call. = FALSE)
  }
  if (qr(fixed_effects)$rank != ncol(fixed_effects)) {
    stop("Null-model fixed effects are rank deficient.", call. = FALSE)
  }

  decomposition <- eigen((grm + t(grm)) / 2, symmetric = TRUE)
  tolerance <- max(1, max(abs(decomposition$values))) * 1e-8
  if (any(decomposition$values < -tolerance)) {
    stop("The genomic relationship matrix is not positive semidefinite.", call. = FALSE)
  }
  eigenvalues <- pmax(decomposition$values, 0)
  rotated_y <- as.numeric(crossprod(decomposition$vectors, y))
  rotated_x <- crossprod(decomposition$vectors, fixed_effects)
  n <- length(y)

  negative_log_likelihood <- function(log_delta) {
    delta <- exp(log_delta)
    inverse_variance <- 1 / (eigenvalues + delta)
    weighted_x <- rotated_x * sqrt(inverse_variance)
    weighted_y <- rotated_y * sqrt(inverse_variance)
    beta <- tryCatch(
      qr.solve(weighted_x, weighted_y, tol = 1e-10),
      error = function(error) NULL
    )
    if (is.null(beta)) return(Inf)
    residual <- rotated_y - as.numeric(rotated_x %*% beta)
    quadratic <- sum(residual^2 * inverse_variance)
    if (!is.finite(quadratic) || quadratic <= 0) return(Inf)
    log_likelihood <- -0.5 * (
      n * (log(2 * pi) + 1 + log(quadratic / n)) +
        sum(log(eigenvalues + delta))
    )
    -log_likelihood
  }

  optimum <- stats::optimize(
    negative_log_likelihood,
    interval = log(c(1e-5, 1e5)),
    tol = 1e-8
  )
  log_likelihood <- -optimum$objective
  parameter_count <- ncol(fixed_effects) + 2L
  list(
    log_likelihood = log_likelihood,
    bic = -2 * log_likelihood + parameter_count * log(n),
    delta = exp(optimum$minimum),
    parameter_count = parameter_count
  )
}

select_trait_gwas_pc_counts <- function(
    dataset,
    structure,
    pc_screen,
    selection_config = list()) {
  validate_genomic_dataset(dataset)
  if (is.null(dataset$phenotypes) || !nrow(dataset$phenotypes)) {
    stop("Phenotypes are required for trait-specific null-model selection.", call. = FALSE)
  }
  phenotype_key <- paste(dataset$phenotypes$sample_id, dataset$phenotypes$trait, sep = "\r")
  if (anyDuplicated(phenotype_key)) {
    stop("Each sample/trait phenotype combination must be unique.", call. = FALSE)
  }
  fixed_pc_cap <- as.integer(selection_config$max_fixed_pcs %||% 10L)
  if (length(fixed_pc_cap) != 1L || is.na(fixed_pc_cap) || fixed_pc_cap < 0L) {
    stop("max_fixed_pcs must be one non-negative integer.", call. = FALSE)
  }
  candidate_max <- min(as.integer(pc_screen$supported_n_pcs), fixed_pc_cap)
  if (candidate_max < 0L) stop("Invalid PC screen result.", call. = FALSE)
  available <- length(grep("^PC[0-9]+$", names(structure$scores)))
  candidate_max <- min(candidate_max, available)
  traits <- sort(unique(as.character(dataset$phenotypes$trait)))
  rows <- list()
  row_index <- 0L

  for (trait_name in traits) {
    phenotype <- dataset$phenotypes[
      dataset$phenotypes$trait == trait_name & is.finite(dataset$phenotypes$value),
      c("sample_id", "value"),
      drop = FALSE
    ]
    score_rows <- match(phenotype$sample_id, structure$scores$sample_id)
    if (anyNA(score_rows)) stop("A phenotype sample is absent from PCA scores.", call. = FALSE)
    grm_rows <- match(phenotype$sample_id, rownames(structure$grm))
    if (anyNA(grm_rows)) stop("A phenotype sample is absent from the GRM.", call. = FALSE)
    y <- as.numeric(phenotype$value)
    trait_max <- min(candidate_max, length(y) - 4L)
    if (trait_max < 0L) stop("Too few observations for null-model PC selection.", call. = FALSE)

    for (n_pcs in 0:trait_max) {
      fixed_effects <- matrix(1, nrow = length(y), ncol = 1L)
      if (n_pcs > 0L) {
        pc_names <- paste0("PC", seq_len(n_pcs))
        fixed_effects <- cbind(
          intercept = 1,
          as.matrix(structure$scores[score_rows, pc_names, drop = FALSE])
        )
      } else {
        colnames(fixed_effects) <- "intercept"
      }
      fit <- profile_null_mixed_model(
        y,
        fixed_effects,
        structure$grm[grm_rows, grm_rows, drop = FALSE]
      )
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        trait = trait_name,
        n_observations = length(y),
        n_pcs = n_pcs,
        log_likelihood = fit$log_likelihood,
        bic = fit$bic,
        delta = fit$delta,
        stringsAsFactors = FALSE
      )
    }
  }

  comparison <- do.call(rbind, rows)
  comparison$selected <- FALSE
  selected_rows <- integer()
  for (trait_name in traits) {
    indices <- which(comparison$trait == trait_name)
    selected_index <- indices[order(comparison$bic[indices], comparison$n_pcs[indices])[[1L]]]
    comparison$selected[[selected_index]] <- TRUE
    selected_rows <- c(selected_rows, selected_index)
  }
  selected <- comparison[selected_rows, c("trait", "n_observations", "n_pcs", "bic", "delta")]
  names(selected)[names(selected) == "n_pcs"] <- "selected_n_pcs"
  rownames(selected) <- NULL
  selected$bic_margin_to_next <- vapply(selected$trait, function(trait_name) {
    trait_bic <- sort(comparison$bic[comparison$trait == trait_name])
    if (length(trait_bic) < 2L) Inf else trait_bic[[2L]] - trait_bic[[1L]]
  }, numeric(1))
  selected$zero_pc_selected <- selected$selected_n_pcs == 0L
  selected$upper_cap_selected <- selected$selected_n_pcs == candidate_max
  selected$near_tie_delta_bic_lt_2 <- selected$bic_margin_to_next < 2

  list(
    selected_by_trait = selected,
    comparison = comparison,
    provenance = list(
      phenotype_values_used = TRUE,
      marker_association_statistics_used = FALSE,
      method = "minimum BIC from a pre-association genomic null mixed model",
      null_model = "phenotype ~ leading LD-pruned PCs + genomic random effect (VanRaden GRM)",
      variance_fit = "maximum likelihood profile over delta = residual/genomic variance",
      maximum_fixed_pcs = fixed_pc_cap,
      candidate_screen_saturated = isTRUE(pc_screen$selection_saturated),
      final_observed_gwas_inflation_used_for_selection = FALSE,
      downstream_requirement = "freeze per-trait PC count before SNP tests; validate calibration and power in simulation"
    )
  )
}

balanced_pc_partition <- function(score_matrix, k, seed, max_iterations = 50L) {
  n <- nrow(score_matrix)
  if (k < 2L || k > n) stop("Invalid number of clusters.", call. = FALSE)
  initial <- with_seed_preserved(
    seed,
    stats::kmeans(score_matrix, centers = k, nstart = 100L, iter.max = 100L)
  )
  centers <- initial$centers
  capacities <- rep(n %/% k, k)
  capacities[seq_len(n %% k)] <- capacities[seq_len(n %% k)] + 1L
  previous <- rep.int(NA_integer_, n)

  for (iteration in seq_len(as.integer(max_iterations))) {
    squared_distance <- vapply(seq_len(k), function(cluster_index) {
      rowSums((score_matrix - matrix(
        centers[cluster_index, ],
        nrow = n,
        ncol = ncol(score_matrix),
        byrow = TRUE
      ))^2)
    }, numeric(n))
    preference_order <- t(apply(squared_distance, 1L, order))
    best <- squared_distance[cbind(seq_len(n), preference_order[, 1L])]
    second <- squared_distance[cbind(seq_len(n), preference_order[, 2L])]
    margin <- second - best
    sample_order <- order(-margin, rownames(score_matrix))

    assignment <- integer(n)
    remaining <- capacities
    for (sample_index in sample_order) {
      preferred <- preference_order[sample_index, ]
      chosen <- preferred[remaining[preferred] > 0L][[1L]]
      assignment[[sample_index]] <- chosen
      remaining[[chosen]] <- remaining[[chosen]] - 1L
    }
    if (identical(assignment, previous)) break
    previous <- assignment
    centers <- vapply(seq_len(k), function(cluster_index) {
      colMeans(score_matrix[assignment == cluster_index, , drop = FALSE])
    }, numeric(ncol(score_matrix)))
    centers <- t(centers)
  }

  center_order <- do.call(order, as.data.frame(centers, check.names = FALSE))
  fold_by_cluster <- setNames(seq_len(k), center_order)
  as.integer(unname(fold_by_cluster[as.character(assignment)]))
}

make_genomic_blocked_folds <- function(
    structure,
    sample_ids,
    k = 5L,
    repeats = 1L,
    seed = 20260930L,
    n_pcs = 10L) {
  assert_unique_nonmissing(sample_ids, "sample_ids")
  k <- as.integer(k)
  repeats <- as.integer(repeats)
  if (k < 2L || k > length(sample_ids)) stop("Invalid number of folds.", call. = FALSE)
  if (repeats < 1L) stop("repeats must be positive.", call. = FALSE)

  score_rows <- match(sample_ids, structure$scores$sample_id)
  if (anyNA(score_rows)) {
    stop("Some fold sample IDs are absent from population-structure scores.", call. = FALSE)
  }
  pc_columns <- grep("^PC[0-9]+$", names(structure$scores), value = TRUE)
  n_pcs <- min(as.integer(n_pcs), length(pc_columns))
  if (n_pcs < 2L) stop("Genomic blocking requires at least two PCs.", call. = FALSE)
  score_matrix <- as.matrix(structure$scores[score_rows, pc_columns[seq_len(n_pcs)], drop = FALSE])
  rownames(score_matrix) <- sample_ids

  pieces <- vector("list", repeats)
  for (repeat_index in seq_len(repeats)) {
    current_seed <- as.integer(seed + repeat_index - 1L)
    fold <- balanced_pc_partition(score_matrix, k = k, seed = current_seed)
    pieces[[repeat_index]] <- data.frame(
      sample_id = sample_ids,
      scheme = "genomic_blocked",
      repeat_id = repeat_index,
      fold = as.integer(fold),
      seed = current_seed,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, pieces)
}

validate_fold_registry <- function(folds, sample_ids, k) {
  required <- c("sample_id", "scheme", "repeat_id", "fold", "seed")
  if (!is.data.frame(folds) || !all(required %in% names(folds))) {
    stop("folds is missing required columns.", call. = FALSE)
  }
  assert_unique_nonmissing(sample_ids, "sample_ids")
  combinations <- unique(folds[, c("scheme", "repeat_id"), drop = FALSE])
  for (index in seq_len(nrow(combinations))) {
    current <- folds[
      folds$scheme == combinations$scheme[[index]] &
        folds$repeat_id == combinations$repeat_id[[index]],
      ,
      drop = FALSE
    ]
    if (anyDuplicated(current$sample_id) ||
        !setequal(as.character(current$sample_id), as.character(sample_ids))) {
      stop("Each scheme/repeat must assign every sample exactly once.", call. = FALSE)
    }
    if (!setequal(current$fold, seq_len(k))) {
      stop("Each scheme/repeat must contain every requested fold.", call. = FALSE)
    }
    for (fold_index in seq_len(k)) {
      split <- fold_split_ids(current, combinations$repeat_id[[index]], fold_index)
      assert_no_split_overlap(split$train, split$test)
    }
  }
  invisible(TRUE)
}

audit_fold_relatedness <- function(folds, structure) {
  fold_sizes <- stats::aggregate(
    sample_id ~ scheme + repeat_id + fold,
    data = folds,
    FUN = length
  )
  names(fold_sizes)[names(fold_sizes) == "sample_id"] <- "n_test"
  rows <- vector("list", nrow(fold_sizes))
  for (index in seq_len(nrow(fold_sizes))) {
    current <- folds[
      folds$scheme == fold_sizes$scheme[[index]] &
        folds$repeat_id == fold_sizes$repeat_id[[index]],
      ,
      drop = FALSE
    ]
    split <- fold_split_ids(current, fold_sizes$repeat_id[[index]], fold_sizes$fold[[index]])
    between <- structure$grm[split$test, split$train, drop = FALSE]
    closest_training_relative <- apply(between, 1L, max)
    rows[[index]] <- data.frame(
      scheme = fold_sizes$scheme[[index]],
      repeat_id = fold_sizes$repeat_id[[index]],
      fold = fold_sizes$fold[[index]],
      mean_max_test_train_grm = mean(closest_training_relative),
      median_max_test_train_grm = stats::median(closest_training_relative),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

select_cv_blocking_pc_count <- function(
    structure,
    sample_ids,
    k,
    repeats,
    seed,
    pc_candidates) {
  available <- length(grep("^PC[0-9]+$", names(structure$scores)))
  candidates <- sort(unique(as.integer(pc_candidates)))
  candidates <- candidates[candidates >= 2L & candidates <= available]
  if (!length(candidates)) stop("No valid CV PC candidates remain.", call. = FALSE)
  folds_by_candidate <- vector("list", length(candidates))
  audit_rows <- vector("list", length(candidates))

  for (index in seq_along(candidates)) {
    folds_by_candidate[[index]] <- make_genomic_blocked_folds(
      structure,
      sample_ids,
      k = k,
      repeats = repeats,
      seed = seed,
      n_pcs = candidates[[index]]
    )
    audit <- audit_fold_relatedness(folds_by_candidate[[index]], structure)
    sizes <- table(folds_by_candidate[[index]]$repeat_id, folds_by_candidate[[index]]$fold)
    audit_rows[[index]] <- data.frame(
      n_pcs = candidates[[index]],
      mean_fold_mean_max_test_train_grm = mean(audit$mean_max_test_train_grm),
      worst_fold_mean_max_test_train_grm = max(audit$mean_max_test_train_grm),
      min_test_n = min(sizes),
      max_test_n = max(sizes),
      stringsAsFactors = FALSE
    )
  }
  comparison <- do.call(rbind, audit_rows)
  selected_index <- order(
    comparison$mean_fold_mean_max_test_train_grm,
    comparison$worst_fold_mean_max_test_train_grm,
    comparison$n_pcs
  )[[1L]]
  comparison$selected <- seq_len(nrow(comparison)) == selected_index
  list(
    selected_n_pcs = comparison$n_pcs[[selected_index]],
    comparison = comparison,
    folds = folds_by_candidate[[selected_index]]
  )
}

make_outer_fold_registry <- function(
    dataset,
    structure,
    resampling,
    seed = 20260930L,
    pc_candidates = c(2L, 3L, 5L, 10L)) {
  validate_genomic_dataset(dataset)
  phenotype_ids <- unique(as.character(dataset$phenotypes$sample_id))
  sample_ids <- dataset$samples$sample_id[dataset$samples$sample_id %in% phenotype_ids]
  assert_unique_nonmissing(sample_ids, "phenotyped sample IDs")

  k <- as.integer(resampling$outer_folds)
  repeats <- as.integer(resampling$outer_repeats_development)
  random <- make_repeated_random_folds(
    sample_ids,
    k = k,
    repeats = repeats,
    seed = seed
  )
  cv_selection <- select_cv_blocking_pc_count(
    structure,
    sample_ids,
    k = k,
    repeats = repeats,
    seed = seed + 10000L,
    pc_candidates = pc_candidates
  )
  blocked <- cv_selection$folds
  folds <- rbind(random, blocked)
  rownames(folds) <- NULL
  sample_rows <- match(folds$sample_id, dataset$samples$sample_id)
  folds$breeding_material_type <- dataset$samples$breeding_material_type[sample_rows]
  folds$breeding_material_type_source <-
    dataset$samples$breeding_material_type_source[sample_rows]
  validate_fold_registry(folds, sample_ids, k)

  fold_sizes <- stats::aggregate(
    sample_id ~ scheme + repeat_id + fold,
    data = folds,
    FUN = length
  )
  names(fold_sizes)[names(fold_sizes) == "sample_id"] <- "n_test"

  trait_membership <- unique(dataset$phenotypes[, c("sample_id", "trait"), drop = FALSE])
  trait_fold_counts <- merge(folds, trait_membership, by = "sample_id", all.x = FALSE, all.y = FALSE)
  trait_fold_counts <- stats::aggregate(
    sample_id ~ trait + scheme + repeat_id + fold,
    data = trait_fold_counts,
    FUN = length
  )
  names(trait_fold_counts)[names(trait_fold_counts) == "sample_id"] <- "n_test"

  list(
    folds = folds,
    fold_sizes = fold_sizes,
    trait_fold_counts = trait_fold_counts,
    relatedness_audit = audit_fold_relatedness(folds, structure),
    cv_pc_selection = cv_selection$comparison,
    provenance = list(
      generated_once = TRUE,
      phenotype_values_used = FALSE,
      eligibility_uses_nonmissing_phenotype_ids = TRUE,
      blocked_group_interpretation = "genomic similarity stress-test blocks, not biological families",
      blocked_algorithm = paste0(
        "capacity-balanced clustering; genotype-only GRM criterion selected ",
        cv_selection$selected_n_pcs, " PCs"
      ),
      blocked_pc_candidates = sort(unique(as.integer(pc_candidates))),
      blocked_selected_n_pcs = cv_selection$selected_n_pcs,
      k = k,
      repeats = repeats,
      random_seed = as.integer(seed),
      blocked_seed = as.integer(seed + 10000L)
    )
  )
}

draw_pca_inspection <- function(structure, registry) {
  pc_columns <- paste0("PC", 1:3)
  if (!all(pc_columns %in% names(structure$scores))) {
    stop("The first three PCs are required for the inspection figure.", call. = FALSE)
  }
  blocked <- registry$folds[
    registry$folds$scheme == "genomic_blocked" & registry$folds$repeat_id == 1L,
    c("sample_id", "fold"),
    drop = FALSE
  ]
  fold <- blocked$fold[match(structure$scores$sample_id, blocked$sample_id)]
  palette <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9")
  if (max(blocked$fold) > length(palette)) {
    palette <- grDevices::hcl.colors(max(blocked$fold), palette = "Dark 3")
  }
  point_colour <- rep("#8A8A8A", length(fold))
  assigned <- !is.na(fold)
  point_colour[assigned] <- palette[fold[assigned]]
  point_shape <- ifelse(assigned, 21L, 4L)
  variance <- setNames(structure$variance$proportion, structure$variance$component)
  pairs <- list(c("PC1", "PC2"), c("PC1", "PC3"), c("PC2", "PC3"))

  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(1L, 3L), mar = c(4.2, 4.3, 2.4, 1), oma = c(3, 0, 3, 0))
  for (index in seq_along(pairs)) {
    x_name <- pairs[[index]][[1L]]
    y_name <- pairs[[index]][[2L]]
    graphics::plot(
      structure$scores[[x_name]],
      structure$scores[[y_name]],
      type = "n",
      xlab = sprintf("%s (%.1f%%)", x_name, 100 * variance[[x_name]]),
      ylab = sprintf("%s (%.1f%%)", y_name, 100 * variance[[y_name]]),
      main = paste(x_name, "vs", y_name),
      bty = "l"
    )
    graphics::abline(h = 0, v = 0, col = "#D9D9D9", lty = 3)
    graphics::points(
      structure$scores[[x_name]],
      structure$scores[[y_name]],
      pch = point_shape,
      bg = point_colour,
      col = ifelse(assigned, "white", "#666666"),
      cex = 0.85,
      lwd = 0.5
    )
    if (index == 3L) {
      n_folds <- max(blocked$fold)
      legend_labels <- paste("CV block", seq_len(n_folds))
      legend_colours <- palette[seq_len(n_folds)]
      legend_border <- rep("white", n_folds)
      legend_shape <- rep(21L, n_folds)
      if (any(!assigned)) {
        legend_labels <- c(legend_labels, "No phenotype/CV fold")
        legend_colours <- c(legend_colours, "#8A8A8A")
        legend_border <- c(legend_border, "#666666")
        legend_shape <- c(legend_shape, 4L)
      }
      graphics::legend(
        "topright",
        legend = legend_labels,
        pt.bg = legend_colours,
        col = legend_border,
        pch = legend_shape,
        pt.cex = 0.9,
        cex = 0.72,
        bty = "n"
      )
    }
  }
  graphics::mtext(
    "LD-pruned genotype PCA: first three components",
    outer = TRUE,
    side = 3,
    line = 1,
    font = 2
  )
  graphics::mtext(
    "Colours show genotype-only CV stress-test blocks; they are not biological groups.",
    outer = TRUE,
    side = 1,
    line = 1,
    cex = 0.82
  )
  invisible(NULL)
}

write_pca_inspection_figure <- function(
    structure,
    registry,
    pdf_path = file.path("results", "figures", "tassel_tutorial_pca_first3.pdf"),
    png_path = file.path("results", "figures", "tassel_tutorial_pca_first3.png")) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(pdf_path, width = 12, height = 4.5, family = "sans", useDingbats = FALSE)
  tryCatch(draw_pca_inspection(structure, registry), finally = grDevices::dev.off())
  grDevices::png(png_path, width = 2400, height = 900, res = 200, type = "windows")
  tryCatch(draw_pca_inspection(structure, registry), finally = grDevices::dev.off())
  c(pdf_path, png_path)
}

write_population_structure_outputs <- function(
    structure,
    registry,
    gwas_pc_screen,
    gwas_pc_selection,
    processed_directory = file.path("data", "processed", "tassel_tutorial"),
    table_directory = file.path("results", "tables")) {
  dir.create(processed_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
  output_paths <- c(
    structure = file.path(processed_directory, "population_structure.rds"),
    folds_rds = file.path(processed_directory, "outer_folds.rds"),
    pca_scores = file.path(table_directory, "tassel_tutorial_pca_scores.tsv"),
    pca_variance = file.path(table_directory, "tassel_tutorial_pca_variance.tsv"),
    grm = file.path(table_directory, "tassel_tutorial_grm.tsv"),
    relatedness_summary = file.path(table_directory, "tassel_tutorial_relatedness_summary.tsv"),
    outer_folds = file.path(table_directory, "tassel_tutorial_outer_folds.tsv"),
    fold_sizes = file.path(table_directory, "tassel_tutorial_fold_sizes.tsv"),
    trait_fold_counts = file.path(table_directory, "tassel_tutorial_trait_fold_counts.tsv"),
    fold_relatedness = file.path(table_directory, "tassel_tutorial_fold_relatedness.tsv"),
    pca_marker_audit = file.path(table_directory, "tassel_tutorial_pca_marker_audit.tsv"),
    gwas_pc_parallel_analysis = file.path(table_directory, "tassel_tutorial_gwas_pc_parallel_analysis.tsv"),
    gwas_pc_selection = file.path(table_directory, "tassel_tutorial_gwas_pc_selection.tsv"),
    gwas_pc_selection_summary = file.path(table_directory, "tassel_tutorial_gwas_pc_selection_summary.tsv"),
    cv_pc_selection = file.path(table_directory, "tassel_tutorial_cv_pc_selection.tsv")
  )

  saveRDS(structure, output_paths[["structure"]], version = 3)
  saveRDS(registry, output_paths[["folds_rds"]], version = 3)
  grm_table <- data.frame(
    sample_id = rownames(structure$grm),
    structure$grm,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  tables <- list(
    pca_scores = structure$scores,
    pca_variance = structure$variance,
    grm = grm_table,
    relatedness_summary = structure$relatedness_summary,
    outer_folds = registry$folds,
    fold_sizes = registry$fold_sizes,
    trait_fold_counts = registry$trait_fold_counts,
    fold_relatedness = registry$relatedness_audit,
    pca_marker_audit = structure$pca_marker_audit,
    gwas_pc_parallel_analysis = gwas_pc_screen$comparison,
    gwas_pc_selection = gwas_pc_selection$comparison,
    gwas_pc_selection_summary = gwas_pc_selection$selected_by_trait,
    cv_pc_selection = registry$cv_pc_selection
  )
  for (table_name in names(tables)) {
    utils::write.table(
      tables[[table_name]],
      output_paths[[table_name]],
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  unname(output_paths)
}
