breeder_check_model_matrix <- function(z, y = NULL) {
  if (!is.matrix(z) || !is.numeric(z) || nrow(z) < 3L || ncol(z) < 1L ||
      any(!is.finite(z)) || is.null(rownames(z)) || is.null(colnames(z))) {
    stop("Model input must be a complete numeric genotype matrix with IDs.",
         call. = FALSE)
  }
  assert_unique_nonmissing(rownames(z), "sample IDs")
  assert_unique_nonmissing(colnames(z), "marker IDs")
  if (!is.null(y) && (length(y) != nrow(z) || any(!is.finite(y)))) {
    stop("Training phenotype must be finite and aligned by sample ID.",
         call. = FALSE)
  }
  invisible(TRUE)
}

breeder_squared_distances <- function(x, reference) {
  x_norm <- rowSums(x^2)
  reference_norm <- rowSums(reference^2)
  distances <- outer(x_norm, reference_norm, "+") -
    2 * tcrossprod(x, reference)
  pmax(distances, 0)
}

fit_gaussian_rkhs <- function(z_train, y_train, gamma = 1, lambda = 1) {
  breeder_check_model_matrix(z_train, y_train)
  if (length(gamma) != 1L || !is.finite(gamma) || gamma <= 0 ||
      length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
    stop("RKHS gamma and lambda must be positive finite scalars.", call. = FALSE)
  }
  distances <- breeder_squared_distances(z_train, z_train)
  upper <- distances[upper.tri(distances)]
  positive <- upper[upper > 1e-12]
  if (!length(positive)) {
    stop("RKHS training lines have no positive genetic distance.", call. = FALSE)
  }
  distance_scale <- stats::median(positive)
  kernel <- exp(-gamma * distances / distance_scale)
  covariance <- kernel + diag(lambda, nrow(kernel))
  inverse_one <- solve(covariance, rep(1, nrow(kernel)))
  intercept <- sum(inverse_one * y_train) / sum(inverse_one)
  alpha <- solve(covariance, y_train - intercept)
  list(
    marker_ids = colnames(z_train),
    training_ids = rownames(z_train),
    z_train = z_train,
    gamma = as.numeric(gamma),
    lambda = as.numeric(lambda),
    distance_scale = distance_scale,
    intercept = intercept,
    alpha = as.numeric(alpha),
    fitted_on_training_only = TRUE
  )
}

predict_gaussian_rkhs <- function(model, z_test) {
  if (!is.matrix(z_test) || !is.numeric(z_test) ||
      !identical(colnames(z_test), model$marker_ids) ||
      any(!is.finite(z_test))) {
    stop("RKHS prediction matrix must have the fitted markers in order.",
         call. = FALSE)
  }
  distances <- breeder_squared_distances(z_test, model$z_train)
  kernel <- exp(-model$gamma * distances / model$distance_scale)
  predicted <- as.numeric(model$intercept + kernel %*% model$alpha)
  if (any(!is.finite(predicted))) {
    stop("RKHS produced non-finite predictions.", call. = FALSE)
  }
  predicted
}

fit_bayesc <- function(z_train, y_train, n_iter = 1500L,
                       burn_in = 500L, thin = 5L, prob_in = 0.05,
                       seed = 1L) {
  breeder_check_model_matrix(z_train, y_train)
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("Project-local package 'BGLR' is required for BayesC.", call. = FALSE)
  }
  counts <- c(n_iter, burn_in, thin, seed)
  if (anyNA(counts) || any(!is.finite(counts)) ||
      any(counts != as.integer(counts)) || n_iter < 20L ||
      burn_in < 1L || burn_in >= n_iter || thin < 1L ||
      length(prob_in) != 1L || !is.finite(prob_in) ||
      prob_in <= 0 || prob_in >= 1) {
    stop("BayesC iteration, seed, or inclusion-prior settings are invalid.",
         call. = FALSE)
  }
  r_user <- Sys.getenv("R_USER", unset = "")
  if (!identical(basename(r_user), ".r_user") ||
      !file.exists(file.path(dirname(r_user), "AGENTS.md"))) {
    stop("BayesC must run through the project-local R launcher.",
         call. = FALSE)
  }
  scratch <- file.path(dirname(r_user), ".tmp", "breeder_bayesc")
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  save_prefix <- tempfile(pattern = "chain_", tmpdir = scratch)
  on.exit({
    generated <- list.files(scratch, full.names = TRUE)
    generated <- generated[startsWith(basename(generated), basename(save_prefix))]
    if (length(generated)) file.remove(generated)
  }, add = TRUE)
  fit <- with_seed_preserved(seed, BGLR::BGLR(
    y = as.numeric(y_train),
    ETA = list(list(X = z_train, model = "BayesC", probIn = prob_in)),
    nIter = as.integer(n_iter), burnIn = as.integer(burn_in),
    thin = as.integer(thin), saveAt = save_prefix,
    verbose = FALSE, rmExistingFiles = TRUE
  ))
  effects <- as.numeric(fit$ETA[[1L]]$b)
  intercept <- as.numeric(fit$mu)
  if (length(effects) != ncol(z_train) || any(!is.finite(effects)) ||
      length(intercept) != 1L || !is.finite(intercept)) {
    stop("BayesC returned invalid posterior mean coefficients.", call. = FALSE)
  }
  names(effects) <- colnames(z_train)
  list(
    marker_ids = colnames(z_train),
    effects = effects,
    intercept = intercept,
    n_iter = as.integer(n_iter), burn_in = as.integer(burn_in),
    thin = as.integer(thin), prob_in = as.numeric(prob_in),
    seed = as.integer(seed), bglr_version = as.character(utils::packageVersion("BGLR")),
    fitted_on_training_only = TRUE
  )
}

predict_bayesc <- function(model, z_test) {
  if (!is.matrix(z_test) || !is.numeric(z_test) ||
      !identical(colnames(z_test), model$marker_ids) ||
      any(!is.finite(z_test))) {
    stop("BayesC prediction matrix must have the fitted markers in order.",
         call. = FALSE)
  }
  predicted <- as.numeric(model$intercept + z_test %*% model$effects)
  if (any(!is.finite(predicted))) {
    stop("BayesC produced non-finite predictions.", call. = FALSE)
  }
  predicted
}
