`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

assert_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be one non-missing logical value.", call. = FALSE)
  }
  invisible(TRUE)
}

assert_unique_nonmissing <- function(x, name) {
  if (length(x) == 0L || anyNA(x) || any(!nzchar(as.character(x)))) {
    stop(name, " must be non-empty and contain no missing/blank values.", call. = FALSE)
  }
  if (anyDuplicated(x)) {
    stop(name, " contains duplicate values.", call. = FALSE)
  }
  invisible(TRUE)
}

with_seed_preserved <- function(seed, code) {
  if (length(seed) != 1L || is.na(seed)) {
    stop("seed must be one non-missing value.", call. = FALSE)
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

