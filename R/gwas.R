selected_gwas_pc_count <- function(pc_selection, trait) {
  if (!is.list(pc_selection) || !is.data.frame(pc_selection$selected_by_trait)) {
    stop("pc_selection must contain selected_by_trait.", call. = FALSE)
  }
  row <- pc_selection$selected_by_trait[
    pc_selection$selected_by_trait$trait == trait,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) {
    stop("Exactly one frozen GWAS PC count is required for trait: ", trait, call. = FALSE)
  }
  n_pcs <- as.integer(row$selected_n_pcs[[1L]])
  if (is.na(n_pcs) || n_pcs < 0L) stop("The frozen GWAS PC count is invalid.", call. = FALSE)
  n_pcs
}

prepare_gapit_numeric_inputs <- function(dataset, structure, pc_selection, trait) {
  validate_genomic_dataset(dataset)
  require_marker_map(dataset, "GWAS")
  trait <- as.character(trait)
  if (length(trait) != 1L || is.na(trait) || !nzchar(trait)) {
    stop("trait must be one non-empty name.", call. = FALSE)
  }
  phenotype <- dataset$phenotypes[
    dataset$phenotypes$trait == trait & is.finite(dataset$phenotypes$value),
    c("sample_id", "value"),
    drop = FALSE
  ]
  if (nrow(phenotype) < 5L) stop("Too few finite observations for GWAS.", call. = FALSE)
  if (anyDuplicated(phenotype$sample_id)) {
    stop("GWAS phenotype sample IDs must be unique within a trait.", call. = FALSE)
  }
  sample_ids <- rownames(dataset$genotype)
  sample_ids <- sample_ids[sample_ids %in% phenotype$sample_id]
  phenotype_rows <- match(sample_ids, phenotype$sample_id)
  genotype_rows <- match(sample_ids, rownames(dataset$genotype))
  if (anyNA(phenotype_rows) || anyNA(genotype_rows)) {
    stop("GWAS sample matching failed.", call. = FALSE)
  }

  genotype <- dataset$genotype[genotype_rows, , drop = FALSE]
  missing_before <- sum(is.na(genotype))
  marker_means <- colMeans(genotype, na.rm = TRUE)
  if (any(!is.finite(marker_means))) {
    stop("A GWAS marker has no callable genotype among phenotyped samples.", call. = FALSE)
  }
  missing <- which(is.na(genotype), arr.ind = TRUE)
  if (nrow(missing)) genotype[missing] <- marker_means[missing[, "col"]]
  if (anyNA(genotype) || any(genotype < 0 | genotype > 2)) {
    stop("GWAS numeric genotype must be complete and remain within [0, 2].", call. = FALSE)
  }

  n_pcs <- selected_gwas_pc_count(pc_selection, trait)
  score_rows <- match(sample_ids, structure$scores$sample_id)
  grm_rows <- match(sample_ids, rownames(structure$grm))
  if (anyNA(score_rows) || anyNA(grm_rows)) {
    stop("GWAS samples are absent from the stored structure object.", call. = FALSE)
  }
  pc_names <- if (n_pcs) paste0("PC", seq_len(n_pcs)) else character()
  if (!all(pc_names %in% names(structure$scores))) {
    stop("The stored structure object lacks a selected GWAS PC.", call. = FALSE)
  }

  y <- data.frame(Taxa = sample_ids, stringsAsFactors = FALSE)
  y[[trait]] <- as.numeric(phenotype$value[phenotype_rows])
  gd <- data.frame(
    Taxa = sample_ids,
    as.data.frame(genotype, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  gm <- data.frame(
    Name = as.character(dataset$markers$marker_id),
    Chromosome = as.character(dataset$markers$chromosome),
    Position = as.numeric(dataset$markers$position_bp),
    stringsAsFactors = FALSE
  )
  cv <- NULL
  if (n_pcs > 0L) {
    cv <- data.frame(
      Taxa = sample_ids,
      structure$scores[score_rows, pc_names, drop = FALSE],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  grm <- structure$grm[grm_rows, grm_rows, drop = FALSE]
  dimnames(grm) <- list(sample_ids, sample_ids)
  ki <- data.frame(
    Taxa = sample_ids,
    as.data.frame(grm, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  sample_rows <- match(sample_ids, dataset$samples$sample_id)
  breeding_material <- dataset$samples[
    sample_rows,
    c("sample_id", "breeding_material_type", "breeding_material_type_source"),
    drop = FALSE
  ]

  output <- list(
    trait = trait,
    Y = y,
    GD = gd,
    GM = gm,
    KI = ki,
    CV = cv,
    breeding_material = breeding_material,
    audit = data.frame(
      trait = trait,
      n_observations = length(sample_ids),
      n_markers = ncol(genotype),
      selected_n_pcs = n_pcs,
      genotype_missing_before_imputation = missing_before,
      genotype_missing_after_imputation = sum(is.na(genotype)),
      imputation = "trait-sample marker mean; genotype-only",
      stringsAsFactors = FALSE
    ),
    provenance = list(
      interface = "GAPIT 4.1 numeric Y/GD/GM/KI/CV",
      sample_matching = "explicit sample_id intersection, returned in genotype order",
      marker_matching = "GM Name exactly equals GD marker columns",
      pc_count_frozen_before_marker_tests = TRUE,
      observed_gwas_inflation_used_for_pc_selection = FALSE,
      genotype_imputation_used_phenotype_values = FALSE
    )
  )
  validate_gapit_numeric_inputs(output)
  output
}

validate_gapit_numeric_inputs <- function(inputs) {
  required <- c("trait", "Y", "GD", "GM", "KI", "CV", "breeding_material", "audit")
  if (!is.list(inputs) || !all(required %in% names(inputs))) {
    stop("GAPIT input bundle is incomplete.", call. = FALSE)
  }
  if (!is.data.frame(inputs$Y) || ncol(inputs$Y) != 2L) {
    stop("Y must contain Taxa and exactly one trait.", call. = FALSE)
  }
  sample_ids <- as.character(inputs$Y[[1L]])
  assert_unique_nonmissing(sample_ids, "GAPIT phenotype Taxa")
  if (any(!is.finite(as.numeric(inputs$Y[[2L]])))) {
    stop("GAPIT phenotype values must be finite.", call. = FALSE)
  }
  if (!is.data.frame(inputs$GD) || ncol(inputs$GD) < 2L ||
      !identical(as.character(inputs$GD[[1L]]), sample_ids)) {
    stop("GD Taxa must exactly match Y Taxa and order.", call. = FALSE)
  }
  marker_ids <- colnames(inputs$GD)[-1L]
  assert_unique_nonmissing(marker_ids, "GAPIT numeric marker IDs")
  dosage <- as.matrix(inputs$GD[, -1L, drop = FALSE])
  storage.mode(dosage) <- "double"
  if (anyNA(dosage) || any(!is.finite(dosage)) || any(dosage < 0 | dosage > 2)) {
    stop("GD marker dosage must be complete, finite, and in [0, 2].", call. = FALSE)
  }
  if (!is.data.frame(inputs$GM) || ncol(inputs$GM) < 3L ||
      !identical(as.character(inputs$GM[[1L]]), marker_ids)) {
    stop("GM marker IDs must exactly match GD marker columns and order.", call. = FALSE)
  }
  if (anyNA(inputs$GM[[2L]]) || any(!is.finite(as.numeric(inputs$GM[[3L]]))) ||
      any(as.numeric(inputs$GM[[3L]]) <= 0)) {
    stop("GM chromosome and position fields must be valid.", call. = FALSE)
  }
  if (!is.data.frame(inputs$KI) || nrow(inputs$KI) != length(sample_ids) ||
      ncol(inputs$KI) != length(sample_ids) + 1L ||
      !identical(as.character(inputs$KI[[1L]]), sample_ids) ||
      !identical(colnames(inputs$KI)[-1L], sample_ids)) {
    stop("KI must contain an exactly aligned Taxa column and square relationship matrix.",
         call. = FALSE)
  }
  kinship <- as.matrix(inputs$KI[, -1L, drop = FALSE])
  storage.mode(kinship) <- "double"
  if (anyNA(kinship) || max(abs(kinship - t(kinship))) > 1e-8) {
    stop("KI must be complete and symmetric.", call. = FALSE)
  }
  if (!is.null(inputs$CV)) {
    if (!is.data.frame(inputs$CV) || ncol(inputs$CV) < 2L ||
        !identical(as.character(inputs$CV[[1L]]), sample_ids) ||
        anyNA(inputs$CV[, -1L, drop = FALSE])) {
      stop("CV must contain aligned Taxa and complete selected PCs.", call. = FALSE)
    }
  }
  if (!is.data.frame(inputs$breeding_material) ||
      !identical(as.character(inputs$breeding_material$sample_id), sample_ids)) {
    stop("Breeding-material metadata must align exactly to GWAS samples.", call. = FALSE)
  }
  invisible(TRUE)
}

gapit_backend_status <- function(minimum_version = "4.1.0") {
  installed <- requireNamespace("GAPIT", quietly = TRUE)
  version <- if (installed) as.character(utils::packageVersion("GAPIT")) else NA_character_
  exported <- installed && "GAPIT" %in% getNamespaceExports("GAPIT")
  compatible <- installed && exported &&
    utils::compareVersion(version, as.character(minimum_version)) >= 0L
  data.frame(
    backend = "GAPIT",
    installed = installed,
    version = version,
    minimum_version = as.character(minimum_version),
    gapit_exported = exported,
    compatible = compatible,
    stringsAsFactors = FALSE
  )
}

neutralize_external_command_calls <- function(expression) {
  lines <- deparse(expression, width.cutoff = 500L)
  external <- grepl(
    "^[[:space:]]*(base::)?(system|system2|shell|shell[.]exec)[(]",
    lines
  )
  if (!any(external)) return(list(expression = expression, count = 0L))
  if (any(!grepl("[)][[:space:]]*$", lines[external]))) {
    stop("A GAPIT external-command call spans multiple deparsed lines; refusing to patch.",
         call. = FALSE)
  }
  lines[external] <- "invisible(0L)"
  parsed <- parse(text = paste(lines, collapse = "\n"))
  if (length(parsed) != 1L) stop("Failed to parse guarded GAPIT function body.", call. = FALSE)
  list(expression = parsed[[1L]], count = sum(external))
}

with_gapit_external_command_guard <- function(code) {
  if (!requireNamespace("GAPIT", quietly = TRUE)) {
    stop("GAPIT must be installed before applying its external-command guard.", call. = FALSE)
  }
  namespace <- asNamespace("GAPIT")
  target_names <- intersect(c("GAPIT.Bus", "GAPIT.GS", "GAPIT.Main"), ls(namespace, all.names = TRUE))
  if (length(target_names) != 3L) {
    stop("The GAPIT 4.1 namespace does not contain the expected guarded functions.", call. = FALSE)
  }
  originals <- lapply(target_names, get, envir = namespace, inherits = FALSE)
  names(originals) <- target_names
  restore <- function() {
    for (name in target_names) {
      was_locked <- bindingIsLocked(name, namespace)
      if (was_locked) unlockBinding(name, namespace)
      assign(name, originals[[name]], envir = namespace)
      if (was_locked) lockBinding(name, namespace)
    }
  }
  on.exit(restore(), add = TRUE)
  neutralized_count <- 0L
  for (name in target_names) {
    replacement <- originals[[name]]
    guarded <- neutralize_external_command_calls(body(replacement))
    body(replacement) <- guarded$expression
    neutralized_count <- neutralized_count + guarded$count
    environment(replacement) <- namespace
    was_locked <- bindingIsLocked(name, namespace)
    if (was_locked) unlockBinding(name, namespace)
    assign(name, replacement, envir = namespace)
    if (was_locked) lockBinding(name, namespace)
  }
  value <- eval(substitute(code), envir = parent.frame())
  attr(value, "gapit_external_command_calls_neutralized") <- neutralized_count
  value
}

normalize_gwas_results <- function(raw, marker_map, method, engine) {
  if (is.matrix(raw)) raw <- as.data.frame(raw, stringsAsFactors = FALSE)
  if (!is.data.frame(raw) || !nrow(raw)) stop("GWAS engine returned no marker results.", call. = FALSE)
  normalized_names <- gsub("[^a-z0-9]", "", tolower(names(raw)))
  find_column <- function(aliases, required = TRUE) {
    hit <- which(normalized_names %in% aliases)
    if (!length(hit)) {
      if (required) stop("GWAS result lacks a required column: ", aliases[[1L]], call. = FALSE)
      return(NA_integer_)
    }
    hit[[1L]]
  }
  marker_column <- find_column(c("snp", "marker", "markerid", "name"))
  p_column <- find_column(c("pvalue", "p", "pvalues"))
  effect_column <- find_column(c("effect", "snpeffect", "beta"), required = FALSE)
  se_column <- find_column(c("se", "stderr", "standarderror"), required = FALSE)
  maf_column <- find_column(c("maf", "minorallelefrequency"), required = FALSE)
  nobs_column <- find_column(c("nobs", "n"), required = FALSE)
  marker_ids <- as.character(raw[[marker_column]])
  assert_unique_nonmissing(marker_ids, "GWAS result marker IDs")
  if (!setequal(marker_ids, marker_map$Name)) {
    stop("GWAS result marker IDs do not exactly match the tested marker map.", call. = FALSE)
  }
  rows <- match(marker_map$Name, marker_ids)
  p_value <- suppressWarnings(as.numeric(raw[[p_column]][rows]))
  if (any(!is.na(p_value) & (!is.finite(p_value) | p_value < 0 | p_value > 1))) {
    stop("GWAS p-values must lie in [0, 1] or be missing.", call. = FALSE)
  }
  optional_numeric <- function(column) {
    if (is.na(column)) rep(NA_real_, nrow(marker_map)) else
      suppressWarnings(as.numeric(raw[[column]][rows]))
  }
  data.frame(
    method = as.character(method),
    engine = as.character(engine),
    marker_id = marker_map$Name,
    chromosome = marker_map$Chromosome,
    position_bp = as.numeric(marker_map$Position),
    p_value = p_value,
    effect = optional_numeric(effect_column),
    standard_error = optional_numeric(se_column),
    maf = optional_numeric(maf_column),
    n_observations = optional_numeric(nobs_column),
    stringsAsFactors = FALSE
  )
}

add_gwas_diagnostics <- function(results, alpha = 0.05, fdr = 0.05) {
  finite <- is.finite(results$p_value)
  n_tested <- sum(finite)
  if (!n_tested) stop("No finite GWAS p-values are available.", call. = FALSE)
  bonferroni <- as.numeric(alpha) / n_tested
  adjusted <- rep(NA_real_, nrow(results))
  adjusted[finite] <- stats::p.adjust(results$p_value[finite], method = "BH")
  chisq <- stats::qchisq(
    1 - pmin(pmax(results$p_value[finite], .Machine$double.xmin), 1 - .Machine$double.eps),
    df = 1
  )
  lambda_gc <- stats::median(chisq) / stats::qchisq(0.5, df = 1)
  results$bh_fdr <- adjusted
  results$bonferroni_threshold <- bonferroni
  results$significant_bonferroni <- finite & results$p_value <= bonferroni
  results$significant_fdr <- finite & adjusted <= as.numeric(fdr)
  list(
    results = results,
    diagnostics = list(
      n_tested = n_tested,
      bonferroni_threshold = bonferroni,
      fdr_threshold = as.numeric(fdr),
      lambda_gc = lambda_gc,
      n_bonferroni = sum(results$significant_bonferroni),
      n_fdr = sum(results$significant_fdr)
    )
  )
}

run_gapit_models <- function(
    inputs,
    models = c("MLM", "FarmCPU", "BLINK"),
    minimum_version = "4.1.0",
    alpha = 0.05,
    fdr = 0.05,
    output_directory = file.path("results", "models", "gapit"),
    engine = NULL) {
  validate_gapit_numeric_inputs(inputs)
  canonical <- c(MLM = "MLM", FARMCPU = "FarmCPU", BLINK = "BLINK")
  keys <- toupper(as.character(models))
  if (!length(keys) || any(!keys %in% names(canonical))) {
    stop("GAPIT models must be selected from MLM, FarmCPU, and BLINK.", call. = FALSE)
  }
  models <- unname(canonical[keys])
  status <- gapit_backend_status(minimum_version)
  injected <- !is.null(engine)
  if (!injected) {
    if (!isTRUE(status$compatible[[1L]])) {
      stop(
        "Compatible GAPIT is unavailable. Install project-pinned GAPIT >= ",
        minimum_version,
        " before running the official backend.",
        call. = FALSE
      )
    }
    engine <- getExportedValue("GAPIT", "GAPIT")
  }
  if (!is.function(engine)) stop("engine must be a function.", call. = FALSE)
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  old_directory <- getwd()
  on.exit(setwd(old_directory), add = TRUE)
  setwd(output_directory)

  result_rows <- list()
  summary_rows <- list()
  for (index in seq_along(models)) {
    model <- models[[index]]
    recorded_warnings <- character()
    elapsed <- system.time({
      execute_engine <- function() engine(
          Y = inputs$Y,
          GD = inputs$GD,
          GM = inputs$GM,
          KI = if (model == "MLM") inputs$KI else NULL,
          CV = inputs$CV,
          model = model,
          PCA.total = 0,
          file.output = FALSE,
          Multiple_analysis = FALSE,
          Geno.View.output = FALSE,
          PCA.View.output = FALSE,
          Phenotype.View = FALSE,
          Inter.Plot = FALSE,
          Random.model = FALSE
        )
      fit <- withCallingHandlers(
        if (injected) execute_engine() else
          with_gapit_external_command_guard(execute_engine()),
        warning = function(condition) {
          recorded_warnings <<- c(recorded_warnings, conditionMessage(condition))
          invokeRestart("muffleWarning")
        }
      )
    })[["elapsed"]]
    raw <- if (is.list(fit) && !is.null(fit$GWAS)) fit$GWAS else fit
    normalized <- normalize_gwas_results(
      raw,
      inputs$GM,
      model,
      if (injected) "injected_test_double" else paste0("GAPIT_", status$version[[1L]])
    )
    diagnostic <- add_gwas_diagnostics(normalized, alpha, fdr)
    diagnostic$results$trait <- inputs$trait
    diagnostic$results$selected_n_pcs <- inputs$audit$selected_n_pcs[[1L]]
    result_rows[[index]] <- diagnostic$results
    summary_rows[[index]] <- data.frame(
      trait = inputs$trait,
      method = model,
      engine = unique(diagnostic$results$engine),
      selected_n_pcs = inputs$audit$selected_n_pcs[[1L]],
      n_observations = nrow(inputs$Y),
      n_markers = diagnostic$diagnostics$n_tested,
      lambda_gc = diagnostic$diagnostics$lambda_gc,
      bonferroni_threshold = diagnostic$diagnostics$bonferroni_threshold,
      n_bonferroni = diagnostic$diagnostics$n_bonferroni,
      n_fdr = diagnostic$diagnostics$n_fdr,
      elapsed_seconds = elapsed,
      warning_count = length(recorded_warnings),
      warnings = paste(unique(recorded_warnings), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }
  list(
    results = do.call(rbind, result_rows),
    summary = do.call(rbind, summary_rows),
    backend_status = status,
    provenance = list(
      official_gapit_executed = !injected,
      injected_test_double = injected,
      external_commands_neutralized = !injected,
      exact_pc_count_supplied_as_external_CV = TRUE,
      GAPIT_internal_PCA_disabled = TRUE
    )
  )
}

run_reference_p3d_mlm <- function(inputs) {
  validate_gapit_numeric_inputs(inputs)
  y <- as.numeric(inputs$Y[[2L]])
  genotype <- as.matrix(inputs$GD[, -1L, drop = FALSE])
  storage.mode(genotype) <- "double"
  relationship <- as.matrix(inputs$KI[, -1L, drop = FALSE])
  storage.mode(relationship) <- "double"
  fixed <- matrix(1, nrow = length(y), ncol = 1L)
  colnames(fixed) <- "intercept"
  if (!is.null(inputs$CV)) {
    fixed <- cbind(fixed, as.matrix(inputs$CV[, -1L, drop = FALSE]))
  }
  null_fit <- profile_null_mixed_model(y, fixed, relationship)
  decomposition <- eigen((relationship + t(relationship)) / 2, symmetric = TRUE)
  tolerance <- max(1, max(abs(decomposition$values))) * 1e-8
  if (any(decomposition$values < -tolerance)) {
    stop("The GWAS relationship matrix is not positive semidefinite.", call. = FALSE)
  }
  eigenvalues <- pmax(decomposition$values, 0)
  inverse_root <- 1 / sqrt(eigenvalues + null_fit$delta)
  whiten <- function(value) {
    sweep(crossprod(decomposition$vectors, value), 1L, inverse_root, FUN = "*")
  }
  whitened_fixed <- whiten(fixed)
  whitened_y <- as.numeric(whiten(matrix(y, ncol = 1L)))
  whitened_genotype <- whiten(genotype)
  fixed_qr <- qr(whitened_fixed, tol = 1e-10)
  if (fixed_qr$rank != ncol(fixed)) stop("GWAS fixed effects are rank deficient.", call. = FALSE)
  q <- qr.Q(fixed_qr)
  residual_y <- whitened_y - as.numeric(q %*% crossprod(q, whitened_y))
  residual_genotype <- whitened_genotype - q %*% crossprod(q, whitened_genotype)
  denominator <- colSums(residual_genotype^2)
  numerator <- as.numeric(crossprod(residual_y, residual_genotype))
  residual_df <- length(y) - fixed_qr$rank - 1L
  if (residual_df < 1L) stop("Too few residual degrees of freedom for GWAS.", call. = FALSE)
  effect <- numerator / denominator
  null_sse <- sum(residual_y^2)
  marker_sse <- pmax(null_sse - numerator^2 / denominator, 0)
  standard_error <- sqrt((marker_sse / residual_df) / denominator)
  statistic <- effect / standard_error
  p_value <- 2 * stats::pt(abs(statistic), df = residual_df, lower.tail = FALSE)
  unusable <- !is.finite(denominator) | denominator <= 1e-12
  effect[unusable] <- NA_real_
  standard_error[unusable] <- NA_real_
  p_value[unusable] <- NA_real_
  maf <- pmin(colMeans(genotype) / 2, 1 - colMeans(genotype) / 2)
  list(
    GWAS = data.frame(
      SNP = inputs$GM$Name,
      Chromosome = inputs$GM$Chromosome,
      Position = inputs$GM$Position,
      P.value = p_value,
      effect = effect,
      SE = standard_error,
      maf = maf,
      nobs = length(y),
      stringsAsFactors = FALSE
    ),
    delta = null_fit$delta,
    residual_df = residual_df,
    null_model = null_fit
  )
}

run_reference_gwas_development <- function(
    dataset,
    structure,
    pc_selection,
    gwas_config = list()) {
  traits <- sort(unique(as.character(dataset$phenotypes$trait)))
  alpha <- as.numeric(gwas_config$alpha %||% 0.05)
  fdr <- as.numeric(gwas_config$fdr %||% 0.05)
  minimum_version <- as.character(gwas_config$gapit_minimum_version %||% "4.1.0")
  result_rows <- list()
  summary_rows <- list()
  audit_rows <- list()
  for (index in seq_along(traits)) {
    inputs <- prepare_gapit_numeric_inputs(dataset, structure, pc_selection, traits[[index]])
    elapsed <- system.time(fit <- run_reference_p3d_mlm(inputs))[["elapsed"]]
    normalized <- normalize_gwas_results(
      fit$GWAS,
      inputs$GM,
      method = "MLM",
      engine = "internal_P3D_reference_v1"
    )
    diagnostic <- add_gwas_diagnostics(normalized, alpha, fdr)
    diagnostic$results$trait <- inputs$trait
    diagnostic$results$selected_n_pcs <- inputs$audit$selected_n_pcs[[1L]]
    result_rows[[index]] <- diagnostic$results
    material <- unique(inputs$breeding_material$breeding_material_type)
    summary_rows[[index]] <- data.frame(
      trait = inputs$trait,
      method = "MLM",
      engine = "internal_P3D_reference_v1",
      selected_n_pcs = inputs$audit$selected_n_pcs[[1L]],
      n_observations = nrow(inputs$Y),
      n_markers = diagnostic$diagnostics$n_tested,
      delta = fit$delta,
      lambda_gc = diagnostic$diagnostics$lambda_gc,
      bonferroni_threshold = diagnostic$diagnostics$bonferroni_threshold,
      n_bonferroni = diagnostic$diagnostics$n_bonferroni,
      n_fdr = diagnostic$diagnostics$n_fdr,
      elapsed_seconds = elapsed,
      breeding_material_type = paste(sort(material), collapse = ";"),
      stringsAsFactors = FALSE
    )
    audit_rows[[index]] <- inputs$audit
  }
  list(
    results = do.call(rbind, result_rows),
    summary = do.call(rbind, summary_rows),
    input_audit = do.call(rbind, audit_rows),
    backend_status = gapit_backend_status(minimum_version),
    provenance = list(
      development_engine = "internal P3D-style MLM reference implementation",
      official_gapit_executed = FALSE,
      purpose = "end-to-end adapter, identity, diagnostic, and plotting validation",
      primary_production_models = c("GAPIT MLM", "GAPIT FarmCPU", "GAPIT BLINK"),
      variance_ratio_estimated_once_under_null = TRUE,
      pc_count_frozen_before_marker_tests = TRUE,
      observed_gwas_inflation_used_for_pc_selection = FALSE,
      locked_simulation_tasks_analyzed = FALSE
    )
  )
}

run_gapit_development <- function(dataset, structure, pc_selection, gwas_config) {
  traits <- sort(unique(as.character(dataset$phenotypes$trait)))
  models <- unlist(gwas_config$primary_models %||% c("MLM", "FarmCPU", "BLINK"))
  minimum_version <- as.character(gwas_config$gapit_minimum_version %||% "4.1.0")
  status <- gapit_backend_status(minimum_version)
  if (!isTRUE(status$compatible[[1L]])) {
    stop("The official GAPIT development target requires compatible GAPIT.", call. = FALSE)
  }
  pieces <- vector("list", length(traits))
  audits <- vector("list", length(traits))
  for (index in seq_along(traits)) {
    inputs <- prepare_gapit_numeric_inputs(dataset, structure, pc_selection, traits[[index]])
    trait_directory <- file.path(
      "results", "models", "gapit_development", gsub("[^A-Za-z0-9_.-]", "_", traits[[index]])
    )
    pieces[[index]] <- run_gapit_models(
      inputs,
      models = models,
      minimum_version = minimum_version,
      alpha = gwas_config$alpha %||% 0.05,
      fdr = gwas_config$fdr %||% 0.05,
      output_directory = trait_directory
    )
    material <- unique(inputs$breeding_material$breeding_material_type)
    pieces[[index]]$summary$breeding_material_type <- paste(sort(material), collapse = ";")
    audits[[index]] <- inputs$audit
  }
  list(
    results = do.call(rbind, lapply(pieces, `[[`, "results")),
    summary = do.call(rbind, lapply(pieces, `[[`, "summary")),
    input_audit = do.call(rbind, audits),
    backend_status = status,
    provenance = list(
      official_gapit_executed = TRUE,
      external_commands_neutralized = TRUE,
      primary_models = models,
      pc_count_frozen_before_marker_tests = TRUE,
      observed_gwas_inflation_used_for_pc_selection = FALSE,
      locked_simulation_tasks_analyzed = FALSE
    )
  )
}

draw_gwas_trait_panels <- function(results, trait, method = NULL, alpha = 0.05) {
  selected <- results$trait == trait & is.finite(results$p_value)
  if (!is.null(method)) selected <- selected & results$method == method
  current <- results[selected, , drop = FALSE]
  if (!nrow(current)) stop("No finite GWAS results for trait: ", trait, call. = FALSE)
  label <- if (is.null(method)) trait else paste(trait, method, sep = " — ")
  chromosome <- as.character(current$chromosome)
  chromosome_numeric <- suppressWarnings(as.numeric(chromosome))
  chromosome_levels <- unique(chromosome[order(
    is.na(chromosome_numeric), chromosome_numeric, chromosome
  )])
  x <- numeric(nrow(current))
  ticks <- numeric(length(chromosome_levels))
  offset <- 0
  for (index in seq_along(chromosome_levels)) {
    rows <- which(chromosome == chromosome_levels[[index]])
    local <- current$position_bp[rows] - min(current$position_bp[rows])
    x[rows] <- local + offset
    ticks[[index]] <- mean(range(x[rows]))
    offset <- max(x[rows]) + max(1e6, 0.02 * max(current$position_bp))
  }
  y <- -log10(pmax(current$p_value, .Machine$double.xmin))
  colors <- rep(c("#2F5597", "#70AD47"), length.out = length(chromosome_levels))
  point_colors <- colors[match(chromosome, chromosome_levels)]
  graphics::plot(
    x, y,
    pch = 16, cex = 0.55, col = point_colors,
    xaxt = "n", xlab = "Chromosome", ylab = expression(-log[10](italic(p))),
    main = paste(label, "Manhattan")
  )
  graphics::axis(1, at = ticks, labels = chromosome_levels, cex.axis = 0.8)
  graphics::abline(
    h = -log10(as.numeric(alpha) / nrow(current)),
    col = "#C00000", lty = 2, lwd = 1.2
  )
  graphics::box()

  observed <- sort(current$p_value)
  expected <- stats::ppoints(length(observed))
  observed_log <- -log10(pmax(observed, .Machine$double.xmin))
  expected_log <- -log10(expected)
  limit <- max(c(observed_log, expected_log))
  graphics::plot(
    expected_log, observed_log,
    pch = 16, cex = 0.6, col = "#2F5597",
    xlim = c(0, limit), ylim = c(0, limit),
    xlab = expression(Expected~~-log[10](italic(p))),
    ylab = expression(Observed~~-log[10](italic(p))),
    main = paste(label, "QQ")
  )
  graphics::abline(0, 1, col = "#C00000", lty = 2)
  chisq <- stats::qchisq(
    1 - pmin(pmax(current$p_value, .Machine$double.xmin), 1 - .Machine$double.eps),
    df = 1
  )
  lambda <- stats::median(chisq) / stats::qchisq(0.5, df = 1)
  graphics::legend(
    "topleft", legend = bquote(lambda[GC] == .(format(round(lambda, 3), nsmall = 3))),
    bty = "n", cex = 0.85
  )
  invisible(TRUE)
}

write_gwas_diagnostic_figures <- function(
    results,
    pdf_path = file.path("results", "figures", "tassel_tutorial_reference_gwas.pdf"),
    png_path = file.path("results", "figures", "tassel_tutorial_reference_gwas.png")) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  if (!is.null(png_path)) {
    dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  }
  groups <- unique(results[, c("trait", "method"), drop = FALSE])
  grDevices::pdf(pdf_path, width = 11, height = 4.5, onefile = TRUE)
  for (index in seq_len(nrow(groups))) {
    graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.2, 1.2))
    draw_gwas_trait_panels(results, groups$trait[[index]], groups$method[[index]])
  }
  grDevices::dev.off()
  if (is.null(png_path)) return(pdf_path)
  grDevices::png(
    png_path,
    width = 2200,
    height = max(900, 800 * nrow(groups)),
    res = 180
  )
  graphics::par(mfrow = c(nrow(groups), 2), mar = c(4.2, 4.2, 3.2, 1.2))
  for (index in seq_len(nrow(groups))) {
    draw_gwas_trait_panels(results, groups$trait[[index]], groups$method[[index]])
  }
  grDevices::dev.off()
  c(pdf_path, png_path)
}

write_reference_gwas_outputs <- function(
    gwas,
    processed_path = file.path("data", "processed", "tassel_tutorial", "reference_gwas.rds"),
    table_directory = file.path("results", "tables"),
    figure_directory = file.path("results", "figures")) {
  dir.create(dirname(processed_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    rds = processed_path,
    results = file.path(table_directory, "tassel_tutorial_reference_gwas_results.tsv"),
    summary = file.path(table_directory, "tassel_tutorial_reference_gwas_summary.tsv"),
    input_audit = file.path(table_directory, "tassel_tutorial_gwas_input_audit.tsv")
  )
  saveRDS(gwas, paths[["rds"]], version = 3)
  tables <- list(
    results = gwas$results,
    summary = gwas$summary,
    input_audit = gwas$input_audit
  )
  for (name in names(tables)) {
    utils::write.table(
      tables[[name]], paths[[name]], sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  figures <- write_gwas_diagnostic_figures(
    gwas$results,
    pdf_path = file.path(figure_directory, "tassel_tutorial_reference_gwas.pdf"),
    png_path = file.path(figure_directory, "tassel_tutorial_reference_gwas.png")
  )
  unname(c(paths, figures))
}

write_gapit_development_outputs <- function(
    gwas,
    processed_path = file.path("data", "processed", "tassel_tutorial", "gapit_development.rds"),
    table_directory = file.path("results", "tables"),
    figure_directory = file.path("results", "figures")) {
  dir.create(dirname(processed_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    rds = processed_path,
    results = file.path(table_directory, "tassel_tutorial_gapit_development_results.tsv"),
    summary = file.path(table_directory, "tassel_tutorial_gapit_development_summary.tsv"),
    input_audit = file.path(table_directory, "tassel_tutorial_gapit_input_audit.tsv"),
    backend_status = file.path(table_directory, "gwas_backend_status.tsv")
  )
  saveRDS(gwas, paths[["rds"]], version = 3)
  tables <- list(
    results = gwas$results,
    summary = gwas$summary,
    input_audit = gwas$input_audit,
    backend_status = gwas$backend_status
  )
  for (name in names(tables)) {
    utils::write.table(
      tables[[name]], paths[[name]], sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  figure <- write_gwas_diagnostic_figures(
    gwas$results,
    pdf_path = file.path(figure_directory, "tassel_tutorial_gapit_development.pdf"),
    png_path = NULL
  )
  unname(c(paths, figure))
}
