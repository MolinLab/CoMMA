breeding_followup_roles <- function(design, discovery_ids, metadata) {
  if (is.character(design)) design <- breeding_delimited_table(design)
  fields <- c("sample_id", "role", "shuffle_stratum")
  if (!is.data.frame(design) || !setequal(names(design), fields) || anyDuplicated(names(design))) {
    stop("Supply sample_id, role and shuffle_stratum in the follow-up design.")
  }
  design <- as.data.frame(lapply(design[fields], as.character))
  assert_unique_nonmissing(design$sample_id, "Follow-up sample IDs")
  if (!setequal(unique(design$role), c("confirmation_fit", "explanation")) ||
      anyNA(design$shuffle_stratum[design$role == "confirmation_fit"]) ||
      any(!nzchar(trimws(design$shuffle_stratum[design$role == "confirmation_fit"])))) {
    stop("Supply both follow-up roles and an explicit shuffle stratum for every confirmation-fit line.")
  }
  if (length(intersect(discovery_ids, design$sample_id))) stop("Discovery and follow-up sample IDs overlap.")
  design <- design[order(design$sample_id, method = "radix"), , drop = FALSE]
  roles <- rbind(data.frame(sample_id = discovery_ids, role = "discovery", shuffle_stratum = NA_character_), design)
  if ("fingerprint_group" %in% names(metadata)) {
    group <- metadata$fingerprint_group[match(roles$sample_id, metadata$sample_id)]
    if (anyNA(group) || any(!nzchar(group))) stop("Supply a fingerprint group for each line when using fingerprint groups.")
    if (anyDuplicated(unique(data.frame(group, role = roles$role))$group)) {
      stop("A genotype fingerprint group spans more than one follow-up role.")
    }
  }
  roles
}

breeding_followup_response <- function(input, trait, allowed, required, label) {
  y <- breeding_input_phenotypes(input, trait)
  if (!is.data.frame(y) || !all(c("sample_id", "value") %in% names(y))) stop("Invalid ", label, " phenotype table.")
  if ("trait" %in% names(y)) y <- y[which(y$trait == trait), , drop = FALSE]
  if (!all(y$sample_id %in% allowed)) stop("The ", label, " phenotype table contains IDs outside its assigned role.")
  gbs_training_response(y, required)
}

breeding_prepare_conditional_inputs <- function(assessment, phenotypes, sample_metadata, trait, output,
    shuffle_columns, explanation_fraction = .25, seed = 20261112L) {
  output_path <- gbs_project_path(output, must_exist = FALSE)
  if (file.exists(output_path)) stop("Follow-up inputs already exist. Preserve them and choose a new directory.")
  receipt <- readRDS(gbs_project_path(file.path(assessment, "assessment_completion.rds")))
  discovery <- receipt$training_ids
  holdout <- receipt$assessment_ids
  if (!is.character(discovery) || !is.character(holdout) || length(intersect(discovery, holdout))) {
    stop("The saved assessment must identify disjoint training and holdout samples.")
  }
  assert_unique_nonmissing(discovery, "Saved discovery sample IDs")
  assert_unique_nonmissing(holdout, "Saved holdout sample IDs")
  metadata <- breeding_input_metadata(sample_metadata)
  if (!is.character(shuffle_columns) || !length(shuffle_columns) || anyDuplicated(shuffle_columns) ||
      !all(shuffle_columns %in% setdiff(names(metadata), c("sample_id", "fingerprint_group", "assessment_eligible")))) {
    stop("Choose explicit shuffle-stratum columns from the sample metadata.")
  }
  panel <- breeding_metadata(metadata, holdout)
  strata <- lapply(panel[shuffle_columns], as.character)
  if (any(vapply(strata, function(x) anyNA(x) || any(!nzchar(trimws(x))), logical(1L)))) {
    stop("Every registered holdout line needs a value in the chosen shuffle-stratum columns.")
  }
  labels <- if (length(strata) == 1L) strata[[1L]] else
    do.call(paste, c(unname(lapply(strata, encodeString, quote = '"')), list(sep = " / ")))
  split <- breeding_assessment_split(panel, explanation_fraction, seed)
  design <- data.frame(sample_id = split$sample_id,
    role = ifelse(split$role == "assessment", "explanation", "confirmation_fit"),
    shuffle_stratum = labels[match(split$sample_id, panel$sample_id)])
  breeding_followup_roles(design, discovery, metadata)
  confirmation <- design$sample_id[design$role == "confirmation_fit"]
  explanation <- design$sample_id[design$role == "explanation"]
  if (length(confirmation) < 30L || length(explanation) < 8L) {
    stop("Allocate at least 30 confirmation-fit and eight explanation lines before marker-level missingness.")
  }

  # Split by genotype metadata before loading phenotype values. Write only the assigned responses.
  responses <- if (is.character(phenotypes)) breeding_delimited_table(phenotypes) else phenotypes
  responses <- breeding_input_name(responses, "sample_id", "Taxa")
  if (!"sample_id" %in% names(responses)) stop("Phenotypes require a sample_id or Taxa column.")
  responses <- breeding_input_trait_table(responses[responses$sample_id %in% c(discovery, confirmation), , drop = FALSE],
                                         trait)$phenotypes
  discovery_y <- gbs_training_response(responses[responses$sample_id %in% discovery, , drop = FALSE], discovery)
  confirmation_y <- gbs_training_response(responses[responses$sample_id %in% confirmation, , drop = FALSE], confirmation)
  if (!dir.create(output_path, recursive = TRUE, showWarnings = FALSE)) stop("Could not create follow-up inputs.")
  files <- file.path(output, c("design.tsv", "discovery_phenotypes.tsv", "confirmation_phenotypes.tsv"))
  data.table::fwrite(design, files[1L], sep = "\t")
  phenotype_table <- function(y) data.frame(sample_id = names(y), trait = trait, value = unname(y))
  data.table::fwrite(phenotype_table(discovery_y), files[2L], sep = "\t")
  data.table::fwrite(phenotype_table(confirmation_y), files[3L], sep = "\t")
  counts <- as.data.frame(table(shuffle_stratum = design$shuffle_stratum, role = design$role), responseName = "lines")
  data.table::fwrite(counts, file.path(output, "role_counts.tsv"), sep = "\t")
  prepared <- list(trait = trait, seed = seed, shuffle_columns = shuffle_columns, analysis = "exploratory",
    explanation_fraction_requested = explanation_fraction, n_discovery = length(discovery),
    n_confirmation = length(confirmation), n_explanation = length(explanation), files = files,
    file_sha256 = vapply(files, digest::digest, character(1L), file = TRUE, algo = "sha256"))
  gbs_parallel_atomic(prepared, file.path(output, "preparation_completion.rds"))
  invisible(prepared)
}

breeding_followup_controls <- function(marker, map, raw) {
  i <- match(marker, map$marker_id)
  local <- map$marker_id[map$chromosome == map$chromosome[i] &
    abs(map$position_bp - map$position_bp[i]) <= 1e7 & map$marker_id != marker]
  a <- raw[, marker]
  ld <- vapply(local, function(id) {
    b <- raw[, id]
    called <- is.finite(a) & is.finite(b)
    if (sum(called) < 20L) return(NA_real_)
    suppressWarnings(stats::cor(a[called], b[called])^2)
  }, numeric(1L))
  passing <- which(ld >= .20)
  head(local[passing[order(-ld[passing], local[passing], method = "radix")]], 2L)
}

breeding_followup_scores <- function(train, query, y, blocks, permutations, seed) {
  stopifnot(identical(rownames(train), names(y)), length(blocks) == length(y),
            is.null(query) || !any(rownames(query) %in% names(y)))
  positions <- split(seq_along(blocks), match(blocks, unique(blocks)))
  index <- with_seed_preserved(seed, vapply(seq_len(permutations), function(i) {
    shuffled <- seq_along(blocks)
    for (at in positions) shuffled[at] <- at[sample.int(length(at))]
    shuffled
  }, integer(length(blocks))))
  p_value <- function(observed, null) (1 + sum(null >= observed - 1e-12 * (1 + abs(observed)))) / (length(null) + 1)
  centered <- train[, 1L] - stats::ave(train[, 1L], blocks)
  observed <- abs(sum(centered * y))
  null <- abs(drop(crossprod(centered, matrix(y[index], nrow = length(y)))))
  result <- list(association_score = observed, association_p = p_value(observed, null),
                 shap_mean_abs = NA_real_, shap_p = NA_real_, shap_fits = 0L, shap_replay_error = NA_real_,
                 shap_replay_max_fraction_of_tolerance = NA_real_)
  if (is.null(query)) return(result)

  settings <- data.frame(eta = .1, max_depth = 2L, min_child_weight = 3, subsample = 1,
                         colsample_bytree = 1, lambda = 1, alpha = 0, gamma = 0)
  params <- xgboost_parameters(settings, list(internal_threads = 1L, tree_method = "hist"), seed)
  query_matrix <- xgboost::xgb.DMatrix(query, nthread = 1L)
  max_error <- 0
  max_fraction <- 0
  fit <- function(response, shuffle) {
    context <- sprintf("marker %s, %s", colnames(train)[1L],
                       if (shuffle == 0L) "observed fit" else paste("shuffle", shuffle))
    data <- xgboost::xgb.DMatrix(train, label = unname(response), nthread = 1L)
    model <- xgboost::xgb.train(params, data, nrounds = 30L, verbose = 0L)
    phi <- stats::predict(model, query_matrix, predcontrib = TRUE, approxcontrib = FALSE)
    predicted <- as.numeric(stats::predict(model, query_matrix))
    if (!identical(dim(phi), c(nrow(query), ncol(train) + 1L)) ||
        (!is.null(colnames(phi)) && !identical(colnames(phi)[seq_len(ncol(train))], colnames(train))) ||
        length(predicted) != nrow(query) || any(!is.finite(phi)) || any(!is.finite(predicted))) {
      stop("Invalid TreeSHAP dimensions, feature order or values (", context, ").", call. = FALSE)
    }
    audit <- breeding_native_shap_audit(phi, predicted)
    fraction <- audit$absolute_residual / audit$allowed_residual
    worst <- which.max(fraction)
    if (!audit$scale_aware_pass[worst]) {
      stop(sprintf("TreeSHAP reconstruction failed (%s, query %s): error %.9g, permitted %.9g.",
        context, rownames(query)[worst], audit$absolute_residual[worst], audit$allowed_residual[worst]), call. = FALSE)
    }
    max_error <<- max(max_error, audit$absolute_residual)
    max_fraction <<- max(max_fraction, fraction[worst])
    mean(abs(phi[, 1L]))
  }
  observed <- fit(y, 0L)
  null <- vapply(seq_len(ncol(index)), function(i) fit(y[index[, i]], i), numeric(1L))
  result$shap_mean_abs <- observed
  result$shap_p <- p_value(observed, null)
  result$shap_fits <- permutations + 1L
  result$shap_replay_error <- max_error
  result$shap_replay_max_fraction_of_tolerance <- max_fraction
  result
}

breeding_followup_blocks <- function(strata, controls) {
  if (!length(strata)) return(integer())
  key <- cbind(data.frame(stratum = unname(strata)), as.data.frame(controls))
  do.call(interaction, c(unname(as.list(key)), list(drop = TRUE)))
}

breeding_followup_informative <- function(by_block) {
  vapply(by_block, function(x) length(x) >= 4L && length(unique(x)) > 1L, logical(1L))
}

breeding_diagnose_conditional_support <- function(followup, model, query, output, markers = NULL, input_memory_gib = 2) {
  destination <- gbs_project_path(output, must_exist = FALSE)
  if (file.exists(destination)) stop("Choose a new conditional support output directory.")
  receipt <- readRDS(gbs_project_path(file.path(followup, "conditional_completion.rds")))
  if (!identical(receipt$status, "conditional_followup_completed")) stop("Supply a completed conditional follow-up.")
  results <- breeding_conditional_plot_data(file.path(followup, "conditional_results.tsv"))
  roles <- breeding_delimited_table(file.path(followup, "line_roles.tsv"))
  saved <- breeding_delimited_table(file.path(followup, "shuffle_blocks.tsv"), allow_empty = TRUE)
  fitted <- breeding_load_model(model)
  sha <- function(x) digest::digest(x, algo = "sha256")
  if (!all(c("local_ld_controls", "other_discovery_lead", "n_shuffle_blocks") %in% names(results)) ||
      !all(c("sample_id", "role", "shuffle_stratum") %in% names(roles)) ||
      !all(c("marker_id", "sample_id", "shuffle_stratum", "shuffle_block", "informative_block") %in% names(saved))) {
    stop("The saved conditional tables are missing required diagnostic columns.")
  }
  if (nrow(results) != receipt$family_size ||
      !identical(sha(fitted$request), receipt$input_sha256$discovery_model_request)) {
    stop("Candidate family or saved discovery fit differs from the completed follow-up.")
  }
  assert_unique_nonmissing(roles$sample_id, "Saved follow-up sample IDs")
  confirmation <- roles$sample_id[roles$role == "confirmation_fit"]
  explanation <- roles$sample_id[roles$role == "explanation"]
  strata <- setNames(roles$shuffle_stratum, roles$sample_id)
  leads <- results$lead_marker_id
  controls <- lapply(results$local_ld_controls, function(x) if (is.na(x)) character() else strsplit(x, ",", fixed = TRUE)[[1L]])
  other <- results$other_discovery_lead
  required <- unique(c(leads, unlist(controls, use.names = FALSE), other[!is.na(other)]))
  map <- gbs_preprocessing_map(fitted$preprocessing$markers, required)
  ids <- c(confirmation, explanation)
  if (is.character(query)) {
    calls <- breeding_read_genotypes(query, markers, sample_ids = ids,
      maximum_memory_gib = input_memory_gib, reference_markers = map)$genotype
  } else {
    gbs_validate_calls(query)
    gbs_preprocessing_map(breeding_input_markers(markers), required, map)
    if (!all(ids %in% rownames(query)) || !all(required %in% colnames(query))) stop("Missing follow-up sample or marker IDs.")
    calls <- query[ids, required, drop = FALSE]
  }
  if (!identical(sha(calls), receipt$input_sha256$followup_calls) ||
      !identical(sha(map), receipt$input_sha256$marker_map)) stop("Original follow-up genotype calls or map differ.")
  if (anyDuplicated(saved[c("marker_id", "sample_id")]) || !all(saved$marker_id %in% leads)) {
    stop("Saved shuffle-block candidate or sample IDs are inconsistent.")
  }

  summaries <- vector("list", length(leads))
  block_tables <- vector("list", length(leads))
  for (i in seq_along(leads)) {
    lead <- leads[i]
    block_markers <- unique(c(controls[[i]], other[i][!is.na(other[i])]))
    raw <- calls[confirmation, c(lead, block_markers), drop = FALSE]
    lead_calls <- raw[, lead]
    n_lead_homozygous <- sum(lead_calls %in% c(0, 2))
    complete <- raw[rowSums(is.na(raw) | raw == 1) == 0L, , drop = FALSE]
    callable_ids <- rownames(complete)
    genotype <- complete[, lead]
    stratum <- strata[callable_ids]
    stages <- list(strata_only = character(), strata_plus_other_lead = other[i][!is.na(other[i])],
                   strata_plus_local_ld = controls[[i]], original_full_blocks = block_markers)
    summaries[[i]] <- do.call(rbind, lapply(names(stages), function(stage) {
      blocks <- breeding_followup_blocks(stratum, complete[, stages[[stage]], drop = FALSE])
      by_block <- split(genotype, blocks)
      informative <- breeding_followup_informative(by_block)
      eligible <- informative[as.integer(blocks)]
      if (stage == "original_full_blocks") {
        observed <- saved[which(saved$marker_id == lead), , drop = FALSE]
        observed <- observed[match(callable_ids, observed$sample_id), , drop = FALSE]
        if (sum(saved$marker_id == lead) != length(callable_ids) ||
            !identical(as.character(observed$sample_id), as.character(callable_ids)) ||
            !identical(as.integer(observed$shuffle_block), as.integer(blocks)) ||
            !identical(as.logical(observed$informative_block), unname(eligible)) ||
            !identical(as.character(observed$shuffle_stratum), as.character(stratum)) ||
            length(by_block) != as.integer(results$n_shuffle_blocks[i]) ||
            length(genotype) != results$n_confirmation[i] || sum(eligible) != results$n_informative[i]) {
          stop("Reconstructed shuffle-block evidence differs for ", lead, ".")
        }
        block_tables[[i]] <<- data.frame(marker_id = rep(lead, length(by_block)), shuffle_block = seq_along(by_block),
          n_callable = unname(lengths(by_block)), n_0 = vapply(by_block, function(x) sum(x == 0), integer(1L)),
          n_2 = vapply(by_block, function(x) sum(x == 2), integer(1L)), informative_block = unname(informative))
      }
      data.frame(priority = results$priority[i], marker_id = lead, stage = stage,
        n_registered = length(lead_calls), n_lead_missing = sum(is.na(lead_calls)),
        n_lead_heterozygous = sum(lead_calls == 1, na.rm = TRUE), n_lead_0 = sum(lead_calls == 0, na.rm = TRUE),
        n_lead_2 = sum(lead_calls == 2, na.rm = TRUE), n_excluded_for_control_calls = n_lead_homozygous - length(genotype),
        n_callable = length(genotype),
        n_callable_0 = sum(genotype == 0), n_callable_2 = sum(genotype == 2),
        n_blocks = length(by_block), n_informative_blocks = sum(informative), n_informative = sum(eligible),
        n_informative_0 = sum(genotype[eligible] == 0), n_informative_2 = sum(genotype[eligible] == 2),
        informative_fraction = if (length(genotype)) mean(eligible) else NA_real_)
    }))
  }
  if (!dir.create(destination, recursive = TRUE, showWarnings = FALSE)) stop("Could not create conditional support output.")
  summary_path <- file.path(output, "candidate_support.tsv")
  block_path <- file.path(output, "genotype_counts_by_block.tsv")
  data.table::fwrite(do.call(rbind, summaries), summary_path, sep = "\t", na = "NA")
  data.table::fwrite(do.call(rbind, block_tables), block_path, sep = "\t", na = "NA")
  invisible(list(candidate_support = summary_path, block_counts = block_path))
}

breeding_run_conditional_followup <- function(model, query, design, discovery_phenotypes,
    confirmation_phenotypes, trait, output, markers = NULL, sample_metadata = NULL, max_candidates = 10L,
    permutations = 999L, shap = FALSE, seed = 20261112L, input_memory_gib = 2) {
  models <- breeding_gwas_shortlist_methods
  breeding_check_gwas_candidates(max_candidates, models)
  whole <- function(x, minimum) is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= minimum && x < .Machine$integer.max && x == floor(x)
  if (!whole(max_candidates, 1L) || !whole(permutations, 1L) || !whole(seed, 0L) ||
      !is.logical(shap) || length(shap) != 1L || is.na(shap) ||
      !is.character(trait) || length(trait) != 1L || is.na(trait) || !nzchar(trait)) {
    stop("Supply valid candidate and permutation counts, seed, SHAP switch and one trait name.")
  }
  if (shap && !requireNamespace("xgboost", quietly = TRUE)) stop("Install XGBoost to run the optional SHAP test.")
  output_path <- gbs_project_path(output, must_exist = FALSE)
  if (file.exists(output_path)) stop("Follow-up output exists. Choose a new directory.")
  fitted <- breeding_load_model(model)
  prep <- fitted$preprocessing
  all_discovery <- prep$line_audit$sample_id
  stopifnot(all(prep$training_ids %in% all_discovery))
  metadata <- breeding_input_metadata(sample_metadata)
  roles <- breeding_followup_roles(design, all_discovery, metadata)
  confirmation <- roles$sample_id[roles$role == "confirmation_fit"]
  explanation <- roles$sample_id[roles$role == "explanation"]
  if (is.null(metadata)) metadata <- data.frame(sample_id = prep$training_ids)
  metadata <- breeding_metadata(metadata, prep$training_ids)
  discovery_y <- breeding_followup_response(discovery_phenotypes, trait, all_discovery, prep$training_ids, "discovery")
  # Reject incompatible input files before starting the discovery scans.
  if (is.character(confirmation_phenotypes)) gbs_project_path(confirmation_phenotypes)
  if (is.character(query)) gbs_project_path(query)
  work_gib <- as.double(length(confirmation)) * permutations * 20 / 1024^3
  model_gib <- as.double(utils::object.size(fitted)) / 1024^3
  gbs_parallel_budget(prep$config$resources, model_gib, .25 + work_gib, 1L)

  if (!dir.create(output_path, recursive = TRUE, showWarnings = FALSE)) stop("Could not create follow-up output.")
  started <- proc.time()[["elapsed"]]
  data.table::fwrite(roles, file.path(output, "line_roles.tsv"), sep = "\t", na = "NA")
  scans <- with_seed_preserved(seed, breeding_run_gwas(prep,
    data.frame(sample_id = names(discovery_y), value = unname(discovery_y)),
    metadata, trait, file.path(output, "discovery_GWAS"), n_pcs = prep$config$GWAS$n_pcs,
    include_source_group = prep$config$GWAS$include_source_group, models = models))
  raw <- prep$imputer$reference
  stopifnot(identical(rownames(raw), prep$training_ids), all(prep$analysis_ids %in% colnames(raw)))
  shortlist <- breeding_gwas_clump_regions(breeding_gwas_rank_regions(scans$results, prep$markers), raw, max_candidates)
  family_size <- nrow(shortlist)
  if (!family_size) stop("The discovery scans produced no eligible candidate SNPs.")
  if (1 / (permutations + 1) > .05 / family_size) stop("Increase permutations to resolve the family-wise 5% threshold.")
  leads <- shortlist$lead_marker_id
  controls <- lapply(leads, breeding_followup_controls, map = prep$markers, raw = raw)
  other_leads <- vapply(leads, function(id) head(c(setdiff(leads, id), NA_character_), 1L), character(1L))
  shortlist$local_ld_controls <- vapply(controls, paste, character(1L), collapse = ",")
  shortlist$other_discovery_lead <- other_leads
  data.table::fwrite(shortlist, file.path(output, "discovery_shortlist.tsv"), sep = "\t", na = "NA")

  required <- unique(c(leads, unlist(controls, use.names = FALSE)))
  map <- gbs_preprocessing_map(prep$markers, required)
  followup_ids <- c(confirmation, explanation)
  if (is.character(query)) {
    query <- breeding_read_genotypes(query, markers, sample_ids = followup_ids,
      maximum_memory_gib = input_memory_gib, reference_markers = map)$genotype
  } else {
    gbs_validate_calls(query)
    gbs_preprocessing_map(breeding_input_markers(markers), required, map)
    if (!all(followup_ids %in% rownames(query)) || !all(required %in% colnames(query))) {
      stop("The query matrix is missing follow-up sample or marker IDs.")
    }
    query <- query[followup_ids, required, drop = FALSE]
  }
  y <- breeding_followup_response(confirmation_phenotypes, trait, confirmation, confirmation, "confirmation-fit")
  strata <- setNames(roles$shuffle_stratum, roles$sample_id)[confirmation]
  details <- vector("list", family_size)
  results <- lapply(seq_along(leads), function(i) {
    lead <- leads[i]
    predictors <- c(lead, controls[[i]])
    block_markers <- unique(c(controls[[i]], stats::na.omit(other_leads[i])))
    training <- query[confirmation, unique(c(lead, block_markers)), drop = FALSE]
    testing <- query[explanation, predictors, drop = FALSE]
    # Use original homozygous calls, as in the audited matched benchmark.
    training <- training[rowSums(is.na(training) | training == 1) == 0L, , drop = FALSE]
    testing <- testing[rowSums(is.na(testing) | testing == 1) == 0L, , drop = FALSE]
    ids <- rownames(training)
    blocks <- breeding_followup_blocks(strata[ids], training[, block_markers, drop = FALSE])
    by_block <- split(training[, lead], blocks)
    informative <- breeding_followup_informative(by_block)
    n_informative <- sum(lengths(by_block)[informative])
    details[[i]] <<- data.frame(marker_id = rep(lead, length(ids)), sample_id = as.character(ids),
      shuffle_stratum = unname(strata[ids]), shuffle_block = as.integer(blocks),
      informative_block = unname(informative[as.integer(blocks)]))
    reason <- if (length(ids) < 30L) "too_few_callable_confirmation_lines" else
      if (n_informative < .4 * length(ids)) "insufficient_within_block_genotype_variation" else "tested"
    shap_status <- if (!shap) "disabled" else if (reason != "tested") reason else
      if (nrow(testing) < 8L) "too_few_callable_explanation_lines" else "tested"
    candidate_seed <- as.integer((as.double(seed) + 100 * i) %% .Machine$integer.max)
    scores <- if (reason == "tested") breeding_followup_scores(training[, predictors, drop = FALSE],
      if (shap_status == "tested") testing else NULL, y[ids], blocks, as.integer(permutations), candidate_seed) else
      list(association_score = NA_real_, association_p = NA_real_, shap_mean_abs = NA_real_,
           shap_p = NA_real_, shap_fits = 0L, shap_replay_error = NA_real_,
           shap_replay_max_fraction_of_tolerance = NA_real_)
    cbind(shortlist[i, , drop = FALSE], data.frame(analysis = "exploratory", n_confirmation = length(ids),
      n_explanation = nrow(testing), n_shuffle_blocks = length(by_block), n_informative = n_informative,
      association_status = reason, shap_status = shap_status, seed = candidate_seed), as.data.frame(scores))
  })
  results <- do.call(rbind, results)
  results$association_p_bonferroni <- pmin(1, results$association_p * family_size)
  results$shap_p_bonferroni <- pmin(1, results$shap_p * family_size)
  result_file <- file.path(output, "conditional_results.tsv")
  data.table::fwrite(results, result_file, sep = "\t", na = "NA")
  data.table::fwrite(do.call(rbind, details), file.path(output, "shuffle_blocks.tsv"), sep = "\t", na = "NA")
  receipt <- list(status = "conditional_followup_completed", analysis = "exploratory", trait = trait,
    results_file = result_file, max_candidates = max_candidates, family_size = family_size,
    permutations = permutations, seed = seed, shap = shap, workers = 1L, threads = 1L,
    input_sha256 = lapply(list(discovery_response = discovery_y, confirmation_response = y,
      discovery_model_request = fitted$request, followup_calls = query, line_roles = roles, marker_map = map),
      digest::digest, algo = "sha256"),
    gapit_version = as.character(utils::packageVersion("GAPIT")),
    R_version = R.version.string, platform = R.version$platform,
    xgboost_version = if (shap) as.character(utils::packageVersion("xgboost")) else NULL,
    elapsed_seconds = unname(proc.time()[["elapsed"]] - started))
  gbs_parallel_atomic(receipt, file.path(output, "conditional_completion.rds"))
  invisible(receipt)
}
