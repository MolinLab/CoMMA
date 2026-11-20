new_genomic_dataset <- function(genotype, samples, markers, phenotypes = NULL,
                                provenance = list()) {
  object <- structure(
    list(
      genotype = genotype,
      samples = samples,
      markers = markers,
      phenotypes = phenotypes,
      provenance = provenance
    ),
    class = "genomic_dataset"
  )
  validate_genomic_dataset(object)
  object
}

validate_genomic_dataset <- function(x, allow_missing_genotype = TRUE) {
  if (!inherits(x, "genomic_dataset")) {
    stop("Object must inherit from 'genomic_dataset'.", call. = FALSE)
  }
  if (!is.matrix(x$genotype) || !is.numeric(x$genotype)) {
    stop("genotype must be a numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(x$genotype)) || is.null(colnames(x$genotype))) {
    stop("genotype must have sample row names and marker column names.", call. = FALSE)
  }
  assert_unique_nonmissing(rownames(x$genotype), "genotype sample IDs")
  assert_unique_nonmissing(colnames(x$genotype), "genotype marker IDs")

  observed <- x$genotype[!is.na(x$genotype)]
  fractional_allowed <- isTRUE(x$provenance$fractional_dosage_allowed)
  outside_range <- length(observed) &&
    any(!is.finite(observed) | observed < 0 | observed > 2)
  non_integer <- length(observed) &&
    any(abs(observed - round(observed)) > 1e-8)
  if (outside_range || (!fractional_allowed && non_integer)) {
    expected <- if (fractional_allowed) "continuous dosage range [0,2]/NA" else "dosage coding 0/1/2/NA"
    stop("genotype contains values outside ", expected, ".", call. = FALSE)
  }
  if (!allow_missing_genotype && anyNA(x$genotype)) {
    stop("Missing genotype values are not allowed at this stage.", call. = FALSE)
  }

  required_sample <- c(
    "sample_id", "breeding_material_type", "breeding_material_type_source"
  )
  required_marker <- c("marker_id", "chromosome", "position_bp", "ref", "alt")
  if (!is.data.frame(x$samples) || !all(required_sample %in% names(x$samples))) {
    stop(
      "samples must contain sample_id, breeding_material_type, and breeding_material_type_source.",
      call. = FALSE
    )
  }
  if (!is.data.frame(x$markers) || !all(required_marker %in% names(x$markers))) {
    stop("markers is missing canonical marker columns.", call. = FALSE)
  }
  assert_unique_nonmissing(x$samples$sample_id, "samples$sample_id")
  if (anyNA(x$samples$breeding_material_type) ||
      any(!nzchar(trimws(as.character(x$samples$breeding_material_type))))) {
    stop("samples$breeding_material_type must be present; use 'unknown' rather than guessing.", call. = FALSE)
  }
  if (anyNA(x$samples$breeding_material_type_source) ||
      any(!nzchar(trimws(as.character(x$samples$breeding_material_type_source))))) {
    stop("samples$breeding_material_type_source must document the metadata origin.", call. = FALSE)
  }
  assert_unique_nonmissing(x$markers$marker_id, "markers$marker_id")

  if (!identical(as.character(x$samples$sample_id), rownames(x$genotype))) {
    stop("samples$sample_id must exactly match genotype row names and order.", call. = FALSE)
  }
  if (!identical(as.character(x$markers$marker_id), colnames(x$genotype))) {
    stop("markers$marker_id must exactly match genotype column names and order.", call. = FALSE)
  }
  marker_map_available <- !identical(x$provenance$marker_map_available, FALSE)
  if (marker_map_available) {
    if (anyNA(x$markers$chromosome) || any(!nzchar(as.character(x$markers$chromosome)))) {
      stop("Marker chromosome values must be present.", call. = FALSE)
    }
    if (anyNA(x$markers$position_bp) || any(x$markers$position_bp <= 0)) {
      stop("Marker positions must be positive and non-missing.", call. = FALSE)
    }
  } else {
    has_chromosome <- !is.na(x$markers$chromosome) & nzchar(as.character(x$markers$chromosome))
    if (any(has_chromosome) || any(!is.na(x$markers$position_bp))) {
      stop(
        "Marker-map provenance says unavailable, so chromosome and position_bp must remain missing.",
        call. = FALSE
      )
    }
  }

  if (!is.null(x$phenotypes)) {
    required_pheno <- c("sample_id", "trait", "value")
    if (!is.data.frame(x$phenotypes) || !all(required_pheno %in% names(x$phenotypes))) {
      stop("phenotypes must be a data.frame with sample_id, trait, and value.", call. = FALSE)
    }
    unknown <- setdiff(unique(as.character(x$phenotypes$sample_id)), x$samples$sample_id)
    if (length(unknown)) {
      stop("Phenotypes contain sample IDs absent from genotype: ",
           paste(utils::head(unknown, 10L), collapse = ", "), call. = FALSE)
    }
  }
  invisible(TRUE)
}

require_marker_map <- function(dataset, context = "This analysis") {
  validate_genomic_dataset(dataset)
  if (identical(dataset$provenance$marker_map_available, FALSE)) {
    stop(
      context,
      " requires verified chromosome and base-pair positions, but this dataset is explicitly unmapped.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

match_phenotype_to_genotype <- function(genotype_ids, phenotype,
                                        require_all_genotypes = FALSE) {
  if (!is.data.frame(phenotype) || !"sample_id" %in% names(phenotype)) {
    stop("phenotype must be a data.frame containing sample_id.", call. = FALSE)
  }
  assert_unique_nonmissing(genotype_ids, "genotype_ids")
  assert_unique_nonmissing(phenotype$sample_id, "phenotype$sample_id")

  genotype_ids <- as.character(genotype_ids)
  phenotype_ids <- as.character(phenotype$sample_id)
  matched_ids <- genotype_ids[genotype_ids %in% phenotype_ids]
  genotype_only <- setdiff(genotype_ids, phenotype_ids)
  phenotype_only <- setdiff(phenotype_ids, genotype_ids)

  if (!length(matched_ids)) {
    stop("No genotype and phenotype sample IDs match.", call. = FALSE)
  }
  if (isTRUE(require_all_genotypes) && length(genotype_only)) {
    stop("Some genotype samples lack phenotypes.", call. = FALSE)
  }

  list(
    phenotype = phenotype[match(matched_ids, phenotype_ids), , drop = FALSE],
    audit = data.frame(
      category = c("genotype", "phenotype", "matched", "genotype_only", "phenotype_only"),
      n = c(length(genotype_ids), length(phenotype_ids), length(matched_ids),
            length(genotype_only), length(phenotype_only)),
      stringsAsFactors = FALSE
    ),
    matched_ids = matched_ids,
    genotype_only = genotype_only,
    phenotype_only = phenotype_only
  )
}

marker_qc_metrics <- function(genotype) {
  if (!is.matrix(genotype) || !is.numeric(genotype)) {
    stop("genotype must be a numeric matrix.", call. = FALSE)
  }
  nonmissing_n <- colSums(!is.na(genotype))
  allele_count <- colSums(genotype, na.rm = TRUE)
  alt_frequency <- allele_count / (2 * nonmissing_n)
  maf <- pmin(alt_frequency, 1 - alt_frequency)
  data.frame(
    marker_id = colnames(genotype),
    missingness = colMeans(is.na(genotype)),
    alt_frequency = alt_frequency,
    maf = maf,
    monomorphic = nonmissing_n == 0L | maf == 0,
    stringsAsFactors = FALSE
  )
}

sample_qc_metrics <- function(genotype) {
  data.frame(
    sample_id = rownames(genotype),
    missingness = rowMeans(is.na(genotype)),
    heterozygosity = rowMeans(genotype == 1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
