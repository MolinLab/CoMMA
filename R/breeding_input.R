breeding_delimited_table <- function(path, dosage = FALSE, allow_empty = FALSE) {
  filename <- gbs_project_path(path)
  suffix <- sub("[.]gz$", "", tolower(path))
  if (!grepl("[.](csv|tsv|txt)$", suffix)) stop("Use a CSV or tab-separated TSV/TXT file.", call. = FALSE)
  sep <- if (endsWith(suffix, ".csv")) "," else "\t"
  read_table <- function(nrows, classes) {
    con <- if (endsWith(tolower(path), ".gz")) gzfile(filename, "rt") else file(filename, "rt")
    on.exit(close(con))
    utils::read.table(con, header = TRUE, sep = sep, nrows = nrows, colClasses = classes,
      row.names = NULL, check.names = FALSE, comment.char = "", quote = "\"", fill = FALSE,
      na.strings = c("", "NA"), stringsAsFactors = FALSE)
  }
  header <- names(read_table(0L, "character"))
  if (anyDuplicated(header) || any(!nzchar(header))) stop("Input column names must be unique and non-empty.")
  if (dosage && (length(header) < 2L || !tolower(header[1L]) %in% c("sample_id", "taxa"))) {
    stop("Start the dosage file with sample_id or Taxa, followed by one column per SNP.", call. = FALSE)
  }
  classes <- if (dosage) c("character", rep("numeric", length(header) - 1L)) else "character"
  x <- read_table(-1L, classes)
  if (!identical(names(x), header) || (!allow_empty && !nrow(x))) stop("Input rows must match the header and contain data.")
  x
}

breeding_input_name <- function(table, name, aliases = character()) {
  i <- which(tolower(names(table)) %in% tolower(c(name, aliases)))
  if (length(i) > 1L) stop("Input has more than one column for ", name, ".", call. = FALSE)
  if (length(i)) names(table)[i] <- name
  table
}

breeding_input_numeric <- function(x, label) {
  y <- suppressWarnings(as.numeric(x))
  if (any(!is.na(x) & !is.finite(y))) stop(label, " must contain finite numeric values or NA.", call. = FALSE)
  y
}

breeding_input_markers <- function(markers) {
  if (!is.character(markers)) return(markers)
  x <- breeding_delimited_table(markers)
  x <- breeding_input_name(x, "marker_id", c("SNP", "Name"))
  x <- breeding_input_name(x, "chromosome", "chr")
  x <- breeding_input_name(x, "position_bp", c("Position", "pos"))
  x <- breeding_input_name(x, "alleles")
  x <- breeding_input_name(x, "ref")
  x <- breeding_input_name(x, "alt")
  if (!all(c("marker_id", "chromosome", "position_bp") %in% names(x))) {
    stop("The marker map needs marker ID, chromosome and base-pair position columns.", call. = FALSE)
  }
  x$position_bp <- breeding_input_numeric(x$position_bp, "Marker positions")
  if (anyNA(x$position_bp) || any(x$position_bp <= 0 | x$position_bp != floor(x$position_bp)) ||
      anyNA(x$chromosome) || any(!nzchar(x$chromosome))) stop("Supply mapped markers with positive whole base-pair positions.")
  for (field in intersect(c("alleles", "ref", "alt"), names(x))) x[[field]] <- toupper(x[[field]])
  assert_unique_nonmissing(x$marker_id, "Marker IDs")
  x
}

breeding_input_phenotypes <- function(phenotypes, trait) {
  if (!is.character(phenotypes)) return(phenotypes)
  x <- breeding_input_name(breeding_delimited_table(phenotypes), "sample_id", "Taxa")
  if (!"sample_id" %in% names(x)) stop("The phenotype file needs a sample_id or Taxa column.", call. = FALSE)
  if (!all(c("trait", "value") %in% names(x))) {
    if (!trait %in% names(x)) stop("The phenotype file has no column named ", trait, ".", call. = FALSE)
    x <- data.frame(sample_id = x$sample_id, trait = trait, value = x[[trait]], stringsAsFactors = FALSE)
  }
  x$value <- breeding_input_numeric(x$value, "Phenotypes")
  x
}

breeding_input_trait_table <- function(phenotypes, traits = NULL) {
  x <- if (is.character(phenotypes)) breeding_delimited_table(phenotypes) else phenotypes
  if (!is.data.frame(x) || !nrow(x) || anyDuplicated(names(x))) stop("Supply one phenotype table with unique column names.")
  x <- breeding_input_name(x, "sample_id", "Taxa")
  if (!"sample_id" %in% names(x)) stop("The phenotype table needs a sample_id or Taxa column.")
  x$sample_id <- as.character(x$sample_id)
  long <- all(c("trait", "value") %in% names(x))
  available <- if (long) unique(as.character(x$trait)) else setdiff(names(x), "sample_id")
  assert_unique_nonmissing(available, "Available phenotype traits")
  if (is.null(traits)) traits <- available
  if (!is.character(traits) || !all(traits %in% available)) stop("Choose traits present in the phenotype table.")
  assert_unique_nonmissing(traits, "Selected phenotype traits")
  if (long) {
    x <- x[x$trait %in% traits, c("sample_id", "trait", "value"), drop = FALSE]
    x$trait <- as.character(x$trait)
  } else {
    x <- do.call(rbind, lapply(traits, function(trait) data.frame(sample_id = x$sample_id, trait = trait,
      value = as.character(x[[trait]]))))
  }
  if (is.factor(x$value)) x$value <- as.character(x$value)
  x$value <- breeding_input_numeric(x$value, "Phenotypes")
  list(phenotypes = x, traits = traits)
}

breeding_input_metadata <- function(metadata) {
  if (!is.character(metadata)) return(metadata)
  x <- breeding_input_name(breeding_delimited_table(metadata), "sample_id", "Taxa")
  breeding_input_name(x, "source_group")
}

breeding_read_genotypes <- function(path, markers = NULL, sample_ids = NULL, maximum_memory_gib = 2,
                                    reference_markers = NULL) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) stop("Supply one genotype file path.")
  if (!grepl("[.](csv|tsv|txt)([.]gz)?$", tolower(path))) {
    stop("Choose a HapMap (.hmp.txt or .hmp.txt.gz), CSV, or tab-separated TSV/TXT genotype file.", call. = FALSE)
  }
  if (!is.numeric(maximum_memory_gib) || length(maximum_memory_gib) != 1L ||
      !is.finite(maximum_memory_gib) || maximum_memory_gib <= 0) stop("Provide a positive HapMap import memory budget in GiB.")
  markers <- breeding_input_markers(markers)
  if (grepl("[.]hmp[.]txt([.]gz)?$", tolower(path))) {
    required <- reference_markers
    if (!is.null(markers)) {
      ids <- if (is.null(required)) markers$marker_id else required$marker_id
      supplied <- gbs_preprocessing_map(markers, ids, required)
      if (is.null(required)) required <- supplied
    }
    return(read_gbs_hapmap(path, sample_ids, required, maximum_memory_gib))
  }
  if (is.null(markers)) stop("Dosage CSV or TSV requires a separate marker map with explicit alleles.", call. = FALSE)
  x <- breeding_delimited_table(path, dosage = TRUE)
  ids <- x[[1L]]
  assert_unique_nonmissing(ids, "Genotype sample IDs")
  if (!is.null(sample_ids)) {
    assert_unique_nonmissing(sample_ids, "Requested genotype sample IDs")
    i <- match(sample_ids, ids)
    if (!length(i) || anyNA(i)) stop("Requested sample IDs are absent from the dosage file.")
    x <- x[i, , drop = FALSE]
    ids <- ids[i]
  }
  columns <- if (is.null(reference_markers)) names(x)[-1L] else reference_markers$marker_id
  if (any(!columns %in% names(x)[-1L])) stop("The query dosage file is missing fitted predictor SNPs.")
  map <- gbs_preprocessing_map(markers, columns, reference_markers)
  calls <- as.matrix(x[columns])
  rownames(calls) <- ids
  gbs_validate_calls(calls)
  list(genotype = calls, markers = map,
       exclusions = data.frame(marker_id = character(), reason = character()),
       provenance = list(source_path = path, raw_call_import = TRUE))
}

breeding_read_query <- function(query, markers, model, maximum_memory_gib) {
  if (!is.character(query) || is.matrix(query)) {
    return(list(genotype = query, markers = breeding_input_markers(markers)))
  }
  if (!inherits(model, "breeding_trait_model")) stop("A saved breeding trait model is required for query file import.")
  breeding_read_genotypes(query, markers, maximum_memory_gib = maximum_memory_gib,
                         reference_markers = model$preprocessing$imputer$markers)
}
