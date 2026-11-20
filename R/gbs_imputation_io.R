gbs_project_path <- function(path, must_exist = TRUE) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) ||
      grepl("^[A-Za-z]:", path) || startsWith(path, "/") ||
      startsWith(path, intToUtf8(92L))) stop("Use a project-relative path.")
  path <- gsub(intToUtf8(92L), "/", path, fixed = TRUE)
  if (any(strsplit(path, "/", fixed = TRUE)[[1L]] == "..")) stop("Parent traversal is forbidden.")
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  absolute <- file.path(root, path)
  parent <- if (file.exists(absolute) || dir.exists(absolute)) absolute else dirname(absolute)
  while (!file.exists(parent) && !dir.exists(parent)) parent <- dirname(parent)
  resolved <- normalizePath(parent, winslash = "/", mustWork = TRUE)
  inside <- identical(tolower(resolved), tolower(root)) ||
    startsWith(tolower(resolved), paste0(tolower(root), "/"))
  if (!inside) stop("The resolved path leaves the project.")
  if (must_exist && !file.exists(absolute)) stop("Input is absent: ", path)
  absolute
}

gbs_decode_hapmap_calls <- function(calls, alleles) {
  pair <- strsplit(toupper(alleles), "/", fixed = TRUE)[[1L]]
  if (length(pair) != 2L || any(!pair %in% c("A", "C", "G", "T")) ||
      pair[[1L]] == pair[[2L]]) stop("A biallelic allele definition is required.")
  meanings <- c(A = "AA", C = "CC", G = "GG", T = "TT", R = "AG", Y = "CT",
    S = "CG", W = "AT", K = "GT", M = "AC")
  calls <- toupper(trimws(calls))
  missing <- is.na(calls) | calls %in% c("", "N", "NN", "NA", ".", "-", "--", "?", "00")
  translated <- unname(meanings[calls])
  double <- !missing & nchar(calls) == 2L
  translated[double] <- calls[double]
  a <- substr(translated, 1L, 1L)
  b <- substr(translated, 2L, 2L)
  invalid <- !missing & (is.na(translated) | !a %in% pair | !b %in% pair)
  dosage <- as.numeric(a == pair[[2L]]) + as.numeric(b == pair[[2L]])
  dosage[missing | invalid] <- NA_real_
  list(dosage = dosage, invalid = invalid)
}

read_gbs_hapmap <- function(path, sample_ids = NULL, reference_markers = NULL,
                            maximum_memory_gib = 2) {
  absolute <- gbs_project_path(path)
  con <- if (endsWith(tolower(path), ".gz")) gzfile(absolute, "rt") else file(absolute, "rt")
  on.exit(close(con))
  header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1L]]
  if (length(header) < 12L || header[[1L]] != "rs#" || match("QCcode", header) != 11L) {
    stop("Expected a standard HapMap header with eleven metadata columns.")
  }
  ids <- header[-seq_len(11L)]
  if (anyDuplicated(ids) || any(!nzchar(ids))) stop("Invalid HapMap sample IDs.")
  if (is.null(sample_ids)) sample_ids <- ids
  if (!length(sample_ids) || anyNA(sample_ids) || anyDuplicated(sample_ids) ||
      anyNA(match(sample_ids, ids))) stop("Sample selection contains absent/duplicate IDs.")
  selected <- match(sample_ids, ids) + 11L
  if (!is.null(reference_markers) && (anyDuplicated(reference_markers$marker_id) ||
      !"alleles" %in% names(reference_markers))) stop("Fitted HapMap allele metadata is required.")
  columns <- metadata <- exclusions <- list()
  seen <- new.env(parent = emptyenv(), hash = TRUE)
  record_exclusion <- function(id, reason) {
    exclusions[[length(exclusions) + 1L]] <<- data.frame(marker_id = id, reason = reason)
  }
  repeat {
    lines <- readLines(con, n = 32L)
    if (!length(lines)) break
    for (line in lines) {
      prefix <- strsplit(substr(line, 1L, 512L), "\t", fixed = TRUE)[[1L]]
      if (length(prefix) < 12L) stop("Truncated HapMap row.")
      id <- prefix[[1L]]
      if (!is.null(reference_markers) && !id %in% reference_markers$marker_id) next
      if (exists(id, envir = seen, inherits = FALSE)) stop("Duplicate HapMap marker: ", id)
      assign(id, TRUE, envir = seen)
      chr <- prefix[[3L]]
      position <- suppressWarnings(as.numeric(prefix[[4L]]))
      if (is.na(position) || !is.finite(position) || position <= 0 ||
          chr %in% c("", "UN", "NA", "0")) {
        if (!is.null(reference_markers)) stop("A fitted marker is unplaced in the query: ", id)
        record_exclusion(id, "unplaced_or_invalid_position")
        next
      }
      alleles <- toupper(prefix[[2L]])
      if (!is.null(reference_markers)) {
        ref <- reference_markers[match(id, reference_markers$marker_id), , drop = FALSE]
        if (ref$chromosome != chr || ref$position_bp != position) stop("Query marker coordinates disagree: ", id)
        declared <- strsplit(alleles, "/", fixed = TRUE)[[1L]]
        canonical <- strsplit(ref$alleles, "/", fixed = TRUE)[[1L]]
        if (alleles != "NA" && any(!declared %in% canonical)) stop("Query allele declaration conflicts: ", id)
        if ("strand" %in% names(ref) && ref$strand %in% c("+", "-") &&
            prefix[[5L]] %in% c("+", "-") && ref$strand != prefix[[5L]]) stop("Query strand differs: ", id)
        alleles <- ref$alleles
      }
      if (!grepl("^[ACGT]/[ACGT]$", alleles) || substr(alleles, 1L, 1L) == substr(alleles, 3L, 3L)) {
        record_exclusion(id, "non_biallelic_allele_definition")
        next
      }
      proposed_cells <- as.double(length(sample_ids)) * (length(columns) + 1L)
      if (proposed_cells * 8 * 5 / 1024^3 + 0.10 > maximum_memory_gib) {
        stop("Streaming import reaches its memory budget. Select fewer line IDs or use a server budget; no markers were silently truncated.")
      }
      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(fields) != length(header)) stop("HapMap row width differs from header: ", id)
      decoded <- gbs_decode_hapmap_calls(fields[selected], alleles)
      if (any(decoded$invalid)) {
        if (!is.null(reference_markers)) stop("Query calls conflict with fitted marker alleles: ", id)
        record_exclusion(id, "call_outside_declared_allele_pair")
        next
      }
      columns[[id]] <- decoded$dosage
      metadata[[id]] <- data.frame(marker_id = id, alleles = alleles,
        chromosome = chr, position_bp = position, strand = prefix[[5L]],
        assembly = prefix[[6L]], dosage_allele = substr(alleles, 3L, 3L),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(columns)) stop("No valid genotype markers remain.")
  genotype <- do.call(cbind, columns)
  rownames(genotype) <- sample_ids
  markers <- do.call(rbind, metadata)
  rownames(markers) <- NULL
  if (!is.null(reference_markers)) {
    if (any(!reference_markers$marker_id %in% colnames(genotype))) stop("Query is missing fitted markers.")
    genotype <- genotype[, reference_markers$marker_id, drop = FALSE]
    markers <- markers[match(reference_markers$marker_id, markers$marker_id), , drop = FALSE]
  }
  gbs_validate_calls(genotype)
  list(genotype = genotype, markers = markers,
    exclusions = if (length(exclusions)) do.call(rbind, exclusions) else
      data.frame(marker_id = character(), reason = character()),
    provenance = list(source_path = path, source_sha256 = digest::digest(absolute, algo = "sha256", file = TRUE),
      source_genotypes_modified = FALSE, raw_call_import = TRUE))
}

read_gbs_tool_input <- function(path, sample_ids_file = NULL, reference_markers = NULL,
                                maximum_memory_gib = 2) {
  absolute <- gbs_project_path(path)
  samples <- if (is.null(sample_ids_file)) NULL else {
    ids <- trimws(readLines(gbs_project_path(sample_ids_file)))
    if (any(!nzchar(ids)) || anyDuplicated(ids)) stop("Sample-ID file has blanks or duplicates.")
    ids
  }
  if (grepl("[.]hmp[.]txt([.]gz)?$", tolower(path))) {
    return(read_gbs_hapmap(path, samples, reference_markers, maximum_memory_gib))
  }
  if (!endsWith(tolower(path), ".rds")) stop("Supported inputs are HapMap (.hmp.txt[.gz]) or canonical genotype RDS.")
  object <- readRDS(absolute)
  if (!is.list(object) || is.null(object$genotype) || is.null(object$markers)) {
    stop("RDS input must contain genotype and marker metadata.")
  }
  if (!"alleles" %in% names(object$markers) && all(c("ref", "alt") %in% names(object$markers))) {
    object$markers$alleles <- paste(object$markers$ref, object$markers$alt, sep = "/")
  }
  if (!all(c("marker_id", "chromosome", "position_bp", "alleles") %in% names(object$markers)) ||
      anyDuplicated(object$markers$marker_id) ||
      anyNA(match(colnames(object$genotype), object$markers$marker_id))) {
    stop("Numeric RDS needs explicit unique marker IDs, positions and allele coding.")
  }
  gbs_validate_calls(object$genotype)
  if (!is.null(samples)) {
    if (anyNA(match(samples, rownames(object$genotype)))) stop("Requested sample IDs are absent.")
    object$genotype <- object$genotype[samples, , drop = FALSE]
  }
  if (!is.null(reference_markers)) {
    if (anyNA(match(reference_markers$marker_id, colnames(object$genotype)))) stop("Query lacks fitted markers.")
    qmap <- object$markers[match(reference_markers$marker_id, object$markers$marker_id), , drop = FALSE]
    if (anyNA(qmap$marker_id) || any(qmap$chromosome != reference_markers$chromosome) ||
        any(qmap$position_bp != reference_markers$position_bp)) stop("RDS query marker metadata disagree.")
    if ("alleles" %in% names(reference_markers) &&
        (!"alleles" %in% names(qmap) || any(qmap$alleles != reference_markers$alleles))) {
      stop("Numeric RDS query allele coding differs; harmonize from original calls first.")
    }
    object$genotype <- object$genotype[, reference_markers$marker_id, drop = FALSE]
    object$markers <- qmap
  }
  if (as.double(length(object$genotype)) * 8 * 5 / 1024^3 + 0.10 > maximum_memory_gib) {
    stop("Loaded RDS exceeds the configured processing budget; use a smaller prepared input.")
  }
  object
}
