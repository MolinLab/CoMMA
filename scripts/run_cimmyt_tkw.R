args <- commandArgs(trailingOnly = TRUE)
phases <- c("check", "core", "explain", "imputation", "report", "conditional_prepare", "conditional",
            "conditional_plot", "conditional_support", "export_imputed")
if (!length(args) || !args[1L] %in% phases ||
    length(args) > switch(args[1L], conditional = 6L, conditional_prepare = 5L, conditional_plot = 4L,
      conditional_support = 4L, report = 5L, export_imputed = 4L, 3L)) {
  stop("Usage: R --vanilla --slave -f scripts/run_cimmyt_tkw.R --args check|core|explain|imputation|report [input] [results/TKW_full] [new_report_directory (report only)]\n",
       "Prepare candidate inputs: --args conditional_prepare input results/TKW_full shuffle_column[,shuffle_column] [explanation_fraction]\n",
       "Run candidate follow-up: --args conditional [input] [results/TKW_full] [max_candidates] [association|shap] [new_followup_directory]\n",
       "Plot saved follow-up: --args conditional_plot [input] [results/TKW_full] [followup_directory]\n",
       "Diagnose genotype support: --args conditional_support [input] [results/TKW_full] [followup_directory]\n",
       "Combine figures: --args report [input] [results/TKW_full] [new_report_directory] [followup_directory]\n",
       "Export imputed dosages: --args export_imputed [input] [results/TKW_full] [final|assessment]")
}
phase <- args[1L]
library_pointer <- file.path("..", ".renv", "breeding_library.txt")
if (identical(Sys.info()[["sysname"]], "Linux") && file.exists(library_pointer)) {
  pinned_library <- trimws(readLines(library_pointer, n = 1L, warn = FALSE))
  if (length(pinned_library) == 1L && nzchar(pinned_library) && dir.exists(pinned_library)) {
    .libPaths(c(pinned_library, .libPaths()))
  }
}
if (phase %in% c("check", "core", "conditional")) {
  gapit <- tryCatch(loadNamespace("GAPIT"), error = function(err) {
    stop("GAPIT failed to load: ", conditionMessage(err), "\nR: ", R.version.string,
      "\nLibrary search paths: ", paste(.libPaths(), collapse = ", "),
      "\nLinux library pointer: ", library_pointer, call. = FALSE)
  })
  cat(R.version.string, "\nGAPIT ", as.character(getNamespaceVersion(gapit)), " from ",
    getNamespaceInfo(gapit, "path"), "\n", sep = "")
  if (!"GAPIT" %in% getNamespaceExports(gapit) || utils::compareVersion(as.character(getNamespaceVersion(gapit)), "4.1.0") < 0L) {
    stop("This workflow requires a GAPIT installation compatible with version 4.1.0.", call. = FALSE)
  }
  if (phase == "check") quit(save = "no")
}
input <- if (length(args) >= 2L) args[2L] else "input"
output <- if (length(args) >= 3L) args[3L] else "results/TKW_full"
reference <- file.path(input, "reference.hmp.txt.gz")
phenotypes <- file.path(input, "phenotypes.csv")
metadata <- file.path(input, "sample_metadata.tsv")
learners <- c("gblup_reml", "elastic_net", "random_forest", "xgboost", "gaussian_rkhs")
source("breeding_pipeline.R")

if (phase == "core") {
  panel <- utils::read.delim(metadata, colClasses = "character")
  if (!all(c("genetic_stratum", "fingerprint_group", "assessment_eligible") %in% names(panel))) {
    stop("Use the full TKW v2 input, which records the genetic strata and lines reserved for assessment.")
  }
  result <- run_breeding_trait(reference, phenotypes = phenotypes, trait = "TKW", output = output,
    gp = learners, gwas = c("MLM", "FarmCPU", "BLINK"), sample_metadata = metadata,
    config = "config/breeding_cimmyt_tkw.yml", gwas_candidates = 10L, input_memory_gib = 8,
    assessment_fraction = .25, unit = "g")
  cat("Completed TKW core. Registered holdout lines:", length(result$assessment$assessment_ids), "\n")
} else if (phase %in% c("explain", "imputation")) {
  receipt <- readRDS(file.path(output, "assessment", "assessment_completion.rds"))
  query <- read_breeding_genotypes(reference, sample_ids = receipt$assessment_ids, maximum_memory_gib = 8)
  model <- file.path(output, "assessment", "GP")
  if (phase == "explain") {
    explain_breeding_models(model, query$genotype, query$markers, learners,
      file.path(output, "explanations"), query_rows = 24L, background_rows = 8L, pairs = 16L,
      region_width_bp = 10000000, trait = "TKW", unit = "g")
  } else {
    validate_breeding_imputation(model, query$genotype, query$markers,
      output = file.path(output, "imputation_validation"), maximum_masks = 50000L)
  }
} else if (phase == "conditional_prepare") {
  if (length(args) < 4L || !nzchar(args[4L])) {
    stop("Choose a shuffle column explicitly, for example: --args conditional_prepare input results/TKW_full source_group")
  }
  fraction <- if (length(args) == 5L) suppressWarnings(as.numeric(args[5L])) else .25
  prepared <- prepare_breeding_conditional_inputs(assessment = file.path(output, "assessment"),
    phenotypes = phenotypes, sample_metadata = metadata, trait = "TKW", output = file.path(input, "conditional_followup"),
    shuffle_columns = trimws(strsplit(args[4L], ",", fixed = TRUE)[[1L]]), explanation_fraction = fraction)
  cat("Wrote follow-up inputs for", prepared$n_discovery, "discovery,", prepared$n_confirmation,
      "confirmation-fit and", prepared$n_explanation, "explanation lines.\n")
} else if (phase == "export_imputed") {
  panel <- if (length(args) == 4L) args[4L] else "final"
  directory <- switch(panel, final = "GP", assessment = file.path("assessment", "GP"),
                      stop("Choose final or assessment for export_imputed."))
  exported <- export_breeding_imputed(file.path(output, directory, "preprocessor.rds"), file.path(output, "imputed", panel))
  cat("Exported", exported$rows, "lines by", exported$markers, "markers to", exported$matrix_file, "\n")
} else if (phase == "conditional") {
  inputs <- file.path(input, "conditional_followup")
  design <- utils::read.delim(file.path(inputs, "design.tsv"), colClasses = "character")
  assessment <- readRDS(file.path(output, "assessment", "assessment_completion.rds"))
  if (!"sample_id" %in% names(design) || !all(design$sample_id %in% assessment$assessment_ids)) {
    stop("Use registered assessment IDs for this exploratory TKW follow-up.")
  }
  candidates <- if (length(args) >= 4L) suppressWarnings(as.numeric(args[4L])) else 10L
  mode <- if (length(args) >= 5L) args[5L] else "association"
  if (!mode %in% c("association", "shap")) stop("Choose association or shap for the conditional phase.")
  target <- if (length(args) == 6L) args[6L] else file.path(output, "conditional_followup")
  result <- run_breeding_conditional_followup(model = file.path(output, "assessment", "GP"),
    query = reference, design = design, discovery_phenotypes = file.path(inputs, "discovery_phenotypes.tsv"),
    confirmation_phenotypes = file.path(inputs, "confirmation_phenotypes.tsv"), trait = "TKW", output = target,
    sample_metadata = metadata, max_candidates = candidates, shap = mode == "shap", input_memory_gib = 8)
  cat("Completed exploratory TKW candidate follow-up:", result$results_file, "\n")
  plot_breeding_report(file.path(target, "plots"), conditional_results = result$results_file, trait = "TKW")
} else if (phase == "conditional_support") {
  target <- if (length(args) == 4L) args[4L] else file.path(output, "conditional_followup")
  diagnostic <- diagnose_breeding_conditional_support(followup = target, model = file.path(output, "assessment", "GP"),
    query = reference, output = file.path(target, "support"), input_memory_gib = 8)
  cat("Conditional genotype support:", diagnostic$candidate_support, "\n")
} else if (phase == "conditional_plot") {
  target <- if (length(args) == 4L) args[4L] else file.path(output, "conditional_followup")
  figure <- plot_breeding_report(file.path(target, "plots"),
    conditional_results = file.path(target, "conditional_results.tsv"), trait = "TKW")
  cat("Conditional follow-up figure:", figure$figures, "\n")
} else {
  explanations <- file.path(output, "explanations")
  imputation <- file.path(output, "imputation_validation", "accuracy_summary.tsv")
  report_output <- if (length(args) >= 4L) args[4L] else file.path(output, "report")
  followup <- if (length(args) == 5L) args[5L] else file.path(output, "conditional_followup")
  conditional <- file.path(followup, "conditional_results.tsv")
  plot_breeding_report(report_output, trait_output = output, observed = phenotypes, unit = "g",
    explanation_output = if (file.exists(file.path(explanations, "explanation_completion.rds"))) explanations,
    imputation_accuracy = if (file.exists(imputation)) imputation,
    conditional_results = if (length(args) == 5L || file.exists(conditional)) conditional)
}
