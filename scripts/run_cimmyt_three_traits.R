args <- commandArgs(trailingOnly = TRUE)
phases <- c("core", "explain", "imputation", "report", "export_imputed")
if (!length(args) || !args[1L] %in% phases || length(args) > 4L ||
    (length(args) == 4L && !args[1L] %in% c("core", "explain"))) {
  stop("Usage: --args core|explain|imputation|report|export_imputed [input] [results/CIMMYT_three_traits] [core config or explanation learners]")
}
phase <- args[1L]
input <- if (length(args) >= 2L) args[2L] else "input"
output <- if (length(args) >= 3L) args[3L] else "results/CIMMYT_three_traits"
reference <- file.path(input, "reference.hmp.txt.gz")
phenotypes <- file.path(input, "phenotypes.csv")
metadata <- file.path(input, "sample_metadata.tsv")
learners <- c("gblup_reml", "elastic_net", "random_forest", "xgboost", "gaussian_rkhs")
units <- c(TKW = "g", grain_protein = "%", irrigated_yield = "t/ha")
source("breeding_pipeline.R")

if (phase == "core") {
  config <- if (length(args) == 4L) args[4L] else "config/breeding_cimmyt_tkw.yml"
  if (!requireNamespace("GAPIT", quietly = TRUE) ||
      utils::compareVersion(as.character(utils::packageVersion("GAPIT")), "4.1.0") < 0L) {
    stop("Load the pinned GAPIT library before starting the three-trait core run.")
  }
  panel <- utils::read.delim(metadata, colClasses = "character")
  if (!all(c("genetic_stratum", "fingerprint_group", "assessment_eligible") %in% names(panel))) {
    stop("Use the three-trait CIMMYT input with its saved cohort and genotype-group metadata.")
  }
  result <- run_breeding_traits(reference, phenotypes = phenotypes, traits = names(units), output = output,
    gp = learners, gwas = c("MLM", "FarmCPU", "BLINK"), sample_metadata = metadata,
    config = config, gwas_candidates = 10L, input_memory_gib = 8,
    assessment_fraction = .25, unit = units)
  cat("Completed core batch. Timing:", result$timing_file, "\n")
} else {
  completion <- readRDS(file.path(output, "multi_trait_completion.rds"))
  manifest <- utils::read.delim(completion$manifest_file, colClasses = "character")
  if (!identical(completion$status, "complete") || !identical(manifest$trait, names(units))) {
    stop("Supply a completed CIMMYT three-trait core run.")
  }
  if (phase == "explain" && length(args) == 4L) {
    chosen <- trimws(strsplit(args[4L], ",", fixed = TRUE)[[1L]])
    if (!length(chosen) || anyDuplicated(chosen) || !all(chosen %in% learners)) stop("Choose valid explanation learner IDs.")
    learners <- chosen
  }
  timing <- file.path(output, "optional_timing.tsv")
  timed <- function(step, trait, code) {
    started <- proc.time()[["elapsed"]]
    status <- "incomplete"
    on.exit(data.table::fwrite(data.frame(phase = phase, step = step, trait = trait,
      elapsed_seconds = proc.time()[["elapsed"]] - started, status = status), timing,
      sep = "\t", append = file.exists(timing), na = "NA"))
    value <- force(code)
    status <- "complete"
    value
  }
  timed("phase_total", NA_character_, {
    imported <- if (phase %in% c("explain", "imputation")) timed("genotype_import", NA_character_,
      read_breeding_genotypes(reference, maximum_memory_gib = 8))
    for (i in seq_len(nrow(manifest))) {
      trait <- manifest$trait[i]
      folder <- manifest$output[i]
      unit <- manifest$unit[i]
      cat(phase, trait, "\n")
      timed("trait_total", trait, {
        if (phase %in% c("explain", "imputation")) {
          model <- file.path(folder, "assessment", "GP")
          ids <- completion$traits[[trait]]$assessment$assessment_ids
          query <- imported$genotype[ids, , drop = FALSE]
          if (phase == "explain") {
            explain_breeding_models(model, query, imported$markers, learners, file.path(folder, "explanations"),
              query_rows = 24L, background_rows = 8L, pairs = 16L, region_width_bp = 10000000, trait = trait, unit = unit)
          } else {
            validate_breeding_imputation(model, query, imported$markers, file.path(folder, "imputation_validation"),
              maximum_masks = 50000L)
          }
        } else if (phase == "report") {
          explanations <- file.path(folder, "explanations")
          imputation <- file.path(folder, "imputation_validation", "accuracy_summary.tsv")
          conditional <- file.path(folder, "conditional_followup", "conditional_results.tsv")
          plot_breeding_report(file.path(folder, "report"), trait_output = folder, observed = phenotypes, unit = unit,
            explanation_output = if (file.exists(file.path(explanations, "explanation_completion.rds"))) explanations,
            imputation_accuracy = if (file.exists(imputation)) imputation,
            conditional_results = if (file.exists(conditional)) conditional)
        } else {
          export_breeding_imputed(file.path(folder, "GP", "preprocessor.rds"), file.path(folder, "imputed", "final"))
        }
      })
    }
    invisible(NULL)
  })
  cat("Completed", phase, "for all three traits. Timing:", timing, "\n")
}
