args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 3:4) stop("Start the workflow with Run_pipeline.R in RStudio or run_pipeline.ps1 in PowerShell.")
source("breeding_pipeline.R")
cfg <- yaml::read_yaml(args[1L])
if (length(args) == 4L) cfg$output <- args[4L]
if (file.exists(cfg$output)) stop("Choose a new output folder in the run configuration.")
workers <- as.integer(args[2L])
memory <- as.numeric(args[3L])
if (is.na(workers) || workers < 1L || !is.finite(memory) || memory < 1) stop("Invalid worker or memory setting.")
if (any(!vapply(list(cfg$validate_imputation, cfg$export_imputed), function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}, logical(1L)))) stop("Set validate_imputation and export_imputed to true or false.")
explain <- unlist(cfg$explain_learners, use.names = FALSE)
grouping <- if (is.null(cfg$explain_grouping)) "coordinate" else cfg$explain_grouping
if (!identical(grouping, "coordinate") && !identical(grouping, "ld")) stop("Set explain_grouping to coordinate or ld.")
ld <- if (is.null(cfg$explain_ld)) list() else cfg$explain_ld
gp <- as.character(unlist(cfg$gp, use.names = FALSE))
gwas <- as.character(unlist(cfg$gwas, use.names = FALSE))
if (anyDuplicated(explain) || !all(explain %in% gp)) stop("Choose explanation learners from the GP learners.")
if ((length(explain) || cfg$validate_imputation) && is.null(cfg$assessment_fraction)) {
  stop("Set assessment_fraction to run regional SHAP or masked-call validation.")
}
if (!length(gp) && !length(gwas)) stop("Select at least one GP or GWAS method.")
model_config <- yaml::read_yaml(cfg$model_config)
model_config$resources <- list(global_cpu_budget = workers, reserved_cpus = 0L, workers = workers,
  worker_threads = 1L, reported_ram_gb = ceiling(memory * 1024^3 / 1e9), memory_budget_gib = memory)
cat("Worker cap:", workers, "| memory budget:", memory, "GiB\n")
started <- proc.time()[["elapsed"]]
input <- read_breeding_genotypes(cfg$reference, markers = cfg$markers, maximum_memory_gib = memory)
import_elapsed <- proc.time()[["elapsed"]] - started

batch <- run_breeding_traits(input$genotype, input$markers, cfg$phenotypes, cfg$output,
  traits = unlist(cfg$traits, use.names = FALSE), gp = gp, gwas = gwas,
  sample_metadata = cfg$sample_metadata, config = model_config, gwas_candidates = cfg$gwas_candidates,
  assessment_fraction = cfg$assessment_fraction,
  unit = if (is.null(cfg$unit)) "trait units" else unlist(cfg$unit))
manifest <- utils::read.delim(batch$manifest_file, colClasses = "character")
timing <- file.path(cfg$output, "first_run_timing.tsv")
record <- function(trait, step, elapsed, status) {
  data.table::fwrite(data.frame(trait, step, elapsed_seconds = elapsed, status), timing,
                     sep = "\t", append = file.exists(timing), na = "NA")
}
record(NA_character_, "genotype_import", import_elapsed, "complete")
timed <- function(trait, step, expr) {
  start <- proc.time()[["elapsed"]]
  status <- "incomplete"
  on.exit(record(trait, step, proc.time()[["elapsed"]] - start, status))
  value <- force(expr)
  status <- "complete"
  value
}
for (i in seq_len(nrow(manifest))) {
  trait <- manifest$trait[i]
  folder <- manifest$output[i]
  unit <- manifest$unit[i]
  model <- file.path(folder, "assessment", "GP")
  explanation <- imputation <- NULL
  if (length(explain) || cfg$validate_imputation) {
    ids <- batch$traits[[trait]]$assessment$assessment_ids
    holdout <- input$genotype[ids, , drop = FALSE]
    if (length(explain)) {
      explanation <- file.path(folder, "explanations")
      timed(trait, "explain", explain_breeding_models(model, holdout, input$markers, explain,
        explanation, trait = trait, unit = unit, plots = FALSE, grouping = grouping, ld_options = ld))
    }
    if (cfg$validate_imputation) {
      validation <- timed(trait, "imputation_validation", validate_breeding_imputation(model, holdout,
        input$markers, output = file.path(folder, "imputation_validation"), plots = FALSE))
      imputation <- validation$accuracy_file
    }
  }
  if (cfg$export_imputed) {
    prep <- if (length(gp)) "GP" else "preprocessing"
    timed(trait, "export_imputed", export_breeding_imputed(file.path(folder, prep, "preprocessor.rds"),
      file.path(folder, "imputed_final")))
  }
  timed(trait, "report", plot_breeding_report(file.path(folder, "report"), trait_output = folder,
    observed = cfg$phenotypes, explanation_output = explanation, imputation_accuracy = imputation, unit = unit))
}
record(NA_character_, "workflow_total", proc.time()[["elapsed"]] - started, "complete")
writeLines(capture.output(utils::sessionInfo()), file.path(cfg$output, "session_info.txt"))
cat("Finished. Reports are in the trait folders listed in", batch$manifest_file, "\n")
