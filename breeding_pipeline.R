# Source from the project/release root after setting up its R library.
local({
  destination <- parent.env(environment()); root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  backend <- new.env(parent = parent.env(globalenv()))
  sys.source(file.path(root, "R/gbs_parallel.R"), envir = backend)
  sys.source(file.path(root, "R/breeding_pipeline.R"), envir = backend)
  for (path in setdiff(backend$breeding_pipeline_sources(), c("R/gbs_parallel.R", "R/breeding_pipeline.R"))) sys.source(file.path(root, path), envir = backend)
  for (path in c("R/breeding_interpretation.R", "R/breeding_group_shapley.R",
                 "R/breeding_coordinate_shapley.R", "R/breeding_ld_shapley.R", "R/breeding_gwas_shortlist.R", "R/breeding_input.R",
                 "R/breeding_chromosome_plot.R", "R/breeding_report.R", "R/breeding_release.R",
                 "R/breeding_xgb_precision.R", "R/breeding_conditional_followup.R")) {
    sys.source(file.path(root, path), envir = backend)
  }
  run <- function(fun, arguments) {
    old_wd <- getwd(); variables <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS")
    old_env <- Sys.getenv(variables, unset = NA_character_); old_threads <- data.table::getDTthreads()
    on.exit({ setwd(old_wd); data.table::setDTthreads(old_threads)
      for (i in seq_along(variables)) if (is.na(old_env[i])) Sys.unsetenv(variables[i]) else do.call(Sys.setenv, setNames(list(old_env[i]), variables[i]))
    }, add = TRUE)
    setwd(root); data.table::setDTthreads(1L)
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1")
    do.call(fun, arguments)
  }
  assign("fit_breeding_models", function(...) run(backend$breeding_fit_models, list(...)), envir = destination)
  assign("prepare_breeding_panel", function(...) run(backend$breeding_prepare_panel, list(...)), envir = destination)
  assign("read_breeding_genotypes", function(path, markers = NULL, sample_ids = NULL, maximum_memory_gib = 2) {
    run(backend$breeding_read_genotypes, list(path, markers, sample_ids, maximum_memory_gib))
  }, envir = destination)
  assign("predict_breeding_models", function(model, query, markers = NULL, sample_metadata = NULL, output,
                                             resources = NULL, learners = NULL, input_memory_gib = 2) {
    run(function() {
      file_query <- is.character(query) && !is.matrix(query)
      if (is.character(model) && (dir.exists(backend$gbs_project_path(model)) || !is.null(learners) || file_query)) {
        model <- backend$breeding_load_model(model)
      }
      if (!is.null(learners)) {
        if (!inherits(model, "breeding_trait_model")) stop("A saved breeding trait model is required.")
        selected <- backend$breeding_release_choices(learners, names(model$models),
                                                     "saved learners", allow_empty = FALSE)
        model$models <- model$models[selected]
      }
      imported <- backend$breeding_read_query(query, markers, model, input_memory_gib)
      metadata <- backend$breeding_input_metadata(sample_metadata)
      if (is.null(metadata)) metadata <- data.frame(sample_id = rownames(imported$genotype))
      backend$breeding_predict_models(model, imported$genotype, imported$markers, metadata, output, resources)
    }, list())
  }, envir = destination)
  assign("score_breeding_predictions", function(predictions, phenotype, assessment_ids = NULL,
                                                 selection_fraction = .10, direction = "higher") {
    run(function() {
      input <- backend$breeding_prediction_evidence(predictions, assessment_ids)
      if (length(intersect(input$ids, input$training_ids))) stop("Holdout IDs overlap the prediction model training panel.")
      backend$breeding_score_predictions(input$predictions, phenotype, input$ids, selection_fraction, direction)
    }, list())
  }, envir = destination)
  assign("run_breeding_gwas", function(...) run(backend$breeding_run_release_gwas, list(...)), envir = destination)
  assign("run_breeding_conditional_followup", function(...) run(backend$breeding_run_conditional_followup, list(...)), envir = destination)
  assign("diagnose_breeding_conditional_support", function(...) run(backend$breeding_diagnose_conditional_support, list(...)), envir = destination)
  assign("prepare_breeding_conditional_inputs", function(...) run(backend$breeding_prepare_conditional_inputs, list(...)), envir = destination)
  assign("export_breeding_imputed", function(...) run(backend$breeding_export_imputed, list(...)), envir = destination)
  assign("run_breeding_trait", function(...) run(backend$breeding_run_trait, list(...)), envir = destination)
  assign("run_breeding_traits", function(...) run(backend$breeding_run_traits, list(...)), envir = destination)
  assign("load_breeding_model", function(path) run(backend$breeding_load_model, list(path)), envir = destination)
  assign("explain_breeding_models", function(...) run(backend$breeding_explain_release_models, list(...)), envir = destination)
  assign("plot_breeding_report", function(...) run(backend$breeding_plot_report, list(...)), envir = destination)
  assign("validate_breeding_imputation", function(...) run(backend$breeding_validate_imputation, list(...)), envir = destination)
})
