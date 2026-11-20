# Open Breeding_pipeline.Rproj, adjust these settings, then click Source.

workers <- 2L
memory_gib <- 2
run_config <- "config/breeding_run.yml"
output_folder <- NULL  # Use the configured folder, or set a new one for this run.

local({
  if (.Platform$OS.type != "windows" || R.version$arch != "x86_64" ||
      getRversion() < "4.6.0" || getRversion() >= "4.7.0") {
    stop("Select 64-bit R 4.6.x in RStudio for this first-run launcher.", call. = FALSE)
  }
  if (!file.exists("breeding_pipeline.R") || !file.exists("scripts/run_breeding_home.R")) {
    stop("Open Breeding_pipeline.Rproj and run this script from that project.", call. = FALSE)
  }
  if (!is.numeric(workers) || length(workers) != 1L || !is.finite(workers) ||
      workers < 1 || workers > 1024 || workers != floor(workers)) {
    stop("Set workers to a whole number between 1 and 1024.", call. = FALSE)
  }
  if (!is.numeric(memory_gib) || length(memory_gib) != 1L || !is.finite(memory_gib) || memory_gib < 1) {
    stop("Set memory_gib to a number of at least 1.", call. = FALSE)
  }
  if (!is.character(run_config) || length(run_config) != 1L || is.na(run_config) || !file.exists(run_config)) {
    stop("Set run_config to an existing YAML file.", call. = FALSE)
  }
  if (!is.null(output_folder) && (!is.character(output_folder) || length(output_folder) != 1L ||
      is.na(output_folder) || !nzchar(output_folder))) {
    stop("Set output_folder to NULL or a new project-relative folder.", call. = FALSE)
  }

  project <- normalizePath(".", winslash = "/")
  scratch <- file.path(project, ".tmp")
  user_dir <- file.path(project, ".r_user")
  library_dir <- file.path(project, ".breeding-library", "R-4.6")
  for (folder in c(scratch, user_dir, library_dir)) {
    if (!dir.exists(folder) && !dir.create(folder, recursive = TRUE)) stop("Could not create ", folder)
  }
  executable <- file.path(R.home(), "bin", "R.exe")
  # Fresh R processes share a project-local library with the PowerShell launcher.
  # Leave the active RStudio session's environment and loaded packages alone.
  settings <- c(TEMP = scratch, TMP = scratch, TMPDIR = scratch, R_USER = user_dir,
                R_LIBS_USER = library_dir, R_LIBS_SITE = library_dir, R_LIBS = "",
                GWAS_R_EXECUTABLE = executable, OMP_NUM_THREADS = "1",
                OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1")
  previous <- Sys.getenv(names(settings), unset = NA_character_)
  on.exit({
    Sys.unsetenv(names(previous)[is.na(previous)])
    do.call(Sys.setenv, as.list(previous[!is.na(previous)]))
  })
  do.call(Sys.setenv, as.list(settings))

  run <- function(script, stage, arguments = character()) {
    options <- c("--vanilla", "--slave", "-f", shQuote(script))
    if (length(arguments)) options <- c(options, "--args", shQuote(arguments))
    cat("\n", stage, "\n", sep = "")
    status <- system2(executable, options)
    if (!identical(as.integer(status), 0L)) {
      stop(stage, " failed (exit ", status, "). See the messages above.", call. = FALSE)
    }
  }
  run("scripts/setup_breeding_home.R", "Checking R packages")
  run("scripts/run_breeding_home.R", "Running the breeding pipeline",
      c(run_config, as.integer(workers), format(memory_gib, scientific = FALSE, trim = TRUE, decimal.mark = "."),
        output_folder))
})
