# Install runtime dependencies in the library assigned by run_pipeline.ps1.
if (.Platform$OS.type != "windows" || R.version$arch != "x86_64" ||
    getRversion() < "4.6.0" || getRversion() >= "4.7.0") {
  stop("The Windows first-run installer requires 64-bit R 4.6.x.", call. = FALSE)
}
library_dir <- normalizePath(".breeding-library/R-4.6", winslash = "/", mustWork = TRUE)
if (!identical(library_dir, normalizePath(Sys.getenv("R_LIBS_USER"), winslash = "/", mustWork = TRUE))) {
  stop("Start the dependency setup with Run_pipeline.R in RStudio or run_pipeline.ps1 in PowerShell.", call. = FALSE)
}
.libPaths(c(library_dir, .Library), include.site = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 1L,
        timeout = max(300, getOption("timeout")), BioC_mirror = "https://bioconductor.org",
        install.packages.compile.from.source = "never")

cran <- c("data.table", "digest", "yaml", "ggplot2", "glmnet", "ranger", "xgboost", "readxl",
          "ape", "bigmemory", "EMMREML", "genetics", "gplots", "htmltools", "magrittr", "lme4",
          "plotly", "RcppArmadillo", "scatterplot3d", "snowfall")
bioc <- c("multtest", "snpStats")
installed <- function(packages) packages %in% rownames(utils::installed.packages(lib.loc = c(library_dir, .Library)))
missing <- cran[!installed(cran)]
if (length(missing)) utils::install.packages(missing, lib = library_dir, type = "binary", dependencies = NA)
missing <- bioc[!installed(bioc)]
if (length(missing)) utils::install.packages(missing, lib = library_dir, type = "binary", dependencies = NA,
  repos = c(BioCsoft = "https://bioconductor.org/packages/3.23/bioc", getOption("repos")))

prerequisites <- c(cran, bioc, "compiler", "grid", "MASS")
unavailable <- prerequisites[!vapply(prerequisites, requireNamespace, logical(1L), quietly = TRUE)]
if (length(unavailable)) stop("Package setup incomplete. These namespaces could not load: ",
                              paste(unavailable, collapse = ", "), call. = FALSE)

gapit_path <- find.package("GAPIT", lib.loc = library_dir, quiet = TRUE)
if (!length(gapit_path) || utils::packageVersion("GAPIT", lib.loc = library_dir) != "4.1.0") {
  source_file <- "vendor/GAPIT_4.1.0.tar.gz"
  if (!file.exists(source_file)) {
    source_file <- file.path(tempdir(), "GAPIT_4.1.0.tar.gz")
    url <- "https://api.github.com/repos/jiabowang/GAPIT/tarball/a3ffee737ac2ddc567dd2fc3a5fe3014d6bf64a9"
    utils::download.file(url, source_file, mode = "wb")
  }
  expected <- "3520262e4e8bc61a1e7ee6d7291d7aeb63b289149fa8c12c5f289257e1631b6d"
  if (!identical(digest::digest(file = source_file, algo = "sha256"), expected)) {
    stop("GAPIT source checksum differs. Use the vendor archive in the Windows starter bundle.", call. = FALSE)
  }
  utils::install.packages(source_file, lib = library_dir, repos = NULL, type = "source", dependencies = FALSE)
}
if (!requireNamespace("GAPIT", quietly = TRUE) || utils::packageVersion("GAPIT", lib.loc = library_dir) != "4.1.0" ||
    !"GAPIT" %in% getNamespaceExports("GAPIT")) stop("GAPIT 4.1.0 failed to load.", call. = FALSE)
cat("Breeding dependencies ready in", library_dir, "\n")
