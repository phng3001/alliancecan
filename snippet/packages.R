# Ensure BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Define packages (can be empty, length 1, or many)
bioc_pkgs <- c(
  "DESeq2",
  "apeglm",
  "clusterProfiler"
)

cran_pkgs <- c(
  "pheatmap",
  "ashr",
  "ggrepel",
  "tidyverse"
)

# Helper: Install missing packages once
install_if_missing <- function(pkgs, install_fun, pkg_type = "CRAN") {

  # Handle NULL, empty, or invalid input
  if (is.null(pkgs) || length(pkgs) == 0) {
    message(sprintf("No %s packages specified. Skipping.", pkg_type))
    return(invisible(NULL))
  }

  pkgs <- as.character(pkgs)
  pkgs <- pkgs[nzchar(pkgs)]  # remove empty strings

  if (length(pkgs) == 0) {
    message(sprintf("No valid %s packages specified. Skipping.", pkg_type))
    return(invisible(NULL))
  }

  installed <- rownames(installed.packages())
  to_install <- setdiff(pkgs, installed)

  if (length(to_install) > 0) {
    message(sprintf(
      "Installing %s package(s): %s",
      pkg_type,
      paste(to_install, collapse = ", ")
    ))
    install_fun(to_install)
  } else {
    message(sprintf("All %s packages are already installed.", pkg_type))
  }

  invisible(NULL)
}

# Install Bioconductor packages
install_if_missing(
  bioc_pkgs,
  function(pkgs) BiocManager::install(pkgs, update = FALSE, ask = FALSE),
  pkg_type = "Bioconductor"
)

# Install CRAN packages
install_if_missing(
  cran_pkgs,
  install.packages,
  pkg_type = "CRAN"
)

# Load all packages
# all_pkgs <- c(bioc_pkgs, cran_pkgs)
# invisible(lapply(all_pkgs, library, character.only = TRUE))
