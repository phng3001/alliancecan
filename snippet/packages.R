# Ensure BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Define packages
bioc_pkgs <- c(
  "DESeq2",
  "apeglm",
  "clusterProfiler",
)

cran_pkgs <- c(
  "pheatmap",
  "ashr",
  "ggrepel",
  "tidyverse",
)

# Helper: Install missing packages once
install_if_missing <- function(pkgs, install_fun, pkg_type = "CRAN") {
  installed <- rownames(installed.packages())
  to_install <- setdiff(pkgs, installed)
  
  if (length(to_install) > 0) {
    message(sprintf("Installing %s package(s): %s", pkg_type, paste(to_install, collapse = ", ")))
    install_fun(to_install)
  } else {
    message(sprintf("All %s packages are already installed.", pkg_type))
  }
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
