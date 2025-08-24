# Installing gbrsR

!!! **NPZ LOADING BACKENDS** !!!

To load NPZ files, install at least one backend at runtime:

- reticulate (with Python + NumPy), or
- RcppCNPy

The package prefers reticulate if available; otherwise it falls back to RcppCNPy. If neither is installed, NPZ loading will error with guidance.

## From GitHub

### Option 1: Using devtools (Recommended)

```r
# Install devtools if you don't have it
if (!require(devtools)) install.packages("devtools")

# Install from GitHub
devtools::install_github("churchill-lab/gbrsR")
```

### Option 2: Using remotes

```r
# Install remotes if you don't have it
if (!require(remotes)) install.packages("remotes")

# Install from GitHub
remotes::install_github("churchill-lab/gbrsR")
```

## From Source

If you have the source code:

```r
# Build the package
system("R CMD build .")

# Install the built package
install.packages("gbrsR_1.0.0.tar.gz", repos = NULL, type = "source")
```

## Verify Installation

```r
library(gbrsR)
# Should load without errors
```

## Dependencies

### R Packages
- Required: ggplot2, dplyr, tidyr, stringr, readr, magrittr
- NPZ backends (install at least one): reticulate, RcppCNPy
- Suggested: testthat, knitr, rmarkdown, extrafont, systemfonts

### If using reticulate backend
- Python 3.6+ installed and accessible from R
- NumPy installed (`pip install numpy` or `conda install numpy`)
