# Installing gbrsR

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

The package requires these R packages:
- ggplot2
- dplyr  
- tidyr
- stringr
- reticulate
- readr
- magrittr

These will be installed automatically when you install gbrsR.
