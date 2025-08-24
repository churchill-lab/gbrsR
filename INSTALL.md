# Installing gbrsR

!!! **CRITICAL DEPENDENCY WARNING** !!!

**This package REQUIRES Python and the `reticulate` R package to function!**

**Without these, you will get errors and the package will not work.**

**See the detailed requirements below.**

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

### R Packages (Required)
The package requires these R packages:
- ggplot2
- dplyr  
- tidyr
- stringr
- **reticulate** !!! **CRITICAL**
- readr
- magrittr

**Note**: These R packages will be installed automatically when you install gbrsR.

### Python Requirements (MANDATORY)
!!! **THIS PACKAGE WILL NOT WORK WITHOUT PYTHON!**

You **MUST** have:
- **Python 3.6+** installed and accessible from R
- **NumPy package** (`pip install numpy` or `conda install numpy`)

**If you get errors about `reticulate` or Python, install these first!**
