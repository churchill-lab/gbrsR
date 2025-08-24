# gbrsR: Genome-Based Reconstruction of Strains Visualization

[![R-CMD-check](https://github.com/churchill-lab/gbrsR/workflows/R-CMD-check/badge.svg)](https://github.com/churchill-lab/gbrsR/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/gbrsR)](https://cran.r-project.org/package=gbrsR)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive R package for visualizing genome reconstructions from GBRS (Genome-Based Reconstruction of Strains) data. Supports NPZ files containing genotypic probability matrices and generates publication-ready genome plots with configurable parameters.

---

🚨 **CRITICAL DEPENDENCY WARNING** 🚨

**This package REQUIRES Python and the `reticulate` R package to function!**

**Without these, you will get errors like:**
- `Error: Please install reticulate`
- `FileNotFoundError: No such file or directory`

**See the [Installation](#installation) section below for setup instructions.**

---

## Features

- **NPZ File Support**: Load and process NPZ files containing genotypic probability data
- **Founder Strain Analysis**: Analyze founder strain contributions across the genome
- **Recombination Detection**: Count and visualize recombination events per chromosome
- **Publication Quality**: Generate high-quality plots suitable for publication
- **Flexible Configuration**: Customize plot appearance, spacing, and styling
- **Multiple Formats**: Export plots as PDF, PNG, SVG, or TIFF
- **Compact Mode**: Pre-configured settings for space-constrained layouts

## Installation

⚠️ **CRITICAL REQUIREMENT**: This package **REQUIRES** `reticulate` and Python to function!

### Prerequisites (MANDATORY)

```r
# Install required packages (reticulate is ESSENTIAL)
install.packages(c("ggplot2", "dplyr", "tidyr", "stringr", "reticulate", "readr", "magrittr"))

# Install suggested packages (optional but recommended)
install.packages(c("testthat", "knitr", "rmarkdown"))
```

### From GitHub (Recommended)

```r
# Install devtools if you don't have it
if (!require(devtools)) install.packages("devtools")

# Install from GitHub
devtools::install_github("churchill-lab/gbrsR")
```

### Alternative: Using remotes

```r
# Install remotes if you don't have it
if (!require(remotes)) install.packages("remotes")

# Install from GitHub
remotes::install_github("churchill-lab/gbrsR")
```

### Python Requirements (MANDATORY)

🚨 **THIS PACKAGE WILL NOT WORK WITHOUT PYTHON!**

This package uses `reticulate` to interface with Python. You **MUST** have:
- **Python 3.6+** installed and accessible from R
- **NumPy package** (`pip install numpy` or `conda install numpy`)
- **The `reticulate` R package** (installed above)

**If you get errors about `reticulate` or Python, install these first!**

## Troubleshooting

### Common Error: "Please install reticulate"

If you see this error:
```
Error in load_npz_file(npz_file) : Please install reticulate
```

**SOLUTION**: Install the required packages first:
```r
# Install reticulate and other dependencies
install.packages(c("reticulate", "ggplot2", "dplyr", "tidyr", "stringr", "readr", "magrittr"))

# Then reinstall gbrsR
devtools::install_github("churchill-lab/gbrsR")
```

### Common Error: "FileNotFoundError" or "No such file or directory"

If you see Python errors about missing files:
```
FileNotFoundError: [Errno 2] No such file or directory
```

**SOLUTION**: Ensure Python and NumPy are installed:
```bash
# Check Python version
python3 --version

# Install NumPy
pip3 install numpy
# OR if using conda
conda install numpy
```

## Quick Start

```r
library(gbrsR)

# Load and plot a genome reconstruction
plot <- plot_gbrs_genome("path/to/your/file.npz")

# Save the plot
save_genome_plot(plot, "genome_plot.pdf")
```

## Usage Examples

### Basic Genome Plot

```r
# Create a basic genome plot
plot <- plot_gbrs_genome("example.npz", sample_name = "DO140")
save_genome_plot(plot, "DO140_genome.pdf")
```

### Custom Configuration

```r
# Create a compact configuration for publication
config <- configure_plot(
    font_size_recomb = 6,
    chrom_spacing = 0.2,
    bar_height = 0.1,
    founder_gap = 0.02,
    font_size_title = 16,
    compact_mode = TRUE
)

# Generate plot with custom settings
plot <- plot_gbrs_genome("example.npz", config = config)
save_genome_plot(plot, "compact_genome.pdf", width = 14, height = 10)
```

### Pre-configured Compact Layout

```r
# Use built-in compact configuration
config <- create_compact_config()
plot <- plot_gbrs_genome("example.npz", config = config)
save_genome_plot(plot, "publication_ready.pdf")
```

### Return Both Plot and Data

```r
# Get both the plot and processed data
result <- plot_gbrs_genome("example.npz", return_data = TRUE)

# Access components
plot <- result$plot
data <- result$data
recombinations <- result$recombinations

# Save plot
save_genome_plot(plot, "genome_with_data.pdf")
```

## Input Data Format

The package expects NPZ files containing:
- **Genotypic probability matrices** for each chromosome
- **Diplotype probabilities** (e.g., AA, AB, AC, etc.)
- **Marker positions** across the genome

## Configuration Options

The package provides extensive customization through the `configure_plot()` function:

- **Layout**: `bar_height`, `founder_gap`, `chrom_spacing`
- **Fonts**: `font_size_title`, `font_size_axis`, `font_size_tick`
- **Colors**: `background_color`, `panel_border`, `grid_alpha`
- **Spacing**: `y_scale_factor`, `plot_margin`, `title_margin`

## Dependencies

- **Required**: ggplot2, dplyr, tidyr, stringr, reticulate, readr, magrittr
- **Suggested**: testthat, knitr, rmarkdown, extrafont, systemfonts

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

