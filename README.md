# gbrsR: Genome-Based Reconstruction of Strains Visualization

[![R-CMD-check](https://github.com/churchill-lab/gbrsR/workflows/R-CMD-check/badge.svg)](https://github.com/churchill-lab/gbrsR/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/gbrsR)](https://cran.r-project.org/package=gbrsR)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive R package for visualizing genome reconstructions from GBRS (Genome-Based Reconstruction of Strains) data. Supports NPZ files containing genotypic probability matrices and generates publication-ready genome plots with configurable parameters.

## Features

- **NPZ File Support**: Load and process NPZ files containing genotypic probability data
- **Founder Strain Analysis**: Analyze founder strain contributions across the genome
- **Recombination Detection**: Count and visualize recombination events per chromosome
- **Publication Quality**: Generate high-quality plots suitable for publication
- **Flexible Configuration**: Customize plot appearance, spacing, and styling
- **Multiple Formats**: Export plots as PDF, PNG, SVG, or TIFF
- **Compact Mode**: Pre-configured settings for space-constrained layouts

## Installation

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
- **Suggested**: testthat, knitr, rmarkdown, svglite

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use gbrsR in your research, please cite:

```
@software{gbrsR,
  title = {gbrsR: Genome-Based Reconstruction of Strains Visualization},
  author = {Matthew Vincent},
  year = {2024},
  url = {https://github.com/churchill-lab/gbrsR}
}
```

## Support

- **Issues**: [GitHub Issues](https://github.com/churchill-lab/gbrsR/issues)
- **Documentation**: [Package Documentation](https://churchill-lab.github.io/gbrsR/)
- **Email**: mvincent@jax.org

## Acknowledgments

- Diversity Outbred (DO) and Collaborative Cross (CC) mouse research community
- R and tidyverse development teams
- Python/numpy community for NPZ file format support
