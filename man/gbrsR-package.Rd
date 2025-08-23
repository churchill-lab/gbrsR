\name{gbrsR-package}
\alias{gbrsR-package}
\alias{gbrsR}
\docType{package}
\title{Genome-Based Reconstruction of Strains Visualization}
\description{
A comprehensive R package for visualizing genome reconstructions from GBRS 
(Genome-Based Reconstruction of Strains) data. Supports NPZ files containing 
genotypic probability matrices and generates publication-ready genome plots 
with configurable parameters.
}
\details{
This package provides functions to:
\itemize{
  \item Load and process NPZ files containing genotypic probability data
  \item Analyze founder strain contributions across the genome
  \item Count recombination events per chromosome
  \item Generate publication-quality genome plots with customizable styling
  \item Save plots in various formats (PDF, PNG, SVG)
}

The main workflow functions are:
\itemize{
  \item \code{\link{plot_gbrs_genome}}: Main function for creating genome plots
  \item \code{\link{load_npz_file}}: Load NPZ files with smart fallback
  \item \code{\link{process_genotype_data}}: Process genotype probabilities
  \item \code{\link{create_genome_plot}}: Generate publication-quality plots
  \item \code{\link{save_genome_plot}}: Save plots in various formats
}

Configuration and customization:
\itemize{
  \item \code{\link{configure_plot}}: Customize plot parameters
  \item \code{\link{create_compact_config}}: Pre-configured compact layout
  \item \code{\link{DEFAULT_CONFIG}}: Default configuration parameters
  \item \code{\link{DEFAULT_FOUNDER_COLORS}}: Standard founder color schemes
}

For more information, see the vignette:
\code{vignette("introduction", package = "gbrsR")}
}
\author{
\strong{Maintainer}: Matthew Vincent \email{mvincent@jax.org}
}
\references{
For more information on GBRS and genome reconstruction:
\itemize{
  \item Diversity Outbred (DO) mouse studies
  \item Collaborative Cross (CC) mouse studies
  \item Genome-based strain reconstruction methods
}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/churchill-lab/gbrsR}
\item Report bugs at \url{https://github.com/churchill-lab/gbrsR/issues}
}
}
\examples{
\dontrun{
# Basic usage
library(gbrsR)
plot <- plot_gbrs_genome("example.npz")
save_genome_plot(plot, "genome_plot.pdf")

# With custom configuration
config <- configure_plot(compact_mode = TRUE)
plot <- plot_gbrs_genome("example.npz", config = config)
}
}
