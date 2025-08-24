#' GBRS Mosaic Plotter
#'
#' @description
#' A comprehensive R package for visualizing genome reconstructions from GBRS
#' (Genome-Based Reconstruction of Strains) data. Supports NPZ files containing
#' genotypic probability matrices and generates publication-ready genome plots
#' with configurable parameters.
#'
#' @details
#' This package provides functions to:
#' \itemize{
#'   \item Load and process NPZ files containing genotypic probability data
#'   \item Analyze founder strain contributions across the genome
#'   \item Count recombination events per chromosome
#'   \item Generate publication-quality genome plots with customizable styling
#'   \item Save plots in various formats (PDF, PNG, SVG)
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import readr
#' @importFrom magrittr %>%
#'
#' @name "_PACKAGE"
#' @title gbrsR: GBRS for R
#' @description
#' A comprehensive R package for analyzing and visualizing genetic recombination
#' data from Diversity Outbred (DO) and Collaborative Cross (CC) mouse studies.
#' Provides tools for loading NPZ genotype probability files, processing
#' recombination data, and creating publication-ready genome plots.
#'
#' @section CRITICAL DEPENDENCY WARNING:
#' This package REQUIRES Python and the reticulate R package to function!
#' Without these dependencies, you will get errors like:
#' - "Error: Please install reticulate"
#' - "FileNotFoundError: No such file or directory"
#'
#' See the README for complete installation instructions.
#'
#' @details
#' The package includes functions for:
#' \itemize{
#'   \item Loading and processing NPZ genotype probability files
#'   \item Analyzing recombination patterns across chromosomes
#'   \item Creating customizable genome plots with founder strain information
#'   \item Generating compact and publication-ready visualizations
#' }
#'
#' @seealso
#' \code{\link{plot_gbrs_genome}} for the main plotting function
#' \code{\link{load_npz_file}} for loading NPZ files
#' \code{\link{configure_plot}} for customizing plot appearance
NULL
