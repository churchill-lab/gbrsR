#' Load NPZ genotype probability files (reticulate backend)
#'
#' @description
#' Loads NPZ files using Python/NumPy via the reticulate package and converts
#' them into a unified R data frame of diplotype probabilities with metadata.
#'
#' @param npz_file Path to the NPZ file to load
#' @param founders Character vector of founder strain names (default: DEFAULT_FOUNDERS)
#'
#' @return Data frame with columns for each diplotype and metadata columns
#'
#' @details
#' This function performs several key operations:
#' \itemize{
#'   \item{Loads NPZ file using Python/NumPy via reticulate}
#'   \item{Generates all possible founder strain combinations (diplotypes)}
#'   \item{Converts probability matrices to R data frames}
#'   \item{Adds chromosome and marker index information}
#' }
#'
#' @section CRITICAL DEPENDENCY WARNING:
#' This function REQUIRES Python and the reticulate R package to function!
#' Without these dependencies, you will get errors like:
#' - 'Error: Please install reticulate'
#' - 'FileNotFoundError: No such file or directory'
#'
#' See the README for complete installation instructions.
#'
#' @examples
#' \dontrun{
#' # Load a standard NPZ file via reticulate
#' data <- load_npz_file_reticulate('example.npz')
#'
#' # Load with custom founders
#' data <- load_npz_file_reticulate('example.npz', founders = c('A', 'B', 'C'))
#' }
load_npz_file_reticulate <- function(npz_file,
                                    founders = DEFAULT_FOUNDERS) {
    # generate diplotypes
    diplotypes <- generate_diplotypes(founders)

    # load data
    message('Loading NPZ file: ', basename(npz_file), '...')

    # Expand file path to handle tilde and relative paths
    npz_file <- normalizePath(npz_file, mustWork = FALSE)

    if (!file.exists(npz_file)) {
        stop(
            'NPZ file not found: ', npz_file, '\n',
            'Please check the file path and ensure the file exists.'
        )
    }

    # Try to load NPZ file using reticulate (Python/numpy)
    if (requireNamespace('reticulate', quietly = TRUE)) {
        tryCatch(
            {
                message('Attempting to load with reticulate (Python/numpy)...')
                np <- reticulate::import('numpy')
                npz <- np$load(npz_file, allow_pickle = TRUE)
                chrom_names <- as.character(npz$files)

                # Extract data for each chromosome
                data_list <- list()
                for (chrom in chrom_names) {
                    data_list[[chrom]] <- npz$f[[chrom]]
                }

                message('Successfully loaded with reticulate')
                loaded_data <- list(
                    data = data_list,
                    chrom_names = chrom_names,
                    method = 'reticulate',
                    success = TRUE
                )
            },
            error = function(e) {
                message('Reticulate failed: ', e$message)
                message('Run `reticulate::py_last_error()` for details.')
                stop(
                    'Failed to load NPZ file with reticulate. Please check:\n',
                    '1. File format is valid NPZ\n',
                    '2. File is not corrupted\n',
                    '3. Python dependencies are properly installed'
                )
            }
        )
    } else {
        stop(
            'Failed to load NPZ file. Please ensure:\n',
            '1. reticulate R package is installed: install.packages("reticulate")\n',
            '2. Python is installed and accessible\n',
            '3. NumPy is installed: pip install numpy\n',
            'See README for complete installation instructions.'
        )
    }

    chrom_names <- loaded_data$chrom_names
    data_list <- loaded_data$data

    message('NPZ file loaded successfully using ', loaded_data$method, ': ', basename(npz_file))
    message('  - Found chromosomes: ', paste(chrom_names, collapse = ', '))
    message('  - Generated diplotypes: ', paste(diplotypes, collapse = ', '))

    # process each chromosome
    all_df <- list()
    for (chrom in chrom_names) {
        mat <- data_list[[chrom]]
        n_markers <- ncol(mat)

        # create data frame for this chromosome
        df_chr <- as.data.frame(t(mat)) # transpose to match expected format
        colnames(df_chr) <- diplotypes
        df_chr$chromosome <- chrom
        df_chr$marker_index <- 0:(n_markers - 1)

        all_df[[chrom]] <- df_chr
    }

    # combine all chromosomes
    df <- dplyr::bind_rows(all_df)

    message('Total markers loaded: ', nrow(df))
    return(df)
}


#' Load NPZ genotype probability files (RcppCNPy backend)
#'
#' @description
#' Loads NPZ files by unzipping the archive and reading each contained .npy
#' array using RcppCNPy, then converts them into a unified R data frame of
#' diplotype probabilities with metadata.
#'
#' @param npz_file Path to the NPZ file to load
#' @param founders Character vector of founder strain names (default: DEFAULT_FOUNDERS)
#' @param cleanup Logical; remove temporary extracted files (default: TRUE)
#'
#' @return Data frame with columns for each diplotype and metadata columns
#'
#' @details
#' This function:
#' \itemize{
#'   \item Unzips the .npz archive to a temporary directory
#'   \item Reads each chromosome's .npy array via \code{RcppCNPy::npyLoad}
#'   \item Generates diplotype column names from \code{founders}
#'   \item Returns a combined data frame with one row per marker per chromosome
#' }
#'
#' @section Dependencies:
#' Requires the \code{RcppCNPy} package. No Python installation is needed.
#'
#' @examples
#' \dontrun{
#' # Load a standard NPZ file via RcppCNPy
#' data <- load_npz_file_rcpp_cnpy('example.npz')
#'
#' # Load with custom founders
#' data <- load_npz_file_rcpp_cnpy('example.npz', founders = c('A', 'B', 'C'))
#' }
load_npz_file_rcpp_cnpy <- function(npz_file,
                                    founders = DEFAULT_FOUNDERS,
                                    cleanup = TRUE) {

    if (!requireNamespace('RcppCNPy', quietly = TRUE)) {
        stop('Please install RcppCNPy: install.packages("RcppCNPy")')
    }

    # load data
    message('Loading NPZ file: ', basename(npz_file), '...')

    npz_file <- normalizePath(npz_file, mustWork = FALSE)
    if (!file.exists(npz_file)) {
        stop('NPZ file not found: ', npz_file)
    }

    diplotypes <- gbrsR::generate_diplotypes(founders)

    # Extract .npz (zip) contents to a temporary directory
    exdir <- file.path(tempdir(), paste0('rcppcnpy_', tools::file_path_sans_ext(basename(npz_file)), '_', as.integer(Sys.time())))

    dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
    message('Attempting to load with RcppCNPy...')
    utils::unzip(npz_file, exdir = exdir)

    # Collect all .npy files (keys inside the npz)
    npy_files <- list.files(exdir, pattern = '\\.npy$', full.names = TRUE, recursive = TRUE)
    if (length(npy_files) == 0) {
        if (cleanup && dir.exists(exdir)) unlink(exdir, recursive = TRUE, force = TRUE)
        stop('No .npy entries found after unzipping NPZ: ', npz_file)
    }
    message('Successfully loaded with RcppCNPy')

    # Map chromosome name from file basename (strip .npy)
    chrom_names <- sub('\\.npy$', '', basename(npy_files))

    message('NPZ file loaded successfully using RcppCNPy: ', basename(npz_file))
    message('  - Found chromosomes: ', paste(chrom_names, collapse = ', '))
    message('  - Generated diplotypes: ', paste(diplotypes, collapse = ', '))

    # Load each .npy and build per-chromosome data frames identical to gbrsR::load_npz_file
    all_df <- vector('list', length(npy_files))
    names(all_df) <- chrom_names

    for (i in seq_along(npy_files)) {
        chrom <- chrom_names[i]
        # Use dotranspose=FALSE to match numpy orientation (validated against reticulate)
        mat <- RcppCNPy::npyLoad(npy_files[i], type = 'numeric', dotranspose = FALSE)

        if (!is.matrix(mat)) {
            mat <- as.matrix(mat)
        }

        n_markers <- ncol(mat)
        df_chr <- as.data.frame(t(mat))
        colnames(df_chr) <- diplotypes
        df_chr$chromosome <- chrom
        df_chr$marker_index <- 0:(n_markers - 1)

        all_df[[i]] <- df_chr
    }

    # Combine
    df <- dplyr::bind_rows(all_df)

    # Optional cleanup of extracted directory
    if (cleanup && dir.exists(exdir)) unlink(exdir, recursive = TRUE, force = TRUE)

    message('Total markers loaded: ', nrow(df))
    return(df)
}

#' Load NPZ genotype probability files (with backend fallback)
#'
#' @description
#' High-level loader for NPZ genotype probability matrices that returns a
#' unified R data frame. It prefers the reticulate (Python/NumPy) backend when
#' available; if \code{reticulate} is not installed, it falls back to the
#' RcppCNPy-based loader. If neither backend is available, it errors.
#'
#' @param npz_file Path to the NPZ file to load
#' @param founders Character vector of founder strain names (default: DEFAULT_FOUNDERS)
#' @param cleanup Logical; remove temporary extracted files (default: TRUE)
#'
#' @return Data frame with diplotype probability columns and metadata columns
#'         \code{chromosome} and \code{marker_index}.
#'
#' @details
#' Behavior:
#' \enumerate{
#'   \item If \code{reticulate} is installed, uses the reticulate backend
#'   \item Else if \code{RcppCNPy} is installed, uses the RcppCNPy backend
#'   \item Else, throws an informative error with install guidance
#' }
#'
#' @examples
#' \dontrun{
#' # Load a standard NPZ file (uses reticulate if available, else RcppCNPy)
#' data <- load_npz_file('example.npz')
#'
#' # Load with custom founders
#' data <- load_npz_file('example.npz', founders = c('A', 'B', 'C'))
#' }
#'
#' @export
load_npz_file <- function(npz_file,
                          founders = DEFAULT_FOUNDERS,
                          cleanup = TRUE) {
    if (requireNamespace('reticulate', quietly = TRUE)) {
        return(load_npz_file_reticulate(npz_file, founders))
    } else if (requireNamespace('RcppCNPy', quietly = TRUE)) {
        return(load_npz_file_rcpp_cnpy(npz_file, founders, cleanup))
    } else {
        stop('Please install reticulate or RcppCNPy to load NPZ files')
    }
}
