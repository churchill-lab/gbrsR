#' Load NPZ genotype probability files
#'
#' @description
#' Loads NPZ files containing genotypic probability matrices and converts them
#' to R data frames. This function handles the critical dependency on Python
#' and NumPy through the reticulate package.
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
#' - "Error: Please install reticulate"
#' - "FileNotFoundError: No such file or directory"
#'
#' See the README for complete installation instructions.
#'
#' @examples
#' \dontrun{
#' # Load a standard NPZ file
#' data <- load_npz_file("example.npz")
#'
#' # Load with custom founders
#' data <- load_npz_file("example.npz", founders = c("A", "B", "C"))
#' }
#'
#' @export
load_npz_file <- function(npz_file,
                          founders = DEFAULT_FOUNDERS) {
    # generate diplotypes
    diplotypes <- generate_diplotypes(founders)

    # load data
    message("Loading NPZ file: ", basename(npz_file), "...")

    # Expand file path to handle tilde and relative paths
    npz_file <- normalizePath(npz_file, mustWork = FALSE)

    if (!file.exists(npz_file)) {
        stop(
            "NPZ file not found: ", npz_file, "\n",
            "Please check the file path and ensure the file exists."
        )
    }

    # Try to load NPZ file using reticulate (Python/numpy)
    if (requireNamespace("reticulate", quietly = TRUE)) {
        tryCatch(
            {
                message("Attempting to load with reticulate (Python/numpy)...")
                np <- reticulate::import("numpy")
                npz <- np$load(npz_file, allow_pickle = TRUE)
                chrom_names <- as.character(npz$files)

                # Extract data for each chromosome
                data_list <- list()
                for (chrom in chrom_names) {
                    data_list[[chrom]] <- npz$f[[chrom]]
                }

                message("Successfully loaded with reticulate")
                loaded_data <- list(
                    data = data_list,
                    chrom_names = chrom_names,
                    method = "reticulate",
                    success = TRUE
                )
            },
            error = function(e) {
                message("Reticulate failed: ", e$message)
                message("Run `reticulate::py_last_error()` for details.")
                stop(
                    "Failed to load NPZ file with reticulate. Please check:\n",
                    "1. File format is valid NPZ\n",
                    "2. File is not corrupted\n",
                    "3. Python dependencies are properly installed"
                )
            }
        )
    } else {
        stop(
            "Failed to load NPZ file. Please ensure:\n",
            '1. reticulate R package is installed: install.packages("reticulate")\n',
            "2. Python is installed and accessible\n",
            "3. NumPy is installed: pip install numpy\n",
            "See README for complete installation instructions."
        )
    }

    chrom_names <- loaded_data$chrom_names
    data_list <- loaded_data$data

    message("NPZ file loaded successfully using ", loaded_data$method, ": ", basename(npz_file))
    message("  - Found chromosomes: ", paste(chrom_names, collapse = ", "))
    message("  - Generated diplotypes: ", paste(diplotypes, collapse = ", "))

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

    message("Total markers loaded: ", nrow(df))
    return(df)
}
