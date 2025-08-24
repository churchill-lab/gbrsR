#' Generate all possible founder strain combinations
#'
#' @description
#' Creates all possible diploid combinations from founder strain names.
#' This is used internally by load_npz_file to generate column names.
#'
#' @param founders Character vector of founder strain names
#'
#' @return Character vector of all possible diplotype combinations
#'
#' @examples
#' generate_diplotypes(c('A', 'B', 'C'))
#' # Returns: c('AA', 'AB', 'AC', 'BA', 'BB', 'BC', 'CA', 'CB', 'CC')
#'
#' @export
generate_diplotypes <- function(founders = DEFAULT_FOUNDERS) {
    diplotypes <- c()
    for (i in seq_along(founders)) {
        for (j in i:length(founders)) {
            diplotypes <- c(diplotypes, paste0(founders[i], founders[j]))
        }
    }
    return(diplotypes)
}

#' Process genotype data to identify founder contributions
#'
#' @description
#' Analyzes the genotype probability data to determine the most likely
#' founder strain at each genomic position. This function converts probability
#' matrices into discrete founder assignments for visualization.
#'
#' @param df Input data frame with diplotype probability columns
#' @param founders Character vector of founder strain names (default: DEFAULT_FOUNDERS)
#' @param config Configuration list (default: DEFAULT_CONFIG)
#' @param chrom_lengths Named vector of chromosome lengths (default: DEFAULT_CHROM_LENGTHS)
#'
#' @return Processed data frame with founder assignments including:
#' \describe{
#'   \item{max_diplotype}{Most likely diplotype at each position}
#'   \item{founder1}{Second founder in the diplotype}
#'   \item{founder2}{First founder in the diplotype}
#'   \item{y_base}{Vertical position for plotting}
#'   \item{cM}{Centimorgan position}
#' }
#'
#' @details
#' This function performs several key operations:
#' \itemize{
#'   \item{Identifies diplotype probability columns in the data}
#'   \item{Determines the most likely diplotype at each marker}
#'   \item{Extracts individual founder assignments}
#'   \item{Orders chromosomes for proper plotting}
#'   \item{Calculates vertical positions and cM coordinates}
#' }
#'
#' @examples
#' \dontrun{
#' # Process genotype data
#' processed_data <- process_genotype_data(loaded_data)
#'
#' # Process with custom configuration
#' config <- configure_plot(compact_mode = TRUE)
#' processed_data <- process_genotype_data(loaded_data, config = config)
#'
#' # Process with custom chromosome ordering
#' custom_chroms <- c('1' = 1000000, '2' = 2000000, 'X' = 1500000)
#' processed_data <- process_genotype_data(loaded_data, chrom_lengths = custom_chroms)
#' }
#'
#' @export
process_genotype_data <- function(df,
                                  founders = DEFAULT_FOUNDERS,
                                  config = DEFAULT_CONFIG,
                                  chrom_lengths = DEFAULT_CHROM_LENGTHS) {
    # identify diplotype columns
    diplotype_cols <- grep('^[A-H][A-H]$', names(df), value = TRUE)

    if (length(diplotype_cols) == 0) {
        stop(
            'No diplotype columns found in data frame.\n',
            'Expected columns with names like AA, AB, AC, etc.\n',
            'Please check your data format and founder strain names.'
        )
    }

    message('Found ', length(diplotype_cols), ' diplotype columns')

    # for each marker, find the diplotype with the highest probability
    df <- df %>%
        dplyr::mutate(
            max_diplotype = diplotype_cols[apply(dplyr::select(df, dplyr::all_of(diplotype_cols)), 1, which.max)],
            founder2 = substr(.data$max_diplotype, 1, 1), # first founder in diplotype
            founder1 = substr(.data$max_diplotype, 2, 2) # second founder in diplotype
        )

    # order chromosomes for plotting (natural sort)
    chrom_order <- names(chrom_lengths)[names(chrom_lengths) %in% unique(df$chromosome)]
    chrom_order <- chrom_order[order(as.numeric(gsub('[^0-9]', '', chrom_order)))]
    df$chromosome <- factor(df$chromosome, levels = chrom_order)

    message('Chromosomes ordered for plotting: ', paste(chrom_order, collapse = ', '))

    # assign vertical positions for chromosomes
    num_chrs <- length(chrom_order)

    # validate configuration parameters
    if (config$bar_height <= 0) {
        warning('bar_height must be positive, setting to 0.1')
        config$bar_height <- 0.1
    }
    if (config$founder_gap < 0) {
        warning('founder_gap cannot be negative, setting to 0')
        config$founder_gap <- 0
    }

    # use provided chrom_spacing or calculate based on bar heights and gaps
    total_chrom_height <- config$bar_height + config$founder_gap + config$bar_height

    # for compact layouts, allow override of tight spacing if explicitly requested
    if (is.null(config$chrom_spacing)) {
        # auto-calculate spacing to accommodate bars and gaps
        chrom_spacing <- total_chrom_height + 0.5
    } else {
        # if user explicitly sets a very small chrom_spacing, use it for compact layouts
        if (config$chrom_spacing < total_chrom_height + 0.05 || config$compact_mode) {
            if (config$compact_mode) {
                message('Compact mode enabled: using minimal spacing')
            } else {
                message('Creating compact layout with minimal spacing')
            }
            chrom_spacing <- config$chrom_spacing
        } else {
            # use provided spacing, but prevent overlap
            min_required_spacing <- total_chrom_height + 0.1
            chrom_spacing <- max(config$chrom_spacing, min_required_spacing)

            if (config$chrom_spacing < min_required_spacing) {
                message('Warning: chrom_spacing (', config$chrom_spacing, ') too small, using minimum (', min_required_spacing, ')')
            }
        }
    }

    # apply y-axis scaling factor with safety bounds
    # ensure we have enough spacing for all chromosomes and positive values
    min_spacing <- max(total_chrom_height + 0.1, 0.5)
    scaled_spacing <- chrom_spacing * config$y_scale_factor
    final_spacing <- max(scaled_spacing, min_spacing)

    if (scaled_spacing < min_spacing) {
        message('Warning: scaled spacing too small, using minimum spacing')
    }

    # calculate chromosome positions with safety checks
    # ensure a reasonable position and use positive spacing
    start_pos <- max(num_chrs * final_spacing, num_chrs * 0.5)
    spacing_step <- max(final_spacing, 0.5)

    chr_ypos <- stats::setNames(seq(start_pos, spacing_step, by = -spacing_step), chrom_order)

    # ensure we have exactly the right number of positions
    if (length(chr_ypos) != length(chrom_order)) {
        # Fallback: create positions manually
        message('Creating fallback chromosome positions')
        chr_ypos <- stats::setNames(seq(num_chrs, 1), chrom_order)
    }

    # Final validation
    if (length(chr_ypos) != length(chrom_order)) {
        stop(
            'Critical error: Cannot create valid chromosome positions.\n',
            'chr_ypos length: ', length(chr_ypos), ', chrom_order length: ', length(chrom_order), '\n',
            'This may indicate corrupted data or configuration issues.\n',
            'Please check your input data and try again.'
        )
    }

    df <- df %>% dplyr::mutate(y_base = chr_ypos[as.character(.data$chromosome)])

    # calculate cM positions
    df <- df %>%
        dplyr::group_by(.data$chromosome) %>%
        dplyr::mutate(cM = .data$marker_index * config$grid_width) %>%
        dplyr::ungroup()

    message('Genotypic data processing completed successfully')
    return(df)
}

#' Count recombination events per chromosome (internal function)
#'
#' @description
#' Identifies and counts recombination events by detecting changes in
#' founder strain assignments between adjacent genomic positions. This
#' function is essential for quantifying genetic diversity and mapping
#' recombination hotspots.
#'
#' @section Internal Function:
#' This function is internal and not exported. Advanced users can access it
#' with \code{gbrsR:::count_recombinations}.
#'
#' @param df Processed data frame from process_genotype_data()
#'
#' @return Data frame with recombination counts containing:
#' \describe{
#'   \item{chromosome}{Chromosome identifier}
#'   \item{num_recomb}{Number of recombination events detected}
#' }
#'
#' @details
#' Recombination events are detected by comparing the diplotype at each
#' marker position with the previous position. A change in diplotype
#' indicates a recombination event has occurred.
#'
#' @examples
#' \dontrun{
#' # Count recombinations in processed data
#' recomb_data <- count_recombinations(processed_data)
#'
#' # View recombination summary
#' print(recomb_data)
#' }
#'
count_recombinations <- function(df) {
    recomb_counts <- df %>%
        dplyr::group_by(.data$chromosome) %>%
        dplyr::summarize(
            num_recomb = sum(dplyr::lag(.data$max_diplotype, default = dplyr::first(.data$max_diplotype)) != .data$max_diplotype),
            .groups = 'drop'
        )

    message('Recombination counts calculated for ', nrow(recomb_counts), ' chromosomes successfully')

    # Print recombination summary
    message('\nRecombination Summary')
    message('---------------------')
    total_recomb <- sum(recomb_counts$num_recomb)
    message('Total recombinations: ', total_recomb)
    message('Per chromosome breakdown:')
    for (i in seq_along(recomb_counts$chromosome)) {
        chr <- recomb_counts$chromosome[i]
        nrec <- recomb_counts$num_recomb[i]
        message('  - Chromosome ', chr, ': ', nrec, ' recombinations')
    }

    return(recomb_counts)
}

#' Prepare data for plotting (internal function)
#'
#' @description
#' Transforms the processed data into long format for ggplot2 visualization
#' with proper positioning. This function converts the wide-format data
#' into a long format suitable for creating stacked bar plots.
#'
#' @section Internal Function:
#' This function is internal and not exported. Advanced users can access it
#' with \code{gbrsR:::prepare_plot_data}.
#'
#' @param df Processed data frame from process_genotype_data()
#' @param config Configuration list (default: DEFAULT_CONFIG)
#'
#' @return Long-format data frame for plotting containing:
#' \describe{
#'   \item{chromosome}{Chromosome identifier}
#'   \item{marker_index}{Marker position index}
#'   \item{max_diplotype}{Most likely diplotype}
#'   \item{hap}{Haplotype identifier ('founder1' or 'founder2')}
#'   \item{founder}{Founder strain assignment}
#'   \item{y}{Vertical position for plotting}
#'   \item{yend}{End position for bar height}
#'   \item{cM}{Centimorgan position}
#' }
#'
#' @details
#' This function performs the following transformations:
#' \itemize{
#'   \item{Converts wide format (founder1, founder2 columns) to long format}
#'   \item{Calculates vertical positions for stacked bars}
#'   \item{Determines bar heights and positions}
#'   \item{Scales cM positions for proper x-axis representation}
#' }
#'
#' @examples
#' \dontrun{
#' # Prepare data for plotting
#' plot_data <- prepare_plot_data(processed_data)
#'
#' # Use custom configuration
#' config <- configure_plot(bar_height = 0.5)
#' plot_data <- prepare_plot_data(processed_data, config = config)
#' }
#'
prepare_plot_data <- function(df, config = DEFAULT_CONFIG) {
    # Create long format for plotting
    df_long <- df %>%
        tidyr::pivot_longer(cols = c('founder1', 'founder2'), names_to = 'hap', values_to = 'founder') %>%
        dplyr::mutate(
            y = ifelse(.data$hap == 'founder1', .data$y_base, .data$y_base + config$bar_height + config$founder_gap),
            yend = .data$y + config$bar_height,
            cM = .data$marker_index * (config$max_cM / max(df$marker_index))
        )

    message('Data transformed to long format for plotting successfully')
    return(df_long)
}
