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
#' @import reticulate
#' @import readr
#' @importFrom magrittr %>%
#'
#' @name "_PACKAGE"
#' @title gbrsR: Genome Browser for Recombination Studies
#' @description
#' A comprehensive R package for analyzing and visualizing genetic recombination
#' data from Diversity Outbred (DO) and Collaborative Cross (CC) mouse studies.
#' Provides tools for loading NPZ genotype probability files, processing
#' recombination data, and creating publication-ready genome plots.
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
#' \code{\link{smart_npz_loader}} for loading NPZ files
#' \code{\link{configure_plot}} for customizing plot appearance
NULL

#' Default configuration for genome plots
#'
#' @description
#' Default configuration parameters for genome plot generation including layout,
#' font sizes, colors, and spacing options.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{bar_height}{Height of each founder bar (default: 0.6)}
#'   \item{founder_gap}{Gap between founder1 and founder2 bars (default: 0.3)}
#'   \item{chrom_spacing}{Vertical spacing between chromosomes (default: 1.2)}
#'   \item{font_size_title}{Title font size (default: 24)}
#'   \item{font_size_axis}{Axis label font size (default: 20)}
#'   \item{font_size_tick}{Tick label font size (default: 18)}
#'   \item{font_size_legend}{Legend font size (default: 16)}
#'   \item{font_size_recomb}{Recombination count font size (default: 6)}
#'   \item{font_family}{Font family for all text (default: 'sans')}
#'   \item{background_color}{Background color (default: 'white')}
#'   \item{panel_border}{Show panel borders (default: TRUE)}
#'   \item{grid_alpha}{Grid line transparency (default: 0.3)}
#'   \item{tile_width_factor}{Width factor for tiles (default: 1.0)}
#'   \item{legend_position}{Legend position (default: 'right')}
#'   \item{show_grid}{Show grid lines (default: TRUE)}
#'   \item{y_scale_factor}{Y-axis scaling factor (default: 1.1)}
#'   \item{compact_mode}{Enable compact mode (default: FALSE)}
#'   \item{max_cM}{Maximum cM for x-axis (default: 90)}
#'   \item{grid_width}{Grid width for cM calculation (default: 0.01)}
#'   \item{plot_margin}{Plot margins (default: c(1.5, 1.5, 1.5, 1.5))}
#'   \item{title_margin}{Title spacing (default: c(0, 0, 0.5, 0))}
#'   \item{axis_line_color}{Axis line color (default: 'black')}
#'   \item{axis_line_width}{Axis line thickness (default: 0.8)}
#' }
#'
#' @export
DEFAULT_CONFIG <- list(
    # layout parameters
    bar_height = 0.6,           # height of each founder bar
    founder_gap = 0.3,          # gap between founder1 and founder2 bars
    chrom_spacing = 1.2,        # vertical spacing between chromosomes

    # font sizes
    font_size_title = 24,       # title font size
    font_size_axis = 20,        # axis label font size
    font_size_tick = 18,        # tick label font size
    font_size_legend = 16,      # legend font size
    font_size_recomb = 6,       # recombination count font size

    # font family
    font_family = 'sans',       # font family for all text

    # background and color settings
    background_color = 'white', # clean white background for publication
    panel_border = TRUE,        # show panel borders for definition
    grid_alpha = 0.3,           # grid line transparency (subtle but visible)

    # plot dimensions and scaling
    tile_width_factor = 1.0,    # width factor for tiles (1.0 = optimal)
    legend_position = 'right',  # legend position: 'right', 'bottom', 'none'
    show_grid = TRUE,           # show subtle grid lines for reference
    y_scale_factor = 1.1,       # slightly increased y-axis scaling for breathing room
    compact_mode = FALSE,       # enable compact mode for tight spacing

    # data parameters
    max_cM = 90,                # maximum cM for x-axis
    grid_width = 0.01,          # grid width for cM calculation

    # Color and theme enhancements
    plot_margin = c(1.5, 1.5, 1.5, 1.5), # plot margins
    title_margin = c(0, 0, 0.5, 0),      # title spacing
    axis_line_color = 'black',           # clear axis lines
    axis_line_width = 0.8                # appropriate axis line thickness
)

#' Default founder strains for DO/CC mice
#'
#' @description
#' Standard founder strain names used in Diversity Outbred (DO) and 
#' Collaborative Cross (CC) mouse studies.
#'
#' @format A character vector with 8 founder strains:
#' \describe{
#'   \item{A}{A/J}
#'   \item{B}{C57BL/6J}
#'   \item{C}{129S1/SvImJ}
#'   \item{D}{NOD/ShiLtJ}
#'   \item{E}{NZO/H1LtJ}
#'   \item{F}{CAST/EiJ}
#'   \item{G}{PWK/PhJ}
#'   \item{H}{WSB/EiJ}
#' }
#'
#' @export
DEFAULT_FOUNDERS <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')

#' Default founder colors for DO/CC mice visualization
#'
#' @description
#' Standard color scheme for founder strains in DO/CC mouse genome plots.
#' These colors provide good contrast and are commonly used in the field.
#'
#' @format A named character vector with hex color codes:
#' \describe{
#'   \item{A}{#F0F000 (Yellow)}
#'   \item{B}{#808080 (Gray)}
#'   \item{C}{#F08080 (Light Red)}
#'   \item{D}{#1010F0 (Blue)}
#'   \item{E}{#00A0F0 (Light Blue)}
#'   \item{F}{#00A000 (Green)}
#'   \item{G}{#F00000 (Red)}
#'   \item{H}{#9000E0 (Purple)}
#' }
#'
#' @export
DEFAULT_FOUNDER_COLORS <- c(
    A = '#F0F000', # A/J
    B = '#808080', # C57BL/6J
    C = '#F08080', # 129S1/SvImJ
    D = '#1010F0', # NOD/ShiLtJ
    E = '#00A0F0', # NZO/H1LtJ
    F = '#00A000', # CAST/EiJ
    G = '#F00000', # PWK/PhJ
    H = '#9000E0'  # WSB/EiJ
)

#' Alternative founder colors for DO/CC mice visualization
#'
#' @description
#' Alternative color scheme for founder strains with improved accessibility
#' and publication quality. These colors are optimized for colorblind-friendly
#' visualization.
#'
#' @format A named character vector with hex color codes:
#' \describe{
#'   \item{A}{#F0E442 (Yellow)}
#'   \item{B}{#555555 (Dark Gray)}
#'   \item{C}{#E69F00 (Orange)}
#'   \item{D}{#0072B2 (Blue)}
#'   \item{E}{#56B4E9 (Light Blue)}
#'   \item{F}{#009E73 (Green)}
#'   \item{G}{#D55E00 (Red-Orange)}
#'   \item{H}{#CC79A7 (Pink)}
#' }
#'
#' @export
DEFAULT_FOUNDER_COLORS_NEW <- c(
    A = '#F0E442', # A/J
    B = '#555555', # C57BL/6J
    C = '#E69F00', # 129S1/SvImJ
    D = '#0072B2', # NOD/ShiLtJ
    E = '#56B4E9', # NZO/H1LtJ
    F = '#009E73', # CAST/EiJ
    G = '#D55E00', # PWK/PhJ
    H = '#CC79A7'  # WSB/EiJ
)


#' Mouse chromosome lengths in base pairs
#'
#' @description
#' Reference chromosome lengths for Mus musculus (house mouse) in base pairs.
#' These values are used for proper scaling and positioning in genome plots.
#'
#' @format A named numeric vector with chromosome lengths:
#' \describe{
#'   \item{1-19}{Autosomal chromosomes with lengths ranging from ~61 to ~195 Mb}
#'   \item{X}{Sex chromosome with length ~169 Mb}
#' }
#'
#' @source Mouse genome assembly data
#' @export
DEFAULT_CHROM_LENGTHS <- c(
    '1' = 195154279, '2' = 181755017, '3' = 159745316, '4' = 156860686, '5' = 151758149,
    '6' = 149588044, '7' = 144995196, '8' = 130127694, '9' = 124359700, '10' = 130530862,
    '11' = 121973369, '12' = 120092757, '13' = 120883175, '14' = 125139656, '15' = 104073951,
    '16' = 98008968, '17' = 95294699, '18' = 90720763, '19' = 61420004, 'X' = 169476592
)

#' Generate diplotypes from founder strains
#'
#' @description
#' Creates all possible diplotype combinations from the founder strains.
#' This is essential for interpreting the genotype probability matrices.
#'
#' @param founders Character vector of founder strain names (default: DEFAULT_FOUNDERS)
#'
#' @return Character vector of diplotypes (e.g., "AA", "AB", "AC", etc.)
#'
#' @examples
#' # Generate diplotypes for standard DO/CC founders
#' diplotypes <- generate_diplotypes()
#' head(diplotypes)
#'
#' # Generate diplotypes for custom founders
#' custom_diplotypes <- generate_diplotypes(c("X", "Y", "Z"))
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


#' Smart NPZ file loader with fallback support
#' 
#' Attempts to load NPZ files using reticulate (Python/numpy).
#' 
#' @param npz_file Path to NPZ file
#' @return List containing loaded data and method used
#' Smart NPZ file loader with fallback options
#'
#' @description
#' Attempts to load NPZ files using different methods with intelligent
#' fallback strategies. First tries reticulate (Python/numpy), then
#' alternative methods if available.
#'
#' @param npz_file Path to the NPZ file to load
#'
#' @return A list containing the loaded data with chromosome names as keys
#'
#' @details
#' This function provides a robust way to load NPZ files containing
#' genotypic probability matrices. It handles various file formats and
#' provides informative error messages if loading fails.
#'
#' @examples
#' \dontrun{
#' # Load an NPZ file
#' data <- smart_npz_loader("path/to/file.npz")
#' }
#'
#' @export
smart_npz_loader <- function(npz_file) {
    if (!file.exists(npz_file)) {
        stop('NPZ file does not exist: ', npz_file)
    }
    
    # try reticulate (Python/numpy)
    if (requireNamespace('reticulate', quietly = TRUE)) {
        tryCatch({
            message('Attempting to load with reticulate (Python/numpy)...')
            np <- reticulate::import('numpy')
            npz <- np$load(npz_file, allow_pickle = TRUE)
            chrom_names <- as.character(npz$files)

            # extract data for each chromosome
            data_list <- list()
            for (chrom in chrom_names) {
                data_list[[chrom]] <- npz$f[[chrom]]
            }

            message('Successfully loaded with reticulate')
            return(list(
                data = data_list,
                chrom_names = chrom_names,
                method = 'reticulate',
                success = TRUE
            ))
        }, error = function(e) {
            message('Reticulate failed: ', e$message)
        })
    }
    
    stop('Please install reticulate')
}

#' Load and process NPZ file
#'
#' @description
#' Loads genotype probability data from NPZ files and prepares it for analysis.
#' This function handles the data loading step with informative messaging and
#' converts the data to a standardized R data frame format.
#'
#' @param npz_file Path to NPZ file
#' @param founders Character vector of founder strain names (default: DEFAULT_FOUNDERS)
#' @param chrom_lengths Named vector of chromosome lengths (default: DEFAULT_CHROM_LENGTHS)
#'
#' @return Processed data frame with the following columns:
#' \itemize{
#'   \item{chromosome}{Chromosome identifier}
#'   \item{marker_index}{Marker position index (0-based)}
#'   \item{diplotype columns}{Probability columns for each diplotype combination}
#' }
#'
#' @details
#' This function is the main entry point for loading NPZ files. It:
#' \itemize{
#'   \item{Generates all possible diplotype combinations from the founder strains}
#'   \item{Loads the NPZ file using smart_npz_loader}
#'   \item{Processes each chromosome's probability matrix}
#'   \item{Combines all chromosomes into a single data frame}
#' }
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
                          founders = DEFAULT_FOUNDERS, 
                          chrom_lengths = DEFAULT_CHROM_LENGTHS) {
    
    # generate diplotypes
    diplotypes <- generate_diplotypes(founders)
    
    # load data
    message('Loading NPZ file: ', basename(npz_file))
    loaded_data <- smart_npz_loader(npz_file)
    
    if (!loaded_data$success) {
        stop('Failed to load NPZ file')
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
        df_chr <- as.data.frame(t(mat))  # transpose to match expected format
        colnames(df_chr) <- diplotypes
        df_chr$chromosome <- chrom
        df_chr$marker_index <- 0:(n_markers-1)
        
        all_df[[chrom]] <- df_chr
    }
    
    # combine all chromosomes
    df <- dplyr::bind_rows(all_df)
    
    message('Total markers loaded: ', nrow(df))
    return(df)
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
#'
#' @return Processed data frame with founder assignments including:
#' \itemize{
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
#' }
#'
#' @export
process_genotype_data <- function(df, founders = DEFAULT_FOUNDERS, config = DEFAULT_CONFIG) {

    # identify diplotype columns
    diplotype_cols <- grep('^[A-H][A-H]$', names(df), value = TRUE)

    if (length(diplotype_cols) == 0) {
        stop('No diplotype columns found in data frame')
    }

    message('Found ', length(diplotype_cols), ' diplotype columns')

    # for each marker, find the diplotype with the highest probability
    df <- df %>%
        mutate(
            max_diplotype = diplotype_cols[apply(select(., all_of(diplotype_cols)), 1, which.max)],
            founder2 = substr(max_diplotype, 1, 1), # first founder in diplotype
            founder1 = substr(max_diplotype, 2, 2)  # second founder in diplotype
        )

    # order chromosomes for plotting (natural sort)
    chrom_order <- names(DEFAULT_CHROM_LENGTHS)[names(DEFAULT_CHROM_LENGTHS) %in% unique(df$chromosome)]
    chrom_order <- chrom_order[order(as.numeric(gsub('[^0-9]', '', chrom_order)))]
    df$chromosome <- factor(df$chromosome, levels = chrom_order)

    message('Chromosomes ordered: ', paste(chrom_order, collapse = ', '))

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

    chr_ypos <- setNames(seq(start_pos, spacing_step, by = -spacing_step), chrom_order)

    # ensure we have exactly the right number of positions
    if (length(chr_ypos) != length(chrom_order)) {
        # Fallback: create positions manually
        message('Creating fallback chromosome positions')
        chr_ypos <- setNames(seq(num_chrs, 1), chrom_order)
    }

    # Final validation
    if (length(chr_ypos) != length(chrom_order)) {
        stop('Critical error: Cannot create valid chromosome positions. chr_ypos length: ', length(chr_ypos), ', chrom_order length: ', length(chrom_order))
    }

    df <- df %>% dplyr::mutate(y_base = chr_ypos[as.character(chromosome)])

    # calculate cM positions
    df <- df %>%
        dplyr::group_by(chromosome) %>%
        dplyr::mutate(cM = marker_index * config$grid_width) %>%
        dplyr::ungroup()

    message('Genotypic data processing completed')
    return(df)
}

#' Count recombination events per chromosome
#'
#' @description
#' Identifies and counts recombination events by detecting changes in
#' founder strain assignments between adjacent genomic positions. This
#' function is essential for quantifying genetic diversity and mapping
#' recombination hotspots.
#'
#' @param df Processed data frame from process_genotype_data()
#'
#' @return Data frame with recombination counts containing:
#' \itemize{
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
#' @export
count_recombinations <- function(df) {

    recomb_counts <- df %>%
        dplyr::group_by(chromosome) %>%
        dplyr::summarize(
            num_recomb = sum(dplyr::lag(max_diplotype, default = dplyr::first(max_diplotype)) != max_diplotype),
            .groups = 'drop'
        )

    message('Recombination counts calculated for ', nrow(recomb_counts), ' chromosomes')
    return(recomb_counts)
}

#' Prepare data for plotting
#'
#' @description
#' Transforms the processed data into long format for ggplot2 visualization
#' with proper positioning. This function converts the wide-format data
#' into a long format suitable for creating stacked bar plots.
#'
#' @param df Processed data frame from process_genotype_data()
#' @param config Configuration list (default: DEFAULT_CONFIG)
#'
#' @return Long-format data frame for plotting containing:
#' \itemize{
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
#' @export
prepare_plot_data <- function(df, config = DEFAULT_CONFIG) {

    # Create long format for plotting
    df_long <- df %>%
        tidyr::pivot_longer(cols = c('founder1', 'founder2'), names_to = 'hap', values_to = 'founder') %>%
        dplyr::mutate(
            y = ifelse(hap == 'founder1', y_base, y_base + config$bar_height + config$founder_gap),
            yend = y + config$bar_height,
            cM = marker_index * (config$max_cM / max(df$marker_index))
        )

    message('Data transformed to long format for plotting')
    return(df_long)
}

#' Create enhanced genome reconstruction plot
#'
#' @description
#' Generates a publication-quality genome plot with improved aesthetics,
#' better spacing, and professional visual design. This function creates
#' the core visualization showing founder strain contributions across
#' the genome.
#'
#' @param df_long Long-format data frame from prepare_plot_data()
#' @param recomb_counts Recombination counts data frame from count_recombinations()
#' @param config Configuration list (default: DEFAULT_CONFIG)
#' @param founder_colors Named vector of founder colors (default: DEFAULT_FOUNDER_COLORS)
#' @param sample_name Sample name for plot title (default: 'Sample')
#'
#' @return ggplot object representing the genome reconstruction plot
#'
#' @details
#' This function creates a comprehensive genome visualization including:
#' \itemize{
#'   \item{Stacked bar plots showing founder contributions}
#'   \item{Chromosome organization with proper spacing}
#'   \item{Recombination event annotations}
#'   \item{Professional styling and color schemes}
#'   \item{Configurable layout and appearance}
#' }
#'
#' @examples
#' \dontrun{
#' # Create a basic genome plot
#' plot <- create_genome_plot(plot_data, recomb_data)
#' 
#' # Create plot with custom colors and configuration
#' config <- configure_plot(compact_mode = TRUE)
#' plot <- create_genome_plot(plot_data, recomb_data, config, 
#'                           founder_colors = DEFAULT_FOUNDER_COLORS_NEW)
#' }
#'
#' @export
create_genome_plot <- function(df_long,
                               recomb_counts,
                               config = DEFAULT_CONFIG,
                               founder_colors = DEFAULT_FOUNDER_COLORS,
                               sample_name = 'Sample') {
    # get chromosome positions
    chr_ypos <- df_long %>%
        dplyr::group_by(chromosome) %>%
        dplyr::summarize(y_base = dplyr::first(y_base), .groups = 'drop') %>%
        dplyr::pull(y_base, name = chromosome)

    # create plot
    p <- ggplot(df_long, aes(x = cM, y = y, fill = founder)) +
        geom_tile(
            aes(
                height = config$bar_height,
                width = (config$max_cM / max(marker_index)) * config$tile_width_factor
            ),
            show.legend = FALSE,
            color = NA, # remove tile borders for cleaner look
            alpha = 0.9 # slight transparency for depth
        ) + 
        # color scheme
        scale_fill_manual(values = founder_colors) +
        # improved y-axis with better spacing
        scale_y_continuous(
            breaks = chr_ypos + config$bar_height + config$founder_gap / 2,
            labels = names(chr_ypos),
            expand = c(0, 0),
            limits = c(0, max(chr_ypos) + config$bar_height + config$founder_gap + config$bar_height + 0.5)
        ) +
        # coordinate system
        coord_cartesian(ylim = c(0, max(chr_ypos) + config$bar_height + config$founder_gap + config$bar_height + 0.5)) +
        # labels
        labs(
            title = paste0('Genome Reconstruction: ', sample_name, '\n(Total ', sum(recomb_counts$num_recomb), ' Recombinations)'),
            x = 'Genetic Distance (cM)',
            y = 'Chromosome',
            fill = 'Founder Strain'
        ) +
        # theme
        theme_minimal(base_size = config$font_size_tick, base_family = config$font_family) +
        theme(
            # panel styling
            panel.background = element_rect(fill = config$background_color, color = NA),
            plot.background = element_rect(fill = config$background_color, color = NA),

            # grid system
            panel.grid.major = if (config$show_grid) element_line(color = 'gray90') else element_blank(),
            panel.grid.minor = element_blank(),

            # panel borders
            panel.border = if (config$panel_border) element_rect(color = config$axis_line_color, fill = NA, linewidth = config$axis_line_width) else element_blank(),

            # axis styling
            axis.line = element_line(color = config$axis_line_color, linewidth = config$axis_line_width),
            axis.text.y = element_text(size = config$font_size_tick, family = config$font_family, color = 'black'),
            axis.text.x = element_text(size = config$font_size_tick, family = config$font_family, color = 'black'),
            axis.title = element_text(size = config$font_size_axis, family = config$font_family, color = 'black', face = 'bold'),

            # title styling
            plot.title = element_text(
                size = config$font_size_title,
                face = 'bold',
                family = config$font_family,
                margin = margin(config$title_margin[1], config$title_margin[2], config$title_margin[3], config$title_margin[4], 'cm'),
                hjust = 0.5,
                color = 'black'
            ),

            # legend styling
            legend.position = if (config$legend_position == 'none') 'none' else config$legend_position,
            legend.text = element_text(size = config$font_size_legend, family = config$font_family),
            legend.title = element_text(size = config$font_size_legend, family = config$font_family, face = 'bold'),

            # margins and spacing
            plot.margin = margin(config$plot_margin[1], config$plot_margin[2], config$plot_margin[3], config$plot_margin[4], 'cm'),

            # remove unnecessary elements for cleaner look
            legend.background = element_blank(),
            legend.box.background = element_blank()
        )

    # add recombination counts with better positioning
    for (i in seq_along(names(chr_ypos))) {
        chr <- names(chr_ypos)[i]
        # center the text in the gap between the two bars
        ypos <- chr_ypos[chr] + config$bar_height + config$founder_gap / 2
        nrec <- recomb_counts$num_recomb[recomb_counts$chromosome == chr]

        # find the max cM for this chromosome
        chr_data <- df_long %>% filter(chromosome == chr, hap == 'founder2')
        xmax <- max(chr_data$cM, na.rm = TRUE)

        p <- p + annotate(
            'text',
            x = xmax + 1,
            y = ypos,
            label = paste0('(', nrec, ')'),
            hjust = 0,
            size = config$font_size_recomb,
            family = config$font_family,
            color = 'darkgray', # subtle color for recombination counts
            fontface = 'italic' # italic for secondary information
        )
    }

    message('Genome plot created successfully')
    return(p)
}

#' Main function to generate enhanced genome reconstruction plot from NPZ file
#'
#' This is the primary function that orchestrates the entire genome visualization
#' pipeline with enhanced aesthetics and informative progression messages.
#'
#' @param npz_file Path to NPZ file
#' @param founders Character vector of founder strain names
#' @param founder_colors Named vector of founder colors
#' @param chrom_lengths Named vector of chromosome lengths
#' Main function to create GBRS genome plots
#'
#' @description
#' This is the main function for creating genome reconstruction plots from NPZ files.
#' It orchestrates the entire workflow from data loading to final plot generation,
#' providing a simple interface for users while maintaining flexibility.
#'
#' @param npz_file Path to NPZ file containing genotypic probability data
#' @param founders Character vector of founder strain names (default: DEFAULT_FOUNDERS)
#' @param founder_colors Named vector of founder colors (default: DEFAULT_FOUNDER_COLORS)
#' @param chrom_lengths Named vector of chromosome lengths (default: DEFAULT_CHROM_LENGTHS)
#' @param config Configuration list (default: DEFAULT_CONFIG)
#' @param sample_name Sample name for plot title (default: extracted from filename)
#' @param return_data Logical, whether to return processed data along with plot (default: FALSE)
#'
#' @return Either a ggplot object (default) or a list containing:
#' \itemize{
#'   \item{plot}{The generated genome plot}
#'   \item{data}{Processed genotype data}
#'   \item{data_long}{Long-format data for plotting}
#'   \item{recombinations}{Recombination count data}
#'   \item{config}{Configuration used}
#' }
#'
#' @details
#' This function performs a complete workflow:
#' \itemize{
#'   \item{Loads and validates NPZ file data}
#'   \item{Processes genotype probabilities into founder assignments}
#'   \item{Counts recombination events per chromosome}
#'   \item{Prepares data for visualization}
#'   \item{Generates publication-quality plots}
#'   \item{Provides comprehensive progress reporting}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot <- plot_gbrs_genome("example.npz")
#' 
#' # With custom configuration
#' config <- configure_plot(compact_mode = TRUE, bar_height = 0.3)
#' plot <- plot_gbrs_genome("example.npz", config = config)
#' 
#' # Return both plot and data
#' result <- plot_gbrs_genome("example.npz", return_data = TRUE)
#' plot <- result$plot
#' data <- result$data
#' }
#'
#' @export
plot_gbrs_genome <- function(npz_file,
                             founders = DEFAULT_FOUNDERS,
                             founder_colors = DEFAULT_FOUNDER_COLORS,
                             chrom_lengths = DEFAULT_CHROM_LENGTHS,
                             config = DEFAULT_CONFIG,
                             sample_name = NULL,
                             return_data = FALSE) {

    message('Starting genome plot...')

    # Input validation
    if (missing(npz_file) || is.null(npz_file)) {
        stop('NPZ file path is required')
    }

    if (is.null(sample_name)) {
        sample_name <- tools::file_path_sans_ext(basename(npz_file))
    }

    message('Input parameters:')
    message('  - NPZ file: ', npz_file)
    message('  - Sample: ', sample_name)
    message('  - Founders: ', paste(founders, collapse = ', '))
    message('  - Configuration: ', length(config), ' parameters')

    # Step 1: Load and process data
    message('\nStep 1/5: Loading NPZ file...')
    df <- load_npz_file(npz_file, founders, chrom_lengths)

    # Step 2: Process genotype data
    message('\nStep 2/5: Processing genotypic data...')
    df <- process_genotype_data(df, founders, config)

    # Step 3: Count recombinations
    message('\nStep 3/5: Counting recombination events...')
    recomb_counts <- count_recombinations(df)

    # Step 4: Prepare plot data
    message('\nStep 4/5: Preparing plot data...')
    df_long <- prepare_plot_data(df, config)

    # Step 5: Create enhanced plot
    message('\nStep 5/5: Creating enhanced genome plot...')
    p <- create_genome_plot(df_long, recomb_counts, config, founder_colors, sample_name)

    # recombination summary
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

    message('\nDone!')

    # Return results
    if (return_data) {
        return(list(
            plot = p,
            data = df,
            data_long = df_long,
            recombinations = recomb_counts,
            config = config
        ))
    } else {
        return(p)
    }
}

#' Configuration function with validation and documentation
#'
#' @description
#' Provides an improved interface for customizing plot parameters with
#' comprehensive validation and helpful documentation. This function allows
#' users to easily modify plot settings while ensuring parameter validity.
#'
#' @param ... Named parameters to update in configuration (e.g., bar_height = 0.5)
#'
#' @return Updated configuration list with validated parameters
#'
#' @details
#' This function provides several benefits:
#' \itemize{
#'   \item{Parameter validation to prevent invalid settings}
#'   \item{Clear feedback on what parameters were updated}
#'   \item{Warnings for invalid parameter names}
#'   \item{Maintains all default values for unspecified parameters}
#' }
#'
#' @examples
#' \dontrun{
#' # Modify a single parameter
#' config <- configure_plot(bar_height = 0.3)
#' 
#' # Modify multiple parameters
#' config <- configure_plot(
#'     bar_height = 0.3,
#'     compact_mode = TRUE,
#'     font_size_title = 18
#' )
#' 
#' # View current configuration
#' print(config)
#' }
#'
#' @export
configure_plot <- function(...) {
    config <- DEFAULT_CONFIG
    updates <- list(...)

    # validate parameter names
    valid_params <- names(DEFAULT_CONFIG)
    invalid_params <- setdiff(names(updates), valid_params)
    if (length(invalid_params) > 0) {
        warning('Invalid parameters ignored: ', paste(invalid_params, collapse = ', '))
        message('Valid parameters are: ', paste(valid_params, collapse = ', '))
    }

    # update valid parameters
    for (param in intersect(names(updates), valid_params)) {
        old_value <- config[[param]]
        config[[param]] <- updates[[param]]
        message('Updated ', param, ': ', old_value, ' → ', updates[[param]])
    }

    return(config)
}

#' Enhanced plot saving with optimal dimensions and quality
#'
#' Saves plots with improved default settings and automatic dimension
#' calculation for optimal visual presentation.
#'
#' @param plot ggplot object
#' @param filename Output filename
#' Save genome plot to file
#'
#' @description
#' Saves genome plots to various file formats with intelligent sizing
#' and comprehensive output options. This function provides a convenient
#' way to export plots for publication or further analysis.
#'
#' @param plot ggplot object to save
#' @param filename Output filename (with or without extension)
#' @param width Figure width in inches (default: 20)
#' @param height Figure height in inches (default: auto-calculated)
#' @param dpi Resolution for raster formats (default: 300)
#' @param format Output format: 'pdf', 'png', 'svg', or 'tiff' (default: 'pdf')
#'
#' @return Invisibly returns the filename of the saved plot
#'
#' @details
#' This function provides several features:
#' \itemize{
#'   \item{Auto-calculates height based on number of chromosomes}
#'   \item{Automatically adds file extensions if missing}
#'   \item{Creates output directories if they don't exist}
#'   \item{Supports multiple output formats}
#'   \item{Provides detailed feedback on the save operation}
#' }
#'
#' @examples
#' \dontrun{
#' # Save as PDF (default)
#' save_genome_plot(plot, "genome_plot.pdf")
#' 
#' # Save as PNG with custom dimensions
#' save_genome_plot(plot, "genome_plot.png", width = 16, height = 12)
#' 
#' # Save as SVG for vector graphics
#' save_genome_plot(plot, "genome_plot", format = "svg")
#' }
#'
#' @export
save_genome_plot <- function(plot,
                             filename,
                             width = 20, 
                             height = NULL,
                             dpi = 300,
                             format = 'pdf') {
    message('Saving genome plot...')

    # auto-calculate height if not specified
    if (is.null(height)) {
        # try to extract number of chromosomes from plot data
        plot_data <- plot$data
        if (!is.null(plot_data) && 'chromosome' %in% names(plot_data)) {
            num_chroms <- length(unique(plot_data$chromosome))
            height <- max(10, num_chroms * 1.8) # Increased height per chromosome
        } else {
            # fallback: use default height
            height <- 14
        }
    }

    # ensure height is numeric
    if (!is.numeric(height)) {
        warning('Height is not numeric, using default value')
        height <- 14
    }

    # ensure file extension matches format
    if (!grepl(paste0('\\.', format, '$'), filename)) {
        filename <- paste0(filename, '.', format)
    }

    # create output directory if it doesn't exist
    output_dir <- dirname(filename)
    if (output_dir != '.' && !dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # save plot
    ggsave(
        filename, 
        plot,
        width = width, 
        height = height, 
        dpi = dpi,
        limitsize = FALSE # allow larger dimensions
    ) 

    message('Plot saved successfully:')
    message('  - File: ', filename)
    message('  - Dimensions: ', width, ' × ', height, ' inches')
    message('  - Format: ', format)
    message('  - Resolution: ', dpi, ' DPI')
}

#' List available fonts on the system
#'
#' @description
#' Lists available fonts on the system using various methods. This function
#' tries multiple approaches to find available fonts, providing fallbacks
#' for different system configurations.
#'
#' @return Character vector of available font names
#'
#' @details
#' This function attempts to list fonts using:
#' \itemize{
#'   \item{extrafont package (if available)}
#'   \item{systemfonts package (if available)}
#'   \item{Fallback to common system fonts}
#' }
#'
#' @examples
#' \dontrun{
#' # List available fonts
#' fonts <- list_available_fonts()
#' head(fonts)
#' 
#' # Check if specific fonts are available
#' "Arial" %in% list_available_fonts()
#' }
#'
#' @export
list_available_fonts <- function() {
    # Check if extrafont is available
    if (requireNamespace('extrafont', quietly = TRUE)) {
        message('Using extrafont to list fonts...')
        extrafont::font_import(prompt = FALSE)
        fonts <- extrafont::fonts()
        return(fonts)
    }

    # Check if systemfonts is available
    if (requireNamespace('systemfonts', quietly = TRUE)) {
        message('Using systemfonts to list fonts...')
        fonts <- systemfonts::system_fonts()$family
        return(unique(fonts))
    }

    # Fallback: common fonts that are usually available
    message('Using fallback method - common fonts:')
    common_fonts <- c(
        'sans', 'serif', 'mono',
        'Arial', 'Times', 'Courier',
        'Helvetica', 'Georgia', 'Verdana',
        'Palatino', 'Bookman', 'Avant Garde'
    )
    return(common_fonts)
}

#' Create compact configuration for publication
#'
#' @description
#' Provides a pre-configured compact layout suitable for smaller publication sections
#' with tight spacing and smaller elements. This function creates an optimized
#' configuration for space-constrained visualizations.
#'
#' @return Compact configuration list with the following modifications:
#' \itemize{
#'   \item{bar_height: 0.15 (smaller bars)}
#'   \item{founder_gap: 0.05 (minimal gap)}
#'   \item{chrom_spacing: 0.25 (tight spacing)}
#'   \item{compact_mode: TRUE (enabled)}
#'   \item{Reduced font sizes for all text elements}
#'   \item{Simplified theme (no grid, no borders)}
#'   \item{Reduced margins and spacing}
#' }
#'
#' @details
#' This configuration is ideal for:
#' \itemize{
#'   \item{Publication figures with space constraints}
#'   \item{Multi-panel layouts}
#'   \item{Poster presentations}
#'   \item{Compact reports}
#' }
#'
#' @examples
#' \dontrun{
#' # Use compact configuration
#' config <- create_compact_config()
#' plot <- plot_gbrs_genome("example.npz", config = config)
#' 
#' # Modify compact configuration further
#' config <- create_compact_config()
#' config <- configure_plot(config, font_size_title = 14)
#' }
#'
#' @export
create_compact_config <- function() {
    config <- DEFAULT_CONFIG

    # compact layout settings
    config$bar_height <- 0.15       # smaller bars
    config$founder_gap <- 0.05      # minimal gap between haplotypes
    config$chrom_spacing <- 0.25    # tight chromosome spacing
    config$compact_mode <- TRUE     # enable compact mode

    # smaller fonts for compact layout
    config$font_size_title <- 18    # smaller title
    config$font_size_axis <- 16     # smaller axis labels
    config$font_size_tick <- 14     # smaller tick labels
    config$font_size_legend <- 12   # smaller legend
    config$font_size_recomb <- 5    # smaller recombination counts

    # reduced margins and spacing
    config$plot_margin <- c(1.0, 1.0, 1.0, 1.0)
    config$title_margin <- c(0, 0, 0.3, 0)

    # simplified theme for compact display
    config$show_grid <- FALSE       # no grid lines
    config$panel_border <- FALSE    # no panel borders

    message('Compact configuration created for publication use')
    return(config)
}


