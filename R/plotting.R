#' Create enhanced genome reconstruction plot
#'
#' @description
#' Generates a publication-quality genome plot with improved aesthetics,
#' better spacing, and professional visual design. This function creates
#' the core visualization showing founder strain contributions across
#' the genome.
#'
#' @param df Processed data frame from process_genotype_data()
#' @param config Configuration list (default: DEFAULT_CONFIG)
#' @param founder_colors Named vector of founder colors (default: DEFAULT_FOUNDER_COLORS)
#' @param sample_name Sample name for plot title (default: 'Sample')
#'
#' @return ggplot object representing the genome reconstruction plot
#'
#' @details
#' This function creates a comprehensive genome visualization including:
#' \itemize{
#'   \item{Automatic internal processing of genotype data}
#'   \item{Stacked bar plots showing founder contributions}
#'   \item{Chromosome organization with proper spacing}
#'   \item{Recombination event annotations}
#'   \item{Professional styling and color schemes}
#'   \item{Configurable layout and appearance}
#' }
#'
#' @section Internal Processing:
#' This function automatically handles:
#' \itemize{
#'   \item{Recombination counting via \code{count_recombinations()}}
#'   \item{Data preparation for plotting via \code{prepare_plot_data()}}
#'   \item{All visualization generation}
#' }
#'
#' @examples
#' \dontrun{
#' # Create a basic genome plot from processed data
#' plot <- create_genome_plot(processed_data)
#'
#' # Create plot with custom colors and configuration
#' config <- configure_plot(compact_mode = TRUE)
#' plot <- create_genome_plot(processed_data, config,
#'     founder_colors = DEFAULT_FOUNDER_COLORS_NEW
#' )
#' }
#'
#' @export
create_genome_plot <- function(df,
                               config = DEFAULT_CONFIG,
                               founder_colors = DEFAULT_FOUNDER_COLORS,
                               sample_name = 'Sample') {
    # Internal processing - count recombinations and prepare plot data
    recomb_counts <- count_recombinations(df)
    df_long <- prepare_plot_data(df, config)

    message('Internal data processing completed for plotting')
    # get chromosome positions
    chr_ypos <- df_long %>%
        dplyr::group_by(.data$chromosome) %>%
        dplyr::summarize(y_base = dplyr::first(.data$y_base), .groups = 'drop') %>%
        dplyr::pull(.data$y_base, name = .data$chromosome)

    # create plot
    p <- ggplot(df_long, aes(x = .data$cM, y = .data$y, fill = .data$founder)) +
        geom_tile(
            aes(
                height = config$bar_height,
                width = (config$max_cM / max(.data$marker_index)) * config$tile_width_factor
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
        chr_data <- df_long %>% dplyr::filter(.data$chromosome == chr, .data$hap == 'founder2')
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

#' Main function to create GBRS genome plots
#'
#' @description
#' This is the main function for creating genome reconstruction plots from NPZ files.
#' It orchestrates the entire workflow from data loading to final plot generation,
#' providing a simple interface for users while maintaining flexibility.
#'
#' @section CRITICAL DEPENDENCY WARNING:
#' This function REQUIRES the reticulate R package and Python with NumPy installed.
#' Without these dependencies, the function will fail with errors like:
#' - 'Error: Please install reticulate'
#' - 'FileNotFoundError: No such file or directory'
#'
#' See the README for complete installation instructions.
#'
#' @param npz_file Path to NPZ file containing genotypic probability data
#' @param founders Character vector of founder strain names (default: DEFAULT_FOUNDERS)
#' @param founder_colors Named vector of founder colors (default: DEFAULT_FOUNDER_COLORS)
#' @param chrom_lengths Named vector of chromosome lengths for ordering (default: DEFAULT_CHROM_LENGTHS)
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
#'   \item{Generates publication-quality plots with automatic internal processing}
#'   \item{Provides comprehensive progress reporting}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot <- plot_gbrs_genome('example.npz')
#'
#' # With custom configuration
#' config <- configure_plot(compact_mode = TRUE, bar_height = 0.3)
#' plot <- plot_gbrs_genome('example.npz', config = config)
#'
#' # With custom chromosome ordering
#' custom_chroms <- c('1' = 1000000, '2' = 2000000, 'X' = 1500000)
#' plot <- plot_gbrs_genome('example.npz', chrom_lengths = custom_chroms)
#'
#' # Return both plot and data
#' result <- plot_gbrs_genome('example.npz', return_data = TRUE)
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
    message('Starting genome plot generation...')

    # Input validation
    if (missing(npz_file) || is.null(npz_file)) {
        stop(
            'NPZ file path is required.\n',
            'Please provide the path to your NPZ file containing genotype data.'
        )
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
    message('\nStep 1/4: Loading NPZ file and processing data...')
    df <- load_npz_file(npz_file, founders)

    # Step 2: Process genotype data
    message('\nStep 2/4: Processing genotypic data and calculating positions...')
    df <- process_genotype_data(df, founders, config, chrom_lengths)

    # Step 3: Create enhanced plot
    message('\nStep 3/3: Creating enhanced genome plot with styling...')
    p <- create_genome_plot(df, config, founder_colors, sample_name)

    message('\nGenome plot generation completed successfully!')

    # Return results
    if (return_data) {
        # Get recombination counts and processed data for return
        recomb_counts <- count_recombinations(df)
        df_long <- prepare_plot_data(df, config)

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

#' Save genome plot to file
#'
#' @description
#' Saves genome plots to various file formats with automatic dimension
#' calculation and comprehensive error handling. This function provides
#' flexible output options for publication and presentation use.
#'
#' @param plot ggplot object to save
#' @param filename Output file path (with or without extension)
#' @param width Plot width in inches (default: 20)
#' @param height Plot height in inches (auto-calculated if NULL)
#' @param dpi Resolution in dots per inch (default: 300)
#' @param format Output format: 'pdf', 'png', 'svg', 'tiff', 'jpeg' (default: 'pdf')
#'
#' @return Invisibly returns the filename for chaining
#'
#' @details
#' This function provides several features:
#' \itemize{
#'   \item{Automatic height calculation based on chromosome count}
#'   \item{Multiple output format support}
#'   \item{Automatic file extension handling}
#'   \item{Directory creation if needed}
#'   \item{Comprehensive progress reporting}
#' }
#'
#' @examples
#' \dontrun{
#' # Save as PDF (default)
#' save_genome_plot(plot, 'genome_plot.pdf')
#'
#' # Save as PNG with custom dimensions
#' save_genome_plot(plot, 'genome_plot.png', width = 16, height = 12)
#'
#' # Save as SVG for vector graphics
#' save_genome_plot(plot, 'genome_plot', format = 'svg')
#' }
#'
#' @export
save_genome_plot <- function(plot,
                             filename,
                             width = 20,
                             height = NULL,
                             dpi = 300,
                             format = 'pdf') {
    message('Saving genome plot to file...')

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
    message('  - Dimensions: ', width, ' x ', height, ' inches')
    message('  - Format: ', format)
    message('  - Resolution: ', dpi, ' DPI')
}
