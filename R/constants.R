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
    bar_height = 0.6,                     # height of each founder bar
    founder_gap = 0.3,                    # gap between founder1 and founder2 bars
    chrom_spacing = 1.2,                  # vertical spacing between chromosomes

    # font sizes
    font_size_title = 24,                 # title font size
    font_size_axis = 20,                  # axis label font size
    font_size_tick = 18,                  # tick label font size
    font_size_legend = 16,                # legend font size
    font_size_recomb = 6,                 # recombination count font size

    # font family
    font_family = 'sans',                 # font family for all text

    # background and color settings
    background_color = 'white',           # clean white background for publication
    panel_border = TRUE,                  # show panel borders for definition
    grid_alpha = 0.3,                     # grid line transparency (subtle but visible)

    # plot dimensions and scaling
    tile_width_factor = 1.0,              # width factor for tiles (1.0 = optimal)
    legend_position = 'right',            # legend position: 'right', 'bottom', 'none'
    show_grid = TRUE,                     # show subtle grid lines for reference
    y_scale_factor = 1.1,                 # slightly increased y-axis scaling for breathing room
    compact_mode = FALSE,                 # enable compact mode for tight spacing

    # data parameters
    max_cM = 90,                          # maximum cM for x-axis
    grid_width = 0.01,                    # grid width for cM calculation

    # Color and theme enhancements
    plot_margin = c(1.5, 1.5, 1.5, 1.5),  # plot margins
    title_margin = c(0, 0, 0.5, 0),       # title spacing
    axis_line_color = 'black',            # clear axis lines
    axis_line_width = 0.8                 # appropriate axis line thickness
)


#' Compact configuration for genome plots
#'
#' @description
#' Compact configuration parameters for genome plot generation including layout,
#' font sizes, colors, and spacing options.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{bar_height}{Height of each founder bar (default: 0.15)}
#'   \item{founder_gap}{Gap between founder1 and founder2 bars (default: 0.05)}
#'   \item{chrom_spacing}{Vertical spacing between chromosomes (default: 0.25)}
#'   \item{font_size_title}{Title font size (default: 18)}
#'   \item{font_size_axis}{Axis label font size (default: 16)}
#'   \item{font_size_tick}{Tick label font size (default: 14)}
#'   \item{font_size_legend}{Legend font size (default: 12)}
#'   \item{font_size_recomb}{Recombination count font size (default: 5)}
#'   \item{font_family}{Font family for all text (default: 'sans')}
#'   \item{background_color}{Background color (default: 'white')}
#'   \item{panel_border}{Show panel borders (default: FALSE)}
#'   \item{grid_alpha}{Grid line transparency (default: 0.3)}
#'   \item{tile_width_factor}{Width factor for tiles (default: 1.0)}
#'   \item{legend_position}{Legend position (default: 'right')}
#'   \item{show_grid}{Show grid lines (default: FALSE)}
#'   \item{y_scale_factor}{Y-axis scaling factor (default: 1.1)}
#'   \item{compact_mode}{Enable compact mode (default: TRUE)}
#'   \item{max_cM}{Maximum cM for x-axis (default: 90)}
#'   \item{grid_width}{Grid width for cM calculation (default: 0.01)}
#'   \item{plot_margin}{Plot margins (default: c(1.0, 1.0, 1.0, 1.0))}
#'   \item{title_margin}{Title spacing (default: c(0, 0, 0.3, 0))}
#'   \item{axis_line_color}{Axis line color (default: 'black')}
#'   \item{axis_line_width}{Axis line thickness (default: 0.8)}
#' }
#'
#' @export
COMPACT_CONFIG <- list(
    # layout parameters
    bar_height = 0.15,                    # height of each founder bar
    founder_gap = 0.05,                   # gap between founder1 and founder2 bars
    chrom_spacing = 0.25,                 # vertical spacing between chromosomes

    # font sizes
    font_size_title = 18,                 # title font size
    font_size_axis = 16,                  # axis label font size
    font_size_tick = 14,                  # tick label font size
    font_size_legend = 12,                # legend font size
    font_size_recomb = 5,                 # recombination count font size

    # font family
    font_family = 'sans',                 # font family for all text

    # background and color settings
    background_color = 'white',           # clean white background for publication
    panel_border = FALSE,                 # show panel borders for definition
    grid_alpha = 0.3,                     # grid line transparency (subtle but visible)

    # plot dimensions and scaling
    tile_width_factor = 1.0,              # width factor for tiles (1.0 = optimal)
    legend_position = 'right',            # legend position: 'right', 'bottom', 'none'
    show_grid = FALSE,                    # show subtle grid lines for reference
    y_scale_factor = 1.1,                 # slightly increased y-axis scaling for breathing room
    compact_mode = TRUE,                  # enable compact mode for tight spacing

    # data parameters
    max_cM = 90,                          # maximum cM for x-axis
    grid_width = 0.01,                    # grid width for cM calculation

    # Color and theme enhancements
    plot_margin = c(1.0, 1.0, 1.0, 1.0),  # plot margins
    title_margin = c(0, 0, 0.3, 0),       # title spacing
    axis_line_color = 'black',            # clear axis lines
    axis_line_width = 0.8                 # appropriate axis line thickness
)


#' Default founder strains for DO/CC mice
#'
#' @description
#' Standard founder strain names used in Diversity Outbred (DO) and
#' Collaborative Cross (CC) mouse studies.
#'
#' @format A character vector with 8 founder strain identifiers:
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
#' Standard color scheme for founder strains in genome plots.
#' Colors are optimized for publication and accessibility.
#'
#' @format A named character vector with hex color codes:
#' \describe{
#'   \item{A}{#F0F000 - Yellow (A/J)}
#'   \item{B}{#808080 - Gray (C57BL/6J)}
#'   \item{C}{#F08080 - Light Red (129S1/SvImJ)}
#'   \item{D}{#1010F0 - Blue (NOD/ShiLtJ)}
#'   \item{E}{#00A0F0 - Light Blue (NZO/H1LtJ)}
#'   \item{F}{#00A000 - Green (CAST/EiJ)}
#'   \item{G}{#F00000 - Red (PWK/PhJ)}
#'   \item{H}{#9000E0 - Purple (WSB/EiJ)}
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
    H = '#9000E0' # WSB/EiJ
)

#' Alternative founder colors for DO/CC mice visualization
#'
#' @description
#' Alternative color scheme optimized for colorblind accessibility
#' and modern visualization standards.
#'
#' @format A named character vector with hex color codes:
#' \describe{
#'   \item{A}{#F0E442 - Yellow (A/J)}
#'   \item{B}{#555555 - Dark Gray (C57BL/6J)}
#'   \item{C}{#E69F00 - Orange (129S1/SvImJ)}
#'   \item{D}{#0072B2 - Blue (NOD/ShiLtJ)}
#'   \item{E}{#56B4E9 - Light Blue (NZO/H1LtJ)}
#'   \item{F}{#009E73 - Green (CAST/EiJ)}
#'   \item{G}{#D55E00 - Red-Orange (PWK/PhJ)}
#'   \item{H}{#CC79A7 - Pink (WSB/EiJ)}
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
    H = '#CC79A7' # WSB/EiJ
)

#' Default chromosome lengths for mouse genome
#'
#' @description
#' Standard chromosome lengths in base pairs for the mouse genome.
#' Used for chromosome ordering and positioning in plots.
#'
#' @format A named numeric vector with chromosome lengths:
#' \describe{
#'   \item{1-19}{Autosomal chromosomes with lengths in base pairs}
#'   \item{X}{X chromosome length}
#' }
#'
#' @export
DEFAULT_CHROM_LENGTHS <- c(
    '1' = 195154279,
    '2' = 181755017,
    '3' = 159745316,
    '4' = 156860686,
    '5' = 151758149,
    '6' = 149588044,
    '7' = 144995196,
    '8' = 130127694,
    '9' = 124359700,
    '10' = 130530862,
    '11' = 121973369,
    '12' = 120092757,
    '13' = 120883175,
    '14' = 125139656,
    '15' = 104073951,
    '16' = 98008968,
    '17' = 95294699,
    '18' = 90720763,
    '19' = 61420004,
     'X' = 169476592
)
