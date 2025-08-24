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
        message('Updated ', param, ': ', old_value, ' -> ', updates[[param]])
    }

    return(config)
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
#' 'Arial' %in% list_available_fonts()
#' }
#'
#' @export
list_available_fonts <- function() {
    #cCheck if extrafont is available
    if (requireNamespace('extrafont', quietly = TRUE)) {
        message('Using extrafont to list fonts...')
        extrafont::font_import(prompt = FALSE)
        fonts <- extrafont::fonts()
        return(fonts)
    }

    # check if systemfonts is available
    if (requireNamespace('systemfonts', quietly = TRUE)) {
        message('Using systemfonts to list fonts...')
        fonts <- systemfonts::system_fonts()$family
        return(unique(fonts))
    }

    # fallback: common fonts that are usually available
    message('Using fallback method - common fonts:')
    common_fonts <- c(
        'sans', 'serif', 'mono',
        'Arial', 'Times', 'Courier',
        'Helvetica', 'Georgia', 'Verdana',
        'Palatino', 'Bookman', 'Avant Garde'
    )
    return(common_fonts)
}


