#' Plot cumulative distributions
#'
#' @param data data to plot
#' @param unit unit of sizes
#'
#' @return plot(s)
#' @export
psa_plot_cumulative <- function(data, unit = "phi") {
  # Set theme
  ggplot2::theme_set(theme(panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black")))

  if (unit == "phi") {
    plot_cumulative <- data |>
      ggplot2::ggplot(aes(size.phi, abundance_cum)) +
      geom_line() +
      scale_y_continuous(name = "Cumulative abundance (%)", limits = c(0, 100), labels = scales::label_number(accuracy = 0.01)) +
      scale_x_continuous(name = "Size (phi)")
  } else if (unit == "um") {
    plot_cumulative <- data |>
      ggplot2::ggplot(aes(size.um, abundance_cum)) +
      geom_line() +
      scale_y_continuous(name = "Cumulative abundance (%)", limits = c(0, 100), labels = scales::label_number(accuracy = 0.01)) +
      scale_x_log10(name = "Size (µm)", labels = scales::label_number())
  }
}

#' Plot particle size distributions
#'
#' @param data data to plot.
#' @param unit unit of sizes.
#'
#' @return particle size distribution plot.
#' @export
psa_plot_distribution <- function(data, unit = "phi") {
  # Set theme
  ggplot2::theme_set(theme(panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black")))

  if (unit == "phi") {
    plot_cumulative <- data |>
      ggplot2::ggplot(aes(size.phi, abundance)) +
      geom_line() +
      scale_y_continuous(name = "Abundance (%)", labels = scales::label_number(accuracy = 0.01)) +
      scale_x_continuous(name = "Size (phi)")
  } else if (unit == "um") {
    plot_cumulative <- data |>
      ggplot2::ggplot(aes(size.um, abundance)) +
      geom_line() +
      scale_y_continuous(name = "Abundance (%)", labels = scales::label_number(accuracy = 0.01)) +
      scale_x_log10(name = "Size (µm)", labels = scales::label_number())
  }
}

#' Plot ternary diagrams
#'
#' @param data data.
#' @param triangle which triangle to use.
#'
#' @return ternary diagram of choice.
#' @export
psa_plot_ternary <- function(data, triangle){

}

#' Plot variable plot types
#'
#' @param data data
#' @param type which plot
#' @param ... optional arguments
#'
#' @return plot.
#' @export
psa_plot <- function(data, type, ...) {
  # Checks

  # Access the additional arguments
  dots <- list(...)

  # Get data from nested tibble
  data <- data$data

  if (type == "cumulative") {
    plots <- purrr::map(data, \(x) psa_plot_cumulative(x, unit = dots$unit))
  } else if (type =="distribution") {
    plots <- purrr::map(data, \(x) psa_plot_distribution(x, unit = dots$unit))
  } else if (type == "ternary"){
    plots <- purrr::map(data, \(x) psa_plot_ternary(x, triangle = dots$triangle))
  }

  # Return plots
  return(plots)
}
