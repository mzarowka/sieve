#' Convert units between metric and phi
#'
#' @param size vector (numeric) of sizes.
#' @param unit either "µm" or "phi" (character).
#'
#' @return a vector of converted sizes.
#' @export
psa_convert_units <- function(size, unit) {
  if (!is.numeric(size)) {
    rlang::abort(cli::format_error("Size is not numeric."))
  }

  if (!(unit %in% c("um", "phi"))) {
    rlang::abort(cli::format_error("Size unit not recognized."))
  }

  # Convert to phi
  if (unit == "um") {
    size_phi <- -log2(size)

    # Convert to metric
  } else if (unit == "phi") {
    size_um <- 2^(-size)
  }
}

#' Find class midpoint
#'
#' @param size vector of sizes.
#'
#' @return a vector of midpoint sizes.
#' @export
psa_find_midpoint <- function(size) {
  midpoint <- dplyr::lag(size + dplyr::lead(size)) / 2
}

#' Prepare and nest particle size data
#'
#' Add secondary size unit, calculate midpoints and cumulative abundance. Size class must come in either µm or phi.
#' All abundance data must come as percent (%).
#'
#' @param data a tibble with size and abundances.
#' @param unit original unit of measurements. One of "um" (default) or "phi".
#'
#' @return A nested tibble with data prepared for further analyses with \pkg{psa}.
#' @export
#'
#' @examples psa_prepare_data(multiple_sample, unit = "um")
psa_prepare_data <- function(data, unit = "um") {
  # Check size unit
  if (!(unit %in% c("um", "phi"))) {
    rlang::abort(cli::format_error("Size unit not recognized."))
  }

  # Check if all data are numeric
  if (!(all(sapply(data, FUN = is.numeric)))) {
    rlang::abort(cli::format_error("There are non-numeric variables."))
  }

  if (unit == "um") {
    # Calculate size variables
    size <- data |>
      dplyr::mutate(
        size.mm = size / 1000,
        size.um = size,
        size.phi = psa_convert_units(size.mm, unit = "um"),
        midpoint.mm = psa_find_midpoint(size.mm),
        midpoint.um = midpoint.mm * 1000,
        midpoint.phi = psa_find_midpoint(size.phi),
        .after = size.mm,
        .keep = "none"
      )

  } else if (unit == "phi") {
    # Calculate size variables
    size <- data |>
      dplyr::mutate(
        size.phi = size,
        size.mm = psa_convert_units(size.phi, unit = "phi"),
        size.um = size.mm * 1000,
        midpoint.mm = psa_find_midpoint(size.mm),
        midpoint.um = midpoint.mm * 1000,
        midpoint.phi = fpsa_ind_midpoint(size.phi),
        .after = size,
        .keep = "none"
      )
  }

  # Keep only abundance data
  sample <- data |>
    # Select only sample
    dplyr::select(-size)

  # Bind back
  data <- dplyr::bind_cols(size, sample) |>
    # Arrange ascending
    dplyr::arrange(size) |>
    # Pivot longer
    tidyr::pivot_longer(-c(contains("size"), contains("midpoint")),
                        names_to = "sample",
                        values_to = "abundance"
    ) |>
    # Nest by sample
    tidyr::nest(.by = sample) |>
    # Add cumulative abundance
    dplyr::mutate(data = purrr::map(data, \(x) dplyr::mutate(x, dplyr::across(abundance, \(x) cumsum(x), .names = "{.col}_cum.p"))))

  # Return data
  return(data)
}
