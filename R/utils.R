################################################################### Prepare data

# #' Convert units between metric and phi
#'
#' @param size vector (numeric) of sizes.
#' @param unit either "µm" or "phi" (character).
#'
#' @return a vector of converted sizes.
#' @export
convert_units <- function(size, unit) {
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
find_midpoint <- function(size) {
  midpoint <- dplyr::lag(size + dplyr::lead(size)) / 2
}

# Check tibble: arrange, size as a first column, all numeric

#' Prepare sieve tibble
#'
#' @param data a tibble with size and abundances.
#' @param unit original unit of measurements. One of "um" or "phi".
#'
#' @return a tibble prepared for further analyses with \pkg{sieve}.
#' @export
prepare_tibble <- function(data, unit) {
  if (!(unit %in% c("um", "phi"))) {
    rlang::abort(cli::format_error("Size unit not recognized"))
  }

  # Calculate sum
  sample_sum <- sum(data$sample)

  if (unit == "um") {
    data <- data |>
      dplyr::mutate(
        size_mm = size / 1000,
        size_um = size,
        size_phi = convert_units(size_mm, unit = "um"),
        sample_g = sample,
        sample_p = (sample_g * 100) / sample_sum,
        cum_retained_p = cumsum(sample_p),
        midpoint_mm = find_midpoint(size_mm),
        midpoint_um = midpoint_mm * 1000,
        midpoint_phi = find_midpoint(size_phi),
        .after = size_mm,
        .keep = "none"
      )

    return(data)
  } else if (unit == "phi") {
    data <- data |>
      dplyr::mutate(
        size_phi = size,
        size_mm = convert_units(size_phi, unit = "phi"),
        size_um = size_mm * 1000,
        sample_g = sample,
        sample_p = (sample_g * 100) / sample_sum,
        cum_retained_p = cumsum(sample_p),
        midpoint_mm = find_midpoint(size_mm),
        midpoint_um = midpoint_mm * 1000,
        midpoint_phi = find_midpoint(size_phi),
        .after = size,
        .keep = "none"
      )
  }
  return(data)
}
