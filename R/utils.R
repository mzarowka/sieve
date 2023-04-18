################################################################### Prepare data

# #' Convert units between metric and phi
#'
#' @param size vector (numeric) of sizes.
#' @param unit either "µm" or "phi" (character).
#'
#' @return
#' @export
#'
#' @examples
convert_units <- function(size, unit) {
  # Convert to phi
  if (unit == "um") {
    size_phi <- -log2(size)

    # Convert to metric
  } else if (unit == "phi") {
    size_um <- 2^(-size)
  }
}

# Find class midpoint
find_midpoint <- function(size) {
  midpoint <- dplyr::lag(size + dplyr::lead(size)) / 2
}

# Prepare tibble
prepare_tibble <- function(data, unit) {
  abundance_sum <- sum(data$abundance)

  if (unit == "um") {
    data <- data |>
      dplyr::mutate(
        size_mm = size / 1000,
        size_um = size,
        size_phi = convert_units(size_mm, unit = "um"),
        abundance_g = abundance,
        abundance_p = (abundance_g * 100) / abundance_sum,
        cum_retained_p = cumsum(abundance_p),
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
        abundance_g = abundance,
        abundance_p = (abundance_g * 100) / abundance_sum,
        cum_retained_p = cumsum(abundance_p),
        midpoint_mm = find_midpoint(size_mm),
        midpoint_um = midpoint_mm * 1000,
        midpoint_phi = find_midpoint(size_phi),
        .after = size,
        .keep = "none"
      )
  }
  return(data)
}
