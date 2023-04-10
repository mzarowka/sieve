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
  if (unit == "um") {
    data <- data |>
      dplyr::mutate(
        size_mm = size / 1000,
        size_phi = convert_units(size_mm, unit = "um"),
        midpoint_mm = find_midpoint(size_mm),
        midpoint_phi = find_midpoint(size_phi),
        .after = size_mm
      )

    return(data)
  } else if (unit == "phi") {
    data <- data |>
      dplyr::mutate(
        size_phi = size,
        size_mm = convert_units(size_phi, unit = "phi"),
        midpoint_mm = find_midpoint(size_mm),
        midpoint_phi = find_midpoint(size_phi),
        .after = size
      )
  }
  return(data)
}
