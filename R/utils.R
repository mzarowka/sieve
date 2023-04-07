################################################################### Prepare data

# Convert units between metric and phi, and vice versa
convert_units <- function(size, class) {
  # Convert to phi
  if (class == "um") {
    size_phi <- -log2(size)

    # Convert to metric
  } else if (class == "phi") {
    size_um <- 2^(-size)
  }
}

# Find class midpoint
find_midpoint <- function(size) {
  midpoint <- dplyr::lag(size + dplyr::lead(size)) / 2
}

# Prepare tibble
prepare_tibble <- function(data, class) {
  if (class == "um") {
    data <- data |>
      dplyr::mutate(
        size_mm = size / 1000,
        size_phi = convert_units(size_mm, class = "um"),
        midpoint_mm = find_midpoint(size_mm),
        midpoint_phi = find_midpoint(size_phi),
        .after = size_mm
      )

    return(data)
  } else if (class == "phi") {
    data <- data |>
      dplyr::mutate(
        size_phi = size,
        size_mm = convert_units(size_phi, class = "phi"),
        midpoint_mm = find_midpoint(size_mm),
        midpoint_phi = find_midpoint(size_phi),
        .after = size
      )
  }
  return(data)
}
