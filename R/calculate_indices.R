################################################## Arithmetic method of moments

# Function to calculate indices according to the arithmetic method of moments
calculate_indices_amm <- function(data) {

}

# Function to calculate mean size according to the arithmetic method of moments
calculate_mean_amm <- function(data) {
  data <- data |>
    dplyr::mutate()
  mean_amm <- sum(data$abundance) * sum(data$size) / 100
}

