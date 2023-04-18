#### Arithmetic method of moments ##############################################
calculate_indices_amm <- function(data) {
  # Calculate first moment intermediates
  fm <- data$midpoint_um * data$abundance_p

  # Calculate first moment
  m_1 <- (sum(fm[fm > 0], na.rm = TRUE) + sum(fm[fm < 0], na.rm = TRUE)) / 100

  # Calculate second moment intermediates
  fm2 <- data$abundance_p * (data$midpoint_um - m_1)^2

  # Calculate second moment
  m_2 <- sqrt((sum(fm2[fm2 > 0], na.rm = TRUE) + sum(fm2[fm2 < 0], na.rm = TRUE)) / 100)

  # Calculate third moment intermediates
  fm3 <- data$abundance_p * (data$midpoint_um - m_1)^3

  # Calculate third moment
  m_3 <- (sum(fm3[fm3 > 0], na.rm = TRUE) + sum(fm3[fm3 < 0], na.rm = TRUE)) / ((100 * (m_2^3)))

  # Calculate fourth moment intermediates
  fm4 <- data$abundance_p * (data$midpoint_um - m_1)^4

  # Calculate fourth moment
  m_4 <- (sum(fm4[fm4 > 0], na.rm = TRUE) + sum(fm4[fm4 < 0], na.rm = TRUE)) / ((100 * (m_2^4)))

  indices <- dplyr::tibble(
    mean = m_1,
    standard_deviation = m_2,
    skewness = m_3,
    kurtosis = m_4,
    kurtosis_3 = m_4 - 3
  ) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) round(x, digits = 3)))
}



#### Logarithmic method of moments #############################################
calculate_indices_lmm <- function(data) {
  # Calculate first moment intermediates
  fm <- data$midpoint_phi * data$abundance_p

  # Calculate first moment
  m_1 <- (sum(fm[fm > 0], na.rm = TRUE) + sum(fm[fm < 0], na.rm = TRUE)) / 100

  # Calculate second moment intermediates
  fm2 <- data$abundance_p * (data$midpoint_phi - m_1)^2

  # Calculate second moment
  m_2 <- sqrt((sum(fm2[fm2 > 0], na.rm = TRUE) + sum(fm2[fm2 < 0], na.rm = TRUE)) / 100)

  # Calculate third moment intermediates
  fm3 <- data$abundance_p * (data$midpoint_phi - m_1)^3

  # Calculate third moment
  m_3 <- (sum(fm3[fm3 > 0], na.rm = TRUE) + sum(fm3[fm3 < 0], na.rm = TRUE)) / ((100 * (m_2^3)))

  # Calculate fourth moment intermediates
  fm4 <- data$abundance_p * (data$midpoint_phi - m_1)^4

  # Calculate fourth moment
  m_4 <- (sum(fm4[fm4 > 0], na.rm = TRUE) + sum(fm4[fm4 < 0], na.rm = TRUE)) / ((100 * (m_2^4)))

  indices <- dplyr::tibble(
    mean = m_1,
    standard_deviation = m_2,
    skewness = m_3,
    kurtosis = m_4,
  ) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) round(x, digits = 3)))
}

#### Geometric method of moments ###############################################
calculate_indices_gmm <- function(data) {
  # Calculate log
  log_m <- log10((2^-(data$midpoint_phi)) * 1000)

  # Calculate first moment intermediates
  fm <- log_m * data$abundance_p

  # Calculate first moment
  m_1 <- 10^((sum(fm[fm > 0], na.rm = TRUE) + sum(fm[fm < 0], na.rm = TRUE)) / 100)

  # Calculate second moment intermediates
  fm2 <- data$abundance_p * (log_m - log10(fm))^2

  # Calculate second moment
  m_2 <- 10^(sqrt((sum(fm2[fm2 > 0], na.rm = TRUE) + sum(fm2[fm2 < 0], na.rm = TRUE)) / 100))

  # Calculate third moment intermediates
  fm3 <- data$abundance_p * ((log_m - log10(fm))^3)

  # Calculate third moment
  m_3 <- (sum(fm3[fm3 > 0], na.rm = TRUE) + sum(fm3[fm3 < 0], na.rm = TRUE)) / (100 * (log10(m_2))^3)

  # Calculate fourth moment intermediates
  fm4 <- data$abundance_p * ((log_m - log10(fm))^4)

  # Calculate fourth moment
  m_4 <- (sum(fm4[fm4 > 0], na.rm = TRUE) + sum(fm4[fm4 < 0], na.rm = TRUE)) / (100 * (log10(m_2))^4)

  indices <- dplyr::tibble(
    mean = m_1,
    standard_deviation = m_2,
    skewness = m_3,
    kurtosis = m_4,
    kurtosis_3 = m_4 - 3
  ) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) round(x, digits = 3)))
}

# Calculate indices
calculate_indices <- function(data, method) {
  if (method == "amm") {
    indices <- calculate_indices_amm(data)
  } else if (method == "gmm") {
    indices <- calculate_indices_gmm(data)
  } else if (method == "lmm") {
    indices <- calculate_indices_lmm(data)
  }
}
