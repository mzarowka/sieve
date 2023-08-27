#### Arithmetic method of moments ##############################################
calculate_indices_amm <- function(data) {
  # Calculate first moment intermediates
  fm <- data$midpoint.um * data$abundance

  # Calculate first moment
  m_1 <- (sum(fm[fm > 0], na.rm = TRUE) + sum(fm[fm < 0], na.rm = TRUE)) / 100

  # Calculate second moment intermediates
  fm2 <- data$abundance * (data$midpoint.um - m_1)^2

  # Calculate second moment
  m_2 <- sqrt((sum(fm2[fm2 > 0], na.rm = TRUE) + sum(fm2[fm2 < 0], na.rm = TRUE)) / 100)

  # Calculate third moment intermediates
  fm3 <- data$abundance * (data$midpoint.um - m_1)^3

  # Calculate third moment
  m_3 <- (sum(fm3[fm3 > 0], na.rm = TRUE) + sum(fm3[fm3 < 0], na.rm = TRUE)) / ((100 * (m_2^3)))

  # Calculate fourth moment intermediates
  fm4 <- data$abundance * (data$midpoint.um - m_1)^4

  # Calculate fourth moment
  m_4 <- (sum(fm4[fm4 > 0], na.rm = TRUE) + sum(fm4[fm4 < 0], na.rm = TRUE)) / ((100 * (m_2^4)))

  indices <- dplyr::tibble(
    mean_amm = m_1,
    standard_deviation_amm = m_2,
    skewness_amm = m_3,
    kurtosis_amm = m_4,
    kurtosis_3_amm = m_4 - 3
  ) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) round(x, digits = 3)))
}



#### Logarithmic method of moments #############################################
calculate_indices_lmm <- function(data) {
  # Calculate first moment intermediates
  fm <- data$midpoint.phi * data$abundance

  # Calculate first moment
  m_1 <- (sum(fm[fm > 0], na.rm = TRUE) + sum(fm[fm < 0], na.rm = TRUE)) / 100

  # Calculate second moment intermediates
  fm2 <- data$abundance * (data$midpoint.phi - m_1)^2

  # Calculate second moment
  m_2 <- sqrt((sum(fm2[fm2 > 0], na.rm = TRUE) + sum(fm2[fm2 < 0], na.rm = TRUE)) / 100)

  # Calculate third moment intermediates
  fm3 <- data$abundance * (data$midpoint.phi - m_1)^3

  # Calculate third moment
  m_3 <- (sum(fm3[fm3 > 0], na.rm = TRUE) + sum(fm3[fm3 < 0], na.rm = TRUE)) / ((100 * (m_2^3)))

  # Calculate fourth moment intermediates
  fm4 <- data$abundance * (data$midpoint.phi - m_1)^4

  # Calculate fourth moment
  m_4 <- (sum(fm4[fm4 > 0], na.rm = TRUE) + sum(fm4[fm4 < 0], na.rm = TRUE)) / ((100 * (m_2^4)))

  indices <- dplyr::tibble(
    mean_lmm = m_1,
    standard_deviation_lmm = m_2,
    skewness_lmm = m_3,
    kurtosis_lmm = m_4,
  ) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) round(x, digits = 3)))
}

#### Geometric method of moments ###############################################
calculate_indices_gmm <- function(data) {
  # Calculate log
  log_m <- log10((2^-(data$midpoint.phi)) * 1000)

  # Calculate first moment intermediates
  fm <- log_m * data$abundance

  # Calculate first moment
  m_1 <- 10^((sum(fm[fm > 0], na.rm = TRUE) + sum(fm[fm < 0], na.rm = TRUE)) / 100)

  # Calculate second moment intermediates
  fm2 <- data$abundance * (log_m - log10(fm))^2

  # Calculate second moment
  m_2 <- 10^(sqrt((sum(fm2[fm2 > 0], na.rm = TRUE) + sum(fm2[fm2 < 0], na.rm = TRUE)) / 100))

  # Calculate third moment intermediates
  fm3 <- data$abundance * ((log_m - log10(fm))^3)

  # Calculate third moment
  m_3 <- (sum(fm3[fm3 > 0], na.rm = TRUE) + sum(fm3[fm3 < 0], na.rm = TRUE)) / (100 * (log10(m_2))^3)

  # Calculate fourth moment intermediates
  fm4 <- data$abundance * ((log_m - log10(fm))^4)

  # Calculate fourth moment
  m_4 <- (sum(fm4[fm4 > 0], na.rm = TRUE) + sum(fm4[fm4 < 0], na.rm = TRUE)) / (100 * (log10(m_2))^4)

  indices <- dplyr::tibble(
    mean_gmm = m_1,
    standard_deviation_gmm = m_2,
    skewness_gmm = m_3,
    kurtosis_gmm = m_4,
    kurtosis_3_gmm = m_4 - 3
  ) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) round(x, digits = 3)))
}

#' Calculate base size indices
#'
#' @param data a tibble prepared with \code{\link{prepare_tibble}}
#' @param method one of arithmetic method of moments ("amm"), logarithmic method of moments ("lmm") or Geometric method of moments ("gmm") or all (default).
#'
#' @return a data frame with calculated indices.
#' @export
#'
psa_calculate_indices <- function(data, method = "all") {
  # Check method
  if (!(method %in% c("all", "amm", "lmm", "gmm"))) {
    rlang::abort(cli::format_error("Method not recognized"))

    # Calculate for all
  } else if (method == "all") {
    # Calculate indices by type
    indices_amm <- calculate_indices_amm(data)
    indices_lmm <- calculate_indices_lmm(data)
    indices_gmm <- calculate_indices_gmm(data)

    # Bind all types
    indices <- dplyr::bind_cols(indices_amm, indices_lmm, indices_gmm)

  } else if (method == "amm") {
    indices <- calculate_indices_amm(data)
  } else if (method == "gmm") {
    indices <- calculate_indices_gmm(data)
  } else if (method == "lmm") {
    indices <- calculate_indices_lmm(data)
  }
}

#' Get all indices
#'
#' @param data nested data prepared with \code{\link{prepare_tibble}}
#' @param method one of arithmetic method of moments ("amm"), logarithmic method of moments ("lmm") or Geometric method of moments ("gmm") or all (default).
#'
#' @return a data frame with calculated indices.
#' @export
#'
#' @examples
psa_get_indices <- function(data, method = "all"){
  # Check method
  if (!(method %in% c("all", "amm", "lmm", "gmm"))) {
    rlang::abort(cli::format_error("Method not recognized"))
  }

  # Calculate indices
  indices <- data |>
    # Calculate indices
    dplyr::mutate(indices = purrr::map(data, \(x) psa_calculate_indices(x))) |>
    # Drop data
    dplyr::select(-data) |>
    # Unnest
    tidyr::unnest(cols = indices)
}
