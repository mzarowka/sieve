#' Calculate statistics using arithmetic method of moments
#'
#' @param data data prepared with \code{\link{psa_prepare_data}}
#'
#' @return a tibble with calculated statistics.
#' @export
psa_calculate_indices_amm <- function(data) {
  # Calculate first moment intermediates
  fm <- data$midpoint.um * data$abundance

  # Calculate first moment
  m_1 <- (sum(fm[fm >= 0], na.rm = TRUE) + sum(fm[fm < 0], na.rm = TRUE)) / 100

  # Calculate second moment intermediates
  fm2 <- data$abundance * (data$midpoint.um - m_1)^2

  # Calculate second moment
  m_2 <- sqrt((sum(fm2[fm2 >= 0], na.rm = TRUE) + sum(fm2[fm2 < 0], na.rm = TRUE)) / 100)

  # Calculate third moment intermediates
  fm3 <- data$abundance * (data$midpoint.um - m_1)^3

  # Calculate third moment
  m_3 <- (sum(fm3[fm3 >= 0], na.rm = TRUE) + sum(fm3[fm3 < 0], na.rm = TRUE)) / ((100 * (m_2^3)))

  # Calculate fourth moment intermediates
  fm4 <- data$abundance * (data$midpoint.um - m_1)^4

  # Calculate fourth moment
  m_4 <- (sum(fm4[fm4 >= 0], na.rm = TRUE) + sum(fm4[fm4 < 0], na.rm = TRUE)) / ((100 * (m_2^4)))

  # Create tibble with indices
  indices <- dplyr::tibble(
    mean_amm = m_1,
    standard_deviation_amm = m_2,
    skewness_amm = m_3,
    kurtosis_amm = m_4
    # kurtosis_3_amm = m_4 - 3
  ) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) round(x, digits = 3)))
}

#' Calculate statistics using logarithmic method of moments
#'
#' @param data data prepared with \code{\link{psa_prepare_data}}
#'
#' @return a tibble with calculated statistics.
#' @export
psa_calculate_indices_lmm <- function(data) {
  # Calculate first moment intermediates
  fm <- data$midpoint.phi * data$abundance

  # Calculate first moment
  m_1 <- (sum(fm[fm >= 0], na.rm = TRUE) + sum(fm[fm < 0], na.rm = TRUE)) / 100

  # Calculate second moment intermediates
  fm2 <- data$abundance * (data$midpoint.phi - m_1)^2

  # Calculate second moment
  m_2 <- sqrt((sum(fm2[fm2 >= 0], na.rm = TRUE) + sum(fm2[fm2 < 0], na.rm = TRUE)) / 100)

  # Calculate third moment intermediates
  fm3 <- data$abundance * (data$midpoint.phi - m_1)^3

  # Calculate third moment
  m_3 <- (sum(fm3[fm3 >= 0], na.rm = TRUE) + sum(fm3[fm3 < 0], na.rm = TRUE)) / ((100 * (m_2^3)))

  # Calculate fourth moment intermediates
  fm4 <- data$abundance * (data$midpoint.phi - m_1)^4

  # Calculate fourth moment
  m_4 <- (sum(fm4[fm4 >= 0], na.rm = TRUE) + sum(fm4[fm4 < 0], na.rm = TRUE)) / ((100 * (m_2^4)))

  # Create tibble with indices
  indices <- dplyr::tibble(
    mean_lmm = m_1,
    standard_deviation_lmm = m_2,
    skewness_lmm = m_3,
    kurtosis_lmm = m_4,
  ) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) round(x, digits = 3)))
}

#' Calculate statistics using geometric method of moments
#'
#' @param data data prepared with \code{\link{psa_prepare_data}}
#'
#' @return a tibble with calculated statistics.
#' @export
psa_calculate_indices_gmm <- function(data) {
  # Calculate log
  log_m <- log((2^-data$midpoint.phi) * 1000)

  # Calculate first moment intermediates
  fm <- log_m * data$abundance

  # Calculate first moment
  m_1 <- exp((sum(fm[fm >= 0], na.rm = TRUE) + sum(fm[fm < 0], na.rm = TRUE)) / 100)

  # Calculate second moment intermediates
  fm2 <- data$abundance * (log_m - log(m_1))^2

  # Calculate second moment
  m_2 <- exp(sqrt((sum(fm2[fm2 >= 0], na.rm = TRUE) + sum(fm2[fm2 < 0], na.rm = TRUE)) / 100))

  # Calculate third moment intermediates
  fm3 <- data$abundance * ((log_m - log(m_1))^3)

  # Calculate third moment
  m_3 <- (sum(fm3[fm3 >= 0], na.rm = TRUE) + sum(fm3[fm3 < 0], na.rm = TRUE)) / (100 * (log(m_2))^3)

  # Calculate fourth moment intermediates
  fm4 <- data$abundance * ((log_m - log(m_1))^4)

  # Calculate fourth moment
  m_4 <- (sum(fm4[fm4 >= 0], na.rm = TRUE) + sum(fm4[fm4 < 0], na.rm = TRUE)) / (100 * (log(m_2))^4)

  # Create tibble with indices
  indices <- dplyr::tibble(
    mean_gmm = m_1,
    standard_deviation_gmm = m_2,
    skewness_gmm = m_3,
    kurtosis_gmm = m_4
    # kurtosis_3_gmm = m_4 - 3
  ) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) round(x, digits = 3)))
}

#' Prepare data for Folk and Ward methods
#'
#' @param data data prepared with \code{\link{psa_prepare_data}}
#'
#' @return a tibble with calculated statistics.
#' @export
psa_get_quantiles <- function(data) {
  # Set vector of quantiles
  qunatiles <- c(0.1, 0.5, 10.0, 16.0, 25.0, 50.0, 75.0, 84.0, 90.0, 95.0, 99.0)

  # Get respective sizes
  data <- data |>
    # Select variables
    dplyr::select(size.um, size.phi, abundance_cum) |>
    # Add rows of quantiles
    dplyr::add_row(abundance_cum = qunatiles) |>
    # Arrange by abundance and size
    dplyr::arrange(abundance_cum, size.phi) |>
    # Interpolate sizes
    dplyr::mutate(dplyr::across(dplyr::contains("size"), \(x) zoo::na.approx(x, rule = 2))) |>
    # Filter by quantiles
    dplyr::filter(abundance_cum %in% qunatiles)

  # Return data
  return(data)

}

#' Calculate statistics using Folk and Ward logarithmic method
#'
#' @param data data prepared with \code{\link{psa_fw_prepare}}
#'
#' @return a tibble with calculated statistics.
#' @export
psa_calculate_indices_lfw <- function(data) {
  # Check if good data here

  # Prepare data for Folk and Ward
  data <- psa_get_quantiles(data)

  # Get vector of sizes and quantiles
  sizes <- data |>
    # Select phi and abundances
    dplyr::select(abundance_cum, size.phi) |>
    # Deframe
    tibble::deframe()

  # Calculate indices
  indices <- tibble::tibble(mean_lfw = (sizes["16"] + sizes["50"] + sizes["84"]) / 3,
                            standard_deviation_lfw = ((sizes["84"] - sizes["16"]) / 4) + ((sizes["95"] - sizes["0.5"]) / 6.6),
                            skewness_lfw = ((sizes["16"] + sizes["84"] - 2 * sizes["50"]) / (2 * (sizes["84"] - sizes["16"]))) + ((sizes["0.5"] + sizes["95"] - 2 * sizes["50"]) / (2 * (sizes["95"] - sizes["0.5"]))),
                            kurtosis_lfw = (sizes["95"] - sizes["0.5"]) / (2.44 * (sizes["75"] - sizes["25"]))
                            ) |>
    # All to numeric
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) as.numeric(x)))

  # Return indices
  return(indices)
}

#' Calculate statistics using Folk and Ward geometric method
#'
#' @param data data prepared with \code{\link{psa_prepare_data}}
#'
#' @return a tibble with calculated statistics.
#' @export
psa_calculate_indices_gfw <- function(data) {
  # Check if good data here

  # Prepare data for Folk and Ward
  data <- psa_get_quantiles(data)

  # Get vector of sizes and quantiles
  sizes <- data |>
    # Select phi and abundances
    dplyr::select(abundance_cum, size.um) |>
    # Deframe
    tibble::deframe()

  # Calculate indices
  indices <- tibble::tibble(mean_gfw = exp((log(sizes["16"]) + log(sizes["50"]) + log(sizes["84"])) / 3),
                            standard_deviation_gfw = exp(((log(sizes["16"]) - log(sizes["84"])) / 4) + ((log(sizes["0.5"]) - log(sizes["95"])) / 6.6)),
                            skewness_gfw = (log(sizes["16"]) + log(sizes["84"]) - 2 * (log(sizes["50"]))) / (2 * (log(sizes["84"]) - log(sizes["16"]))) + ((log(sizes["0.5"]) + log(sizes["95"]) - (2 * (log(sizes["50"])))) / (2 * (log(sizes["95"]) - log(sizes["0.5"])))),
                            kurtosis_gfw = (log(sizes["0.5"]) - log(sizes["95"])) / (2.44 * (log(sizes["25"]) - log(sizes["75"])))
                            ) |>
    # All to numeric
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) as.numeric(x)))

  # Return indices
  return(indices)
}

#' Calculate qunatile based phi diameter statistics
#'
#' @param data data prepared with \code{\link{psa_prepare_data}}
#'
#' @return a tibble with calculated statistics.
#' @export
psa_calculate_qunatiles_phi <- function(data) {
  # Check if good data here

  # Prepare data for Folk and Ward
  data <- psa_get_quantiles(data)

  # Get vector of sizes and quantiles
  sizes <- data |>
    # Select phi and abundances
    dplyr::select(abundance_cum, size.phi) |>
    # Deframe
    tibble::deframe()

  # Calculate indices
  indices <- tibble::tibble(d01.phi = sizes["0.1"],
                            d05.phi = sizes["0.5"],
                            d10.phi = sizes["10"],
                            d50.phi = sizes["50"],
                            d90.phi = sizes["90"],
                            ratio_d90d10.phi = sizes["90"] / sizes["10"],
                            diff_d90d10.phi = sizes["90"] - sizes["10"],
                            ratio_d75d25.phi = sizes["75"] / sizes["25"],
                            diff_d75d25.phi = sizes["75"] - sizes["25"]
                            ) |>
    # All to numeric
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) as.numeric(x)))

  # Return indices
  return(indices)
}

#' Calculate qunatile based metric diameter statistics
#'
#' @param data data prepared with \code{\link{psa_prepare_data}}
#'
#' @return a tibble with calculated statistics.
#' @export
psa_calculate_qunatiles_um <- function(data) {
  # Check if good data here

  # Prepare data for Folk and Ward
  data <- psa_get_quantiles(data)

  # Get vector of sizes and quantiles
  sizes <- data |>
    # Select phi and abundances
    dplyr::select(abundance_cum, size.phi) |>
    # Deframe
    tibble::deframe()

  # Calculate indices
  indices <- tibble::tibble(d01.um = psa_convert_units(sizes["0.1"], unit = "phi"),
                            d05.um = psa_convert_units(sizes["0.5"], unit = "phi"),
                            d10.um = psa_convert_units(sizes["10"], unit = "phi"),
                            d50.um = psa_convert_units(sizes["50"], unit = "phi"),
                            d90.um = psa_convert_units(sizes["90"], unit = "phi"),
                            ratio_d90d10.um = psa_convert_units(sizes["90"], unit = "phi") / psa_convert_units(sizes["10"], unit = "phi"),
                            diff_d90d10.um = psa_convert_units(sizes["90"], unit = "phi") - psa_convert_units(sizes["10"], unit = "phi"),
                            ratio_d75d25.um = psa_convert_units(sizes["75"], unit = "phi") / psa_convert_units(sizes["25"], unit = "phi"),
                            diff_d75d25.um = psa_convert_units(sizes["75"], unit = "phi") - psa_convert_units(sizes["25"], unit = "phi")
  ) |>
    # All to numeric
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) as.numeric(x)))

  # Return indices
  return(indices)
}

#' Calculate base size indices
#'
#' @param data a tibble prepared with \code{\link{psa_prepare_data}}
#' @param method character, one of: arithmetic method of moments ("amm"), logarithmic method of moments ("lmm"), geometric method of moments ("gmm"), Folk and Ward logarithmic ("lfw"), Folk and Ward geometric ("gfw"), metric quantiles ("umq"), phi quantiles ("phiq") or all ("all" - default).
#'
#' @return a data frame with calculated indices.
#' @export
#'
psa_calculate_indices <- function(data, method) {
  # Check method
  if (!(method %in% c("all", "amm", "lmm", "gmm", "lfw", "gfw", "umq", "phiq"))) {
    rlang::abort(cli::format_error("Method not recognized"))

    # Calculate for all
  } else if (method == "all") {
    # Calculate indices by type
    indices_amm <- psa_calculate_indices_amm(data)
    indices_gmm <- psa_calculate_indices_gmm(data)
    indices_lmm <- psa_calculate_indices_lmm(data)
    indices_gfw <- psa_calculate_indices_gfw(data)
    indices_lfw <- psa_calculate_indices_lfw(data)
    indices_umq <- psa_calculate_qunatiles_um(data)
    indices_phiq <- psa_calculate_qunatiles_phi(data)

    # Bind all types
    indices <- dplyr::bind_cols(indices_amm, indices_gmm, indices_lmm, indices_gfw, indices_lfw, indices_umq, indices_phiq)

  } else if (method == "amm") {
    indices <- psa_calculate_indices_amm(data)
  } else if (method == "gmm") {
    indices <- psa_calculate_indices_gmm(data)
  } else if (method == "lmm") {
    indices <- psa_calculate_indices_lmm(data)
  } else if (method == "lfw") {
    indices <- psa_calculate_indices_lfw(data)
  } else if (method == "gfw") {
    indices <- psa_calculate_indices_gfw(data)
  } else if (method == "umq") {
    indices <- psa_calculate_qunatiles_um(data)
  } else if (method == "phiq") {
    indices <- psa_calculate_qunatiles_phi(data)
  }

  # Return indices
  return(indices)
}

#' Get all indices
#'
#' @param data nested data prepared with \code{\link{prepare_tibble}}
#' @param method character, one of: arithmetic method of moments ("amm"), logarithmic method of moments ("lmm"), geometric method of moments ("gmm"), Folk and Ward logarithmic ("lfw"), Folk and Ward geometric ("gfw"), metric quantiles ("umq"), phi quantiles ("phiq") or all ("all" - default).
#'
#' @return a data frame with calculated indices.
#' @export
psa_get_indices <- function(data, method = "all"){
  # Check method
  if (!(method %in% c("all", "amm", "lmm", "gmm", "lfw", "gfw", "umq", "phiq"))) {
    rlang::abort(cli::format_error("Method not recognized"))
  }

  # Calculate indices
  indices <- data |>
    # Calculate indices
    dplyr::mutate(indices = purrr::map(data, \(x) psa_calculate_indices(x, method = method))) |>
    # Drop data
    dplyr::select(-data) |>
    # Unnest
    tidyr::unnest(cols = indices)
}
