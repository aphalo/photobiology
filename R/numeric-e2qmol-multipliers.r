#' Calculate energy to quantum (mol) multipliers
#'
#' Multipliers as a function of wavelength, for converting from energy to
#' photon (quantum) molar units.
#'
#' @param w.length numeric Vector of wavelengths (nm)
#'
#' @return A numeric vector of multipliers
#'
#' @export
#' @examples
#' with(sun.data, e2qmol_multipliers(w.length))
#'
#' @family quantity conversion functions
#'
e2qmol_multipliers <- function(w.length) {
  e2quantum_multipliers(w.length, molar = TRUE)
}
