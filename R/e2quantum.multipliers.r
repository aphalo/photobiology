#' Calculate energy to quantum multipliers
#'
#' Gives multipliers as a function of wavelength, for converting from energy to
#' photon (quantum) units (number of photons as default, or moles of photons).
#'
#' @param w.length numeric Vector of wavelengths (nm)
#' @param molar logical Flag indicating whether output should be in moles or
#'   numbers
#'
#' @return A numeric vector of multipliers
#'
#' @export
#' @examples
#' with(sun.data, e2quantum_multipliers(w.length))
#' with(sun.data, e2quantum_multipliers(w.length, molar = TRUE))
#'
#' @family quantity conversion functions
#'
e2quantum_multipliers <- function(w.length, molar = FALSE){
  # E = hc / w.length energy of one photon
  # converting (energy) irradiance (I) to photon irradiance (Q): I / E = Q
  h <- 6.62607015e-34 # Plank's constant (Js)
  c <- 2.99792458e8 # speed of light in vacuum (m/s), so we convert m/s to nm/s
  Na <- 6.02214129e23 # Avogadro's number, photons per mol

  # compute the constants to be used
  hc <- h * c * 1e9
  if (molar) hc <- hc * Na
  # do vectorization if needed
  return(w.length  / hc)
}
