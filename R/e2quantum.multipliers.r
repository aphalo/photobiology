#' Gives multipliers as a function of wavelength, for converting from energy to photon 
#' (quantum) units (number of photons as default, or moles of photons).
#'
#' @usage e2quantum_multipliers(w.length, molar=FALSE)
#' 
#' @param w.length numeric array of wavelengths (nm)
#' @param molar logical indicating whether output should be in moles or numbers
#' 
#' @return a numeric array of multipliers
#' @keywords manip misc
#' @export
#' @examples
#' with(sun.data, e2quantum_multipliers(w.length))
#' with(sun.data, e2quantum_multipliers(w.length, molar=TRUE))
#' 
e2quantum_multipliers <- function(w.length, molar=FALSE){
  # E = hc / w.length energy of one photon
  # coverting (energy) irradiance (I) to photon irradiance (Q): I / E = Q
  h = 6.626e-34 # Plank's contsant (Js)
  c = 2.998e8 # speed of light in vaccum (m/s), so we need to convert wavelength to metres
  Na = 6.02214129e23 # Avogadro's number, photons per mol
  
  hc = h * c * 1e9
  if (molar) hc = hc * Na
  return(w.length  / hc)
}
