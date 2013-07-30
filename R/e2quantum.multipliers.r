#' Gives multipliers as a function of wavelength, for converting from energy to photon 
#' (quantum) units.
#'
#' @param w.length numeric array of wavelengths (nm)
#' 
#' @return a numeric array of multipliers
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' with(sun.data, e2quantum_multipliers(w.length))
#' 
e2quantum_multipliers <- function(w.length){
  # E = hc / w.length energy of one photon
  # covertine (energy) irradiance (I) to photon flux (Q): I / E = Q
  h = 6.626e-34 # Plank's contsant (Js)
  c = 2.998e8 # speed of light in vaccum (m/s), so we need to convert wavelength to metres
  hc = h * c
  return((w.length * 1e-9) / hc)
}
