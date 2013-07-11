#' Gives multipliers as a function of wavelength, for converting from energy to photon 
#' (quantum) molar unoits.
#'
#' @param w.length numeric array of wavelengths (nm)
#' 
#' @return a numeric array of multipliers
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' with(sun.data, e2qmol_multipliers(w.length))
#' 
e2qmol_multipliers <- function(w.length){
  h = 6.626e-34 # Plank's contsant (Js)
  c = 2.998e8 # speed of light in vaccum (m/s), so we need to convert wavelength to metres
  Na = 6.022e23 # Avogadro's number, photons per mol
  return((w.length * 1e-9) / (h * c * Na))
}
