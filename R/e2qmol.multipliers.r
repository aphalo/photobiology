#' Gives multipliers as a function of wavelength, for converting from energy to photon 
#' (quantum) molar units.
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
  Na = 6.02214129e23 # Avogadro's number, photons per mol
  return(e2quantum_multipliers(w.length) / Na)
}
