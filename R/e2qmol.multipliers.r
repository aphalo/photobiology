#' Gives multipliers as a function of wavelength, for converting from energy to photon 
#' (quantum) molar units.
#'
#' @usage e2qmol_multipliers(w.length)
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
  return(e2quantum_multipliers(w.length, molar=TRUE))
}
