#' Convert spectral energy irradiance into spectral photon irradiance
#'
#' For example an spectrum [W m-2 nm-1] is converted into a spectrum [s-1 m-2 nm-1]
#'
#' @usage as_quantum(w.length, s.e.irrad)
#' 
#' @param w.length numeric array of wavelength (nm)
#' @param s.e.irrad a numeric array of spectral (energy) irradiances
#'
#' @return a numeric array of spectral photon irradiances
#' @export
#' @keywords manip misc
#' @examples
#' data(sun.data)
#' with(sun.data, as_quantum(w.length, s.e.irrad))

as_quantum <- function(w.length, s.e.irrad){
  return(s.e.irrad * e2quantum_multipliers(w.length))
}
