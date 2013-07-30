#' Convert spectral photon irradiance into spectral energy irradiance
#'
#' For example an spectrum [mol s-1 m-2 nm-1] is converted into a spectrum [W m-2 nm-1] 
#'
#' @param w.length numeric array of wavelength (nm)
#' @param s.qmol.irrad a numeric array of spectral photon irradiances
#'
#' @return a numeric array of spectral (energy) irradiances
#' @export
#' @keywords manip misc
#' @examples
#' data(sun.data)
#' with(sun.data, as_energy(w.length, s.q.irrad))
#' 

as_energy <- function(w.length, s.qmol.irrad){
  return(s.qmol.irrad / e2qmol_multipliers(w.length))
}
