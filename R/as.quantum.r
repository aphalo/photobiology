#' Convert spectral energy irradiance into spectral photon irradiance
#'
#' Convert spectral energy irradiance [W m-2 nm-1] into spectral photon
#' irradiance expressed as number of photons [s-1 m-2 nm-1]
#'
#' @param w.length numeric vector of wavelengths (nm).
#' @param s.e.irrad numeric vector of spectral (energy) irradiance values.
#'
#' @return A numeric vector of spectral photon irradiances.
#' @export
#'
#' @examples
#' with(sun.data, as_quantum(w.length, s.e.irrad))
#'
#' @family quantity conversion functions
#'
as_quantum <- function(w.length, s.e.irrad){
  return(s.e.irrad * e2quantum_multipliers(w.length))
}
