#' Convert spectral energy irradiance into spectral photon irradiance
#'
#' Convert spectral energy irradiance [\eqn{W\,m^{-2}\,nm^{-1}}{W m-2 nm-1}]
#' into a spectral photon irradiance expressed in number of molds of photons
#' [\eqn{mol\,s^{-1}\,m^{-2}\,nm^{-1}}{mol s-1 m-2 nm-1}].
#'
#' @param w.length numeric vector of wavelengths (nm).
#' @param s.e.irrad numeric vector of spectral (energy) irradiance values.
#'
#' @return a numeric vector of spectral photon irradiances.
#' @export
#'
#' @examples
#' with(sun.data, as_quantum_mol(w.length, s.e.irrad))
#'
#' @family low-level functions operating on numeric vectors.
#'
as_quantum_mol <- function(w.length, s.e.irrad) {
  s.e.irrad * e2qmol_multipliers(w.length)
}
