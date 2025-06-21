#' Convert spectral photon irradiance into spectral energy irradiance
#'
#' Convert a spectral photon irradiance [\eqn{mol\,s^{-1}\,m^{-2}\,nm^{-1}}{mol
#' s-1 m-2 nm-1}] into a spectral energy irradiance [\eqn{W\,m^{-2}\,nm^{-1}}{W
#' m-2 nm-1}].
#'
#' @param w.length numeric vector of wavelengths [\eqn{nm}]).
#' @param s.qmol.irrad numeric vector of spectral photon irradiance values.
#'
#' @return A numeric vector of spectral (energy) irradiances.
#'
#' @export
#'
#' @examples
#' with(sun.spct, as_energy(w.length, s.q.irrad))
#'
#' @family low-level functions operating on numeric vectors.
#'
as_energy <- function(w.length, s.qmol.irrad) {
  s.qmol.irrad / e2qmol_multipliers(w.length)
}
