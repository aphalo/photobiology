#' Calculate (energy) irradiance from spectral irradiance
#'
#' Energy irradiance for a waveband from a radiation spectrum, optionally
#' applying a "biological spectral weighting function" or BSWF.
#'
#' @param w.length numeric vector of wavelength [\eqn{nm}].
#' @param s.irrad numeric vector of spectral irradiances in
#'   [\eqn{W\,m^{-2}\,nm^{-1}}{W m-2 nm-1}] or
#'   [\eqn{mol\,s^{-1}\,sm^{-2}\,nm^{-1}}{mol s-1 m-2 nm-1}] as indicated by the
#'   argument pased to \code{unit.in}.
#' @param w.band waveband.
#' @param unit.in character Allowed values \code{"energy"}, and \code{"photon"},
#'   or its alias \code{"quantum"}.
#' @param check.spectrum logical Flag indicating whether to sanity check input
#'   data, default is \code{TRUE}.
#' @param use.cached.mult logical Flag indicating whether multiplier values
#'   should be cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#'
#' @return A single numeric value with no change in scale factor:
#'   [\eqn{W\,m^{-2}}{W m-2}].
#'
#' @export
#' @examples
#' with(sun.data, energy_irradiance(w.length, s.e.irrad))
#' with(sun.data, energy_irradiance(w.length, s.e.irrad, new_waveband(400,700)))
#'
#' @family low-level functions operating on numeric vectors.
#'
energy_irradiance <-
  function(w.length, s.irrad,
           w.band = NULL,
           unit.in = "energy",
           check.spectrum = TRUE,
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges", default = NULL) ) {
    return(irradiance(w.length = w.length, s.irrad = s.irrad, w.band = w.band,
                      unit.out = "energy", unit.in = unit.in,
                      check.spectrum = check.spectrum,
                      use.cached.mult = use.cached.mult,
                      use.hinges = use.hinges))
  }
