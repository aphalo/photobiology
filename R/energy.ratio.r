#' Energy:energy ratio
#'
#' Energy irradiance ratio between two wavebands for a radiation spectrum.
#'
#' @param w.length numeric vector of wavelengths (nm).
#' @param s.irrad numeric vector of spectral (energy) irradiances (W m-2
#'   nm-1).
#' @param w.band.num waveband object used to compute the numerator of the ratio.
#' @param w.band.denom waveband object used to compute the denominator of the ratio.
#' @param unit.in character Allowed values "energy", and "photon", or its alias
#'   "quantum".
#' @param check.spectrum logical Flag indicating whether to sanity check input
#'   data, default is TRUE.
#' @param use.cached.mult logical Flag indicating whether multiplier values
#'   should be cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#'
#' @note The default for both \code{w.band} parameters is a waveband covering
#'   the whole range of \code{w.length}.
#'
#' @return a single numeric value giving the unitless ratio.
#'
#' @export
#' @examples
#' with(sun.data,
#'      energy_ratio(w.length, s.e.irrad, new_waveband(400,500), new_waveband(400,700)))
#'
#' @family low-level functions operating on numeric vectors.
#'
energy_ratio <- function(w.length, s.irrad,
                         w.band.num = NULL, w.band.denom = NULL,
                         unit.in = "energy",
                         check.spectrum = TRUE,
                         use.cached.mult = FALSE,
                         use.hinges = NULL) {
  return(waveband_ratio(w.length, s.irrad, w.band.num, w.band.denom,
                        unit.out.num = "energy", unit.out.denom = "energy",
                        unit.in = unit.in,
                        check.spectrum = check.spectrum,
                        use.cached.mult = use.cached.mult,
                        use.hinges = use.hinges))
}
