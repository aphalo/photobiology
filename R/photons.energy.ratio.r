#' Photon:energy ratio
#'
#' This function gives the photons:energy ratio between for one given waveband
#' of a radiation spectrum.
#'
#' @param w.length numeric vector of wavelength (nm).
#' @param s.irrad numeric vector of spectral irradiances in
#'   [\eqn{W\,m^{-2}\,nm^{-1}}{W m-2 nm-1}] or
#'   [\eqn{mol\,s^{-1}\,sm^{-2}\,nm^{-1}}{mol s-1 m-2 nm-1}] as indicated by the
#'   argument pased to \code{unit.in}.
#' @param w.band waveband object.
#' @param unit.in character Allowed values \code{"energy"}, and \code{"photon"},
#'   or its alias \code{"quantum"}.
#' @param check.spectrum logical Flag telling whether to sanity check input
#'   data, default is TRUE.
#' @param use.cached.mult logical Flag telling whether multiplier values should
#'   be cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#'
#' @note The default for the \code{w.band} parameter is a waveband covering
#'   the whole range of \code{w.length}.
#'
#' @return A single numeric value giving the ratio moles-photons per Joule.
#'
#' @export
#' @examples
#' # photons:energy ratio
#' with(sun.data, photons_energy_ratio(w.length, s.e.irrad, new_waveband(400,500)))
#' # photons:energy ratio for whole spectrum
#' with(sun.data, photons_energy_ratio(w.length, s.e.irrad))
#'
#' @family low-level functions operating on numeric vectors.
#'
photons_energy_ratio <-
  function(w.length, s.irrad,
           w.band = NULL,
           unit.in = "energy",
           check.spectrum = TRUE,
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges", default = NULL) ){
    return(waveband_ratio(w.length, s.irrad, w.band, w.band,
                          unit.out.num = "photon", unit.out.denom = "energy",
                          unit.in = unit.in,
                          check.spectrum = check.spectrum,
                          use.cached.mult = use.cached.mult,
                          use.hinges = use.hinges))
  }
