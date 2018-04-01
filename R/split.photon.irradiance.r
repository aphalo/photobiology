#' Photon irradiance for split spectrum regions
#'
#' This function returns the photon irradiance for a series of contiguous
#' wavebands from a radiation spectrum. The returned values can be either
#' absolute or relative to their sum.
#'
#' @param w.length numeric vector of wavelengths (nm).
#' @param s.irrad numeric vector of spectral (energy or photon) irradiance values
#'   (W m-2 nm-1).
#' @param cut.w.length numeric vector of wavelengths (nm).
#' @param unit.in character Allowed values "energy", and "photon", or its alias
#'   "quantum".
#' @param scale a character A string indicating the scale used for the returned
#'   values ("absolute", "relative", "percent").
#' @param check.spectrum logical Flag indicating whether to sanity check input
#'   data, default is TRUE.
#' @param use.cached.mult logical Flag indicating whether multiplier values
#'   should be cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#'
#' @return a numeric vector of photon irradiances with no change in scale factor:
#'   [W m-2 nm-1] -> [mol s-1 m-2], [mol s-1 m-2 nm-1] -> [mol s-1 m-2] or
#'   relative values (fraction of one based on photon units) if scale
#'   = "relative" or scale = "percent".
#'
#' @export
#' @examples
#' with(sun.data,
#'      split_photon_irradiance(w.length, s.e.irrad,
#'              cut.w.length = c(300, 400, 500, 600, 700)))
#' with(sun.data,
#'      split_photon_irradiance(w.length, s.e.irrad))
#'
#' @note The last three parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you set \code{check.spectrum=FALSE}
#'   then you should call \code{\link{check_spectrum}} at least once for your
#'   spectrum before using any of the other functions. If you will use
#'   repeatedly the same SWFs on many spectra measured at exactly the same
#'   wavelengths you may obtain some speed up by setting
#'   \code{use.cached.mult=TRUE}. However, be aware that you are responsible for
#'   ensuring that the wavelengths are the same in each call, as the only test
#'   done is for the length of the \code{w.length} vector.
#'
#' @family low-level functions operating on numeric vectors.
#'
split_photon_irradiance <-
  function(w.length, s.irrad,
           cut.w.length = range(w.length),
           unit.in = "energy",
           scale = "absolute",
           check.spectrum = TRUE,
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges", default = NULL) )
  {
    return(split_irradiance(w.length, s.irrad,
                            cut.w.length = cut.w.length,
                            unit.out = "photon",
                            unit.in = unit.in,
                            scale = scale,
                            check.spectrum = check.spectrum,
                            use.cached.mult = use.cached.mult,
                            use.hinges = use.hinges))
  }
