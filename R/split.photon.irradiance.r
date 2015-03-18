#' Split a spectrum into contiguous bands and calculate photon irradiance
#' from spectral (energy) or photon irradiance.
#'
#' This function returns the photon irradiance for a series of contiguous wavebands
#' from a radiation spectrum. The returned values can be either absolute or relative to their
#' sum.
#'
#' @usage split_photon_irradiance(w.length, s.irrad, cut.w.length=range(w.length),
#'                                unit.in="energy", scale = "absolute",
#'                                check.spectrum=TRUE, use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
#'                use.hinges=getOption("photobiology.use.hinges", default=NULL) )
#'
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral (energy) irradiances (W m-2 nm-1)
#' @param cut.w.length numeric array of wavelengths (nm)
#' @param unit.in character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param scale a character string indicating the scale used for the returned values ("absolute", "relative", "percent")
#' @param check.spectrum logical indicating whether to sanity check input data, default is TRUE
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @return a numeric array of photon irradiances with no change in scale factor: [W m-2 nm-1] -> [mol s-1 m-2]
#' or relative values (fraction of one) if scale = "relative" or scale = "percent"
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' with(sun.data,
#'      split_photon_irradiance(w.length, s.e.irrad,
#'              cut.w.length=c(300, 400, 500, 600, 700)))
#' with(sun.data,
#'      split_photon_irradiance(w.length, s.e.irrad,
#'              cut.w.length=c(200, 400, 500, 600, 900)))
#' with(sun.data,
#'      split_photon_irradiance(w.length, s.e.irrad,
#'              cut.w.length=c(300, 400, 500)))
#' with(sun.data,
#'      split_photon_irradiance(w.length, s.e.irrad,
#'              cut.w.length=c(300, 400)))
#' with(sun.data,
#'      split_photon_irradiance(w.length, s.e.irrad,
#'              cut.w.length=c(300, 400, 300)))
#' with(sun.data,
#'      split_photon_irradiance(w.length, s.e.irrad,
#'              cut.w.length=c(100, 200)))
#' with(sun.data,
#'      split_photon_irradiance(w.length, s.e.irrad,
#'              cut.w.length=c(1000, 1200)))
#' with(sun.data,
#'      split_photon_irradiance(w.length, s.e.irrad,
#'              cut.w.length=c(300)))
#' with(sun.data,
#'      split_photon_irradiance(w.length, s.e.irrad))
#'
#' @note The last three parameters control speed optimizations. The defaults should be suitable
#' in mosts cases. If you set \code{check.spectrum=FALSE} then you should call \code{check_spectrum()}
#' at least once for your spectrum before using any of the other functions. If you will use repeatedly
#' the same SWFs on many spectra measured at exactly the same wavelengths you may obtain some speed up
#' by setting \code{use.cached.mult=TRUE}. However, be aware that you are responsible for ensuring
#' that the wavelengths are the same in each call, as the only test done is for the length of the
#' \code{w.length} vector.

split_photon_irradiance <- function(w.length, s.irrad, cut.w.length=range(w.length), unit.in="energy", scale="absolute",
                 check.spectrum=TRUE, use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
                 use.hinges=getOption("photobiology.use.hinges", default=NULL) )
{
    return(split_irradiance(w.length, s.irrad, cut.w.length=cut.w.length, unit.out="photon", unit.in=unit.in, scale=scale,
                            check.spectrum=check.spectrum, use.cached.mult=use.cached.mult, use.hinges=use.hinges))
}
