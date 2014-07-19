#' Calculate energy irradiance from spectral irradiance.
#'
#' This function returns the energy irradiance for a given
#' waveband of a light source spectrum.
#'
#' @usage e_irrad_spct(spct, w.band=NULL,
#' use.cached.mult=FALSE, use.hinges=NULL)
#'
#' @param spct an object of class "source.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @return a single numeric value with no change in scale factor: [W m-2 nm-1] -> [W m-2]
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.spct)
#' e_irrad_spct(sun.spct, new_waveband(400,700))
#'
#' @note The last two parameters control speed optimizations. The defaults should be suitable
#' in mosts cases. If you will use repeatedly
#' the same SWFs on many spectra measured at exactly the same wavelengths you may obtain some speed up
#' by setting \code{use.cached.mult=TRUE}. However, be aware that you are responsible for ensuring
#' that the wavelengths are the same in each call, as the only test done is for the length of the
#' \code{w.length} vector.

e_irrad_spct <-
  function(spct, w.band=NULL, use.cached.mult=FALSE, use.hinges=NULL){
     return(irrad_spct(spct, w.band=w.band, unit.out="energy",
                    use.cached.mult=use.cached.mult, use.hinges=use.hinges))
  }
