#' Calculate RGB values from spectral irradiance.
#'
#' This function returns the RGB values for a source spectrum.
#'
#' @usage rgb_spct(spct, sens=ciexyzCMF2.data, color.name=NULL)
#'
#' @param spct an object of class "source.spct"
#' @param sens a dataframe with variables w.length, x, y, and z, giving the CC or CMF definition (default is the
#' proposed human CMF according to CIE 2006.)
#' @param color.name character string for naming the rgb color definition
#'
#' @return A color defined using \code{rgb()}. The numeric values of the RGB components can be obtained
#' @keywords manip misc
#' @export
#' @examples
#' rgb_spct(sun.spct)
#'

rgb_spct <-
  function(spct, sens=ciexyzCMF2.data, color.name=NULL){
    if (is(spct, "source.spct")) {
      return(s_e_irrad2rgb(spct$w.length, spct$s.e.irrad, sens=sens, color.name=color.name))
    } else {
      return(NA)
    }
  }
