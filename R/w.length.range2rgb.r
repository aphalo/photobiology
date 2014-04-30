#' Wavelength range to rgb color conversion
#'
#' Calculates rgb values from spectra based on human color matching functions
#'
#' @usage w_length_range2rgb(w.length, sens=ciexyzCMF2.data, color.name=NULL)
#'
#' @param w.length numeric array of wavelengths (nm) of length 2. If longer, its range is used.
#' @param sens a dataframe with variables w.length, x, y, and z, giving the chromaticity definition 
#' @param color.name character string for naming the rgb color definition
#' 
#' @return An atrray of colors defined using \code{rgb()}. The numeric values of the RGB components can be obtained 
#' using function \code{col2rgb()}.
#' 
#' @export
#' @examples
#' col2rgb(w_length_range2rgb(c(500,600)))
#' col2rgb(w_length_range2rgb(550))
#' col2rgb(w_length_range2rgb(500:600))
#' @author Pedro J. Aphalo 

w_length_range2rgb <- function(w.length, sens=ciexyzCMF2.data, color.name=NULL) {
  w.length <- unique(sort(w.length))
  len <- length(w.length)
  if (len < 2 ) {
    warning("Calculating RGB values for monochromatic light.")
    return(w_length2rgb(w.length, sens, color.name))
  } else if (len > 2) {
    warning("Using only extreme wavelength values.")
    w.length <- range(w.length)
  }
  num.values <- max(50L, trunc(w.length[2] - w.length[1]))
  w.length.values <- seq(w.length[1], w.length[2], length.out = num.values)
  s.e.irrad.values <- rep(1.0, length.out = num.values)
  color <-  s_e_irrad2rgb(w.length.values, s.e.irrad.values, sens=sens, 
                                color.name=ifelse(is.null(color.name), 
                              paste(as.character(w.length[1]), "-", as.character(w.length[2]), " nm", sep=""), 
                           color.name))
  return(color)
}
