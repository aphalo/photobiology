#' Wavelength to rgb color conversion
#'
#' Calculates rgb values from spectra based on human color matching functions
#'
#' @usage w_length2rgb(w.length, sens=ciexyz.data, colour.name=NULL)
#'
#' @param w.length numeric array of wavelengths (nm)
#' @param sens a dataframe with variables w.length, x, y, and z, giving the chromaticity definition 
#' @param colour.name character string for naming the rgb colour definition
#' 
#' @return A colour defined using \code{rgb()}. The numeric values of the RGB components can be obtained 
#' using function \code{col2rgb()}.
#' 
#' @export
#' @examples
#' col2rgb(s_e_irrad2rgb(580))
#' 
#' @author Pedro J. Aphalo 

w_length2rgb <- function(w.length, sens=ciexyz.data, colour.name=NULL) {
  len <- length(w.length)
  colors <- NULL
  for (i in 1:len) {
    colors[i] <-  s_e_irrad2rgb(w.length[i], 1.0, sens=sens, 
                                colour.name=ifelse(is.null(colour.name), paste(as.character(w.length[i], "nm")), colour.name))
  }
  return(colors)
}
