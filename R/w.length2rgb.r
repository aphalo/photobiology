#' Wavelength to rgb color conversion
#'
#' Calculates rgb values from spectra based on human color matching functions
#'
#' @usage w_length2rgb(w.length, sens=ciexyzCMF2.spct, color.name=NULL)
#'
#' @param w.length numeric array of wavelengths (nm)
#' @param sens a chroma.spct object with variables w.length, x, y, and z, giving the chromaticity definition
#' @param color.name character string for naming the rgb color definition
#'
#' @return An atrray of colors defined using \code{rgb()}. The numeric values of the RGB components can be obtained
#'   using function \code{col2rgb()}.
#'
#' @export
#' @examples
#' col2rgb(w_length2rgb(580))
#' col2rgb(w_length2rgb(c(400, 500, 600, 700)))
#' col2rgb(w_length2rgb(c(400, 500, 600, 700), color.name=c("a","b","c","d")))
#' col2rgb(w_length2rgb(c(400, 500, 600, 700), color.name=c("a","b")))
#' col2rgb(w_length2rgb(c(400, 500, 600, 700), color.name="a"))
#'
#' @author Pedro J. Aphalo

w_length2rgb <- function(w.length, sens=ciexyzCMF2.spct, color.name=NULL) {
  len.wl <- length(w.length)
  generate.names <- is.null(color.name)
  if (!generate.names) {
    len.col <- length(color.name)
    if (len.col == 1L) {
      color.names <- rep(color.name[1], length.out = len.wl)
    } else if (len.col < len.wl) {
      warning("color.name argument shorter than w.length argument.")
      color.names <- color.name
    } else {
      color.names <- color.name
    }
  } else {
    color.names <-NULL
  }
  colors <- NULL
  for (i in 1:len.wl) {
    colors[i] <-  s_e_irrad2rgb(w.length[i], 1.0, sens=sens)
    if (generate.names) color.names[i] <- paste(as.character(w.length[i]), "nm")
  }
  names(colors) <- color.names
  return(colors)
}
