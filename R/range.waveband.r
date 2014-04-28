#' Wavelength range of a "waveband" object.
#' 
#' A function to extract the wavelength range from objects of class "waveband".
#' 
#' @param x an object of class "waveband"
#' @param ... not used in current version
#' @param na.rm ignored
#' @export
#' 
range.waveband <- function(x, ..., na.rm = FALSE) {
  return(c(x$low, x$high))
}

#' Wavelength range of a "waveband" object.
#' 
#' A function to extract the wavelength range from objects of class "waveband".
#' 
#' @param x an object of class "waveband"
#' @param ... not used in current version
#' @param na.rm ignored
#' @export
#' 
min.waveband <- function(x, ..., na.rm = FALSE) {
  return(x$low)
}

#' Wavelength range of a "waveband" object.
#' 
#' A function to extract the wavelength range from objects of class "waveband".
#' 
#' @param x an object of class "waveband"
#' @param ... not used in current version
#' @param na.rm ignored
#' @export
#' 
max.waveband <- function(x, ..., na.rm = FALSE) {
  return(x$high)
}

#' Wavelength at center of a "waveband" object.
#' 
#' A function to extract the wavelength range from objects of class "waveband".
#' 
#' @param x an object of class "waveband"
#' @param trim ignored
#' @param ... not used in current version
#' @param na.rm ignored
#' @export
#' 
mean.waveband <- function(x, trim = 0, na.rm = FALSE, ...) {
  return((x$low + x$high) / 2)  
}
