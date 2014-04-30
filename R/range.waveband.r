#' Wavelength range of a "waveband" object.
#' 
#' A function that returns the wavelength range from objects of class "waveband".
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
#' A function that returns the wavelength minimum from objects of class "waveband".
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
#' A function that returns the wavelength maximum from objects of class "waveband".
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
#' A function that returns the wavelength at the middle of the wavelength range of
#' objects of class "waveband".
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

#' Generic function
#' 
#' A function that returns the wavelength at the middle of the wavelength range.
#' 
#' @param x an R object
#' @export center_wl center_wl.default
center_wl <- function(x) UseMethod("center_wl", x)

#' Default for generic function
#' 
#' A function that returns the wavelength at the middle of the wavelength range.
#' 
#' @param x an R object
#' @export center_wl center_wl.default
center_wl.default <- function(x) {
  return(NA)
}

#' Wavelength at center of a "waveband" object.
#' 
#' A function that returns the wavelength at the middle of the wavelength range of
#' objects of class "waveband".
#' 
#' @param x an object of class "waveband"
#' @export center_wl.waveband
#' 
center_wl.waveband <- function(x) {
  return((x$low + x$high) / 2)  
}

#' Generic function that returns the color of an object.
#' 
#' A function that returns the equivalent RGB color of an object.
#' 
#' @param x an R object
#' @export color color.default
#' 
color <- function(x) UseMethod("color", x)

#' Default of function that returns Color of an object.
#' 
#' A function that returns the equivalent RGB color of an object.
#' 
#' @param x an R object
#' @export color color.default
#' 
color.default <- function(x) {
  return("black")
}

#' Color at center of a "waveband" object.
#' 
#' A function that returns the equivalent RGB colour of an object of class "waveband".
#' 
#' @param x an object of class "waveband"
#' @export color.waveband
#' 
color.waveband <- function(x) {
  color <- c(w_length_range2rgb(range(x), sens=ciexyzCMF2.data, color.name=paste(labels(x)[1], "CMF")), 
             w_length_range2rgb(range(x), sens=ciexyzCC2.data, color.name=paste(labels(x)[1], "CC")))
  return(color)           
}
