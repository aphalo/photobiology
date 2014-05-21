# range -------------------------------------------------------------------

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

# min ---------------------------------------------------------------------

#' Wavelength minimum of a "waveband" object.
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

# max ---------------------------------------------------------------------

#' Wavelength maximum of a "waveband" object.
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

# midpoint ------------------------------------------------------------------

#' Generic function
#' 
#' A function that returns the wavelength at the center of the wavelength range.
#' 
#' @param x an R object
#' @export midpoint
midpoint <- function(x) UseMethod("midpoint", x)

#' Default for generic function
#' 
#' A function that returns the wavelength at the center of the wavelength range.
#' 
#' @param x an R object
#' @export midpoint.default
midpoint.default <- function(x) {
  return(min(x) + (max(x) - min(x)) / 2)
}

#' Wavelength at center of a "waveband" object.
#' 
#' A function that returns the wavelength at the center of the wavelength range of
#' objects of class "waveband".
#' 
#' @param x an object of class "waveband"
#' @export midpoint.waveband
#' 
midpoint.waveband <- function(x) {
  return(x$low + (x$high - x$low) / 2)
}

# spread ------------------------------------------------------------------

#' Generic function
#' 
#' A function that returns the spread (max(x) - min(x)) for R objects.
#' 
#' @param x an R object
#' @export spread
spread <- function(x) UseMethod("spread", x)

#' Default for generic function
#' 
#' A function that returns the spread (max(x) - min(x)) for objects.
#' 
#' @param x an R object
#' @export spread.default
#' 
spread.default <- function(x) {
  return(max(x) - min(x))
}

#' Wavelength spread (max-min) of a "waveband" object.
#' 
#' A function that returns the wavelength spread from objects of class "waveband".
#' 
#' @param x an object of class "waveband"
#' @export
#' 
spread.waveband <- function(x) {
  return(x$high - x$low)
}

# color -------------------------------------------------------------------

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
