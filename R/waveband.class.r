# labels ------------------------------------------------------------------

#' Name and label of a "waveband" object.
#'
#' A function to obtain the name and label of objects of class "waveband".
#'
#' @param object an object of class "waveband"
#' @param ... not used in current version
#'
#' @export
#'
labels.waveband <- function(object, ...) {
  return(list(label = object$label, name = object$name))
}

# range -------------------------------------------------------------------

#' Wavelength range of a "waveband" object.
#'
#' A function that returns the wavelength range from objects of class "waveband".
#'
#' @param ... not used in current version
#' @param na.rm ignored
#' @export
#'
range.waveband <- function(..., na.rm = FALSE) {
  x <- c(...)
  return(c(x$low, x$high))
}

# min ---------------------------------------------------------------------

#' Wavelength minimum of a "waveband" object.
#'
#' A function that returns the wavelength minimum from objects of class "waveband".
#'
#' @param ... not used in current version
#' @param na.rm ignored
#' @export
#'
min.waveband <- function(..., na.rm = FALSE) {
  x <- c(...)
    return(x$low)
}

# max ---------------------------------------------------------------------

#' Wavelength maximum of a "waveband" object.
#'
#' A function that returns the wavelength maximum from objects of class "waveband".
#'
#' @param ... not used in current version
#' @param na.rm ignored
#' @export
#'
max.waveband <- function(..., na.rm = FALSE) {
  x <- c(...)
  return(x$high)
}

# midpoint ------------------------------------------------------------------

#' Generic function
#'
#' A function that returns the wavelength at the center of the wavelength range.
#'
#' @param x an R object
#' @export midpoint
midpoint <- function(x) UseMethod("midpoint")

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
spread <- function(x) UseMethod("spread")

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
#' @param ... not used in current version
#' @export color
#'
color <- function(x, ...) UseMethod("color")

#' Default of function that returns Color of an object.
#'
#' A function that returns the equivalent RGB color of an object.
#'
#' @param x an R object
#' @param ... not used in current version
#' @export color.default
#'
color.default <- function(x, ...) {
  return("#000000")
}

#' Color at center of a "waveband" object.
#'
#' A function that returns the equivalent RGB colour of an object of class "waveband".
#'
#' @param x an object of class "waveband"
#' @param short.names logical indicating whether to use short or long names for wavebands
#' @param ... not used in current version
#' @export color.waveband
#'
color.waveband <- function(x, short.names=TRUE, ...) {
  idx <- ifelse(short.names, "name", "label")
  color <- c(w_length_range2rgb(range(x), sens=ciexyzCMF2.spct, color.name=paste(labels(x)[[idx]], "CMF")),
             w_length_range2rgb(range(x), sens=ciexyzCC2.spct, color.name=paste(labels(x)[[idx]], "CC")))
  return(color)
}


# normalization -----------------------------------------------------------

#' Normalization of an R object.
#'
#' A generic function that returns the normalization of an R object.
#'
#' @param x an R object
#' @export normalization.default
#'
normalization <- function(x) UseMethod("normalization")

#' Normalization of an R object.
#'
#' A generic function that returns the normalization of an R object.
#'
#' @param x an R object
#' @export normalization
#'
normalization.default <- function(x) {
  return(NA)
}

#' Normalization of a "waveband" object.
#'
#' A function that returns the normalization wavelength of a waveband object.
#'
#' @param x an object of class "waveband"
#' @export normalization.waveband
#'
normalization.waveband <- function(x) {
  return(ifelse(is.null(x$norm), NA, x$norm))
}

# is_effective -----------------------------------------------------------

#' Is an R object "effective".
#'
#' A generic function that returns the normalization of an R object.
#'
#' @param x an R object
#' @export is_effective.default
#'
is_effective <- function(x) UseMethod("is_effective")

#' Is an R object "effective".
#'
#' A generic function that returns the normalization of an R object.
#'
#' @param x an R object
#' @export is_effective
#'
is_effective.default <- function(x) {
  return(NA)
}

#' Is a "waveband" object defining an effective irradiance.
#'
#' A function that returns the normalization wavelength of a waveband object.
#'
#' @param x an object of class "waveband"
#' @export is_effective.waveband
#'
is_effective.waveband <- function(x) {
  return(x$weight != "none")
}