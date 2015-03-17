
# options(datatable.print.topn=1)
# options(datatable.print.nrows=50)


# summary -----------------------------------------------------------------

#' Summary of a "generic.spct" object.
#'
#' A method of generic function summary for objects of class "generic.spct".
#'
#' @param object an object of class "generic.spct" for which a summary is desired
#' @param digits integer, used for number formatting with signif()
#' @param ... additional arguments affecting the summary produced, ignored in current version
#'
#' @export summary.generic.spct
#'
summary.generic.spct <- function(object, digits = 4L, ...) {
  z <- c(
    max.w.length = max(object),
    min.w.length = min(object),
    midpoint.w.length = midpoint(object),
    w.length.step = stepsize(object)[1],
  )
  z <- signif(z, digits)
  class(z) <- c("summary.generic.spct", class(z))
  return(z)
}

#' Print a "summary.generic.spct" object.
#'
#' A function to nicely print objects of class "summary.generic.spct".
#'
#' @param x an object of class "summary.generic.spct"
#' @param ... not used in current version
#'
#' @export print.summary.generic.spct
#'
print.summary.generic.spct <- function(x, ...) {
  time.unit <- attr(x, "time.unit")
  cat("wavelength ranges from", x[["min.w.length"]], "to", x[["max.w.length"]], "nm \n")
  cat("largest wavelength step size is", x[["w.length.step"]], "nm \n")
}

#' Summary of a "source.spct" object.
#'
#' A method of generic function summary for objects of class "source.spct".
#'
#' @param object an object of class "source.spct" for which a summary is desired
#' @param digits integer, used for number formatting with signif()
#' @param ... additional arguments affecting the summary produced, ignored in current version
#'
#' @export summary.source.spct
#'
#' @examples
#' str(summary(sun.spct))
summary.source.spct <- function(object, digits = 4L, ...) {
  time.unit <- getTimeUnit(object)
  bswf.used <- getBSWFUsed(object)
  z <- c(
  max.w.length = max(object),
  min.w.length = min(object),
  midpoint.w.length = midpoint(object),
  w.length.step = stepsize(object)[1],
  max.s.e.irrad = max(object$s.e.irrad),
  min.s.e.irrad = min(object$s.e.irrad),
  e.irrad = as.numeric(e_irrad(object)),
  q.irrad = as.numeric(q_irrad(object))
  )
  z <- signif(z, digits)
  attr(z, "time.unit") <- time.unit
  attr(z, "bswf.used") <- bswf.used
  class(z) <- c("summary.source.spct", class(z))
  return(z)
}

#' Print a "summary.source.spct" object.
#'
#' A function to nicely print objects of class "summary.source.spct".
#'
#' @param x an object of class "summary.source.spct"
#' @param ... not used in current version
#'
#' @export print.summary.source.spct
#'
#' @examples
#' summary(sun.spct)
#' summary(sun.daily.spct)


print.summary.source.spct <- function(x, ...) {
  time.unit <- attr(x, "time.unit")
  bswf.used <- attr(x, "bswf.used")
  cat("wavelength ranges from", x[["min.w.length"]], "to", x[["max.w.length"]], "nm \n")
  cat("largest wavelength step size is", x[["w.length.step"]], "nm \n")
  if (bswf.used != "none") {
    cat("effective irradiances based on BSWF =", bswf.used, "\n")
  }
  if (time.unit == "day") {
    cat("spectral irradiance ranges from", x[["min.s.e.irrad"]] * 1e-3, "to", x[["max.s.e.irrad"]] * 1e-3, "kJ d-1 m-2 nm-1 \n")
    cat("energy irradiance is", x[["e.irrad"]] * 1e-6, "MJ m-2 \n")
    cat("photon irradiance is", x[["q.irrad"]], "mol d-1 m-2 \n")
  } else if (time.unit == "second") {
    cat("spectral irradiance ranges from", x[["min.s.e.irrad"]], "to", x[["max.s.e.irrad"]], "W m-2 nm-1 \n")
    cat("energy irradiance is", x[["e.irrad"]], "W m-2 \n")
    cat("photon irradiance is", x[["q.irrad"]] * 1e6, "umol s-1 m-2\n")
  } else {
    cat("spectral irradiance ranges from", x[["min.s.e.irrad"]], "to", x[["max.s.e.irrad"]], "\n")
    cat("energy irradiance is", x[["e.irrad"]], "\n")
    cat("photon irradiance is", x[["q.irrad"]], "\n")
  }
}

#' Summary of a "filter.spct" object.
#'
#' A method of generic function summary for objects of class "filter.spct".
#'
#' @param object an object of class "filter.spct" for which a summary is desired
#' @param digits integer, used for number formatting with signif()
#' @param ... additional arguments affecting the summary produced, ignored in current version
#'
#' @export summary.filter.spct
#'
summary.filter.spct <- function(object, digits = 4L, ...) {
  Tfr.type <- getTfrType(object)
  z <- c(
    max.w.length = max(object),
    min.w.length = min(object),
    midpoint.w.length = midpoint(object),
    w.length.step = stepsize(object)[1],
    max.Tfr = max(object$Tfr),
    min.Tfr = min(object$Tfr),
    mean.Tfr = as.numeric(integrate_spct(object) / spread(object))
  )
  z <- signif(z, digits)
  attr(z, "Tfr.type") <- Tfr.type
  class(z) <- c("summary.filter.spct", class(z))
  return(z)
}

#' Print a "summary.filter.spct" object.
#'
#' A function to nicely print objects of class "summary.filter.spct".
#'
#' @param x an object of class "summary.filter.spct"
#' @param ... not used in current version
#'
#' @export print.summary.filter.spct
#'
print.summary.filter.spct <- function(x, ...) {
  Tfr.type <- attr(x, "Tfr.type")
  cat("wavelength ranges from", x[["min.w.length"]], "to", x[["max.w.length"]], "nm \n")
  cat("largest wavelength step size is", x[["w.length.step"]], "nm \n")
  cat("Spectral transmittance ranges from", x[["min.Tfr"]], "to", x[["max.Tfr"]], "\n")
  cat("Mean transmittance is", x[["mean.Tfr"]], "\n")
  cat("Quantity is", Tfr.type, "\n")
}

#' Summary of a "reflector.spct" object.
#'
#' A method of generic function summary for objects of class "reflector.spct".
#'
#' @param object an object of class "reflector.spct" for which a summary is desired
#' @param digits integer, used for number formatting with signif()
#' @param ... additional arguments affecting the summary produced, ignored in current version
#'
#' @export summary.reflector.spct
#'
summary.reflector.spct <- function(object, digits = 4L, ...) {
  Rfr.type <- getRfrType(object)
  z <- c(
    max.w.length = max(object),
    min.w.length = min(object),
    midpoint.w.length = midpoint(object),
    w.length.step = stepsize(object)[1],
    max.Rfr = max(object$Rfr),
    min.Rfr = min(object$Rfr),
    mean.Rfr = as.numeric(integrate_spct(object) / spread(object))
  )
  z <- signif(z, digits)
  attr(z, "Rfr.type") <- Rfr.type
  class(z) <- c("summary.reflector.spct", class(z))
  return(z)
}

#' Print a "summary.reflector.spct" object.
#'
#' A function to nicely print objects of class "summary.reflector.spct".
#'
#' @param x an object of class "summary.reflector.spct"
#' @param ... not used in current version
#'
#' @export print.summary.reflector.spct
#'
print.summary.reflector.spct <- function(x, ...) {
  Rfr.type <- attr(x, "Rfr.type")
  cat("wavelength ranges from", x[["min.w.length"]], "to", x[["max.w.length"]], "nm \n")
  cat("largest wavelength step size is", x[["w.length.step"]], "nm \n")
  cat("Spectral reflectance ranges from", x[["min.Rfr"]], "to", x[["max.Rfr"]], "\n")
  cat("Mean reflectance is", x[["mean.Rfr"]], "\n")
  cat("Quantity is", Rfr.type, "\n")
}

#' Summary of a "response.spct" object.
#'
#' A method of generic function summary for objects of class "response.spct".
#'
#' @param object an object of class "response.spct" for which a summary is desired
#' @param digits integer, used for number formatting with signif()
#' @param ... additional arguments affecting the summary produced, ignored in current version
#'
#' @export summary.response.spct
#'
summary.response.spct <- function(object, digits = 4L, ...) {
  time.unit <- getTimeUnit(object)
  z <- c(
    max.w.length = max(object),
    min.w.length = min(object),
    midpoint.w.length = midpoint(object),
    w.length.step = stepsize(object)[1],
    max.response = max(object$s.e.response),
    min.response = min(object$s.e.response),
    total.response = as.numeric(integrate_spct(object)),
    mean.response = as.numeric(integrate_spct(object) / spread(object))
  )
  z <- signif(z, digits)
  class(z) <- c("summary.response.spct", class(z))
  attr(z, "time.unit") <- time.unit
  return(z)
}

#' Print a "summary.response.spct" object.
#'
#' A function to nicely print objects of class "summary.response.spct".
#'
#' @param x an object of class "summary.response.spct"
#' @param ... not used in current version
#'
#' @export print.summary.response.spct
#'
print.summary.response.spct <- function(x, ...) {
  time.unit <- attr(x, "time.unit")
  cat("wavelength ranges from", x[["min.w.length"]], "to", x[["max.w.length"]], "nm \n")
  cat("largest wavelength step size is", x[["w.length.step"]], "nm \n")
  cat("Spectral response ranges from", x[["min.response"]], "to", x[["max.response"]], "nm-1 \n")
  cat("Mean response is", x[["mean.response"]], "nm-1 \n")
  cat("Time unit is", time.unit, "\n")
}

#' Summary of a "chroma.spct" object.
#'
#' A method of generic function summary for objects of class "chroma.spct".
#'
#' @param object an object of class "chroma.spct" for which a summary is desired
#' @param digits integer, used for number formatting with signif()
#' @param ... additional arguments affecting the summary produced, ignored in current version
#'
#' @export summary.chroma.spct
#'
summary.chroma.spct <- function(object, digits = 4L, ...) {
  z <- c(
    max.w.length = max(object),
    min.w.length = min(object),
    midpoint.w.length = midpoint(object),
    w.length.step = stepsize(object)[1],
    x.max = max(object[["x"]]),
    y.max = max(object[["y"]]),
    z.max = max(object[["z"]])
  )
  z <- signif(z, digits)
  class(z) <- c("summary.chroma.spct", class(z))
  return(z)
}

#' Print a "summary.chroma.spct" object.
#'
#' A function to nicely print objects of class "summary.chroma.spct".
#'
#' @param x an object of class "summary.chroma.spct"
#' @param ... not used in current version
#'
#' @export print.summary.chroma.spct
#'
print.summary.chroma.spct <- function(x, ...) {
  time.unit <- attr(x, "time.unit")
  cat("wavelength ranges from", x[["min.w.length"]], "to", x[["max.w.length"]], "nm \n")
  cat("largest wavelength step size is", x[["w.length.step"]], "nm \n")
  cat(paste("maximum (x, y, z) values are (", paste(x[["x.max"]], x[["y.max"]], x[["z.max"]], sep=", "), ")", sep=""), "\n")
}

#' Color of a source.spct object.
#'
#' A function that returns the equivalent RGB colour of an object of class "source.spct".
#'
#' @param x an object of class "source.spct"
#' @param ... not used in current version
#' @export color.source.spct
#'
color.source.spct <- function(x, ...) {
#  x.name <- as.character(substitute(x))
  x.name <- "source"
  q2e(x, byref=TRUE)
  color <- c(s_e_irrad2rgb(x[["w.length"]], x[["s.e.irrad"]], sens=ciexyzCMF2.spct, color.name=paste(x.name, "CMF")),
             s_e_irrad2rgb(x[["w.length"]], x[["s.e.irrad"]], sens=ciexyzCC2.spct, color.name=paste(x.name, "CC")))
  return(color)
}

# w.length summaries ------------------------------------------------------


#' "range" function for spectra
#'
#' Range function for spectra, returning wavelength range.
#'
#' @param ... not used in current version
#' @param na.rm a logical indicating whether missing values should be removed.
#' @export
#'
range.generic.spct <- function(..., na.rm=FALSE) {
  x <- c(...)
  return(range(x[["w.length"]], na.rm=na.rm))
}

#' "max" function for spectra
#'
#' Maximun function for spectra, returning wavelength maximum.
#'
#' @param ... not used in current version
#' @param na.rm a logical indicating whether missing values should be removed.
#' @export
#'
max.generic.spct <- function(..., na.rm=FALSE) {
  x <- c(...)
  return(max(x[["w.length"]], na.rm=na.rm))
}

#' "min" function for spectra
#'
#' Minimun function for spectra, returning wavelength minimum.
#'
#' @param ... not used in current version
#' @param na.rm a logical indicating whether missing values should be removed.
#' @export
#'
min.generic.spct <- function(..., na.rm=FALSE) {
  x <- c(...)
  return(min(x[["w.length"]], na.rm=na.rm))
}

#' Generic function
#'
#' Function that returns the range of step sizes in an object.
#'
#' @param x an R object
#' @param ... not used in current version
#' @export stepsize
stepsize <- function(x, ...) UseMethod("stepsize")

#' Default for generic function
#'
#' Function that returns the range of step sizes in an object.
#'
#' @param x an R object
#' @param ... not used in current version
#' @export stepsize.default
stepsize.default <- function(x, ...) {
  return(range(diff(x)))
}

#' Method for "generic.spct" objects for generic function
#'
#' Function that returns the range of wavelength step sizes in a "generic.spct" object.
#'
#' @param x an R object
#' @param ... not used in current version
#' @export stepsize.generic.spct
#'
#' @examples
#' stepsize(sun.spct)
#'
stepsize.generic.spct <- function(x, ...) {
  range(diff(x[["w.length"]]))
}

#' Labels of a "generic.spct" object.
#'
#' A function to obtain the labels of a spectrum. Currently returns 'names'.
#'
#' @param object an object of generic.spct
#' @param ... not used in current version
#'
#' @export
#'
labels.generic.spct <- function(object, ...) {
  return(names(object))
}
