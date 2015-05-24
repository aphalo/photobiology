# summary -----------------------------------------------------------------

#' Summary of a spectral object
#'
#' Methods of generic function summary for objects of spectral classes.
#'
#' @param object An object of one of the spectral classes for which a summary is
#'   desired
#' @param digits integer Used for number formatting with \code{signif()}
#' @param ... additional arguments affecting the summary produced, ignored in
#'   current version
#'
#' @method summary generic_spct
#'
#' @return A summary object matching the class of \code{object}.
#'
#' @export summary.generic_spct
#'
summary.generic_spct <- function(object, digits = max(3, getOption("digits")-3), ...) {
  z <- c(
    max.w.length = max(object),
    min.w.length = min(object),
    midpoint.w.length = midpoint(object),
    w.length.step = stepsize(object)[1]
  )
  z <- signif(z, digits)
  class(z) <- c("summary.generic_spct", class(z))
  return(z)
}

# @describeIn summary.generic_spct Summary of a "source_spct" object.
#'
#' @param time.unit character or lubridate::duration
#'
#' @export
#' @rdname summary.generic_spct
#'
summary.source_spct <- function(object,
                                digits = max(3, getOption("digits")-3),
                                time.unit = NULL,
                                ...) {
  if (!is.null(time.unit)) {
    object <- convertTimeUnit(object, time.unit = time.unit, byref = FALSE)
  } else {
    time.unit <- getTimeUnit(object)
  }
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
  class(z) <- c("summary.source_spct", class(z))
  return(z)
}

# @describeIn summary.generic_spct Summary of a \code{filter_spct} object.
#'
#' @export
#' @rdname summary.generic_spct
#'
summary.filter_spct <- function(object, digits = max(3, getOption("digits")-3), ...) {
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
  class(z) <- c("summary.filter_spct", class(z))
  return(z)
}

# @describeIn summary.generic_spct Summary of a "reflector_spct" object.
#'
#' @export
#' @rdname summary.generic_spct
#'
summary.reflector_spct <- function(object, digits = max(3, getOption("digits")-3), ...) {
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
  class(z) <- c("summary.reflector_spct", class(z))
  return(z)
}

# @describeIn summary.generic_spct Summary of a "response_spct" object.
#'
#' @export
#' @rdname summary.generic_spct
#'
summary.response_spct <- function(object,
                                  digits = max(3, getOption("digits")-3),
                                  time.unit = NULL,
                                  ...) {
  if (!is.null(time.unit)) {
    object <- convertTimeUnit(object, time.unit = time.unit, byref = FALSE)
  } else {
    time.unit <- getTimeUnit(object)
  }
  z <- c(
    max.w.length = max(object),
    min.w.length = min(object),
    midpoint.w.length = midpoint(object),
    w.length.step = stepsize(object)[1],
    max.response = max(object$s.e.response),
    min.response = min(object$s.e.response),
    total.q_response = as.numeric(q_response(object, quantity = "total")),
    mean.q_response = as.numeric(q_response(object, quantity = "mean")),
    total.e_response = as.numeric(e_response(object, quantity = "total")),
    mean.e_response = as.numeric(e_response(object, quantity = "mean"))
  )
  z <- signif(z, digits)
  class(z) <- c("summary.response_spct", class(z))
  attr(z, "time.unit") <- time.unit
  return(z)
}

# @describeIn summary.generic_spct Summary of a "chroma_spct" object.
#'
#' @export
#' @rdname summary.generic_spct
#'
summary.chroma_spct <- function(object, digits = max(3, getOption("digits")-3), ...) {
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
  class(z) <- c("summary.chroma_spct", class(z))
  return(z)
}


# Print spectral summaries ------------------------------------------------

#' Print a summary object of a spectrum.
#'
#' A function to nicely print objects of classes "summary...spct".
#'
#' @param x An object of one of the summary classes for spectra
#' @param ... not used in current version
#'
#' @method print summary.generic_spct
#'
#' @export print.summary.generic_spct
#'
print.summary.generic_spct <- function(x, ...) {
  time.unit <- attr(x, "time.unit")
  cat("wavelength ranges from", x[["min.w.length"]], "to", x[["max.w.length"]], "nm \n")
  cat("largest wavelength step size is", x[["w.length.step"]], "nm \n")
}

# @describeIn print.summary.generic_spct Print a "summary.source_spct" object.
#'
#' @export
#' @rdname print.summary.generic_spct
#'
print.summary.source_spct <- function(x, ...) {
  time.unit <- attr(x, "time.unit")
  bswf.used <- attr(x, "bswf.used")
  cat("wavelength ranges from", x[["min.w.length"]], "to", x[["max.w.length"]], "nm \n")
  cat("largest wavelength step size is", x[["w.length.step"]], "nm \n")
  if (bswf.used != "none") {
    cat("effective irradiances based on BSWF =", bswf.used, "\n")
  }
  if (time.unit == "day" || time.unit == lubridate::duration(1, "days")) {
    cat("spectral irradiance ranges from", x[["min.s.e.irrad"]] * 1e-3,
        "to", x[["max.s.e.irrad"]] * 1e-3, "kJ d-1 m-2 nm-1 \n")
    cat("energy irradiance is", x[["e.irrad"]] * 1e-6, "MJ d-1 m-2 \n")
    cat("photon irradiance is", x[["q.irrad"]], "mol d-1 m-2 \n")
  } else if (time.unit == "second" || time.unit == lubridate::duration(1, "seconds")) {
    cat("spectral irradiance ranges from", x[["min.s.e.irrad"]],
        "to", x[["max.s.e.irrad"]], "W m-2 nm-1 \n")
    cat("energy irradiance is", x[["e.irrad"]], "W m-2 \n")
    cat("photon irradiance is", x[["q.irrad"]] * 1e6, "umol s-1 m-2\n")
  } else {
    cat("spectral irradiance ranges from", x[["min.s.e.irrad"]],
        "to", x[["max.s.e.irrad"]], "J m-2 nm-1 per", as.character(time.unit), "\n")
    cat("energy irradiance is", x[["e.irrad"]], "J m-2 per", as.character(time.unit), "\n")
    cat("photon irradiance is", x[["q.irrad"]], "mol m-2 per", as.character(time.unit), "\n")
  }
}

# @describeIn print.summary.generic_spct Print a "summary.filter_spct" object.
#'
#' @export
#' @rdname print.summary.generic_spct
#'
print.summary.filter_spct <- function(x, ...) {
  Tfr.type <- attr(x, "Tfr.type")
  cat("wavelength ranges from", x[["min.w.length"]], "to", x[["max.w.length"]], "nm \n")
  cat("largest wavelength step size is", x[["w.length.step"]], "nm \n")
  cat("Spectral transmittance ranges from", x[["min.Tfr"]], "to", x[["max.Tfr"]], "\n")
  cat("Mean transmittance is", x[["mean.Tfr"]], "\n")
  cat("Quantity is", Tfr.type, "\n")
}

# @describeIn print.summary.generic_spct Print a "summary.reflector_spct" object.
#'
#' @export
#' @rdname print.summary.generic_spct
#'
print.summary.reflector_spct <- function(x, ...) {
  Rfr.type <- attr(x, "Rfr.type")
  cat("wavelength ranges from", x[["min.w.length"]], "to", x[["max.w.length"]], "nm \n")
  cat("largest wavelength step size is", x[["w.length.step"]], "nm \n")
  cat("Spectral reflectance ranges from", x[["min.Rfr"]], "to", x[["max.Rfr"]], "\n")
  cat("Mean reflectance is", x[["mean.Rfr"]], "\n")
  cat("Quantity is", Rfr.type, "\n")
}

# @describeIn print.summary.generic_spct Print a "summary.response_spct" object.
#'
#' @export
#' @rdname print.summary.generic_spct
#'
print.summary.response_spct <- function(x, ...) {
  time.unit <- attr(x, "time.unit")
  cat("wavelength ranges from", x[["min.w.length"]], "to", x[["max.w.length"]], "nm \n")
  cat("largest wavelength step size is", x[["w.length.step"]], "nm \n")
  cat("Spectral response ranges from", x[["min.response"]],
      "to", x[["max.response"]], "response units W-1 nm-1 per", as.character(time.unit), "\n")
  cat("Mean response is", x[["mean.e_response"]], "response units W-1 nm-1 per",
      as.character(time.unit), "\n")
  cat("Total response is", x[["total.e_response"]], "response units J-1 per",
      as.character(time.unit), "\n")
  cat("Mean response is", x[["mean.q_response"]], "response units mol-1 nm-1 per",
      as.character(time.unit), "\n")
  cat("Total response is", x[["total.q_response"]], "response units mol-1 per",
      as.character(time.unit), "\n")
}

# @describeIn print.summary.generic_spct Print a "summary.chrome.spct" object.
#'
#' @export
#' @rdname print.summary.generic_spct
#'
print.summary.chroma_spct <- function(x, ...) {
  time.unit <- attr(x, "time.unit")
  cat("wavelength ranges from", x[["min.w.length"]], "to", x[["max.w.length"]], "nm \n")
  cat("largest wavelength step size is", x[["w.length.step"]], "nm \n")
  cat(paste("maximum (x, y, z) values are (",
            paste(x[["x.max"]], x[["y.max"]], x[["z.max"]], sep=", "),
            ")", sep=""), "\n")
}

# Color -------------------------------------------------------------------

#' Color of a source_spct object.
#'
#' A function that returns the equivalent RGB colour of an object of class
#' "source_spct".
#'
#' @param x an object of class "source_spct"
#' @param ... not used in current version
#' @export
#'
color.source_spct <- function(x, ...) {
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
#' @family wavelength summaries
#'
range.generic_spct <- function(..., na.rm=FALSE) {
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
#' @family wavelength summaries
#'
max.generic_spct <- function(..., na.rm=FALSE) {
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
#' @family wavelength summaries
#'
min.generic_spct <- function(..., na.rm=FALSE) {
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
#' @family wavelength summaries
stepsize <- function(x, ...) UseMethod("stepsize")

#' @describeIn stepsize Default function usable on numeric vectors.
#' @export stepsize.default
stepsize.default <- function(x, ...) {
  return(range(diff(x)))
}

#' @describeIn stepsize  Method for "generic_spct" objects for generic function.
#'
#' @export stepsize.generic_spct
#'
#' @examples
#' stepsize(sun.spct)
#'
stepsize.generic_spct <- function(x, ...) {
  range(diff(x[["w.length"]]))
}


# Labels ------------------------------------------------------------------

#' Labels of a "generic_spct" object
#'
#' A function to obtain the labels of a spectrum. Currently returns 'names'.
#'
#' @param object an object of generic_spct
#' @param ... not used in current version
#'
#' @export
#'
labels.generic_spct <- function(object, ...) {
  return(names(object))
}
