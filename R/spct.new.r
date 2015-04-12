

# Constructors ------------------------------------------------------------

#' Create a new source.spct
#'
#' This fucntion can be used to create source.spct objects from numeric vectors.
#'
#' @usage source.spct(w.length, s.e.irrad=NULL, s.q.irrad=NULL,
#'   time.unit=c("second", "day"), bswf.used = c("none", "unknown"),
#'   comment=NULL, strict.range=TRUE)
#'
#' @param w.length numeric vector with wavelengths in nanometres
#' @param s.e.irrad numeric vecror with spectral energy irradiance in W m-2 nm-1
#'   or J d-1 m-2 nm-1
#' @param s.q.irrad numeric vecror with spectral photon irradiance in mol s-1
#'   m-2 nm-1 or mol d-1 m-2 nm-1
#' @param time.unit character string indicating the time unit used for spectral
#'   irradiance or exposure ("second" or "day")
#' @param bswf.used character string indicating the BSWF used, if any, for
#'   spectral effective irradiance or exposure ("second" or "day")
#' @param comment character string to be added as a comment attribute to the
#'   created object
#' @param strict.range logical indicating whether off-range values result in an
#'   error instead of a warning
#'
#' @return An source.spct object
#'
#' @export
#'
#' @family creation of spectral objects functions
#'
source.spct <- function(w.length, s.e.irrad=NULL, s.q.irrad=NULL,
                        time.unit=c("second", "day"), bswf.used = c("none", "unknown"),
                        comment=NULL, strict.range=TRUE) {
  if (is.null(s.q.irrad) && (is.numeric(s.e.irrad))) {
    out.spct <- data.table(w.length, s.e.irrad)
  } else if (is.null(s.e.irrad) && (is.numeric(s.q.irrad))) {
    out.spct <- data.table(w.length, s.q.irrad)
  } else {
    warning("One and only one of s.e.irrad or s.q.irrad should be different from NULL.")
    return(NA)
  }
  if (!is.null(comment)) {
    setattr(out.spct, "comment", comment)
  }
  setSourceSpct(out.spct,
                time.unit = time.unit, bswf.used = bswf.used,
                strict.range = strict.range)
  return(out.spct)
}

#' Create a new response.spct
#'
#' This fucntion can be used to create response.spct objects from numeric
#' vectors.
#'
#' @usage response.spct(w.length, s.e.response=NULL, s.q.response=NULL,
#'   time.unit=c("second", "day"), comment=NULL)
#'
#' @param w.length numeric vector with wavelengths in nanometres
#' @param s.e.response numeric vecror with spectral energy irradiance in W m-2
#'   nm-1 or J d-1 m-2 nm-1
#' @param s.q.response numeric vecror with spectral photon irradiance in mol s-1
#'   m-2 nm-1 or mol d-1 m-2 nm-1
#' @param time.unit character string indicating the time unit used for spectral
#'   irradiance or exposure ("second" or "day")
#' @param comment character string to be added as a comment attribute to the
#'   created object
#'
#' @return A response.spct object
#'
#' @export response.spct
#'
#' @family creation of spectral objects functions
#'
response.spct <- function(w.length, s.e.response=NULL, s.q.response=NULL, time.unit=c("second", "day"), comment=NULL) {
  if (is.null(s.q.response) && (is.numeric(s.e.response))) {
    out.spct <- data.table(w.length, s.e.response)
  } else if (is.null(s.e.response) && (is.numeric(s.q.response))) {
    out.spct <- data.table(w.length, s.q.response)
  } else {
    warning("One and only one of s.e.response or s.q.response should be different from NULL.")
    return(NA)
  }
  if (!is.null(comment)) {
    setattr(out.spct, "comment", comment)
  }
  setResponseSpct(out.spct, time.unit)
  return(out.spct)
}

#' Create a new filter.spct
#'
#' This fucntion can be used to create filter.spct objects from numeric vectors.
#'
#' @usage filter.spct(w.length, Tfr=NULL, Tpc=NULL, A=NULL, Tfr.type=c("total",
#'   "internal"), comment=NULL, strict.range=TRUE)
#'
#' @param w.length numeric vector with wavelengths in nanometres
#' @param Tfr numeric vector with spectral transmittance as fraction of one
#' @param Tpc numeric vector with spectral transmittance as percent values
#' @param A   numeric vector of absorbance values (log10 based)
#' @param Tfr.type character string indicating whether transmittance values are
#'   "total" or "internal" values
#' @param comment character string to be added as a comment attribute to the
#'   created object
#' @param strict.range logical indicating whether off-range values result in an
#'   error instead of a warning
#'
#' @return A filter.spct object
#'
#' @note "internal" transmittance is defined as the transmittance of the
#'   material body itself, while "total" transmittance includes the effects of
#'   surface reflectance on the ammount of light transmitted.
#'
#' @export
#' @family creation of spectral objects functions
#'
#'
filter.spct <- function(w.length, Tfr=NULL, Tpc=NULL, A=NULL, Tfr.type=c("total", "internal"),
                        comment=NULL, strict.range=TRUE) {
  if (is.null(Tpc) && is.null(A) && is.numeric(Tfr)) {
    out.spct <- data.table(w.length, Tfr)
  } else if (is.null(Tfr) && is.null(A) && is.numeric(Tpc)) {
    out.spct <- data.table(w.length, Tpc)
  } else if (is.null(Tpc) && is.null(Tfr) && is.numeric(A)) {
    out.spct <- data.table(w.length, A)
  } else {
    warning("One and only one of Tfr, Tpc or Abs should be different from NULL.")
    return(NA)
  }
  if (!is.null(comment)) {
    setattr(out.spct, "comment", comment)
  }
  setFilterSpct(out.spct, Tfr.type, strict.range = strict.range)
  return(out.spct)
}

#' Create a new reflector.spct
#'
#' This fucntion can be used to create reflector.spct objects from numeric
#' vectors.
#'
#' @usage reflector.spct(w.length, Rfr=NULL, Rpc=NULL, comment=NULL,
#'   strict.range=TRUE)
#'
#' @param w.length numeric vector with wavelengths in nanometres
#' @param Rfr numeric vector with spectral refletance as fraction of one
#' @param Rpc numeric vector with spectral reflectance as percent values
#' @param comment character string to be added as a comment attribute to the
#'   created object
#' @param strict.range logical indicating whether off-range values result in an
#'   error instead of a warning
#'
#' @return A reflector.spct object
#'
#' @note "internal" transmittance is defined as the transmittance of the
#'   material body itself, while "total" transmittance includes the effects of
#'   surface reflectance on the ammount of light transmitted.
#'
#' @export
#'
#' @family creation of spectral objects functions
#'
reflector.spct <- function(w.length, Rfr=NULL, Rpc=NULL, comment=NULL, strict.range=TRUE) {
  if (is.null(Rpc) && is.numeric(Rfr)) {
    out.spct <- data.table(w.length, Rfr)
  } else if (is.null(Rfr) && is.numeric(Rpc)) {
    out.spct <- data.table(w.length, Rpc)
  } else {
    warning("One and only one of Rfr, or Rpc should be different from NULL.")
    return(NA)
  }
  if (!is.null(comment)) {
    setattr(out.spct, "comment", comment)
  }
  setReflectorSpct(out.spct, strict.range=strict.range)
  return(out.spct)
}

#' Create a new object.spct
#'
#' This fucntion can be used to create object.spct objects from numeric vectors.
#' Such objects can be used to describe the optical properties of objects that
#' both reflect and transmit radiation.
#'
#' @usage object.spct(w.length, Rfr=NULL, Tfr=NULL, Tfr.type=c("total",
#'   "internal"), Rfr.type=c("total", "specular"), comment=NULL,
#'   strict.range=TRUE)
#'
#' @param w.length numeric vector with wavelengths in nanometres
#' @param Rfr numeric vector with spectral Reflectance as fraction of one
#' @param Tfr numeric vector with spectral transmittance as fraction of one
#' @param Tfr.type a character string, either "total" or "internal"
#' @param Rfr.type a character string, either "total" or "specular"
#' @param comment character string to be added as a comment attribute to the
#'   created object
#' @param strict.range logical indicating whether off-range values result in an
#'   error instead of a warning
#'
#' @return An object.spct object
#'
#' @note "internal" transmittance is defined as the transmittance of the
#'   material body itself, while "total" transmittance includes the effects of
#'   surface reflectance on the ammount of light transmitted.
#'
#' @export
#'
#' @family creation of spectral objects functions
#'
object.spct <- function(w.length, Rfr=NULL, Tfr=NULL,
                        Tfr.type=c("total", "internal"),
                        Rfr.type=c("total", "specular"),
                        comment=NULL, strict.range=TRUE) {
  out.spct <- data.table(w.length, Rfr, Tfr)
  if (!is.null(comment)) {
    setattr(out.spct, "comment", comment)
  }
  setObjectSpct(out.spct,
                Tfr.type = Tfr.type,
                Rfr.type = Rfr.type,
                strict.range=strict.range)
  return(out.spct)
}


# merge -------------------------------------------------------------------


#' Merge two generic.spct objects
#'
#' Relatively quick merge of two spct objects based on w.length.
#'
#' @param x generic.spct (or derived) objects to be merged
#' @param y generic.spct (or derived) objects to be merged
#' @param by a vector of shared column names in \code{x} and \code{y} to merge on;
#' \code{by} defaults to \code{w.length}.
#' @param ... other arguments passed to \code{merge.data.table}
#'
#' @note if the class of x and y is the same, it is preserved, but
#' if it differs \code{generic.spct} is used for the returned value,
#' except when x and y, are one each of classes reflector.spct and
#' filter.spct in which case an object.spct is returned.
#' In the current implementation only wavelengths values shared
#' by x and y are preserved.
#'
#' @seealso \code{\link[data.table]{merge}}
#'
#' @export
#'
merge.generic.spct <- function(x, y, by = "w.length", ...) {
  if (identical(class.spct(x), class.spct(y))) {
    z <- data.table:::merge.data.table(x, y, by = by, ...)
    setattr(z, "class", class(x))
  } else if (is.filter.spct(x) && is.reflector.spct(y)) {
    xx <- A2T(x, action = "replace", byref = FALSE)
    z <- data.table:::merge.data.table(xx, y, by = "w.length", ...)
    setObjectSpct(z, Tfr.type = getTfrType(x), Rfr.type = getRfrType(y))
  } else if (is.reflector.spct(x) && is.filter.spct(y)) {
    yy <- A2T(y, action = "replace", byref = FALSE)
    z <- data.table:::merge.data.table(xx, yy, by = "w.length", ...)
    setObjectSpct(z, Tfr.type = getTfrType(y), Rfr.type = getRfrType(x))
  } else {
    z <- data.table:::merge.data.table(x, y, by = "w.length", ...)
    setGenericSpct(z)
  }
  new.comment <- paste("Merged spectrum\ncomment(x):\n",
                       comment(x),
                       "\nclass: ",
                       class.spct(x),
                       "\n\ncomment(y):\n",
                       comment(y),
                       "\nclass: ",
                       class.spct(y))
  setattr(z, "comment", new.comment)
  return(z)
}

