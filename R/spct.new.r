

# Constructors ------------------------------------------------------------

#' Create new spectral objects
#'
#' These fucntions can be used to create spectral objects derived from
#' \code{generic_spct}. They take as arguments numeric vectors for the data
#' character scalars for attributes, and a logical flag.
#'
#' @param w.length numeric vector with wavelengths in nanometres
#' @param s.e.irrad numeric vector with spectral energy irradiance in [W m-2
#'   nm-1] or [J d-1 m-2 nm-1]
#' @param s.q.irrad numeric A vector with spectral photon irradiance in [mol s-1
#'   m-2 nm-1] or [mol d-1 m-2 nm-1].
#' @param time.unit character string indicating the time unit used for spectral
#'   irradiance or exposure ("second" or "day").
#' @param bswf.used character A string indicating the BSWF used, if any, for
#'   spectral effective irradiance or exposure ("none" or the name of the BSWF).
#' @param comment character A string to be added as a comment attribute to the
#'   object created.
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning.
#'
#' @return An source_spct object
#'
#' @export
#'
#' @family creation of spectral objects functions
#'
#' @rdname source_spct
#'
#' @note The functions can be used to add only one spectral quantity to a
#'   spectral object. Some of the functions have different arguments, for the
#'   same quantity expressed in different units. An actual parameter can be
#'   supplied to only one of these formal paprameters in a given call to the
#'   fucntion.
#'
source_spct <- function(w.length, s.e.irrad=NULL, s.q.irrad=NULL,
                        time.unit=c("second", "day"), bswf.used = c("none", "unknown"),
                        comment=NULL, strict.range=TRUE) {
  if (is.null(s.q.irrad) && (is.numeric(s.e.irrad))) {
    z <- data.table(w.length, s.e.irrad)
  } else if (is.null(s.e.irrad) && (is.numeric(s.q.irrad))) {
    z <- data.table(w.length, s.q.irrad)
  } else {
    warning("One and only one of s.e.irrad or s.q.irrad should be different from NULL.")
    return(NA)
  }
  if (!is.null(comment)) {
    setattr(z, "comment", comment)
  }
  setSourceSpct(z,
                time.unit = time.unit, bswf.used = bswf.used,
                strict.range = strict.range)
  return(z)
}

#' @rdname source_spct
#'
#' @param cps numeric vector with linearized raw counts expressed per second
#'
#' @export
#'
cps_spct <- function(w.length, cps=NULL, comment=NULL) {
  z <- data.table(w.length = w.length, cps = cps)
  if (!is.null(comment)) {
    setattr(z, "comment", comment)
  }
  setCpsSpct(z)
  return(z)
}

#' @rdname source_spct
#'
#' @param s.e.response numeric vector with spectral energy irradiance in W m-2
#'   nm-1 or J d-1 m-2 nm-1
#' @param s.q.response numeric vector with spectral photon irradiance in mol s-1
#'   m-2 nm-1 or mol d-1 m-2 nm-1
#'
#' @export
#'
response_spct <- function(w.length, s.e.response=NULL, s.q.response=NULL, time.unit=c("second", "day"), comment=NULL) {
  if (is.null(s.q.response) && (is.numeric(s.e.response))) {
    z <- data.table(w.length, s.e.response)
  } else if (is.null(s.e.response) && (is.numeric(s.q.response))) {
    z <- data.table(w.length, s.q.response)
  } else {
    warning("One and only one of s.e.response or s.q.response should be different from NULL.")
    return(NA)
  }
  if (!is.null(comment)) {
    setattr(z, "comment", comment)
  }
  setResponseSpct(z, time.unit)
  return(z)
}

#' @rdname source_spct
#'
#' @param Tfr numeric vector with spectral transmittance as fraction of one
#' @param Tpc numeric vector with spectral transmittance as percent values
#' @param A   numeric vector of absorbance values (log10 based)
#' @param Tfr.type character string indicating whether transmittance values are
#'   "total" or "internal" values
#'
#' @note "internal" transmittance is defined as the transmittance of the
#'   material body itself, while "total" transmittance includes the effects of
#'   surface reflectance on the ammount of light transmitted.
#'
#' @export
#'
filter_spct <- function(w.length, Tfr=NULL, Tpc=NULL, A=NULL, Tfr.type=c("total", "internal"),
                        comment=NULL, strict.range=TRUE) {
  if (is.null(Tpc) && is.null(A) && is.numeric(Tfr)) {
    z <- data.table(w.length, Tfr)
  } else if (is.null(Tfr) && is.null(A) && is.numeric(Tpc)) {
    z <- data.table(w.length, Tpc)
  } else if (is.null(Tpc) && is.null(Tfr) && is.numeric(A)) {
    z <- data.table(w.length, A)
  } else {
    warning("One and only one of Tfr, Tpc or Abs should be different from NULL.")
    return(NA)
  }
  if (!is.null(comment)) {
    setattr(z, "comment", comment)
  }
  setFilterSpct(z, Tfr.type, strict.range = strict.range)
  return(z)
}

#' @rdname source_spct
#'
#' @param Rfr numeric vector with spectral refletance as fraction of one
#' @param Rpc numeric vector with spectral reflectance as percent values
#' @param Rfr.type character A string, either "total" or "specular".
#'
#' @export
#'
reflector_spct <- function(w.length, Rfr=NULL, Rpc=NULL,
                           Rfr.type=c("total", "specular"),
                           comment=NULL, strict.range=TRUE) {
  if (is.null(Rpc) && is.numeric(Rfr)) {
    z <- data.table(w.length, Rfr)
  } else if (is.null(Rfr) && is.numeric(Rpc)) {
    z <- data.table(w.length, Rpc)
  } else {
    warning("One and only one of Rfr, or Rpc should be different from NULL.")
    return(NA)
  }
  if (!is.null(comment)) {
    setattr(z, "comment", comment)
  }
  setReflectorSpct(z, Rfr.type = Rfr.type, strict.range=strict.range)
  return(z)
}

#' @rdname source_spct
#'
#' @export
#'
object_spct <- function(w.length, Rfr=NULL, Tfr=NULL,
                        Tfr.type=c("total", "internal"),
                        Rfr.type=c("total", "specular"),
                        comment=NULL, strict.range=TRUE) {
  z <- data.table(w.length, Rfr, Tfr)
  if (!is.null(comment)) {
    setattr(z, "comment", comment)
  }
  setObjectSpct(z,
                Tfr.type = Tfr.type,
                Rfr.type = Rfr.type,
                strict.range=strict.range)
  return(z)
}


# merge -------------------------------------------------------------------


#' Merge two generic_spct objects
#'
#' Relatively quick merge of two spct objects based on w.length.
#'
#' @param x generic_spct (or derived) objects to be merged
#' @param y generic_spct (or derived) objects to be merged
#' @param by a vector of shared column names in \code{x} and \code{y} to merge on;
#' \code{by} defaults to \code{w.length}.
#' @param ... other arguments passed to \code{merge.data.table}
#'
#' @note if the class of x and y is the same, it is preserved, but
#' if it differs \code{generic_spct} is used for the returned value,
#' except when x and y, are one each of classes reflector_spct and
#' filter_spct in which case an object_spct is returned.
#' In the current implementation only wavelengths values shared
#' by x and y are preserved.
#'
#' @seealso \code{\link[data.table]{merge}}
#'
#' @export
#'
merge.generic_spct <- function(x, y, by = "w.length", ...) {
  if (identical(class_spct(x), class_spct(y))) {
    z <- data.table:::merge.data.table(x, y, by = by, ...)
    setattr(z, "class", class(x))
  } else if (is.filter_spct(x) && is.reflector_spct(y)) {
    xx <- A2T(x, action = "replace", byref = FALSE)
    z <- data.table:::merge.data.table(xx, y, by = "w.length", ...)
    setObjectSpct(z, Tfr.type = getTfrType(x), Rfr.type = getRfrType(y))
  } else if (is.reflector_spct(x) && is.filter_spct(y)) {
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
                       class_spct(x),
                       "\n\ncomment(y):\n",
                       comment(y),
                       "\nclass: ",
                       class_spct(y))
  setattr(z, "comment", new.comment)
  return(z)
}

