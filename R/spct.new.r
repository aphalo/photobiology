#' Create a new source.spct
#'
#' This fucntion can be used to create source.spct objects from numeric vectors.
#'
#' @usage source.spct(w.length, s.e.irrad=NULL, s.q.irrad=NULL, time.unit=c("second", "day"), comment=NULL)
#'
#' @param w.length numeric vector with wavelengths in nanometres
#' @param s.e.irrad numeric vecror with spectral energy irradiance in W m-2 nm-1 or J d-1 m-2 nm-1
#' @param s.q.irrad numeric vecror with spectral photon irradiance in mol s-1 m-2 nm-1 or mol d-1 m-2 nm-1
#' @param time.unit character string indicating the time unit used for spectral irradiance or exposure ("second" or "day")
#' @param comment character string to be added as a comment attribute to the created object
#'
#' @return a source.spct object
#'
#' @export
#'
source.spct <- function(w.length, s.e.irrad=NULL, s.q.irrad=NULL, time.unit=c("second", "day"), comment=NULL) {
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
  setSourceSpct(out.spct, time.unit)
  return(out.spct)
}

#' Create a new filter.spct
#'
#' This fucntion can be used to create source.spct objects from numeric vectors.
#'
#' @usage filter.spct(w.length, Tfr=NULL, Tpc=NULL, A=NULL, Tfr.type=c("total", "internal"), comment=NULL)
#'
#' @param w.length numeric vector with wavelengths in nanometres
#' @param Tfr numeric vector with spectral transmittance as fraction of one
#' @param Tpc numeric vector with spectral transmittance as percent values
#' @param A   numeric vector of absorbance values (log10 based)
#' @param Tfr.type character string indicating whether transmittance values are "total" or "internal" values
#' @param comment character string to be added as a comment attribute to the created object
#'
#' @return a source.spct object
#'
#' @note "internal" transmittance is defined as the transmittance of the material body itself, while
#' "total" transmittance includes the effects of surface reflectance on the ammount of light transmitted.
#'
#' @export
#'
filter.spct <- function(w.length, Tfr=NULL, Tpc=NULL, A=NULL, Tfr.type=c("total", "internal"), comment=NULL) {
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
  setFilterSpct(out.spct, Tfr.type)
  return(out.spct)
}

#' Create a new reflector.spct
#'
#' This fucntion can be used to create source.spct objects from numeric vectors.
#'
#' @usage reflector.spct(w.length, Rfr=NULL, Rpc=NULL, comment=NULL)
#'
#' @param w.length numeric vector with wavelengths in nanometres
#' @param Rfr numeric vector with spectral transmittance as fraction of one
#' @param Rpc numeric vector with spectral transmittance as percent values
#' @param comment character string to be added as a comment attribute to the created object
#'
#' @return a source.spct object
#'
#' @note "internal" transmittance is defined as the transmittance of the material body itself, while
#' "total" transmittance includes the effects of surface reflectance on the ammount of light transmitted.
#'
#' @export
#'
reflector.spct <- function(w.length, Rfr=NULL, Rpc=NULL, comment=NULL) {
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
  setReflectorSpct(out.spct)
  return(out.spct)
}
