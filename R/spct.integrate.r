#' Integrate spectral data.
#'
#' This function gives the result of integrating spectral data over wavelengths.
#'
#' @param spct generic_spct
#'
#' @return One or more numeric values with no change in scale factor: e.g. [W
#'   m-2 nm-1] -> [W m-2]. Each value in the returned vector corresponds to a
#'   variable in the spectral object, except for wavelenght.
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.spct)
#' integrate_spct(sun.spct)
#'
integrate_spct <- function(spct) {
  names.spct <- names(spct)
  names.data <- names.spct[names.spct != "w.length"]
  comment.spct <- comment(spct)
  integrals <- NULL
  for (data.col in names.data) {
    integrals <- c(integrals, integrate_irradiance(spct[["w.length"]], spct[[eval(data.col)]]))
  }
  names(integrals) <- gsub("^s.", x = names.data, replacement = "")
  comment(integrals) <- comment.spct
  return(integrals)
}

#' Average spectral data.
#'
#' This function gives the result of integrating spectral data over
#' wavelengths and dividing the result by the spread or span of the
#' wavelengths.
#'
#' @param spct generic_spct
#'
#' @return One or more numeric values with no change in scale factor: e.g. [W
#'   m-2 nm-1] -> [W m-2 nm-1]. Each value in the returned vector corresponds to a
#'   variable in the spectral object, except for wavelenght.
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.spct)
#' average_spct(sun.spct)
#'
average_spct <- function(spct) {
  return(integrate_spct(spct) / (max(spct) - min(spct)))
}

#' Map a spectrum to new wavelength values.
#'
#' This function gives the result of interpolating spectral data from the original set of
#' wavelengths to a new one.
#'
#' @param spct generic_spct
#' @param w.length.out numeric array of wavelengths (nm)
#' @param fill a value to be assigned to out of range wavelengths
#' @param length.out numeric value
#'
#' @details If \code{length.out} it is a numeric value, then gives the number of rows in the
#' output, if it is \code{NULL}, the values in the numeric vector \code{w.length.out} are used.
#' If both are not \code{NULL} then the range of \code{w.length.out} and \code{length.out} are
#' used to generate a vector of wavelength. A value of \code{NULL} for \code{fill} prevents
#' extrapolation.
#'
#' @note The default \code{fill = NA} fills extrpolated values with NA. Giving NULL as
#' argument for \code{fill} deletes wavelengths outside the input data range from the
#' returned spectrum. A numerical value can be also be provided as fill. This function calls
#' \code{interpolate_spectrum} for each non-wavelength column in the input spectra object.
#'
#' @return A new spectral object of the same class as argument \code{spct}.
#'
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.spct)
#' interpolate_spct(sun.spct, 400:500, NA)
#' interpolate_spct(sun.spct, 400:500, NULL)
#' interpolate_spct(sun.spct, seq(200, 1000, by=0.1), 0)
#' interpolate_spct(sun.spct, c(400,500), length.out=201)
#'
interpolate_spct <- function(spct,
                             w.length.out=NULL,
                             fill=NA,
                             length.out=NULL) {
  stopifnot(is.any_spct(spct))
  if (!is.null(length.out) && (is.na(length.out) || length.out < 1L) ) {
    return(spct[NA])
  }
  if (length(w.length.out) > 1L) {
    step.ratio <- stepsize(spct)[1] / max(diff(w.length.out), na.rm = TRUE)
    if (step.ratio < 1 && length(spct) > 100) {
      # this a temporary kludge as the degree of smoothing is not well tuned and only
      # tested with the solar spectrum.
      warning("Smoothing before interpolation, as w.length.out is more sparse./nIt could be better to use the original data as is.")
      spct <- smooth_spct(spct, method = "supsmu", strength = step.ratio * 1e-2)
    }
  }
  class_spct <- class(spct)
  if (!is.null(length.out)  && length.out == 1L) {
    if (is.null(w.length.out)) {
      w.length.out <- midpoint(spct)
    } else {
      w.length.out <- midpoint(w.length.out)
    }
  }
  if (!is.null(length.out) && length.out > 1) {
    if (is.null(w.length.out) || length(w.length.out) < 2L) {
      w.length.out <- seq(min(spct), max(spct), length.out = length.out)
    } else {
      w.length.out <- seq(min(w.length.out), max(w.length.out), length.out = length.out)
    }
  } else if (is.null(w.length.out)) {
    # nothing to do
    return(spct)
  }
  names.spct <- names(spct)
  names.data <- setdiff(names.spct, "w.length")
  max.spct <- max(spct)
  min.spct <- min(spct)
  max.wl.out <- max(w.length.out)
  min.wl.out <- min(w.length.out)
  if (min.spct > min.wl.out && min.spct < max.wl.out) w.length.out <- c(min.spct, w.length.out)
  if (max.spct < max.wl.out && max.spct > min.wl.out) w.length.out <- c(w.length.out, max.spct)
  if (is.null(fill)) {
    w.length.out <- w.length.out[w.length.out >= min.spct & w.length.out <= max.spct]
    if (length(w.length.out) == 0) {
      return(spct[NA])
    }
  }
  w.length.out <- unique(sort(w.length.out))
  new.spct <- dplyr::data_frame(w.length = w.length.out)

  for (data.col in names.data) {
    temp.values <-  with(spct, get(data.col))
    if (is.numeric(temp.values)) {
      new.values <- interpolate_spectrum(spct$w.length,
                                         temp.values,
                                         w.length.out,
                                         fill)
      new.spct[[data.col]] <- new.values
    }
  }
  if (class_spct[1] == "source_spct") {
    setSourceSpct(new.spct,
                  time.unit = getTimeUnit(spct),
                  bswf.used = getBSWFUsed(spct))
  } else if (class_spct[1] == "filter_spct") {
    setFilterSpct(new.spct,
                  Tfr.type = getTfrType(spct))
  } else if (class_spct[1] == "reflector_spct") {
    setReflectorSpct(new.spct,
                     Rfr.type = getRfrType(spct))
  } else if (class_spct[1] == "object_spct") {
    setObjectSpct(new.spct,
                  Tfr.type = getTfrType(spct),
                  Rfr.type = getRfrType(spct))
  } else if (class_spct[1] == "response_spct") {
    setResponseSpct(new.spct,
                    time.unit = getTimeUnit(spct)
    )
  } else if (class_spct[1] == "chroma_spct") {
    setChromaSpct(new.spct)
  } else if (class_spct[1] == "generic_spct") {
    setGenericSpct(new.spct)
  }
  setNormalized(new.spct, getNormalized(spct))
  setScaled(new.spct, getScaled(spct))
  comment(new.spct) <- comment(spct)
  return(new.spct)
}

#' @rdname interpolate_spct
#'
#' @param mspct an object of class "generic_mspct"
#'
#' @export
#'
interpolate_mspct <- function(mspct,
                             w.length.out=NULL,
                             fill=NA,
                             length.out=NULL) {

  msmsply(mspct = mspct,
          .fun = interpolate_spct,
          w.length.out = w.length.out,
          fill = fill,
          length.out = length.out)
}

#' Map a spectrum to new wavelength values.
#'
#' This function gives the result of interpolating spectral data from the original set of
#' wavelengths to a new one.
#'
#' @param spct generic_spct
#' @param w.length.out numeric array of wavelengths (nm)
#' @param fill a value to be assigned to out of range wavelengths
#' @param length.out numeric value
#'
#' @details If \code{length.out} it is a numeric value, then gives the number of rows in the
#' output, if it is \code{NULL}, the values in the numeric vector \code{w.length.out} are used.
#' If both are not \code{NULL} then the range of \code{w.length.out} and \code{length.out} are
#' used to generate a vector of wavelength. A value of \code{NULL} for \code{fill} prevents
#' extrapolation.
#'
#' @note The default \code{fill = NA} fills extrpolated values with NA. Giving NULL as
#' argument for \code{fill} deletes wavelengths outside the input data range from the
#' returned spectrum. A numerical value can be also be provided as fill. This function calls
#' \code{interpolate_spectrum} for each non-wavelength column in the input spectra object.
#'
#' @return A new spectral object of the same class as argument \code{spct}.
#'
#' @export
#' @examples
#' data(sun.spct)
#' interpolate_wl(sun.spct, 400:500, NA)
#' interpolate_wl(sun.spct, 400:500, NULL)
#' interpolate_wl(sun.spct, seq(200, 1000, by=0.1), 0)
#' interpolate_wl(sun.spct, c(400,500), length.out=201)
#'

#' Clip head and/or tail of a spectrum
#'
#' Clipping of head and tail of a spectrum based on wavelength limits, no
#' interpolation used.
#'
#' @param x an R object
#' @param range a numeric vector of length two, or any other object for which
#'   function \code{range()} will return range of walengths expressed in
#'   nanometres.
#' @param ... not used
#'
#' @return an R object of same class as input, most frequently of a shorter
#'   length, and never longer.
#'
#' @note The condition tested is \code{wl >= range[1] & wl < (range[2] + 1e-13)}.
#'
#' @family trim functions
#' @export
#' @examples
#' interpolate_wl(sun.spct, range = c(400, 500))
#' interpolate_wl(sun.spct, range = c(NA, 500))
#' interpolate_wl(sun.spct, range = c(400, NA))
#'
interpolate_wl <- function(x, range, ...) UseMethod("interpolate_wl")

#' @describeIn interpolate_wl Default for generic function
#'
#' @export
#'
interpolate_wl.default <- function(x, range, ...) {
  warning("'interpolate_wl' is not defined for objects of class ", class(x)[1])
  x
}

#' @describeIn interpolate_wl Clip an object of class "generic_spct" or derived.
#'
#' @export
#'
interpolate_wl.generic_spct <- function(x, range = NULL, ...) {
  if (is.null(range)) {
    return(x)
  }
  guard <- 1e-13
  stopifnot(is.any_spct(x))
  stopifnot(!all(is.na(range)))
  if (is.numeric(range) && length(range) == 2) {
    if (is.na(range[1])) {
      x[x[["w.length"]] < range[2] + guard, ]
    } else if (is.na(range[2])) {
      x[x[["w.length"]] >= range[1], ]
    } else {
      x[x[["w.length"]] >= range[1] & x[["w.length"]] < range[2] + guard, ]
    }
  } else {
    range = range(range)
    x[x[["w.length"]] >= range[1] & x[["w.length"]] < range[2] + guard, ]
  }
}

#' @describeIn interpolate_wl  Clip an object of class "generic_mspct" or derived.
#'
#' @export
#'
interpolate_wl.generic_mspct <- function(x, range = NULL, ...) {
  msmsply(mspct = x,
          .fun = interpolate_wl,
          range = range)
}
