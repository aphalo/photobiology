#' Integrate spectral data.
#'
#' This function gives the result of integrating spectral data over wavelengths.
#'
#' @param spct generic_spct
#'
#' @return One or more numeric values with no change in scale factor: e.g. [W
#'   m-2 nm-1] -> [W m-2]. Each value in the returned vector corresponds to a
#'   variable in the spectral object, except for wavelength.
#'
#' @export
#' @examples
#'
#' integrate_spct(sun.spct)
#'
integrate_spct <- function(spct) {
  names.spct <- names(spct)
  data.cols <- names.spct[names.spct != "w.length"]
  comment.spct <- comment(spct)
  integrals <- NULL
  for (data.col in data.cols) {
    integrals <- c(integrals,
                   integrate_xy(spct[["w.length"]],
                                        spct[[eval(data.col)]]))
  }
  names(integrals) <- gsub("^s.", x = data.cols, replacement = "")
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
#'   variable in the spectral object, except for wavelength.
#'
#' @export
#' @examples
#'
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
#' @param w.length.out numeric vector of wavelengths (nm)
#' @param fill a value to be assigned to out of range wavelengths
#' @param length.out numeric value
#'
#' @details If \code{length.out} it is a numeric value, then gives the number of
#'   rows in the output, if it is \code{NULL}, the values in the numeric vector
#'   \code{w.length.out} are used. If both are not \code{NULL} then the range of
#'   \code{w.length.out} and \code{length.out} are used to generate a vector of
#'   wavelength. A value of \code{NULL} for \code{fill} prevents extrapolation.
#'   If both \code{w.length.out} and \code{length.out} are \code{NULL} the input
#'   is returned as is. If \code{w.length.out} has length equal to zero, zero
#'   rows from the input are returned.
#'
#' @note The default \code{fill = NA} fills extrapolated values with NA. Giving NULL as
#' argument for \code{fill} deletes wavelengths outside the input data range from the
#' returned spectrum. A numerical value can be also be provided as fill. This function calls
#' \code{interpolate_spectrum} for each non-wavelength column in the input spectra object.
#'
#' @return A new spectral object of the same class as argument \code{spct}.
#'
#' @export
#' @examples
#'
#' interpolate_spct(sun.spct, 400:500, NA)
#' interpolate_spct(sun.spct, 400:500, NULL)
#' interpolate_spct(sun.spct, seq(200, 1000, by=0.1), 0)
#' interpolate_spct(sun.spct, c(400,500), length.out=201)
#'
interpolate_spct <- function(spct,
                             w.length.out = NULL,
                             fill = NA,
                             length.out = NULL) {
  stopifnot(is.generic_spct(spct))
  if (length(w.length.out) == 0 && is.null(length.out)) {
    if (is.null(w.length.out)) {
      # with default we return the input
      return(spct)
    } else {
      # with no wavelengths we return a spectrum of length zero
      return(spct[FALSE, ])
    }
  }
  if (!is.null(w.length.out) && any(is.na(w.length.out))) {
    warning("NAs omitted from 'w.length.out'.")
    w.length.out <- stats::na.omit(w.length.out)
  }
  if (is.null(fill)) {
    if (!is.null(w.length.out)) {
      w.length.out <-
        unique(sort(c(ifelse(min(spct) > min(w.length.out), min(spct), min(w.length.out)),
                      w.length.out,
                      ifelse(max(spct) < max(w.length.out), max(spct), max(w.length.out)))))
      w.length.out <- w.length.out[w.length.out >= min(spct) & w.length.out <= max(spct)]
      fill <- NA_real_
    }
  } else if (is.na(fill)) {
    fill <- NA_real_
  }
  names.spct <- names(spct)
  numeric.cols <- names.spct[sapply(spct[ ,names.spct], is.numeric)]
#  other.cols <- setdiff(names.spct, numeric.cols)
  data.cols <- setdiff(numeric.cols, "w.length")
  if (nrow(spct) == 0) {
    new.spct <- tibble::tibble(w.length = w.length.out)
    new.spct[ , data.cols] <- fill
  }
  if (!is.null(length.out)) {
    if (!is.numeric(length.out) ||
        length(length.out) == 0 ||
        is.na(length.out) ||
        length.out == 0) {
      return(spct[FALSE, numeric.cols]) # same columns as in other cases
    } else{
      length.out <- round(length.out, 0)
    }
  }
  if (is.null(w.length.out)) {
    if (is.null(length.out)) {
      return(spct[ , numeric.cols]) # same columns as in other cases
    } else {
      w.length.out <- range(spct)
    }
  }
  if (length(w.length.out) == 0) {
    return(spct[FALSE, numeric.cols]) # same columns as in other cases
  }
  if (length(w.length.out) == 0) {
    # we can get here only if all(is.na(w.length.out))
    out.spct <- spct[1, numeric.cols]
    out.spct[ , ] <- NA_real_
    return(out.spct)
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
  new.spct <- tibble::tibble(w.length = w.length.out)

  for (data.col in data.cols) {
    temp.values <-  with(spct, get(data.col))
    if (is.numeric(temp.values)) {
      new.values <- interpolate_spectrum(spct$w.length,
                                         temp.values,
                                         w.length.out,
                                         fill)
      new.spct[[data.col]] <- new.values
    }
  }
  setGenericSpct(new.spct)
  copy_attributes(spct, new.spct, copy.class = TRUE)
}

#' @rdname interpolate_spct
#'
#' @param mspct an object of class "generic_mspct"
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
interpolate_mspct <- function(mspct,
                              w.length.out = NULL,
                              fill = NA,
                              length.out = NULL,
                              .parallel = FALSE,
                              .paropts = NULL) {

  msmsply(mspct = mspct,
          .fun = interpolate_spct,
          w.length.out = w.length.out,
          fill = fill,
          length.out = length.out,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' Map spectra to new wavelength values.
#'
#' This function returns the result of interpolating spectral data from the original set of
#' wavelengths to a new one.
#'
#' @param x an R object
#' @param w.length.out numeric vector of wavelengths (nm)
#' @param fill a value to be assigned to out of range wavelengths
#' @param length.out numeric value
#' @param ... not used
#'
#' @details If \code{length.out} it is a numeric value, then gives the number of rows in the
#' output, if it is \code{NULL}, the values in the numeric vector \code{w.length.out} are used.
#' If both are not \code{NULL} then the range of \code{w.length.out} and \code{length.out} are
#' used to generate a vector of wavelength. A value of \code{NULL} for \code{fill} prevents
#' extrapolation.
#'
#' @note The default \code{fill = NA} fills extrapolated values with NA. Giving NULL as
#' argument for \code{fill} deletes wavelengths outside the input data range from the
#' returned spectrum. A numerical value can be also be provided as fill. This function calls
#' \code{interpolate_spectrum} for each non-wavelength column in the input spectra object.
#'
#' @return A new spectral object of the same class as argument \code{spct}.
#'
#' @family interpolate functions
#' @export
#' @examples
#' interpolate_wl(sun.spct, 400:500, NA)
#' interpolate_wl(sun.spct, 400:500, NULL)
#' interpolate_wl(sun.spct, seq(200, 1000, by=0.1), 0)
#' interpolate_wl(sun.spct, c(400,500), length.out=201)
#'
interpolate_wl <- function(x,
                           w.length.out,
                           fill,
                           length.out,
                           ...) UseMethod("interpolate_wl")

#' @describeIn interpolate_wl Default for generic function
#'
#' @export
#'
interpolate_wl.default <- function(x,
                                   w.length.out,
                                   fill,
                                   length.out,
                                   ...) {
  stop("'interpolate_wl()' is not defined for objects of class '", class(x)[1], "'.")
}

#' @describeIn interpolate_wl  Interpolate wavelength in an object of class
#'   "generic_spct" or derived.
#'
#' @export
#'
interpolate_wl.generic_spct <- function(x,
                                        w.length.out = NULL,
                                        fill = NA,
                                        length.out = NULL,
                                        ...) {
  interpolate_spct(spct = x,
                   w.length.out = w.length.out,
                   fill = fill,
                   length.out = length.out)
}

#' @describeIn interpolate_wl  Interpolate wavelength in an object of class
#'   "generic_mspct" or derived.
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
interpolate_wl.generic_mspct <- function(x,
                                         w.length.out = NULL,
                                         fill = NA,
                                         length.out = NULL,
                                         ...,
                                         .parallel = FALSE,
                                         .paropts = NULL) {
  interpolate_mspct(mspct = x,
                    w.length.out = w.length.out,
                    fill = fill,
                    length.out = length.out,
                    .parallel = .parallel,
                    .paropts = .paropts)
}
