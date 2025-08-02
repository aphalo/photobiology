#' Integrate spectral data.
#'
#' This function gives the result of integrating spectral data over wavelengths.
#'
#' @param spct generic_spct
#'
#' @return One or more numeric values with no change in scale factor: e.g. [W
#'   m-2 nm-1] -> [W m-2]. Each value in the returned vector corresponds to a
#'   variable in the spectral object, except for wavelength. For non-numeric
#'   variables the returned value is \code{NA}.
#'
#' @export
#' @examples
#'
#' integrate_spct(sun.spct)
#'
integrate_spct <- function(spct) {

  # we look for multiple spectra in long form
  num.spectra <- getMultipleWl(spct)
  if (num.spectra != 1) {
    warning("Integrating ", num.spectra, " spectra as one.")
  }

  names.spct <- names(spct)
  data.cols <- names.spct[names.spct != "w.length"]
  comment.spct <- comment(spct)
  integrals <- NULL
  for (data.col in data.cols) {
    if (is.numeric(spct[[eval(data.col)]])) {
      integrals <- c(integrals,
                     integrate_xy(spct[["w.length"]],
                                  spct[[eval(data.col)]]))
    } else {
      message("Skipping non-numeric column ", data.col)
      integrals <- c(integrals, NA_real_)
    }
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
  return(integrate_spct(spct) / (wl_max(spct) - wl_min(spct)))
}

#' Map a spectrum to new wavelength values.
#'
#' This function gives the result of interpolating spectral data from the original set of
#' wavelengths to a new one.
#'
#' @inheritParams interpolate_spectrum
#' @param spct generic_spct
#' @param length.out integer Length of the wavelength vector in the returned
#'   value. Overrides \code{w.length.out} is not \code{NULL}, respects its
#'   range but overrides the actual values.
#'
#' @inherit interpolate_spectrum details
#'
#' @return A new spectral object of the same class as argument \code{spct} with
#'   a different number of rows than \code{x}, different \code{w.length} values
#'   and new numeric values for spectral data obtained by interpolation.
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
                             length.out = NULL,
                             method = "approx",
                             ...) {

  stopifnot(is.generic_spct(spct))

  # we look for multiple spectra in long form
  num.spectra <- getMultipleWl(spct)
  if (num.spectra != 1) {
    stop("Cannot interpolate ", num.spectra, " spectra in long form")
  }

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
        unique(sort(c(ifelse(wl_min(spct) > min(w.length.out), wl_min(spct), min(w.length.out)),
                      w.length.out,
                      ifelse(wl_max(spct) < max(w.length.out), wl_max(spct), max(w.length.out)))))
      w.length.out <- w.length.out[w.length.out >= wl_min(spct) & w.length.out <= wl_max(spct)]
      fill <- NA_real_
    }
  } else if (is.na(fill)) {
    fill <- NA_real_
  }
  names.spct <- names(spct)
  numeric.cols <- names.spct[sapply(spct, is.numeric)]
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
  if (!is.null(length.out)  && length.out == 1L) {
    if (is.null(w.length.out)) {
      w.length.out <- midpoint(spct)
    } else {
      w.length.out <- midpoint(w.length.out)
    }
  }
  if (!is.null(length.out) && length.out > 1) {
    if (is.null(w.length.out) || length(w.length.out) < 2L) {
      w.length.out <- seq(wl_min(spct), max(spct), length.out = length.out)
    } else {
      w.length.out <- seq(min(w.length.out), max(w.length.out), length.out = length.out)
    }
  } else if (is.null(w.length.out)) {
    # nothing to do
    return(spct)
  }
  max.spct <- wl_max(spct)
  min.spct <- wl_min(spct)
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
      new.values <- interpolate_spectrum(w.length.in = spct[["w.length"]],
                                         s.irrad = temp.values,
                                         w.length.out = w.length.out,
                                         fill = fill,
                                         method = method)
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
                              method = "approx",
                              ...,
                              .parallel = FALSE,
                              .paropts = NULL) {

  msmsply(mspct = mspct,
          .fun = interpolate_spct,
          w.length.out = w.length.out,
          fill = fill,
          length.out = length.out,
          method = method,
          .parallel = .parallel,
          .paropts = .paropts,
          ...)
}

#' Map spectra to new wavelength values.
#'
#' This method returns the result of interpolating spectral data from the
#' original set of wavelengths to a new one.
#'
#' @inheritParams interpolate_spct
#' @param x an R object
#'
#' @inherit interpolate_spct details
#'
#' @return A new spectral object or collection of spectral objects, of the same
#'   class as argument \code{x}. Each spectrum returned with more or fewer rows
#'   than in \code{x}, the requested new \code{w.length} values and new numeric
#'   values for spectral quantities, obtained by interpolation.
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
                           method,
                           ...) UseMethod("interpolate_wl")

#' @describeIn interpolate_wl Default for generic function
#'
#' @export
#'
interpolate_wl.default <- function(x,
                                   w.length.out,
                                   fill,
                                   length.out,
                                   method,
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
                                        method = "approx",
                                        ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(x) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    # call method on the collection
    return(interpolate_wl(x = mspct,
                          w.length.out = w.length.out,
                          fill = fill,
                          length.out = length.out,
                          method = method,
                          ...))
  }

  interpolate_spct(spct = x,
                   w.length.out = w.length.out,
                   fill = fill,
                   length.out = length.out,
                   method = method)
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
                                         method = "approx",
                                         ...,
                                         .parallel = FALSE,
                                         .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  interpolate_mspct(mspct = x,
                    w.length.out = w.length.out,
                    fill = fill,
                    length.out = length.out,
                    method = method,
                    .parallel = .parallel,
                    .paropts = .paropts)
}
