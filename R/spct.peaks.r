#' Find peaks in a spectrum
#'
#' This function finds all peaks (local maxima) in a spectrum, using a user
#' selectable size threshold relative to the tallest peak (global maximum). This
#' a wrapper built on top of function peaks from package splus2R.
#'
#' @param x numeric vector
#' @param ignore_threshold numeric value between 0.0 and 1.0 indicating the size
#'   threshold below which peaks will be ignored.
#' @param span a peak is defined as an element in a sequence which is greater
#'   than all other elements within a window of width span centered at that
#'   element. The default value is 3, meaning that a peak is bigger than both of
#'   its neighbors. Default: 3.
#' @param strict logical flag: if TRUE, an element must be strictly greater than
#'   all other values in its window to be considered a peak. Default: TRUE.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for peaks.
#'
#' @return an object like s.irrad of logical values. Values that are TRUE
#'   correspond to local peaks in the data.
#'
#' @export
#' @examples
#' with(sun.data, w.length[find_peaks(s.e.irrad)])
#'
#' @note This function is a wrapper built on function
#'   \code{\link[splus2R]{peaks}} from \pkg{splus2R} and handles non-finite
#'   (including NA) values differently than \code{peaks}, instead of giving an
#'   error they are replaced with the smallest finite value in \code{x}.
#'
#' @seealso \code{\link[splus2R]{peaks}}
#'
#' @family peaks and valleys functions
#'
find_peaks <-
  function(x,
           ignore_threshold = 0.0,
           span = 3,
           strict = TRUE,
           na.rm = FALSE) {
    if (na.rm) {
      x <- na.omit(x)
    }
    if(is.null(span)) {
      return(x == max(x))
    }
    range_x <- range(x, finite = TRUE)
    min_x <- range_x[1]
    max_x <- range_x[2]
    x <- ifelse(!is.finite(x), min_x, x)
    # the next two lines cater for the case when max_x < 0, which is quite common with logs
    delta <- max_x - min_x
    top_flag <- ignore_threshold > 0.0
    scaled_threshold <- delta * abs(ignore_threshold)
    pks <- splus2R::peaks(x = x, span = span, strict = strict)
    if (abs(ignore_threshold) < 1e-5)
      return(pks)
    if (top_flag) {
      return(ifelse(x - min_x > scaled_threshold, pks , FALSE))
    } else {
      return(ifelse(max_x - x > scaled_threshold, pks , FALSE))
    }
  }

#' Get peaks and valleys in a spectrum
#'
#' These functions find peaks (local maxima) or valleys (local minima) in a
#' spectrum, using a user selectable size threshold relative to the tallest peak
#' (global maximum). This a wrapper built on top of function peaks from package
#' splus2R.
#'
#' @param x numeric
#' @param y numeric
#' @param ignore_threshold numeric Value between 0.0 and 1.0 indicating the
#'   relative size compared to tallest peak or deepest valley of the peaks
#'   to return.
#' @param span numeric A peak is defined as an element in a sequence which is
#'   greater than all other elements within a window of width \code{span}
#'   centered at that element. For example, a value of 3 means that a peak is
#'   bigger than both of its neighbors.
#' @param strict logical Flag: if TRUE, an element must be strictly greater than
#'   all other values in its window to be considered a peak. Default: TRUE.
#' @param x_unit character Vector of texts to be pasted at end of labels built
#'   from x value at peaks.
#' @param x_digits numeric Number of significant digits in wavelength label.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for peaks.
#'
#' @return A data frame with variables w.length and s.irrad with their values at
#'   the peaks or valleys plus a character variable of labels.
#'
#' @export
#' @examples
#' with(sun.spct, get_peaks(w.length, s.e.irrad))
#' with(sun.spct, get_valleys(w.length, s.e.irrad))
#'
#' @family peaks and valleys functions
#'
get_peaks <- function(x,
                      y,
                      ignore_threshold = 0.0,
                      span = 5,
                      strict = TRUE,
                      x_unit = "",
                      x_digits = 3,
                      na.rm = FALSE) {
  stopifnot(length(x) == length(y))
  selector <- find_peaks(x = y,
                         ignore_threshold = ignore_threshold,
                         span = span,
                         strict = strict,
                         na.rm = na.rm)
  if (sum(selector) < 1) {
    return(data.frame(
      x = numeric(0),
      y = numeric(0),
      label = character(0)
    ))
  } else {
    peaks.x <- x[selector]
    peaks.y <- y[selector]
    return(data.frame(
      x = peaks.x,
      y = peaks.y,
      label = paste(as.character(signif(
        x = peaks.x, digits = x_digits
      )), x_unit, sep = "")
    ))
  }
}

#' @rdname get_peaks
#' @export
#'
get_valleys <- function(x, y,
                        ignore_threshold = 0.0,
                        span = 5,
                        strict = TRUE,
                        x_unit = "",
                        x_digits = 3,
                        na.rm = FALSE) {
  xy.data <- get_peaks(x = x, y = -y,
                       ignore_threshold = -ignore_threshold,
                       span = span,
                       strict = strict,
                       x_unit = x_unit,
                       x_digits = x_digits,
                       na.rm = na.rm)
  xy.data$y <- -xy.data$y
  return(xy.data)
}


# peaks -------------------------------------------------------------------

#' Peaks or local maxima
#'
#' Function that returns a subset of an R object with observations corresponding
#' to local maxima.
#'
#' @param x an R object
#' @param ignore_threshold numeric value between 0.0 and 1.0 indicating the
#'   relative size compared to tallest peak threshold below which peaks will be
#'   ignored.
#' @param span a peak is defined as an element in a sequence which is greater
#'   than all other elements within a window of width span centered at that
#'   element. The default value is 3, meaning that a peak is bigger than both of
#'   its neighbors. Default: 3.
#' @param strict logical flag: if TRUE, an element must be strictly greater than
#'   all other values in its window to be considered a peak. Default: TRUE.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for peaks.
#' @param ... ignored
#'
#' @return A subset of \code{x} with rows corresponding to local maxima.
#'
#' @export
#'
#' @examples
#' peaks(sun.spct, span = 50)
#' peaks(sun.spct, span = NULL)
#'
#' @family peaks and valleys functions
#'
peaks <- function(x, span, ignore_threshold, strict, na.rm, ...) UseMethod("peaks")

#' @describeIn peaks Default returning always NA.
#' @export
peaks.default <-
  function(x, span = NA, ignore_threshold = NA, strict = NA, na.rm = FALSE, ...) {
  warning("Method 'peaks' not implemented for objects of class ", class(x)[1])
  x[NA]
}

#' @describeIn peaks Default function usable on numeric vectors.
#' @export
peaks.numeric <-
  function(x, span = 5, ignore_threshold = NA, strict = TRUE, na.rm = FALSE, ...) {
  x[find_peaks(x = x, span = span, strict = strict, na.rm = na.rm)]
}

#' @describeIn peaks  Method for "data.frame" objects.
#'
#' @param var.name Name of column where to look for peaks.
#'
#' @export
#'
peaks.data.frame <-
  function(x, span = 5, ignore_threshold = 0, strict = TRUE, na.rm = FALSE, var.name, ...) {
    if (is.null(var.name)) {
      return(x[NA, ])
    }
    peaks.idx <- find_peaks(x[[var.name]],
                            span = span, ignore_threshold = ignore_threshold,
                            strict = strict)
    x[peaks.idx, ]
  }

#' @describeIn peaks  Method for "generic_spct" objects.
#'
#' @export
#'
peaks.generic_spct <-
  function(x, span = 5, ignore_threshold = 0, strict = TRUE, na.rm = FALSE, var.name = NULL, ...) {
  if (is.null(var.name)) {
  # find target variable
    var.name <- names(x)
    var.name <- subset(var.name, sapply(x, is.numeric))
    var.name <- setdiff(var.name, "w.length")
    if (length(var.name) > 1L) {
      warning("Multiple numeric data columns found, explicit argument to 'var.name' required.")
      return(x[NA, ])
    }
  }
  peaks.idx <- find_peaks(x[[var.name]],
                          span = span, ignore_threshold = ignore_threshold,
                          strict = strict)
  x[peaks.idx, ]
}

#' @describeIn peaks  Method for "source_spct" objects.
#'
#' @param unit.out character One of "energy" or "photon"
#'
#' @export
#'
#' @examples
#' peaks(sun.spct)
#'
peaks.source_spct <-
  function(x, span = 5, ignore_threshold = 0, strict = TRUE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {
    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.irrad"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.irrad"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    peaks.idx <- find_peaks(z[[col.name]],
                            span = span, ignore_threshold = ignore_threshold,
                            strict = strict,
                            na.rm = na.rm)
    z[peaks.idx, ]
  }

#' @describeIn peaks  Method for "response_spct" objects.
#'
#' @export
#'
peaks.response_spct <-
  function(x, span = 5, ignore_threshold = 0.0, strict = TRUE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {
    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.response"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.response"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    peaks.idx <- find_peaks(z[[col.name]],
                            span = span, ignore_threshold = ignore_threshold,
                            strict = strict,
                            na.rm = na.rm)
    z[peaks.idx, ]
  }

#' @describeIn peaks  Method for "filter_spct" objects.
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
peaks.filter_spct <-
  function(x, span = 5, ignore_threshold = 0, strict = TRUE,
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty", default = "transmittance"),
           ...) {
    if (filter.qty == "transmittance") {
      z <- A2T(x, "replace", FALSE)
      col.name <- "Tfr"
    } else if (filter.qty == "absorbance") {
      z <- T2A(x, "replace", FALSE)
      col.name <- "A"
    } else {
      stop("Unrecognized 'filter.qty': ", filter.qty)
    }
    peaks.idx <- find_peaks(z[[col.name]],
                            span = span, ignore_threshold = ignore_threshold,
                            strict = strict,
                            na.rm = na.rm)
    z[peaks.idx, ]
  }

#' @describeIn peaks  Method for "reflector_spct" objects.
#'
#' @export
#'
peaks.reflector_spct <- function(x, span = 5, ignore_threshold = 0, strict = TRUE,
                                 na.rm = FALSE,
                                 ...) {
  peaks.idx <- find_peaks(x[["Rfr"]],
                          span = span, ignore_threshold = ignore_threshold,
                          strict = strict,
                          na.rm = na.rm)
  x[peaks.idx, ]
}

#' @describeIn peaks  Method for "cps_spct" objects.
#'
#' @export
#'
peaks.cps_spct <- function(x, span = 5, ignore_threshold = 0, strict = TRUE,
                           na.rm = FALSE,
                           ...) {
  peaks.idx <- find_peaks(x[["cps"]],
                          span = span, ignore_threshold = ignore_threshold,
                          strict = strict,
                          na.rm = na.rm)
  x[peaks.idx, ]
}

#' @describeIn peaks  Method for "cps_spct" objects.
#'
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
peaks.generic_mspct <- function(x,
                                span = 5,
                                ignore_threshold = 0,
                                strict = TRUE,
                                na.rm = FALSE,
                                ...,
                                .parallel = FALSE,
                                .paropts = NULL) {
  msmsply(x,
          .fun = peaks,
          span = span,
          ignore_threshold = ignore_threshold,
          strict = strict,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
  }

# valleys -------------------------------------------------------------------

#' Valleys or local minima
#'
#' Function that returns a subset of an R object with observations corresponding
#' to local maxima.
#'
#' @param x an R object
#' @param ignore_threshold numeric value between 0.0 and 1.0 indicating the
#'   relative size compared to tallest peak threshold below which valleys will be
#'   ignored.
#' @param span a peak is defined as an element in a sequence which is greater
#'   than all other elements within a window of width span centered at that
#'   element. The default value is 3, meaning that a peak is bigger than both of
#'   its neighbors. Default: 3.
#' @param strict logical flag: if TRUE, an element must be strictly greater than
#'   all other values in its window to be considered a peak. Default: TRUE.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for peaks.
#' @param ... ignored
#'
#' @return A subset of \code{x} with rows corresponding to local minima.
#'
#' @examples
#' valleys(sun.spct, span = 50)
#'
#' @export
#'
#' @family peaks and valleys functions
#'
valleys <- function(x, span, ignore_threshold, strict, ...) UseMethod("valleys")

#' @describeIn valleys Default function usable on numeric vectors.
#' @export
valleys.default <- function(x, span, ignore_threshold, strict, ...) {
  x[NA]
}

#' @describeIn valleys Default returning always NA.
#' @export
valleys.default <-
  function(x, span = NA, ignore_threshold = NA, strict = NA, na.rm = FALSE, ...) {
    warning("Method 'valleys' not implemented for objects of class ", class(x)[1])
    x[NA]
  }

#' @describeIn valleys Default function usable on numeric vectors.
#' @export
valleys.numeric <-
  function(x, span = 5, ignore_threshold, strict = TRUE, na.rm = FALSE, ...) {
    x[find_peaks(x = -x, span = span, strict = strict, na.rm = na.rm)]
  }

#' @describeIn valleys  Method for "data.frame" objects.
#'
#' @param var.name Name of column where to look for peaks.
#'
#' @export
#'
valleys.data.frame <-
  function(x, span = 5, ignore_threshold = 0, strict = TRUE, na.rm = FALSE, var.name, ...) {
    if (is.null(var.name)) {
      return(x[NA, ])
    }
    peaks.idx <- find_peaks(-x[[var.name]],
                            span = span, ignore_threshold = ignore_threshold,
                            strict = strict)
    x[peaks.idx, ]
  }

#' @describeIn valleys  Method for "generic_spct" objects.
#'
#' @export
#'
valleys.generic_spct <-
  function(x, span = 5, ignore_threshold = 0, strict = TRUE, na.rm = FALSE, var.name = NULL, ...) {
    if (is.null(var.name)) {
      # find target variable
      var.name <- names(x)
      var.name <- subset(var.name, sapply(x, is.numeric))
      var.name <- setdiff(var.name, "w.length")
      if (length(var.name) > 1L) {
        warning("Multiple numeric data columns found, explicit argument to 'var.name' required.")
        return(x[NA, ])
      }
    }
    peaks.idx <- find_peaks(-x[[var.name]],
                            span = span, ignore_threshold = ignore_threshold,
                            strict = strict)
    x[peaks.idx, ]
  }

#' @describeIn valleys  Method for "source_spct" objects.
#'
#' @param unit.out character One of "energy" or "photon"
#'
#' @export
#'
#' @examples
#' valleys(sun.spct)
#'
valleys.source_spct <-
  function(x, span = 5, ignore_threshold = 0.0, strict = TRUE, na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {
    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.irrad"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.irrad"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    valleys.idx <- find_peaks(-z[[col.name]],
                          span = span, ignore_threshold = ignore_threshold,
                          strict = strict,
                          na.rm = na.rm)
    z[valleys.idx, ]
  }

#' @describeIn valleys  Method for "response_spct" objects.
#'
#' @export
#'
valleys.response_spct <-
  function(x, span = 5, ignore_threshold = 0.0, strict = TRUE, na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {
    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.response"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.response"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    valleys.idx <- find_peaks(-z[[col.name]],
                            span = span, ignore_threshold = ignore_threshold,
                            strict = strict,
                            na.rm = na.rm)
    z[valleys.idx, ]
  }

#' @describeIn valleys  Method for "filter_spct" objects.
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
valleys.filter_spct <-
  function(x, span = 5, ignore_threshold = 0, strict = TRUE, na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty", default = "transmittance"),
           ...) {
    if (filter.qty == "transmittance") {
      z <- A2T(x, "replace", FALSE)
      col.name <- "Tfr"
    } else if (filter.qty == "absorbance") {
      z <- T2A(x, "replace", FALSE)
      col.name <- "A"
    } else {
      stop("Unrecognized 'filter.qty': ", filter.qty)
    }
    valleys.idx <- find_peaks(-z[[col.name]],
                              span = span, ignore_threshold = ignore_threshold,
                              strict = strict,
                              na.rm = na.rm)
    z[valleys.idx, ]
  }

#' @describeIn valleys  Method for "reflector_spct".
#'
#' @export
#'
valleys.reflector_spct <-
  function(x, span = 5, ignore_threshold = 0, strict = TRUE, na.rm = FALSE, ...) {
  valleys.idx <- find_peaks(-x[["Rfr"]],
                          span = span, ignore_threshold = ignore_threshold,
                          strict = strict,
                          na.rm = na.rm)
  x[valleys.idx, ]
}

#' @describeIn valleys  Method for "cps_spct" objects.
#'
#' @export
#'
valleys.cps_spct <-
  function(x, span = 5, ignore_threshold = 0, strict = TRUE, na.rm = FALSE, ...) {
  valleys.idx <- find_peaks(-x[["cps"]],
                          span = span, ignore_threshold = ignore_threshold,
                          strict = strict,
                          na.rm = na.rm)
  x[valleys.idx, ]
}

#' @describeIn valleys  Method for "generic_mspct" objects.
#'
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
valleys.generic_mspct <- function(x,
                                  span = 5,
                                  ignore_threshold = 0,
                                  strict = TRUE,
                                  na.rm = FALSE,
                                  ...,
                                  .parallel = FALSE,
                                  .paropts = NULL) {
  msmsply(x,
          .fun = valleys,
          span = span,
          ignore_threshold = ignore_threshold,
          strict = strict,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

# find wavelengths for a target y ----------------------------------------------

#' Find wavelength values in a spectrum
#'
#' Find wavelength values corresponding to a target y value in any spectrum. The
#' name of the column of the spectral data to be used to match the target needs
#' to be passed as argument unless the spectrum contains a single numerical
#' variable in addition to "w.length".
#'
#' @param x an R object
#' @param target numeric value indicating the spectral quantity value for which
#'   wavelengths are to be searched and interpolated if need. The character
#'   strings "half.maximum" and "half.range" are also accepted as arguments.
#' @param col.name.x character The name of the column in which to the
#'   independent variable is stored. Defaults to "w.length" for objects of
#'   class \code{"generic_spct"} or derived.
#' @param col.name character The name of the column in which to
#'    search for the target value.
#' @param .fun function A binary comparison function or operator.
#' @param interpolate logical Indicating whether the nearest wavelength value
#'   in \code{x} should be returned or a value calculated by linear
#'   interpolation between wavelength values stradling the target.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for the target.
#'
#' @note This function is used internally by method \code{wls_at_target()}, and
#'   these methods should be preferred in user code and scripts.
#'
#' @return A spectrum object of the same class as \code{x} with fewer rows,
#'   possibly even no rows. If \code{FALSE} is passed to \code{interpolate} a
#'   subset of \code{x} is returned, otherwise a new object of the same class
#'   containing interpolated wavelenths for the \code{target} value is
#'   returned.
#'
#' @export
#'
find_wls <- function(x,
                     target = NULL,
                     col.name.x = NULL,
                     col.name = NULL,
                     .fun = `<=`,
                     interpolate = FALSE,
                     na.rm = FALSE) {
  stopifnot(is.data.frame(x))
  x.class <- class(x)[1]
  if (is.null(target) || is.na(target)) {
    return(x[NULL, ])
  }
  if (is.null(col.name.x)) {
    if (is.any_spct(x)) {
      col.name.x <- "w.length"
    } else {
      warning("Object is not a \"generic spectrum\" explicit argument to 'col.name' required.")
      return(x[NULL, ])
    }
  }
  if (is.null(col.name)) {
    # find target variable
    col.name <- names(x)
    col.name <- subset(col.name, sapply(x, is.numeric))
    col.name <- setdiff(col.name, col.name.x)
    if (length(col.name) > 1L) {
      warning("Multiple numeric data columns found, explicit argument to 'col.name' required.")
      return(x[NULL, ])
    }
  }
  if (na.rm) {
    x <- na.omit(x)
  }
  if (is.character(target)) {
    if (target %in% c("half.maximum", "HM")) {
      target <- max(x[[col.name]]) / 2
    } else if (target %in% c("half.range", "HR")) {
      target <- mean(range(x[[col.name]]))
    } else {
      warning("Unrecognized character string: '", target, "' passed to 'target'", sep = "")
      target <- NA_real_
    }
    if (is.na(target)) {
      return(x[NULL, ])
    }
  }
  # test all rows for the condition
  true.rows <- .fun(x[[col.name]], target)
  # use run length to find transition points
  runs <- rle(true.rows)
  if (length(runs$lengths) < 2) {
    return(do.call(x.class, args = list()))
  }
  # accumulate run lengths to get index positions
  opening.idx <- cumsum(runs$lengths[-length(runs$lengths)])
  closing.idx <- opening.idx + 1L
  if (max(closing.idx) > nrow(x)) {
    closing.idx[length(closing.idx)] <- nrow(x)
  }
  if (interpolate) {
    # do vectorized interpolation to fetch true intersects
    delta.wl <- x[[col.name.x]][closing.idx] - x[[col.name.x]][opening.idx]
    delta.col <- x[[col.name]][closing.idx] - x[[col.name]][opening.idx]
    delta.col.target <- target - x[[col.name]][opening.idx]
    wl.increment <- delta.wl * abs(delta.col.target / delta.col)
    wls <- x[[col.name.x]][opening.idx] + wl.increment

    # return as a "short" spectrum containing only matching wls and target values
    z <- tibble::tibble(wls, target)
    names(z) <- c(col.name.x, col.name)
    if (x.class %in% spct_classes()) {
      z <- do.call(paste("as", x.class, sep = "."), args = list(x = z))
      # we need to copy our private attributes as we are building a new object
      z <- copy_attributes(x, z)
    }
  } else {
    # extract nearest wl value for target
    idxs <- ifelse(abs((x[[col.name]][closing.idx] - target) /
                     (x[[col.name]][closing.idx] - x[[col.name]][opening.idx])) > 0.5,
                   opening.idx,
                   closing.idx)
    # if the target value is close to a peak or valley, we may pick the same idx on both sides of it.
    z <- x[unique(idxs), ]
  }
  z
}

# find wavelengths for a target y ----------------------------------------------

#' Find wavelengths values corresponding to a target spectral value
#'
#' Find wavelength values corresponding to a target spectral value in a spectrum.
#' The name of the column of the spectral data to be used is inferred from the
#' class of \code{x} and the argument passed to \code{unit.out} or
#' \code{filter.qty} or their defaults that depend on R options set.
#'
#' @param x data.frame or spectrum object.
#' @param target numeric value indicating the spectral quantity value for which
#'   wavelengths are to be searched and interpolated if need. The character
#'   string "half.maximum" is also accepted as argument.
#' @param interpolate logical Indicating whether the nearest wavelength value
#'   in \code{x} should be returned or a value calculated by linear
#'   interpolation between wavelength values stradling the target.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for the target.
#' @param ... currently ignored.
#'
#' @return A spectrum object of the same class as \code{x} with fewer rows,
#'   possibly even no rows. If \code{FALSE} is passed to \code{interpolate} a
#'   subset of \code{x} is returned, otherwise a new object of the same class
#'   containing interpolated wavelenths for the \code{target} value is
#'   returned.
#'
#' @note When interpolation is used, only column \code{w.length} and the column
#'   against which the target value was compared are included in the returned
#'   object, otherwise, all columns in \code{x} are returned. We implement
#'   support for \code{data.frame} to simplify the coding of 'ggplot2' stats
#'   using this function.
#'
#' @examples
#' wls_at_target(sun.spct, target = 0.1)
#'
#' @export
#'
#' @family peaks and valleys functions
#'
wls_at_target <- function(x,
                          target = NULL,
                          interpolate = FALSE,
                          na.rm = FALSE,
                          ...) UseMethod("wls_at_target")

#' @describeIn wls_at_target Default returning always an empty object of the
#'   same class as \code{x}.
#' @export
#'
wls_at_target.default <-
  function(x,
           target = NULL,
           interpolate = FALSE,
           na.rm = FALSE,
           ...) {
    warning("Method 'wls_at_target' not implemented for objects of class ", class(x)[1])
    x[NULL]
  }

#' @describeIn wls_at_target Method for "generic_spct" objects.
#'
#' @param col.name character The name of the column in which to search for the
#'   target value.
#'
#' @export
#'
wls_at_target.generic_spct <-
  function(x,
           target = "half.maximum",
           interpolate = FALSE,
           na.rm = FALSE,
           col.name = NULL,
           ...) {
    find_wls(x,
             target = target,
             col.name = col.name,
             interpolate = interpolate,
             na.rm = na.rm,
             ...)
  }

#' @describeIn wls_at_target Method for "source_spct" objects.
#'
#' @param unit.out character One of "energy" or "photon"
#'
#' @export
#'
wls_at_target.source_spct <-
  function(x,
           target = "half.maximum",
           interpolate = FALSE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {
    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.irrad"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.irrad"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    find_wls(x,
             target = target,
             col.name = col.name,
             interpolate = interpolate,
             na.rm = na.rm)
  }

#' @describeIn wls_at_target Method for "response_spct" objects.
#' @export
#'
wls_at_target.response_spct <-
  function(x,
           target = "half.maximum",
           interpolate = FALSE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {
    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.response"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.response"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    find_wls(x,
             target = target,
             col.name = col.name,
             interpolate = interpolate,
             na.rm = na.rm)
  }

#' @describeIn wls_at_target Method for "filter_spct" objects.
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
wls_at_target.filter_spct <-
  function(x,
           target = "half.maximum",
           interpolate = FALSE,
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty", default = "transmittance"),
           ...) {
    if (filter.qty == "transmittance") {
      z <- A2T(x, "replace", FALSE)
      col.name <- "Tfr"
    } else if (filter.qty == "absorbance") {
      z <- T2A(x, "replace", FALSE)
      col.name <- "A"
    } else {
      stop("Unrecognized 'filter.qty': ", filter.qty)
    }
    find_wls(x,
             target = target,
             col.name = col.name,
             interpolate = interpolate,
             na.rm = na.rm)
  }

#' @describeIn wls_at_target Method for "reflector_spct" objects.
#' @export
#'
wls_at_target.reflector_spct <-
  function(x,
           target = "half.maximum",
           interpolate = FALSE,
           na.rm = FALSE,
           ...) {
    find_wls(x,
             target = target,
             col.name = "Rfr",
             interpolate = interpolate,
             na.rm = na.rm)
  }

#' @describeIn wls_at_target Method for "cps_spct" objects.
#'
#' @export
#'
wls_at_target.cps_spct <-
  function(x,
           target = "half.maximum",
           interpolate = FALSE,
           na.rm = FALSE,
           ...) {
    find_wls(x,
             target = target,
             col.name = "cps",
             interpolate = interpolate,
             na.rm = na.rm)
  }

#' @describeIn wls_at_target  Method for "generic_mspct" objects.
#'
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
wls_at_target.generic_mspct <- function(x,
                                        target = "half.maximum",
                                        interpolate = FALSE,
                                        na.rm = FALSE,
                                        ...,
                                        .parallel = FALSE,
                                        .paropts = NULL) {
  msmsply(x,
          .fun = wls_at_target,
          target = target,
          interpolate = interpolate,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

