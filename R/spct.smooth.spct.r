#' Smooth a spectrum
#'
#' These functions implement one original methods and acts as a wrapper for
#' other common R smoothing functions. The advantage of using this function for
#' smoothing spectral objects is that it simplifies the user interface and sets,
#' when needed, defaults suitable for spectral data.
#'
#' @param x an R object.
#' @param method a character string "custom", "lowess", "supsmu".
#' @param strength numeric value to adjust the degree of smoothing. Ignored if
#'   method-specific parameters are passed through \code{...}.
#' @param wl.range any R object on which applying the method \code{range()}
#'   yields a vector of two numeric values, describing a range of wavelengths
#'   (nm). \code{NA} is interpreted as data's own min or max value.
#' @param na.rm	logical A flag indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param ... other parameters passed to the underlying smoothing functions.
#'
#' @return A copy of \code{x} with spectral data values replaced by smoothed
#'   ones.
#'
#' @note Method "custom" is our home-brewed method which applies strong
#'   smoothing to low signal regions of the spectral data, and weaker or no
#'   smoothing to the high signal areas. Values very close to zero are set to
#'   zero with a limit which depends on the local variation. This method is an
#'   ad-hock method suitable for smoothing spectral data obtained with
#'   spectrometers. In the cased of methods "lowess" and "supsmu" the current
#'   function behaves like a wrapper of the functions of the same names from
#'   base R.
#'
#' @export
#'
#' @examples
#'
#' my.spct <- clip_wl(sun.spct, c(400, 500))
#' smooth_spct(my.spct)
#' smooth_spct(my.spct, method = "custom", strength = 1)
#' smooth_spct(my.spct, method = "custom", strength = 4)
#' smooth_spct(my.spct, method = "supsmu", strength = 4)
#'
smooth_spct <- function(x,
                        method,
                        strength,
                        wl.range,
                        ...) UseMethod("smooth_spct")

#' @describeIn smooth_spct Default for generic function
#'
#' @export
#'
smooth_spct.default <- function(x, method, strength, wl.range, ...) {
  warning("'smooth_spct' is not defined for objects of class ", class(x)[1])
  return(x)
}

#' @describeIn smooth_spct Smooth a source spectrum
#'
#' @export
#'
smooth_spct.source_spct <- function(x,
                                    method = "custom",
                                    strength = 1,
                                    wl.range = NULL,
                                    na.rm = FALSE,
                                    ...) {
  num.spectra <- getMultipleWl(x)
  if (num.spectra != 1) {
    warning("Skipping smoothing as object contains ",
            num.spectra, " spectra in long form.")
    return(x)
  }

  stopifnot(strength >= 0)
  if (strength == 0) {
    return(x)
  }

  # if both energy and photon data present, we retain only energy and set a flag
  e.and.q.input <- all(c("s.e.irrad", "s.q.irrad") %in% names(x))
  if (e.and.q.input) {
    x <- q2e(x, action = "replace")
  }

  if ("s.e.irrad" %in% names(x) && anyNA(x[["s.e.irrad"]])) {
    if (na.rm) {
      message("Removing NA values at ", sum(is.na(x[["s.e.irrad"]])), " wavelengths.")
      x <- na.omit(x)
    } else {
      stop("NAs encountered when smoothing.")
    }
  } else if ("s.q.irrad" %in% names(x) && anyNA(x[["s.q.irrad"]])) {
    if (na.rm) {
      message("Removing NA values at ", sum(is.na(x[["s.q.irrad"]])), " wavelengths.")
      x <- na.omit(x)
    } else {
      stop("NAs encountered when smoothing.")
    }
  }

  # Skip checks for intermediate results
  # as intermediate values may be off-range
  prev_state <- disable_check_spct()
  on.exit(set_check_spct(prev_state), add = TRUE)

  # this needs to be after NAs are omitted
  if (!is.null(wl.range)) {
    wl.range <- range(wl.range)
    wl.selector <-
      x[["w.length"]] >= wl.range[1] & x[["w.length"]] <= wl.range[2]
  } else {
    wl.selector <- TRUE
  }

  xx <- x[wl.selector, ]

  if (method == "lowess") {
    span = 1/50 * strength
    if ("s.e.irrad" %in% names(x)) {
      z <- stats::lowess(xx[["w.length"]], xx[["s.e.irrad"]], f = span, ...)
      x[wl.selector, "s.e.irrad"] <- z[["y"]]
    } else if ("s.q.irrad" %in% names(x)) {
      z <- stats::lowess(xx[["w.length"]], xx[["s.q.irrad"]], f = span, ...)
      x[wl.selector, "s.q.irrad"] <- z[["y"]]
    }
    comment.text <- paste("Smoothed using 'lowess', f =", signif(span, 3))
  } else if (method == "supsmu") {
    span = 1/50 * strength
    if ("s.e.irrad" %in% names(x)) {
      z <- stats::supsmu(xx[["w.length"]], xx[["s.e.irrad"]], span = span, ...)
      x[wl.selector, "s.e.irrad"] <- z[["y"]]
    } else if ("s.q.irrad" %in% names(x)) {
      z <- stats::supsmu(xx[["w.length"]], xx[["s.q.irrad"]], span = span, ...)
      x[wl.selector, "s.q.irrad"] <- z[["y"]]
    }
    comment.text <-  paste("Smoothed using 'supsmu', span =", signif(span, 3))
  } else if (method == "custom") {
    smooth.limit <- 1e-3 * strength
    if ("s.e.irrad" %in% names(x)) {
      zero.limit.cnst <- max(x[["s.e.irrad"]]) * 3e-4 * strength
      smooth.threshold <-  max(x[["s.e.irrad"]]) * 5e-2 * strength
      z <- adaptive_smoothing(xx[["w.length"]], xx[["s.e.irrad"]],
                              zero.limit.cnst = zero.limit.cnst,
                              smooth.limit = smooth.limit,
                              smooth.threshold = smooth.threshold,
                                     ...)
      x[wl.selector, "s.e.irrad"] <- z[["y"]]
    } else if ("s.q.irrad" %in% names(x)) {
      zero.limit.cnst <- max(x[["s.q.irrad"]]) * 3e-4 * strength
      smooth.threshold <-  max(x[["s.q.irrad"]]) * 5e-2 * strength
      z <- adaptive_smoothing(xx[["w.length"]], xx[["s.q.irrad"]],
                              zero.limit.cnst = zero.limit.cnst,
                              smooth.limit = smooth.limit,
                              smooth.threshold = smooth.threshold,
                                     ...)
      x[wl.selector, "s.q.irrad"] <- z[["y"]]
    }
    comment.text <- paste("Smoothed using 'custom', smooth.limit =",
                          signif(smooth.limit, 3))
  }

#  out.spct <- copy_attributes(x, out.spct)
  # restore s.q.irrad if needed
  if (e.and.q.input) {
    e2q(x, action = "add", byref = TRUE)
  }

  if (!is.null(comment(x))) {
    comment(x) <- paste(comment.text, "\n\n", comment(x))
  } else {
    comment(x) <- comment.text
  }
  check_spct(x, force = FALSE)
}

#' @describeIn smooth_spct Smooth a filter spectrum
#'
#' @export
#'
smooth_spct.filter_spct <- function(x,
                                    method = "custom",
                                    strength = 1,
                                    wl.range = NULL,
                                    na.rm = FALSE,
                                    ...) {
  num.spectra <- getMultipleWl(x)
  if (num.spectra != 1) {
    warning("Skipping smoothing as object contains ",
            num.spectra, " spectra in long form.")
    return(x)
  }

  stopifnot(strength >= 0)
  if (strength == 0) {
    return(x)
  }

  # if both Tfr and A data present, we retain only Tfr and set a flag
  T.and.A.input <- all(c("Tfr", "A") %in% names(x))

  if ("Tfr" %in% names(x) && any(c("Afr", "A") %in% names(x))) {
    warning("Retaing \"Tfr\", droping \"Afr\" and or  \"A\".")
    x[ , setdiff(names(x), c("Afr", "A"))]
  }
  if ("A" %in% names(x) && "Afr" %in% names(x)) {
    warning("Retaing \"A\", droping \"Afr\".")
    x[ , setdiff(names(x), "Afr")]
  }

  if ("Tfr" %in% names(x) && anyNA(x[["Tfr"]])) {
    if (na.rm) {
      message("Removing NA values at ", sum(is.na(x[["Tfr"]])), " wavelengths.")
      x <- na.omit(x)
    } else {
      stop("NAs encountered when smoothing.")
    }
  } else if ("A" %in% names(x) && anyNA(x[["A"]])) {
    if (na.rm) {
      message("Removing NA values at ", sum(is.na(x[["A"]])), " wavelengths.")
      x <- na.omit(x)
    } else {
      stop("NAs encountered when smoothing.")
    }
  } else if ("Afr" %in% names(x) && anyNA(x[["A"]])) {
    if (na.rm) {
      message("Removing NA values at ", length(is.na(x[["Afr"]])), " wavelengths.")
      x <- na.omit(x)
    } else {
      stop("NAs encountered when smoothing.")
    }
  }

  # Skip checks for intermediate results
  # as intermediate values may be off-range
  prev_state <- disable_check_spct()
  on.exit(set_check_spct(prev_state), add = TRUE)

  # this needs to be after NAs are omitted
  if (!is.null(wl.range)) {
    wl.range <- range(wl.range)
    wl.selector <-
      x[["w.length"]] >= wl.range[1] & x[["w.length"]] <= wl.range[2]
  } else {
    wl.selector <- TRUE
  }

  xx <- x[wl.selector, ]

  if (method == "lowess") {
    span = 1/50 * strength
    if ("Tfr" %in% names(x)) {
      z <- stats::lowess(xx[["w.length"]], xx[["Tfr"]], f = span, ...)
      x[wl.selector, "Tfr"] <- z[["y"]]
    } else if ("A" %in% names(x)) {
      z <- stats::lowess(xx[["w.length"]], xx[["A"]], f = span, ...)
      x[wl.selector, "A"] <- z[["y"]]
    } else if ("Afr" %in% names(x)) {
      z <- stats::lowess(xx[["w.length"]], xx[["Afr"]], f = span, ...)
      x[wl.selector, "Afr"] <- z[["y"]]
    }
    comment.text <- paste("Smoothed using 'lowess', f =", signif(span, 3))
  } else if (method == "supsmu") {
    span = 1/50 * strength
    if ("Tfr" %in% names(x)) {
      z <- stats::supsmu(xx[["w.length"]], xx[["Tfr"]], span = span, ...)
      x[wl.selector, "Tfr"] <- z[["y"]]
    } else if ("A" %in% names(x)) {
      z <- stats::supsmu(xx[["w.length"]], xx[["A"]], span = span, ...)
      x[wl.selector, "A"] <- z[["y"]]
    } else if ("Afr" %in% names(x)) {
      z <- stats::supsmu(xx[["w.length"]], xx[["Afr"]], span = span, ...)
      x[wl.selector, "Afr"] <- z[["y"]]
    }
    comment.text <-  paste("Smoothed using 'supsmu', span =", signif(span, 3))
  } else if (method == "custom") {
    smooth.limit <- 1e-3 * strength
    if ("Tfr" %in% names(x)) {
      zero.limit.cnst <- max(x[["Tfr"]]) * 3e-4 * strength
      smooth.threshold <-  max(x[["Tfr"]]) * 5e-2 * strength
      z <- adaptive_smoothing(xx[["w.length"]], xx[["Tfr"]],
                              zero.limit.cnst = zero.limit.cnst,
                              smooth.limit = smooth.limit,
                              smooth.threshold = smooth.threshold,
                              ...)
      x[wl.selector, "Tfr"] <- z[["y"]]
    } else if ("A" %in% names(x)) {
      zero.limit.cnst <- max(x[["A"]]) * 3e-4 * strength
      smooth.threshold <-  max(x[["A"]]) * 5e-2 * strength
      z <- adaptive_smoothing(xx[["w.length"]], xx[["A"]],
                              zero.limit.cnst = zero.limit.cnst,
                              smooth.limit = smooth.limit,
                              smooth.threshold = smooth.threshold,
                              ...)
      x[wl.selector, "A"] <- z[["y"]]
    } else if ("Afr" %in% names(x)) {
      zero.limit.cnst <- max(x[["Afr"]]) * 3e-4 * strength
      smooth.threshold <-  max(x[["Afr"]]) * 5e-2 * strength
      z <- adaptive_smoothing(xx[["w.length"]], xx[["Afr"]],
                              zero.limit.cnst = zero.limit.cnst,
                              smooth.limit = smooth.limit,
                              smooth.threshold = smooth.threshold,
                              ...)
      x[wl.selector, "Afr"] <- z[["y"]]
    }
    comment.text <- paste("Smoothed using 'custom', smooth.limit =",
                          signif(smooth.limit, 3))
  }

  #  out.spct <- copy_attributes(x, out.spct)
  # restore s.q.irrad if needed
  if (T.and.A.input) {
    T2A(x, action = "add", byref = TRUE)
  }

  if (!is.null(comment(x))) {
    comment(x) <- paste(comment.text, "\n\n", comment(x))
  } else {
    comment(x) <- comment.text
  }
  check_spct(x, force = FALSE)
}

#' @describeIn smooth_spct Smooth a reflector spectrum
#'
#' @export
#'
smooth_spct.reflector_spct <- function(x,
                                       method = "custom",
                                       strength = 1,
                                       wl.range = NULL,
                                       na.rm = FALSE,
                                       ...) {
  num.spectra <- getMultipleWl(x)
  if (num.spectra != 1) {
    warning("Skipping smoothing as object contains ",
            num.spectra, " spectra in long form.")
    return(x)
  }

  stopifnot(strength >= 0)
  if (strength == 0) {
    return(x)
  }

  if (anyNA(x[["Rfr"]])) {
    if (na.rm) {
      message("Removing NA values at ", sum(is.na(x[["Rfr"]])), " wavelengths.")
      x <- na.omit(x)
    } else {
      stop("NAs encountered when smoothing.")
    }
  }

  # Skip checks for intermediate results
  # as intermediate values may be off-range
  prev_state <- disable_check_spct()
  on.exit(set_check_spct(prev_state), add = TRUE)

  # this needs to be after NAs are omitted
  if (!is.null(wl.range)) {
    wl.range <- range(wl.range)
    wl.selector <-
      x[["w.length"]] >= wl.range[1] & x[["w.length"]] <= wl.range[2]
  } else {
    wl.selector <- TRUE
  }

  xx <- x[wl.selector, ]

  if (method == "lowess") {
    span = 1/50 * strength
    z <- stats::lowess(xx[["w.length"]], xx[["Rfr"]], f = span, ...)
    x[wl.selector, "Rfr"] <- z[["y"]]
    comment.text <- paste("Smoothed using 'lowess', f =", signif(span, 3))
  } else if (method == "supsmu") {
    span = 1/50 * strength
    z <- stats::supsmu(xx[["w.length"]], xx[["Rfr"]], span = span, ...)
    x[wl.selector, "Rfr"] <- z[["y"]]
    comment.text <-  paste("Smoothed using 'supsmu', span =", signif(span, 3))
  } else if (method == "custom") {
    smooth.limit <- 1e-3 * strength
    zero.limit.cnst <- max(x[["Rfr"]]) * 3e-4 * strength
    smooth.threshold <-  max(x[["Rfr"]]) * 5e-2 * strength
    z <- adaptive_smoothing(xx[["w.length"]], xx[["Rfr"]],
                            zero.limit.cnst = zero.limit.cnst,
                            smooth.limit = smooth.limit,
                            smooth.threshold = smooth.threshold,
                            ...)
    x[wl.selector, "Rfr"] <- z[["y"]]
    comment.text <- paste("Smoothed using 'custom', smooth.limit =",
                          signif(smooth.limit, 3))
  }

  if (!is.null(comment(x))) {
    comment(x) <- paste(comment.text, "\n\n", comment(x))
  } else {
    comment(x) <- comment.text
  }
  check_spct(x, force = FALSE)
}

#' @describeIn smooth_spct Smooth a response spectrum
#'
#' @export
#'
smooth_spct.response_spct <- function(x,
                                      method = "custom",
                                      strength = 1,
                                      wl.range = NULL,
                                      na.rm = FALSE,
                                      ...) {
  num.spectra <- getMultipleWl(x)
  if (num.spectra != 1) {
    warning("Skipping smoothing as object contains ",
            num.spectra, " spectra in long form.")
    return(x)
  }

  stopifnot(strength >= 0)
  if (strength == 0) {
    return(x)
  }

  # if both energy and photon data present, we retain only energy and set a flag
  e.and.q.input <- all(c("s.e.response", "s.q.response") %in% names(x))
  if (e.and.q.input) {
    x <- q2e(x, action = "replace")
  }

  if ("s.e.response" %in% names(x) && anyNA(x[["s.e.response"]])) {
    if (na.rm) {
      message("Removing NA values at ", sum(is.na(x[["s.e.response"]])), " wavelengths.")
      x <- na.omit(x)
    } else {
      stop("NAs encountered when smoothing.")
    }
  } else if ("s.q.response" %in% names(x) && anyNA(x[["s.q.response"]])) {
    if (na.rm) {
      message("Removing NA values at ", sum(is.na(x[["s.q.response"]])), " wavelengths.")
      x <- na.omit(x)
    } else {
      stop("NAs encountered when smoothing.")
    }
  }

  # Skip checks for intermediate results
  # as intermediate values may be off-range
  prev_state <- disable_check_spct()
  on.exit(set_check_spct(prev_state), add = TRUE)

  # this needs to be after NAs are omitted
  if (!is.null(wl.range)) {
    wl.range <- range(wl.range)
    wl.selector <-
      x[["w.length"]] >= wl.range[1] & x[["w.length"]] <= wl.range[2]
  } else {
    wl.selector <- TRUE
  }

  xx <- x[wl.selector, ]

  if (method == "lowess") {
    span = 1/50 * strength
    if ("s.e.response" %in% names(x)) {
      z <- stats::lowess(xx[["w.length"]], xx[["s.e.response"]], f = span, ...)
      x[wl.selector, "s.e.response"] <- z[["y"]]
    } else if ("s.q.response" %in% names(x)) {
      z <- stats::lowess(xx[["w.length"]], xx[["s.q.response"]], f = span, ...)
      x[wl.selector, "s.q.response"] <- z[["y"]]
    }
    comment.text <- paste("Smoothed using 'lowess', f =", signif(span, 3))
  } else if (method == "supsmu") {
    span = 1/50 * strength
    if ("s.e.response" %in% names(x)) {
      z <- stats::supsmu(xx[["w.length"]], xx[["s.e.response"]], span = span, ...)
      x[wl.selector, "s.e.response"] <- z[["y"]]
    } else if ("s.q.response" %in% names(x)) {
      z <- stats::supsmu(xx[["w.length"]], xx[["s.q.response"]], span = span, ...)
      x[wl.selector, "s.q.response"] <- z[["y"]]
    }
    comment.text <-  paste("Smoothed using 'supsmu', span =", signif(span, 3))
  } else if (method == "custom") {
    smooth.limit <- 1e-3 * strength
    if ("s.e.response" %in% names(x)) {
      zero.limit.cnst <- max(x[["s.e.response"]]) * 3e-4 * strength
      smooth.threshold <-  max(x[["s.e.response"]]) * 5e-2 * strength
      z <- adaptive_smoothing(xx[["w.length"]], xx[["s.e.response"]],
                              zero.limit.cnst = zero.limit.cnst,
                              smooth.limit = smooth.limit,
                              smooth.threshold = smooth.threshold,
                              ...)
      x[wl.selector, "s.e.response"] <- z[["y"]]
    } else if ("s.q.response" %in% names(x)) {
      zero.limit.cnst <- max(x[["s.q.response"]]) * 3e-4 * strength
      smooth.threshold <-  max(x[["s.q.response"]]) * 5e-2 * strength
      z <- adaptive_smoothing(xx[["w.length"]], xx[["s.q.response"]],
                              zero.limit.cnst = zero.limit.cnst,
                              smooth.limit = smooth.limit,
                              smooth.threshold = smooth.threshold,
                              ...)
      x[wl.selector, "s.q.response"] <- z[["y"]]
    }
    comment.text <- paste("Smoothed using 'custom', smooth.limit =",
                          signif(smooth.limit, 3))
  }

  #  out.spct <- copy_attributes(x, out.spct)
  # restore s.q.response if needed
  if (e.and.q.input) {
    e2q(x, action = "add", byref = TRUE)
  }

  if (!is.null(comment(x))) {
    comment(x) <- paste(comment.text, "\n\n", comment(x))
  } else {
    comment(x) <- comment.text
  }
  check_spct(x, force = FALSE)
}


#' @describeIn smooth_spct
#'
#' @export
#'
smooth_spct.generic_mspct <-
  function(x,
           method = "custom",
           strength = 1,
           wl.range = NULL,
           na.rm = FALSE,
           ...) {
    msmsply(x,
            smooth_spct,
            method = method,
            strength = strength,
            wl.range = wl.range,
            na.rm = na.rm,
            ...)
  }

#' Custom smoothing
#'
#' Generic implementation of the custom smoothing method
#'
#' @param x,y numeric vectors of equal length.
#' @param zero.limit numeric vector of length one or of the same leangth as
#'   \code{x} and \code{y}. Smaller values in y are forced to zero.
#' @param smooth.limt numeric mad/median value above which no smoothing is
#'   applied.
#' @param smooth.thershold numeric y value above which no smoothing is applied.
#'
#' @keywords internal
#'
adaptive_smoothing <- function(x, y,
                               zero.limit.cnst = max(y) * 1e-4,
                               smooth.limit = 1e-3,
                               smooth.threshold =  max(y) * 1e-1,
                               ...) {
  stopifnot(!anyNA(x) && !anyNA(y))
  max_response <- max(y)
  # this could be tweaked in many ways...
  if (length(zero.limit.cnst == 1L)) {
    zero_limit <- (zero.limit.cnst * 600) / x # stronger zero cleaning at smaller x
  } else {
    zero_limit <- zero.limit.cnst
  }
  runmad7 <- zoo::rollapply(y, width = 7, FUN = stats::mad, partial = TRUE)
  runmed3 <- stats::runmed(y, 3, endrule = "median")
  runmed7 <- stats::runmed(y, 7, endrule = "median")
  runmed19 <- stats::runmed(y, 19, endrule = "median")
  runmin5 <- zoo::rollapply(y, width = 5, FUN = min, partial = TRUE)
  # runmadmed is similar to coefficient of variation but distribution free
  # we need to avoid division by 0.0 and we use zero_limit / 10 as close enough to zero
  runmadmed <- ifelse(runmad7 < zero.limit.cnst * 1e-1 | runmed7 < zero_limit * 1e-1,
                                         0.0, runmad7/abs(runmed7))
  z <-
    ifelse( (runmed19 < zero_limit) | (runmin5 < zero_limit * 5e-2) | (y < 0),
            0,
            ifelse((y > smooth.threshold) | (runmadmed < smooth.limit),
                   y,
                   ifelse(runmadmed < 2 * smooth.limit,
                          runmed3,
                          ifelse(runmadmed < 4 * smooth.limit,
                                 runmed7, runmed19))))
  response.good <- runmadmed / runmed19 * max(runmed19) < 1.0
  if (anyNA(z)) {
    warning(sum(is.na(z)), " NAs introduced during smoothing")
  }
  num_bad <- sum(!response.good, na.rm=TRUE)
  if (num_bad > length(x) / 20) {
    message(num_bad, " possibly 'bad' values in smoothed spectral response")
  }
  # returned value similar to that of R's lowess() and supsmu()
  list(x = x, y = z)
}
