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
      message("Removing NA values at ", length(is.na(x[["s.e.irrad"]])), " wavelengths.")
      x <- na.omit(x)
    } else {
      stop("NAs encountered when smoothing.")
    }
  } else if ("s.q.irrad" %in% names(x) && anyNA(x[["s.q.irrad"]])) {
    if (na.rm) {
      message("Removing NA values at ", length(is.na(x[["s.q.irrad"]])), " wavelengths.")
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
    zero.limit.cnst <- max(x[["s.e.irrad"]]) * 3e-4 / strength
    smooth.limit <- 1e-3
    smooth.threshold <-  max(x[["s.e.irrad"]]) * 5e-2 / strength
    if ("s.e.irrad" %in% names(x)) {
      z <- adaptive_smoothing(xx[["w.length"]], xx[["s.e.irrad"]],
                              zero.limit.cnst = zero.limit.cnst,
                              smooth.limit = smooth.limit,
                              smooth.threshold = smooth.threshold,
                                     ...)
      x[wl.selector, "s.e.irrad"] <- z[["y"]]
    } else if ("s.q.irrad" %in% names(x)) {
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
            num.spectra, " spectra")
    return(x)
  }
  T.and.A.input <- all(c("Tfr", "A") %in% names(x))

  if ("Tfr" %in% names(x) && anyNA(x[["Tfr"]])) {
    if (na.rm) {
      message("Removing NA values at ", length(is.na(x[["Tfr"]])), " wavelengths.")
      x <- na.omit(A2T(x, action = "replace"))
    } else {
      warning("NAs encountered when smoothing, returning input unchanged.")
      return(x)
    }
  } else if ("A" %in% names(x) && anyNA(x[["A"]])) {
    if (na.rm) {
      message("Removing NA values at ", length(is.na(x[["A"]])), " wavelengths.")
      x <- na.omit(T2A(x, action = "replace"))
    } else {
      warning("NAs encountered when smoothing, returning input unchanged.")
      return(x)
    }
  }

  # Skip checks for intermediate results
  # as intermediate values may be off-range
  prev_state <- disable_check_spct()
  on.exit(set_check_spct(prev_state), add = TRUE)
  if (method == "lowess") {
    span = 1/50 * strength
    if ("Tfr" %in% names(x)) {
      out.spct <- stats::lowess(x[["w.length"]], x[["Tfr"]], f = span, ...)
      names(out.spct) <- c("w.length", "Tfr")
    } else if ("A" %in% names(x)) {
      out.spct <- stats::lowess(x[["w.length"]], x[["A"]], f = span, ...)
      names(out.spct) <- c("w.length", "A")
    }
    setFilterSpct(out.spct, Tfr.type = attr(x, "Tfr.type", exact = TRUE))
    if (T.and.A.input) {
      T2A(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      comment.text <- paste("Smoothed using 'lowess', f =", signif(span, 3), "\n\n", comment(x))
    } else {
      comment.text <- paste("Smoothed using 'lowess', f =", signif(span, 3))
    }
  } else if (method == "supsmu") {
    span = 1/50 * strength
    if ("Tfr" %in% names(x)) {
      out.spct <- stats::supsmu(x[["w.length"]], x[["Tfr"]], span = span, ...)
      names(out.spct) <- c("w.length", "Tfr")
    } else if ("A" %in% names(x)) {
      out.spct <- stats::supsmu(x[["w.length"]], x[["A"]], span = span, ...)
      names(out.spct) <- c("w.length", "A")
    }
    setFilterSpct(out.spct, Tfr.type = attr(x, "Tfr.type", exact = TRUE))
    if (T.and.A.input) {
      T2A(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      comment.text <- paste("Smoothed using 'supsmu', span =", signif(span, 3), "\n\n", comment(x))
    } else {
      comment.text <-  paste("Smoothed using 'supsmu', span =", signif(span, 3))
    }
  } else if (method == "custom") {
    # my own and inefficient method!
    # as the spectrum is already in energy units, we need to normalize thresholds
    out.spct <- x # just to avoid editing the code
    A2T(out.spct, action = "replace", byref = TRUE)
    max_Tfr <- 1
    smoothing_coef <- 1
    smoothing_hi_lim <- max(out.spct[["w.length"]])
    # this could be tweaked in many ways...
    zero_limit_cnst <- max_Tfr * 3e-4 / strength
    out.spct[["zero_limit"]] <-  (zero_limit_cnst * 600) / out.spct[["w.length"]]
    smooth_limit <- 1e-3 * smoothing_coef # just a guess for runmadmed
    smooth_threshold <- 5e-2 * max_Tfr / strength # for Tfr
#    out.spct[["runmad"]] <- caTools::runmad(out.spct[["Tfr"]], 7, endrule="mad")
    out.spct[["runmad"]] <- zoo::rollapply(out.spct[["Tfr"]], width = 7, FUN = stats::mad, partial = TRUE)
    out.spct[["runmed3"]] <- stats::runmed(out.spct[["Tfr"]], 3, endrule="median")
    out.spct[["runmed7"]] <- stats::runmed(out.spct[["Tfr"]], 7, endrule="median")
    out.spct[["runmed19"]] <- stats::runmed(out.spct[["Tfr"]], 19, endrule="median")
#    out.spct[["runmin5"]] <- caTools::runmin(out.spct[["Tfr"]], 5)
    out.spct[["runmin5"]] <- zoo::rollapply(out.spct[["Tfr"]], width = 5, FUN = min, partial = TRUE)
    # we need to avoid division by 0.0 and we use zero_limit / 10 close enough to zero
    out.spct[["runmadmed"]] <- with(out.spct,
                                    ifelse(runmad < zero_limit_cnst * 1e-1 | runmed7 < zero_limit * 1e-1,
                                           0.0, runmad/abs(runmed7)))
    out.spct[["Tfr.sm"]] <- with(out.spct,
               ifelse( (runmed19 < zero_limit) | (runmin5 < zero_limit * 5e-2), 0.0,
                       ifelse((Tfr > smooth_threshold) | (runmadmed < smooth_limit), Tfr,
                              ifelse(runmadmed < 2 * smooth_limit, runmed3,
                                     ifelse(runmadmed < 4 * smooth_limit, runmed7, runmed19)))))
    out.spct[["Tfr"]] <- with(out.spct,
                              ifelse(w.length < smoothing_hi_lim, Tfr.sm, Tfr))
    out.spct[["Tfr.good"]] <- out.spct[["runmadmed"]] / out.spct[["runmed19"]] * max(out.spct[["runmed19"]]) < 1.0
    if (any(is.na(out.spct[["Tfr"]]))) {
      warning(sum(is.na(out.spct[["Tfr"]])), " NAs in returned spectral irradiance")
    }
    num_bad <- sum(!out.spct[["Tfr.good"]], na.rm=TRUE)
    if (num_bad > length(out.spct) / 20) {
      message(num_bad, " possibly 'bad' values in smoothed spectral Tfr")
    }
    out.spct <- out.spct[ , c("w.length", "Tfr")]
    setFilterSpct(out.spct, Tfr.type = attr(x, "Tfr.type", exact = TRUE))
    if (T.and.A.input) {
      T2A(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      comment.text <- paste("Smoothed using 'custom', smooth_limit =",
                                 signif(smooth_limit, 3), "\n\n", comment(x))
    } else {
      comment.text <- paste("Smoothed using 'custom', smooth_limit =", signif(smooth_limit, 3))
    }
  }
  out.spct <- copy_attributes(x, out.spct)
  comment(out.spct) <- comment.text
  check_spct(out.spct, force = FALSE)
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
            num.spectra, " spectra")
    return(x)
  }
  if ("Rfr" %in% names(x) && anyNA(x[["Rfr"]])) {
    if (na.rm) {
      message("Removing NA values at ", length(is.na(x[["Rfr"]])), " wavelengths.")
      x <- na.omit(x)
    } else {
      warning("NAs encountered when smoothing, returning input unchanged.")
      return(x)
    }
  }

  # Skip checks for intermediate results
  # as intermediate values may be off-range
  prev_state <- disable_check_spct()
  on.exit(set_check_spct(prev_state), add = TRUE)
  if (method == "lowess") {
    span = 1/50 * strength
    if ("Rfr" %in% names(x)) {
      out.spct <- stats::lowess(x[["w.length"]], x[["Rfr"]], f = span, ...)
      names(out.spct) <- c("w.length", "Rfr")
    }
    setReflectorSpct(out.spct)
    if (!is.null(comment(x))) {
      comment.text <- paste("Smoothed using 'lowess', f =", signif(span, 3), "\n\n", comment(x))
    } else {
      comment.text <- paste("Smoothed using 'lowess', f =", signif(span, 3))
    }
  } else if (method == "supsmu") {
    span = 1/50 * strength
    if ("Rfr" %in% names(x)) {
      out.spct <- stats::supsmu(x[["w.length"]], x[["Rfr"]], span = span, ...)
      names(out.spct) <- c("w.length", "Rfr")
    }
    setReflectorSpct(out.spct)
    if (!is.null(comment(x))) {
      comment.text <- paste("Smoothed using 'supsmu', span =", signif(span, 3), "\n\n", comment(x))
    } else {
      comment.text <- paste("Smoothed using 'supsmu', span =", signif(span, 3))
    }
  } else if (method == "custom") {
    # my own and inefficient method!
    # as the spectrum is already in energy units, we need to normalize thresholds
    out.spct <- x # we make a working copy
    max_Rfr <- 1
    smoothing_coef <- 1
    smoothing_hi_lim <- max(out.spct[["w.length"]])
    # this could be tweaked in many ways...
    zero_limit_cnst <- max_Rfr * 3e-4 / strength
    out.spct[["zero_limit"]] <-  (zero_limit_cnst * 600) / out.spct[["w.length"]]
    smooth_limit <- 1e-3 * smoothing_coef # just a guess for runmadmed
    smooth_threshold <- 5e-2 * max_Rfr / strength # for Rfr
#    out.spct[["runmad"]] <- caTools::runmad(out.spct[["Rfr"]], 7, endrule="mad")
    out.spct[["runmad"]] <- zoo::rollapply(out.spct[["Rfr"]], width = 7, FUN = stats::mad, partial = TRUE)
    out.spct[["runmed3"]] <- stats::runmed(out.spct[["Rfr"]], 3, endrule="median")
    out.spct[["runmed7"]] <- stats::runmed(out.spct[["Rfr"]], 7, endrule="median")
    out.spct[["runmed19"]] <- stats::runmed(out.spct[["Rfr"]], 19, endrule="median")
#    out.spct[["runmin5"]] <- caTools::runmin(out.spct[["Rfr"]], 5)
    out.spct[["runmin5"]] <- zoo::rollapply(out.spct[["Rfr"]], width = 5, FUN = min, partial = TRUE)
    # we need to avoid division by 0.0 and we use zero_limit / 10 close enough to zero
    out.spct[["runmadmed"]] <- with(out.spct,
                                    ifelse(runmad < zero_limit_cnst * 1e-1 | runmed7 < zero_limit * 1e-1,
                                           0.0, runmad/abs(runmed7)))
    out.spct[["Rfr.sm"]] <- with(out.spct,
               ifelse( (runmed19 < zero_limit) | (runmin5 < zero_limit * 5e-2), 0.0,
                       ifelse((Rfr > smooth_threshold) | (runmadmed < smooth_limit), Rfr,
                              ifelse(runmadmed < 2 * smooth_limit, runmed3,
                                     ifelse(runmadmed < 4 * smooth_limit, runmed7, runmed19)))))
    out.spct[["Rfr"]] <- with(out.spct,
                              ifelse(w.length < smoothing_hi_lim, Rfr.sm, Rfr))
    out.spct[["Rfr.good"]] <- out.spct[["runmadmed"]] / out.spct[["runmed19"]] * max(out.spct[["runmed19"]]) < 1.0
    if (anyNA(out.spct[["Rfr"]])) {
      warning(sum(is.na(out.spct[["Rfr"]])), " NAs in spectral irradiance")
    }
    num_bad <- sum(!out.spct[["Rfr.good"]], na.rm=TRUE)
    if (num_bad > length(out.spct) / 20) {
      message(num_bad, " possibly 'bad' values in smoothed spectral Rfr")
    }
    out.spct <- out.spct[ , c("w.length", "Rfr")]
    setReflectorSpct(out.spct)
    if (!is.null(comment(x))) {
      comment.text <- paste("Smoothed using 'custom', smooth_limit =",
                                 signif(smooth_limit, 3), "\n\n", comment(x))
    } else {
      comment.text <- paste("Smoothed using 'custom', smooth_limit =", signif(smooth_limit, 3))
    }
  }
  out.spct <- copy_attributes(x, out.spct)
  comment(out.spct) <- comment.text
  check_spct(out.spct, force = FALSE)
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
            num.spectra, " spectra")
    return(x)
  }
  e.and.q.input <- all(c("s.e.irrad", "s.q.irrad") %in% names(x))

  if ("s.e.response" %in% names(x) && anyNA(x[["s.e.response"]])) {
    if (na.rm) {
      message("Removing NA values at ", length(is.na(x[["s.e.response"]])), " wavelengths.")
      x <- na.omit(q2e(x, action = "replace"))
    } else {
      warning("NAs encountered when smoothing, returning input unchanged.")
      return(x)
    }
  } else if ("s.q.response" %in% names(x) && anyNA(x[["s.q.response"]])) {
    if (na.rm) {
      message("Removing NA values at ", length(is.na(x[["s.q.response"]])), " wavelengths.")
      x <- na.omit(e2q(x, action = "replace"))
    } else {
      warning("NAs encountered when smoothing, returning input unchanged.")
      return(x)
    }
  }

  # Skip checks for intermediate results
  # as intermediate values may be off-range
  prev_state <- disable_check_spct()
  on.exit(set_check_spct(prev_state), add = TRUE)
  if (method == "lowess") {
    span = 1/50 * strength
    if ("s.e.response" %in% names(x)) {
      out.spct <- stats::lowess(x[["w.length"]], x[["s.e.response"]], f = span, ...)
      names(out.spct) <- c("w.length", "s.e.response")
    } else if ("s.q.response" %in% names(x)) {
      out.spct <- stats::lowess(x[["w.length"]], x[["s.q.response"]], f = span, ...)
      names(out.spct) <- c("w.length", "s.q.response")
    }
    setResponseSpct(out.spct, time.unit = attr(x, "time.unit", exact = TRUE))
    if (e.and.q.input) {
      e2q(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      comment.text <- paste("Smoothed using 'lowess', f =", signif(span, 3), "\n\n", comment(x))
    } else {
      comment.text <- paste("Smoothed using 'lowess', f =", signif(span, 3))
    }
  } else if (method == "supsmu") {
    span = 1/50 * strength
    if ("s.e.response" %in% names(x)) {
      out.spct <- stats::supsmu(x[["w.length"]], x[["s.e.response"]], span = span, ...)
      names(out.spct) <- c("w.length", "s.e.response")
    } else if ("s.q.response" %in% names(x)) {
      out.spct <- stats::supsmu(x[["w.length"]], x[["s.q.response"]], span = span, ...)
      names(out.spct) <- c("w.length", "s.q.response")
    }
    setResponseSpct(out.spct, time.unit = attr(x, "time.unit", exact = TRUE))
    if (e.and.q.input) {
      e2q(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      comment.text <- paste("Smoothed using 'supsmu', span =", signif(span, 3), "\n\n", comment(x))
    } else {
      comment.text <- paste("Smoothed using 'supsmu', span =", signif(span, 3))
    }
  } else if (method == "custom") {
    # my own and inefficient method!
    # as the spectrum is already in energy units, we need to normalize thresholds
    out.spct <- x # we make a working copy
    q2e(out.spct, action = "replace", byref = TRUE)
    max_response <- max(out.spct[["s.e.response"]], na.rm=TRUE)
    smoothing_coef <- 1
    smoothing_hi_lim <- max(out.spct[["w.length"]])
    # this could be tweaked in many ways...
    zero_limit_cnst <- max_response * 3e-4 / strength
    out.spct[["zero_limit"]] <-  (zero_limit_cnst * 600) / out.spct[["w.length"]]
    smooth_limit <- 1e-3 * smoothing_coef # just a guess for runmadmed
    smooth_threshold <- 5e-2 * max_response / strength # for s.e.response
#    out.spct[["runmad"]] <- caTools::runmad(out.spct[["s.e.response"]], 7, endrule="mad")
    out.spct[["runmad"]] <- zoo::rollapply(out.spct[["s.e.response"]], width = 7, FUN = stats::mad, partial = TRUE)
    out.spct[["runmed3"]] <- stats::runmed(out.spct[["s.e.response"]], 3, endrule="median")
    out.spct[["runmed7"]] <- stats::runmed(out.spct[["s.e.response"]], 7, endrule="median")
    out.spct[["runmed19"]] <- stats::runmed(out.spct[["s.e.response"]], 19, endrule="median")
#    out.spct[["runmin5"]] <- caTools::runmin(out.spct[["s.e.response"]], 5)
    out.spct[["runmin5"]] <- zoo::rollapply(out.spct[["s.e.response"]], width = 5, FUN = min, partial = TRUE)
    # we need to avoid division by 0.0 and we use zero_limit / 10 close enough to zero
    out.spct[["runmadmed"]] <- with(out.spct,
                                    ifelse(runmad < zero_limit_cnst * 1e-1 | runmed7 < zero_limit * 1e-1,
                                           0.0, runmad/abs(runmed7)))
    out.spct[["s.e.response.sm"]] <- with(out.spct,
               ifelse( (runmed19 < zero_limit) | (runmin5 < zero_limit * 5e-2), 0.0,
                       ifelse((s.e.response > smooth_threshold) | (runmadmed < smooth_limit), s.e.response,
                              ifelse(runmadmed < 2 * smooth_limit, runmed3,
                                     ifelse(runmadmed < 4 * smooth_limit, runmed7, runmed19)))))
    out.spct[["s.e.response"]] <- with(out.spct,
                                       ifelse(w.length < smoothing_hi_lim, s.e.response.sm, s.e.response))
    out.spct[["s.e.response.good"]] <- out.spct[["runmadmed"]] / out.spct[["runmed19"]] * max(out.spct[["runmed19"]]) < 1.0
    if (anyNA(out.spct[["s.e.response"]])) {
      warning(sum(is.na(out.spct[["s.e.response"]])), " NAs in spectral response")
    }
    num_bad <- sum(!out.spct[["s.e.response.good"]], na.rm=TRUE)
    if (num_bad > length(out.spct) / 20) {
      message(num_bad, " possibly 'bad' values in smoothed spectral response")
    }
    out.spct <- out.spct[ , c("w.length", "s.e.response")]
    setResponseSpct(out.spct, time.unit = attr(x, "time.unit", exact = TRUE))
    if (e.and.q.input) {
      e2q(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      comment.text <- paste("Smoothed using 'custom', smooth_limit =",
                                 signif(smooth_limit, 3), "\n\n", comment(x))
    } else {
      comment.text <- paste("Smoothed using 'custom', smooth_limit =", signif(smooth_limit, 3))
    }
  }
  out.spct <- copy_attributes(x, out.spct)
  comment(out.spct) <- comment.text
  check_spct(out.spct, force = FALSE)
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
  s.e.response.good <- runmadmed / runmed19 * max(runmed19) < 1.0
  if (anyNA(z)) {
    warning(sum(is.na(z)), " NAs in smoothed spectrum")
  }
  num_bad <- sum(!s.e.response.good, na.rm=TRUE)
  if (num_bad > length(x) / 20) {
    message(num_bad, " possibly 'bad' values in smoothed spectral response")
  }
  # returned value similar to that of R's lowess() and supsmu()
  list(x = x, y = z)
}
