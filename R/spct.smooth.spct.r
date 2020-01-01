#' Smooth a spectrum
#'
#' These functions implement one original methods and acts as a wrapper for
#' other common R smoothing functions. The advantage of using this function for
#' smoothing spectral objects is that it simplifies the user interface and sets,
#' when needed, defaults suitable for spectral data.
#'
#' @param x an R object.
#' @param method a character string "custom", "lowess", "supsmu".
#' @param strength numeric value to adjust the degree of smoothing.
#' @param na.rm	logical A value indicating whether NA values should be stripped
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
smooth_spct <- function(x, method, strength, ...) UseMethod("smooth_spct")

#' @describeIn smooth_spct Default for generic function
#'
#' @export
#'
smooth_spct.default <- function(x, method, strength, ...) {
  warning("'smooth_spct' is not defined for objects of class ", class(x)[1])
  return(x)
}

#' @describeIn smooth_spct Smooth a source spectrum
#'
#' @export
#'
smooth_spct.source_spct <- function(x, method = "custom", strength = 1, na.rm = FALSE, ...) {
  num.spectra <- getMultipleWl(x)
  if (num.spectra != 1) {
    warning("Skipping smoothing as object contains ",
            num.spectra, " spectra in long form.")
    return(x)
  }
  e.and.q.input <- all(c("s.e.irrad", "s.q.irrad") %in% names(x))

  if ("s.e.irrad" %in% names(x) && anyNA(x[["s.e.irrad"]])) {
    if (na.rm) {
      message("Removing NA values at ", length(is.na(x[["s.e.irrad"]])), " wavelengths.")
      x <- na.omit(q2e(x, action = "replace"))
    } else {
      warning("NAs encountered when smoothing, returning input unchanged.")
      return(x)
    }
  } else if ("s.q.irrad" %in% names(x) && anyNA(x[["s.q.irrad"]])) {
    if (na.rm) {
      message("Removing NA values at ", length(is.na(x[["s.q.irrad"]])), " wavelengths.")
      x <- na.omit(e2q(x, action = "replace"))
    } else {
      warning("NAs encountered when smoothing, returning input unchanged.")
      return(x)
    }
  }

  # we disable range checks for spectra until end of function
  old.options <- options(photobiology.strict.range = NA_integer_)
  if (method == "lowess") {
    span = 1/50 * strength
    if ("s.e.irrad" %in% names(x)) {
      out.spct <- stats::lowess(x$w.length, x$s.e.irrad, f = span, ...)
      names(out.spct) <- c("w.length", "s.e.irrad")
    } else if ("s.q.irrad" %in% names(x)) {
      out.spct <- stats::lowess(x$w.length, x$s.q.irrad, f = span, ...)
      names(out.spct) <- c("w.length", "s.q.irrad")
    }
    setSourceSpct(out.spct, time.unit = getTimeUnit(x), bswf.used = getBSWFUsed(x))
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
    if ("s.e.irrad" %in% names(x)) {
      out.spct <- stats::supsmu(x$w.length, x$s.e.irrad, span = span, ...)
      names(out.spct) <- c("w.length", "s.e.irrad")
    } else if ("s.q.irrad" %in% names(x)) {
      out.spct <- stats::supsmu(x$w.length, x$s.q.irrad, span = span, ...)
      names(out.spct) <- c("w.length", "s.q.irrad")
    }
    setSourceSpct(out.spct, time.unit = attr(x, "time.unit", exact = TRUE))
    if (e.and.q.input) {
      e2q(out.spct, action = "add", byref = TRUE)
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
    q2e(out.spct, action = "replace", byref = TRUE)
    max_irrad <- max(out.spct[["s.e.irrad"]], na.rm=TRUE)
    smoothing_coef <- 1
    smoothing_hi_lim <- max(out.spct$w.length)
    # this could be tweaked in many ways...
    zero_limit_cnst <- max_irrad * 3e-4 / strength
    out.spct[["zero_limit"]] <-  (zero_limit_cnst * 600) / out.spct[["w.length"]]
    smooth_limit <- 1e-3 * smoothing_coef # just a guess for runmadmed
    smooth_threshold <- 5e-2 * max_irrad / strength # for s.e.irrad
#    out.spct[["runmad"]] <- caTools::runmad(out.spct[["s.e.irrad"]], 7, endrule="mad")
    out.spct[["runmad"]] <- zoo::rollapply(out.spct[["s.e.irrad"]], width = 7, FUN = stats::mad, partial = TRUE)
    out.spct[["runmed3"]] <- stats::runmed(out.spct[["s.e.irrad"]], 3, endrule="median")
    out.spct[["runmed7"]] <- stats::runmed(out.spct[["s.e.irrad"]], 7, endrule="median")
    out.spct[["runmed19"]] <- stats::runmed(out.spct[["s.e.irrad"]], 19, endrule="median")
#    out.spct[["runmin5"]] <- caTools::runmin(out.spct[["s.e.irrad"]], 5)
    out.spct[["runmin5"]] <- zoo::rollapply(out.spct[["s.e.irrad"]], width = 5, FUN = min, partial = TRUE)
    # we need to avoid division by 0.0 and we use zero_limit / 10 close enough to zero
    out.spct[["runmadmed"]] <- with(out.spct,
                                    ifelse(runmad < zero_limit_cnst * 1e-1 | runmed7 < zero_limit * 1e-1,
                                           0.0, runmad/abs(runmed7)))
    out.spct[["s.e.irrad.sm"]] <-
               with(out.spct, ifelse( (runmed19 < zero_limit) | (runmin5 < zero_limit * 5e-2), 0.0,
                       ifelse((s.e.irrad > smooth_threshold) | (runmadmed < smooth_limit), s.e.irrad,
                              ifelse(runmadmed < 2 * smooth_limit, runmed3,
                                     ifelse(runmadmed < 4 * smooth_limit, runmed7, runmed19)))))
    out.spct[["s.e.irrad"]] <- with(out.spct, ifelse(w.length < smoothing_hi_lim, s.e.irrad.sm, s.e.irrad))
    out.spct[["s.e.irrad.good"]] <- out.spct[["runmadmed"]] / out.spct[["runmed19"]] * max(out.spct[["runmed19"]]) < 1.0
    if (with(out.spct, any(is.na(s.e.irrad)))) {
      warning(sum(with(out.spct, any(is.na(s.e.irrad)))), " NAs in spectral irradiance")
    }
    num_bad <- with(out.spct, sum(!s.e.irrad.good, na.rm = TRUE))
    if (num_bad > length(out.spct) / 20) {
      message(num_bad, " possibly 'bad' values in smoothed spectral irradiance")
    }
    out.spct <- out.spct[ , c("w.length", "s.e.irrad")]
    setSourceSpct(out.spct, time.unit = attr(x, "time.unit", exact = TRUE))
    if (e.and.q.input) {
      e2q(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      comment.text <- paste("Smoothed using 'custom', smooth_limit =",
                                 signif(smooth_limit, 3), "\n\n", comment(x))
    } else {
      comment.text <-  paste("Smoothed using 'custom', smooth_limit =",
                                  signif(smooth_limit, 3))
    }
  }
  options(old.options)
  out.spct <- copy_attributes(x, out.spct)
  comment(out.spct) <- comment.text
  check_spct(out.spct)
}

#' @describeIn smooth_spct Smooth a filter spectrum
#'
#' @export
#'
smooth_spct.filter_spct <- function(x, method = "custom", strength = 1, na.rm = FALSE, ...) {
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

  # we disable range checks for spectra until end of function
  old.options <- options(photobiology.strict.range = NA_integer_)
  if (method == "lowess") {
    span = 1/50 * strength
    if ("Tfr" %in% names(x)) {
      out.spct <- stats::lowess(x$w.length, x$Tfr, f = span, ...)
      names(out.spct) <- c("w.length", "Tfr")
    } else if ("A" %in% names(x)) {
      out.spct <- stats::lowess(x$w.length, x$A, f = span, ...)
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
      out.spct <- stats::supsmu(x$w.length, x$Tfr, span = span, ...)
      names(out.spct) <- c("w.length", "Tfr")
    } else if ("A" %in% names(x)) {
      out.spct <- stats::supsmu(x$w.length, x$A, span = span, ...)
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
    smoothing_hi_lim <- max(out.spct$w.length)
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
  options(old.options)
  out.spct <- copy_attributes(x, out.spct)
  comment(out.spct) <- comment.text
  check_spct(out.spct)
}


#' @describeIn smooth_spct Smooth a reflector spectrum
#'
#' @export
#'
smooth_spct.reflector_spct <- function(x, method = "custom", strength = 1, na.rm = FALSE, ...) {
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
  # we disable range checks for spectra until end of function
  old.options <- options(photobiology.strict.range = NA_integer_)
  if (method == "lowess") {
    span = 1/50 * strength
    if ("Rfr" %in% names(x)) {
      out.spct <- stats::lowess(x$w.length, x$Rfr, f = span, ...)
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
      out.spct <- stats::supsmu(x$w.length, x$Rfr, span = span, ...)
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
    smoothing_hi_lim <- max(out.spct$w.length)
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
  options(old.options)
  out.spct <- copy_attributes(x, out.spct)
  comment(out.spct) <- comment.text
  check_spct(out.spct)
}

#' @describeIn smooth_spct Smooth a response spectrum
#'
#' @export
#'
smooth_spct.response_spct <- function(x, method = "custom", strength = 1, na.rm = FALSE, ...) {
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

  # we disable range checks for spectra until end of function
  old.options <- options(photobiology.strict.range = NA_integer_)
  if (method == "lowess") {
    span = 1/50 * strength
    if ("s.e.response" %in% names(x)) {
      out.spct <- stats::lowess(x$w.length, x$s.e.response, f = span, ...)
      names(out.spct) <- c("w.length", "s.e.response")
    } else if ("s.q.response" %in% names(x)) {
      out.spct <- stats::lowess(x$w.length, x$s.q.response, f = span, ...)
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
      out.spct <- stats::supsmu(x$w.length, x$s.e.response, span = span, ...)
      names(out.spct) <- c("w.length", "s.e.response")
    } else if ("s.q.response" %in% names(x)) {
      out.spct <- stats::supsmu(x$w.length, x$s.q.response, span = span, ...)
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
    smoothing_hi_lim <- max(out.spct$w.length)
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
  options(old.options)
  out.spct <- copy_attributes(x, out.spct)
  comment(out.spct) <- comment.text
  check_spct(out.spct)
}

#' @describeIn smooth_spct
#'
#' @export
#'
smooth_spct.generic_mspct <-
  function(x, method = "custom", strength = 1, na.rm = FALSE, ...) {
    msmsply(x, smooth_spct, method = method, strength = strength,
            na.rm = na.rm, ...)
  }
