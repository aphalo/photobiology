#' Generic function
#'
#' Smooth a spectrum.
#'
#' @param x an R object
#' @param method a character string "custom", "lowess", "supsmu"
#' @param strength numeric value to adjust the degree of smoothing
#' @param ... other parameters passed to the underlying smoothing functions
#'
#' @export smooth_spct
#'
smooth_spct <- function(x, method, strength, ...) UseMethod("smooth_spct")

#' Default for generic function
#'
#' Smooth a spectrum.
#'
#' @param x an R object
#' @param method a character string "custom", "lowess", "supsmu"
#' @param strength numeric value to adjust the degree of smoothing
#' @param ... other parameters passed to the underlying smoothing functions
#'
#' @export smooth_spct.default
#'
smooth_spct.default <- function(x, method, strength, ...) {
  return(x)
}

#' Smooth a source spectrum using one of several different methods
#'
#' This function implements one original methods and acts as a wrapper for other
#' common R smoothing functions. The advantage of using this function for smoothing
#' \code{source.spct} objects is that it simplifies the user interface and sets,
#' when needed, defaults suitable for spectral irardiance data.
#'
#' @usage smooth_spct.source.spct(x, method = "custom", strength = 1, ...)
#'
#' @param x a source.spct object
#' @param method a character string "custom", "lowess", "supsmu"
#' @param strength numeric value to adjust the degree of smoothing
#' @param ... other parameters passed to the underlying smoothing functions
#'
#' @return a source.spct object smoothed
#'
#' @export smooth_spct.source.spct
#' @importFrom caTools runmad runmin
#'
#' @note Method "custom" is our home-brewed method which applies strong smoothing
#' to low signal regions of the spectral data, and weaker or no smoothing to the
#' high signal areas. Values very close to zero are set to zero with a limit which
#' depends on the local variation. This method is an ad-hock method suitable for
#' smoothing spectral data obtained with array spectrometers. In the cased of methods
#' "lowess" and "supsmu" the current function behaves like a wrapper of the functions
#' of the same names from base R.
#'
#' @keywords manip misc

smooth_spct.source.spct <- function(x, method = "custom", strength = 1, ...) {
  if (method == "lowess") {
    span = 1/50 * strength
    if ("s.e.irrad" %in% names(x)) {
      out.spct <- lowess(x$w.length, x$s.e.irrad, f = span, ...)
      names(out.spct) <- c("w.length", "s.e.irrad")
    } else if ("s.q.irrad" %in% names(x)) {
      out.spct <- lowess(x$w.length, x$s.q.irrad, f = span, ...)
      names(out.spct) <- c("w.length", "s.q.irrad")
    }
    setSourceSpct(out.spct, time.unit = getTimeUnit(x))
    if (all(c("s.e.irrad", "s.q.irrad") %in% names(x))) {
      e2q(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      setattr(out.spct, "comment", paste("Smoothed using 'lowess', f =", signif(span, 3), "\n\n", comment(x)))
    } else {
      setattr(out.spct, "comment", paste("Smoothed using 'lowess', f =", signif(span, 3)))
    }
    return(out.spct)
  } else if (method == "supsmu") {
    span = 1/50 * strength
    if ("s.e.irrad" %in% names(x)) {
      out.spct <- supsmu(x$w.length, x$s.e.irrad, span = span, ...)
      names(out.spct) <- c("w.length", "s.e.irrad")
    } else if ("s.q.irrad" %in% names(x)) {
      out.spct <- supsmu(x$w.length, x$s.q.irrad, span = span, ...)
      names(out.spct) <- c("w.length", "s.q.irrad")
    }
    setSourceSpct(out.spct, time.unit = attr(x, "time.unit", exact = TRUE))
    if (all(c("s.e.irrad", "s.q.irrad") %in% names(x))) {
      e2q(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      setattr(out.spct, "comment", paste("Smoothed using 'supsmu', span =", signif(span, 3), "\n\n", comment(x)))
    } else {
      setattr(out.spct, "comment", paste("Smoothed using 'supsmu', span =", signif(span, 3)))
    }
    return(out.spct)
  } else if (method == "custom") {
    # my own and inefficient method!
    # as the spectrum is already in energy units, we need to normalize thresholds
    out.spct <- copy(x) # just to avoid editing the code
    q2e(out.spct, action = "replace", byref = TRUE)
    max_irrad <- out.spct[ , max(s.e.irrad, na.rm=TRUE)]
    smoothing_coef <- 1
    smoothing_hi_lim <- max(out.spct$w.length)
    # this could be tweeked in many ways...
    zero_limit_cnst <- max_irrad * 3e-4 / strength
    out.spct[ , zero_limit :=  (zero_limit_cnst * 600) / w.length]
    smooth_limit <- 1e-3 * smoothing_coef # just a guess for runmadmed
    smooth_threshold <- 5e-2 * max_irrad / strength # for s.e.irrad
    out.spct[ , runmad := runmad(s.e.irrad, 7, endrule="mad")]
    out.spct[ , runmed3 := runmed(s.e.irrad, 3, endrule="median")]
    out.spct[ , runmed7 := runmed(s.e.irrad, 7, endrule="median")]
    out.spct[ , runmed19 := runmed(s.e.irrad, 19, endrule="median")]
    out.spct[ , runmin5 := runmin(s.e.irrad, 5)]
    # we need to avoid division by 0.0 and we use zero_limit / 10 close enough to zero
    out.spct[ , runmadmed := ifelse(runmad < zero_limit_cnst * 1e-1 | runmed7 < zero_limit * 1e-1, 0.0, runmad/abs(runmed7))]
    out.spct[ , s.e.irrad.sm :=
               ifelse( (runmed19 < zero_limit) | (runmin5 < zero_limit * 5e-2), 0.0,
                       ifelse((s.e.irrad > smooth_threshold) | (runmadmed < smooth_limit), s.e.irrad,
                              ifelse(runmadmed < 2 * smooth_limit, runmed3,
                                     ifelse(runmadmed < 4 * smooth_limit, runmed7, runmed19))))]
    out.spct[w.length < smoothing_hi_lim, s.e.irrad := s.e.irrad.sm]
    out.spct[ , s.e.irrad.good := runmadmed < 1.0]
    if (out.spct[ , any(is.na(s.e.irrad))]) {
      warning(out.spct[ , sum(is.na(s.e.irrad))], " NAs in spectral irradiance")
    }
    num_bad <- out.spct[ , sum(!s.e.irrad.good, na.rm=TRUE)]
    if (num_bad > 50) {
      warning(num_bad, " 'bad' estimates in spectral irradiance")
    }
    out.spct <- out.spct[ , .(w.length, s.e.irrad)]
    setSourceSpct(out.spct, time.unit = attr(x, "time.unit", exact = TRUE))
    if (all(c("s.e.irrad", "s.q.irrad") %in% names(x))) {
      e2q(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      setattr(out.spct, "comment", paste("Smoothed using 'custom', smooth_limit =", signif(smooth_limit, 3), "\n\n", comment(x)))
    } else {
      setattr(out.spct, "comment", paste("Smoothed using 'custom', smooth_limit =", signif(smooth_limit, 3)))
    }
    return(out.spct)
  }
}

#' Smooth a filter spectrum using one of several different methods
#'
#' This function implements one original methods and acts as a wrapper for other
#' common R smoothing functions. The advantage of using this function for smoothing
#' \code{filter.spct} objects is that it simplifies the user interface and sets,
#' when needed, defaults suitable for spectral irardiance data.
#'
#' @usage smooth_spct.filter.spct(x, method = "custom", strength = 1, ...)
#'
#' @param x a filter.spct object
#' @param method a character string "custom", "lowess", "supsmu"
#' @param strength numeric value to adjust the degree of smoothing
#' @param ... other parameters passed to the underlying smoothing functions
#'
#' @return a filter.spct object smoothed
#'
#' @export smooth_spct.filter.spct
#' @importFrom caTools runmad runmin
#'
#' @note Method "custom" is our home-brewed method which applies strong smoothing
#' to low signal regions of the spectral data, and weaker or no smoothing to the
#' high signal areas. Values very close to zero are set to zero with a limit which
#' depends on the local variation. This method is an ad-hock method suitable for
#' smoothing spectral data obtained with array spectrometers. In the cased of methods
#' "lowess" and "supsmu" the current function behaves like a wrapper of the functions
#' of the same names from base R.
#'
#' @keywords manip misc

smooth_spct.filter.spct <- function(x, method = "custom", strength = 1, ...) {
  if (method == "lowess") {
    span = 1/50 * strength
    if ("Tfr" %in% names(x)) {
      out.spct <- lowess(x$w.length, x$Tfr, f = span, ...)
      names(out.spct) <- c("w.length", "Tfr")
    } else if ("A" %in% names(x)) {
      out.spct <- lowess(x$w.length, x$A, f = span, ...)
      names(out.spct) <- c("w.length", "A")
    }
    setFilterSpct(out.spct, Tfr.type = attr(x, "Tfr.type", exact = TRUE))
    if (all(c("Tfr", "A") %in% names(x))) {
      T2A(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      setattr(out.spct, "comment", paste("Smoothed using 'lowess', f =", signif(span, 3), "\n\n", comment(x)))
    } else {
      setattr(out.spct, "comment", paste("Smoothed using 'lowess', f =", signif(span, 3)))
    }
    return(out.spct)
  } else if (method == "supsmu") {
    span = 1/50 * strength
    if ("Tfr" %in% names(x)) {
      out.spct <- supsmu(x$w.length, x$Tfr, span = span, ...)
      names(out.spct) <- c("w.length", "Tfr")
    } else if ("A" %in% names(x)) {
      out.spct <- supsmu(x$w.length, x$A, span = span, ...)
      names(out.spct) <- c("w.length", "A")
    }
    setFilterSpct(out.spct, Tfr.type = attr(x, "Tfr.type", exact = TRUE))
    if (all(c("Tfr", "A") %in% names(x))) {
      T2A(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      setattr(out.spct, "comment", paste("Smoothed using 'supsmu', span =", signif(span, 3), "\n\n", comment(x)))
    } else {
      setattr(out.spct, "comment", paste("Smoothed using 'supsmu', span =", signif(span, 3)))
    }
    return(out.spct)
  } else if (method == "custom") {
    # my own and inefficient method!
    # as the spectrum is already in energy units, we need to normalize thresholds
    out.spct <- copy(x) # just to avoid editing the code
    A2T(out.spct, action = "replace", byref = TRUE)
    max_Tfr <- 1
    smoothing_coef <- 1
    smoothing_hi_lim <- max(out.spct$w.length)
    # this could be tweeked in many ways...
    zero_limit_cnst <- max_Tfr * 3e-4 / strength
    out.spct[ , zero_limit :=  (zero_limit_cnst * 600) / w.length]
    smooth_limit <- 1e-3 * smoothing_coef # just a guess for runmadmed
    smooth_threshold <- 5e-2 * max_Tfr / strength # for Tfr
    out.spct[ , runmad := runmad(Tfr, 7, endrule="mad")]
    out.spct[ , runmed3 := runmed(Tfr, 3, endrule="median")]
    out.spct[ , runmed7 := runmed(Tfr, 7, endrule="median")]
    out.spct[ , runmed19 := runmed(Tfr, 19, endrule="median")]
    out.spct[ , runmin5 := runmin(Tfr, 5)]
    # we need to avoid division by 0.0 and we use zero_limit / 10 close enough to zero
    out.spct[ , runmadmed := ifelse(runmad < zero_limit_cnst * 1e-1 | runmed7 < zero_limit * 1e-1, 0.0, runmad/abs(runmed7))]
    out.spct[ , Tfr.sm :=
               ifelse( (runmed19 < zero_limit) | (runmin5 < zero_limit * 5e-2), 0.0,
                       ifelse((Tfr > smooth_threshold) | (runmadmed < smooth_limit), Tfr,
                              ifelse(runmadmed < 2 * smooth_limit, runmed3,
                                     ifelse(runmadmed < 4 * smooth_limit, runmed7, runmed19))))]
    out.spct[w.length < smoothing_hi_lim, Tfr := Tfr.sm]
    out.spct[ , Tfr.good := runmadmed < 1.0]
    if (out.spct[ , any(is.na(Tfr))]) {
      warning(out.spct[ , sum(is.na(Tfr))], " NAs in spectral irradiance")
    }
    num_bad <- out.spct[ , sum(!Tfr.good, na.rm=TRUE)]
    if (num_bad > 50) {
      warning(num_bad, " 'bad' estimates in spectral irradiance")
    }
    out.spct <- out.spct[ , .(w.length, Tfr)]
    setFilterSpct(out.spct, Tfr.type = attr(x, "Tfr.type", exact = TRUE))
    if (all(c("Tfr", "A") %in% names(x))) {
      T2A(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      setattr(out.spct, "comment", paste("Smoothed using 'custom', smooth_limit =", signif(smooth_limit, 3), "\n\n", comment(x)))
    } else {
      setattr(out.spct, "comment", paste("Smoothed using 'custom', smooth_limit =", signif(smooth_limit, 3)))
    }
    return(out.spct)
  }
}


#' Smooth a reflector spectrum using one of several different methods
#'
#' This function implements one original methods and acts as a wrapper for other
#' common R smoothing functions. The advantage of using this function for smoothing
#' \code{reflector.spct} objects is that it simplifies the user interface and sets,
#' when needed, defaults suitable for spectral irardiance data.
#'
#' @usage smooth_spct.reflector.spct(x, method = "custom", strength = 1, ...)
#'
#' @param x a reflector.spct object
#' @param method a character string "custom", "lowess", "supsmu"
#' @param strength numeric value to adjust the degree of smoothing
#' @param ... other parameters passed to the underlying smoothing functions
#'
#' @return a reflector.spct object smoothed
#'
#' @export smooth_spct.reflector.spct
#' @importFrom caTools runmad runmin
#'
#' @note Method "custom" is our home-brewed method which applies strong smoothing
#' to low signal regions of the spectral data, and weaker or no smoothing to the
#' high signal areas. Values very close to zero are set to zero with a limit which
#' depends on the local variation. This method is an ad-hock method suitable for
#' smoothing spectral data obtained with array spectrometers. In the cased of methods
#' "lowess" and "supsmu" the current function behaves like a wrapper of the functions
#' of the same names from base R.
#'
#' @keywords manip misc

smooth_spct.reflector.spct <- function(x, method = "custom", strength = 1, ...) {
  if (method == "lowess") {
    span = 1/50 * strength
    if ("Rfr" %in% names(x)) {
      out.spct <- lowess(x$w.length, x$Rfr, f = span, ...)
      names(out.spct) <- c("w.length", "Rfr")
    }
    setReflectorSpct(out.spct)
    if (!is.null(comment(x))) {
      setattr(out.spct, "comment", paste("Smoothed using 'lowess', f =", signif(span, 3), "\n\n", comment(x)))
    } else {
      setattr(out.spct, "comment", paste("Smoothed using 'lowess', f =", signif(span, 3)))
    }
    return(out.spct)
  } else if (method == "supsmu") {
    span = 1/50 * strength
    if ("Rfr" %in% names(x)) {
      out.spct <- supsmu(x$w.length, x$Rfr, span = span, ...)
      names(out.spct) <- c("w.length", "Rfr")
    }
    setReflectorSpct(out.spct)
    if (!is.null(comment(x))) {
      setattr(out.spct, "comment", paste("Smoothed using 'supsmu', span =", signif(span, 3), "\n\n", comment(x)))
    } else {
      setattr(out.spct, "comment", paste("Smoothed using 'supsmu', span =", signif(span, 3)))
    }
    return(out.spct)
  } else if (method == "custom") {
    # my own and inefficient method!
    # as the spectrum is already in energy units, we need to normalize thresholds
    out.spct <- copy(x) # just to avoid editing the code
    max_Rfr <- 1
    smoothing_coef <- 1
    smoothing_hi_lim <- max(out.spct$w.length)
    # this could be tweeked in many ways...
    zero_limit_cnst <- max_irrad * 3e-4 / strength
    out.spct[ , zero_limit :=  (zero_limit_cnst * 600) / w.length]
    smooth_limit <- 1e-3 * smoothing_coef # just a guess for runmadmed
    smooth_threshold <- 5e-2 * max_Rfr / strength # for Rfr
    out.spct[ , runmad := runmad(Rfr, 7, endrule="mad")]
    out.spct[ , runmed3 := runmed(Rfr, 3, endrule="median")]
    out.spct[ , runmed7 := runmed(Rfr, 7, endrule="median")]
    out.spct[ , runmed19 := runmed(Rfr, 19, endrule="median")]
    out.spct[ , runmin5 := runmin(Rfr, 5)]
    # we need to avoid division by 0.0 and we use zero_limit / 10 close enough to zero
    out.spct[ , runmadmed := ifelse(runmad < zero_limit_cnst * 1e-1 | runmed7 < zero_limit * 1e-1, 0.0, runmad/abs(runmed7))]
    out.spct[ , Rfr.sm :=
               ifelse( (runmed19 < zero_limit) | (runmin5 < zero_limit * 5e-2), 0.0,
                       ifelse((Rfr > smooth_threshold) | (runmadmed < smooth_limit), Rfr,
                              ifelse(runmadmed < 2 * smooth_limit, runmed3,
                                     ifelse(runmadmed < 4 * smooth_limit, runmed7, runmed19))))]
    out.spct[w.length < smoothing_hi_lim, Rfr := Rfr.sm]
    out.spct[ , Rfr.good := runmadmed < 1.0]
    if (out.spct[ , any(is.na(Rfr))]) {
      warning(out.spct[ , sum(is.na(Rfr))], " NAs in spectral irradiance")
    }
    num_bad <- out.spct[ , sum(!Rfr.good, na.rm=TRUE)]
    if (num_bad > 50) {
      warning(num_bad, " 'bad' estimates in spectral irradiance")
    }
    out.spct <- out.spct[ , .(w.length, Rfr)]
    setReflectorSpct(out.spct)
    if (all(c("Rfr", "A") %in% names(x))) {
      T2A(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      setattr(out.spct, "comment", paste("Smoothed using 'custom', smooth_limit =", signif(smooth_limit, 3), "\n\n", comment(x)))
    } else {
      setattr(out.spct, "comment", paste("Smoothed using 'custom', smooth_limit =", signif(smooth_limit, 3)))
    }
    return(out.spct)
  }
}

#' Smooth a response spectrum using one of several different methods
#'
#' This function implements one original methods and acts as a wrapper for other
#' common R smoothing functions. The advantage of using this function for smoothing
#' \code{response.spct} objects is that it simplifies the user interface and sets,
#' when needed, defaults suitable for spectral irardiance data.
#'
#' @usage smooth_spct.response.spct(x, method = "custom", strength = 1, ...)
#'
#' @param x a response.spct object
#' @param method a character string "custom", "lowess", "supsmu"
#' @param strength numeric value to adjust the degree of smoothing
#' @param ... other parameters passed to the underlying smoothing functions
#'
#' @return a response.spct object smoothed
#'
#' @export smooth_spct.response.spct
#' @importFrom caTools runmad runmin
#'
#' @note Method "custom" is our home-brewed method which applies strong smoothing
#' to low signal regions of the spectral data, and weaker or no smoothing to the
#' high signal areas. Values very close to zero are set to zero with a limit which
#' depends on the local variation. This method is an ad-hock method suitable for
#' smoothing spectral data obtained with array spectrometers. In the cased of methods
#' "lowess" and "supsmu" the current function behaves like a wrapper of the functions
#' of the same names from base R.
#'
#' @keywords manip misc

smooth_spct.response.spct <- function(x, method = "custom", strength = 1, ...) {
  if (method == "lowess") {
    span = 1/50 * strength
    if ("s.e.response" %in% names(x)) {
      out.spct <- lowess(x$w.length, x$s.e.response, f = span, ...)
      names(out.spct) <- c("w.length", "s.e.response")
    } else if ("s.q.response" %in% names(x)) {
      out.spct <- lowess(x$w.length, x$s.q.response, f = span, ...)
      names(out.spct) <- c("w.length", "s.q.response")
    }
    setResponseSpct(out.spct, time.unit = attr(x, "time.unit", exact = TRUE))
    if (all(c("s.e.response", "s.q.response") %in% names(x))) {
      e2q(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      setattr(out.spct, "comment", paste("Smoothed using 'lowess', f =", signif(span, 3), "\n\n", comment(x)))
    } else {
      setattr(out.spct, "comment", paste("Smoothed using 'lowess', f =", signif(span, 3)))
    }
    return(out.spct)
  } else if (method == "supsmu") {
    span = 1/50 * strength
    if ("s.e.response" %in% names(x)) {
      out.spct <- supsmu(x$w.length, x$s.e.response, span = span, ...)
      names(out.spct) <- c("w.length", "s.e.response")
    } else if ("s.q.response" %in% names(x)) {
      out.spct <- supsmu(x$w.length, x$s.q.response, span = span, ...)
      names(out.spct) <- c("w.length", "s.q.response")
    }
    setResponseSpct(out.spct, time.unit = attr(x, "time.unit", exact = TRUE))
    if (all(c("s.e.response", "s.q.response") %in% names(x))) {
      e2q(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      setattr(out.spct, "comment", paste("Smoothed using 'supsmu', span =", signif(span, 3), "\n\n", comment(x)))
    } else {
      setattr(out.spct, "comment", paste("Smoothed using 'supsmu', span =", signif(span, 3)))
    }
    return(out.spct)
  } else if (method == "custom") {
    # my own and inefficient method!
    # as the spectrum is already in energy units, we need to normalize thresholds
    out.spct <- copy(x) # just to avoid editing the code
    q2e(out.spct, action = "replace", byref = TRUE)
    max_response <- out.spct[ , max(s.e.response, na.rm=TRUE)]
    smoothing_coef <- 1
    smoothing_hi_lim <- max(out.spct$w.length)
    # this could be tweeked in many ways...
    zero_limit_cnst <- max_response * 3e-4 / strength
    out.spct[ , zero_limit :=  (zero_limit_cnst * 600) / w.length]
    smooth_limit <- 1e-3 * smoothing_coef # just a guess for runmadmed
    smooth_threshold <- 5e-2 * max_response / strength # for s.e.response
    out.spct[ , runmad := runmad(s.e.response, 7, endrule="mad")]
    out.spct[ , runmed3 := runmed(s.e.response, 3, endrule="median")]
    out.spct[ , runmed7 := runmed(s.e.response, 7, endrule="median")]
    out.spct[ , runmed19 := runmed(s.e.response, 19, endrule="median")]
    out.spct[ , runmin5 := runmin(s.e.response, 5)]
    # we need to avoid division by 0.0 and we use zero_limit / 10 close enough to zero
    out.spct[ , runmadmed := ifelse(runmad < zero_limit_cnst * 1e-1 | runmed7 < zero_limit * 1e-1, 0.0, runmad/abs(runmed7))]
    out.spct[ , s.e.response.sm :=
               ifelse( (runmed19 < zero_limit) | (runmin5 < zero_limit * 5e-2), 0.0,
                       ifelse((s.e.response > smooth_threshold) | (runmadmed < smooth_limit), s.e.response,
                              ifelse(runmadmed < 2 * smooth_limit, runmed3,
                                     ifelse(runmadmed < 4 * smooth_limit, runmed7, runmed19))))]
    out.spct[w.length < smoothing_hi_lim, s.e.response := s.e.response.sm]
    out.spct[ , s.e.response.good := runmadmed < 1.0]
    if (out.spct[ , any(is.na(s.e.response))]) {
      warning(out.spct[ , sum(is.na(s.e.response))], " NAs in spectral responseiance")
    }
    num_bad <- out.spct[ , sum(!s.e.response.good, na.rm=TRUE)]
    if (num_bad > 50) {
      warning(num_bad, " 'bad' estimates in spectral responseiance")
    }
    out.spct <- out.spct[ , .(w.length, s.e.response)]
    setResponseSpct(out.spct, time.unit = attr(x, "time.unit", exact = TRUE))
    if (all(c("s.e.response", "s.q.response") %in% names(x))) {
      e2q(out.spct, action = "add", byref = TRUE)
    }
    if (!is.null(comment(x))) {
      setattr(out.spct, "comment", paste("Smoothed using 'custom', smooth_limit =", signif(smooth_limit, 3), "\n\n", comment(x)))
    } else {
      setattr(out.spct, "comment", paste("Smoothed using 'custom', smooth_limit =", signif(smooth_limit, 3)))
    }
    return(out.spct)
  }
}