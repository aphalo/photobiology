# response methods --------------------------------------------------------


#' Integrated response
#'
#' Calculate average photon- or energy-based photo-response.
#'
#' @param spct an R object of class "generic_spct".
#' @param w.band waveband or list of waveband objects or a numeric vector of
#'   length two. The waveband(s) determine the region(s) of the spectrum that
#'   are summarized. If a numeric range is supplied a waveband object is
#'   constructed on the fly from it.
#' @param unit.out character Allowed values "energy", and "photon", or its alias
#'   "quantum".
#' @param quantity character string One of "total", "average" or "mean",
#'   "contribution", "contribution.pc", "relative" or "relative.pc".
#' @param time.unit character or lubridate::duration object.
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical Flag telling if wavebands crossing spectral data boundaries
#'   are trimmed or ignored.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @note The parameter \code{use.hinges} controls speed optimization. The
#'   defaults should be suitable in most cases. Only the range of wavelengths
#'   in the wavebands is used and all BSWFs are ignored.
#'
#' @return A named \code{numeric} vector in the case of methods for individual
#'   spectra, with one value for each \code{waveband} passed to parameter
#'   \code{w.band}. A \code{data.frame} in the case of collections of spectra,
#'   containing one column for each \code{waveband} object, an index column with
#'   the names of the spectra, and optionally additional columns with metadata
#'   values retrieved from the attributes of the member spectra.
#'
#'   Whether returned values are expressed in energy-based or photon-based units
#'   depends on \code{unit.out}. By default values are only integrated, but
#'   depending on the argument passed to parameter \code{quantity} they can be
#'   re-expressed as relative fractions or percentages. In the case of vector
#'   output, \code{names} attribute is set to the name of the corresponding
#'   waveband unless a named list is supplied in which case the names of the
#'   list members are used.
#'
#' @export
#' @family response functions
#'
response <- function(spct, w.band, unit.out, quantity, time.unit, scale.factor, wb.trim, use.hinges, ...) UseMethod("response")

#' @describeIn response Default for generic function
#'
#' @export
#'
response.default <- function(spct, w.band, unit.out, quantity, time.unit, scale.factor, wb.trim, use.hinges, ...) {
  warning("'response' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn response Method for response spectra.
#'
#' @export
#'
response.response_spct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           quantity = "total",
           time.unit = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL), ... ) {
    resp_spct(spct = spct,
              w.band = w.band,
              unit.out = unit.out,
              quantity = quantity,
              time.unit = time.unit,
              scale.factor = scale.factor,
              wb.trim = wb.trim,
              use.hinges = use.hinges )
  }

#' Calculate response from spectral response
#'
#' This function returns the mean response for a given waveband and a response
#' spectrum.
#'
#' @param spct an object of class response_spct".
#' @param w.band waveband or list of waveband objects or a numeric vector of
#'   length two. The waveband(s) determine the region(s) of the spectrum that
#'   are summarized. If a numeric range is supplied a waveband object is
#'   constructed on the fly from it.
#' @param unit.out character with allowed values "energy", and "photon", or its
#'   alias "quantum".
#' @param quantity character string One of "total", "average" or "mean",
#'   "contribution", "contribution.pc", "relative" or "relative.pc".
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @return a single numeric value expressed either as a fraction of one or a
#'   percentage, or a vector of the same length as the list of \code{waveband}
#'   objects.
#' @keywords internal
#'
#' @note The parameter \code{use.hinges} controls speed optimization. The
#'   defaults should be suitable in most cases. Only the range of wavelengths
#'   in the wavebands is used and all BSWFs are ignored.
#'
#' @keywords internal
#'
resp_spct <-
  function(spct,
           w.band,
           unit.out,
           quantity,
           time.unit,
           scale.factor,
           wb.trim,
           use.hinges,
           ...) {
    num.spectra <- getMultipleWl(spct)
    if (num.spectra != 1) {
      warning("Skipping response calculation as object contains ",
              num.spectra, " spectra")
      return(NA_real_)
    }
    if (is_normalized(spct)) {
      warning("The spectral data has been normalized, making impossible to calculate absorbance")
      return(NA_real_)
    }
    if (is_scaled(spct)) {
      warning("Summary calculated from rescaled data")
    }
    # makes "quantum" synonym for "photon" without changes to other code
    if (unit.out == "quantum") {
      unit.out <- "photon"
    }

    data.time.unit <- getTimeUnit(spct, force.duration = lubridate::is.duration(time.unit))

    if (!is.null(time.unit) && time.unit != data.time.unit) {
      if (!lubridate::is.duration(time.unit) && !is.character(time.unit)) {
        message("converting 'time.unit' ", time.unit, " into a lubridate::duration")
        time.unit <- lubridate::as.duration(time.unit)
      }
      spct <- convertTimeUnit(spct, time.unit = time.unit, byref = FALSE)
    } else {
      time.unit <- data.time.unit
    }

    if (unit.out == "photon") {
      spct <- e2q(spct)
      spct <- spct[ , c("w.length", "s.q.response")]
    } else if (unit.out == "energy") {
      spct <- q2e(spct)
      spct <- spct[ , c("w.length", "s.e.response")]
    } else {
      stop("Invalid 'unit.out'")
    }

    # if the waveband is undefined then use all data
    if (length(w.band) == 0) {
      w.band <- waveband(spct)
    }
    if (is.numeric(w.band)) {
      w.band <- waveband(w.band)
    }
    if (is.waveband(w.band)) {
      # if the argument is a single w.band, we enclose it in a list
      # so that the for loop works as expected. This is a bit of a
      # kludge but it let's us avoid treating it as a special case
      w.band <- list(w.band)
    }
    w.band <- trim_waveband(w.band = w.band, range = spct, trim = wb.trim)

    # if the w.band includes 'hinges' we insert them,
    # but if not, we decide whether to insert hinges or not
    # hinges or not based of the wavelength resolution of the
    # spectrum. This can produce small errors for high
    # spectral resolution data, but speed up the calculations.
    if (is.null(use.hinges)) {
      use.hinges <- auto_hinges(spct[["w.length"]])
    }

    # we collect all hinges and insert them in one go
    # this may alter very slightly the returned values
    # but improves calculation speed
    if (use.hinges) {
      all.hinges <- NULL
      for (wb in w.band) {
        if (!is.null(wb$hinges) && length(wb$hinges) > 0) {
          all.hinges <- c(all.hinges, wb$hinges)
        }
      }
      if (!is.null(all.hinges)) {
        spct <- insert_spct_hinges(spct, all.hinges)
      }
    }

    # we prepare labels for output
    wb_name <- names(w.band)
    no_names_flag <- is.null(wb_name)
    if (no_names_flag) {
      wb_name <- character(length(w.band))
    }

    # we iterate through the list of wavebands
    response <- double(length(w.band))
    i <- 0
    for (wb in w.band) {
      i <- i + 1
      # we get names from wb if needed
      if (no_names_flag) {
        if (is_effective(wb)) {
          warning("Using only wavelength range from a weighted waveband object.")
          wb_name[i] <- paste("range", as.character(signif(min(wb), 4)),
                              as.character(signif(max(wb), 4)), sep = ".")
        } else {
          wb_name[i] <- wb$name
        }
      }
      # we calculate the integrated response.
      response[i] <- integrate_spct(trim_spct(spct, wb, use.hinges = FALSE))
    }
    if (quantity %in% c("contribution", "contribution.pc")) {
      if (any(sapply(w.band, is_effective))) {
        warning("'quantity '", quantity,
                "' not supported when using BSWFs, returning 'total' instead")
        quantity <- "total"
      } else {
        total <- resp_spct(spct,
                           w.band = NULL,
                           unit.out = unit.out,
                           quantity = "total",
                           time.unit = time.unit,
                           scale.factor = scale.factor,
                           wb.trim = FALSE,
                           use.hinges = use.hinges)
        response <- response / total
        if (quantity == "contribution.pc") {
          response <- response * 1e2
        }
      }
    } else if (quantity %in% c("relative", "relative.pc")) {
      if (any(sapply(w.band, is_effective))) {
        warning("'quantity '", quantity,
                "' not supported when using BSWFs, returning 'total' instead")
        quantity <- "total"
      } else {
        total <- sum(response)
        response <- response / total
        if (quantity == "relative.pc") {
          response <- response * 1e2
        }
      }
    } else if (quantity %in% c("average", "mean")) {
      response <- response / sapply(w.band, wl_expanse)
    } else if (quantity != "total") {
      warning("'quantity '", quantity, "' is invalid, returning 'total' instead")
      quantity <- "total"
    }

    if (length(response) == 0) {
      response <- NA
      names(response) <- "out of range"
    } else {
      names(response) <- paste(names(response), wb_name)
    }

    if (length(scale.factor)  == 1L ||
        length(scale.factor) == length(w.band)) {
      if (any(abs(log10(scale.factor) %% 1) > 1e-5)) {
        warning("Scale factor is not decimal!")
      }
      response <- response * scale.factor
    } else {
      stop("'scale.factor' must be of length = 1 or of same length as 'w.band'.")
    }

    attr(response, "time.unit") <- getTimeUnit(spct)
    attr(response, "radiation.unit") <- paste(unit.out, "response", quantity)
    response
  }

# e_response methods --------------------------------------------------------

#' Energy-based photo-response
#'
#' This function returns the mean, total, or contribution of response for each
#' waveband and a response spectrum.
#'
#' @param spct an R object.
#' @param w.band waveband or list of waveband objects or a numeric vector of
#'   length two. The waveband(s) determine the region(s) of the spectrum that
#'   are summarized. If a numeric range is supplied a waveband object is
#'   constructed on the fly from it.
#' @param quantity character string One of "total", "average" or "mean",
#'   "contribution", "contribution.pc", "relative" or "relative.pc".
#' @param time.unit character or lubridate::duration object.
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @return A named \code{numeric} vector in the case of methods for individual
#'   spectra, with one value for each \code{waveband} passed to parameter
#'   \code{w.band}. A \code{data.frame} in the case of collections of spectra,
#'   containing one column for each \code{waveband} object, an index column with
#'   the names of the spectra, and optionally additional columns with metadata
#'   values retrieved from the attributes of the member spectra.
#'
#'   By default values are only integrated, but depending on the argument passed
#'   to parameter \code{quantity} they can be re-expressed as relative fractions
#'   or percentages. In the case of vector output, \code{names} attribute is set
#'   to the name of the corresponding waveband unless a named list is supplied
#'   in which case the names of the list members are used.
#'
#' @export
#' @examples
#' e_response(ccd.spct, new_waveband(200,300))
#' e_response(photodiode.spct)
#'
#' @note The parameter \code{use.hinges} controls speed optimization. The
#'   defaults should be suitable in most cases. Only the range of wavelengths
#'   in the wavebands is used and all BSWFs are ignored.
#'
#' @family response functions
#'
e_response <- function(spct, w.band, quantity, time.unit, scale.factor, wb.trim, use.hinges, ...) UseMethod("e_response")

#' @describeIn e_response Default method for generic function
#'
#' @export
#'
e_response.default <- function(spct, w.band, quantity, time.unit, scale.factor, wb.trim, use.hinges, ...) {
  warning("'e_response' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn e_response Method for response spectra.
#'
#' @export
#'
e_response.response_spct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL), ...) {
    resp_spct(spct = spct,
              w.band = w.band,
              unit.out = "energy",
              quantity = quantity,
              time.unit = time.unit,
              scale.factor = scale.factor,
              wb.trim = wb.trim,
              use.hinges = use.hinges )
  }

# q_response methods --------------------------------------------------------

##' Photon-based photo-response
#'
#' This function returns the mean response for a given
#' waveband and a response spectrum.
#'
#' @param spct an R object.
#' @param w.band waveband or list of waveband objects or a numeric vector of
#'   length two. The waveband(s) determine the region(s) of the spectrum that
#'   are summarized. If a numeric range is supplied a waveband object is
#'   constructed on the fly from it.
#' @param quantity character string One of "total", "average" or "mean",
#'   "contribution", "contribution.pc", "relative" or "relative.pc".
#' @param time.unit character or lubridate::duration object.
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @return A named \code{numeric} vector in the case of methods for individual
#'   spectra, with one value for each \code{waveband} passed to parameter
#'   \code{w.band}. A \code{data.frame} in the case of collections of spectra,
#'   containing one column for each \code{waveband} object, an index column with
#'   the names of the spectra, and optionally additional columns with metadata
#'   values retrieved from the attributes of the member spectra.
#'
#'   By default values are only integrated, but depending on the argument passed
#'   to parameter \code{quantity} they can be re-expressed as relative fractions
#'   or percentages. In the case of vector output, \code{names} attribute is set
#'   to the name of the corresponding waveband unless a named list is supplied
#'   in which case the names of the list members are used.
#'
#' @export
#' @examples
#' q_response(ccd.spct, new_waveband(200,300))
#' q_response(photodiode.spct)
#'
#' @note The parameter \code{use.hinges} controls speed optimization. The
#'   defaults should be suitable in most cases. Only the range of wavelengths
#'   in the wavebands is used and all BSWFs are ignored.
#'
#' @family response functions
#'
q_response <- function(spct,
                       w.band,
                       quantity,
                       time.unit,
                       scale.factor,
                       wb.trim,
                       use.hinges,
                       ...) UseMethod("q_response")

#' @describeIn q_response Default method for generic function
#'
#' @export
#'
q_response.default <- function(spct, w.band, quantity, time.unit, scale.factor, wb.trim, use.hinges, ...) {
  warning("'q_response' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn q_response Method for response spectra.
#'
#' @export
#'
q_response.response_spct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL), ... ) {
    resp_spct(spct = spct,
              w.band = w.band,
              unit.out = "photon",
              quantity = quantity,
              time.unit = time.unit,
              scale.factor = scale.factor,
              wb.trim = wb.trim,
              use.hinges = use.hinges )
  }

# response_mspct methods -----------------------------------------------

#' @describeIn response Calculates response from a \code{response_mspct}
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax for \code{attr2tb} passed as is to formal parameter \code{col.names}.
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
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
response.response_mspct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           quantity = "total",
           time.unit = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {
    z <-
      msdply(
        mspct = spct,
        .fun = response,
        w.band = w.band,
        unit.out = unit.out,
        quantity = quantity,
        time.unit = time.unit,
        scale.factor = scale.factor,
        wb.trim = wb.trim,
        use.hinges = use.hinges,
        idx = idx,
        col.names = names(w.band),
        .parallel = .parallel,
        .paropts = .paropts
      )
    add_attr2tb(tb = z,
                mspct = spct,
                col.names = attr2tb,
                idx = idx)
  }

#' @describeIn q_response Calculates photon (quantum) response from a
#'   \code{response_mspct}
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax for \code{attr2tb} passed as is to formal parameter \code{col.names}.
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
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
q_response.response_mspct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {
    z <-
      msdply(
        mspct = spct,
        .fun = q_response,
        w.band = w.band,
        quantity = quantity,
        time.unit = time.unit,
        scale.factor = scale.factor,
        wb.trim = wb.trim,
        use.hinges = use.hinges,
        idx = idx,
        col.names = names(w.band),
        .parallel = .parallel,
        .paropts = .paropts
      )
    add_attr2tb(tb = z,
                mspct = spct,
                col.names = attr2tb,
                idx = idx)
  }

#' @describeIn e_response Calculates energy response from a
#'   \code{response_mspct}
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax for \code{attr2tb} passed as is to formal parameter \code{col.names}.
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
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
e_response.response_mspct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {
    z <-
      msdply(
        mspct = spct,
        .fun = e_response,
        w.band = w.band,
        quantity = quantity,
        time.unit = time.unit,
        scale.factor = scale.factor,
        wb.trim = wb.trim,
        use.hinges = use.hinges,
        idx = idx,
        col.names = names(w.band),
        .parallel = .parallel,
        .paropts = .paropts
      )
    add_attr2tb(tb = z,
                mspct = spct,
                col.names = attr2tb,
                idx = idx)
  }
