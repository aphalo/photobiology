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
#' @param unit.out character Allowed values \code{"energy"}, and
#'   \code{"photon"}, or its alias \code{"quantum"}.
#' @param quantity character string One of \code{"average"} or \code{"mean"},
#'   \code{"total"}, \code{"contribution"}, \code{"contribution.pc"},
#'   \code{"relative"} or \code{"relative.pc"}.
#' @param time.unit character or lubridate::duration object.
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if \code{TRUE} wavebands crossing spectral data
#'   boundaries are trimmed, if \code{FALSE}, they are discarded.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param naming character one of \code{"long"}, \code{"default"},
#'   \code{"short"} or \code{"none"}. Used to select the type of names to assign
#'   to returned value.
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
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           naming = "default",
           ... ) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(response(spct = mspct,
                      w.band = w.band,
                      unit.out = unit.out,
                      quantity = quantity,
                      time.unit = time.unit,
                      scale.factor = scale.factor,
                      wb.trim = wb.trim,
                      use.hinges = use.hinges,
                      naming = naming,
                      ...))
    }

    resp_spct(spct = spct,
              w.band = w.band,
              unit.out = unit.out,
              quantity = quantity,
              time.unit = time.unit,
              scale.factor = scale.factor,
              wb.trim = wb.trim,
              use.hinges = use.hinges,
              naming = naming)
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
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
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
           allow.scaled = !quantity %in% c("average", "mean", "total"),
           naming,
           ...) {
    if (unit.out == "quantum") {
      unit.out <- "photon"
    }

    if (unit.out == "photon") {
      spct <- e2q(spct, action = "replace")
    } else if (unit.out == "energy") {
      spct <- q2e(spct, action = "replace")
    } else {
      stop("Unrecognized value", unit.out, " for unit.out")
    }

    if (quantity == "total") {
      summary.name <- switch(unit.out,
                             photon = "R[/q]",
                             energy = "R[/e]")
    } else if (quantity %in% c("average", "mean")) {
      summary.name <- switch(unit.out,
                             photon = "R(wl)[/q]",
                             energy = "R(wl)[/e]")
    } else if (quantity %in% c("contribution", "contribution.pc")) {
      summary.name <- switch(unit.out,
                             photon = "R/Rtot[/q]",
                             energy = "R/Rtot[/e]")
    } else if (quantity %in% c("relative", "relative.pc")) {
      summary.name <- switch(unit.out,
                             photon = "R/Rsum[/q]",
                             energy = "R/Rsum[/e]")
    } else {
      stop("Unrecognized 'quantity' : \"", quantity, "\"")
    }

    if (!allow.scaled && is_normalized(spct)) {
      warning("The spectral data has been normalized, making impossible to calculate response")
      return(NA_real_)
    }
    if (!allow.scaled && is_scaled(spct) && quantity == "total") {
      warning("Summary calculated from rescaled data")
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

    # "source_spct" objects are not guaranteed to contain spectral irradiance
    # expressed in t for he needed type of units.
    if (unit.out == "energy") {
      q2e(spct, byref = TRUE)
      w.length <- spct[["w.length"]]
      s.irrad <- spct[["s.e.irrad"]]
    } else if (unit.out == "photon") {
      e2q(spct, byref = TRUE)
      w.length <- spct[["w.length"]]
      s.irrad <- spct[["s.q.irrad"]]
    } else {
      stop("Unrecognized value for unit.out")
    }

    if (is.numeric(w.band)) {
      # range of wavelengths
      w.band <- waveband(w.band)
    }
    if (length(w.band) == 0) {
      # whole range of spectrum
      w.band <- waveband(spct)
    }
    if (is.waveband(w.band)) {
      # if the argument is a single w.band, we enclose it in a list
      # so that it can be handled below as a normal case.
      w.band <- list(w.band)
    }

    # convert effective wavebands
    flag.effective <- FALSE
    for (i in seq_along(w.band)) {
      if (is_effective(w.band[[i]])) {
        flag.effective <- TRUE
        w.band[[i]] <- waveband(range(w.band[[i]]))
      }
    }
    if (flag.effective) {
      warning("Using only wavelength range for effective waveband(s)")
    }
    # we trim the wavebands so that they are within the range of spct
    w.band <- trim_waveband(w.band = w.band, range = spct, trim = wb.trim)
    # if the elements of the list are named we collect them
    wb.number <- length(w.band) # number of wavebands in list
    wb.name <- names(w.band) # their names in the list
    # if no names returned, we fill the vector with "".
    if (is.null(wb.name)) {
      wb.name <- character(wb.number)
    }

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
        if (!is.null(wb[["hinges"]]) && length(wb[["hinges"]]) > 0) {
          all.hinges <- c(all.hinges, wb[["hinges"]])
        }
      }
      if (!is.null(all.hinges)) {
        spct <- insert_spct_hinges(spct, all.hinges)
      }
    }

    # we iterate through the list of wavebands
    response <- double(length(w.band))
    i <- 0
    for (wb in w.band) {
      i <- i + 1
      # we get names from wb if needed
      if (wb.name[i] == "") {
        if (naming == "short") {
          wb.name[i] <- labels(wb)[["label"]] # short name
        } else {
          wb.name[i] <- labels(wb)[["name"]] # full name
        }
      }
      # we calculate the integrated response.
      response[i] <- integrate_spct(trim_spct(spct, wb, use.hinges = FALSE))
    }
    if (quantity %in% c("contribution", "contribution.pc")) {
      total <- resp_spct(spct,
                         w.band = NULL,
                         unit.out = unit.out,
                         quantity = "total",
                         time.unit = time.unit,
                         scale.factor = scale.factor,
                         wb.trim = FALSE,
                         use.hinges = use.hinges,
                         naming = naming)
      response <- response / total
      if (quantity == "contribution.pc") {
        response <- response * 1e2
      }
    } else if (quantity %in% c("relative", "relative.pc")) {
      total <- sum(response)
      response <- response / total
      if (quantity == "relative.pc") {
        response <- response * 1e2
      }
    } else if (quantity %in% c("average", "mean")) {
      response <- response / sapply(w.band, wl_expanse)
    } else if (quantity != "total") {
      warning("'quantity '", quantity, "' is invalid, returning 'total' instead")
      quantity <- "total"
    }

    if (length(response) == 0) {
      response <- NA_real_
      names(response) <- "out of range"
    } else if (naming %in% c("long", "default")) {
      names(response) <- paste(summary.name, wb.name, sep = "_")
    } else if (naming == "short") {
      names(response) <- wb.name
    } else if (naming != "none") {
      warning("Argument to 'naming' unrecognized, assuming \"none\".")
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
    attr(response, "radiation.unit") <- paste(quantity, unit.out, "response")
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
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
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
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           naming = "default",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(e_response(spct = mspct,
                        w.band = w.band,
                        quantity = quantity,
                        time.unit = time.unit,
                        scale.factor = scale.factor,
                        wb.trim = wb.trim,
                        use.hinges = use.hinges,
                        naming = naming,
                        ...))
    }

    resp_spct(spct = spct,
              w.band = w.band,
              unit.out = "energy",
              quantity = quantity,
              time.unit = time.unit,
              scale.factor = scale.factor,
              wb.trim = wb.trim,
              use.hinges = use.hinges,
              naming = naming)
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
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
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
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           naming = "default",
           ... ) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(q_response(spct = mspct,
                       w.band = w.band,
                       quantity = quantity,
                       time.unit = time.unit,
                       scale.factor = scale.factor,
                       wb.trim = wb.trim,
                       use.hinges = use.hinges,
                       naming = naming,
                       ...))
    }

    resp_spct(spct = spct,
              w.band = w.band,
              unit.out = "photon",
              quantity = quantity,
              time.unit = time.unit,
              scale.factor = scale.factor,
              wb.trim = wb.trim,
              use.hinges = use.hinges,
              naming = naming)
  }

# response_mspct methods -----------------------------------------------

#' @describeIn response Calculates response from a \code{response_mspct}
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax for \code{attr2tb} passed as is to formal parameter \code{col.names}.
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
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
           naming = "default",
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {

    spct <- subset2mspct(spct) # expand long form spectra within collection

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
        naming = naming,
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
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
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
           naming = "default",
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {

    spct <- subset2mspct(spct) # expand long form spectra within collection

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
        naming = naming,
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
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
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
           naming = "default",
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {

    spct <- subset2mspct(spct) # expand long form spectra within collection

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
        naming = naming,
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
