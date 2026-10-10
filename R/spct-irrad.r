
# irradiance --------------------------------------------------------------

#' Irradiance
#'
#' This function returns the irradiance for a given waveband of a light source
#' spectrum.
#'
#' @param spct an R object.
#' @param w.band waveband or list of waveband objects The waveband(s) determine
#'   the region(s) of the spectrum that are summarized.
#' @param unit.out character Allowed values \code{"energy"}, and
#'   \code{"photon"}, or its alias \code{"quantum"}.
#' @param quantity character string One of "total", "average" or "mean",
#'   "contribution", "contribution.pc", "relative" or "relative.pc".
#' @param time.unit character or lubridate::duration object.
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if \code{TRUE} wavebands crossing spectral data
#'   boundaries are trimmed, if \code{FALSE}, they are discarded.
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands. If NULL, default is chosen based on data.
#' @param allow.scaled logical indicating whether scaled or normalized spectra
#'   as argument to \code{spct} trigger an error.
#' @param naming character one of \code{"long"}, \code{"default"},
#'   \code{"short"} or \code{"none"}. Used to select the type of names to assign
#'   to returned value.
#' @param return.tb logical Flag forcing a tibble to be always returned, even
#'   for a single spectrum as argumnet to \code{spct}. The default is
#'   \code{FALSE} for backwards compatibility.
#' @param ... other arguments (possibly ignored)
#'
#' @note Formal parameter \code{allow.scaled} is used internally for calculation
#'   of ratios, as rescaling and normalization do not invalidate the calculation
#'   of ratios.
#'
#' @return A named \code{numeric} vector in the case of a \code{_spct} object
#'   containing a single spectrum and \code{return.tb = FALSE}. The vector has
#'   one member one value for each \code{waveband} passed to parameter
#'   \code{w.band}. In all other cases a \code{tibble}, containing one column
#'   for each \code{waveband} object, an index column with the names of the
#'   spectra, and optionally additional columns with metadata values retrieved
#'   from the attributes of the member spectra.
#'
#'   If \code{naming = "long"} the names generated reflect both quantity and
#'   waveband, if \code{naming = "short"}, names are based only on the wavebands,
#'   and if \code{naming = "none"} the returned vector has no names.
#'
#'   By default values are only integrated, but depending on the argument passed
#'   to parameter \code{quantity} they can be re-expressed as relative fractions
#'   or percentages. In the case of vector output, \code{names} attribute is set
#'   to the name of the corresponding waveband unless a named list is supplied
#'   in which case the names of the list members are used. The \code{time.unit}
#'   attribute is copied from the spectrum object to the output. Units are as
#'   follows: If time.unit is second, [W m-2 nm-1] -> [mol s-1 m-2] or [W m-2
#'   nm-1] -> [W m-2] If time.unit is day, [J d-1 m-2 nm-1] -> [mol d-1 m-2] or
#'   [J d-1 m-2 nm-1] -> [J m-2]
#'
#' @export
#' @examples
#' irrad(sun.spct, waveband(c(400,700)))
#' irrad(sun.spct, waveband(c(400,700)), "energy")
#' irrad(sun.spct, waveband(c(400,700)), "photon")
#' irrad(sun.spct, split_bands(c(400,700), length.out = 3))
#' irrad(sun.spct, split_bands(c(400,700), length.out = 3), quantity = "total")
#' irrad(sun.spct, split_bands(c(400,700), length.out = 3), quantity = "average")
#' irrad(sun.spct, split_bands(c(400,700), length.out = 3), quantity = "relative")
#' irrad(sun.spct, split_bands(c(400,700), length.out = 3), quantity = "relative.pc")
#' irrad(sun.spct, split_bands(c(400,700), length.out = 3), quantity = "contribution")
#' irrad(sun.spct, split_bands(c(400,700), length.out = 3), quantity = "contribution.pc")
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @family irradiance functions
#'
irrad <- function(spct, w.band, unit.out, quantity, time.unit, scale.factor, wb.trim,
                  use.cached.mult, use.hinges, allow.scaled, ...) UseMethod("irrad")

#' @rdname irrad
#'
#' @export
#'
irrad.default <- function(spct, w.band, unit.out, quantity, time.unit, scale.factor, wb.trim,
                          use.cached.mult, use.hinges, allow.scaled, ...) {
  warning("'irrad' is not defined for objects of class ", class(spct)[1])
  return(NA_real_)
}

#' @rdname irrad
#'
#' @method irrad source_spct
#' @export
#'
irrad.source_spct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           quantity = "total",
           time.unit = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = NULL,
           allow.scaled = !quantity %in% c("average", "mean", "total"),
           naming = "default",
           return.tb = FALSE,
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # compute in place
      idx.var.name <- getIdFactor(spct)
      idx.var <- spct[[idx.var.name]] # not faster than inline
      if (is.factor(idx.var)) {
        idx.levels <- levels(idx.var)
      } else {
        idx.levels <- unique(idx.var)
      }
      # do conversion in one go and delete values not used
      if (unit.out == "energy") {
        spct <- q2e(spct, action = "replace")
      } else {
        spct <- e2q(spct, action = "replace")
      }
      z <- list()
      for (idx in idx.levels) {
        target.rows <- which(idx.var == idx) # a lot faster than logical
        temp.spct <- spct[target.rows, ]
        z[[idx]] <- irrad_spct(spct = temp.spct,
                               w.band = w.band,
                               unit.out = unit.out,
                               quantity = quantity,
                               time.unit = time.unit,
                               scale.factor = scale.factor,
                               wb.trim = wb.trim,
                               use.cached.mult = use.cached.mult,
                               use.hinges = use.hinges,
                               allow.scaled = allow.scaled,
                               naming = naming,
                               return.tb = TRUE,
                               ...)
      }
      z <- dplyr::bind_rows(z)
      z[[idx.var.name]] <- idx.levels
      z[["when.measured"]] <-
        as.POSIXct(unlist(when_measured(spct), use.names = FALSE),
                   tz = "UTC",
                   origin = lubridate::origin)
      attr(z, "time.unit") <- getTimeUnit(spct)
      if (is_effective(spct)) {
        attr(z, "radiation.unit") <-
          paste(unit.out, "irradiance", quantity, "effective:", getBSWFUsed(spct))
      } else {
        attr(z, "radiation.unit") <- paste(quantity, unit.out, "irradiance")
      }

      return(z)
    }

    if (unit.out == "quantum") {
      unit.out <- "photon"
    }

    if (quantity == "total") {
      summary.name <- switch(unit.out,
                             photon = "Q",
                             energy = "E")
    } else if (quantity %in% c("average", "mean")) {
      summary.name <- switch(unit.out,
                             photon = "Q(wl)",
                             energy = "E(wl)")
    } else if (quantity %in% c("contribution", "contribution.pc")) {
      summary.name <- switch(unit.out,
                             photon = "Q/Qtot",
                             energy = "E/Etot")
    } else if (quantity %in% c("relative", "relative.pc")) {
      summary.name <- switch(unit.out,
                             photon = "Q/Qsum",
                             energy = "E/Esum")
    } else {
      stop("Unrecognized 'quantity' : \"", quantity, "\"")
    }

    if (!allow.scaled && is_normalized(spct)) {
      warning("The spectral data have been normalized, ",
              "preventing calculation of irradiance. ",
              "'allow.scaled = TRUE' disables this test.")
      return(NA_real_)
    }
    if (!allow.scaled && is_scaled(spct)) {
      warning("The spectral data have been scaled, ",
              "preventing calculation of irradiance. ",
              "'allow.scaled = TRUE' disables this test.")
      return(NA_real_)
    }

    data.time.unit <-
      getTimeUnit(spct, force.duration = lubridate::is.duration(time.unit))

    if (!is.null(time.unit) && time.unit != data.time.unit) {
      if (!lubridate::is.duration(time.unit) && !is.character(time.unit)) {
        message("converting 'time.unit' ", time.unit, " into a lubridate::duration")
        time.unit <- lubridate::as.duration(time.unit)
      }
      spct <- convertTimeUnit(spct, time.unit = time.unit, byref = FALSE)
    } else {
      time.unit <- data.time.unit
    }

    if (is.null(unit.out) || is.na(unit.out)) {
      warning("'unit.out' set to an invalid value")
      return(NA_real_)
    }

    # "source_spct" objects are not guaranteed to contain spectral irradiance
    # expressed in the needed units.
    if (unit.out == "energy") {
      q2e(spct, byref = TRUE)
      w.length <- spct[["w.length"]]
      s.irrad <- spct[["s.e.irrad"]]
    } else if (unit.out == "photon") {
      e2q(spct, byref = TRUE)
      w.length <- spct[["w.length"]]
      s.irrad <- spct[["s.q.irrad"]]
    } else {
      stop("Unrecognized value", unit.out, " for unit.out")
    }

    if (length(w.band) == 0) {
      # whole range of spectrum
      w.band <- waveband(spct)
    }
    if (is.numeric(w.band)) {
      # range of wavelengths
      w.band <- waveband(w.band)
    }
    if (is.waveband(w.band)) {
      # if the argument is a single w.band, we enclose it in a list
      # so that it can be handled below as a normal case.
      w.band <- list(w.band)
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

    # hinges
    if (is.null(use.hinges)) {
       use.hinges <- auto_hinges(w.length)
    }
    # we collect all hinges and insert them in one go
    if (use.hinges) {
      all.hinges <- NULL
      for (wb in w.band) {
        all.hinges <- c(all.hinges, wb[["hinges"]])
      }
      if (anyNA(all.hinges)) {
        warning("Missing hinges encountered and removed!")
        all.hinges <- na.omit(all.hinges)
      }
      lst <- l_insert_hinges(w.length, s.irrad, all.hinges)
      w.length <- lst[["x"]]
      s.irrad <- lst[["y"]]
    }

    # We iterate through the list of wavebands collecting the irradiances,
    # and waveband names.
    irrad <- numeric(wb.number)
    i <- 0L
    is.effective.spectrum <- is_effective(spct)
    for (wb in w.band) {
      i <- i + 1L
      # get names from wb if needed
      if (wb.name[i] == "") {
        if (naming == "short") {
          wb.name[i] <- labels(wb)[["label"]] # short name
        } else {
          wb.name[i] <- labels(wb)[["name"]] # full name
        }
      }
      # check for NA wavebands
      if (is.na(wb[["low"]]) || is.na(wb[["high"]])) {
        irrad[i] <- NA_real_
        next
      }
      if (is.effective.spectrum && is_effective(wb)) {
        warning("Effective spectral irradiance is not compatible with a BSWF: ",
                wb.name[i])
        irrad[i] <- NA_real_
      } else {
        if (is.effective.spectrum) {
          wb.name[i] <- paste(getBSWFUsed(spct), "*", wb.name[i], sep = "")
        }
        wl.selector <- which(w.length >= min(wb) & w.length <= max(wb))
        if (wl.selector[1] > 1) {
          wl.selector <- c(wl.selector[1] - 1, wl.selector)
        }
        if (wl.selector[length(wl.selector)] < length(w.length)) {
          wl.selector <- c(wl.selector,  wl.selector[length(wl.selector)] + 1)
        }

        # calculate the multipliers
        mult <- calc_multipliers(w.length = w.length[wl.selector],
                                 w.band = wb,
                                 unit.out = unit.out,
                                 unit.in = unit.out,
                                 use.cached.mult = use.cached.mult)
        # calculate weighted spectral irradiance
        # the ifelse is needed to override NAs in spectral data for regions
        # where mult == 0
          irrad[i] <- integrate_xy(w.length[wl.selector],
                                   ifelse(mult == 0, 0, s.irrad[wl.selector] * mult))
      }
    }

    if (quantity %in% c("contribution", "contribution.pc")) {
      if (any(sapply(w.band, is_effective))) {
        warning("'quantity '", quantity,
                "' not supported when using BSWFs, returning 'total' instead")
        quantity <- "total"
      } else {
        # recursive call
        total <- irrad_spct(spct, w.band = NULL,
                            unit.out = unit.out,
                            quantity = "total",
                            time.unit = time.unit,
                            use.cached.mult = use.cached.mult,
                            wb.trim = wb.trim,
                            use.hinges = use.hinges,
                            naming = naming)
        irrad <- irrad / total
        if (quantity == "contribution.pc") {
          irrad <- irrad * 1e2
        }
      }
    } else if (quantity %in% c("relative", "relative.pc")) {
      if (any(sapply(w.band, is_effective))) {
        warning("'quantity '", quantity,
                "' not supported when using BSWFs, returning 'total' instead")
        quantity <- "total"
      } else {
        total <- sum(irrad)
        irrad <- irrad / total
        if (quantity == "relative.pc") {
          irrad <- irrad * 1e2
        }
      }
    } else if (quantity %in% c("average", "mean") ) {
      irrad <- irrad / sapply(w.band, wl_expanse)
    } else if (quantity != "total") {
      warning("'quantity '", quantity, "' is invalid, returning 'total' instead")
      quantity <- "total"
    }

    if (length(irrad) == 0) {
      irrad <- NA_real_
      names(irrad) <- "out of range or NAs in waveband"
    } else if (naming %in% c("long", "default")) {
      names(irrad) <- paste(summary.name, wb.name, sep = "_")
    } else if (naming == "short") {
      names(irrad) <- wb.name
    } else if (naming != "none") {
      warning("Argument to 'naming' unrecognized, assuming \"none\".")
    }

    if (length(scale.factor)  == 1L ||
        length(scale.factor) == length(w.band)) {
      if (any(abs(log10(scale.factor) %% 1) > 1e-5)) {
        warning("Scale factor is not decimal!")
      }
      irrad <- irrad * scale.factor
    } else {
      stop("'scale.factor' must be of length = 1 or of same length as 'w.band'.")
    }

    if (return.tb) {
      irrad <- tibble::as_tibble_row(irrad, .name_repair = "minimal")
    } else {
      attr(irrad, "time.unit") <- getTimeUnit(spct)
      if (is_effective(spct)) {
        attr(irrad, "radiation.unit") <-
          paste(unit.out, "irradiance", quantity, "effective:", getBSWFUsed(spct))
      } else {
        attr(irrad, "radiation.unit") <- paste(quantity, unit.out, "irradiance")
      }
    }

    irrad
  }

#' @keywords internal
irrad_spct <- irrad.source_spct

# energy irradiance -------------------------------------------------------

#' Energy irradiance
#'
#' Energy irradiance for one or more wavebands of a light source spectrum.
#'
#' @param spct an R object.
#' @param w.band a list of \code{waveband} objects or a \code{waveband} object.
#' @param quantity character string One of "total", "average" or "mean",
#'   "contribution", "contribution.pc", "relative" or "relative.pc".
#' @param time.unit character or lubridate::duration object.
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded.
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param allow.scaled logical indicating whether scaled or normalized spectra
#'   as argument to spct are flagged as an error.
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param return.tb logical Flag forcing a tibble to be always returned, even
#'   for a single spectrum as argumnet to \code{spct}. The default is
#'   \code{FALSE} for backwards compatibility.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @export
#'
#' @examples
#' e_irrad(sun.spct, waveband(c(400,700)))
#' e_irrad(sun.spct, split_bands(c(400,700), length.out = 3))
#' e_irrad(sun.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "total")
#' e_irrad(sun.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "average")
#' e_irrad(sun.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "relative")
#' e_irrad(sun.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "relative.pc")
#' e_irrad(sun.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "contribution")
#' e_irrad(sun.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "contribution.pc")
#'
#' @return A named \code{numeric} vector in the case of a \code{_spct} object
#'   containing a single spectrum and \code{return.tb = FALSE}. The vector has
#'   one member one value for each \code{waveband} passed to parameter
#'   \code{w.band}. In all other cases a \code{tibble}, containing one column
#'   for each \code{waveband} object, an index column with the names of the
#'   spectra, and optionally additional columns with metadata values retrieved
#'   from the attributes of the member spectra.
#'
#'   By default values are only integrated, but depending on the argument passed
#'   to parameter \code{quantity} they can be re-expressed as relative fractions
#'   or percentages. In the case of vector output, \code{names} attribute is set
#'   to the name of the corresponding waveband unless a named list is supplied
#'   in which case the names of the list members are used. The time.unit
#'   attribute is copied from the spectrum object to the output. Units are as
#'   follows: If units are absolute and time.unit is second, [W m-2 nm-1] -> [W
#'   m-2] If time.unit is day, [J d-1 m-2 nm-1] -> [J m-2]; if units are
#'   relative, fraction of one or percent.
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @family irradiance functions
#'
e_irrad <- function(spct, w.band,
                    quantity, time.unit, scale.factor, wb.trim,
                    use.cached.mult, use.hinges, allow.scaled,
                    ...) UseMethod("e_irrad")

#' @rdname e_irrad
#'
#' @export
#'
e_irrad.default <- function(spct, w.band,
                            quantity, time.unit, scale.factor, wb.trim,
                            use.cached.mult, use.hinges, allow.scaled, ...) {
  warning("'e_irrad' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @rdname e_irrad
#'
#' @export
#'
e_irrad.source_spct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = NULL,
           allow.scaled = !quantity  %in% c("average", "mean", "total"),
           naming = "default",
           return.tb = FALSE,
           ...) {

    irrad_spct(spct, w.band = w.band, unit.out = "energy",
               scale.factor = scale.factor,
               quantity = quantity,
               time.unit = time.unit, wb.trim = wb.trim,
               use.cached.mult = use.cached.mult, use.hinges = use.hinges,
               allow.scaled = allow.scaled,
               naming = naming,
               return.tb = return.tb)
  }

# photon irradiance -------------------------------------------------------

#' Photon irradiance
#'
#' Photon irradiance (i.e. quantum irradiance) for one or more wavebands of a
#' light source spectrum.
#'
#' @param spct an R object.
#' @param w.band a list of \code{waveband} objects or a \code{waveband} object.
#' @param quantity character string One of "total", "average" or "mean",
#'   "contribution", "contribution.pc", "relative" or "relative.pc".
#' @param time.unit character or lubridate::duration object.
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded.
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param allow.scaled logical indicating whether scaled or normalized spectra
#'   as argument to spct are flagged as an error.
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param return.tb logical Flag forcing a tibble to be always returned, even
#'   for a single spectrum as argumnet to \code{spct}. The default is
#'   \code{FALSE} for backwards compatibility.
#' @param ... other arguments (possibly ignored).
#'
#' @export
#'
#' @examples
#' q_irrad(sun.spct, waveband(c(400,700)))
#' q_irrad(sun.spct, split_bands(c(400,700), length.out = 3))
#' q_irrad(sun.spct, split_bands(c(400,700), length.out = 3), quantity = "total")
#' q_irrad(sun.spct, split_bands(c(400,700), length.out = 3), quantity = "average")
#' q_irrad(sun.spct, split_bands(c(400,700), length.out = 3), quantity = "relative")
#' q_irrad(sun.spct, split_bands(c(400,700), length.out = 3), quantity = "relative.pc")
#' q_irrad(sun.spct, split_bands(c(400,700), length.out = 3), quantity = "contribution")
#' q_irrad(sun.spct, split_bands(c(400,700), length.out = 3), quantity = "contribution.pc")
#'
#' @return A named \code{numeric} vector in the case of a \code{_spct} object
#'   containing a single spectrum and \code{return.tb = FALSE}. The vector has
#'   one member one value for each \code{waveband} passed to parameter
#'   \code{w.band}. In all other cases a \code{tibble}, containing one column
#'   for each \code{waveband} object, an index column with the names of the
#'   spectra, and optionally additional columns with metadata values retrieved
#'   from the attributes of the member spectra.
#'
#'   By default values are only integrated, but depending on the argument passed
#'   to parameter \code{quantity} they can be re-expressed as relative fractions
#'   or percentages. In the case of vector output, \code{names} attribute is set
#'   to the name of the corresponding waveband unless a named list is supplied
#'   in which case the names of the list members are used. The time.unit
#'   attribute is copied from the spectrum object to the output. Units are as
#'   follows: If time.unit is second, [W m-2 nm-1] -> [mol s-1 m-2] If time.unit
#'   is day, [J d-1 m-2 nm-1] -> [mol d-1 m-2]
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @export
#' @family irradiance functions
q_irrad <- function(spct, w.band,
                    quantity, time.unit, scale.factor, wb.trim,
                    use.cached.mult, use.hinges, allow.scaled, ...) UseMethod("q_irrad")

#' @rdname q_irrad
#'
#' @export
#'
q_irrad.default <- function(spct, w.band,
                            quantity, time.unit, scale.factor, wb.trim,
                            use.cached.mult, use.hinges, allow.scaled, ...) {
  warning("'q_irrad' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @rdname q_irrad
#'
#' @export
#'
q_irrad.source_spct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = NULL,
           allow.scaled = !quantity  %in% c("average", "mean", "total"),
           naming = "default",
           return.tb = FALSE,
           ...) {

    irrad_spct(spct, w.band = w.band, unit.out = "photon", quantity = quantity,
               time.unit = time.unit,
               scale.factor = scale.factor,
               wb.trim = wb.trim,
               use.cached.mult = use.cached.mult, use.hinges = use.hinges,
               allow.scaled = allow.scaled,
               naming = naming,
               return.tb = return.tb)
  }


# fluence -----------------------------------------------------------------

#' Fluence
#'
#' Energy or photon fluence for one or more wavebands of a light source spectrum
#' and a duration of exposure.
#'
#' @param spct an R object.
#' @param w.band a list of \code{waveband} objects or a \code{waveband} object.
#' @param unit.out character string with allowed values "energy", and "photon",
#'   or its alias "quantum".
#' @param exposure.time lubridate::duration object.
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded.
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param allow.scaled logical indicating whether scaled or normalized spectra
#'   as argument to spct are flagged as an error.
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @export
#'
#' @examples
#' library(lubridate)
#' fluence(sun.spct,
#'         w.band = waveband(c(400,700)),
#'         exposure.time = lubridate::duration(3, "minutes") )
#'
#' @return One numeric value for each waveband with no change in scale factor,
#'   with name attribute set to the name of each waveband unless a named list is
#'   supplied in which case the names of the list elements are used. The
#'   time.unit attribute is copied from the spectrum object to the output. Units
#'   are as follows: If time.unit is second, [W m-2 nm-1] -> [mol s-1 m-2] If
#'   time.unit is day, [J d-1 m-2 nm-1] -> [mol d-1 m-2]
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @export
#' @family irradiance functions
fluence <- function(spct, w.band, unit.out, exposure.time, scale.factor, wb.trim,
                    use.cached.mult, use.hinges, allow.scaled, ...) UseMethod("fluence")

#' @rdname fluence
#'
#' @export
#'
fluence.default <- function(spct, w.band, unit.out, exposure.time, scale.factor,
                            wb.trim, use.cached.mult, use.hinges, allow.scaled, ...) {
  warning("'fluence' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @rdname fluence
#'
#' @export
#'
fluence.source_spct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           exposure.time,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = NULL,
           allow.scaled = FALSE,
           naming = "default",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(fluence(spct = mspct,
                     w.band = w.band,
                     unit.out = unit.out,
                     exposure.time = exposure.time,
                     scale.factor = scale.factor,
                     wb.trim = wb.trim,
                     use.cached.mult = use.cached.mult,
                     use.hinges = use.hinges,
                     allow.scaled = allow.scaled,
                     naming = naming,
                     ...))
    }

    if (!lubridate::is.duration(exposure.time) &&
        !lubridate::is.period(exposure.time) &&
        !is.numeric(exposure.time) ) {
      stop("Invalid value ", exposure.time, " for 'exposure.time'")
    } else if (is.na(exposure.time)) {
      return.value <- NA_real_
    } else {
      return.value <-
        irrad_spct(spct, w.band = w.band, unit.out = unit.out, quantity = "total",
                   time.unit = exposure.time,
                   scale.factor = scale.factor,
                   wb.trim = wb.trim,
                   use.cached.mult = use.cached.mult, use.hinges = use.hinges,
                   allow.scaled = allow.scaled,
                   naming = naming)
    }
    if (unit.out %in% c("photon", "quantum")) {
      attr(return.value, "radiation.unit") <- "photon fluence (mol m-2)"
    } else if (unit.out == "energy") {
      attr(return.value, "radiation.unit") <- "energy fluence (J m-2)"
    }
    attr(return.value, "exposure.duration") <- exposure.time
    attr(return.value, "time.unit") <- NULL
    return.value
  }


# photon fluence ----------------------------------------------------------

#' Photon fluence
#'
#' Photon irradiance (i.e. quantum irradiance) for one or more waveband of a
#' light source spectrum.
#'
#' @param spct an R object.
#' @param w.band a list of \code{waveband} objects or a \code{waveband} object
#' @param exposure.time lubridate::duration object.
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded.
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param allow.scaled logical indicating whether scaled or normalized spectra
#'   as argument to spct are flagged as an error.
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param ... other arguments (possibly ignored).
#'
#' @examples
#' library(lubridate)
#' q_fluence(sun.spct,
#'           w.band = waveband(c(400,700)),
#'           exposure.time = lubridate::duration(3, "minutes") )
#'
#' @return One numeric value for each waveband with no change in scale factor,
#'   with name attribute set to the name of each waveband unless a named list is
#'   supplied in which case the names of the list elements are used. The
#'   exposure.time is copied from the spectrum object to the output as an attribute.
#'   Units are as follows: moles of photons per exposure.
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @export
#'
#' @family irradiance functions
q_fluence <- function(spct, w.band, exposure.time, scale.factor, wb.trim, use.cached.mult,
                      use.hinges, allow.scaled, ...) UseMethod("q_fluence")

#' @rdname q_fluence
#'
#' @export
#'
q_fluence.default <- function(spct, w.band, exposure.time, scale.factor, wb.trim,
                              use.cached.mult, use.hinges, allow.scaled, ...) {
  warning("'q_fluence' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @rdname q_fluence
#'
#' @export
#'
q_fluence.source_spct <-
  function(spct, w.band = NULL,
           exposure.time,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = NULL,
           allow.scaled = FALSE,
           naming = "default",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(q_fluence(spct = mspct,
                       w.band = w.band,
                       exposure.time = exposure.time,
                       scale.factor = scale.factor,
                       wb.trim = wb.trim,
                       use.cached.mult = use.cached.mult,
                       use.hinges = use.hinges,
                       allow.scaled = allow.scaled,
                       naming = naming,
                       ...))
    }

    if (!lubridate::is.duration(exposure.time) &&
        !lubridate::is.period(exposure.time) &&
        !is.numeric(exposure.time) ) {
      stop("Invalid value ", exposure.time, " for 'exposure.time'")
    } else if (is.na(exposure.time)) {
      return.value <- NA_real_
    } else {
      return.value <-
        irrad_spct(spct, w.band = w.band, unit.out = "photon", quantity = "total",
                   time.unit = exposure.time,
                   scale.factor = scale.factor,
                   wb.trim = wb.trim,
                   use.cached.mult = use.cached.mult, use.hinges = use.hinges,
                   allow.scaled = allow.scaled,
                   naming = naming)
    }
    attr(return.value, "radiation.unit") <- "photon fluence (mol m-2)"
    attr(return.value, "exposure.duration") <- exposure.time
    attr(return.value, "time.unit") <- NULL
    return.value
  }


# energy fluence ----------------------------------------------------------

#' Energy fluence
#'
#' Energy fluence for one or more wavebands of a light source spectrum and a
#' duration of the exposure.
#'
#' @param spct an R object
#' @param w.band a list of \code{waveband} objects or a \code{waveband} object
#' @param exposure.time lubridate::duration object.
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param allow.scaled logical indicating whether scaled or normalized spectra
#'   as argument to spct are flagged as an error
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param ... other arguments (possibly ignored)
#'
#' @examples
#' library(lubridate)
#' e_fluence(sun.spct, w.band = waveband(c(400,700)),
#'           exposure.time = lubridate::duration(3, "minutes") )
#'
#' @return One numeric value for each waveband with no change in scale factor,
#'   with name attribute set to the name of each waveband unless a named list is
#'   supplied in which case the names of the list elements are used. The
#'   exposure.time is copied to the output as an attribute. Units are as
#'   follows: (J) joules per exposure.
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @export
#' @family irradiance functions
e_fluence <- function(spct, w.band, exposure.time, scale.factor, wb.trim, use.cached.mult,
                      use.hinges, allow.scaled, ...) UseMethod("e_fluence")

#' @rdname e_fluence
#'
#' @export
#'
e_fluence.default <- function(spct, w.band, exposure.time, scale.factor, wb.trim, use.cached.mult,
                              use.hinges, allow.scaled, ...) {
  warning("'e_fluence' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @rdname e_fluence
#'
#' @export
#'
e_fluence.source_spct <-
  function(spct, w.band = NULL,
           exposure.time,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = NULL,
           allow.scaled = FALSE,
           naming = "default",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(q_fluence(spct = mspct,
                       w.band = w.band,
                       exposure.time = exposure.time,
                       scale.factor = scale.factor,
                       wb.trim = wb.trim,
                       use.cached.mult = use.cached.mult,
                       use.hinges = use.hinges,
                       allow.scaled = allow.scaled,
                       naming = naming,
                       ...))
    }

    if (!lubridate::is.duration(exposure.time) &&
        !lubridate::is.period(exposure.time) &&
        !is.numeric(exposure.time) ) {
      stop("Invalid value ", exposure.time, " for 'exposure.time'")
    } else if (is.na(exposure.time)) {
      return.value <- NA_real_
    } else {
      return.value <-
        irrad_spct(spct, w.band = w.band, unit.out = "energy", quantity = "total",
                   time.unit = exposure.time,
                   scale.factor = scale.factor,
                   wb.trim = wb.trim,
                   use.cached.mult = use.cached.mult, use.hinges = use.hinges,
                   allow.scaled = allow.scaled,
                   naming = naming)
    }
    attr(return.value, "radiation.unit") <- "energy fluence (J m-2)"
    attr(return.value, "exposure.duration") <- exposure.time
    attr(return.value, "time.unit") <- NULL
    return.value
  }

# source_mspct methods -----------------------------------------------

#' @rdname irrad
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax
#'   for \code{attr2tb} passed as is to formal parameter \code{col.names}.
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach.
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
irrad.source_mspct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           quantity = "total",
           time.unit = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = NULL,
           allow.scaled = !quantity  %in% c("average", "mean", "total"),
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
        .fun = irrad.source_spct,
        w.band = w.band,
        quantity = quantity,
        time.unit,
        unit.out = unit.out,
        scale.factor = scale.factor,
        wb.trim = wb.trim,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
        allow.scaled = allow.scaled,
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

#' @rdname q_irrad
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
q_irrad.source_mspct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = NULL,
           allow.scaled = !quantity  %in% c("average", "mean", "total"),
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
        .fun = q_irrad.source_spct,
        w.band = w.band,
        quantity = quantity,
        time.unit = time.unit,
        scale.factor = scale.factor,
        wb.trim = wb.trim,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
        allow.scaled = allow.scaled,
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

#' @rdname e_irrad
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
e_irrad.source_mspct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = NULL,
           allow.scaled = !quantity  %in% c("average", "mean", "total"),
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
        .fun = e_irrad.source_spct,
        w.band = w.band,
        quantity = quantity,
        time.unit = time.unit,
        scale.factor = scale.factor,
        wb.trim = wb.trim,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
        allow.scaled = allow.scaled,
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

#' @rdname fluence
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
fluence.source_mspct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           exposure.time,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = NULL,
           allow.scaled = FALSE,
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
        .fun = fluence.source_spct,
        w.band = w.band,
        unit.out = unit.out,
        exposure.time = exposure.time,
        scale.factor = scale.factor,
        wb.trim = wb.trim,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
        allow.scaled = allow.scaled,
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

#' @rdname e_fluence
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
e_fluence.source_mspct <-
  function(spct, w.band = NULL,
           exposure.time,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = NULL,
           allow.scaled = FALSE,
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {

    spct <- subset2mspct(spct) # expand long form spectra within collection

    z <-
      msdply(
        mspct = spct,
        .fun = e_fluence.source_spct,
        w.band = w.band,
        exposure.time = exposure.time,
        scale.factor = scale.factor,
        wb.trim = wb.trim,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
        allow.scaled = allow.scaled,
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

#' @rdname q_fluence
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
q_fluence.source_mspct <-
  function(spct, w.band = NULL,
           exposure.time,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges = NULL,
           allow.scaled = FALSE,
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
        .fun = q_fluence.source_spct,
        w.band = w.band,
        exposure.time = exposure.time,
        scale.factor = scale.factor,
        wb.trim = wb.trim,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
        allow.scaled = allow.scaled,
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
