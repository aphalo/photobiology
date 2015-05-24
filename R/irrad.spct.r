
# irradiance --------------------------------------------------------------

#' Calculate irradiance from spectral irradiance.
#'
#' This function returns the irradiance for a given waveband of a light source
#' spectrum.
#'
#' @param spct an R object
#' @param w.band waveband or list of waveband objects The waveband(s) determine
#'   the region(s) of the spectrum that are summarized.
#' @param unit.out character string with allowed values "energy", and "photon",
#'   or its alias "quantum"
#' @param quantity character string
#' @param time.unit character or lubridate::duration
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#' @param allow.scaled logical indicating whether scaled or normalized spectra
#'   as argument to spct are flagged as an error
#'
#' @note Formal parameter \code{allow.scaled} is used internally for calculation
#'   of ratios, as rescaling and normalization do not invalidate the calculation
#'   of ratios.
#'
#' @return One numeric value for each waveband with no change in scale factor,
#'   with name attribute set to the name of each waveband unless a named list is
#'   supplied in which case the names of the list elements are used. The
#'   time.unit attribute is copied from the spectrum object to the output. Units
#'   are as follows: If time.unit is second, [W m-2 nm-1] -> [mol s-1 m-2] or [W
#'   m-2 nm-1] -> [W m-2] If time.unit is day, [J d-1 m-2 nm-1] -> [mol d-1 m-2]
#'   or [J d-1 m-2 nm-1] -> [J m-2]
#'
#' @keywords manip misc
#' @export
#' @examples
#' irrad(sun.spct, new_waveband(400,700), "photon")
#' irrad(sun.spct, new_waveband(400,700), "energy")
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in mosts cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @aliases irrad
#'
#' @family irradiance functions
#'
irrad <- function(spct, w.band, unit.out, quantity, time.unit, wb.trim,
                  use.cached.mult, use.hinges, allow.scaled) UseMethod("irrad")

#' @describeIn irrad Default for generic function
#'
#' @export
#'
irrad.default <- function(spct, w.band, unit.out, quantity, time.unit, wb.trim,
                          use.cached.mult, use.hinges, allow.scaled) {
  warning("'irrad' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn irrad  Calculates irradiance from a \code{source_spct}
#'   object.
#'
#' @method irrad source_spct
#' @export
#'
irrad.source_spct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           quantity="total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL),
           allow.scaled = FALSE){
    # we have a default, but we check for invalid arguments
    if (!allow.scaled && (is_normalized(spct) || is_scaled(spct))) {
      warning("The espectral data has been normalized or scaled, making impossible to calculate irradiance")
      return(NA)
    }

    if (!is.null(time.unit)) {
      spct <- convertTimeUnit(spct, time.unit = time.unit, byref = FALSE)
    } else {
      time.unit <- getTimeUnit(spct)
    }

    if (identical(attr(spct, ".data.table.locked"), TRUE)) {
      spct_x <- copy(spct)
    } else {
      spct_x <- spct
    }
    if (is.null(unit.out) || is.na(unit.out)){
      warning("'unit.out' set to an invalid value")
      return(NA)
    }
    if (unit.out == "quantum") {
      unit.out <- "photon"
    }
    if (is.null(w.band)) {
      w.band <- waveband(spct_x)
    }
    if (is.waveband(w.band)) {
      # if the argument is a single w.band, we enclose it in a list
      # so that the for loop works as expected.This is a bit of a
      # cludge but lets us avoid treating it as a special case
      w.band <- list(w.band)
    }
    w.band <- trim_waveband(w.band=w.band, range=spct_x, trim=wb.trim)
    # we check if the list elements are named, if not we set a flag
    # and an empty vector that will be later filled in with data from
    # the waveband definitions.
    wb.number <- length(w.band) # number of wavebands in list
    wb.name <- names(w.band) # their names in the list
    if (is.null(wb.name)) {
      wb.name <- character(wb.number)
    }
    # if the w.band includes 'hinges' we insert them
    # choose whether to use hinges or not
    # if the user has specified its value, we leave it alone
    # but if it was not requested, we decide whether to use
    # it or not based of the wavelength resolution of the
    # spectrum. This will produce small errors for high
    # spectral resolution data, and speed up the calculations
    # a lot in such cases
    if (is.null(use.hinges)) {
      use.hinges <- stepsize(spct_x)[2] > getOption("photobiology.auto.hinges.limit", default = 0.5) # nm
    }
    # we collect all hinges and insert them in one go
    # this may alter a little the returned values
    # but should be faster
    if (use.hinges) {
      all.hinges <- NULL
      for (wb in w.band) {
        if (!is.null(wb$hinges) & length(wb$hinges)>0) {
          all.hinges <- c(all.hinges, wb$hinges)
        }
      }
      if (!is.null(all.hinges)) {
        spct_x <- insert_spct_hinges(spct_x, all.hinges)
      }
    }

    # "source_spct" objects are not guaranteed to contain spectral irradiance
    # expressed in the needed type of scale, if the needed one is missing
    # we add the missing it.
    # As spectra are passed by reference the changes propagate to the argument
    if (unit.out == "energy") {
      q2e(spct_x, byref=TRUE)
    } else if (unit.out == "photon") {
      e2q(spct_x, byref=TRUE)
    } else {
      stop("Unrecognized value for unit.out")
    }
    unit.in <- unit.out

    # We iterate through the list of wavebands collecting the integrated irradiances,
    # possibly weighted depending on the waveband definition
    irrad <- numeric(wb.number)
    i <- 0L
    is.effective.spectrum <- is_effective(spct)
    for (wb in w.band) {
      i <- i + 1L
      # get names from wb if needed
      if (wb.name[i] == "") {
        wb.name[i] <- wb$name
      }
      if (is.effective.spectrum && is_effective(wb)) {
        warning("Effective spectral irradiance is not compatible with a BSWF: ", wb.name[i])
        irrad[i] <- NA
      } else {
        if (is.effective.spectrum) {
          wb.name[i] <- paste(getBSWFUsed(spct), "*", wb.name[i])
        }
        # calculate the multipliers
        mult <- calc_multipliers(w.length=spct_x$w.length, w.band=wb, unit.out=unit.out,
                                 unit.in=unit.in, use.cached.mult=use.cached.mult)
        # calculate weighted spectral irradiance
        # the ifelse is needed to overrride NAs in spectral data for regions where mult == 0
        if (unit.out == "energy") {
          irr <- with(spct_x, integrate_irradiance(w.length, ifelse(mult == 0, 0, s.e.irrad * mult)))
        } else {
          irr <- with(spct_x, integrate_irradiance(w.length, ifelse(mult == 0, 0, s.q.irrad * mult)))
        }
        irrad[i] <- irr
      }
    }
    if (quantity %in% c("contribution", "contribution.pc")) {
      if (any(sapply(w.band, is_effective))) {
        warning("'quantity '", quantity, "' not supported when using BSWFs, returning 'total' instead")
        quantity <- "total"
      } else {
        total <- irrad_spct(spct_x, w.band=NULL, unit.out=unit.out,
                            quantity="total", use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE), use.hinges=FALSE)
        irrad <- irrad / total
        if (quantity == "contribution.pc") {
          irrad <- irrad * 1e2
        }
      }
    } else if (quantity %in% c("relative", "relative.pc")) {
      if (any(sapply(w.band, is_effective))) {
        warning("'quantity '", quantity, "' not supported when using BSWFs, returning 'total' instead")
        quantity <- "total"
      } else {
        total <- sum(irrad)
        irrad <- irrad / total
        if (quantity == "relative.pc") {
          irrad <- irrad * 1e2
        }
      }
    } else if (quantity == "average") {
      irrad <- irrad / sapply(w.band, spread)
    } else if (quantity != "total") {
      warning("'quantity '", quantity, "' is invalid, returning 'total' instead")
      quantity <- "total"
    }
    if (length(irrad) == 0) {
      irrad <- NA
      names(irrad) <- "out of range"
    }
    names(irrad) <- paste(names(irrad), wb.name)
    setattr(irrad, "time.unit", getTimeUnit(spct_x))
    if (is_effective(spct_x)) {
      setattr(irrad, "radiation.unit",
              paste(unit.out, "irradiance", quantity, "effective:", getBSWFUsed(spct_x)))
    } else {
      setattr(irrad, "radiation.unit", paste(unit.out, "irradiance", quantity))
    }
    return(irrad)
  }

#' @keywords internal
irrad_spct <- irrad.source_spct

# energy irradiance -------------------------------------------------------


#' Calculate energy irradiance from spectral irradiance.
#'
#' This function returns the energy irradiance for a given waveband of a light
#' source spectrum.
#'
#' @param spct an R object
#' @param w.band a list of \code{waveband} objects or a \code{waveband} object
#' @param quantity character string
#' @param time.unit character or lubridate::duration
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#' @param allow.scaled logical indicating whether scaled or normalized spectra
#'   as argument to spct are flagged as an error
#'
#' @keywords manip misc
#'
#' @export
#'
#' @examples
#' e_irrad(sun.spct, new_waveband(400,700))
#'
#' @return One numeric value for each waveband with no change in scale factor,
#'   with name attribute set to the name of each waveband unless a named list is
#'   supplied in which case the names of the list elements are used. The
#'   time.unit attribute is copied from the spectrum object to the output. Units
#'   are as follows: If time.unit is second, [W m-2 nm-1] -> [W m-2] If
#'   time.unit is day, [J d-1 m-2 nm-1] -> [J m-2]
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in mosts cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @family irradiance functions
#'
e_irrad <- function(spct, w.band, quantity, time.unit, wb.trim, use.cached.mult, use.hinges, allow.scaled) UseMethod("e_irrad")

#' @describeIn e_irrad Default for generic function
#'
#' @export
#'
e_irrad.default <- function(spct, w.band, quantity, time.unit, wb.trim, use.cached.mult, use.hinges, allow.scaled) {
  warning("'e_irrad' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn e_irrad  Calculates energy irradiance from a \code{source_spct}
#'   object.
#'
#' @export
#'
e_irrad.source_spct <-
  function(spct, w.band=NULL,
           quantity="total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL),
           allow.scaled = FALSE ) {
    irrad_spct(spct, w.band=w.band, unit.out="energy", quantity=quantity,
               time.unit = time.unit, wb.trim=wb.trim,
               use.cached.mult=use.cached.mult, use.hinges=use.hinges,
               allow.scaled = allow.scaled)
  }

# photon irradiance -------------------------------------------------------


#' Calculate photon irradiance from spectral irradiance.
#'
#' This function returns the photon irradiance (or quantum irradiance) for a
#' given waveband of a light source spectrum.
#'
#' @param spct an R object
#' @param w.band a list of \code{waveband} objects or a \code{waveband} object
#' @param quantity character string
#' @param time.unit character or lubridate::duration
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#' @param allow.scaled logical indicating whether scaled or normalized spectra
#'   as argument to spct are flagged as an error
#'
#' @keywords manip misc
#'
#' @export
#'
#' @examples
#' q_irrad(sun.spct, new_waveband(400,700))
#'
#' @return One numeric value for each waveband with no change in scale factor,
#'   with name attribute set to the name of each waveband unless a named list is
#'   supplied in which case the names of the list elements are used. The
#'   time.unit attribute is copied from the spectrum object to the output. Units
#'   are as follows: If time.unit is second, [W m-2 nm-1] -> [mol s-1 m-2] If
#'   time.unit is day, [J d-1 m-2 nm-1] -> [mol d-1 m-2]
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in mosts cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @export
#' @family irradiance functions
q_irrad <- function(spct, w.band, quantity, time.unit, wb.trim, use.cached.mult, use.hinges, allow.scaled) UseMethod("q_irrad")

#' @describeIn q_irrad Default for generic function
#'
#' @export
#'
q_irrad.default <- function(spct, w.band, quantity, time.unit, wb.trim, use.cached.mult, use.hinges, allow.scaled) {
  warning("'q_irrad' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn q_irrad  Calculates photon irradiance from a \code{source_spct}
#'   object.
#'
#' @export
#'
q_irrad.source_spct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL),
           allow.scaled = FALSE ) {
    irrad_spct(spct, w.band=w.band, unit.out="photon", quantity=quantity,
               time.unit = time.unit, wb.trim=wb.trim,
               use.cached.mult=use.cached.mult, use.hinges=use.hinges,
               allow.scaled = allow.scaled)
  }


# fluence -----------------------------------------------------------------

#' Calculate energy or photon fluence from spectral irradiance
#'
#' This function returns the energy or photon fluence for a given waveband of a
#' light source spectrum and the duration of the exposure.
#'
#' @param spct an R object
#' @param w.band a list of \code{waveband} objects or a \code{waveband} object
#' @param unit.out character string with allowed values "energy", and "photon",
#'   or its alias "quantum"
#' @param exposure.time lubridate::duration
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#' @param allow.scaled logical indicating whether scaled or normalized spectra
#'   as argument to spct are flagged as an error
#'
#' @keywords manip misc
#'
#' @export
#'
#' @examples
#' library(lubridate)
#' fluence(sun.spct,
#'         w.band = new_waveband(400,700),
#'         exposure.time = duration(3, "minutes") )
#'
#' @return One numeric value for each waveband with no change in scale factor,
#'   with name attribute set to the name of each waveband unless a named list is
#'   supplied in which case the names of the list elements are used. The
#'   time.unit attribute is copied from the spectrum object to the output. Units
#'   are as follows: If time.unit is second, [W m-2 nm-1] -> [mol s-1 m-2] If
#'   time.unit is day, [J d-1 m-2 nm-1] -> [mol d-1 m-2]
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in mosts cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @export
#' @family irradiance functions
fluence <- function(spct, w.band, unit.out, exposure.time, wb.trim,
                    use.cached.mult, use.hinges, allow.scaled) UseMethod("fluence")

#' @describeIn fluence Default for generic function
#'
#' @export
#'
fluence.default <- function(spct, w.band, unit.out, exposure.time,
                            wb.trim, use.cached.mult, use.hinges, allow.scaled) {
  warning("'fluence' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn fluence  Calculate photon fluence from a \code{source_spct}
#'   object and the duration of the exposure
#'
#' @export
#'
fluence.source_spct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           exposure.time = NA,
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL),
           allow.scaled = FALSE ) {
    if (!lubridate::is.duration(exposure.time)) {
      exposure.time <- lubridate::as.duration(exposure.time)
    }
    return.value <-
      irrad_spct(spct, w.band=w.band, unit.out=unit.out, quantity="total",
                 time.unit = exposure.time, wb.trim=wb.trim,
                 use.cached.mult=use.cached.mult, use.hinges=use.hinges,
                 allow.scaled = allow.scaled)
    if (unit.out %in% c("photon", "quantum")) {
      setattr(return.value, "radiation.unit", "photon fluence (mol m-2)")
    } else if (unit.out == "energy") {
      setattr(return.value, "radiation.unit", "energy fluence (J m-2)")
    }
    setattr(return.value, "exposure.duration", exposure.time)
    setattr(return.value, "time.unit", NULL)
    return(return.value)
  }


# photon fluence ----------------------------------------------------------

#' Calculate photon fluence from spectral irradiance.
#'
#' This function returns the photon irradiance (or quantum irradiance) for a
#' given waveband of a light source spectrum.
#'
#' @param spct an R object
#' @param w.band a list of \code{waveband} objects or a \code{waveband} object
#' @param exposure.time lubridate::duration
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#' @param allow.scaled logical indicating whether scaled or normalized spectra
#'   as argument to spct are flagged as an error
#'
#' @keywords manip misc
#'
#' @export
#'
#' @examples
#' library(lubridate)
#' q_fluence(sun.spct,
#'           w.band = new_waveband(400,700),
#'           exposure.time = duration(3, "minutes") )
#'
#' @return One numeric value for each waveband with no change in scale factor,
#'   with name attribute set to the name of each waveband unless a named list is
#'   supplied in which case the names of the list elements are used. The
#'   exposure.time is copied from the spectrum object to the output as an attibute.
#'   Units are as follows: moles of photons per exposure.
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in mosts cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @export
#' @family irradiance functions
q_fluence <- function(spct, w.band, exposure.time, wb.trim, use.cached.mult,
                      use.hinges, allow.scaled) UseMethod("q_fluence")

#' @describeIn q_fluence Default for generic function
#'
#' @export
#'
q_fluence.default <- function(spct, w.band, exposure.time, wb.trim,
                              use.cached.mult, use.hinges, allow.scaled) {
  warning("'q_fluence' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn q_fluence  Calculate photon fluence from a \code{source_spct}
#'   object and the duration of the exposure
#'
#' @export
#'
q_fluence.source_spct <-
  function(spct, w.band = NULL,
           exposure.time = NA,
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL),
           allow.scaled = FALSE ) {
    if (!lubridate::is.duration(exposure.time)) {
      exposure.time <- lubridate::as.duration(exposure.time)
    }

    return.value <-
      irrad_spct(spct, w.band=w.band, unit.out="photon", quantity="total",
                 time.unit = exposure.time, wb.trim=wb.trim,
                 use.cached.mult=use.cached.mult, use.hinges=use.hinges,
                 allow.scaled = allow.scaled)
    setattr(return.value, "radiation.unit", "photon fluence (mol m-2)")
    setattr(return.value, "exposure.duration", exposure.time)
    setattr(return.value, "time.unit", NULL)
    return(return.value)
  }


# energy fluence ----------------------------------------------------------

#' Calculate energy fluence from spectral irradiance.
#'
#' This function returns the energy flurnce for a given waveband of a light
#' source spectrum given the duration of the exposure.
#'
#' @param spct an R object
#' @param w.band a list of \code{waveband} objects or a \code{waveband} object
#' @param exposure.time lubridate::duration
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#' @param allow.scaled logical indicating whether scaled or normalized spectra
#'   as argument to spct are flagged as an error
#'
#' @keywords manip misc
#'
#' @export
#'
#' @examples
#' library(lubridate)
#' e_fluence(sun.spct, w.band = new_waveband(400,700),
#'           exposure.time = duration(3, "minutes") )
#'
#' @return One numeric value for each waveband with no change in scale factor,
#'   with name attribute set to the name of each waveband unless a named list is
#'   supplied in which case the names of the list elements are used. The
#'   exposure.time is copied to the output as an attribute. Units are as
#'   follows: (J) joules per exposure.
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in mosts cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @export
#' @family irradiance functions
e_fluence <- function(spct, w.band, exposure.time, wb.trim, use.cached.mult,
                      use.hinges, allow.scaled) UseMethod("e_fluence")

#' @describeIn e_fluence Default for generic function
#'
#' @export
#'
e_fluence.default <- function(spct, w.band, exposure.time, wb.trim, use.cached.mult,
                              use.hinges, allow.scaled) {
  warning("'e_fluence' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn e_fluence  Calculate energy fluence from a \code{source_spct}
#'   object and the duration of the exposure.
#'
#' @export
#'
e_fluence.source_spct <-
  function(spct, w.band = NULL,
           exposure.time = NA,
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL),
           allow.scaled = FALSE ) {
    if (!lubridate::is.duration(exposure.time)) {
      exposure.time <- lubridate::as.duration(exposure.time)
    }
  return.value <-
      irrad_spct(spct, w.band=w.band, unit.out="energy", quantity="total",
                 time.unit = exposure.time, wb.trim=wb.trim,
                 use.cached.mult=use.cached.mult, use.hinges=use.hinges,
                 allow.scaled = allow.scaled)
  setattr(return.value, "radiation.unit", "energy fluence (J m-2)")
  setattr(return.value, "exposure.duration", exposure.time)
  setattr(return.value, "time.unit", NULL)
  return(return.value)
  }


