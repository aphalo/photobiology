
# irradiance --------------------------------------------------------------

#' Calculate irradiance from spectral irradiance.
#'
#' This function returns the irradiance for a given
#' waveband of a light source spectrum.
#'
#' @usage irrad_spct(spct, w.band=NULL,
#'                   unit.out=getOption("photobiology.radiation.unit", default="energy"),
#'                   quantity="total", wb.trim=NULL, use.cached.mult=FALSE, use.hinges=NULL, allow.scaled = FALSE)
#'
#' @param spct an object of class "source.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @param allow.scaled logical indicating whether rescaled or normalized spectra as argument to spct are flagged as an error
#'
#' @note Formal parameter allow scaled is mainly used internally for calculation of ratios, as rescaling
#' and normalization do not invalidate the calcualtion of ratios.
#'
#' @return One numeric value for each waveband with no change in scale factor, with name attribute set to
#' the name of each waveband unless a named list is supplied in which case the names of the list elements are
#' used. The time.unit attribute is copied from the spectrum object to the output. Units are as follows:
#' If time.unit is second, [W m-2 nm-1] -> [mol s-1 m-2] or [W m-2 nm-1] -> [W m-2]
#' If time.unit is day, [J d-1 m-2 nm-1] -> [mol d-1 m-2] or [J d-1 m-2 nm-1] -> [J m-2]
#'
#' @keywords manip misc
#' @export
#' @examples
#' irrad_spct(sun.spct, new_waveband(400,700), "photon")
#' irrad_spct(sun.spct, new_waveband(400,700), "energy")
#'
#' @note The last two parameters control speed optimizations. The defaults should be suitable
#' in mosts cases. If you will use repeatedly
#' the same SWFs on many spectra measured at exactly the same wavelengths you may obtain some speed up
#' by setting \code{use.cached.mult=TRUE}. However, be aware that you are responsible for ensuring
#' that the wavelengths are the same in each call, as the only test done is for the length of the
#' \code{w.length} vector.
#'
#' @aliases irrad.source.spct irrad_spct

irrad_spct <-
  function(spct, w.band=NULL, unit.out=getOption("photobiology.radiation.unit", default="energy"),
           quantity="total", wb.trim=NULL, use.cached.mult=FALSE, use.hinges=NULL, allow.scaled = FALSE){
    # we have a default, but we check for invalid arguments
    if (!allow.scaled && (is.normalized(spct) || is.rescaled(spct))) {
      warning("The espectral data has been normalized or rescaled, making impossible to calculate irradiance")
      return(NA)
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
      length.wl <- length(spct_x$w.length)
      use.hinges <- (spct_x$w.length[length.wl] - spct_x$w.length[1]) / length.wl > 1.1
      # we use 1.1 nm as performance degradation by using hinges is very significant
      # in the current version.
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

    # "source.spct" objects are not guaranteed to contain spectral irradiance
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
    for (wb in w.band) {
      i <- i + 1L
      # get names from wb if needed
      if (wb.name[i] == "") {
        wb.name[i] <- wb$name
      }
      # calculate the multipliers
      mult <- calc_multipliers(w.length=spct_x$w.length, w.band=wb, unit.out=unit.out,
                               unit.in=unit.in, use.cached.mult=use.cached.mult)
      # calculate weighted spectral irradiance
      if (unit.out == "energy") {
        irr <- with(spct_x, integrate_irradiance(w.length, s.e.irrad * mult))
      } else {
        irr <- with(spct_x, integrate_irradiance(w.length, s.q.irrad * mult))
      }
      irrad[i] <- irr
    }
    if (quantity %in% c("contribution", "contribution.pc")) {
      if (any(sapply(w.band, is_effective))) {
        warning("'quantity '", quantity, "' not supported when using BSWFs, returning 'total' instead")
        quantity <- "total"
      } else {
        total <- irrad_spct(spct_x, w.band=NULL, unit.out=unit.out,
                            quantity="total", use.cached.mult=FALSE, use.hinges=FALSE)
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
    setattr(irrad, "radiation.unit", paste(unit.out, "irradiance", quantity))
    return(irrad)
  }

#' @method irrad source.spct
#' @export
irrad.source.spct <- irrad_spct

#' Default for generic function
#'
#' Calculate energy or photon irradiance.
#'
#' @param spct an object of class "generic.spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @export
#'
irrad.default <- function(spct, w.band, unit.out, quantity, wb.trim, use.cached.mult, use.hinges) {
  return(NA)
}

#' Generic function
#'
#' Calculate energy or photon irradiance.
#'
#' @param spct an R object of class "generic.spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @export
#'
irrad <- function(spct, w.band, unit.out, quantity, wb.trim, use.cached.mult, use.hinges) UseMethod("irrad")

# energy irradiance -------------------------------------------------------


#' Calculate energy irradiance from spectral irradiance.
#'
#' This function returns the energy irradiance for a given
#' waveband of a light source spectrum.
#'
#' @usage e_irrad.source.spct(spct, w.band=NULL,
#'                quantity="total", wb.trim=NULL, use.cached.mult=FALSE, use.hinges=NULL)
#'
#' @param spct an object of class "source.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE,
#'        they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @keywords manip misc
#'
#' @export
#'
#' @examples
#' e_irrad(sun.spct, new_waveband(400,700))
#'
#' @return One numeric value for each waveband with no change in scale factor, with name attribute set to
#' the name of each waveband unless a named list is supplied in which case the names of the list elements are
#' used. The time.unit attribute is copied from the spectrum object to the output. Units are as follows:
#' If time.unit is second, [W m-2 nm-1] -> [W m-2]
#' If time.unit is day, [J d-1 m-2 nm-1] -> [J m-2]
#'
#' @note The last two parameters control speed optimizations. The defaults should be suitable
#' in mosts cases. If you will use repeatedly
#' the same SWFs on many spectra measured at exactly the same wavelengths you may obtain some speed up
#' by setting \code{use.cached.mult=TRUE}. However, be aware that you are responsible for ensuring
#' that the wavelengths are the same in each call, as the only test done is for the length of the
#' \code{w.length} vector.
#'
e_irrad.source.spct <-
  function(spct, w.band=NULL,
           quantity="total", wb.trim=NULL, use.cached.mult=FALSE, use.hinges=NULL){
    irrad_spct(spct, w.band=w.band, unit.out="energy", quantity=quantity, wb.trim=wb.trim,
                      use.cached.mult=use.cached.mult, use.hinges=use.hinges)
  }

# photon irradiance -------------------------------------------------------


#' Calculate photon irradiance from spectral irradiance.
#'
#' This function returns the photon irradiance for a given
#' waveband of a light source spectrum.
#'
#' @usage q_irrad.source.spct(spct, w.band=NULL,
#'                            quantity="total", wb.trim=NULL,
#'                            use.cached.mult=FALSE, use.hinges=NULL)
#'
#' @param spct an object of class "source.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @keywords manip misc
#'
#' @export
#'
#' @note The three functions are at the moment identical. The "_spct" versions are deprecated
#'
#' @examples
#' q_irrad(sun.spct, new_waveband(400,700))
#'
#' @return One numeric value for each waveband with no change in scale factor, with name attribute set to
#' the name of each waveband unless a named list is supplied in which case the names of the list elements are
#' used. The time.unit attribute is copied from the spectrum object to the output. Units are as follows:
#' If time.unit is second, [W m-2 nm-1] -> [mol s-1 m-2]
#' If time.unit is day, [J d-1 m-2 nm-1] -> [mol d-1 m-2]
#'
#' @note The last two parameters control speed optimizations. The defaults should be suitable
#' in mosts cases. If you will use repeatedly
#' the same SWFs on many spectra measured at exactly the same wavelengths you may obtain some speed up
#' by setting \code{use.cached.mult=TRUE}. However, be aware that you are responsible for ensuring
#' that the wavelengths are the same in each call, as the only test done is for the length of the
#' \code{w.length} vector.
#'
#' @export

q_irrad.source.spct <-
  function(spct, w.band=NULL,
           quantity="total",
           wb.trim=NULL,
           use.cached.mult=FALSE,
           use.hinges=NULL){
    irrad_spct(spct, w.band=w.band, unit.out="photon", quantity=quantity, wb.trim=wb.trim,
                      use.cached.mult=use.cached.mult, use.hinges=use.hinges)
  }

#' Generic function
#'
#' Calculate (energy) irradiance.
#'
#' @param spct an R object of class "generic.spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @export
#'
e_irrad <- function(spct, w.band, quantity, wb.trim, use.cached.mult, use.hinges) UseMethod("e_irrad")

#' Generic function
#'
#' Calculate photon irradiance.
#'
#' @param spct an R object of class "generic.spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @export
#'
q_irrad <- function(spct, w.band, quantity, wb.trim, use.cached.mult, use.hinges) UseMethod("q_irrad")

#' Default for generic function
#'
#' Calculate (energy) irradiance.
#'
#' @param spct an object of class "generic.spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @export
#'
e_irrad.default <- function(spct, w.band, quantity, wb.trim, use.cached.mult, use.hinges) {
  return(NA)
}

#' Default for generic function
#'
#' Calculate photon irradiance.
#'
#' @param spct an object of class "generic.spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @export
#'
q_irrad.default <- function(spct, w.band, quantity, wb.trim, use.cached.mult, use.hinges) {
  return(NA)
}
