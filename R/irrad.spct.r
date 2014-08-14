
# irradiance --------------------------------------------------------------

#' Calculate irradiance from spectral irradiance.
#'
#' This function returns the irradiance for a given
#' waveband of a light source spectrum.
#'
#' @usage irrad(spct, w.band=NULL, unit.out=NULL,
#' use.cached.mult=FALSE, use.hinges=NULL)
#'
#' @usage irrad_spct(spct, w.band=NULL, unit.out=NULL,
#' use.cached.mult=FALSE, use.hinges=NULL)
#'
#' @param spct an object of class "source.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @return One numeric value for each waveband with no change in scale factor, with name attribute set to
#' the name of each waveband unless a named list is supplied in which case the names of the list elements are
#' used. The time.unit attribute is copied from the spectrum object to the output. Units are as follows:
#' If time.unit is second, [W m-2 nm-1] -> [mol s-1 m-2] or [W m-2 nm-1] -> [W m-2]
#' If time.unit is day, [J d-1 m-2 nm-1] -> [mol d-1 m-2] or [J d-1 m-2 nm-1] -> [J m-2]
#'
#' @keywords manip misc
#' @export irrad irrad_spct
#' @examples
#' data(sun.spct)
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
#' @aliases irrad irrad_spct

irrad <-
  function(spct, w.band=NULL, unit.out=NULL,
           use.cached.mult=FALSE, use.hinges=NULL){
    # what output? seems safer to not have a default here
    if (is.null(unit.out) || is.na(unit.out)){
      warning("'unit.out' has no default value")
      return(NA)
    }
    if (unit.out == "quantum") {
      unit.out <- "photon"
    }
    if (is.null(w.band)) {
      w.band <- new_waveband(min(spct), max(spct) + 1e-4)
    }
    if (is(w.band, "waveband")) {
      # if the argument is a single w.band, we enclose it in a list
      # so that the for loop works as expected.This is a bit of a
      # cludge but lets us avoid treating it as a special case
      w.band <- list(w.band)
    }
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
      length.wl <- length(spct$w.length)
      use.hinges <- (spct$w.length[length.wl] - spct$w.length[1]) / length.wl > 1.1
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
        spct <- insert_spct_hinges(spct, all.hinges)
      }
    }

    # "source.spct" objects are not guaranteed to contain spectral irradiance
    # expressed in the needed type of scale, if the needed one is missing
    # we add the missing it.
    # As spectra are passed by reference the changes propagate to the argument
    if (unit.out == "energy") {
      q2e(spct, byref=TRUE)
    } else if (unit.out == "photon") {
      e2q(spct, byref=TRUE)
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
      mult <- calc_multipliers(w.length=spct$w.length, w.band=wb, unit.out=unit.out,
                               unit.in=unit.in, use.cached.mult=use.cached.mult)
      # calculate weighted spectral irradiance
      if (unit.out == "energy") {
        irr <- with(spct, integrate_irradiance(w.length, s.e.irrad * mult))
      } else {
        irr <- with(spct, integrate_irradiance(w.length, s.q.irrad * mult))
      }
      irrad[i] <- irr
    }

    names(irrad) <- wb.name
    attr(irrad, "time.unit") <- attr(spct, "time.unit", exact=TRUE)
    return(irrad)
  }

irrad_spct <- irrad


# energy irradiance -------------------------------------------------------


#' Calculate energy irradiance from spectral irradiance.
#'
#' This function returns the energy irradiance for a given
#' waveband of a light source spectrum.
#'
#' @usage e_irrad(spct, w.band=NULL,
#' use.cached.mult=FALSE, use.hinges=NULL)
#'
#' @usage e_irrad_spct(spct, w.band=NULL,
#' use.cached.mult=FALSE, use.hinges=NULL)
#'
#' @param spct an object of class "source.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @return a single numeric value with no change in scale factor: [W m-2 nm-1] -> [W m-2]
#' @keywords manip misc
#' @export e_irrad e_irrad_spct
#' @examples
#' data(sun.spct)
#' e_irrad_spct(sun.spct, new_waveband(400,700))
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
#' @aliases e_irrad e_irrad_spect

e_irrad <-
  function(spct, w.band=NULL, use.cached.mult=FALSE, use.hinges=NULL){
    return(irrad(spct, w.band=w.band, unit.out="energy",
                      use.cached.mult=use.cached.mult, use.hinges=use.hinges))
  }

e_irrad_spct <- e_irrad


# photon irradiance -------------------------------------------------------


#' Calculate photon irradiance from spectral irradiance.
#'
#' This function returns the photon irradiance for a given
#' waveband of a light source spectrum.
#'
#' @usage q_irrad(spct, w.band=NULL,
#' use.cached.mult=FALSE, use.hinges=NULL)
#'
#' @usage q_irrad_spct(spct, w.band=NULL,
#' use.cached.mult=FALSE, use.hinges=NULL)
#'
#' @param spct an object of class "source.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @return a single numeric value with no change in scale factor: [W m-2 nm-1] -> [mol s-1 m-2]
#' @keywords manip misc
#' @export q_irrad q_irrad_spct
#'
#' @note The three functions are at the moment identical. The "_spct" versions are deprecated
#' @examples
#' data(sun.spct)
#' q_irrad_spct(sun.spct, new_waveband(400,700))
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
#' @name q_irrad
#' @aliases q_irrad q_irrad_spct

q_irrad <-
  function(spct, w.band=NULL, use.cached.mult=FALSE, use.hinges=NULL){
    return(irrad(spct, w.band=w.band, unit.out="photon",
                      use.cached.mult=use.cached.mult, use.hinges=use.hinges))
  }

q_irrad_spct <- q_irrad

