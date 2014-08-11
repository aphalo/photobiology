#' Calculate irradiance from spectral irradiance.
#'
#' This function returns the irradiance for a given
#' waveband of a light source spectrum.
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
#' @return a single numeric value with no change in scale factor: [W m-2 nm-1] -> [mol s-1 m-2]
#' @keywords manip misc
#' @export
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

irrad_spct <-
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
    return(irrad)
  }
