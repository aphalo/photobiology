#' Photon or energy irradiance from spectral energy or photon irradiance.
#'
#' Energy or photon irradiance for one or more wavebands of a radiation
#' spectrum.
#'
#' @param w.length numeric Vector of wavelength (nm).
#' @param s.irrad numeric vector of spectral (energy) irradiances (W m-2 nm-1).
#' @param w.band waveband or list of waveband objects The waveband(s) determine
#'   the region(s) of the spectrum that are summarized.
#' @param unit.out character Allowed values "energy", and "photon", or its alias
#'   "quantum".
#' @param unit.in character Allowed values "energy", and "photon", or its alias
#'   "quantum".
#' @param check.spectrum logical Flag indicating whether to sanity check input
#'   data, default is TRUE.
#' @param use.cached.mult logical Flag indicating whether multiplier values
#'   should be cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#'
#' @return A single numeric value or a vector of numeric values with no change
#'   in scale factor: [W m-2 nm-1] -> [mol s-1 m-2]
#'
#' @export
#' @examples
#'
#' with(sun.data, irradiance(w.length, s.e.irrad, new_waveband(400,700), "photon"))
#' @note The last three parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you set \code{check.spectrum=FALSE}
#'   then you should call \code{check_spectrum()} at least once for your
#'   spectrum before using any of the other functions. If you will use
#'   repeatedly the same SWFs on many spectra measured at exactly the same
#'   wavelengths you may obtain some speed up by setting
#'   \code{use.cached.mult=TRUE}. However, be aware that you are responsible for
#'   ensuring that the wavelengths are the same in each call, as the only test
#'   done is for the length of the \code{w.length} vector. The is no reason for
#'   setting \code{use.cpp.code=FALSE} other than for testing the improvement in
#'   speed, or in cases where there is no suitable C++ compiler for building the
#'   package.
#'
#' @family low-level functions operating on numeric vectors.
#'
irradiance <-
  function(w.length, s.irrad, w.band = NULL, unit.out = NULL, unit.in = "energy",
           check.spectrum = TRUE,
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges", default = NULL) ){
    # what output? seems safer to not have a default here
    if (is.null(unit.out)){
      warning("'unit.out' has no default value")
      return(NA_real_)
    }
    # make code a bit simpler further down
    if (unit.in=="quantum") {unit.in <- "photon"}
    # sanity check for spectral data and abort if check fails
    if (check.spectrum && !check_spectrum(w.length, s.irrad)) {
      return(NA_real_)
    }
    # if the waveband is undefined then use all data
    if (length(w.band) == 0){
#      w.band <- new_waveband(min(w.length), max(w.length))
      w.band <- new_waveband(min(w.length), max(w.length) + 1e-12)
      # we need to add a small number as the test is "<"
      # this affects significantly the result only when no hinges are used
    }
    if (is.waveband(w.band)) {
      # if the argument is a single w.band, we enclose it in a list
      # so that the for loop works as expected. This is a bit of a
      # kludge but let's us avoid treating it as a special case
      w.band <- list(w.band)
    }
    # if the w.band includes 'hinges' we insert them
    # depending on the value of use.hinges. If no explicit argument
    # was supplied, we decide based on the step size for w.length
    if (is.null(use.hinges)) {
      use.hinges <- auto_hinges(w.length)
    }
    # we collect all hinges and insert them in one go.
    # This may alter very slightly the returned values
    # for overlapping wavebands but should be faster.
    if (use.hinges) {
      all.hinges <- NULL
      for (wb in w.band) {
        if (!is.null(wb[["hinges"]]) & length(wb[["hinges"]]) > 0) {
          all.hinges <- c(all.hinges, wb[["hinges"]])
        }
      }
      if (!is.null(all.hinges)) {
        new.data <- l_insert_hinges(x = w.length, y = s.irrad, all.hinges)
        w.length <- new.data[["x"]]
        s.irrad <- new.data[["y"]]
      }
    }
    wb_name <- names(w.band)
    no_names_flag <- is.null(wb_name)
    if (no_names_flag) wb_name <- character(length(w.band))
    irrad <- numeric(length(w.band))

    i <- 0
    for (wb in w.band) {
      i <- i + 1
      # get names from wb if needed
      if (no_names_flag) wb_name[i] <- wb[["name"]]
      # calculate the multipliers
      mult <- calc_multipliers(w.length = w.length, w.band = wb,
                               unit.out = unit.out, unit.in = unit.in,
                               use.cached.mult = use.cached.mult)
      # calculate weighted spectral irradiance
      irr <- integrate_xy(w.length, s.irrad * mult)
      irrad[i] <- irr
    }
    names(irrad) <- wb_name
    return(irrad)
  }
