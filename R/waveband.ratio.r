#' Photon or energy ratio
#'
#' This function gives the (energy or photon) irradiance ratio between two given
#' wavebands of a radiation spectrum.
#'
#' @param w.length numeric Vector of wavelengths [\eqn{nm}].
#' @param s.irrad numeric vector of spectral irradiances in
#'   [\eqn{W\,m^{-2}\,nm^{-1}}{W m-2 nm-1}] or
#'   [\eqn{mol\,s^{-1}\,sm^{-2}\,nm^{-1}}{mol s-1 m-2 nm-1}] as indicated by the
#'   argument pased to \code{unit.in}.
#' @param w.band.num,w.band.denom waveband objects used to compute the numerator
#'   and denominator of the ratio.
#' @param unit.out.num,unit.out.denom character Base of expression used to
#'   compute the numerator and denominator of the ratio. Allowed values
#'   \code{"energy"}, and \code{"photon"}, or its alias \code{"quantum"}.
#' @param unit.in character Allowed values \code{"energy"}, and \code{"photon"},
#'   or its alias \code{"quantum"}.
#' @param check.spectrum logical Flag indicating whether to sanity check input
#'   data, default is TRUE.
#' @param use.cached.mult logical Flag indicating whether multiplier values
#'   should be cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#'
#' @note The default for both \code{w.band} parameters is a waveband covering
#'   the whole range of \code{w.length}. From version 0.9.19 onwards use of this
#'   default does not trigger a warning, but instead is used silently.
#'
#' @return a single numeric value giving the ratio
#'
#' @export
#' @examples
#' # photon:photon ratio
#' with(sun.data,
#'      waveband_ratio(w.length, s.e.irrad,
#'                     new_waveband(400,500),
#'                     new_waveband(400,700), "photon"))
#' # energy:energy ratio
#' with(sun.data,
#'      waveband_ratio(w.length, s.e.irrad,
#'                     new_waveband(400,500),
#'                     new_waveband(400,700), "energy"))
#' # energy:photon ratio
#' with(sun.data,
#'      waveband_ratio(w.length, s.e.irrad,
#'                     new_waveband(400,700),
#'                     new_waveband(400,700),
#'                     "energy", "photon"))
#' # photon:photon ratio waveband : whole spectrum
#' with(sun.data,
#'      waveband_ratio(w.length, s.e.irrad,
#'                     new_waveband(400,500),
#'                     unit.out.num="photon"))
#' # photon:photon ratio of whole spectrum should be equal to 1.0
#' with(sun.data,
#'      waveband_ratio(w.length, s.e.irrad,
#'      unit.out.num="photon"))
#'
waveband_ratio <-
  function(w.length, s.irrad,
           w.band.num = NULL, w.band.denom = NULL,
           unit.out.num = NULL, unit.out.denom = unit.out.num,
           unit.in = "energy",
           check.spectrum = TRUE,
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges",
                                  default = NULL)) {
    # We duplicate code from irradiance() here to avoid repeated checks
    # and calculations on the same data
    #
    # what output? seems safer to not have a default here
    if (is.null(unit.out.num) || is.null(unit.out.denom)) {
      warning("'unit.out.num' has no default value")
      return(NA_real_)
    }
    # make code a bit simpler further down
    if (unit.in == "quantum") {unit.in <- "photon"}
    # sanity check for wavelengths
    if (check.spectrum && !check_spectrum(w.length, s.irrad)) {
      return(NA_real_)
    }
    # if the waveband for numerator is undefined then use
    # the whole wavelength range of the spectrum for numerator
    if (is.null(w.band.num)) {
      w.band.num <- new_waveband(min(w.length),max(w.length))
    }
    # if the waveband for denominator is undefined then use
    # the whole wavelength range of the spectrum for denominator
    if (is.null(w.band.denom)) {
      w.band.denom <- new_waveband(min(w.length),max(w.length))
    }
    # choose whether to use hinges or not
    # if the user has specified its value, we leave it alone
    # but if it was not requested, we decide whether to use
    # it or not based of the wavelength resolution of the
    # spectrum.
    if (is.null(use.hinges)) {
      use.hinges <- auto_hinges(w.length)
    }
    # if the w.band.num and/or w.band.denom include 'hinges' we insert them.
    # it is o.k. to have hinges unsorted!
    # in new_waveband() NULL hinges are replaced with numeric(0)
    if (use.hinges) {
      merged.hinges <- c(w.band.denom[["hinges"]], w.band.num[["hinges"]])
      if (length(merged.hinges) > 0) {
        new.data <- l_insert_hinges(x = w.length, y = s.irrad, merged.hinges)
        w.length <- new.data[["x"]]
        s.irrad <- new.data[["y"]]
      }
    }
    # calculate the multipliers
    mult.num <- calc_multipliers(w.length, w.band.num,
                                 unit.out.num, unit.in,
                                 use.cached.mult = use.cached.mult)
    mult.denom <- calc_multipliers(w.length, w.band.denom,
                                   unit.out.denom, unit.in,
                                   use.cached.mult = use.cached.mult)

    # calculate weighted spectral irradiance
    irrad.num <- integrate_xy(w.length, s.irrad * mult.num)
    irrad.denom <- integrate_xy(w.length, s.irrad * mult.denom)

    return(irrad.num / irrad.denom)
  }
