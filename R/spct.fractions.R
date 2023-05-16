#' Photon:photon fraction
#'
#' This function returns the photon fraction for a given pair of wavebands of a
#' light source spectrum.
#'
#' @param spct an object of class "source_spct".
#' @param w.band.num waveband object or a list of waveband objects used to
#'   compute the numerator(s) and denominator(s) of the fraction(s).
#' @param w.band.denom waveband object or a list of waveband objects used to
#'   compute the denominator(s) of the fraction(s).
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param name.tag character Used to tag the name of the returned values.
#' @param ... other arguments (possibly ignored)
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector of adimensional values giving a photon fraction between integrated
#'   photon irradiances for pairs of wavebands, with name attribute set to
#'   the name of the wavebands unless a named list of wavebands is supplied in
#'   which case the names of the list elements are used, with "(q:q)" appended.
#'   A \code{data.frame} in the case of collections of spectra, containing one
#'   column for each fraction definition, an index column with the names of the
#'   spectra, and optionally additional columns with metadata values retrieved
#'   from the attributes of the member spectra.
#'
#'   Fraction definitions are "assembled" from the arguments passed to
#'   \code{w.band.num} and \code{w.band.denom}. If both arguments are of equal
#'   length, then the wavebands are paired to obtain as many fractions as the
#'   number of wavebands in each list. Recycling for wavebands takes place when
#'   the number of denominator and numerator wavebands differ.
#'
#' @export
#' @examples
#' q_fraction(sun.spct, new_waveband(400,500), new_waveband(400,700))
#'
#' @note The last two parameters control speed
#'   optimizations. The defaults should be suitable in most cases. If you will
#'   use repeatedly the same SWFs on many spectra measured at exactly the same
#'   wavelengths you may obtain some speed up by setting
#'   \code{use.cached.mult=TRUE}. However, be aware that you are responsible for
#'   ensuring that the wavelengths are the same in each call, as the only test
#'   done is for the length of the \code{w.length} vector.
#'
#' @family photon and energy ratio functions
#'
q_fraction <- function(spct, w.band.num, w.band.denom, scale.factor, wb.trim,
                    use.cached.mult, use.hinges, ...) UseMethod("q_fraction")

#' @describeIn q_fraction Default for generic function
#'
#' @export
#'
q_fraction.default <- function(spct, w.band.num, w.band.denom, scale.factor, wb.trim,
                            use.cached.mult, use.hinges, ...) {
  warning("'q_fraction' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn q_fraction Method for \code{source_spct} objects
#'
#' @export
#'
q_fraction.source_spct <-
  function(spct,
           w.band.num = NULL,
           w.band.denom = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = NULL,
           naming = "short",
           name.tag = ifelse(naming != "none", "[q:q]", ""),
           ... ) {

    irrads <- two_irrads(spct = spct,
                         w.band.num = w.band.num,
                         w.band.denom = w.band.denom,
                         unit.out.num = "photon",
                         unit.out.denom = "photon",
                         quantity = "total",
                         wb.trim = wb.trim,
                         use.cached.mult = use.cached.mult,
                         use.hinges = use.hinges,
                         naming = naming,
                         ...)

    q.irrad.num <- irrads[["irrad.num"]]
    q.irrad.denom <- irrads[["irrad.denom"]]
    fraction <- q.irrad.num / (q.irrad.denom + q.irrad.num) * scale.factor
    names(fraction) <- paste(names(q.irrad.num), ":", names(q.irrad.denom), name.tag, sep = "")
    attr(fraction, "time.unit") <- NULL
    attr(fraction, "radiation.unit") <- "q:q fraction"
    return(fraction)
  }

#' Energy:energy fraction
#'
#' This function returns the energy fraction for a given pair of wavebands of a
#' light source spectrum.
#'
#' @param spct source_spct
#' @param w.band.num waveband object or a list of waveband objects used to
#'   compute the numerator(s) and denominator(s) of the fraction(s).
#' @param w.band.denom waveband object or a list of waveband objects used to
#'   compute the denominator(s) of the fraction(s).
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical Flag telling whether multiplier values should be
#'   cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param name.tag character Used to tag the name of the returned values.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector of adimensional values giving a energy fraction between
#'   integrated energy irradiances for pairs of wavebands, with name attribute
#'   set to the name of the wavebands unless a named list of wavebands is
#'   supplied in which case the names of the list elements are used, with
#'   "(e:e)" appended. A \code{data.frame} in the case of collections of
#'   spectra, containing one column for each fraction definition, an index column
#'   with the names of the spectra, and optionally additional columns with
#'   metadata values retrieved from the attributes of the member spectra.
#'
#'   Fraction definitions are "assembled" from the arguments passed to
#'   \code{w.band.num} and \code{w.band.denom}. If both arguments are of equal
#'   length, then the wavebands are paired to obtain as many fractions as the
#'   number of wavebands in each list. Recycling for wavebands takes place when
#'   the number of denominator and numerator wavebands differ.
#'
#' @export
#' @examples
#' e_fraction(sun.spct, new_waveband(400,700), new_waveband(400,500))
#'
#' @note Recycling for wavebands takes place when the number of denominator and
#'   denominator wavebands differ. The last two parameters control speed
#'   optimizations. The defaults should be suitable in most cases. If you will
#'   use repeatedly the same SWFs on many spectra measured at exactly the same
#'   wavelengths you may obtain some speed up by setting
#'   \code{use.cached.mult=TRUE}. However, be aware that you are responsible for
#'   ensuring that the wavelengths are the same in each call, as the only test
#'   done is for the length of the \code{w.length} vector.
#'
#' @family photon and energy fraction functions
#'
e_fraction <- function(spct, w.band.num, w.band.denom, scale.factor, wb.trim,
                    use.cached.mult, use.hinges, ...) UseMethod("e_fraction")

#' @describeIn e_fraction Default for generic function
#'
#' @export
#'
e_fraction.default <- function(spct, w.band.num, w.band.denom, scale.factor, wb.trim,
                            use.cached.mult, use.hinges, ...) {
  warning("'e_fraction' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn e_fraction Method for \code{source_spct} objects
#'
#' @export
#'
e_fraction.source_spct <-
  function(spct,
           w.band.num = NULL, w.band.denom = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = NULL,
           naming = "short",
           name.tag = ifelse(naming != "none", "[e:e]", ""),
           ...) {

    irrads <- two_irrads(spct = spct,
                         w.band.num = w.band.num,
                         w.band.denom = w.band.denom,
                         unit.out.num = "energy",
                         unit.out.denom = "energy",
                         quantity = "total",
                         wb.trim = wb.trim,
                         use.cached.mult = use.cached.mult,
                         use.hinges = use.hinges,
                         naming = naming,
                         ...)

    e.irrad.num <- irrads[["irrad.num"]]
    e.irrad.denom <- irrads[["irrad.denom"]]
    fraction <- e.irrad.num / (e.irrad.denom + e.irrad.num) * scale.factor
    names(fraction) <- paste(names(e.irrad.num), ":", names(e.irrad.denom), name.tag, sep="")
    attr(fraction, "time.unit") <- NULL
    attr(fraction, "radiation.unit") <- "e:e fraction"
    return(fraction)
  }

