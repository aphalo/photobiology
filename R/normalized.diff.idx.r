#' Calculate a normalized index.
#'
#' This method returns a normalized difference index value for an arbitrary
#' pair of wavebands. There are many such indexes in use, such as NDVI
#' (normalized difference vegetation index), NDWI (normalized difference water
#' index), NDMI (normalized difference moisture index), etc., the only
#' difference among then is in the wavebands used.
#'
#' @param spct an R object
#' @param plus.w.band,minus.w.band waveband objects The wavebands determine the
#'   regions of the spectrum used in the calculations.
#' @param f function used for integration taking spct as first argument and a
#'   list of wavebands as second argument.
#' @param ... additional arguments passed to f
#'
#' @return A named numeric value for the index, or a tibble depending on whether
#'   a spectrum or a collection of spectra is passed as first argument. If
#'   the wavelength range of \code{spct} does not fully overlap with both
#'   wavebands \code{NA} is silently returned.
#'
#' @export
#'
#' @details \code{f} is most frequently \code{\link{reflectance}}, but also
#'   \code{\link{transmittance}}, or even \code{\link{absorbance}},
#'   \code{\link{response}}, \code{\link{irradiance}} or a user-defined function
#'   can be used if there is a good reason for it. In every case \code{spct}
#'   should be of the class expected by \code{f}. When using two wavebands of
#'   different widths do consider passing to \code{f} a suitable \code{quantity}
#'   argument, for example to compare averages rather than integrals. Wavebands
#'   can describe weighting functions if desired.
#'
#'   \deqn{\mathrm{NDxI} = \frac{f(s, wb_\mathrm{plus}) - f(s, wb_\mathrm{minus})}{f(s, wb_\mathrm{plus}) + f(s, wb_\mathrm{minus})}}
#'
#' @note Some NDxI indexes are directly based on satellite instrument data, such
#'   as those in the Landsat satellites. To simulate such indexes using spectral
#'   reflectande as input, constructors of \code{waveband} definitions from package
#'   'photobiologyWavebands' can be useful.
#'
#' @seealso \code{\link{Rfr_normdiff}}
#'
#' @export
#'
normalized_diff_ind <-
  function(spct, plus.w.band, minus.w.band, f, ...) {
    UseMethod("normalized_diff_ind")
  }

#' @rdname normalized_diff_ind
#'
#' @note \code{normalised_diff_ind()} is a synonym for \code{normalized_diff_ind()}.
#'
#' @export
#'
normalised_diff_ind <- normalized_diff_ind

#' @rdname normalized_diff_ind
#'
#' @note \code{NDxI()} is a shorthand for \code{normalized_diff_ind()}.
#'
#' @export
#'
NDxI <- normalized_diff_ind

#' @describeIn normalized_diff_ind default
#'
#' @export
#'
normalized_diff_ind.default <-
  function(spct, plus.w.band, minus.w.band, f, ...) {
    warning("'normalized_diff_ind' is not defined for objects of class ",
            class(spct)[1])
    return(spct)
  }

#' @describeIn normalized_diff_ind
#'
#' @export
#'
normalized_diff_ind.generic_spct <- function(spct, plus.w.band, minus.w.band, f, ...) {
  # check that spectral data fully covers both wavebands
  min.wl.bands <- min(wl_min(plus.w.band), wl_min(minus.w.band))
  max.wl.bands <- max(wl_max(plus.w.band), wl_max(minus.w.band))
  if (wl_min(spct) > min.wl.bands || wl_max(spct) < max.wl.bands) {
    NA_real_
  } else {
    x <- as.numeric(f(spct, list(plus.w.band, minus.w.band), ...))
    z <- (x[1] - x[2]) / (x[1] + x[2])
    name <- paste("NDI ", as.character(substitute(f)), " [",
                  sub("range.", "", labels(plus.w.band)[["label"]]), "] - [",
                  sub("range.", "", labels(minus.w.band)[["label"]]), "]",
                  sep = "")
    names(z) <- name
    z
  }
}

#' @describeIn normalized_diff_ind
#'
#' @export
#'
normalized_diff_ind.generic_mspct <- function(spct, plus.w.band, minus.w.band, f, ...) {
  msdply(mspct = spct,
         plus.w.band = plus.w.band,
         minus.w.band = minus.w.band,
         f = f,
         ...)
}

#' reflectance:reflectance ratio
#'
#' This function returns the reflectance ratio for a given pair of wavebands of a
#' light reflector spectrum.
#'
#' @param spct an object of class "reflector_spct".
#' @param w.band.num waveband object or a list of waveband objects used to
#'   compute the numerator(s) and denominator(s) of the ratio(s).
#' @param w.band.denom waveband object or a list of waveband objects used to
#'   compute the denominator(s) of the ratio(s).
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param quantity character One of "total", "average" or "mean".
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param name.tag character Used to tag the name of the returned values.
#' @param ... other arguments (possibly ignored)
#'
#' @details With the default \code{quantity = "total"} the ratio is based on
#'   two \strong{photon irradiances}, one computed for each waveband.
#'
#'   \deqn{\frac{\mathrm{Rfr}(s, wb_\mathrm{num})}{\mathrm{Rfr}(s, wb_\mathrm{denom})}}
#'
#' If the argument is set to \code{quantity = "mean"} or
#'  \code{quantity = "average"} the ratio is based on
#'   two \strong{mean spectral photon irradiances}, one computed for each waveband.
#'
#'   \deqn{\frac{\bar{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{num})}{\bar{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{denom}))}}
#'
#' Only if the wavelength expanse of the two wavebands is the same, these two
#' ratios are numerically identical.
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector with name attribute set. The name is based on the name of the
#'   wavebands unless a named list of wavebands is supplied in which case the
#'   names of the list elements are used. "[Rfr:Rfr]" is appended if \code{quantity
#'   = "total"} and "[Rfr(wl):Rfr(wl)]" if \code{quantity = "mean"} or
#'   \code{quantity = "average"}.
#'
#'   A \code{data.frame} is returned in the case of collections of spectra,
#'   containing one column for each fraction definition, an index column with
#'   the names of the spectra, and optionally additional columns with metadata
#'   values retrieved from the attributes of the member spectra.
#'
#'   Fraction definitions are "assembled" from the arguments passed to
#'   \code{w.band.num} and \code{w.band.denom}. If both arguments are lists of
#'   waveband definitions, with an equal number of members, then the wavebands
#'   are paired to obtain as many fractions as the number of wavebands in each
#'   list. Recycling for wavebands takes place when the number of denominator
#'   and numerator wavebands differ.
#'
#' @export
#' @examples
#' Rfr_ratio(white_body.spct, new_waveband(400,500), new_waveband(400,700))
#'
#' @note The last two parameters control speed
#'   optimizations. The defaults should be suitable in most cases. If you will
#'   use repeatedly the same SWFs on many spectra measured at exactly the same
#'   wavelengths you may obtain some speed up by setting
#'   \code{use.cached.mult=TRUE}. However, be aware that you are responsible for
#'   ensuring that the wavelengths are the same in each call, as the only test
#'   done is for the length of the \code{w.length} vector.
#'
#' @family Reflectance ratio functions
#'
Rfr_ratio <- function(spct,
                      w.band.num,
                      w.band.denom,
                      scale.factor,
                      wb.trim,
                      use.cached.mult,
                      use.hinges,
                      ...) UseMethod("Rfr_fraction")

#' @describeIn Rfr_ratio Default for generic function
#'
#' @export
#'
Rfr_ratio.default <- function(spct,
                              w.band.num,
                              w.band.denom,
                              scale.factor,
                              wb.trim,
                              use.cached.mult,
                              use.hinges,
                              ...) {
  warning("'Rfr_ratio' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn Rfr_ratio Method for \code{reflector_spct} objects
#'
#' @export
#'
Rfr_ratio.reflector_spct <-
  function(spct,
           w.band.num = NULL,
           w.band.denom = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = NULL,
           quantity = "total",
           naming = "short",
           name.tag = NULL,
           ... ) {

    if (is.null(name.tag) && naming != "none") {
      if (quantity  == "total") {
        name.tag <- "[Rfr:Rfr]"
      } else {
        name.tag <- "[Rfr(wl):Rfr(wl)]"
      }
    }

    reflectances <-
      two_reflectances(spct = spct,
                       w.band.1 = w.band.num,
                       w.band.2 = w.band.denom,
                       quantity = quantity,
                       wb.trim = wb.trim,
                       use.cached.mult = use.cached.mult,
                       use.hinges = use.hinges,
                       naming = naming)

    Rfr.num <- reflectances[["Rfr.1"]]
    Rfr.denom <- reflectances[["Rfr.2"]]
    fraction <- Rfr.num / Rfr.denom * scale.factor
    names(fraction) <-
      paste(names(Rfr.num), ":", names(Rfr.num), name.tag, sep = "")
    setRfrType(fraction, getRfrType(spct))
    if (quantity == "total") {
      attr(fraction, "radiation.unit") <- "Rfr:Rfr fraction"
    } else {
      attr(fraction, "radiation.unit") <- "Rfr(wl):Rfr(wl) fraction"
    }
    return(fraction)
  }

#' reflectance:reflectance fraction
#'
#' This function returns the reflectance fraction for a given pair of wavebands of a
#' light reflector spectrum.
#'
#' @param spct an object of class "reflector_spct".
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
#' @param quantity character One of "total", "average" or "mean".
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param name.tag character Used to tag the name of the returned values.
#' @param ... other arguments (possibly ignored)
#'
#' @details With the default \code{quantity = "total"} the fraction is based on
#'   two \strong{photon irradiances}, one computed for each waveband.
#'
#'   \deqn{\frac{\mathrm{Rfr}(s, wb_\mathrm{num})}{\mathrm{Rfr}(s, wb_\mathrm{denom}) + \mathrm{Rfr}(s, wb_\mathrm{num})}}
#'
#' If the argument is set to \code{quantity = "mean"} or
#'  \code{quantity = "average"} the ratio is based on
#'   two \strong{mean spectral photon irradiances}, one computed for each waveband.
#'
#'   \deqn{\frac{\bar{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{num})}{\bar{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{denom}) + \bar{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{num})}}
#'
#' Only if the wavelength expanse of the two wavebands is the same, these two
#' ratios are numerically identical.
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector with name attribute set. The name is based on the name of the
#'   wavebands unless a named list of wavebands is supplied in which case the
#'   names of the list elements are used. "[Rfr:Rfr]" is appended if \code{quantity
#'   = "total"} and "[Rfr(wl):Rfr(wl)]" if \code{quantity = "mean"} or
#'   \code{quantity = "average"}.
#'
#'   A \code{data.frame} is returned in the case of collections of spectra,
#'   containing one column for each fraction definition, an index column with
#'   the names of the spectra, and optionally additional columns with metadata
#'   values retrieved from the attributes of the member spectra.
#'
#'   Fraction definitions are "assembled" from the arguments passed to
#'   \code{w.band.num} and \code{w.band.denom}. If both arguments are lists of
#'   waveband definitions, with an equal number of members, then the wavebands
#'   are paired to obtain as many fractions as the number of wavebands in each
#'   list. Recycling for wavebands takes place when the number of denominator
#'   and numerator wavebands differ.
#'
#' @export
#' @examples
#' Rfr_fraction(white_body.spct, new_waveband(400,500), new_waveband(400,700))
#'
#' @note The last two parameters control speed
#'   optimizations. The defaults should be suitable in most cases. If you will
#'   use repeatedly the same SWFs on many spectra measured at exactly the same
#'   wavelengths you may obtain some speed up by setting
#'   \code{use.cached.mult=TRUE}. However, be aware that you are responsible for
#'   ensuring that the wavelengths are the same in each call, as the only test
#'   done is for the length of the \code{w.length} vector.
#'
#' @family Reflectance ratio functions
#'
Rfr_fraction <- function(spct,
                         w.band.num,
                         w.band.denom,
                         scale.factor,
                         wb.trim,
                         use.cached.mult,
                         use.hinges,
                         ...) UseMethod("Rfr_fraction")

#' @describeIn Rfr_fraction Default for generic function
#'
#' @export
#'
Rfr_fraction.default <- function(spct,
                                 w.band.num,
                                 w.band.denom,
                                 scale.factor,
                                 wb.trim,
                                 use.cached.mult,
                                 use.hinges,
                                 ...) {
  warning("'Rfr_fraction' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn Rfr_fraction Method for \code{reflector_spct} objects
#'
#' @export
#'
Rfr_fraction.reflector_spct <-
  function(spct,
           w.band.num = NULL,
           w.band.denom = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = NULL,
           quantity = "total",
           naming = "short",
           name.tag = NULL,
           ... ) {

    if (is.null(name.tag) && naming != "none") {
      if (quantity  == "total") {
        name.tag <- "[Rfr:Rfr]"
      } else {
        name.tag <- "[Rfr(wl):Rfr(wl)]"
      }
    }

    reflectances <-
      two_reflectances(spct = spct,
                       w.band.1 = w.band.num,
                       w.band.2 = w.band.denom,
                       quantity = quantity,
                       wb.trim = wb.trim,
                       use.cached.mult = use.cached.mult,
                       use.hinges = use.hinges,
                       naming = naming)

    Rfr.num <- reflectances[["Rfr.1"]]
    Rfr.denom <- reflectances[["Rfr.2"]]
    fraction <- Rfr.num / (Rfr.denom + Rfr.num) * scale.factor
    names(fraction) <- paste(names(Rfr.num), ":(",
                             names(Rfr.num), "+", names(Rfr.denom), ")",
                             name.tag, sep = "")
    setRfrType(fraction, getRfrType(spct))
    if (quantity == "total") {
      attr(fraction, "radiation.unit") <- "Rfr:Rfr fraction"
    } else {
      attr(fraction, "radiation.unit") <- "Rfr(wl):Rfr(wl) fraction"
    }
    return(fraction)
  }


#' reflectance:reflectance NDI
#'
#' This function returns the reflectance normalized difference index for a given
#' pair of wavebands of a light reflector spectrum.
#'
#' @param spct an object of class "reflector_spct".
#' @param w.band.plus,w.band.minus waveband object(s) or a list(s) of waveband
#'   objects used to compute the additive and subtractive reflectance terms of
#'   the normalized difference index.
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param quantity character One of "total", "average" or "mean".
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param name.tag character Used to tag the name of the returned values.
#' @param ... other arguments (possibly ignored)
#'
#' @details With the default \code{quantity = "total"} the fraction is based on
#'   two \strong{photon reflectances}, one computed for each waveband.
#'
#'   \deqn{\frac{\mathrm{Rfr}(s, wb_\mathrm{plus}) - \mathrm{Rfr}(s, wb_\mathrm{minus})}{\mathrm{Rfr}(s, wb_\mathrm{plus}) + \mathrm{Rfr}(s, wb_\mathrm{minus})}}
#'
#' If the argument is set to \code{quantity = "mean"} or
#'  \code{quantity = "average"} the ratio is based on
#'   two \strong{mean spectral photon reflectances}, one computed for each waveband.
#'
#'   \deqn{\frac{\bar{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{plus}) - \bar{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{minus})}{\bar{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{plus}) + \bar{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{minus})}}
#'
#' Only if the wavelength expanse of the two wavebands is the same, these two
#' ratios are numerically identical.
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector with name attribute set. The name is based on the name of the
#'   wavebands unless a named list of wavebands is supplied in which case the
#'   names of the list elements are used. "[Rfr:Rfr]" is appended if \code{quantity
#'   = "total"} and "[Rfr(wl):Rfr(wl)]" if \code{quantity = "mean"} or
#'   \code{quantity = "average"}.
#'
#'   A \code{data.frame} is returned in the case of collections of spectra,
#'   containing one column for each fraction definition, an index column with
#'   the names of the spectra, and optionally additional columns with metadata
#'   values retrieved from the attributes of the member spectra.
#'
#'   Fraction definitions are "assembled" from the arguments passed to
#'   \code{w.band.num} and \code{w.band.denom}. If both arguments are lists of
#'   waveband definitions, with an equal number of members, then the wavebands
#'   are paired to obtain as many fractions as the number of wavebands in each
#'   list. Recycling for wavebands takes place when the number of denominator
#'   and numerator wavebands differ.
#'
#' @export
#' @examples
#' Rfr_normdiff(white_body.spct, new_waveband(400,500), new_waveband(400,700))
#'
#' @note The last two parameters control speed
#'   optimizations. The defaults should be suitable in most cases. If you will
#'   use repeatedly the same SWFs on many spectra measured at exactly the same
#'   wavelengths you may obtain some speed up by setting
#'   \code{use.cached.mult=TRUE}. However, be aware that you are responsible for
#'   ensuring that the wavelengths are the same in each call, as the only test
#'   done is for the length of the \code{w.length} vector.
#'
#' @family Reflectance ratio functions
#' @seealso \code{\link{normalized_diff_ind}}, accepts different summary
#' functions.
#'
Rfr_normdiff <- function(spct,
                         w.band.plus,
                         w.band.minus,
                         scale.factor,
                         wb.trim,
                         use.cached.mult,
                         use.hinges,
                         ...) UseMethod("Rfr_fraction")

#' @describeIn Rfr_normdiff Default for generic function
#'
#' @export
#'
Rfr_normdiff.default <- function(spct,
                                 w.band.plus,
                                 w.band.minus,
                                 scale.factor,
                                 wb.trim,
                                 use.cached.mult,
                                 use.hinges,
                                 ...) {
  warning("'Rfr_fraction' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn Rfr_normdiff Method for \code{reflector_spct} objects
#'
#' @export
#'
Rfr_normdiff.reflector_spct <-
  function(spct,
           w.band.plus = NULL,
           w.band.minus = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = NULL,
           quantity = "total",
           naming = "short",
           name.tag = NULL,
           ... ) {

    if (is.null(name.tag) && naming != "none") {
      if (quantity  == "total") {
        name.tag <- "[Rfr:Rfr]"
      } else {
        name.tag <- "[Rfr(wl):Rfr(wl)]"
      }
    }

    reflectances <-
      two_reflectances(spct = spct,
                       w.band.1 = w.band.plus,
                       w.band.2 = w.band.minus,
                       quantity = quantity,
                       wb.trim = wb.trim,
                       use.cached.mult = use.cached.mult,
                       use.hinges = use.hinges,
                       naming = naming)

    Rfr.plus <- reflectances[["Rfr.1"]]
    Rfr.minus <- reflectances[["Rfr.2"]]
    Rfr.NDI <- (Rfr.plus - Rfr.minus) / (Rfr.plus + Rfr.minus) * scale.factor
    names(Rfr.NDI) <- paste("(",
                            names(Rfr.plus), "-", names(Rfr.minus),
                            "):(",
                            names(Rfr.plus), "+", names(Rfr.minus), ")",
                            name.tag, sep = "")
    setRfrType(Rfr.NDI, getRfrType(spct))
    if (quantity == "total") {
      attr(Rfr.NDI, "radiation.unit") <- "Rfr:Rfr NDI"
    } else {
      attr(Rfr.NDI, "radiation.unit") <- "Rfr(wl):Rfr(wl) NDI"
    }
    return(Rfr.NDI)
  }


#' Compute two reflectances for ratio, fraction or index
#'
#' Internal function that computes the two reflectances needed to compute
#' various waveband ratios and fractions.
#'
#' @details See \code{\link{reflectance}} for details on the reflectance calculations.
#'
#' @param spct an object of class "reflector_spct" or "object_spct".
#' @param w.band.1,w.band.2 waveband objects or lists of waveband objects
#'   used to compute the numerator(s) and denominator(s) of the ratio(s). The
#'   waveband(s) determine the region(s) of the spectrum that are summarized.
#' @param quantity character string One of "total", "average" or "mean".
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#'
#' @keywords internal
#'
# This function is extremely simple but ensures consistency and avoids repetition
# It is used to define ratios, fractions and NDIs.
#
two_reflectances <- function(spct,
                             w.band.1,
                             w.band.2,
                             quantity,
                             wb.trim,
                             use.cached.mult,
                             use.hinges,
                             naming) {

  # we look for multiple spectra in long form
  num.spectra <- getMultipleWl(spct)
  if (num.spectra > 1) {
    message("Object contains ", num.spectra, " spectra in long form")
    # convert to a collection of spectra
    mspct <- subset2mspct(x = spct,
                          idx.var = getIdFactor(spct),
                          drop.idx = FALSE)
    # call method on the collection
    return(two_reflectances(spct = mspct,
                            w.band.1 = w.band.1,
                            w.band.2 = w.band.2,
                            quantity = quantity,
                            wb.trim = wb.trim,
                            use.cached.mult = use.cached.mult,
                            use.hinges = use.hinges,
                            naming = naming))
  }

  stopifnot("Unsupported argument passed to 'quantity'" =
              quantity %in% c("total", "average", "mean"))

  Rfr.1 <- reflectance(spct,
                       w.band = w.band.1,
                       quantity = quantity,
                       scale.factor = 1,
                       wb.trim = wb.trim,
                       use.cached.mult = use.cached.mult,
                       use.hinges = use.hinges,
                       allow.scaled = TRUE,
                       naming = naming)

  Rfr.2 <- reflectance(spct,
                       w.band = w.band.2,
                       quantity = quantity,
                       scale.factor = 1,
                       wb.trim = wb.trim,
                       use.cached.mult = use.cached.mult,
                       use.hinges = use.hinges,
                       allow.scaled = TRUE,
                       naming = naming)

  list(Rfr.1 = Rfr.1, Rfr.2 = Rfr.2)
}

