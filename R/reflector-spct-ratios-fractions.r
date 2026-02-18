# Rfr_ratio() -------------------------------------------------------------

#' reflectance:reflectance ratio
#'
#' This function returns the reflectance ratio for a given pair of wavebands
#' of a reflector spectrum.
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
#' @details With the default \code{quantity = "mean"} or
#'  \code{quantity = "average"} the ratio is based on
#'   two \strong{mean spectral reflectance}, one computed for each waveband.
#'
#'   \deqn{\frac{\overline{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{num})}{\overline{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{denom}))}}
#'
#' If the argument is set to \code{quantity = "total"} the ratio is based on
#'   two \strong{integrated reflectance}, one computed for each waveband.
#'
#'   \deqn{\frac{\mathrm{Rfr}(s, wb_\mathrm{num})}{\mathrm{Rfr}(s, wb_\mathrm{denom})}}
#'
#' Only if the wavelength expanse of the two wavebands is the same, these two
#' ratios are numerically identical.
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector with name attribute set. The name is based on the name of the
#'   wavebands unless a named list of wavebands is supplied in which case the
#'   names of the list elements are used. "[Rfr:Rfr]" is appended if
#'   \code{quantity = "total"} and "[Rfr(wl):Rfr(wl)]" if
#'   \code{quantity = "mean"} or \code{quantity = "average"}.
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
#' Rfr_ratio(Ler_leaf_rflt.spct,
#'           waveband(c(400,500), wb.name = "Blue"),
#'           waveband(c(600,700), wb.name = "Red"))
#' Rfr_ratio(Ler_leaf_rflt.spct,
#'           waveband(c(400,500), wb.name = "Blue"),
#'           waveband(c(600,700), wb.name = "Red"),
#'           quantity = "total")
#' Rfr_ratio(Ler_leaf_rflt.spct,
#'           waveband(c(400,500), wb.name = "Blue"),
#'           waveband(c(600,700), wb.name = "Red"),
#'           quantity = "mean")
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
                      ...) UseMethod("Rfr_ratio")

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
  NA_real_
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
           quantity = "mean",
           naming = "short",
           name.tag = NULL,
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(Rfr_ratio(spct = mspct,
                       w.band.num = w.band.num,
                       w.band.denom = w.band.denom,
                       scale.factor = scale.factor,
                       wb.trim = wb.trim,
                       use.cached.mult = use.cached.mult,
                       use.hinges = use.hinges,
                       quantity = quantity,
                       naming = naming,
                       name.tag = name.tag,
                       ...))
    }

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
    ratio <- Rfr.num / Rfr.denom * scale.factor
    names(ratio) <-
      paste(names(Rfr.num), ":", names(Rfr.denom), name.tag, sep = "")
    attr(ratio, "Rfr.type") <- getRfrType(spct)
    if (quantity == "total") {
      attr(ratio, "radiation.unit") <- "Rfr:Rfr ratio"
    } else {
      attr(ratio, "radiation.unit") <- "Rfr(wl):Rfr(wl) ratio"
    }
    return(ratio)
  }

#' @describeIn Rfr_ratio Calculates Rfr:Rfr from a \code{reflector_mspct}
#'   object.
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax
#'   for \code{attr2tb} passed as is to formal parameter \code{col.names}.
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
Rfr_ratio.reflector_mspct <-
  function(spct,
           w.band.num = NULL,
           w.band.denom = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = NULL,
           quantity = "mean",
           naming = "short",
           name.tag = NULL,
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {
    if (naming == "none") {
      # need names for columns
      naming <- "short"
    }

    spct <- subset2mspct(spct) # expand long form spectra within collection

    z <-
      msdply(
        mspct = spct,
        .fun = Rfr_ratio.reflector_spct,
        w.band.num = w.band.num,
        w.band.denom = w.band.denom,
        wb.trim = wb.trim,
        scale.factor = scale.factor,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
        quantity = quantity,
        naming = naming,
        name.tag = name.tag,
        idx = idx,
        .parallel = .parallel,
        .paropts = .paropts
      )
    add_attr2tb(tb = z,
                mspct = spct,
                col.names = attr2tb,
                idx = idx)
  }

# Rfr_fraction() ----------------------------------------------------------

#' reflectance:reflectance fraction
#'
#' This function returns the reflectance fraction for a given pair of wavebands
#' of a reflector spectrum.
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
#' @details With the default \code{quantity = "mean"} or \code{quantity =
#'   "average"} the ratio is based on two \strong{mean spectral reflectance},
#'   one computed for each waveband.
#'
#'   \deqn{\frac{\overline{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{num})}{\overline{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{denom}) + \overline{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{num})}}
#'
#' If the argument is set to \code{quantity = "total"} the fraction is based on
#'   two \strong{integrated reflectance}, one computed for each waveband.
#'
#'   \deqn{\frac{\mathrm{Rfr}(s, wb_\mathrm{num})}{\mathrm{Rfr}(s, wb_\mathrm{denom}) + \mathrm{Rfr}(s, wb_\mathrm{num})}}
#'
#' Only if the wavelength expanse of the two wavebands is the same, these two
#' ratios are numerically identical.
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector with name attribute set. The name is based on the name of the
#'   wavebands unless a named list of wavebands is supplied in which case the
#'   names of the list elements are used. "[Rfr:Rfr]" is appended if
#'   \code{quantity = "total"} and "[Rfr(wl):Rfr(wl)]" if
#'   \code{quantity = "mean"} or \code{quantity = "average"}.
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
#' Rfr_fraction(Ler_leaf_rflt.spct,
#'              waveband(c(400,500), wb.name = "Blue"),
#'              waveband(c(600,700), wb.name = "Red"))
#' Rfr_fraction(Ler_leaf_rflt.spct,
#'              waveband(c(400,500), wb.name = "Blue"),
#'              waveband(c(600,700), wb.name = "Red"),
#'              quantity = "total")
#' Rfr_fraction(Ler_leaf_rflt.spct,
#'              waveband(c(400,500), wb.name = "Blue"),
#'              waveband(c(600,700), wb.name = "Red"),
#'              quantity = "mean")
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
  NA_real_
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
           quantity = "mean",
           naming = "short",
           name.tag = NULL,
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(Rfr_ratio(spct = mspct,
                       w.band.num = w.band.num,
                       w.band.denom = w.band.denom,
                       scale.factor = scale.factor,
                       wb.trim = wb.trim,
                       use.cached.mult = use.cached.mult,
                       use.hinges = use.hinges,
                       quantity = quantity,
                       naming = naming,
                       name.tag = name.tag,
                       ...))
    }

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
    attr(fraction, "Rfr.type") <- getRfrType(spct)
    if (quantity == "total") {
      attr(fraction, "radiation.unit") <- "Rfr:Rfr fraction"
    } else {
      attr(fraction, "radiation.unit") <- "Rfr(wl):Rfr(wl) fraction"
    }
    return(fraction)
  }

#' @describeIn Rfr_fraction Calculates Rfr:Rfr from a \code{reflector_mspct}
#'   object.
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax
#'   for \code{attr2tb} passed as is to formal parameter \code{col.names}.
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
Rfr_fraction.reflector_mspct <-
  function(spct,
           w.band.num = NULL,
           w.band.denom = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = NULL,
           quantity = "mean",
           naming = "short",
           name.tag = NULL,
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {
    if (naming == "none") {
      # need names for columns
      naming <- "short"
    }

    spct <- subset2mspct(spct) # expand long form spectra within collection

    z <-
      msdply(
        mspct = spct,
        .fun = Rfr_fraction.reflector_spct,
        w.band.num = w.band.num,
        w.band.denom = w.band.denom,
        wb.trim = wb.trim,
        scale.factor = scale.factor,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
        quantity = quantity,
        naming = naming,
        name.tag = name.tag,
        idx = idx,
        .parallel = .parallel,
        .paropts = .paropts
      )
    add_attr2tb(tb = z,
                mspct = spct,
                col.names = attr2tb,
                idx = idx)
  }


# Rfr_normdiff() ----------------------------------------------------------

#' reflectance:reflectance normalised difference
#'
#' This function returns the reflectance normalized difference index for a given
#' pair of wavebands of a reflector spectrum.
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
#' @details With the default \code{quantity = "mean"} or
#'   \code{quantity = "average"} the ratio is based on two values of
#'   \strong{mean spectral photon reflectance}, one computed for each waveband.
#'
#'   \deqn{\frac{\overline{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{plus}) - \overline{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{minus})}{\overline{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{plus}) + \overline{\mathrm{Rfr}_\lambda}(s, wb_\mathrm{minus})}}
#'
#' If the argument is set to \code{quantity = "total"} the fraction is based on
#'   two \strong{photon reflectances}, one computed for each waveband.
#'
#'   \deqn{\frac{\mathrm{Rfr}(s, wb_\mathrm{plus}) - \mathrm{Rfr}(s, wb_\mathrm{minus})}{\mathrm{Rfr}(s, wb_\mathrm{plus}) + \mathrm{Rfr}(s, wb_\mathrm{minus})}}
#'
#' Only if the wavelength expanse of the two wavebands is the same, these two
#' ratios are numerically identical.
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector with name attribute set. The name is based on the name of the
#'   wavebands unless a named list of wavebands is supplied in which case the
#'   names of the list elements are used. "[Rfr:Rfr]" is appended if
#'   \code{quantity = "total"} and "[Rfr(wl):Rfr(wl)]" if
#'   \code{quantity = "mean"} or \code{quantity = "average"}.
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
#' Rfr_normdiff(Ler_leaf_rflt.spct,
#'              waveband(c(400,500), wb.name = "Blue"),
#'              waveband(c(600,700), wb.name = "Red"))
#' Rfr_normdiff(Ler_leaf_rflt.spct,
#'              waveband(c(400,500), wb.name = "Blue"),
#'              waveband(c(600,700), wb.name = "Red"),
#'              quantity = "total")
#' Rfr_normdiff(Ler_leaf_rflt.spct,
#'              waveband(c(400,500), wb.name = "Blue"),
#'              waveband(c(600,700), wb.name = "Red"),
#'              quantity = "mean")
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult =T RUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
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
                         ...) UseMethod("Rfr_normdiff")

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
  NA_real_
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
           quantity = "mean",
           naming = "short",
           name.tag = NULL,
           ... ) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(Rfr_normdiff(spct = mspct,
                          w.band.plus = w.band.plus,
                          w.band.minus = w.band.minus,
                          scale.factor = scale.factor,
                          wb.trim = wb.trim,
                          use.cached.mult = use.cached.mult,
                          use.hinges = use.hinges,
                          quantity = quantity,
                          naming = naming,
                          name.tag = name.tag,
                          ...))
    }

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
    Rfr.normdiff <-
      (Rfr.plus - Rfr.minus) / (Rfr.plus + Rfr.minus) * scale.factor
    names(Rfr.normdiff) <-
      paste("(",
            names(Rfr.plus), "-", names(Rfr.minus), "):(",
            names(Rfr.plus), "+", names(Rfr.minus), ")",
            name.tag, sep = "")
    attr(Rfr.normdiff, "Rfr.type") <- getRfrType(spct)
    if (quantity == "total") {
      attr(Rfr.normdiff, "radiation.unit") <- "Rfr:Rfr normdiff"
    } else {
      attr(Rfr.normdiff, "radiation.unit") <- "Rfr(wl):Rfr(wl) normdiff"
    }
    return(Rfr.normdiff)
  }

#' @describeIn Rfr_normdiff Calculates Rfr:Rfr from a \code{reflector_mspct}
#'   object.
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax
#'   for \code{attr2tb} passed as is to formal parameter \code{col.names}.
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
Rfr_normdiff.reflector_mspct <-
  function(spct,
           w.band.plus = NULL,
           w.band.minus = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = NULL,
           quantity = "mean",
           naming = "short",
           name.tag = NULL,
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {
    if (naming == "none") {
      # need names for columns
      naming <- "short"
    }

    spct <- subset2mspct(spct) # expand long form spectra within collection

    z <-
      msdply(
        mspct = spct,
        .fun = Rfr_normdiff.reflector_spct,
        w.band.plus = w.band.plus,
        w.band.minus = w.band.minus,
        wb.trim = wb.trim,
        scale.factor = scale.factor,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
        quantity = quantity,
        naming = naming,
        name.tag = name.tag,
        idx = idx,
        .parallel = .parallel,
        .paropts = .paropts
      )
    add_attr2tb(tb = z,
                mspct = spct,
                col.names = attr2tb,
                idx = idx)
  }

# internal utility function --------------------------------------------------

#' Compute two reflectances for ratio, fraction or normalised difference
#'
#' Internal function that computes the two reflectances needed to compute
#' various waveband ratios and fractions.
#'
#' @details See \code{\link{reflectance}} for details on the reflectance
#'    calculations.
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
# This function is extremely simple but ensures consistency and avoids
# repetition. It is used to define ratios, fractions and NDIs.
#
two_reflectances <- function(spct,
                             w.band.1,
                             w.band.2,
                             quantity,
                             wb.trim,
                             use.cached.mult,
                             use.hinges,
                             naming) {

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
