# Tfr_ratio() -------------------------------------------------------------

#' transmittance:transmittance ratio
#'
#' Transmittance ratio for a given pair of wavebands of a
#' filter spectrum.
#'
#' @param spct an object of class "filter_spct".
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
#' @details With the default \code{quantity = "mean"} or \code{quantity =
#'   "average"} the ratio is based on two \strong{mean spectral transmittance},
#'   one computed for each waveband.
#'
#'   \deqn{\frac{\overline{\mathrm{Tfr}_\lambda}(s, wb_\mathrm{num})}{\overline{\mathrm{Tfr}_\lambda}(s, wb_\mathrm{denom}))}}
#'
#' If the argument is set to \code{quantity = "total"} the ratio is based on
#'   two \strong{integrated transmittance}, one computed for each waveband.
#'
#'   \deqn{\frac{\mathrm{Tfr}(s, wb_\mathrm{num})}{\mathrm{Tfr}(s, wb_\mathrm{denom})}}
#'
#' Only if the wavelength expanse of the two wavebands is the same, these two
#' ratios are numerically identical.
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector with name attribute set. The name is based on the name of the
#'   wavebands unless a named list of wavebands is supplied in which case the
#'   names of the list elements are used. "[Tfr:Tfr]" is appended if
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
#' Tfr_ratio(Ler_leaf_rflt.spct,
#'           waveband(c(400,500), wb.name = "Blue"),
#'           waveband(c(600,700), wb.name = "Red"))
#' Tfr_ratio(Ler_leaf_rflt.spct,
#'           waveband(c(400,500), wb.name = "Blue"),
#'           waveband(c(600,700), wb.name = "Red"),
#'           quantity = "total")
#' Tfr_ratio(Ler_leaf_rflt.spct,
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
#' @family transmittance ratio functions
#'
Tfr_ratio <- function(spct,
                      w.band.num,
                      w.band.denom,
                      scale.factor,
                      wb.trim,
                      use.cached.mult,
                      use.hinges,
                      ...) UseMethod("Tfr_ratio")

#' @describeIn Tfr_ratio Default for generic function
#'
#' @export
#'
Tfr_ratio.default <- function(spct,
                              w.band.num,
                              w.band.denom,
                              scale.factor,
                              wb.trim,
                              use.cached.mult,
                              use.hinges,
                              ...) {
  warning("'Tfr_ratio' is not defined for objects of class ", class(spct)[1])
  NA_real_
}

#' @describeIn Tfr_ratio Method for \code{filter_spct} objects
#'
#' @export
#'
Tfr_ratio.filter_spct <-
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
           ... ) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(Tfr_ratio(spct = mspct,
                       w.band.num =  w.band.num,
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
        name.tag <- "[Tfr:Tfr]"
      } else {
        name.tag <- "[Tfr(wl):Tfr(wl)]"
      }
    }

    transmittances <-
      two_transmittances(spct = spct,
                         w.band.1 = w.band.num,
                         w.band.2 = w.band.denom,
                         quantity = quantity,
                         wb.trim = wb.trim,
                         use.cached.mult = use.cached.mult,
                         use.hinges = use.hinges,
                         naming = naming)

    Tfr.num <- transmittances[["Tfr.1"]]
    Tfr.denom <- transmittances[["Tfr.2"]]
    ratio <- Tfr.num / Tfr.denom * scale.factor
    names(ratio) <-
      paste(names(Tfr.num), ":", names(Tfr.denom), name.tag, sep = "")
    attr(ratio, "Tfr.type") <- getTfrType(spct)
    if (quantity == "total") {
      attr(ratio, "radiation.unit") <- "Tfr:Tfr ratio"
    } else {
      attr(ratio, "radiation.unit") <- "Tfr(wl):Tfr(wl) ratio"
    }
    return(ratio)
  }

#' @describeIn Tfr_ratio Calculates Tfr:Tfr from a \code{filter_mspct}
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
Tfr_ratio.filter_mspct <-
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
        .fun = Tfr_ratio.filter_spct,
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

# Tfr_fraction() ----------------------------------------------------------

#' transmittance:transmittance fraction
#'
#' Transmittance fraction for a given pair of wavebands of a
#' filter spectrum.
#'
#' @param spct an object of class "filter_spct".
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
#'   "average"} the ratio is based on two \strong{mean spectral transmittance},
#'   one computed for each waveband.
#'
#'   \deqn{\frac{\overline{\mathrm{Tfr}_\lambda}(s, wb_\mathrm{num})}{\overline{\mathrm{Tfr}_\lambda}(s, wb_\mathrm{denom}) + \overline{\mathrm{Tfr}_\lambda}(s, wb_\mathrm{num})}}
#'
#' If the argument is set to \code{quantity = "total"} the fraction is based on
#'   two \strong{integrated transmittance}, one computed for each waveband.
#'
#'   \deqn{\frac{\mathrm{Tfr}(s, wb_\mathrm{num})}{\mathrm{Tfr}(s, wb_\mathrm{denom}) + \mathrm{Tfr}(s, wb_\mathrm{num})}}
#'
#' Only if the wavelength expanse of the two wavebands is the same, these two
#' ratios are numerically identical.
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector with name attribute set. The name is based on the name of the
#'   wavebands unless a named list of wavebands is supplied in which case the
#'   names of the list elements are used. "[Tfr:Tfr]" is appended if \code{quantity
#'   = "total"} and "[Tfr(wl):Tfr(wl)]" if \code{quantity = "mean"} or
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
#' Tfr_fraction(Ler_leaf_rflt.spct,
#'              waveband(c(400,500), wb.name = "Blue"),
#'              waveband(c(600,700), wb.name = "Red"))
#' Tfr_fraction(Ler_leaf_rflt.spct,
#'              waveband(c(400,500), wb.name = "Blue"),
#'              waveband(c(600,700), wb.name = "Red"),
#'              quantity = "total")
#' Tfr_fraction(Ler_leaf_rflt.spct,
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
#' @family transmittance ratio functions
#'
Tfr_fraction <- function(spct,
                         w.band.num,
                         w.band.denom,
                         scale.factor,
                         wb.trim,
                         use.cached.mult,
                         use.hinges,
                         ...) UseMethod("Tfr_fraction")

#' @describeIn Tfr_fraction Default for generic function
#'
#' @export
#'
Tfr_fraction.default <- function(spct,
                                 w.band.num,
                                 w.band.denom,
                                 scale.factor,
                                 wb.trim,
                                 use.cached.mult,
                                 use.hinges,
                                 ...) {
  warning("'Tfr_fraction' is not defined for objects of class ", class(spct)[1])
  NA_real_
}

#' @describeIn Tfr_fraction Method for \code{filter_spct} objects
#'
#' @export
#'
Tfr_fraction.filter_spct <-
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
      return(Tfr_fraction(spct = mspct,
                          w.band.num =  w.band.num,
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
        name.tag <- "[Tfr:Tfr]"
      } else {
        name.tag <- "[Tfr(wl):Tfr(wl)]"
      }
    }

    transmittances <-
      two_transmittances(spct = spct,
                         w.band.1 = w.band.num,
                         w.band.2 = w.band.denom,
                         quantity = quantity,
                         wb.trim = wb.trim,
                         use.cached.mult = use.cached.mult,
                         use.hinges = use.hinges,
                         naming = naming)

    Tfr.num <- transmittances[["Tfr.1"]]
    Tfr.denom <- transmittances[["Tfr.2"]]
    fraction <- Tfr.num / (Tfr.denom + Tfr.num) * scale.factor
    names(fraction) <- paste(names(Tfr.num), ":(",
                             names(Tfr.num), "+", names(Tfr.denom), ")",
                             name.tag, sep = "")
    attr(fraction, "Tfr.type") <- getTfrType(spct)
    if (quantity == "total") {
      attr(fraction, "radiation.unit") <- "Tfr:Tfr fraction"
    } else {
      attr(fraction, "radiation.unit") <- "Tfr(wl):Tfr(wl) fraction"
    }
    return(fraction)
  }

#' @describeIn Tfr_fraction Calculates Tfr:Tfr from a \code{filter_mspct}
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
Tfr_fraction.filter_mspct <-
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
        .fun = Tfr_fraction.filter_spct,
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

# Tfr_normdiff() ----------------------------------------------------------

#' transmittance:transmittance normalised difference
#'
#' Transmittance normalized difference index for a given
#' pair of wavebands computed from a filter spectrum.
#'
#' @param spct an object of class "filter_spct".
#' @param w.band.plus,w.band.minus waveband object(s) or a list(s) of waveband
#'   objects used to compute the additive and subtractive transmittance terms of
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
#'  \code{quantity = "average"} the ratio is based on
#'  two \strong{mean spectral photon transmittances}, one computed for each
#'  waveband.
#'
#'   \deqn{\frac{\overline{\mathrm{Tfr}_\lambda}(s, wb_\mathrm{plus}) - \overline{\mathrm{Tfr}_\lambda}(s, wb_\mathrm{minus})}{\overline{\mathrm{Tfr}_\lambda}(s, wb_\mathrm{plus}) + \overline{\mathrm{Tfr}_\lambda}(s, wb_\mathrm{minus})}}
#'
#' If the argument is set to \code{quantity = "total"} the fraction is based on
#'   two \strong{photon transmittances}, one computed for each waveband.
#'
#'   \deqn{\frac{\mathrm{Tfr}(s, wb_\mathrm{plus}) - \mathrm{Tfr}(s, wb_\mathrm{minus})}{\mathrm{Tfr}(s, wb_\mathrm{plus}) + \mathrm{Tfr}(s, wb_\mathrm{minus})}}
#'
#' Only if the wavelength expanse of the two wavebands is the same, these two
#' ratios are numerically identical.
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector with name attribute set. The name is based on the name of the
#'   wavebands unless a named list of wavebands is supplied in which case the
#'   names of the list elements are used. "[Tfr:Tfr]" is appended if
#'   \code{quantity= "total"} and "[Tfr(wl):Tfr(wl)]" if
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
#' Tfr_normdiff(Ler_leaf_rflt.spct,
#'              waveband(c(400,500), wb.name = "Blue"),
#'              waveband(c(600,700), wb.name = "Red"))
#' Tfr_normdiff(Ler_leaf_rflt.spct,
#'              waveband(c(400,500), wb.name = "Blue"),
#'              waveband(c(600,700), wb.name = "Red"),
#'              quantity = "total")
#' Tfr_normdiff(Ler_leaf_rflt.spct,
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
#' @family transmittance ratio functions
#' @seealso \code{\link{normalized_diff_ind}}, accepts different summary
#' functions.
#'
Tfr_normdiff <- function(spct,
                         w.band.plus,
                         w.band.minus,
                         scale.factor,
                         wb.trim,
                         use.cached.mult,
                         use.hinges,
                         ...) UseMethod("Tfr_normdiff")

#' @describeIn Tfr_normdiff Default for generic function
#'
#' @export
#'
Tfr_normdiff.default <- function(spct,
                                 w.band.plus,
                                 w.band.minus,
                                 scale.factor,
                                 wb.trim,
                                 use.cached.mult,
                                 use.hinges,
                                 ...) {
  warning("'Tfr_fraction' is not defined for objects of class ", class(spct)[1])
  NA_real_
}

#' @describeIn Tfr_normdiff Method for \code{filter_spct} objects
#'
#' @export
#'
Tfr_normdiff.filter_spct <-
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
      return(Tfr_normdiff(spct = mspct,
                          w.band.plus =  w.band.plus,
                          w.band.minus = w.band.minus,
                          scale.factor = scale.factor,
                          wb.trim = wb.trim,
                          use.cached.mult = use.cached.mult,
                          use.hinges = use.hinges,
                          quantity = quantity,
                          naming = naming,
                          name.tag = name.tag,
                          ... ))
    }

    if (is.null(name.tag) && naming != "none") {
      if (quantity  == "total") {
        name.tag <- "[Tfr:Tfr]"
      } else {
        name.tag <- "[Tfr(wl):Tfr(wl)]"
      }
    }

    transmittances <-
      two_transmittances(spct = spct,
                         w.band.1 = w.band.plus,
                         w.band.2 = w.band.minus,
                         quantity = quantity,
                         wb.trim = wb.trim,
                         use.cached.mult = use.cached.mult,
                         use.hinges = use.hinges,
                         naming = naming)

    Tfr.plus <- transmittances[["Tfr.1"]]
    Tfr.minus <- transmittances[["Tfr.2"]]
    Tfr.normdiff <-
      (Tfr.plus - Tfr.minus) / (Tfr.plus + Tfr.minus) * scale.factor
    names(Tfr.normdiff) <- paste("(",
                                 names(Tfr.plus), "-", names(Tfr.minus),
                                 "):(",
                                 names(Tfr.plus), "+", names(Tfr.minus), ")",
                                 name.tag, sep = "")
    attr(Tfr.normdiff, "Tfr.type") <- getTfrType(spct)
    if (quantity == "total") {
      attr(Tfr.normdiff, "radiation.unit") <- "Tfr:Tfr normdiff"
    } else {
      attr(Tfr.normdiff, "radiation.unit") <- "Tfr(wl):Tfr(wl) normdiff"
    }
    return(Tfr.normdiff)
  }

#' @describeIn Tfr_normdiff Calculates Tfr:Tfr from a \code{filter_mspct}
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
Tfr_normdiff.filter_mspct <-
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
        .fun = Tfr_normdiff.filter_spct,
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

#' Compute two transmittances for ratio, fraction or normalised difference
#'
#' Internal function that computes the two transmittances needed to compute
#' various waveband ratios and fractions.
#'
#' @details # This function is extremely simple but ensures consistency and
#'   avoids repetition. It is used to define ratios, fractions and NDIs.
#' @seealso See \code{\link{transmittance}} for details on the transmittance
#'   calculations.
#'
#' @param spct an object of class "filter_spct" or "object_spct".
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
two_transmittances <- function(spct,
                               w.band.1,
                               w.band.2,
                               quantity,
                               wb.trim,
                               use.cached.mult,
                               use.hinges,
                               naming) {

  stopifnot("Unsupported argument passed to 'quantity'" =
              quantity %in% c("total", "average", "mean"))

  Tfr.1 <- transmittance(spct,
                         w.band = w.band.1,
                         quantity = quantity,
                         scale.factor = 1,
                         wb.trim = wb.trim,
                         use.cached.mult = use.cached.mult,
                         use.hinges = use.hinges,
                         allow.scaled = TRUE,
                         naming = naming)

  Tfr.2 <- transmittance(spct,
                         w.band = w.band.2,
                         quantity = quantity,
                         scale.factor = 1,
                         wb.trim = wb.trim,
                         use.cached.mult = use.cached.mult,
                         use.hinges = use.hinges,
                         allow.scaled = TRUE,
                         naming = naming)

  list(Tfr.1 = Tfr.1, Tfr.2 = Tfr.2)
}
