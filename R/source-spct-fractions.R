# q_fraction() ---------------------------------------------------------------

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
#' @param quantity character One of "total", "average" or "mean".
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param name.tag character Used to tag the name of the returned values.
#' @param ... other arguments (possibly ignored)
#'
#' @details With the default \code{quantity = "total"} the fraction is based on
#'   two \strong{photon irradiances}, one computed for each waveband.
#'
#'   \deqn{\frac{Q(s, wb_\mathrm{num})}{Q(s, wb_\mathrm{denom}) + Q(s, wb_\mathrm{num})}}
#'
#' If the argument is set to \code{quantity = "mean"} or
#'  \code{quantity = "average"} the ratio is based on two
#'  \strong{mean spectral photon irradiances}, one computed for each waveband.
#'
#'   \deqn{\frac{\overline{Q_\lambda}(s, wb_\mathrm{num})}{\overline{Q_\lambda}(s, wb_\mathrm{denom}) + \overline{Q_\lambda}(s, wb_\mathrm{num})}}
#'
#' Only if the wavelength expanse of the two wavebands is the same, these two
#' ratios are numerically identical.
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector with name attribute set. The name is based on the name of the
#'   wavebands unless a named list of wavebands is supplied in which case the
#'   names of the list elements are used. "[q:q]" is appended if \code{quantity
#'   = "total"} and "[q(wl):q(wl)]" if \code{quantity = "mean"} or
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
q_fraction <- function(spct,
                       w.band.num,
                       w.band.denom,
                       scale.factor,
                       wb.trim,
                       use.cached.mult,
                       use.hinges,
                       ...) UseMethod("q_fraction")

#' @describeIn q_fraction Default for generic function
#'
#' @export
#'
q_fraction.default <- function(spct,
                               w.band.num,
                               w.band.denom,
                               scale.factor,
                               wb.trim,
                               use.cached.mult,
                               use.hinges,
                               ...) {
  warning("'q_fraction' is not defined for objects of class ", class(spct)[1])
  NA_real_
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
           quantity = "total",
           naming = "short",
           name.tag = NULL,
           ...) {

    # we look for multiple spectra in long form
    num.spectra <- getMultipleWl(spct)
    if (num.spectra > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(q_fraction(spct = mspct,
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
        name.tag <- "[q:q]"
      } else {
        name.tag <- "[q(wl):q(wl)]"
      }
    }

    irrads <- two_irrads(spct = spct,
                         w.band.num = w.band.num,
                         w.band.denom = w.band.denom,
                         unit.out.num = "photon",
                         unit.out.denom = "photon",
                         quantity = quantity,
                         wb.trim = wb.trim,
                         use.cached.mult = use.cached.mult,
                         use.hinges = use.hinges,
                         naming = naming)

    q.irrad.num <- irrads[["irrad.num"]]
    q.irrad.denom <- irrads[["irrad.denom"]]
    fraction <- q.irrad.num / (q.irrad.denom + q.irrad.num) * scale.factor
    names(fraction) <- paste(names(q.irrad.num), ":(",
                             names(q.irrad.num), "+", names(q.irrad.denom), ")",
                             name.tag, sep = "")
    attr(fraction, "time.unit") <- NULL
    if (quantity == "total") {
      attr(fraction, "radiation.unit") <- "q:q fraction"
    } else {
      attr(fraction, "radiation.unit") <- "q(wl):q(wl) fraction"
    }
    return(fraction)
  }

#' @describeIn q_fraction Calculates photon:photon from a \code{source_mspct}
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
q_fraction.source_mspct <-
  function(spct,
           w.band.num = NULL,
           w.band.denom = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = NULL,
           quantity = "total",
           naming = "short",
           name.tag = ifelse(naming != "none", "[q:q]", ""),
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
        .fun = q_fraction.source_spct,
        w.band.num = w.band.num,
        w.band.denom = w.band.denom,
        quantity = quantity,
        wb.trim = wb.trim,
        scale.factor = scale.factor,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
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

# e_fraction() ---------------------------------------------------------------

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
#' @param use.cached.mult logical Flag telling whether multiplier values should
#'   be cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param quantity character One of "total", "average" or "mean".
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param name.tag character Used to tag the name of the returned values.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @details With the default \code{quantity = "total"} the fraction is based on
#'   two \strong{energy irradiances}, one computed for each waveband.
#'
#'   \deqn{\frac{E(s, wb_\mathrm{num})}{E(s, wb_\mathrm{denom}) + E(s, wb_\mathrm{num})}}
#'
#' If the argument is set to \code{quantity = "mean"} or
#'  \code{quantity = "average"} the ratio is based on two
#'  \strong{mean spectral energy irradiances}, one computed for each waveband.
#'
#'   \deqn{\frac{\overline{Q_\lambda}(s, wb_\mathrm{num})}{\overline{Q_\lambda}(s, wb_\mathrm{denom}) + \overline{Q_\lambda}(s, wb_\mathrm{num})}}
#'
#' Only if the wavelength expanse of the two wavebands is the same, these two
#' ratios are numerically identical.
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector with name attribute set. The name is based on the name of the
#'   wavebands unless a named list of wavebands is supplied in which case the
#'   names of the list elements are used. "[e:e]" is appended if \code{quantity
#'   = "total"} and "[e(wl):e(wl)]" if \code{quantity = "mean"} or
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
#' @family photon and energy ratio functions
#'
e_fraction <- function(spct,
                       w.band.num,
                       w.band.denom,
                       scale.factor,
                       wb.trim,
                       use.cached.mult,
                       use.hinges,
                       ...) UseMethod("e_fraction")

#' @describeIn e_fraction Default for generic function
#'
#' @export
#'
e_fraction.default <- function(spct,
                               w.band.num,
                               w.band.denom,
                               scale.factor,
                               wb.trim,
                               use.cached.mult,
                               use.hinges,
                               ...) {
  warning("'e_fraction' is not defined for objects of class ", class(spct)[1])
  NA_real_
}

#' @describeIn e_fraction Method for \code{source_spct} objects
#'
#' @export
#'
e_fraction.source_spct <-
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
           ...) {

    # we look for multiple spectra in long form
    num.spectra <- getMultipleWl(spct)
    if (num.spectra > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(e_fraction(spct = mspct,
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
        name.tag <- "[e:e]"
      } else {
        name.tag <- "[e(wl):e(wl)]"
      }
    }

    irrads <- two_irrads(spct = spct,
                         w.band.num = w.band.num,
                         w.band.denom = w.band.denom,
                         unit.out.num = "energy",
                         unit.out.denom = "energy",
                         quantity = quantity,
                         wb.trim = wb.trim,
                         use.cached.mult = use.cached.mult,
                         use.hinges = use.hinges,
                         naming = naming)

    e.irrad.num <- irrads[["irrad.num"]]
    e.irrad.denom <- irrads[["irrad.denom"]]
    fraction <- e.irrad.num / (e.irrad.denom + e.irrad.num) * scale.factor
    names(fraction) <- paste(names(e.irrad.num), ":(",
                             names(e.irrad.num), "+", names(e.irrad.denom), ")",
                             name.tag, sep = "")
    attr(fraction, "time.unit") <- NULL
    if (quantity == "total") {
      attr(fraction, "radiation.unit") <- "e:e fraction"
    } else {
      attr(fraction, "radiation.unit") <- "e(wl):e(wl) fraction"
    }
    return(fraction)
  }

#' @describeIn e_fraction Calculates energy:energy fraction from a
#'   \code{source_mspct} object.
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax
#'   for \code{attr2tb} passed as is to formal parameter \code{col.names}.
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach.
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
e_fraction.source_mspct <-
  function(spct,
           w.band.num = NULL,
           w.band.denom = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = NULL,
           quantity = "total",
           naming = "short",
           name.tag = ifelse(naming != "none", "[e:e]", ""),
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
        .fun = e_fraction.source_spct,
        w.band.num = w.band.num,
        w.band.denom = w.band.denom,
        quantity = quantity,
        wb.trim = wb.trim,
        scale.factor = scale.factor,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
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
