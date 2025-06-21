# q_ratio() ---------------------------------------------------------------

#' Photon:photon ratio
#'
#' This function returns the photon ratio for a given pair of wavebands of a
#' light source spectrum.
#'
#' @param spct an object of class "source_spct".
#' @param w.band.num waveband object or a list of waveband objects used to
#'   compute the numerator(s) of the ratio(s).
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
#'   two photon irradiances, one computed for each waveband.
#'
#'   \deqn{\frac{Q(s, wb_\mathrm{num})}{Q(s, wb_\mathrm{denom})}}
#'
#' If the argument is set to \code{quantity = "mean"} or
#'  \code{quantity = "average"} the ratio is based on
#'   two mean spectral photon irradiances, one computed for each waveband.
#'
#'   \deqn{\frac{\overline{Q_\lambda}(s, wb_\mathrm{num})}{\overline{Q_\lambda}(s, wb_\mathrm{denom})}}
#'
#'   Ratios based on totals and means are numerically identical only if the
#'   wavelength expanse of the two wavebands is the same.
#'
#'   Fraction definitions are "assembled" from the arguments passed to
#'   \code{w.band.num} and \code{w.band.denom}. If both arguments are lists of
#'   waveband definitions, with an equal number of members, then the wavebands
#'   are paired to obtain as many fractions as the number of wavebands in each
#'   list. Recycling for wavebands takes place when the number of denominator
#'   and numerator wavebands differ.
#'
#'   The last two parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
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
#' @export
#' @examples
#' q_ratio(sun.spct,
#'         waveband(c(400,500), wb.name = "Blue"),
#'         waveband(c(400,700), wb.name = "White"))
#'
#' @section Performance: As this method accepts spectra as its input, it
#'   computes irradiances before computing the ratios. If you need to compute
#'   both ratios and irradiances from several hundreds or thousands of spectra,
#'   computing the ratios from previously computed irradiances avoids their
#'   repeated computation. A less dramatic, but still important, increase in
#'   performance is available when computing in the same function call ratios
#'   that share the same denominator.
#'
#' @family photon and energy ratio functions
#'
q_ratio <- function(spct,
                    w.band.num,
                    w.band.denom,
                    scale.factor,
                    wb.trim,
                    use.cached.mult,
                    use.hinges,
                    ...) UseMethod("q_ratio")

#' @describeIn q_ratio Default for generic function
#'
#' @export
#'
q_ratio.default <- function(spct,
                            w.band.num,
                            w.band.denom,
                            scale.factor,
                            wb.trim,
                            use.cached.mult,
                            use.hinges,
                            ...) {
  warning("'q_ratio' is not defined for objects of class ", class(spct)[1])
  NA_real_
}

#' @describeIn q_ratio Method for \code{source_spct} objects
#'
#' @export
#'
q_ratio.source_spct <-
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
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(q_ratio(spct = mspct,
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
    ratio <- q.irrad.num / q.irrad.denom * scale.factor
    names(ratio) <-
      paste(names(q.irrad.num), ":", names(q.irrad.denom), name.tag, sep = "")
    attr(ratio, "time.unit") <- NULL
    if (quantity == "total") {
      attr(ratio, "radiation.unit") <- "q:q ratio"
    } else {
      attr(ratio, "radiation.unit") <- "q(wl):q(wl) ratio"
    }
    return(ratio)
  }

#' @describeIn q_ratio Calculates photon:photon from a \code{source_mspct}
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
q_ratio.source_mspct <-
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
        .fun = q_ratio,
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

# e_ratio() ---------------------------------------------------------------

#' Energy:energy ratio
#'
#' This function returns the photon ratio for a given pair of wavebands of a
#' light source spectrum.
#'
#' @param spct source_spct
#' @param w.band.num waveband object or a list of waveband objects used to
#'   compute the numerator(s) of the ratio(s).
#' @param w.band.denom waveband object or a list of waveband objects used to
#'   compute the denominator(s) of the ratio(s).
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
#' @details With the default \code{quantity = "total"} the ratio is based on
#'   two energy irradiances, one computed for each waveband.
#'
#'   \deqn{\frac{I(s, wb_\mathrm{num})}{I(s, wb_\mathrm{denom})}}
#'
#' If the argument is set to \code{quantity = "mean"} or
#'  \code{quantity = "average"} the ratio is based on
#'   two mean spectral photon irradiances, one computed for each waveband.
#'
#'   \deqn{\frac{\overline{I_\lambda}(s, wb_\mathrm{num})}{\overline{I_\lambda}(s, wb_\mathrm{denom})}}
#'
#' Only if the wavelength expanse of the two wavebands is the same, these two
#' ratios are numerically identical.
#'
#'   Fraction definitions are "assembled" from the arguments passed to
#'   \code{w.band.num} and \code{w.band.denom}. If both arguments are lists of
#'   waveband definitions, with an equal number of members, then the wavebands
#'   are paired to obtain as many fractions as the number of wavebands in each
#'   list. Recycling for wavebands takes place when the number of denominator
#'   and numerator wavebands differ.
#'
#'   The last two parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
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
#' @export
#' @examples
#' e_ratio(sun.spct,
#'         waveband(c(400,500), wb.name = "Blue"),
#'         waveband(c(400,700), wb.name = "White"))
#'
#' @section Performance: As this method accepts spectra as its input, it
#'   computes irradiances before computing the ratios. If you need to compute
#'   both ratios and irradiances from several hundreds or thousands of spectra,
#'   computing the ratios from previously computed irradiances avoids their
#'   repeated computation. A less dramatic, but still important, increase in
#'   performance is available when computing in the same function call ratios
#'   that share the same denominator.
#'
#' @family photon and energy ratio functions
#'
e_ratio <- function(spct,
                    w.band.num,
                    w.band.denom,
                    scale.factor,
                    wb.trim,
                    use.cached.mult,
                    use.hinges, ...) UseMethod("e_ratio")

#' @describeIn e_ratio Default for generic function
#'
#' @export
#'
e_ratio.default <- function(spct,
                            w.band.num,
                            w.band.denom,
                            scale.factor,
                            wb.trim,
                            use.cached.mult,
                            use.hinges,
                            ...) {
  warning("'e_ratio' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn e_ratio Method for \code{source_spct} objects
#'
#' @export
#'
e_ratio.source_spct <-
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
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(e_ratio(spct = mspct,
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
                         quantity = quantity,
                         unit.out.num = "energy",
                         unit.out.denom = "energy",
                         wb.trim = wb.trim,
                         use.cached.mult = use.cached.mult,
                         use.hinges = use.hinges,
                         naming = naming)

    e.irrad.num <- irrads[["irrad.num"]]
    e.irrad.denom <- irrads[["irrad.denom"]]
    ratio <- e.irrad.num / e.irrad.denom * scale.factor
    names(ratio) <-
      paste(names(e.irrad.num), ":", names(e.irrad.denom), name.tag, sep = "")
    attr(ratio, "time.unit") <- NULL
    if (quantity == "total") {
      attr(ratio, "radiation.unit") <- "e:e ratio"
    } else {
      attr(ratio, "radiation.unit") <- "e(wl):e(wl) ratio"
    }
    return(ratio)
  }

#' @describeIn e_ratio Calculates energy:energy ratio from a \code{source_mspct}
#'   object.
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
e_ratio.source_mspct <-
  function(spct,
           w.band.num = NULL, w.band.denom = NULL,
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
        .fun = e_ratio,
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

# qe_ratio() --------------------------------------------------------------

#' Photon:energy ratio
#'
#' This function returns the photon to energy ratio for each waveband of a light
#' source spectrum.
#'
#' @param spct source_spct.
#' @param w.band waveband or list of waveband objects.
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded.
#' @param use.cached.mult logical Flag telling whether multiplier values should
#'   be cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param name.tag character Used to tag the name of the returned values.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @details The ratio is based on one photon irrandiance and one energy
#'   irradiance, both computed for the same waveband.
#'
#'   \deqn{\frac{Q(s, wb)}{I(s, wb)}}
#'
#'   The last two parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @return Computed values are ratios between photon irradiance and energy
#'   irradiance for a given waveband. A named \code{numeric} vector in the case
#'   of methods for individual spectra, with one value for each \code{waveband}
#'   passed to parameter \code{w.band}. A \code{data.frame} in the case of
#'   collections of spectra, containing one column for each \code{waveband}
#'   object, an index column with the names of the spectra, and optionally
#'   additional columns with metadata values retrieved from the attributes of
#'   the member spectra.
#'
#'   By default values are only integrated, but depending on the argument passed
#'   to parameter \code{quantity} they can be re-expressed as relative fractions
#'   or percentages. In the case of vector output, \code{names} attribute is set
#'   to the name of the corresponding waveband unless a named list is supplied
#'   in which case the names of the list members are used, with "[q:e]"
#'   prepended. Units are [mol J-1].
#'
#' @export
#' @examples
#' qe_ratio(sun.spct,
#'          waveband(c(400,700), wb.name = "White")) # mol J-1
#' qe_ratio(sun.spct,
#'          waveband(c(400,700), wb.name = "White"),
#'          scale.factor = 1e6) # umol J-1
#'
#' @section Performance: As this method accepts spectra as its input, it
#'   computes irradiances before computing the ratios. If you need to compute
#'   both ratios and irradiances from several hundreds or thousands of spectra,
#'   computing the ratios from previously computed irradiances avoids their
#'   repeated computation. A less dramatic, but still important, increase in
#'   performance is available when computing in the same function call ratios
#'   that share the same denominator.
#'
#' @family photon and energy ratio functions
#'
qe_ratio <- function(spct,
                     w.band,
                     scale.factor,
                     wb.trim,
                     use.cached.mult,
                     use.hinges,
                     ...) UseMethod("qe_ratio")

#' @describeIn qe_ratio Default for generic function
#'
#' @export
#'
qe_ratio.default <- function(spct,
                             w.band,
                             scale.factor,
                             wb.trim,
                             use.cached.mult,
                             use.hinges,
                             ...) {
  warning("'qe_ratio' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn qe_ratio Method for \code{source_spct} objects
#'
#' @export
#'
qe_ratio.source_spct <-
  function(spct,
           w.band = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = NULL,
           naming = "short",
           name.tag = ifelse(naming != "none", "[q:e]", ""),
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(qe_ratio(spct = mspct,
                      w.band = w.band,
                      scale.factor = scale.factor,
                      wb.trim = wb.trim,
                      use.cached.mult = use.cached.mult,
                      use.hinges = use.hinges,
                      naming = naming,
                      name.tag = name.tag,
                      ...))
    }

    irrads <- two_irrads(spct = spct,
                         w.band.num = w.band,
                         w.band.denom = w.band,
                         unit.out.num = "photon",
                         unit.out.denom = "energy",
                         quantity = "total",
                         wb.trim = wb.trim,
                         use.cached.mult = use.cached.mult,
                         use.hinges = use.hinges,
                         naming = naming)

    q.irrad <- irrads[["irrad.num"]]
    e.irrad <- irrads[["irrad.denom"]]

    ratio <- q.irrad / e.irrad * scale.factor
    names(ratio) <- paste(names(q.irrad), name.tag, sep = "")
    attr(ratio, "time.unit") <- NULL
    attr(ratio, "radiation.unit") <- "q:e ratio"
    return(ratio)
  }

#' @describeIn qe_ratio Calculates photon:energy ratio from a
#'   \code{source_mspct} object.
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
qe_ratio.source_mspct <-
  function(spct, w.band=NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = NULL,
           naming = "short",
           name.tag = ifelse(naming != "none", "[q:e]", ""),
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
        spct,
        .fun = qe_ratio,
        w.band = w.band,
        wb.trim = wb.trim,
        scale.factor = scale.factor,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
        naming = naming,
        name.tag = name.tag,
        idx = idx,
        col.names = names(w.band),
        .parallel = .parallel,
        .paropts = .paropts
      )
    add_attr2tb(tb = z,
                mspct = spct,
                col.names = attr2tb,
                idx = idx)
  }

# eq_ratio() --------------------------------------------------------------

#' Energy:photon ratio
#'
#' This function returns the energy to mole of photons ratio for each waveband
#' and a light source spectrum.
#'
#' @param spct source_spct.
#' @param w.band waveband or list of waveband objects.
#' @param scale.factor numeric vector of length 1, or length equal to that of
#'   \code{w.band}. Numeric multiplier applied to returned values.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded.
#' @param use.cached.mult logical Flag telling whether multiplier values should
#'   be cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param name.tag character Used to tag the name of the returned values.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @details The ratio is based on one photon irradiance and one energy
#'   irradiance, both computed for the same waveband.
#'
#'   \deqn{\frac{I(s, wb)}{Q(s, wb)}}
#'
#'   The last two parameters control speed optimizations. The defaults should be
#'   suitable in most cases. If you will use repeatedly the same SWFs on many
#'   spectra measured at exactly the same wavelengths you may obtain some speed
#'   up by setting \code{use.cached.mult=TRUE}. However, be aware that you are
#'   responsible for ensuring that the wavelengths are the same in each call, as
#'   the only test done is for the length of the \code{w.length} vector.#'
#'   @return Computed values are ratios between energy irradiance and photon
#'   irradiance for a given waveband. A named \code{numeric} vector in the case
#'   of methods for individual spectra, with one value for each \code{waveband}
#'   passed to parameter \code{w.band}. A \code{data.frame} in the case of
#'   collections of spectra, containing one column for each \code{waveband}
#'   object, an index column with the names of the spectra, and optionally
#'   additional columns with metadata values retrieved from the attributes of
#'   the member spectra.
#'
#'   By default values are only integrated, but depending on the argument passed
#'   to parameter \code{quantity} they can be re-expressed as relative fractions
#'   or percentages. In the case of vector output, \code{names} attribute is set
#'   to the name of the corresponding waveband unless a named list is supplied
#'   in which case the names of the list members are used, with "[e:q]"
#'   prepended. Units [J mol-1].
#'
#' @return Computed values are ratios between energy irradiance and photon
#'   irradiance for a given waveband. A named \code{numeric} vector in the case
#'   of methods for individual spectra, with one value for each \code{waveband}
#'   passed to parameter \code{w.band}. A \code{data.frame} in the case of
#'   multiple spectra, containing one column with ratios for each \code{waveband}
#'   object, an index column with the names of the spectra, and optionally
#'   additional columns with metadata values retrieved from the attributes of
#'   the member spectra.
#'
#'   By default values are only integrated, but depending on the argument passed
#'   to parameter \code{quantity} they are expressed as relative fractions
#'   or percentages. In the case of vector output, \code{names} attribute is set
#'   to the name of the corresponding waveband unless a named list is supplied
#'   in which case the names of the list members are used, with "[e:q]" prepended.
#'   Units [mol J-1].
#'
#' @export
#'
#' @examples
#' eq_ratio(sun.spct,
#'          waveband(c(400,700), wb.name = "White")) # J mol-1
#' eq_ratio(sun.spct,
#'          waveband(c(400,700), wb.name = "White"),
#'          scale.factor = 1e-6) # J umol-1
#'
#' @section Performance: As this method accepts spectra as its input, it
#'   computes irradiances before computing the ratios. If you need to compute
#'   both ratios and irradiances from several hundreds or thousands of spectra,
#'   computing the ratios from previously computed irradiances avoids their
#'   repeated computation. A less dramatic, but still important, increase in
#'   performance is available when computing in the same function call ratios
#'   that share the same denominator.
#'
#' @family photon and energy ratio functions
#'
eq_ratio <- function(spct,
                     w.band,
                     scale.factor,
                     wb.trim,
                     use.cached.mult,
                     use.hinges,
                     ...) UseMethod("eq_ratio")

#' @describeIn eq_ratio Default for generic function
#'
#' @export
#'
eq_ratio.default <- function(spct,
                             w.band,
                             scale.factor,
                             wb.trim,
                             use.cached.mult,
                             use.hinges,
                             ...) {
  warning("'eq_ratio' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn eq_ratio Method for \code{source_spct} objects
#'
#' @export
#'
eq_ratio.source_spct <-
  function(spct,
           w.band = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges  = NULL,
           naming = "short",
           name.tag = ifelse(naming != "none", "[e:q]", ""),
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(eq_ratio(spct = mspct,
                      w.band = w.band,
                      scale.factor = scale.factor,
                      wb.trim = wb.trim,
                      use.cached.mult = use.cached.mult,
                      use.hinges = use.hinges,
                      naming = naming,
                      name.tag = name.tag,
                      ...))
    }

    ratio <- scale.factor /
      qe_ratio(spct = spct, w.band = w.band, wb.trim = wb.trim,
               use.cached.mult = use.cached.mult, use.hinges = use.hinges)
    names(ratio) <- gsub("q:e", "e:q", names(ratio), fixed = TRUE )
    attr(ratio, "time.unit") <- NULL
    attr(ratio, "radiation.unit") <- "e:q ratio"
    return(ratio)
  }

#' @describeIn eq_ratio Calculates energy:photon from a \code{source_mspct}
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
eq_ratio.source_mspct <-
  function(spct, w.band = NULL,
           scale.factor = 1,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = NULL,
           naming = "short",
           name.tag = ifelse(naming != "none", "[e:q]", ""),
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
        .fun = eq_ratio,
        w.band = w.band,
        wb.trim = wb.trim,
        scale.factor = scale.factor,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
        naming = naming,
        name.tag = name.tag,
        idx = idx,
        col.names = names(w.band),
        .parallel = .parallel,
        .paropts = .paropts
      )
    add_attr2tb(tb = z,
                mspct = spct,
                col.names = attr2tb,
                idx = idx)
  }

# internal utility function --------------------------------------------------

#' Compute two irrads for ratio, fraction or normalised difference
#'
#' Internal function that computes the two irradiances needed to compute
#' various waveband ratios and fractions.
#'
#' @details See \code{\link{irrad}} for details on the irradiance calculations.
#'
#' @param spct an object of class "source_spct".
#' @param w.band.num,w.band.denom waveband objects or lists of waveband objects
#'   used to compute the numerator(s) and denominator(s) of the ratio(s). The
#'   waveband(s) determine the region(s) of the spectrum that are summarized.
#' @param unit.out.num,unit.out.denom character Allowed values "energy", and
#'   "photon", or its alias "quantum".
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
# repetition. It is used to define ratios and fractions
#
two_irrads <- function(spct,
                       w.band.num,
                       w.band.denom,
                       unit.out.num,
                       unit.out.denom,
                       quantity,
                       wb.trim,
                       use.cached.mult,
                       use.hinges,
                       naming) {

  stopifnot("Unsupported argument passed to 'quantity'" =
              quantity %in% c("total", "average", "mean"))

  irrad.num <- irrad(spct,
                     w.band = w.band.num,
                     unit.out = unit.out.num,
                     quantity = quantity,
                     scale.factor = 1,
                     wb.trim = wb.trim,
                     use.cached.mult = use.cached.mult,
                     use.hinges = use.hinges,
                     allow.scaled = TRUE,
                     naming = naming)

  irrad.denom <- irrad(spct,
                       w.band = w.band.denom,
                       unit.out = unit.out.denom,
                       quantity = quantity,
                       scale.factor = 1,
                       wb.trim = wb.trim,
                       use.cached.mult = use.cached.mult,
                       use.hinges = use.hinges,
                       allow.scaled = TRUE,
                       naming = naming)

  list(irrad.num = irrad.num, irrad.denom = irrad.denom)
}
