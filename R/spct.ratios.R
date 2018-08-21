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
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param ... other arguments (possibly ignored)
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector of adimensional values giving a photon ratio between integrated
#'   photon irradiances for pairs of wavebands, with name attribute set to
#'   the name of the wavebands unless a named list of wavebands is supplied in
#'   which case the names of the list elements are used, with "(q:q)" appended.
#'   A \code{data.frame} in the case of collections of spectra, containing one
#'   column for each ratio definition, an index column with the names of the
#'   spectra, and optionally additional columns with metadata values retrieved
#'   from the attributes of the member spectra.
#'
#'   Ratio definitions are "assembled" from the arguments passed to
#'   \code{w.band.num} and \code{w.band.denom}. If both arguments are of equal
#'   length, then the wavebands are paired to obtain as many ratios as the
#'   number of wavebands in each list. Recycling for wavebands takes place when
#'   the number of denominator and numerator wavebands differ.
#'
#' @export
#' @examples
#' q_ratio(sun.spct, new_waveband(400,500), new_waveband(400,700))
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
q_ratio <- function(spct, w.band.num, w.band.denom, wb.trim,
                  use.cached.mult, use.hinges, ...) UseMethod("q_ratio")

#' @describeIn q_ratio Default for generic function
#'
#' @export
#'
q_ratio.default <- function(spct, w.band.num, w.band.denom, wb.trim,
                            use.cached.mult, use.hinges, ...) {
  warning("'q_ratio' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn q_ratio Method for \code{source_spct} objects
#'
#' @export
#'
q_ratio.source_spct <-
  function(spct,
           w.band.num = NULL, w.band.denom = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges"), ... ) {
    q.irrad.num <- irrad_spct(spct, w.band = w.band.num,
                              unit.out = "photon", quantity = "total",
                              wb.trim = wb.trim,
                              use.cached.mult = use.cached.mult,
                              use.hinges = use.hinges,
                              allow.scaled = TRUE)
    q.irrad.denom <- irrad_spct(spct, w.band = w.band.denom,
                                unit.out = "photon", quantity = "total",
                                wb.trim = wb.trim,
                                use.cached.mult = use.cached.mult,
                                use.hinges = use.hinges,
                                allow.scaled = TRUE)
    ratio <- q.irrad.num / q.irrad.denom
    names(ratio) <- paste(names(q.irrad.num), ":", names(q.irrad.denom), "(q:q)", sep = "")
    attr(ratio, "time.unit") <- NULL
    attr(ratio, "radiation.unit") <- "q:q ratio"
    return(ratio)
  }

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
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical Flag telling whether multiplier values should be
#'   cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @return In the case of methods for individual spectra, a \code{numeric}
#'   vector of adimensional values giving a energy ratio between
#'   integrated energy irradiances for pairs of wavebands, with name attribute
#'   set to the name of the wavebands unless a named list of wavebands is
#'   supplied in which case the names of the list elements are used, with
#'   "(e:e)" appended. A \code{data.frame} in the case of collections of
#'   spectra, containing one column for each ratio definition, an index column
#'   with the names of the spectra, and optionally additional columns with
#'   metadata values retrieved from the attributes of the member spectra.
#'
#'   Ratio definitions are "assembled" from the arguments passed to
#'   \code{w.band.num} and \code{w.band.denom}. If both arguments are of equal
#'   length, then the wavebands are paired to obtain as many ratios as the
#'   number of wavebands in each list. Recycling for wavebands takes place when
#'   the number of denominator and numerator wavebands differ.
#'
#' @export e_ratio
#' @examples
#' e_ratio(sun.spct, new_waveband(400,500), new_waveband(400,700))
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
e_ratio <- function(spct, w.band.num, w.band.denom, wb.trim,
                    use.cached.mult, use.hinges, ...) UseMethod("e_ratio")

#' @describeIn e_ratio Default for generic function
#'
#' @export
#'
e_ratio.default <- function(spct, w.band.num, w.band.denom, wb.trim,
                            use.cached.mult, use.hinges, ...) {
  warning("'e_ratio' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn e_ratio Method for \code{source_spct} objects
#'
#' @export
#'
e_ratio.source_spct <-
  function(spct,
           w.band.num = NULL, w.band.denom = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges=getOption("photobiology.use.hinges"), ...) {
    e.irrad.num <- irrad_spct(spct, w.band = w.band.num, unit.out = "energy", quantity = "total",
                              wb.trim = wb.trim,
                              use.cached.mult = use.cached.mult, use.hinges = use.hinges,
                              allow.scaled=TRUE)
    e.irrad.denom <- irrad_spct(spct, w.band = w.band.denom, unit.out = "energy", quantity = "total",
                                wb.trim = wb.trim,
                                use.cached.mult = use.cached.mult, use.hinges = use.hinges,
                                allow.scaled = TRUE)
    ratio <- e.irrad.num / e.irrad.denom
    names(ratio) <- paste(names(e.irrad.num), ":", names(e.irrad.denom), "(e:e)", sep="")
    attr(ratio, "time.unit") <- NULL
    attr(ratio, "radiation.unit") <- "e:e ratio"
    return(ratio)
  }

#' Photon:energy ratio
#'
#' This function returns the photon to energy ratio for each waveband of a light
#' source spectrum.
#'
#' @param spct source_spct.
#' @param w.band waveband or list of waveband objects.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded.
#' @param use.cached.mult logical Flag telling whether multiplier values should be
#'   cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param ... other arguments (possibly used by derived methods).
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
#'   in which case the names of the list members are used, with "q:e" prepended.
#'   Units [mol J-1].
#'
#' @export
#' @examples
#' qe_ratio(sun.spct, new_waveband(400,700))
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @family photon and energy ratio functions
#'
qe_ratio <- function(spct, w.band, wb.trim,
                     use.cached.mult, use.hinges, ...) UseMethod("qe_ratio")

#' @describeIn qe_ratio Default for generic function
#'
#' @export
#'
qe_ratio.default <- function(spct, w.band, wb.trim,
                             use.cached.mult, use.hinges, ...) {
  warning("'qe_ratio' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn qe_ratio Method for \code{source_spct} objects
#'
#' @export
#'
qe_ratio.source_spct <-
  function(spct, w.band = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges"), ...) {
    q.irrad <- irrad_spct(spct, w.band=w.band, unit.out = "photon",
                          quantity ="total",
                          wb.trim = wb.trim,
                          use.cached.mult = use.cached.mult,
                          use.hinges = use.hinges,
                          allow.scaled = TRUE)
    e.irrad <- irrad_spct(spct, w.band=w.band, unit.out = "energy",
                          quantity = "total",
                          wb.trim = wb.trim,
                          use.cached.mult = use.cached.mult,
                          use.hinges = use.hinges,
                          allow.scaled = TRUE)
    ratio <- q.irrad / e.irrad
    names(ratio) <- paste("q:e(", names(q.irrad), ")", sep = "")
    attr(ratio, "time.unit") <- NULL
    attr(ratio, "radiation.unit") <- "q:e ratio"
    return(ratio)
  }

#' Energy:photon ratio
#'
#' This function returns the energy to mole of photons ratio for each waveband and a
#' light source spectrum.
#'
#' @param spct source_spct.
#' @param w.band waveband or list of waveband objects.
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded.
#' @param use.cached.mult logical Flag telling whether multiplier values should
#'   be cached between calls.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @return Computed values are ratios between energy irradiance and photon
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
#'   in which case the names of the list members are used, with "e:q" prepended.
#'   Units [J mol-1].
#'
#' @export
#' @examples
#' eq_ratio(sun.spct, new_waveband(400,700))
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in most cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @family photon and energy ratio functions
#'
eq_ratio <- function(spct, w.band, wb.trim,
                     use.cached.mult, use.hinges, ...) UseMethod("eq_ratio")

#' @describeIn eq_ratio Default for generic function
#'
#' @export
#'
eq_ratio.default <- function(spct, w.band, wb.trim,
                             use.cached.mult, use.hinges, ...) {
  warning("'eq_ratio' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn eq_ratio Method for \code{source_spct} objects
#'
#' @export
#'
eq_ratio.source_spct <-
  function(spct, w.band=NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult =FALSE,
           use.hinges  = getOption("photobiology.use.hinges"), ...) {
    ratio <- 1 / qe_ratio(spct = spct, w.band = w.band, wb.trim = wb.trim,
                          use.cached.mult = use.cached.mult, use.hinges = use.hinges)
    names(ratio) <- gsub("q:e", "e:q", names(ratio), fixed = TRUE )
    attr(ratio, "time.unit") <- NULL
    attr(ratio, "radiation.unit") <- "e:q ratio"
    return(ratio)
  }

# source_mspct methods ----------------------------------------------------

#' @describeIn q_ratio Calculates photon:photon from a \code{source_mspct}
#'   object.
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax for \code{attr2tb} passed as is to formal parameter \code{col.names}.
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
           w.band.num = NULL, w.band.denom = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges"),
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {
    z <-
      msdply(
        mspct = spct,
        .fun = q_ratio,
        w.band.num = w.band.num,
        w.band.denom = w.band.denom,
        wb.trim = wb.trim,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
        idx = idx,
        .parallel = .parallel,
        .paropts = .paropts
      )
    add_attr2tb(tb = z,
                mspct = spct,
                col.names = attr2tb,
                idx = idx)
  }

#' @describeIn e_ratio Calculates energy:energy ratio from a \code{source_mspct}
#'   object.
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax for \code{attr2tb} passed as is to formal parameter \code{col.names}.
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
e_ratio.source_mspct <-
  function(spct,
           w.band.num = NULL, w.band.denom = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges"),
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {
    z <-
      msdply(
        mspct = spct,
        .fun = e_ratio,
        w.band.num = w.band.num,
        w.band.denom = w.band.denom,
        wb.trim = wb.trim,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
        idx = idx,
        .parallel = .parallel,
        .paropts = .paropts
      )
    add_attr2tb(tb = z,
                mspct = spct,
                col.names = attr2tb,
                idx = idx)
  }

#' @describeIn eq_ratio Calculates energy:photon from a \code{source_mspct}
#'   object.
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax for \code{attr2tb} passed as is to formal parameter \code{col.names}.
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
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges"),
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {
    z <-
      msdply(
        mspct = spct,
        .fun = eq_ratio,
        w.band = w.band,
        wb.trim = wb.trim,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
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

#' @describeIn qe_ratio Calculates photon:energy ratio from a
#'   \code{source_mspct} object.
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax for \code{attr2tb} passed as is to formal parameter \code{col.names}.
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
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges=getOption("photobiology.use.hinges"),
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {
    z <-
      msdply(
        spct,
        .fun = qe_ratio,
        w.band = w.band,
        wb.trim = wb.trim,
        use.cached.mult = use.cached.mult,
        use.hinges = use.hinges,
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
