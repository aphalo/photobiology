#' Reflectance
#'
#' Function to calculate the mean, total, or other summary of reflectance for
#' spectral data stored in a \code{reflector_spct} or in an \code{object_spct}.
#'
#' @param spct an R object
#' @param w.band waveband or list of waveband objects or a numeric vector of
#'   length two. The waveband(s) determine the region(s) of the spectrum that
#'   are summarized. If a numeric range is supplied a waveband object is
#'   constructed on the fly from it.
#' @param quantity character string One of \code{"average"} or \code{"mean"},
#'   \code{"total"}, \code{"contribution"}, \code{"contribution.pc"},
#'   \code{"relative"} or \code{"relative.pc"}.
#' @param wb.trim logical if \code{TRUE} wavebands crossing spectral data
#'   boundaries are trimmed, if \code{FALSE}, they are discarded.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param naming character one of \code{"long"}, \code{"default"},
#'   \code{"short"} or \code{"none"}. Used to select the type of names to assign
#'   to returned value.
#' @param ... other arguments
#'
#' @note The \code{use.hinges} parameter controls speed optimization. The
#'   defaults should be suitable in most cases. Only the range of wavelengths
#'   in the wavebands is used and all BSWFs are ignored.
#'
#' @return A named \code{numeric} vector in the case of methods for individual
#'   spectra, with one value for each \code{waveband} passed to parameter
#'   \code{w.band}. A \code{data.frame} in the case of collections of spectra,
#'   containing one column for each \code{waveband} object, an index column with
#'   the names of the spectra, and optionally additional columns with metadata
#'   values retrieved from the attributes of the member spectra.
#'
#'   By default values are only integrated, but depending on the argument passed
#'   to parameter \code{quantity} they can be re-expressed as relative fractions
#'   or percentages. In the case of vector output, \code{names} attribute is set
#'   to the name of the corresponding waveband unless a named list is supplied
#'   in which case the names of the list members are used.
#'
#' @examples
#' reflectance(black_body.spct, waveband(c(400,700)))
#' reflectance(white_body.spct, waveband(c(400,700)))
#'
#' @export
#'
reflectance <- function(spct, w.band, quantity, wb.trim, use.hinges, ...) UseMethod("reflectance")

#' @describeIn reflectance Default for generic function
#'
#' @export
#'
reflectance.default <- function(spct, w.band, quantity, wb.trim, use.hinges, ...) {
  warning("'reflectance' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn reflectance Specialization for reflector_spct
#'
#' @export
#'
reflectance.reflector_spct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges = NULL,
           naming = "default",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(reflectance(spct = mspct,
                         w.band = w.band,
                         quantity = quantity,
                         wb.trim = wb.trim,
                         use.hinges = use.hinges,
                         naming = naming,
                         ...))
    }

    reflectance_spct(spct = spct, w.band = w.band,
                     quantity = quantity,
                     wb.trim = wb.trim,
                     use.hinges = use.hinges,
                     naming = naming)
  }

#' @describeIn reflectance Specialization for object_spct
#'
#' @export
#'
reflectance.object_spct <-
  function(spct,
           w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges = NULL,
           naming = "default",
           ... ) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(reflectance(spct = mspct,
                         w.band = w.band,
                         quantity = quantity,
                         wb.trim = wb.trim,
                         use.hinges = use.hinges,
                         naming = naming,
                         ...))
    }

    reflectance_spct(spct = spct, w.band = w.band,
                     quantity = quantity,
                     wb.trim = wb.trim,
                     use.hinges = use.hinges,
                     naming = naming)
  }

#' Calculate reflectance from spectral reflectance
#'
#' This function returns the mean reflectance for a given waveband and a
#' reflectance spectrum.
#'
#' @param spct an object of class generic_spct"
#' @param w.band waveband or list of waveband objects or a numeric vector of
#'   length two. The waveband(s) determine the region(s) of the spectrum that
#'   are summarized. If a numeric range is supplied a waveband object is
#'   constructed on the fly from it.
#' @param quantity character string One of "total", "average" or "mean",
#'   "contribution", "contribution.pc", "relative" or "relative.pc"
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @return A single numeric value expressed as a fraction of one
#' @keywords internal
#'
reflectance_spct <-
  function(spct,
           w.band,
           quantity,
           wb.trim,
           use.hinges,
           naming,
           ...){

    summary.name <-
      switch(quantity,
             total = "Rfr",
             average = "Rfr(wl)",
             mean = "Rfr(wl)",
             contribution = "Rfr/Rfrtot",
             contribution.pc = "Rfr/Rfrtot[%]",
             relative = "Rfr/Rfrsum",
             relative.pc = "Rfr/Rfrsum[%]",
             stop("Unrecognized 'quantity' : \"", quantity, "\"")
      )

    if (is_normalized(spct)) {
      warning("The spectral data has been normalized,",
              "making impossible to calculate reflectance")
      return(NA_real_)
    }
    if (is_scaled(spct)) {
      warning("Reflectance calculated from rescaled data")
    }

    if (is.object_spct(spct)) {
      spct <- as.reflector_spct(spct)
    }
    spct <- spct[ , c("w.length", "Rfr")]

    if (length(w.band) == 0) {
      # whole range of spectrum
      w.band <- waveband(spct)
    }
    if (is.numeric(w.band)) {
      w.band <- waveband(w.band)
    }
    if (is.waveband(w.band)) {
      # if the argument is a single w.band, we enclose it in a list
      # so that it can be handled below as a normal case.
      w.band <- list(w.band)
    }
    # we trim the wavebands so that they are within the range of spct
    w.band <- trim_waveband(w.band = w.band, range = spct, trim = wb.trim)
    # if the elements of the list are named we collect them
    wb.number <- length(w.band) # number of wavebands in list
    wb.name <- names(w.band) # their names in the list
    # if no names returned, we fill the vector with "".
    if (is.null(wb.name)) {
      wb.name <- character(wb.number)
    }

    # hinges
    if (is.null(use.hinges)) {
      use.hinges <- auto_hinges(spct[["w.length"]])
    }
    # we collect all hinges and insert them in one go
    if (use.hinges) {
      all.hinges <- NULL
      for (wb in w.band) {
        all.hinges <- c(all.hinges, wb[["hinges"]])
      }
      if (!is.null(all.hinges)) {
        spct <- insert_spct_hinges(spct, all.hinges)
      }
    }

    # We iterate through the list of wavebands collecting the transmittances,
    # and waveband names.
    reflectance <- numeric(length(w.band))
    i <- 0L
    for (wb in w.band) {
      i <- i + 1L
      # weighting functions are not meaningful
      if (is_effective(wb)) {
        warning("Using wavelength range from a weighted waveband object.")
        wb <- waveband(wl_range(wb))
      }
      # we get names from wb if needed
      if (wb.name[i] == "") {
        if (naming == "short") {
          wb.name[i] <- labels(wb)[["label"]] # short name
        } else {
          wb.name[i] <- labels(wb)[["name"]] # full name
        }
      }

      # we calculate the average reflectance
      reflectance[i] <- integrate_spct(trim_spct(spct, wb, use.hinges = FALSE))
    }

    if (quantity %in% c("contribution", "contribution.pc")) {
      total <- reflectance_spct(spct,
                                w.band = NULL,
                                wb.trim = wb.trim,
                                quantity = "total",
                                use.hinges = use.hinges,
                                naming = naming)
      reflectance <- reflectance / total
      if (quantity == "contribution.pc") {
        reflectance <- reflectance * 1e2
      }
    } else if (quantity %in% c("relative", "relative.pc")) {
      total <- sum(reflectance)
      reflectance <- reflectance / total
      if (quantity == "relative.pc") {
        reflectance <- reflectance * 1e2
      }
    } else if (quantity %in% c("average", "mean")) {
      reflectance <- reflectance / sapply(w.band, wl_expanse)
    } else if (quantity == "total") {
    } else if (quantity != "total") {
      warning("'quantity '", quantity, "' is invalid, returning 'total' instead")
      quantity <- "total"
    }

    if (length(reflectance) == 0) {
      reflectance <- NA_real_
      names(reflectance) <- "out of range"
    } else if (naming %in% c("long", "default")) {
      names(reflectance) <- paste(summary.name, wb.name, sep = "_")
    } else if (naming == "short") {
      names(reflectance) <- wb.name
    } else if (naming != "none") {
      warning("Argument to 'naming' unrecognized, assuming \"none\".")
    }

    attr(reflectance, "Rfr.type") <- getRfrType(spct)
    attr(reflectance, "radiation.unit") <- paste("reflectance", quantity)
    reflectance
  }

# reflector_mspct methods -----------------------------------------------

#' @describeIn reflectance Calculates reflectance from a \code{reflector_mspct}
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
reflectance.reflector_mspct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = NULL,
           naming = "default",
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {

    spct <- subset2mspct(spct) # expand long form spectra within collection

    z <-
      msdply(
        mspct = spct,
        .fun = reflectance,
        w.band = w.band,
        quantity = quantity,
        wb.trim = wb.trim,
        use.hinges = use.hinges,
        naming = naming,
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

# object_mspct methods -----------------------------------------------

#' @describeIn reflectance Calculates reflectance from a \code{object_mspct}
#'
#' @export
#'
reflectance.object_mspct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = NULL,
           naming = "default",
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {

    spct <- subset2mspct(spct) # expand long form spectra within collection

    z <-
      msdply(
        mspct = spct,
        .fun = reflectance,
        w.band = w.band,
        quantity = quantity,
        wb.trim = wb.trim,
        use.hinges = use.hinges,
        naming = naming,
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
