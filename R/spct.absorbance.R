#' Absorbance
#'
#' Function to calculate the mean, total, or other summary of absorbance for
#' spectral data stored in a \code{filter_spct} or in an \code{object_spct}.
#'
#' @param spct an R object.
#' @param w.band waveband or list of waveband objects or a numeric vector of
#'   length two. The waveband(s) determine the region(s) of the spectrum that
#'   are summarized. If a numeric range is supplied a waveband object is
#'   constructed on the fly from it.
#' @param quantity character string One of "average" or "mean", "total",
#'   "contribution", "contribution.pc", "relative" or "relative.pc".
#' @param wb.trim logical Flag indicating if wavebands crossing spectral data
#'   boundaries are trimmed or ignored.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param ... other arguments (possibly used by derived methods).
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
#' @note The \code{use.hinges} parameter controls speed optimization. The
#'   defaults should be suitable in most cases. Only the range of wavelengths in
#'   the wavebands is used and all BSWFs are ignored.
#'
#' @examples
#' absorbance(polyester.spct, new_waveband(400,700))
#' absorbance(yellow_gel.spct, new_waveband(400,700))
#' absorbance(yellow_gel.spct, split_bands(c(400,700), length.out = 3))
#' absorbance(yellow_gel.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "average")
#' absorbance(yellow_gel.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "total")
#' absorbance(yellow_gel.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "relative")
#' absorbance(yellow_gel.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "relative.pc")
#' absorbance(yellow_gel.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "contribution")
#' absorbance(yellow_gel.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "contribution.pc")
#'
#' @export
#'
absorbance <- function(spct, w.band, quantity, wb.trim, use.hinges, ...) UseMethod("absorbance")

#' @describeIn absorbance Default for generic function
#'
#' @export
#'
absorbance.default <- function(spct, w.band, quantity, wb.trim, use.hinges, ...) {
  warning("'absorbance' is not defined for objects of class ", class(spct)[1])
  return(NA_real_)
}

#' @describeIn absorbance Specialization for filter spectra
#'
#' @export
#'
absorbance.filter_spct <-
  function(spct,
           w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = NULL,
           naming = "default",
           ...) {
    absorbance_spct(spct = spct,
                    w.band = w.band,
                    quantity = quantity,
                    wb.trim = wb.trim,
                    use.hinges = use.hinges,
                    naming = naming)
  }

#' @describeIn absorbance Specialization for object spectra
#'
#' @export
#'
absorbance.object_spct <-
  function(spct,
           w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = NULL,
           naming = "default",
           ...) {
    spct <- as.filter_spct(spct)
    absorbance_spct(spct,
                    w.band = w.band,
                    quantity = quantity,
                    wb.trim = wb.trim,
                    use.hinges = use.hinges,
                    naming = naming)
  }

#' Calculate absorbance from spectral absorbance.
#'
#' This function returns the mean absorbance for a given
#' waveband of a absorbance spectrum.
#'
#' @param spct filter_spct
#' @param w.band waveband or list of waveband objects or a numeric vector of
#'   length two. The waveband(s) determine the region(s) of the spectrum that
#'   are summarized. If a numeric range is supplied a waveband object is
#'   constructed on the fly from it.
#' @param quantity character string One of "average" or "mean", "total",
#'   "contribution", "contribution.pc", "relative" or "relative.pc"
#' @param wb.trim logical Flag if wavebands crossing spectral data boundaries
#'   are trimmed or ignored
#' @param use.hinges logical Flag indicating whether to use hinges to reduce
#'   interpolation errors
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @keywords internal
#'
absorbance_spct <-
  function(spct,
           w.band,
           quantity,
           wb.trim,
           use.hinges,
           naming,
           ...) {

    # we look for multiple spectra in long form
    num.spectra <- getMultipleWl(spct)
    if (num.spectra != 1) {
      message("Object contains ", num.spectra, " spectra in long form")
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(absorbance(spct = mspct,
                        w.band = w.band,
                        quantity = quantity,
                        wb.trim = wb.trim,
                        use.hinges = use.hinges,
                        naming,
                        ...))
    }

    summary.name <-
      switch(quantity,
             total = "A",
             average = "A(wl)",
             mean = "A(wl)",
             contribution = "A/Atot",
             contribution.pc = "A/Atot[%]",
             relative = "A/Asum",
             relative.pc = "A/Asum[%]",
             stop("Unrecognized 'quantity' : \"", quantity, "\"")
      )

    if (is_normalized(spct)) {
      warning("The spectral data has been normalized,",
              "making impossible to calculate absorbance")
      return(NA_real_)
    }
    if (is_scaled(spct)) {
      warning("Summary calculated from rescaled data")
    }
    spct <- T2A(spct, action = "replace", byref = FALSE)
    spct <- spct[ , c("w.length", "A")]

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

    # We iterate through the list of wavebands collecting the absorbances,
    # and waveband names.
    absorbance <- numeric(length(w.band))
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

      # we calculate the average absorbance.
      absorbance[i] <- integrate_spct(trim_spct(spct, wb, use.hinges = FALSE))
    }

    if (quantity %in% c("contribution", "contribution.pc")) {
      total <- absorbance_spct(spct,
                               w.band = NULL,
                               quantity = "total",
                               wb.trim = wb.trim,
                               use.hinges = use.hinges,
                               naming = naming)
      absorbance <- absorbance / total
      if (quantity == "contribution.pc") {
        absorbance <- absorbance * 1e2
      }
    } else if (quantity %in% c("relative", "relative.pc")) {
      total <- sum(absorbance)
      absorbance <- absorbance / total
      if (quantity == "relative.pc") {
        absorbance <- absorbance * 1e2
      }
    } else if (quantity %in% c("average", "mean")) {
      absorbance <- absorbance / sapply(w.band, wl_expanse)
    }

    if (length(absorbance) == 0) {
      absorbance <- NA_real_
      names(absorbance) <- "out of range"
    } else if (naming %in% c("long", "default")) {
      names(absorbance) <- paste(summary.name, wb.name, sep = "_")
    } else if (naming == "short") {
      names(absorbance) <- wb.name
    } else if (naming != "none") {
      warning("Argument to 'naming' unrecognized, assuming \"none\".")
    }

    attr(absorbance, "Tfr.type") <- getTfrType(spct)
    attr(absorbance, "radiation.unit") <- paste("absorbance", quantity)

    absorbance
  }

# filter_mspct methods -----------------------------------------------

#' @describeIn absorbance Calculates absorbance from a \code{filter_mspct}
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax for \code{attr2tb} passed as is to formal parameter \code{col.names}.
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#'
#' @export
#'
absorbance.filter_mspct <-
  function(spct, w.band=NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = NULL,
           naming = "default",
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {
    z <-
      msdply(
        mspct = spct,
        .fun = absorbance,
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

#' @describeIn absorbance Calculates absorbance from a \code{object_mspct}
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
absorbance.object_mspct <-
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
    z <-
      msdply(
        mspct = spct,
        .fun = absorbance,
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
