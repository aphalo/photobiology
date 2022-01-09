#' Absorptance
#'
#' Function to calculate the mean, total, or other summary of absorptance for
#' spectral data stored in a \code{filter_spct} or in an \code{object_spct}.
#' Absorptance is a different quantity than absorbance.
#'
#' @param spct an R object.
#' @param w.band waveband or list of waveband objects or a numeric vector of
#'   length two. The waveband(s) determine the region(s) of the spectrum that
#'   are summarized. If a numeric range is supplied a waveband object is
#'   constructed on the fly from it.
#' @param quantity character string One of "average" or "mean", "total",
#'   "contribution", "contribution.pc", "relative" or "relative.pc".
#' @param wb.trim logical Flag, if TRUE wavebands crossing spectral data
#'   boundaries are trimmed and otherwise ignored.
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   spectral data before integration so as to reduce interpolation errors at
#'   the boundaries of the wavebands.
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param ... other arguments (possibly used by derived methods).
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
#' absorptance(black_body.spct, new_waveband(400,500))
#' absorptance(white_body.spct, new_waveband(300,400))
#' absorptance(black_body.spct, split_bands(c(400,700), length.out = 3))
#' absorptance(black_body.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "average")
#' absorptance(black_body.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "total")
#' absorptance(black_body.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "relative")
#' absorptance(black_body.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "relative.pc")
#' absorptance(black_body.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "contribution")
#' absorptance(black_body.spct, split_bands(c(400,700), length.out = 3),
#'         quantity = "contribution.pc")
#'
#' @export
#'
absorptance <- function(spct, w.band, quantity, wb.trim, use.hinges, ...) UseMethod("absorptance")

#' @describeIn absorptance Default for generic function
#'
#' @export
#'
absorptance.default <- function(spct, w.band, quantity, wb.trim, use.hinges, ...) {
  warning("'absorptance' is not defined for objects of class ", class(spct)[1])
  return(NA_real_)
}

#' @describeIn absorptance Specialization for filter spectra
#'
#' @export
#'
absorptance.filter_spct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = NULL,
           naming = "default",
           ... ) {
    if (getTfrType(spct) != "internal") {
      warning("Internal absorptance cannot be calculated from total transmittance alone")
      return(NA)
    } else {
      absorptance_spct(spct = spct,
                       w.band = w.band,
                       quantity = quantity,
                       wb.trim = wb.trim,
                       use.hinges = use.hinges,
                       naming = naming)
    }
  }

#' @describeIn absorptance Specialization for object spectra
#'
#' @export
#'
absorptance.object_spct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = NULL,
           naming = "default",
           ...)  {
    absorptance_spct(spct = spct,
                     w.band = w.band,
                     quantity = quantity,
                     wb.trim = wb.trim,
                     use.hinges = use.hinges,
                     naming = naming)
  }

#' Calculate absorptance from spectral absorptance.
#'
#' This function returns the summary absorptance for a given waveband of a
#' \code{object_spct} object
#'
#' @param spct object_spct
#' @param w.band waveband or list of waveband objects or a numeric vector of
#'   length two. The waveband(s) determine the region(s) of the spectrum that
#'   are summarized. If a numeric range is supplied a waveband object is
#'   constructed on the fly from it.
#' @param quantity character string One of "average" or "mean", "contribution",
#'   "contribution.pc", "relative" or "relative.pc"
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#' @param naming character one of "long", "default", "short" or "none". Used to
#'   select the type of names to assign to returned value.
#' @param ... other arguments (possibly used by derived methods).
#'
#' @keywords internal
#'
absorptance_spct <-
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
      return(absorptance(spct = mspct,
                         w.band = w.band,
                         quantity = quantity,
                         wb.trim = wb.trim,
                         use.hinges = use.hinges,
                         naming,
                         ...))
    }

    summary.name <-
      switch(quantity,
             total = "Afr",
             average = "Afr(wl)",
             mean = "Afr(wl)",
             contribution = "Afr/Afrtot",
             contribution.pc = "Afr/Afrtot[%]",
             relative = "Afr/Afrsum",
             relative.pc = "Afr/Afrsum[%]",
             stop("Unrecognized 'quantity' : \"", quantity, "\"")
      )

    if (is_normalized(spct)) {
      warning("The spectral data has been normalized,",
              "making impossible to calculate absorptance")
      return(NA_real_)
    }
    if (is_scaled(spct)) {
      warning("Summary calculated from rescaled data")
    }

    # we calculate absorptance
    Tfr.type <- getTfrType(spct)
    Rfr.type <- getRfrType(spct)
    if (is.filter_spct(spct) && Tfr.type == "internal") {
      Rfr.type <- "unknown" # otherwise NA would require special handling
      A2T(spct, action = "add", byref = TRUE)
      temp.spct <- tibble::tibble(w.length = spct[["w.length"]],
                                     Afr = 1 - spct[["Tfr"]])
    } else if (Tfr.type == "total" && Rfr.type == "total") {
      temp.spct <- tibble::tibble(w.length = spct[["w.length"]],
                               Afr = 1 - spct[["Tfr"]] - spct[["Rfr"]])
     } else if (Tfr.type == "internal" && Rfr.type == "total") {
      temp.spct <- tibble::tibble(w.length = spct[["w.length"]],
                                   Afr = (1 - spct[["Tfr"]]) * (1 - spct[["Rfr"]]))
    } else if (Tfr.type == "unknown" || Rfr.type == "unknown") {
      warning("'unknown' Tfr.type or Rfr.type, skipping absorptance calculation")
      absorptance <- NA
      attr(absorptance, "radiation.unit") <- paste("absorptance", quantity)
      return(absorptance)
    } else if (Rfr.type == "specular") {
      warning("'specular' Rfr.type, skipping absorptance calculation")
      absorptance <- NA
      attr(absorptance, "radiation.unit") <- paste("absorptance", quantity)
      return(absorptance)
    } else {
      stop("Failed assertion with Tfr.type: ", Tfr.type, "and Rfr.type: ", Rfr.type)
    }
    temp.spct <- setGenericSpct(temp.spct)

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
      use.hinges <- auto_hinges(temp.spct[["w.length"]])
    }
    # we collect all hinges and insert them in one go
    if (use.hinges) {
      all.hinges <- NULL
      for (wb in w.band) {
        all.hinges <- c(all.hinges, wb[["hinges"]])
      }
      if (!is.null(all.hinges)) {
        temp.spct <- insert_spct_hinges(temp.spct, all.hinges)
      }
    }

    # We iterate through the list of wavebands collecting the absorptances,
    # and waveband names.
    absorptance <- numeric(length(w.band))
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

      # we calculate the average absorptance.
      absorptance[i] <-
        integrate_spct(trim_spct(temp.spct, wb,
                                 use.hinges = FALSE))
    }

    if (quantity %in% c("contribution", "contribution.pc")) {
      total <- absorptance_spct(spct, w.band = NULL,
                                quantity = "total",
                                wb.trim = wb.trim,
                                use.hinges = use.hinges,
                                naming = naming)
      absorptance <- absorptance / total
      if (quantity == "contribution.pc") {
        absorptance <- absorptance * 1e2
      }
    } else if (quantity %in% c("relative", "relative.pc")) {
      total <- sum(absorptance)
      absorptance <- absorptance / total
      if (quantity == "relative.pc") {
        absorptance <- absorptance * 1e2
      }
    } else if (quantity %in% c("average", "mean")) {
      absorptance <- absorptance / sapply(w.band, wl_expanse)
    }

    if (length(absorptance) == 0) {
      absorptance <- NA_real_
      names(absorptance) <- "out of range"
    } else if (naming %in% c("long", "default")) {
      names(absorptance) <- paste(summary.name, wb.name, sep = "_")
    } else if (naming == "short") {
      names(absorptance) <- wb.name
    } else if (naming != "none") {
      warning("Argument to 'naming' unrecognized, assuming \"none\".")
    }

    attr(absorptance, "radiation.unit") <- paste("absorptance", quantity)
    return(absorptance)
  }

# filter_mspct methods -----------------------------------------------

#' @describeIn absorptance Calculates absorptance from a \code{filter_mspct}
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax for \code{attr2tb} passed as is to formal parameter \code{col.names}.
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#'
#' @export
#'
absorptance.filter_mspct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = NULL,
           naming = "default",
           ...,
           attr2tb = NULL,
           idx = "spct.idx" ) {
    z <-
      msdply(
        mspct = spct,
        .fun = absorptance,
        w.band = w.band,
        quantity = quantity,
        wb.trim = wb.trim,
        use.hinges = use.hinges,
        naming = naming,
        idx = idx,
        col.names = names(w.band)
      )
    add_attr2tb(tb = z,
                mspct = spct,
                col.names = attr2tb,
                idx = idx)
  }

# object_mspct methods -----------------------------------------------

#' @describeIn absorptance Calculates absorptance from a \code{object_mspct}
#'
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
absorptance.object_mspct <-
  function(spct, w.band=NULL,
           quantity="average",
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
        .fun = absorptance,
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
