#' Calculate transmittance from spectral transmittance.
#'
#' This function returns the summary transmittance for given wavebands from a
#' filter or object spectrum.
#'
#' @usage transmittance(spct, w.band, pc.out, quantity, wb.trim, use.hinges)
#'
#' @param spct an R object
#' @param w.band waveband or list of waveband objects The waveband(s) determine
#'   the region(s) of the spectrum that are summarized.
#' @param pc.out logical Flag indicating whether result should be a percentage or a
#'   fraction of one
#' @param quantity character
#' @param wb.trim logical Flag indicating if wavebands crossing spectral data boundaries
#'   are trimmed or ignored
#' @param use.hinges logical Flag indicating whether to use hinges to reduce
#'   interpolation errors
#'
#' @return A numeric vector with no change in scale factor (expressed as percent
#'   or fraction)
#' @keywords manip misc
#' @export transmittance
#' @examples
#' library(photobiologyFilters)
#' transmittance(polyester.new.spct, new_waveband(400,700), pc.out=TRUE)
#' transmittance(polyester.new.spct, new_waveband(400,700), pc.out=FALSE)
#'
#' @note The \code{use.hinges} parameter controls speed optimization. The
#'   defaults should be suitable in mosts cases. Only the range of wavelengths
#'   in the wavebands is used and all BSWFs are ignored.
#'
#' @export transmittance
#'
transmittance <- function(spct, w.band, pc.out, quantity, wb.trim, use.hinges) UseMethod("transmittance")

#' @describeIn transmittance Default method
#'
#' @export
#'
transmittance.default <- function(spct, w.band, pc.out, quantity, wb.trim, use.hinges) {
  return(NA)
}

#' @describeIn transmittance Method for filter spectra
#'
#' @export
#'
transmittance.filter_spct <-
  function(spct, w.band=NULL,
           pc.out=FALSE,
           quantity="average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL) ) {
    transmittance_spct(spct = spct,
                       w.band = w.band,
                       pc.out = pc.out,
                       quantity = quantity,
                       wb.trim = wb.trim,
                       use.hinges = use.hinges)
  }

#' @describeIn transmittance Method for object spectra
#'
#' @export
#'
transmittance.object_spct <-
  function(spct, w.band=NULL,
           pc.out=FALSE,
           quantity="average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL) ) {
    transmittance_spct(spct = spct,
                       w.band = w.band,
                       pc.out = pc.out,
                       quantity = quantity,
                       wb.trim = wb.trim,
                       use.hinges = use.hinges)
  }

#' Calculate transmittance from spectral transmittance.
#'
#' This function returns the mean transmittance for a given waveband of a
#' transmittance spectrum.
#'
#' @usage transmittance_spct(spct, w.band, pc.out, quantity, wb.trim, use.hinges)
#'
#' @param spct an object of class "generic_spct"
#' @param w.band waveband or list of waveband objects The waveband(s) determine
#'   the region(s) of the spectrum that are summarized.
#' @param pc.out logical Flag indicating whether result should be a percentage or a
#'   fraction of one
#' @param quantity character
#' @param wb.trim logical Flag indicating if wavebands crossing spectral data boundaries
#'   are trimmed or ignored
#' @param use.hinges logical Flag indicating whether to use hinges to reduce
#'   interpolation errors
#'
#' @return a single numeric value
#' @keywords internal
#'
#' @note The last parameter controls speed optimization. The defaults should be
#'   suitable in mosts cases. Only the range of wavelengths in the wavebands is
#'   used and all BSWFs are ignored.

transmittance_spct <-
  function(spct, w.band, pc.out, quantity, wb.trim, use.hinges) {
    if (is.normalized(spct) || is.rescaled(spct)) {
      warning("The espectral data has been normalized or rescaled, making impossible to calculate transmittance")
      return(NA)
    }
    if (is_filter_spct(spct)) {
      spct <- A2T(spct, action="replace", byref=FALSE)
    } else {
      spct <- as_filter_spct(spct)
    }
    Tfr.type <- getTfrType(spct)
    spct <- spct[ , .(w.length, Tfr)] # data.table removes attributes!
    setTfrType(spct, Tfr.type = Tfr.type)
    # if the waveband is undefined then use all data
    if (is.null(w.band)){
      w.band <- waveband(spct)
    }
    if (is(w.band, "waveband")) {
      # if the argument is a single w.band, we enclose it in a list
      # so that the for loop works as expected.This is a bit of a
      # cludge but let's us avoid treating it as a special case
      w.band <- list(w.band)
    }
    w.band <- trim_waveband(w.band=w.band, range=spct, trim=wb.trim)

    # if the w.band includes 'hinges' we insert them
    # choose whether to use hinges or not
    # if the user has specified its value, we leave it alone
    # but if it was not requested, we decide whether to insert
    # hinges or not based of the wavelength resolution of the
    # spectrum. This will produce small errors for high
    # spectral resolution data, and speed up the calculations
    # a lot in such cases
    if (is.null(use.hinges)) {
      use.hinges <- stepsize(spct)[2] > getOption("photobiology.auto.hinges.limit", default = 0.5) # nm
    }

    # we collect all hinges and insert them in one go
    # this may alter a little the returned values
    # but should be faster
    if (use.hinges) {
      all.hinges <- NULL
      for (wb in w.band) {
        if (!is.null(wb$hinges) & length(wb$hinges)>0) {
          all.hinges <- c(all.hinges, wb$hinges)
        }
      }
      if (!is.null(all.hinges)) {
        spct <- insert_spct_hinges(spct, all.hinges)
      }
    }

    # we prepare labels for output
    wb.name <- names(w.band)
    no_names_flag <- is.null(wb.name)
    if (no_names_flag) {
      wb.name <- character(length(w.band))
    }
    # we iterate through the list of wavebands
    transmittance <- numeric(length(w.band))
    i <- 0
    for (wb in w.band) {
      i <- i + 1
      # we get names from wb if needed
      if (no_names_flag) {
        if (is_effective(wb)) {
          warning("Using only wavelength range from a weighted waveband object.")
          wb.name[i] <- paste("range", as.character(signif(min(wb), 4)), as.character(signif(max(wb), 4)), sep=".")
        } else {
          wb.name[i] <- wb$name
        }
      }
      # we calculate the average transmittance.
      transmittance[i] <- integrate_spct(trim_spct(spct, wb, use.hinges=FALSE))
    }

    if (quantity %in% c("contribution", "contribution.pc")) {
      total <- transmittance_spct(spct, w.band=NULL, pc.out=FALSE,
                                quantity="total", use.hinges=FALSE)
      transmittance <- transmittance / total
      if (quantity == "contribution.pc") {
        transmittance <- transmittance * 1e2
      }
    } else if (quantity %in% c("relative", "relative.pc")) {
      total <- sum(transmittance)
      transmittance <- transmittance / total
      if (quantity == "relative.pc") {
        transmittance <- transmittance * 1e2
      }
    } else if (quantity %in% c("average", "mean")) {
      transmittance <- transmittance / sapply(w.band, spread)
      if (pc.out) {
        transmittance <- transmittance * 1e2
        quantity <- paste(quantity, "(%)")
      }
    } else if (quantity == "total") {
      if (pc.out) {
        transmittance <- transmittance * 1e2
        quantity <- paste(quantity, "(%)")
      }
    } else if (quantity != "total") {
      warning("'quantity '", quantity, "' is invalid, returning 'total' instead")
      quantity <- "total"
    }

    if (length(transmittance) == 0) {
      transmittance <- NA
      names(transmittance) <- "out of range"
    }
    names(transmittance) <- paste(names(transmittance), wb.name)
    setattr(transmittance, "Tfr.type", getTfrType(spct))
    setattr(transmittance, "radiation.unit", paste("transmittance", quantity))
    return(transmittance)
  }

