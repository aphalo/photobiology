#' Calculate absorptance from spectral absorptance.
#'
#' This function returns the mean absorptance for a given
#' waveband of a absorptance spectrum.
#'
#' @usage absorptance_spct(spct, w.band=NULL, quantity="average", wb.trim=NULL, use.hinges=NULL)
#' @usage absorptance(spct, w.band=NULL, use.hinges=NULL)
#'
#' @param spct an object of class "object.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @return a single numeric value with no change in scale factor: AU (absorptance units, using log10)
#' @keywords manip misc
#' @export absorptance_spct absorptance.object.spct
#' @examples
#' library(photobiologyFilters)
#' absorptance(polyester.new.spct, new_waveband(400,700))
#'
#' @note The last parameter controls speed optimization. The defaults should be suitable
#' in mosts cases. Only the range of wavelengths in the wavebands is used and all BSWFs are ignored.

absorptance_spct <-
  function(spct, w.band=NULL, quantity="average", wb.trim=NULL, use.hinges=NULL){
    if (is.normalized(spct) || is.rescaled(spct)) {
      warning("The espectral data has been normalized or rescaled, making impossible to calculate absorptance")
      return(NA)
    }
    spct <- copy(spct)
    # we calculate absorptance
    Tfr.type <- getTfrType(spct)
    Rfr.type <- getRfrType(spct)
    if (Tfr.type == "total" && Rfr.type == "total") {
      spct[ , Afr := 1 - Tfr - Rfr]
    } else if (Tfr.type == "internal" && Rfr.type == "total") {
      spct[ , Afr == (1 - Tfr) * (1 - Rfr)]
    } else if (Tfr.type == "unknown" || Rfr.type == "unknown") {
      warning("'unknown' Tfr.type or Rfr.type, skipping absorptance calculation")
      absorptance <- NA
      setattr(absorptance, "Afr.type", "unknown")
      setattr(absorptance, "radiation.unit", paste("absorptance", quantity))
      return(absorptance)
    } else if (Rfr.type == "specular") {
      warning("'specular' Rfr.type, skipping absorptance calculation")
      absorptance <- NA
      setattr(absorptance, "Afr.type", "unknown")
      setattr(absorptance, "radiation.unit", paste("absorptance", quantity))
      return(absorptance)
    } else {
      stop("Failed assertion with Tfr.type: ", Tfr.type, "and Rfr.type: ", Rfr.type)
    }
    spct <- spct[ , .(w.length, Afr)] # data.table removes attributes!
    setGenericSpct(spct)
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
      length.wl <- length(spct$w.length)
      use.hinges <- (spct$w.length[length.wl] - spct$w.length[1]) / length.wl > 0.7 #
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
    #
    # we iterate through the list of wavebands
    absorptance <- numeric(length(w.band))
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
      absorptance[i] <- integrate_spct(trim_spct(spct, wb, use.hinges=FALSE))
    }

    if (quantity %in% c("contribution", "contribution.pc")) {
      total <- absorptance_spct(spct, w.band=NULL,
                                  quantity="total", use.hinges=FALSE)
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
      absorptance <- absorptance / sapply(w.band, spread)
    }
    if (length(absorptance) == 0) {
      absorptance <- NA
      names(absorptance) <- "out of range"
    }
    names(absorptance) <- paste(names(absorptance), wb.name)
    setattr(absorptance, "Afr.type", "total")
    setattr(absorptance, "radiation.unit", paste("absorptance", quantity))
    return(absorptance)
  }

#' Generic function
#'
#' Calculate average absorptance.
#'
#' @param spct an object of class "generic.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param quantity character string
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @export absorptance
#'
absorptance <- function(spct, w.band, quantity, wb.trim, use.hinges) UseMethod("absorptance")

#' Default for generic function
#'
#' Calculate average absorptance.
#'
#' @param spct an object of class "generic.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @export absorptance.default
#'
absorptance.default <- function(spct, w.band, quantity, wb.trim, use.hinges) {
  warning("'absorptance' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' Specialization for object.spct
#'
#' Calculate average absorptance
#'
#' @param spct an object of class "object.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @export absorptance.object.spct
#'
absorptance.object.spct <- absorptance_spct
