#' Calculate absorbance from spectral absorbance.
#'
#' This function returns the mean absorbance for a given
#' waveband of a absorbance spectrum.
#'
#' @usage absorbance_spct(spct, w.band=NULL, quantity="average", wb.trim=NULL, use.hinges=NULL)
#' @usage absorbance(spct, w.band=NULL, use.hinges=NULL)
#'
#' @param spct an object of class "filter.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @return a single numeric value with no change in scale factor: AU (absorbance units, using log10)
#' @keywords manip misc
#' @export absorbance_spct absorbance.filter.spct
#' @examples
#' library(photobiologyFilters)
#' absorbance(polyester.new.spct, new_waveband(400,700))
#'
#' @note The last parameter controls speed optimization. The defaults should be suitable
#' in mosts cases. Only the range of wavelengths in the wavebands is used and all BSWFs are ignored.

absorbance_spct <-
  function(spct, w.band=NULL, quantity="average", wb.trim=NULL, use.hinges=NULL){
    if (is.normalized(spct) || is.rescaled(spct)) {
      warning("The espectral data has been normalized or rescaled, making impossible to calculate absorbance")
      return(NA)
    }
    spct <- T2A(spct, action="replace", byref=FALSE)
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
    absorbance <- numeric(length(w.band))
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
      absorbance[i] <- integrate_spct(trim_spct(spct, wb, use.hinges=FALSE))
    }

    if (quantity %in% c("contribution", "contribution.pc")) {
      total <- absorbance_spct(spct, w.band=NULL,
                                  quantity="total", use.hinges=FALSE)
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
    } else if (quantity == "average") {
      absorbance <- absorbance / sapply(w.band, spread)
    }
    if (length(absorbance) == 0) {
      irrad <- NA
      names(absorbance) <- "out of range"
    }
    names(absorbance) <- paste(names(absorbance), wb.name)
    setattr(absorbance, "Tfr.type", getTfrType(spct))
    setattr(absorbance, "radiation.unit", paste("absorbance", quantity))
    return(absorbance)
  }

#' Generic function
#'
#' Calculate average absorbance.
#'
#' @param spct an object of class "generic.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param quantity character string
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @export absorbance
#'
absorbance <- function(spct, w.band, quantity, wb.trim, use.hinges) UseMethod("absorbance")

#' Default for generic function
#'
#' Calculate average absorbance.
#'
#' @param spct an object of class "generic.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @export absorbance.default
#'
absorbance.default <- function(spct, w.band, quantity, wb.trim, use.hinges) {
  return(NA)
}

#' Specialization for filter.spct
#'
#' Calculate average absorbance
#'
#' @param spct an object of class "filter.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @export absorbance.filter.spct
#'
absorbance.filter.spct <- absorbance_spct
