#' Calculate reflectance from spectral reflectance.
#'
#' This function returns the mean reflectance for a given
#' waveband and a reflectance spectrum.
#'
#' @usage reflectance_spct(spct, w.band=NULL, pc.out=FALSE,
#'                         quantity="average", wb.trim=NULL, use.hinges=NULL)
#'
#' @param spct an object of class generic.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param pc.out a logical indicating whether result should be a percentage or a fraction of one
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @return a single numeric value expressed either as a fraction of one or a percentage
#' @keywords manip misc
#' @export
#' @examples
#' # library(photobiologyReflectors)
#' # reflectance(glass_refl.spct, new_waveband(400,700), pc.out=TRUE)
#' # reflectance(glass_refl.spct, new_waveband(400,700), pc.out=FALSE)
#'
#' @note The last parameter controls speed optimization. The defaults should be suitable
#' in mosts cases. Only the range of wavelengths in the wavebands is used and all BSWFs are ignored.

reflectance_spct <-
  function(spct, w.band=NULL, pc.out=FALSE, quantity="average", wb.trim=NULL, use.hinges=NULL){
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
      use.hinges <- ((max(spct) - min(spct)) / length.wl) > 0.7
    }

    # we collect all hinges and insert them in one go
    # this may alter a little the returned values
    # but should be faster
    if (use.hinges) {
      all.hinges <- NULL
      for (wb in w.band) {
        if (!is.null(wb$hinges) && length(wb$hinges)>0) {
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
    reflectance <- numeric(length(w.band))
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
      # we calculate the average reflectance.
      reflectance[i] <- integrate_spct(trim_spct(spct, wb, use.hinges=FALSE))
    }

   if (quantity %in% c("contribution", "contribution.pc")) {
     total <- reflectance_spct(spct, w.band=NULL, pc.out=FALSE,
                                 quantity="total", use.hinges=FALSE)
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
   } else if (quantity == "average") {
     reflectance <- reflectance / sapply(w.band, spread)
     if (pc.out) {
       reflectance <- reflectance * 1e2
       quantity <- paste(quantity, "(%)")
     }
   } else if (quantity == "total") {
     if (pc.out) {
       reflectance <- reflectance * 1e2
       quantity <- paste(quantity, "(%)")
     }
   } else if (quantity != "total") {
     warning("'quantity '", quantity, "' is invalid, returning 'total' instead")
     quantity <- "total"
   }
   if (length(reflectance) == 0) {
     irrad <- NA
     names(reflectance) <- "out of range"
   }
   names(reflectance) <- paste(names(reflectance), wb.name)
   setattr(reflectance, "time.unit", "none")
   setattr(reflectance, "radiation.unit", paste("reflectance", quantity))
   return(reflectance)
  }


#' Generic function
#'
#' Calculate average reflectance.
#'
#' @param spct an object of class "generic.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param pc.out a logical indicating whether result should be a percentage or a fraction of one
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @export reflectance
#'
reflectance <- function(spct, w.band, pc.out, quantity, wb.trim, use.hinges) UseMethod("reflectance")

#' Default for generic function
#'
#' Calculate average reflectance.
#'
#' @param spct an object of class "generic.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param pc.out a logical indicating whether result should be a percentage or a fraction of one
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @export reflectance.default
#'
reflectance.default <- function(spct, w.band, pc.out, quantity, wb.trim, use.hinges) {
  return(NA)
}

#' Specialization for reflector.spct
#'
#' Calculate average reflectance.
#'
#' @param spct an object of class "reflector.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param pc.out a logical indicating whether result should be a percentage or a fraction of one
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @export reflectance.reflector.spct
#'
reflectance.reflector.spct <- reflectance_spct
