#' Calculate response from spectral response.
#'
#' This function returns the mean response for a given
#' waveband and a response spectrum.
#'
#' @usage response_intg(spct, w.band=NULL, unit.out="energy", use.hinges=NULL)
#'
#' @param spct an object of class response.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @return a single numeric value expressed either as a fraction of one or a percentage, or a
#' vector of the same length as the list of wave.bands.
#' @keywords manip misc
#' @export
#' @examples
#' # library(photobiologySensors)
#' # data(Vital_BW_20.spct)
#' # response(Vital_BW_20.spct, new_waveband(400,700), unit.out="photon")
#' # response(Vital_BW_20.spct, new_waveband(400,700), unit.out="energy")
#'
#' @note The parameter \code{use.hinges} controls speed optimization. The defaults should be suitable
#' in mosts cases. Only the range of wavelengths in the wavebands is used and all BSWFs are ignored.

response_intg <-
  function(spct, w.band=NULL, unit.out="energy", use.hinges=NULL){
    # if the waveband is undefined then use all data
    if (is.null(w.band)){
      #      w.band <- new_waveband(min(w.length), max(w.length))
      w.band <- new_waveband(min(spct), max(spct) + 1e-4)
      # we need to add a small number as the test is "<"
      # this affects signifcantly the result only when no hinges are used
    }
    if (is(w.band, "waveband")) {
      # if the argument is a single w.band, we enclose it in a list
      # so that the for loop works as expected.This is a bit of a
      # cludge but let's us avoid treating it as a special case
      w.band <- list(w.band)
    }

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

    #!!!!
    if (unit.out=="photon") {
      trimmed_spct <- trimmed_spct[ , list(w.length, s.q.response)]
    } else if (unit.out=="energy") {
      trimmed_spct <- trimmed_spct[ , list(w.length, s.e.response)]
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
    wb_name <- names(w.band)
    no_names_flag <- is.null(wb_name)
    if (no_names_flag) {
      wb_name <- character(length(w.band))
    }
    # "response.spct" objects are not guaranteed to contain response
    # data.
    if (!exists("s.e.response", spct, inherits=FALSE) &&
          ! exists("s.q.response", spct, inherits=FALSE)) {
          warning("No light response data found.")
          return(NA)
        }

    # we iterate through the list of wavebands
    response <- numeric(length(w.band))
    i <- 0
    for (wb in w.band) {
      i <- i + 1
      # we get names from wb if needed
      if (no_names_flag) {
        if (is_effective(wb)) {
          warning("Using only wavelength range from a weighted waveband object.")
          wb_name[i] <- paste("range", as.character(signif(min(wb), 4)), as.character(signif(max(wb), 4)), sep=".")
        } else {
          wb_name[i] <- wb$name
        }
      }
      # we calculate the integrated response.
      trimmed_spct <- trim_spct(temp.spct, wb, use.hinges=FALSE)
      response[i] <- integrate_spct(trim_spct(spct, wb, use.hinges=FALSE))
    }

    names(response) <- paste(names(response), wb_name)
    return(response)
  }

#' Generic function
#'
#' Calculate average response.
#'
#' @param spct an R object of class "generic.spct"
#' @param w.band list of waveband objects
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @export response
#'
response <- function(spct, w.band, unit.out, use.hinges) UseMethod("response")

#' Default for generic function
#'
#' Calculate average response.
#'
#' @param spct an object of class "generic.spct"
#' @param w.band list of waveband objects
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @export reflectance.default
#'
response.default <- function(spct, w.band, unit.out, use.hinges) {
  return(NA)
}

#' Specialization for response.spct
#'
#' Calculate average response.
#'
#' @param spct an object of class "response.spct"
#' @param w.band list of waveband objects
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @export response.response.spct
#'
response.response.spct <- response_intg
