#' Calculate response from spectral response.
#'
#' This function returns the mean response for a given
#' waveband and a response spectrum.
#'
#' @usage response_intg(spct, w.band=NULL, unit.out, use.hinges=NULL)
#'
#' @param spct an object of class response.spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @return a single numeric value expressed either as a fraction of one or a percentage, or a
#' vector of the same length as the list of wave.bands.
#' @keywords manip misc
#' @export
#' @examples
#' library(photobiologySensors)
#' response_intg(Vital_BW_20.spct, new_waveband(200,300), unit.out="photon") * 1e-6
#' response_intg(Vital_BW_20.spct, new_waveband(200,300), unit.out="energy")
#' response_intg(Vital_BW_20.spct, unit.out="energy")
#' response_intg(Vital_BW_20.spct, unit.out="photon") * 1e-6
#'
#' @note The parameter \code{use.hinges} controls speed optimization. The defaults should be suitable
#' in mosts cases. Only the range of wavelengths in the wavebands is used and all BSWFs are ignored.

response_intg <-
  function(spct, w.band=NULL, unit.out, use.hinges=NULL){
    # "response.spct" objects are not guaranteed to contain response
    # data.
    if (!exists("s.e.response", spct, inherits=FALSE) &&
          ! exists("s.q.response", spct, inherits=FALSE)) {
      warning("No spectral response data found.")
      return(NA)
    }
    # makes "quantum" synonym for "photon" without changes to other code
    if (unit.out == "quantum") {
      unit.out <- "photon"
    }
    if (unit.out=="photon") {
      spct <- e2q(spct)
      spct <- spct[ , list(w.length, s.q.response)]
    } else if (unit.out=="energy") {
      spct <- q2e(spct)
      spct <- spct[ , list(w.length, s.e.response)]
    } else {
      stop("Invalid 'unit.out'")
    }

    # if the waveband is undefined then use all data
    if (is.null(w.band)){
      w.band <- list(Total=waveband(range(spct), wb.name="Total"))
    }
    if (is.waveband(w.band)) {
      # if the argument is a single w.band, we enclose it in a list
      # so that the for loop works as expected. This is a bit of a
      # cludge but it let's us avoid treating it as a special case
      w.band <- list(w.band)
    }

    # if the w.band includes 'hinges' we insert them,
    # but if not, we decide whether to insert hinges or not
    # hinges or not based of the wavelength resolution of the
    # spectrum. This can produce small errors for high
    # spectral resolution data, but speed up the calculations.
    if (is.null(use.hinges)) {
      length.wl <- length(spct$w.length)
      use.hinges <- ((max(spct) - min(spct)) / length.wl) > 0.7
    }

    # we collect all hinges and insert them in one go
    # this may alter very slightly the returned values
    # but improves calculation speed
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

    # we iterate through the list of wavebands
    response <- double(length(w.band))
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
      response[i] <- integrate_spct(trim_spct(spct, wb, use.hinges=FALSE))
    }

    names(response) <- paste(names(response), wb_name)
    return(response)
  }

#' Calculate energy-based photo-response from spectral response.
#'
#' This function returns the mean response for a given
#' waveband and a response spectrum.
#'
#' @usage e_response.response.spct(spct, w.band=NULL, use.hinges=NULL)
#'
#' @param spct an object of class response.spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @return a single numeric value expressed either as a fraction of one or a percentage, or a
#' vector of the same length as the list of wave.bands.
#' @keywords manip misc
#' @export
#' @examples
#' library(photobiologySensors)
#' e_response(Vital_BW_20.spct, new_waveband(200,300))
#' e_response(Vital_BW_20.spct)
#'
#' @note The parameter \code{use.hinges} controls speed optimization. The defaults should be suitable
#' in mosts cases. Only the range of wavelengths in the wavebands is used and all BSWFs are ignored.

e_response.response.spct <-
  function(spct, w.band=NULL, use.hinges=NULL){
    response_intg(spct=spct, w.band=w.band, unit.out="energy", use.hinges=use.hinges)
  }

#' Calculate photon-based photo-response from spectral response.
#'
#' This function returns the mean response for a given
#' waveband and a response spectrum.
#'
#' @usage q_response.response.spct(spct, w.band=NULL, use.hinges=NULL)
#'
#' @param spct an object of class response.spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @return a single numeric value expressed either as a fraction of one or a percentage, or a
#' vector of the same length as the list of wave.bands.
#' @keywords manip misc
#' @export
#' @examples
#' library(photobiologySensors)
#' q_response(Vital_BW_20.spct, new_waveband(200,300)) * 1e-6
#' q_response(Vital_BW_20.spct) * 1e-6
#'
#' @note The parameter \code{use.hinges} controls speed optimization. The defaults should be suitable
#' in mosts cases. Only the range of wavelengths in the wavebands is used and all BSWFs are ignored.

q_response.response.spct <-
  function(spct, w.band=NULL, use.hinges=NULL){
    response_intg(spct=spct, w.band=w.band, unit.out="photon", use.hinges=use.hinges)
  }

#' Generic function
#'
#' Calculate average energy-based photo-response.
#'
#' @param spct an R object of class "generic.spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @export e_response
#'
e_response <- function(spct, w.band, use.hinges) UseMethod("e_response")

#' Generic function
#'
#' Calculate average photon-based photo-response.
#'
#' @param spct an R object of class "generic.spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @export q_response
#'
q_response <- function(spct, w.band, use.hinges) UseMethod("q_response")

#' Default for generic function
#'
#' Calculate average energy-based photo-response.
#'
#' @param spct an object of class "generic.spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @export reflectance.default
#'
e_response.default <- function(spct, w.band, use.hinges) {
  return(NA)
}

#' Default for generic function
#'
#' Calculate average photon-based photo-response.
#'
#' @param spct an object of class "generic.spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#' @export reflectance.default
#'
q_response.default <- function(spct, w.band, use.hinges) {
  return(NA)
}
