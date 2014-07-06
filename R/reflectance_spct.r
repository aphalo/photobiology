#' Calculate reflectance from spectral reflectance.
#'
#' This function returns the mean reflectance for a given
#' waveband and a reflectance spectrum.
#'
#' @usage reflectance_spct(spct, w.band=NULL, pc.out=FALSE, use.hinges=NULL)
#'
#' @param spct an object of class "generic.spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param pc.out a logical indicating whether result should be a percentage or a fraction of one
#' @param use.hinges logical indicating whether to use hinges to reduce interpolation errors
#'
#' @return a single numeric value expressed either as a fraction of one or a percentage
#' @keywords manip misc
#' @export
#' @examples
#' # reflectance_spct(glass_refl.spct, new_waveband(400,700), pc.out=TRUE)
#' # reflectance_spct(glass_refl.spct, new_waveband(400,700), pc.out=FALSE)
#'
#' @note The last parameter controls speed optimization. The defaults should be suitable
#' in mosts cases. Only the range of wavelengths in the wavebands is used and all BSWFs are ignored.

reflectance_spct <-
  function(spct, w.band=NULL, pc.out=FALSE, use.hinges=NULL){
    # if the waveband is undefined then use all data
    if (is.null(w.band)){
      #      w.band <- new_waveband(min(w.length), max(w.length))
      w.band <- new_waveband(min(spct$w.length), max(spct$w.length) + 1e-4)
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
    wb_name <- names(w.band)
    no_names_flag <- is.null(wb_name)
    if (no_names_flag) {
      wb_name <- character(length(w.band))
    }
    # "filter.spct" objects are not guaranteed to contain transmittance
    # expressed in the needed scale, we add the needed columns and as
    # spectra are passed by reference they propagate to the argument
    if (pc.out) {
      if (with(spct, !exists("Rpc"))) {
        if (with(spct, exists("Rfr"))) {
          spct[ , Rpc := Rfr * 100]
        } else {
          warning("No light source data found.")
          return(NA)
        }
      }
    } else {
      if (with(spct, exists("Rfr"))) {
        if (with(spct, exists("Rpc"))) {
          spct[ , Rfr := Rpc / 100]
        } else {
          warning("No light source data found.")
          return(NA)
        }
      }
    }

    #
    if (pc.out) {
      spct.cols <- spct[ , list(w.length, Rpc)]
    } else {
      spct.cols <- spct[ , list(w.length, Rfr)]
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
          wb_name[i] <- paste("range", as.character(signif(min(wb), 4)), as.character(signif(max(wb), 4)), sep=".")
        } else {
          wb_name[i] <- wb$name
        }
      }
      # we calculate the average transmittance.
      reflectance[i] <- average_spct(trim_spct(spct.cols, wb, use.hinges=FALSE))
    }

    names(reflectance) <- paste(names(reflectance), wb_name)
    return(reflectance)
  }
