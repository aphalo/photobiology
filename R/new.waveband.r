#' Build a "waveband" object that can be used as imput when calculating irradiances.
#'
#' @usage waveband(x, weight=NULL, SWF.e.fun=NULL, SWF.q.fun=NULL, norm=NULL, SWF.norm=NULL, hinges=NULL, wb.name=NULL, wb.label=wb.name)
#'
#' @param x any R object on which applying the function range yields an array of two numeric values, describimg a range of wavelengths (nm)
#' @param weight a character string "SWF" or "BSWF", use NULL (the defalt) to indicate no weighting used
#' when calculating irradiance
#' @param SWF.e.fun a function giving multipliers for a spectral weighting function (energy) as a function of wavelength (nm)
#' @param SWF.q.fun a function giving multipliers for a spectral weighting function (quantum) as a function of wavelength (nm)
#' @param SWF.norm a numeric value giving the native normalization wavelength (nm) used by SWF.e.fun and SWF.q.fun
#' @param norm a single numeric value indicating the wavelength at which the SWF should be normalized
#' to 1.0, in nm. "NULL" means no normalization.
#' @param hinges a numeric array giving the wavelengths at which the s.irrad should be inserted by
#' interpolation, no interpolation is indicated by an empty array (numeric(0)), if NULL then interpolation
#' will take place at both ends of the band.
#' @param wb.name character string giving the name for the waveband defined, default is NULL
#' @param wb.label character string giving the label of the waveband to be used for ploting, default is wb.name
#'
#' @return a list with components low, high, weight, SWF.fun, norm, hinges, name
#' @keywords manip misc
#' @export
#' @exportClass waveband
#' @examples
#' e_irrad(sun.spct, waveband(sun.spct))
#' q_irrad(sun.spct, waveband(c(400,700)))
#'
waveband <- function(x,
                     weight=NULL, SWF.e.fun=NULL, SWF.q.fun=NULL, norm=NULL,
                     SWF.norm=NULL, hinges=NULL, wb.name=NULL, wb.label=wb.name) {
  if (is.generic.spct(x) && is.null(wb.name)) {
    wb.name = "Total"
  }
  x.range <- range(x)
  new_waveband(x.range[1], x.range[2], weight=weight, SWF.e.fun=SWF.e.fun, SWF.q.fun=SWF.q.fun,
               norm=norm, SWF.norm=SWF.norm, hinges=hinges, wb.name=wb.name, wb.label=wb.label)
}

#' Build a "waveband" object that can be used as imput when calculating irradiances.
#'
#' @usage new_waveband(w.low, w.high, weight=NULL, SWF.e.fun=NULL, SWF.q.fun=NULL, norm=NULL, SWF.norm=NULL, hinges=NULL, wb.name=NULL, wb.label=wb.name)
#'
#' @param w.low numeric value, wavelength at the short end of the band (nm)
#' @param w.high numeric value, wavelength at the long end of the band (nm)
#' @param weight a character string "SWF" or "BSWF", use NULL (the defalt) to indicate no weighting used
#' when calculating irradiance
#' @param SWF.e.fun a function giving multipliers for a spectral weighting function (energy) as a function of wavelength (nm)
#' @param SWF.q.fun a function giving multipliers for a spectral weighting function (quantum) as a function of wavelength (nm)
#' @param SWF.norm a numeric value giving the native normalization wavelength (nm) used by SWF.e.fun and SWF.q.fun
#' @param norm a single numeric value indicating the wavelength at which the SWF should be normalized
#' to 1.0, in nm. "NULL" means no normalization.
#' @param hinges a numeric array giving the wavelengths at which the s.irrad should be inserted by
#' interpolation, no interpolation is indicated by an empty array (numeric(0)), if NULL then interpolation
#' will take place at both ends of the band.
#' @param wb.name character string giving the name for the waveband defined, default is NULL
#' @param wb.label character string giving the label of the waveband to be used for ploting, default is wb.name
#'
#' @return a list with components low, high, weight, SWF.fun, norm, hinges, name
#' @keywords manip misc
#' @export
#' @exportClass waveband
#' @examples
#' with(sun.spct, irradiance(w.length, s.e.irrad, new_waveband(400,700), "photon"))
#'
new_waveband <- function(w.low, w.high,
                         weight=NULL, SWF.e.fun=NULL, SWF.q.fun=NULL, norm=NULL,
                         SWF.norm=NULL, hinges=NULL, wb.name=NULL, wb.label=wb.name){
  # we make sure that hinges is not NULL, as this would cause problems elsewhere
  # if we are not using a SWF then we do not need to add hinges as we will be anyway interpolating
  # raw irradiances rather than weighted irradiances
  if (is.null(hinges)) {
    hinges <- c(w.low - 1e-4, w.low, w.high, w.high + 1e-4)
  }
  if (!is.null(weight) && weight!="none") {
    #
    if (!is.null(SWF.e.fun) && is.null(SWF.q.fun)){
      if (!is.null(SWF.norm)){
        SWF.q.fun <- function(w.length){SWF.e.fun(w.length) *  SWF.norm / w.length}
      } else {
        warning("Warning: either both photon and energy SWFs should be supplied, or a value for the
              wavelength at which the function supplied is normalized should be supplied through SWF.norm")
      }
    } else if (!is.null(SWF.q.fun) && is.null(SWF.e.fun)  && !is.null(SWF.norm)){
      if (!is.null(SWF.norm)){
        SWF.e.fun <- function(w.length){SWF.q.fun(w.length) * w.length / SWF.norm}
      } else {
        warning("Warning: either both photon and energy SWFs should be supplied, or a value for the
              wavelength at which the function supplied is normalized should be supplied through SWF.norm")
      }
    } else if (is.null(SWF.e.fun) && is.null(SWF.q.fun)){
      warning("weight != NULL, but no SWFs supplied")
      return(NA)
    }
    if (is.null(wb.name)) wb.name <- paste("range", as.character(signif(w.low, 4)), as.character(signif(w.high, 4)), "wtd", sep=".")
  } else {
    weight <- "none"
    if (is.null(wb.name)) wb.name <- paste("range", as.character(signif(w.low, 4)), as.character(signif(w.high, 4)), sep=".")
  }
  w_band <- list(low = w.low, high = w.high,
                 weight = weight, SWF.e.fun = SWF.e.fun, SWF.q.fun = SWF.q.fun, SWF.norm = SWF.norm,
                 norm = norm, hinges = hinges, name = wb.name, label = wb.label)
  setattr(w_band, "class", c("waveband", class(w_band)))
  return(w_band)
}

#' Build a list of unweighted "waveband" objects that can be used as imput when calculating irradiances.
#'
#' @usage split_bands(x, list.names=NULL, short.names=is.null(wb.names), length.out=NULL)
#'
#' @param x a numeric array of wavelengths to split at (nm), or a range of wavelengths or
#' a generic.spct or a waveband, or a list composed of wavelength ranges (as numeric vectors) for each individual
#' waveband.
#' @param list.names character vector with names for the component wavebands in the returned list
#' @param short.names logical indicating whether to use short or long names for wavebands
#' @param length.out numeric giving the number of regions to split the range into (ignored if w.length is not numeric).
#'
#' @return an un-named list of wabeband objects
#' @keywords manip misc
#' @export
#' @examples
#' split_bands(c(400,500,600))
#' split_bands(list(c(400,500),c(550,650)))
#' split_bands(c(400,500,600), short.names=FALSE)
#' split_bands(c(400,500,600), list.names=c("a","b"))
#' split_bands(c(400,700), length.out=6)
#' split_bands(400:700, length.out=3)
#' split_bands(sun.spct, length.out=10)
#'
#' @note \code{list.names} is used to assign names to the elements of the list, while the waveband objects themselves
#' always retain their \code{wb.label} and \code{wb.name} as generated during their creation.

split_bands <- function(x, list.names=NULL, short.names=is.null(list.names), length.out=NULL) {
  if (is(x, "generic.spct") || is(x, "waveband")) {
    w.length <- range(x)
    w.length[2] <- w.length[2]
  } else if (!is.null(length.out)){
    w.length <- range(x)
  } else {
    w.length <- x
  }
  if (is.numeric(w.length)) {
    unique(sort(w.length))
    wl.len <- length(w.length)
    if (wl.len < 2) {
      warning("At least two wavelength values are needed.")
      return(list())
    } else {
      if (!is.null(length.out)) {
        if (length.out < 1L) {
          return(NA)
        } else {
          wl.len <- length.out + 1
          w.length <- seq(min(w.length), max(w.length), length.out=wl.len)
         }
      } else {
        length.out <- wl.len - 1
      }
      bands.out <- list()
      names <- character()
      for (i in 1:(wl.len - 1)) {
        wb.temp <- new_waveband(w.length[i], w.length[i+1],
                     hinges=NULL, wb.name=NULL)
        names[i] <- labels(wb.temp)[[1]]
        bands.out <- c(bands.out, list(wb.temp))
      }
      num.bands <- wl.len - 1
    }
  } else if (is.list(w.length)) {
    bands.out <- list()
    if (is.null(list.names)) {
      list.names <- names(w.length)
    }
    names <- character()
    i <- 0L
    for (wl.range in w.length) {
      i <- i + 1L
      wb.temp <- new_waveband(min(wl.range), max(wl.range),
                                   hinges=NULL, wb.name=NULL)
      names[i] <- labels(wb.temp)[[1]]
      bands.out <- c(bands.out, list(wb.temp))
    }
    num.bands <- i
  } else {
    warning("Invalid x input.")
    return(NA)
  }
  if (length(bands.out) == 1 && is.null(list.names)) {
    return(bands.out[[1]])
  } else {
    if (short.names) {
      names(bands.out) <- paste("wb", 1:length(bands.out), sep="")
    } else {
      if (!is.null(list.names) && length(list.names) >= num.bands)  {
        names(bands.out) <- list.names[1:num.bands]
      } else {
        names(bands.out) <- names
      }
    }
    return(bands.out)
  }
}

#' Query if it is a waveband
#'
#' Functions to check if an object is waveband.
#'
#' @usage is.waveband(x)
#'
#' @param x any R object
#'
#' @return is.waveband returns TRUE if its argument is a waveband and FALSE otherwise.
#'
#' @export
#'
is.waveband <- function(x) {
  inherits(x, "waveband")
}
