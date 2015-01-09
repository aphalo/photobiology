#' Calculate multipliers for selecting a range of wavelengths and optionally applying a
#' biological spectral weighting function (BSWF) and wavelength normalization.
#'
#' This function gives a set of numeric multipliers that can be used
#' to select a waveband and apply a weight.
#'
#' @usage calc_multipliers(w.length, w.band,
#'                         unit.out="energy", unit.in="energy",
#'                         use.cached.mult=FALSE, fill=0)
#'
#' @param w.length numeric array of wavelength (nm)
#' @param w.band a waveband object (but not a list of waveband objects)
#' @param unit.out a character string: "photon" or "energy", default is "energy"
#' @param unit.in a character string: "photon" or "energy", default is "energy"
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param fill if fill==NA then values returned for wavelengths outside the range of the waveband are set to NA
#'
#' @return a numeric array of multipliers of the same length as \code{w.length}
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' with(sun.data, calc_multipliers(w.length, new_waveband(400,700),"photon"))
#' ## with(sun.data, calc_multipliers(w.length, new_waveband(400,700),"photon", use.cached.mult=TRUE))

calc_multipliers <- function(w.length, w.band, unit.out="energy", unit.in="energy", use.cached.mult=FALSE, fill=0){
  cache.needs.saving <- FALSE
  if (use.cached.mult && !is.null(w.band$name)) {
    # this needs to be changed to something better
    ourEnv <- photobio.cache
    # search for cached multipliers
    cache.name <- paste(w.band$name, as.character(fill), unit.in, unit.out, sep=".")
    if (exists(cache.name, where = ourEnv)) {
      mult <- get(cache.name, envir = ourEnv)
      if (length(w.length) == length(mult)) {
        return(mult)
      } else {
        cache.needs.saving <- TRUE
      }
    } else {
      cache.needs.saving <- TRUE
    }
  }
  mult <- numeric(length(w.length))
  outside.band <- w.length < w.band$low | w.length >= w.band$high
  inside.band <- !outside.band
  mult[outside.band] <- fill
  if (unit.out=="energy"){
    if (unit.in=="energy"){
      mult[inside.band] <- 1.0
    } else if (unit.in=="photon"){
      mult[inside.band] <- 1.0 / e2qmol_multipliers(w.length[inside.band])
    }
    if (!is.null(w.band$weight) && (w.band$weight=="BSWF" || w.band$weight=="SWF")){
      mult[inside.band] <- mult[inside.band] * w.band$SWF.e.fun(w.length[inside.band])
    }
  }
  else if (unit.out=="photon"){
    if (unit.in=="photon"){
      mult[inside.band] <- 1.0
    } else if (unit.in=="energy"){
      mult[inside.band] <- e2qmol_multipliers(w.length[inside.band])
    }
    if (!is.null(w.band$weight) && (w.band$weight=="BSWF" || w.band$weight=="SWF")){
      mult[inside.band] <- mult[inside.band] * w.band$SWF.q.fun(w.length[inside.band])
    }
  }
  if (!is.null(w.band$norm) && (!is.null(w.band$weight))){
    if (w.band$norm >= w.band$low && w.band$norm <= w.band$high){
      norm.divisor <- ifelse(unit.out=="energy", w.band$SWF.e.fun(w.band$norm),  w.band$SWF.q.fun(w.band$norm))
      mult[inside.band] <- mult[inside.band] / norm.divisor
    } else {
      warning("normalization wavelength outside range of SWF")
      return(NA)
    }
  }
  if (use.cached.mult && cache.needs.saving) {
    assign(cache.name, mult, envir = ourEnv)
  }
  return(mult)
}

.onLoad <- function(libname, pkgname) {
  photobio.cache <<- new.env(parent = emptyenv())
}

.onUnload <- function(libpath) {
  rm(photobio.cache, envir = emptyenv())
}
