#' Calculate multipliers for selecting a range of wavelengths and optionally applying a
#' biological spectral weighting function (BSWF) and wavelength normalization.
#'
#' This function gives a set of numeric multipliers that can be used
#' to select a waveband and apply a weight.
#'
#' @usage calc_multipliers(w.length, w.band, unit.out="energy", unit.in="energy")
#' 
#' @param w.length numeric array of wavelength (nm)
#' @param w.band list(low, high, weight, BSWF.fun, norm)
#' @param unit.out a character string: "photon" or "energy", default is "energy"
#' @param unit.in a character string: "photon" or "energy", default is "energy"
#'
#' @return a numeric array of multipliers of the same length as \code{w.length}
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' with(sun.data, calc_multipliers(w.length, new_waveband(400,700),"photon"))

calc_multipliers <- function(w.length, w.band, unit.out="energy", unit.in="energy"){
  mult <- numeric(length(w.length))
  outside.band <- w.length < w.band$low | w.length >= w.band$high
  inside.band <- !outside.band
  mult[outside.band] <- 0.0
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
  return(mult)
}

  