#' Gives values for the Plant Growth BSWF as a function of wavelength
#'
#' This function gives a set of numeric multipliers that can be used
#' as a weight to calculate effective doses and irradiances. The
#' returned values are on quantum based effectiveness relative units.
#'
#' @param w.length numeric array of wavelengths (nm)
#'
#' @return a numeric array of the same length as \code{w.length} with values for the BSWF normalized
#' as in the original source (300 nm)
#' 
#' @references \url{http://uv4growth.dyndns.org/}
#' @keywords misc
#' @examples
#' data(sun.data)
#' with(sun.data, PG.q(w.length))

PG.q.fun <-
function(w.length){
#    band.high <- 390 # (nm)
#    within.band <- logical(length(w.length))
#    within.band <- (w.length >= band.low) & (w.length < band.high)
#    if (!any(within.band)){
#      warning("No wavelengths within range of SWF")
#      return(NA)
#    }
#    above.band <- w.length > band.high
#    uncertain <- w.length < band.low
#   spectral_weigths <- numeric(length(w.length))
    spectral_weigths <- 
      exp(4.688272*exp(-exp(0.1703411*(w.length-307.867)/1.15)) +
      ((390-w.length)/121.7557-4.183832))
#     if (any(above.band)){
#       spectral_weigths[above.band] <- 0
#     }
#     if (any(uncertain)){
#       spectral_weigths[uncertain] <- NA # do not extrapolate too far
#     }
    return(spectral_weigths) 
}

