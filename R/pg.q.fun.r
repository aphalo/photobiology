#' Gives values for the Plant Growth BSWF as a function of wavelength
#'
#' This function gives a set of numeric multipliers that can be used
#' as a weight to calculate effective doses and irradiances. The
#' returned values are on quantum based effectiveness relative units.
#'
#' @param w.length numeric array of wavelengths (nm)
#'
#' @return a numeric array of the same length as \cde{w.length} with values for the BSWF
#' 
#' @references \url{http://uv4growth.dyndns.org/}
#' @keywords misc
#' @examples
#' data(sun.data)
#' with(sun.data, PG.q(w.length))

PG.q <-
function(w.length,w.norm=300){
    band.low <- 250 # (nm)
    band.high <- 390 # (nm)
    wihin.band <- (w.length >= band.low) & (w.length < band.high)
    above.band <- w.length > band.high
    uncertain <- w.length < band.low
    weigths <- numeric(length(w.length))
    weigths[within.band] <- 
      exp(4.688272*exp(-exp(0.1703411*(w.length[within.band]-307.867)/1.15)) +
      ((390-w.length[within.band])/121.7557-4.183832))
    weights[above.band] <- 0
    weights[uncertain] <- NA # do not extrapolate too far
    return(weights) 
}

