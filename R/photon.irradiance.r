#' Calculate photon irradiance from spectral irradiance
#'
#' This function gives the photon irradiance for a given
#' waveband of a radiation spectrum, optionally applies
#' a BSWF.
#' 
#' @usage photon_irradiance(w.length, s.irrad, w.band=NULL, unit.in="energy",
#' check.spectrum=TRUE, use.cached.mult=FALSE, use.cpp.code=TRUE)
#'
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral irradiances, by default as energy (W m-2 nm-1)
#' @param w.band list (see \code{\link{new_waveband}} for details)
#' @param unit.in a character string: "photon" or "energy", default is "energy"
#' @param check.spectrum logical indicating whether to sanity check input data, default is TRUE
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' @param use.cpp.code logical indicating whether to use compiled C++ function for integartion
#' 
#' @return a single numeric value with no change in scale factor: [W m-2 nm-1] -> [W m-2]
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' 
#' with(sun.data, photon_irradiance(w.length, s.e.irrad))
#' with(sun.data, photon_irradiance(w.length, s.e.irrad, new_waveband(400,700)))

photon_irradiance <- 
  function(w.length, s.irrad, w.band=NULL, unit.in="energy", 
           check.spectrum=TRUE, use.cached.mult=FALSE, use.cpp.code=TRUE){
    return(irradiance(w.length=w.length, s.irrad=s.irrad, w.band=w.band, 
                      unit.out="photon", unit.in=unit.in, 
                      check.spectrum=check.spectrum, use.cached.mult=use.cached.mult, use.cpp.code=use.cpp.code))
  }
