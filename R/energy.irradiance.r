#' Calculate (energy) irradiance from spectral irradiance.
#'
#' This function gives the energy irradiance for a given
#' waveband of a radiation spectrum, optionally applies
#' a BSWF.
#'
#' @param w.length numeric array of wavelength (nm)
#' @param s.e.irrad numeric array of spectral (energy) irradiances (W m-2 nm-1)
#' @param w.band list with elements 'lo' and 'hi' giving the boundaries of the waveband (nm)
#' 
#' @return a single numeric value with no change in scale factor: [W m-2 nm-1] -> [W m-2]
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' with(sun.data, energy.irradiance(w.length, s.e.irrad)) # the whole spectrum
#' with(sun.data, energy.irradiance(w.length, s.e.irrad, PAR)) # photosyntheti energy irradiance
#' with(sun.data, energy.irradiance(w.length, s.e.irrad, waveband(400,700))) # idem
#' with(sun.data, energy.irradiance(w.length, s.e.irrad, CIE)) # effective UV irradiance

energy.irradiance <- 
  function(w.length, s.irrad, w.band=NULL, unit.in="energy"){
    return(irradiance(w.length, s.irrad, w.band, unit.out="energy", unit.in))
  }
