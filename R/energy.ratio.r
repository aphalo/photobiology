#' Calculate energy ratio from spectral (energy or photon) irradiance.
#'
#' This function gives the energy ratio between two given
#' wavebands of a radiation spectrum.
#'
#' @usage energy_ratio(w.length, s.irrad, w.band.num=NULL, w.band.denom=NULL, unit.in="energy")
#'
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral (energy) irradiances (W m-2 nm-1)
#' @param w.band.num list with elements 'lo' and 'hi' giving the boundaries of the waveband (nm)
#' @param w.band.denom list with elements 'lo' and 'hi' giving the boundaries of the waveband (nm)
#' @param unit.in character string with allowed values "energy", and "photon", or its alias "quantum"
#' 
#' @return a single numeric value giving the unitless ratio 
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' # energy:energy ratio
#' with(sun.data, energy_ratio(w.length, s.e.irrad, new_waveband(400,500), new_waveband(400,700)))
#' # energy:energy ratio waveband : whole spectrum
#' with(sun.data, energy_ratio(w.length, s.e.irrad, new_waveband(400,500)))
#' # energy:energy ratio of whole spectrum should be equal to 1.0
#' with(sun.data, energy_ratio(w.length, s.e.irrad))

energy_ratio <- function(w.length, s.irrad, 
                           w.band.num=NULL, w.band.denom=NULL, 
                           unit.in="energy"){
  return(waveband_ratio(w.length, s.irrad, w.band.num, w.band.denom, 
                        unit.out.num="energy", unit.out.denom="energy", 
                        unit.in=unit.in))
}