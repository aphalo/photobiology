#' Calculate moles of photons per Joule ratio from spectral (energy or photon) irradiance.
#'
#' This function gives the photons:energy ratio between for one given
#' waveband of a radiation spectrum.
#'
#' @usage photons_energy_ratio(w.length, s.irrad, w.band=NULL, unit.in="energy")
#' 
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral (energy) irradiances (W m-2 nm-1)
#' @param w.band list with elements 'lo' and 'hi' giving the boundaries of the waveband (nm)
#' @param unit.in character string with allowed values "energy", and "photon", or its alias "quantum"
#' 
#' @return a single numeric value giving the ratio [moles photons] per Joule.
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' # photons:energy ratio
#' with(sun.data, photons_energy_ratio(w.length, s.e.irrad, new_waveband(400,500)))
#' # photons:energy ratio for whole spectrum
#' with(sun.data, photons_energy_ratio(w.length, s.e.irrad))

photons_energy_ratio <- function(w.length, s.irrad, 
                           w.band=NULL, 
                           unit.in="energy"){
  return(waveband_ratio(w.length, s.irrad, w.band, w.band, 
                        unit.out.num="photon", unit.out.denom="energy", 
                        unit.in=unit.in))
}