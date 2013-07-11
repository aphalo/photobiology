#' Calculate photon irradiance from spectral irradiance
#'
#' This function gives the photon irradiance for a given
#' waveband of a radiation spectrum, optionally applies
#' a BSWF.
#'
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral irradiances, by default as energy (W m-2 nm-1)
#' @param w.band list (see \code{\link{new_waveband}} for details)
#' @param unit.in a character string: "photon" or "energy", default is "energy"
#' 
#' @return a single numeric value with no change in scale factor: [W m-2 nm-1] -> [W m-2]
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' with(sun.data, photon_irradiance(w.length, s.e.irrad)) # the whole spectrum
#' with(sun.data, photon_irradiance(w.length, s.e.irrad, PAR())) # photosynthetic photon irradiance
#' with(sun.data, photon_irradiance(w.length, s.e.irrad, waveband(400,700))) # idem
#' with(sun.data, photon_irradiance(w.length, s.e.irrad, PG())) # effective UV photon irradiance

photon_irradiance <- 
  function(w.length, s.irrad, w.band=NULL, unit.in="energy"){
    return(irradiance(w.length, s.irrad, w.band, unit.out="photon", unit.in))
  }
