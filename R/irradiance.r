#' Calculate photon (quantum) irradiance from spectral (energy) irradiance.
#'
#' This function gives the energy irradiance for a given
#' waveband of a radiation spectrum.
#'
#' @usage irradiance(w.length, s.irrad, w.band=NULL, unit.out=NULL, unit.in="energy", 
#' check.spectrum=TRUE, use.cached.mult=FALSE, use.cpp.code=TRUE)
#' 
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral (energy) irradiances (W m-2 nm-1)
#' @param w.band list with elements 'lo' and 'hi' giving the boundaries of the waveband (nm)
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param unit.in character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param check.spectrum logical indicating whether to sanity check input data, default is TRUE
#' @param use.cached.mult logical indicating whether multiplier values should be cached between calls
#' 
#' @return a single numeric value with no change in scale factor: [W m-2 nm-1] -> [mol s-1 m-2]
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' with(sun.data, irradiance(w.length, s.e.irrad, new_waveband(400,700), "photon"))

irradiance <- 
  function(w.length, s.irrad, w.band=NULL, unit.out=NULL, unit.in="energy", 
           check.spectrum=TRUE, use.cached.mult=FALSE, use.cpp.code=TRUE){
    # what output? seems safer to not have a default here
    if (is.null(unit.out)){
      warning("'unit.out' has no default value")
      return(NA)
    }
    # make code a bit simpler further down
    if (unit.in=="quantum") {unit.in <- "photon"}
    # sanity check for spectral data and abort if check fails
    if (check.spectrum && !check_spectrum(w.length, s.irrad)) {
      return(NA)
    } 
    # if the waveband is undefined then use all data
    if (is.null(w.band)){
      w.band <- new_waveband(min(w.length),max(w.length))
    }
    # if the w.band includes 'hinges' we insert them
    if (!is.null(w.band$hinges)){
      new.data <- insert_hinges(w.length, s.irrad, w.band$hinges)
      w.length <- new.data$w.length
      s.irrad <- new.data$s.irrad
    }
    # calculate the multipliers
    mult <- calc_multipliers(w.length=w.length, w.band=w.band, unit.out=unit.out, 
                             unit.in=unit.in, use.cached.mult=use.cached.mult)
    
    # calculate weighted spectral irradiance
    if (use.cpp.code) {
      irrad <- integrate_irradianceC(w.length, s.irrad * mult)
    } else {
      irrad <- integrate_irradianceR(w.length, s.irrad * mult)
    }
    
    return(irrad)
  }
