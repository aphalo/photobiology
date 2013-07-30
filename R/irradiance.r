#' Calculate photon (quantum) irradiance from spectral (energy) irradiance.
#'
#' This function gives the energy irradiance for a given
#' waveband of a radiation spectrum.
#'
#' @usage irradiance(w.length, s.irrad, w.band=NULL, unit.out=NULL, unit.in="energy")
#' 
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral (energy) irradiances (W m-2 nm-1)
#' @param w.band list with elements 'lo' and 'hi' giving the boundaries of the waveband (nm)
#' @param unit.out character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param unit.in character string with allowed values "energy", and "photon", or its alias "quantum"
#' 
#' @return a single numeric value with no change in scale factor: [W m-2 nm-1] -> [mol s-1 m-2]
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' with(sun.data, irradiance(w.length, s.e.irrad, new_waveband(400,700), "photon"))

irradiance <- 
  function(w.length, s.irrad, w.band=NULL, unit.out=NULL, unit.in="energy"){
    # what output? seems safer to not have a default here
    if (is.null(unit.out)){
      warning("'unit.out' has no default value")
      return(NA)
    }
    # make code a bit simpler further down
    if (unit.in=="quantum") {unit.in <- "photon"}
    # sanity check for wavelengths
    if (is.unsorted(w.length, strictly=TRUE)) {
      warning("Error: wavelengths should be sorted in ascending order")
      return(NA)
    }
    if (length(w.length) != length(s.irrad)){
      warning("Error: wavelengths vector and s.e.irrad vector should have same length")
      return(NA)
    }
    # check for NAs
    if (any(is.na(w.length)|is.na(s.irrad))){
      warning("Error: at least one NA value in wavelengths vector and/or s.e.irrad vector")
      return(NA)
    }
    # warn if w.length values are not reasonable
    if (min(w.length < 200.0) || max(w.length > 1000.0)){
      warning("Warning: wavelength values should be in nm\n data contains values < 200 nm and/or > 1000 nm")
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
    mult <- calc_multipliers(w.length, w.band, unit.out, unit.in)
    
    # calculate weighted spectral irradiance
    irrad <- integrate_irradiance(w.length, s.irrad * mult)
    
    return(irrad)
  }
