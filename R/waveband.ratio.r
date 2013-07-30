#' Calculate photon (quantum) or energy ratio from spectral (energy or photon) irradiance.
#'
#' This function gives the (energy or photon) irradiance ratio between two given
#' wavebands of a radiation spectrum.
#'
#' @usage waveband_ratio(w.length, s.irrad, w.band.num=NULL, w.band.denom=NULL, 
#' unit.out.num=NULL, unit.out.denom=unit.out.num, unit.in="energy")
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral (energy) irradiances (W m-2 nm-1)
#' @param w.band.num list with elements 'lo' and 'hi' giving the boundaries of the waveband (nm)
#' @param w.band.denom list with elements 'lo' and 'hi' giving the boundaries of the waveband (nm)
#' @param unit.out.num character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param unit.out.denom character string with allowed values "energy", and "photon", or its alias "quantum"
#' @param unit.in character string with allowed values "energy", and "photon", or its alias "quantum"
#' 
#' @return a single numeric value giving the unitless ratio 
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' # photon:photon ratio
#' with(sun.data, waveband_ratio(w.length, s.e.irrad, new_waveband(400,500), new_waveband(400,700), "photon"))
#' # energy:energy ratio
#' with(sun.data, waveband_ratio(w.length, s.e.irrad, new_waveband(400,500), new_waveband(400,700), "energy"))
#' # energy:photon ratio
#' with(sun.data, waveband_ratio(w.length, s.e.irrad, new_waveband(400,700), new_waveband(400,700), "energy", "photon"))
#' # photon:photon ratio waveband : whole spectrum
#' with(sun.data, waveband_ratio(w.length, s.e.irrad, new_waveband(400,500), unit.out.num="photon"))
#' # photon:photon ratio of whole spectrum should be equal to 1.0
#' with(sun.data, waveband_ratio(w.length, s.e.irrad, unit.out.num="photon"))

waveband_ratio <- function(w.length, s.irrad, 
                             w.band.num=NULL, w.band.denom=NULL, 
                             unit.out.num=NULL, unit.out.denom=unit.out.num, 
                             unit.in="energy"){
    # We duplicate code from irradiance() here to avoid repeated checks
    # and calculations on the same data
    #
    # what output? seems safer to not have a default here
    if (is.null(unit.out.num) || is.null(unit.out.denom)){
      warning("'unit.out.num' has no default value")
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
      warning("Error: at least one NA value in wavelengths vector and/or s.irrad vector")
      return(NA)
    }
    # warn if w.length values are not reasonable
    if (min(w.length < 200.0) || max(w.length > 1000.0)){
      warning("Warning: wavelength values should be in nm\n data contains values < 200 nm and/or > 1000 nm")
    }
    # if the waveband for numerator is undefined then use
    # the whole wavelength range of the spectrum for numerator
    if (is.null(w.band.num)){
      w.band.num <- new_waveband(min(w.length),max(w.length))
      warning("'w.band.num' not supplied, using whole range of data instead.")
    }
    # if the waveband for denominator is undefined then use
    # the whole wavelength range of the spectrum for denominator
    if (is.null(w.band.denom)){
      w.band.denom <- new_waveband(min(w.length),max(w.length))
      warning("'w.band.denom' not supplied, using whole range of data instead.")
    }
    # if the w.band.num and/or w.band.denom include 'hinges' we insert them
    # it is o.k. to have hinges unsorted!
    merged.hinges <- c(ifelse(is.null(w.band.denom$hinges), numeric(0), w.band.denom$hinges), 
                       ifelse(is.null(w.band.num$hinges), numeric(0), w.band.num$hinges))
    if (length(merged.hinges) > 0){
      new.data <- insert_hinges(w.length, s.irrad,merged.hinges)
      w.length <- new.data$w.length
      s.e.irrad <- new.data$s.irrad
    }
    # calculate the multipliers
    mult.num <- calc_multipliers(w.length, w.band.num, unit.out.num, unit.in)
    mult.denom <- calc_multipliers(w.length, w.band.denom, unit.out.denom, unit.in)
    
    # calculate weighted spectral irradiance
    irrad.num <- integrate_irradiance(w.length, s.irrad * mult.num)
    irrad.denom <- integrate_irradiance(w.length, s.irrad * mult.denom)
    
    return(irrad.num / irrad.denom)
  }
