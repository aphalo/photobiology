#' Build a "waveband" object that can be used as imput when calculating irradiances.
#'
#' @param w.low numeric value, wavelength at the short end of the band (nm)
#' @param w.high numeric value, wavelength at the long end of the band (nm)
#' @param weight a character string "SWF" or "BSWF", use NULL (the defalt) to indicate no weighting used
#' when calculating irradiance
#' @param SWF.fun a function giving multipliers for a spectral weighting function as a function of wavelength (nm)
#' @param norm a single numeric value indicating the wavelength at which the SWF should be normalized
#' to 1.0, in nm. "NULL" means no normalization.
#' @param hinges a numeric array giving the wavelengths at which the s.irrad should be inserted by
#' interpolation, no interpolation is indicated by an empty array (numeric(0))
#' 
#' @return a list with components low, high, weight, SWF.fun, norm, hinges
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' with(sun.data, irradiance(w.length, s.e.irrad, new_waveband(400,700), "photon"))
#' 
new_waveband <- function(w.low, w.high, weight=NULL, SWF.fun=NULL, norm=NULL, hinges=NULL){
  return(list(low=w.low, high=w.high, weight=weight, SWF.fun=SWF.fun, norm=norm, hinges=hinges))
} 