#' Build a "waveband" object that can be used as imput when calculating irradiances.
#'
#' @usage new_waveband(w.low, w.high, weight=NULL, SWF.e.fun=NULL, SWF.q.fun=NULL, norm=NULL, SWF.norm=NULL, hinges=c(w.low-0.01,w.low,w.high-0.01,w.high), wb.name=NULL)
#' 
#' @param w.low numeric value, wavelength at the short end of the band (nm)
#' @param w.high numeric value, wavelength at the long end of the band (nm)
#' @param weight a character string "SWF" or "BSWF", use NULL (the defalt) to indicate no weighting used
#' when calculating irradiance
#' @param SWF.e.fun a function giving multipliers for a spectral weighting function (energy) as a function of wavelength (nm)
#' @param SWF.q.fun a function giving multipliers for a spectral weighting function (quantum) as a function of wavelength (nm)
#' @param SWF.norm a numeric value giving the native normalization wavelength (nm) used by SWF.e.fun and SWF.q.fun
#' @param norm a single numeric value indicating the wavelength at which the SWF should be normalized
#' to 1.0, in nm. "NULL" means no normalization.
#' @param hinges a numeric array giving the wavelengths at which the s.irrad should be inserted by
#' interpolation, no interpolation is indicated by an empty array (numeric(0))
#' @param wb.name character string giving the name for the waveband defined, default is "no.name"
#' 
#' @return a list with components low, high, weight, SWF.fun, norm, hinges, name
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' with(sun.data, irradiance(w.length, s.e.irrad, new_waveband(400,700), "photon"))
#' 
new_waveband <- function(w.low, w.high, 
                         weight=NULL, SWF.e.fun=NULL, SWF.q.fun=NULL, norm=NULL, 
                         SWF.norm=NULL, hinges=c(w.low-0.01,w.low,w.high-0.01,w.high), wb.name=NULL){
  if (!is.null(weight)) {
    # 
    if (!is.null(SWF.e.fun) && is.null(SWF.q.fun)){
      if (!is.null(SWF.norm)){
        SWF.q.fun <- function(w.length,s.e.irrad){SWF.e.fun(w.length) *  SWF.norm / w.length}
      } else {
        warning("Warning: either both photon and energy SWFs should be supplied, or a value for the 
              wavelength at which the function supplied is normalized should be supplied through SWF.norm")
      }
    } else if (!is.null(SWF.q.fun) && is.null(SWF.e.fun)  && !is.null(SWF.norm)){
      if (!is.null(SWF.norm)){
        SWF.e.fun <- function(w.length,s.e.irrad){SWF.q.fun(w.length) * w.length / SWF.norm}
      } else {
        warning("Warning: either both photon and energy SWFs should be supplied, or a value for the 
              wavelength at which the function supplied is normalized should be supplied through SWF.norm")
      }
    } else if (is.null(SWF.e.fun) && is.null(SWF.q.fun)){
      warning("weight != NULL, but no SWFs supplied")
      return(NA)
    }
  }
  return(list(low=w.low, high=w.high, weight=weight, SWF.e.fun=SWF.e.fun, SWF.q.fun=SWF.q.fun, SWF.norm=SWF.norm, norm=norm, hinges=hinges, name=wb.name))
} 