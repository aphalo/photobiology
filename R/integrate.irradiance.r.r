#' Gives irradiance from spectral irradiance.
#'
#' This function gives the result of integrating spectral irradiance over
#' wavelengths.
#' 
#' @usage integrate_irradianceR(w.length, s.irrad)
#' 
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral irradiances
#' 
#' @return a single numeric value with no change in scale factor: e.g. [W m-2 nm-1] -> [W m-2]
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' with(sun.data, integrate_irradianceR(w.length, s.e.irrad))

integrate_irradianceR <-
function(w.length, s.irrad){
    irrad <- 0.0
    for (i in 1:(length(w.length)-1)){
      irrad <- irrad + 
        ((s.irrad[i] + s.irrad[i+1]) / 2 * (w.length[i+1] - w.length[i]))
    }
    return(irrad)
}
