#' Definition of PAR waveband
#' 
#' Photosythetically active radiation (400...700 nm), no weighting
#' applied.
#' 
#' @usage PAR()
#' @export
#' @seealso \code{\link{photon_irradiance}} and \code{\link{energy_irradiance}}
#' @examples
#' PAR()
PAR <- function(){
  return(new_waveband(400,700,hinges=c(399.99,400,699.99,700)))
}