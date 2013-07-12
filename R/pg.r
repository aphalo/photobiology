#' Definition of PG weighted waveband
#' 
#' Plant growth BSWF
#' 
#' @param norm normalization wavelength (nm)
#' @usage PG(norm=NULL)
#' @references
#' Flint, S. and Caldwell M. M. (2003)
#' 
#' @export
#' @seealso \code{\link{photon_irradiance}} and \code{\link{energy_irradiance}}
#' @examples
#' PG()
#' PG(300)

PG <- function(norm=NULL) {
  new_waveband(w.low=250,w.high=390,weight="SWF",SWF.fun=PG.q.fun,norm=norm)
}