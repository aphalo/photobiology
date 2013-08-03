#' Quantification of ultraviolet and visible radiation for photobiology
#' 
#' Calculation of quantities relevant to the effects of radiation on different 
#' organisms and biological processes from spectral data. The package is designed 
#' so that it is easy for the user to create new quantification functions.
#' 
#' @docType package
#' @keywords misc
#' @name photobiology-package
#' @author Pedro J. Apahalo
#' @details
#' \tabular{ll}{
#' Package: \tab photobiology\cr
#' Type: \tab Package\cr
#' Version: \tab 0.2.1\cr
#' Date: \tab 2013-08-03\cr
#' License: \tab GPL (>2.0)\cr
#' }
#' The most important functions in the package are \code{\link{energy_irradiance}}, 
#' \code{\link{photon_irradiance}}, and \code{\link{new_waveband}}. The first two, 
#' are used to obtain energy and photon irradiances from spectral data. The third 
#' function is used to define how to calculate new quantities.
#' @references
#' Aphalo, P. J., Albert, A., Björn, L. O., McLeod, A. R., Robson, T. M., 
#' Rosenqvist, E. (Eds.). (2012). Beyond the Visible: A handbook of best 
#' practice in plant UV photobiology (1st ed., p. xxx + 174). 
#' Helsinki: University of Helsinki, Department of Biosciences, 
#' Division of Plant Biology. ISBN 978-952-10-8363-1 (PDF), 
#' 978-952-10-8362-4 (paperback). Open access PDF download available at 
#' http://hdl.handle.net/10138/37558
#' @note When released, this package will replace the package 
#' \code{\link[UVcalc:UVcalc-package]{UVcalc}}. 
#' @examples
#' data(sun.data)
#' with(sun.data, photon_irradiance(w.length, s.e.irrad)) # the whole spectrum
#' with(sun.data, photon_irradiance(w.length, s.e.irrad, new_waveband(400,700))) # idem
#' 
#' with(sun.data, energy_irradiance(w.length, s.e.irrad)) # the whole spectrum
#' with(sun.data, energy_irradiance(w.length, s.e.irrad, new_waveband(400,700))) # idem
NULL