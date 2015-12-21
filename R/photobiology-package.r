#' Base Photobiological Calculations
#'
#' Calculation of quantities relevant to the effects of radiation on different
#' organisms and biological processes from spectral data. The package is
#' designed so that it is easy for the user to create new quantification
#' functions.
#'
#' @docType package
#' @keywords misc
#' @name photobiology-package
#' @aliases photobiology
#' @author Pedro J. Aphalo
#' @details
#' \tabular{ll}{
#' Package: \tab photobiology\cr
#' Type: \tab Package\cr
#' Version: \tab 0.9.1\cr
#' Date: \tab 2015-12-21\cr
#' License: \tab GPL (>= 3.0)\cr
#' URL: \tab \url{http://www.r4photobiology.info},\cr
#' \tab \url{https://bitbucket.org/aphalo/photobiology}\cr
#' BugReports: \tab \url{https://bitbucket.org/aphalo/photobiology}\cr
#' }
#' This package is the core of a suite of packages for photobiological
#' data analysis and plotting. The accompanying packages are data and
#' definitions that are to a large extent application-area specific
#' while the functions in this package are widely useful in photobiology
#' and radiation quantification in geophysics and meteorology.
#'
#' @references
#' Aphalo, P. J., Albert, A., Björn, L. O., McLeod, A. R., Robson, T. M.,
#' Rosenqvist, E. (Eds.). (2012). Beyond the Visible: A handbook of best
#' practice in plant UV photobiology (1st ed., p. xxx + 174).
#' Helsinki: University of Helsinki, Department of Biosciences,
#' Division of Plant Biology. ISBN 978-952-10-8363-1 (PDF),
#' 978-952-10-8362-4 (paperback). Open access PDF download available at
#' http://hdl.handle.net/10138/37558
#'
#' @note This package is still under development, but current functionality
#'  is mostly stable.
#'
#' @importFrom Rcpp evalCpp
#'
#' @examples
#' # irradiance of the whole spectrum
#' irrad(sun.spct)
#' # photon irradiance 400 nm to 700 nm
#' q_irrad(sun.spct, waveband(c(400,700)))
#' # energy irradiance 400 nm to 700 nm
#' e_irrad(sun.spct, waveband(c(400,700)))
#' # simulating the effect of a filter on solar irradiance
#' e_irrad(sun.spct * yellow_gel.spct, waveband(c(400,500)))
#' e_irrad(sun.spct * yellow_gel.spct, waveband(c(500,700)))
#' # daylength
#' sunrise_time(lubridate::today(tzone = "EET"), tz = "EET",
#'              lat = 60, lon = 25, unit.out = "hour")
#' day_length(lubridate::today(tzone = "EET"), tz = "EET",
#'               lat = 60, lon = 25, unit.out = "hour")
#' # colour as seen by humans
#' color(sun.spct)
#' color(sun.spct * yellow_gel.spct)
#' # filter transmittance
#' transmittance(yellow_gel.spct)
#' transmittance(yellow_gel.spct, waveband(c(400,500)))
#' transmittance(yellow_gel.spct, waveband(c(500,700)))
NULL
