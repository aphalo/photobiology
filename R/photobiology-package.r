#' Quantification of ultraviolet and visible radiation for photobiology
#'
#' Calculation of quantities relevant to the effects of radiation on different
#' organisms and biological processes from spectral data. The package is designed
#' so that it is easy for the user to create new quantification functions.
#'
#' @docType package
#' @keywords misc
#' @name photobiology-package
#' @author Pedro J. Aphalo
#' @details
#' \tabular{ll}{
#' Package: \tab photobiology\cr
#' Type: \tab Package\cr
#' Version: \tab 0.5.9\cr
#' Date: \tab 2015-01-30\cr
#' License: \tab GPL (>= 2.0)\cr
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
#' @note This package is still under development, but is by now quite stable.
#'
#' @import lubridate methods
#'
#' @importFrom Rcpp evalCpp
#' @importFrom data.table data.table tables setkey setkeyv key "key<-" haskey CJ SJ copy
#' @importFrom data.table as.data.table is.data.table test.data.table last like "%like%" between "%between%"
#' @importFrom data.table truelength alloc.col ":="
#' @importFrom data.table setattr setnames setcolorder set setDT setDF
#' @importFrom data.table setorder setorderv
#' @importFrom data.table setNumericRounding getNumericRounding
#' @importFrom data.table chmatch "%chin%" chorder chgroup
#' @importFrom data.table fread
#' @importFrom data.table address
#' @importFrom data.table .SD .N .I .GRP .BY
#' @importClassesFrom data.table data.table
#'
#' @examples
#' data(sun.data)
#' with(sun.data, photon_irradiance(w.length, s.e.irrad)) # the whole spectrum
#' with(sun.data, photon_irradiance(w.length, s.e.irrad, new_waveband(400,700)))
#'
#' with(sun.data, energy_irradiance(w.length, s.e.irrad)) # the whole spectrum
#' with(sun.data, energy_irradiance(w.length, s.e.irrad, new_waveband(400,700)))
#'
NULL
