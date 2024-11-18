#' @details Package \sQuote{photobiology} is at the core of a suite of R packages
#'   supporting computations and plotting relevant to photobiology (described at
#'   \url{https://www.r4photobiology.info/}). Package 'photobiology' has its
#'   main focus in the characterization of the light environment, the
#'   description of optical properties of objects and substances and description
#'   of light responses of organisms and devices used to measure light. The
#'   facilities for spectral data storage and manipulations are widely useful in
#'   photobiology, chemistry, geophysics, radiation climatology and remote
#'   sensing. Astronomical computations for the sun are also implemented. The
#'   design of object classes for spectral data supports reproducibility by
#'   facilitating the consistent use of units and physical quantities and
#'   consistent embedding of metadata. Data are expressed throughout using SI
#'   base units, except for wavelengths which are consistently expressed in
#'   nanometres [\eqn{nm}]. Please see the vignette \emph{0: The R for
#'   Photobiology Suite} for a description of the suite.
#'
#' @references
#' Aphalo, P. J., Albert, A., Björn, L. O., McLeod, A. R., Robson, T. M.,
#' Rosenqvist, E. (Eds.). (2012). \emph{Beyond the Visible: A handbook of best
#' practice in plant UV photobiology} (1st ed., p. xx + 174). Helsinki:
#' University of Helsinki, Department of Biosciences, Division of Plant Biology.
#' ISBN 978-952-10-8363-1 (PDF), 978-952-10-8362-4 (paperback). Open access PDF
#' download available at \doi{10.31885/9789521083631}.
#'
#' Aphalo, Pedro J. (2015) The r4photobiology suite. \emph{UV4Plants Bulletin}, 2015:1,
#' 21-29. \doi{10.19232/uv4pb.2015.1.14}.
#'
#' Maia, R., Eliason, C. M., Bitton, P. P., Doucet, S. M., Shawkey, M. D. (2013)
#' pavo: an R package for the analysis, visualization and organization of
#' spectral data. \emph{Methods in Ecology and Evolution}, 4(10):906-913.
#' \doi{10.1111/2041-210X.12069}.
#'
#' @section Acknowledgements: This work was funded by the Academy of Finland
#'   (decision 252548). COST Action FA9604 \sQuote{UV4Growth} facilitated discussions
#'   and exchanges of ideas that lead to the development of this package. The
#'   contributions of Andy McLeod, Lars Olof Björn, Nigel Paul, Lasse Ylianttila,
#'   T. Matthew Robson and Titta Kotilainen were specially significant.
#'   Tutorials by Hadley Wickham and comments on my presentation at UseR!2015
#'   allowed me to significantly improve the coding and functionality.
#'
#' @importFrom rlang .data
#' @import tibble
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
#' sunrise_time(lubridate::today(tzone = "Europe/Helsinki"), tz = "Europe/Helsinki",
#'              geocode = data.frame(lat = 60, lon = 25),
#'              unit.out = "hour")
#' day_length(lubridate::today(tzone = "Europe/Helsinki"), tz = "Europe/Helsinki",
#'            geocode = data.frame(lat = 60, lon = 25),
#'            unit.out = "hour")
#' # colour as seen by humans
#' color_of(sun.spct)
#' color_of(sun.spct * yellow_gel.spct)
#' # filter transmittance
#' transmittance(yellow_gel.spct)
#' transmittance(yellow_gel.spct, waveband(c(400,500)))
#' transmittance(yellow_gel.spct, waveband(c(500,700)))
"_PACKAGE"
