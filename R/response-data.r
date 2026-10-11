#' @title Spectral response of a GaAsP photodiode
#'
#' @description A dataset containing wavelengths at a 1 nm interval and
#'   spectral response as \eqn{A / (W / nm)} for GaAsP photodiode type
#'   G6262 from Hamamatsu. Data digitized from manufacturer's data sheet.
#'   The value at the peak is 0.19 \eqn{A / W}.
#'
#' @details \itemize{ \item w.length (nm). \item s.e.response
#' (A/W)  }
#'
#' @references
#' Hamamatsu (2011) Datasheet: GaAsP Photodiodes G5645 G5842 G6262. Hamamatsu
#' Photonics KK, Hamamatsu, City.
#' http://www.hamamatsu.com/jp/en/G6262.html.
#' Visited 2017-12-15.
#'
#' @docType data
#' @keywords datasets
#' @format A \code{response_spct} object with 94 rows and 2 variables
#' @family Spectral data examples
#'
#' @examples photodiode.spct
#'
"photodiode.spct"

#' @title Spectral response of a back-thinned CCD image sensor.
#'
#' @description A dataset containing wavelengths at a 1 nm interval and
#'   spectral response as quantum efficiency for CCD sensor type
#'   S11071/S10420 from Hamamatsu (measured without a quartz window). These
#'   vectors are frequently used as sensors in high-UV-sensitivity vector
#'   spectrometers. Data digitized from manufacturer's data sheet.
#'   The original data is expressed as percent quantum efficiency with a value
#'   of 77\% at the peak. The data have been re-expressed as fractions of one.
#'
#' @details \itemize{ \item w.length (nm). \item s.q.response
#' (fractional quantum efficiency)  }
#'
#' @references
#' Hamamatsu (2014) Datasheet: CCD Image Sensors S11071/S10420-01 Series.
#' Hamamatsu Photonics KK, Hamamatsu, City.
#' http://www.hamamatsu.com/jp/en/S11071-1004.html.
#' Visited 2017-12-15.
#'
#' @docType data
#' @keywords datasets
#' @format A \code{response_spct} object with 186 rows and 2 variables
#' @family Spectral data examples
#'
#' @examples
#' ccd.spct
#'
"ccd.spct"

#' @title Spectral response of two light sensors.
#'
#' @description A dataset containing a collection of two spectra.
#'
#' @details The spectra in \code{\link{photodiode.spct}} and
#'   \code{\link{ccd.spct}} stored as a collection in a
#'   \code{\link{response_mspct}} object named \code{response.mspct} with
#'   members \code{photodiode} and \code{ccd}, and and in long form in a
#'   \code{link{response_spct}} object named \code{response.mspct} identified
#'   bit the levels of factor \code{spct.idx}.
#'
#' @seealso \code{\link{photodiode.spct}} and \code{\link{ccd.spct}}.
#'
#' @docType data
#' @keywords datasets
#' @format A \code{response_spct} object with 186 rows and 2 variables
#' @family Spectral data examples
#'
#' @examples
#' two_sensors.mspct
#' two_sensors.spct
#'
"two_sensors.mspct"

#' @rdname two_sensors.mspct
#'
"two_sensors.spct"
