#' White led bulb spectrum
#'
#' @description Datasets containing wavelengths and the
#'   corresponding spectral irradiance data for an Osram warm white led lamp,
#'   and the corresponding raw instrument counts and counts per second data
#'   underlying them.
#'
#' @details
#' \itemize{ \item w.length (nm), range 250 to 900 nm. \item s.e.irrad
#' (W m-2 nm-1)}
#'
#' or
#'
#' \itemize{ \item w.length (nm), range 188 to 1117 nm. \item cps }
#'
#' or
#'
#' \itemize{ \item w.length (nm), range 188 to 1117 nm.
#'           \item counts_1 \item counts_2 \item counts_3 }
#'
#' @docType data
#' @keywords datasets
#' @format A \code{source_spct} object with 1421 rows and 2 columns,
#'   a \code{cps_spct} object with 2068 rows and 2 columns, and
#'   a \code{raw_spct} object with 2068 rows and 4 columns.
#' @family Spectral data examples
#'
#' @examples
#' white_led.source_spct
#'
"white_led.source_spct"

#' @rdname white_led.source_spct
#'
"white_led.cps_spct"

#' @rdname white_led.source_spct
#'
#'
"white_led.raw_spct"
