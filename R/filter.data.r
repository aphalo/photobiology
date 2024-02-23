#' @title Transmittance spectrum of plastic films
#'
#' @description Datasets containing the wavelengths at a 1 nm interval and
#'   fractional total transmittance for a clear polyester film and a yellow
#'   theatrical "gel".
#'
#' @details \itemize{ \item \code{w.length} (nm). \item \code{Tfr}
#' (0..1). \item \code{spct.idx} (names, only in \code{two_filters.spct}).}
#'
#' @docType data
#' @keywords datasets
#' @format A \code{filter_spct} object with 611 rows and 2 variables.
#'   Individually as \code{filter_spct} objects, and together as a collection
#'   stored in a \code{filter_mspct} object and in a long-form
#'   \code{filter_spct} object.
#' @family Spectral data examples
#'
#' @note Package 'photobiologyFilters' contains data sets for hundreds of
#'   optical filters and materials in objects of these same classes, ready to be
#'   used with package 'photobiology'.
#'
#' @examples
#' polyester.spct
#' yellow.spct
#' summary(two_filters.mspct)
#'
"two_filters.spct"

#' @rdname two_filters.spct
#'
"two_filters.mspct"

#' @rdname two_filters.spct
#'
"polyester.spct"

#' @rdname two_filters.spct
#'
"yellow_gel.spct"

#' @title Theoretical spectrum of clear and apaque materials
#'
#' @description Dataset for hypothetical objects with transmittance 1/1
#'   (100\%) and transmittance 0/1 (0\%)
#'
#' @details \itemize{ \item w.length (nm). \item Tfr
#' (0..1)  }
#'
#' @docType data
#' @keywords datasets
#' @format A \code{filter_spct} object with 4 rows and 2 variables
#' @family Spectral data examples
#'
#' @examples
#' clear.spct
#' opaque.spct
#'
"clear.spct"

#' @rdname clear.spct
#'
"opaque.spct"
