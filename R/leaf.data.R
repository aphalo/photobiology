#' @title Green birch leaf reflectance.
#'
#' @description A dataset of spectral reflectance expressed as a fraction of one.
#'
#' @details \itemize{ \item w.length (nm) \item Rfr (0..1) }
#'
#' @docType data
#' @keywords datasets
#' @format A \code{reflector_spct} object with 226 rows and 2 variables
#' @family Spectral data examples
#' @references
#' Aphalo, P. J. & Lehto, T. Effects of light quality on growth and N
#' accumulation in birch seedlings Tree Physiology, 1997, 17, 125-132
#'
#' @examples
#' green_leaf.spct
#'
"green_leaf.spct"

#' @title Green Arabidopsis leaf reflectance and transmittance.
#'
#' @description A dataset of total spectral reflectance and total spectral
#'   transmittance expressed as fractions of one from the upper surface of a
#'   leaf of an Arabidopsis thaliana 'Ler' rosette.
#'
#' @details \itemize{ \item w.length (nm) \item Rfr (0..1) \item Tfr (0..1)}
#'
#' @note Measured with a Jaz spectrometer from Ocean Optics (USA) configured
#'   with a PX Xenon lamp module and Spectroclip double integrating spheres.
#'
#' @docType data
#' @keywords datasets
#' @format Datasets stored as \code{object_spct}, \code{reflector_spct} and
#'    \code{filter_spct} objects, containing transmittance and reflectance
#'    data.
#' @family Spectral data examples
#' @author
#' Aphalo, P. J. & Wang, F (unpublished data)
#'
#' @examples
#' Ler_leaf.spct
#' Ler_leaf_rflt.spct
#'
"Ler_leaf.spct"

#' @rdname Ler_leaf.spct
#'
"Ler_leaf_rflt.spct"

#' @rdname Ler_leaf.spct
#'
"Ler_leaf_trns.spct"

#' @rdname Ler_leaf.spct
#'
"Ler_leaf_trns_i.spct"
