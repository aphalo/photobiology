#' @title Molar spectral attenuation coefficient of water
#'
#' @description A dataset containing the wavelengths at a 2 nm interval and the
#'   corresponding attenuation coefficients.
#'
#' @details \itemize{ \item w.length (nm), range 300 to 800 nm. \item K.mole
#' (cm-1/M)}
#'
#' @author Buiteveld et al. (1994) (original data)
#'
#' @references H. Buiteveld and J. M. H. Hakvoort and M. Donze (1994) "The
#'   optical properties of pure water," in SPIE Proceedings on Ocean Optics XII,
#'   edited by J. S. Jaffe, 2258, 174--183.
#'
#' \url{https://omlc.org/spectra/water/}
#'
#' @docType data
#' @keywords datasets
#' @format A \code{solute_spct} object with 251 rows and 2 variables
#' @family Spectral data examples
#'
#' @examples
#' head(water.spct)
#' summary(water.spct)
#' solute_properties(water.spct)
#' cat(comment(water.spct))
#'
"water.spct"

#' @title Molar spectral attenuation coefficient of phenylalanine
#'
#' @description A dataset containing the wavelengths at a 0.25 nm interval
#'   and the corresponding attenuation coefficients.
#'
#' @details \itemize{ \item w.length (nm), range 222 to 720 nm. \item K.mole
#' (cm-1/M)}
#'
#' @author Du et ql. (original data); Scott Prahl (included data).
#'
#' @references
#' \url{https://omlc.org/spectra/PhotochemCAD/html/073.html}
#'
#' H. Du, R. A. Fuh, J. Li, A. Corkan, J. S. Lindsey, "PhotochemCAD: A
#' computer-aided design and research tool in photochemistry," Photochem.
#' Photobiol., 68, 141-142, 1998.
#'
#' J. M. Dixon, M. Taniguchi and J. S. Lindsey "PhotochemCAD 2. A refined
#' program with accompanying spectral databases for photochemical calculations",
#' Photochem. Photobiol., 81, 212-213, 2005.
#'
#' @docType data
#' @keywords datasets
#' @format A \code{solute_spct} object with 1993 rows and 2 variables
#' @family Spectral data examples
#'
#' @examples
#' head(phenylalanine.spct)
#' summary(phenylalanine.spct)
#' solute_properties(phenylalanine.spct)
#' cat(comment(phenylalanine.spct))
#'
"phenylalanine.spct"
