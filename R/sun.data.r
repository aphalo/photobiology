#' @title Solar spectral irradiance (simulated)
#'
#' @description A dataset containing the wavelengths at a 1 nm interval and the
#'   corresponding spectral (energy) irradiance and spectral photon irradiance.
#'   Values simulated for 22 June 2010, near midday, at Helsinki, under partly
#'   cloudy conditions. The variables are as follows:
#'
#' @details \itemize{ \item w.length (nm), range 293 to 800 nm. \item s.e.irrad
#' (W m-2 nm-1) \item s.q.irrad (mol m-2 nm-1) }
#'
#' @author Anders K. Lindfors (data)
#' @references Lindfors, A.; Heikkilä, A.; Kaurola, J.; Koskela, T. & Lakkala,
#' K. (2009) Reconstruction of Solar Spectral Surface UV Irradiances Using
#' Radiative Transfer Simulations. Photochemistry and Photobiology, 85:
#' 1233-1239
#'
#' @docType data
#' @keywords datasets
#' @format A \code{source_spct} object and a \code{data.frame}, each with 511
#' rows and 3 variables
#' @family Spectral data examples
#'
#' @examples
#' sun.spct
#' summary(sun.spct)
#'
"sun.spct"

#' @rdname sun.spct
#'
"sun.data"

#' Daily solar spectral irradiance (simulated)
#'
#' A dataset containing the wavelengths at a 1 nm interval and the corresponding
#' spectral (energy) irradiance. Values simulated for 2 June 2012, at Helsinki,
#' under clear sky conditions. The variables are as follows:
#'
#' \itemize{ \item w.length (nm), range 290 to 800 nm. \item s.e.irrad (J d-1
#' m-2 nm-1) \item s.q.irrad (mol d-1 m-2 nm-1) }
#'
#' @author Anders K. Lindfors (data)
#' @references Lindfors, A.; Heikkilä, A.; Kaurola, J.; Koskela, T. & Lakkala,
#' K. (2009) Reconstruction of Solar Spectral Surface UV Irradiances Using
#' Radiative Transfer Simulations. Photochemistry and Photobiology, 85:
#' 1233-1239
#'
#' @note The simulations are based on libRadTran using hourly mean global
#' radiation measurements to estimate cloud cover. The simulations were for
#' each hour and the results integrated for the whole day.
#'
#' @section Deprecation!:
#' Objects \code{sun.daily.spct} and \code{sun.daily.data} have been renamed
#' into \code{sun_daily.spct} and \code{sun_daily.data}, for consistency with
#' other data sets in the package. Please, use the new names for new code.
#'
#' @docType data
#' @keywords datasets
#' @format A \code{source_spct} object and a \code{data.frame}, each with 511
#' rows and 3 variables
#' @family Spectral data examples
#'
#' @examples
#' sun.daily.spct
#' summary(sun.daily.spct)
#'
"sun_daily.spct"

#' @rdname sun_daily.spct
#'
"sun_daily.data"

#' @rdname sun_daily.spct
#'
"sun.daily.spct"

#' @rdname sun_daily.spct
#'
"sun.daily.data"

#' Time series of solar spectral irradiance (measured)
#'
#' Two data objects containing containing the same time series of five spectra.
#' Values measured in Viikki, Helsinki, under nearly clear sky in a summer
#' evening.
#'
#' The variables are as follows:
#'
#' \itemize{ \item w.length (nm), range 290 to 1000 nm. \item s.e.irrad (J d-1
#' m-2 nm-1) \item s.q.irrad (mol d-1 m-2 nm-1) }
#'
#' @author Pedro J. Aphalo (data)
#'
#' @docType data
#' @keywords datasets
#' @format A \code{source_spct} object and a \code{source_mspct} object.
#' @family Spectral data examples
#'
#' @examples
#' summary(sun_evening.mspct)
#' colnames(sun_evening.spct)
#'
"sun_evening.spct"

#' @rdname sun_evening.spct
#'
"sun_evening.mspct"
