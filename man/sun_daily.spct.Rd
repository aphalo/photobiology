% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sun.data.r
\docType{data}
\name{sun_daily.spct}
\alias{sun_daily.spct}
\alias{sun_daily.data}
\alias{sun.daily.spct}
\alias{sun.daily.data}
\title{Daily solar spectral irradiance (simulated)}
\format{
A \code{source_spct} object and a \code{data.frame}, each with 511
rows and 3 variables

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 511 rows and 3 columns.

An object of class \code{source_spct} (inherits from \code{generic_spct}, \code{tbl_df}, \code{tbl}, \code{data.frame}) with 522 rows and 3 columns.

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 511 rows and 3 columns.
}
\usage{
sun_daily.spct

sun_daily.data

sun.daily.spct

sun.daily.data
}
\description{
A dataset containing the wavelengths at a 1 nm interval and the corresponding
spectral (energy) irradiance. Values simulated for 2 June 2012, at Helsinki,
under clear sky conditions. The variables are as follows:
}
\details{
\itemize{ \item w.length (nm), range 290 to 800 nm. \item s.e.irrad (J d-1
m-2 nm-1) \item s.q.irrad (mol d-1 m-2 nm-1) }
}
\note{
The simulations are based on libRadTran using hourly mean global
radiation measurements to estimate cloud cover. The simulations were for
each hour and the results integrated for the whole day.
}
\section{Deprecation!}{

Objects \code{sun.daily.spct} and \code{sun.daily.data} have been renamed
into \code{sun_daily.spct} and \code{sun_daily.data}, for consistency with
other data sets in the package. Please, use the new names for new code.
}

\examples{
sun.daily.spct
summary(sun.daily.spct)

}
\references{
Lindfors, A.; Heikkilä, A.; Kaurola, J.; Koskela, T. & Lakkala,
K. (2009) Reconstruction of Solar Spectral Surface UV Irradiances Using
Radiative Transfer Simulations. Photochemistry and Photobiology, 85:
1233-1239
}
\seealso{
Other Spectral data examples: 
\code{\link{A.illuminant.spct}},
\code{\link{D50.illuminant.spct}},
\code{\link{D65.illuminant.spct}},
\code{\link{Ler_leaf.spct}},
\code{\link{black_body.spct}},
\code{\link{ccd.spct}},
\code{\link{clear.spct}},
\code{\link{filter_cps.mspct}},
\code{\link{green_leaf.spct}},
\code{\link{phenylalanine.spct}},
\code{\link{photodiode.spct}},
\code{\link{sun.spct}},
\code{\link{sun_evening.spct}},
\code{\link{two_filters.spct}},
\code{\link{two_sensors.mspct}},
\code{\link{water.spct}},
\code{\link{white_led.source_spct}}
}
\author{
Anders K. Lindfors (data)
}
\concept{Spectral data examples}
\keyword{datasets}
