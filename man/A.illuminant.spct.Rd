% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cieillum.data.r
\docType{data}
\name{A.illuminant.spct}
\alias{A.illuminant.spct}
\title{CIE A illuminant data}
\format{
A source spectrum with 96 rows and 2 variables
}
\usage{
A.illuminant.spct
}
\description{
A dataset containing wavelengths at a 5 nm interval (300 nm to
  830 nm) and the corresponding spectral energy irradiance normalized to 1 at
  560 nm. Spectrum approximates typical, domestic, tungsten-filament lighting
  and 'corresponds' to a black body a 2856 K. CIE standard illuminant A is
  intended to represent typical, domestic, tungsten-filament lighting.
  Original data from \url{http://files.cie.co.at/204.xls} downloaded on
  2014-07-25 The variables are as follows:
}
\details{
\itemize{ \item w.length (nm) \item s.e.irrad (rel. units) }
}
\examples{
A.illuminant.spct

}
\seealso{
Other Spectral data examples: 
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
\code{\link{sun_daily.spct}},
\code{\link{sun_evening.spct}},
\code{\link{two_filters.spct}},
\code{\link{water.spct}},
\code{\link{white_led.source_spct}}
}
\author{
CIE
}
\concept{Spectral data examples}
\keyword{datasets}
