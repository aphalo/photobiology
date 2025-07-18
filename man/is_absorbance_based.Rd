% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spct.classes.r
\name{is_absorbance_based}
\alias{is_absorbance_based}
\alias{is_absorptance_based}
\alias{is_transmittance_based}
\title{Query if a spectrum contains absorbance or transmittance data}
\usage{
is_absorbance_based(x)

is_absorptance_based(x)

is_transmittance_based(x)
}
\arguments{
\item{x}{an R object}
}
\value{
\code{is_absorbance_based} returns a logical value, \code{TRUE}
  if its argument is a \code{filter_spct}
  object that contains spectral absorbance data and \code{FALSE} otherwise,
  but returns \code{NA} for any other R object, including those belonging
  other \code{generic_spct}-derived classes.

\code{is_absorptance_based} returns a logical value, if its argument
  is a \code{filter_spct} object, \code{TRUE} if it contains data as spectral
  absorptance and \code{FALSE} otherwise, but returns \code{NA} for any other
  R object, including those belonging other \code{generic_spct}-derived
  classes.

\code{is_transmittance_based} returns TRUE if its argument is a
  \code{filter_spct} object that contains spectral transmittance data and
  FALSE if it does not contain such data, but returns NA for any other R
  object, including those belonging other \code{generic_spct}-derived
  classes.
}
\description{
Functions to query if an filter spectrum contains spectral absorbance data or
spectral transmittance data.
}
\examples{
is_absorbance_based(polyester.spct)
my.spct <- T2A(polyester.spct)
is.filter_spct(my.spct)
is_absorbance_based(my.spct)

is_absorptance_based(polyester.spct)

is_transmittance_based(polyester.spct)

}
\seealso{
Other query units functions: 
\code{\link{is_mole_based}()},
\code{\link{is_photon_based}()}
}
\concept{query units functions}
