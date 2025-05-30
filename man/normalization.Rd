% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/waveband.class.r
\name{normalization}
\alias{normalization}
\alias{normalization.default}
\alias{normalization.waveband}
\alias{normalization.generic_spct}
\alias{normalization.summary_generic_spct}
\alias{normalization.generic_mspct}
\title{Normalization of an R object}
\usage{
normalization(x)

\method{normalization}{default}(x)

\method{normalization}{waveband}(x)

\method{normalization}{generic_spct}(x)

\method{normalization}{summary_generic_spct}(x)

\method{normalization}{generic_mspct}(x)
}
\arguments{
\item{x}{an R object}
}
\value{
A single numeric value of wavelength [\eqn{nm}] or a list with
  with members.
}
\description{
Normalization wavelength [\eqn{nm}] and other normalization metadata of an R
object, retrieved from the object's attributes.
}
\details{
In the case of wavebands for spectral weighting functions
  (\code{waveband} objects), the normalization wavelength is returned. For
  spectral objects (\code{generic_spct} and derived ), the normalization
  descriptor, a list object, is returned. This list contains in addition to
  the normalization wavelength, the multiplier used and type of normalization
  applied. These metadata makes it possible to "undo" the normalization and
  to "update" the normalization after a transformation, such as conversion to
  a related physical quantity, of the spectral data.
}
\section{Methods (by class)}{
\itemize{
\item \code{normalization(default)}: Default methods.

\item \code{normalization(waveband)}: Normalization of a \code{\link{waveband}} object.

\item \code{normalization(generic_spct)}: Normalization of a \code{\link{generic_spct}} object.

\item \code{normalization(summary_generic_spct)}: Normalization of a \code{\link{summary.generic_spct}} object.

\item \code{normalization(generic_mspct)}: Normalization of a \code{\link{generic_mspct}} object.

}}
\note{
Older versions of the package stored only a subset of the metadata or
  only a flag to indicate that normalization had been applied. For such
  objects some or even all fields in the returned list are set to \code{NA}.
}
\examples{

is_normalized(sun.spct)
normalization(sun.spct)
sun_norm.spct <- normalize(sun.spct)
is_normalized(sun_norm.spct)
normalization(sun_norm.spct)

my_wband <- waveband(c(400,700))
is_normalized(my_wband)
normalization(my_wband)

}
\seealso{
Other waveband attributes: 
\code{\link{is_effective}()},
\code{\link{labels}()}
}
\concept{waveband attributes}
