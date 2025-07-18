% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mspct.sd.r
\name{s_sd}
\alias{s_sd}
\alias{s_sd.default}
\alias{s_sd.generic_spct}
\alias{s_sd.filter_mspct}
\alias{s_sd.source_mspct}
\alias{s_sd.response_mspct}
\alias{s_sd.reflector_mspct}
\alias{s_sd.calibration_mspct}
\alias{s_sd.cps_mspct}
\alias{s_sd.raw_mspct}
\title{Standard Deviation of a collection of spectra}
\usage{
s_sd(x, na.rm, ...)

\method{s_sd}{default}(x, na.rm = FALSE, ...)

\method{s_sd}{generic_spct}(x, na.rm = FALSE, ...)

\method{s_sd}{filter_mspct}(x, na.rm = FALSE, ...)

\method{s_sd}{source_mspct}(x, na.rm = FALSE, ...)

\method{s_sd}{response_mspct}(x, na.rm = FALSE, ...)

\method{s_sd}{reflector_mspct}(x, na.rm = FALSE, ...)

\method{s_sd}{calibration_mspct}(x, na.rm = FALSE, ...)

\method{s_sd}{cps_mspct}(x, na.rm = FALSE, ...)

\method{s_sd}{raw_mspct}(x, na.rm = FALSE, ...)
}
\arguments{
\item{x}{An R object.}

\item{na.rm}{logical A value indicating whether NA values should be stripped
before the computation proceeds.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
If \code{x} is a collection spectral of objects, such as a
  \code{"filter_mspct"} object, the returned object belongs to the same class
  as the members of the collection, such as \code{"filter_spct"}, containing
  the summary spectrum, with variables with names tagged for summaries other
  than mean or median.
}
\description{
Method to compute the "parallel" standard deviation of values across members
of a collection of spectra or of a spectral object containing multiple
spectra in long form.
}
\details{
Method specializations compute the standard deviation at each
  wavelength across a group of spectra stored in an object of one of the
  classes defined in package 'photobiology'. Omission of NAs is done
  separately at each wavelength. Interpolation is not applied, so all spectra
  in \code{x} must share the same set of wavelengths. An error is triggered
  if this condition is nor fulfilled.
}
\note{
Objects of classes \code{raw_spct} and \code{cps_spct} can contain data
  from multiple scans in multiple variables or "columns". The parallel
  summaries' methods accept as arguments objects of these classes only if
  spectra contain data for a single spectrometer scan. In the case of
  \code{cps_spct} objects, a single column can also contain data from
  multiple scans spliced into a single variable.
}
\examples{
s_sd(sun_evening.mspct)

}
\seealso{
See \code{\link[stats]{sd}} for details about \code{sd()} methods
  for other classes.
}
