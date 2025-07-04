% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/join.mspct.R
\name{join_mspct}
\alias{join_mspct}
\alias{join_mspct.default}
\alias{join_mspct.generic_mspct}
\alias{join_mspct.source_mspct}
\alias{join_mspct.response_mspct}
\alias{join_mspct.filter_mspct}
\alias{join_mspct.reflector_mspct}
\alias{join_mspct.object_mspct}
\alias{join_mspct.solute_mspct}
\title{Join all spectra in a collection}
\usage{
join_mspct(x, type, ...)

\method{join_mspct}{default}(x, type = "full", ...)

\method{join_mspct}{generic_mspct}(x, type = "full", col.name, validate.names = TRUE, ...)

\method{join_mspct}{source_mspct}(x, type = "full", unit.out = "energy", validate.names = TRUE, ...)

\method{join_mspct}{response_mspct}(x, type = "full", unit.out = "energy", validate.names = TRUE, ...)

\method{join_mspct}{filter_mspct}(
  x,
  type = "full",
  qty.out = "transmittance",
  validate.names = TRUE,
  ...
)

\method{join_mspct}{reflector_mspct}(x, type = "full", validate.names = TRUE, ...)

\method{join_mspct}{object_mspct}(x, type = "full", qty.out, validate.names = TRUE, ...)

\method{join_mspct}{solute_mspct}(x, type = "full", validate.names = TRUE, ...)
}
\arguments{
\item{x}{generic_mspct object, or an object of a class derived from
\code{generic_mspct}.}

\item{type}{character Type of join: \code{"inner"} (default) or
\code{"full"}. See details for more information.}

\item{...}{ignored (possibly used by derived methods).}

\item{col.name}{character, name of the column in the spectra to be preserved,
in addition to "w.length".}

\item{validate.names}{logical A flag to enable (default) or disable
validation of column names with \code{\link[base]{make.names}}.}

\item{unit.out}{character Allowed values \code{"energy"}, and
\code{"photon"}, or its alias \code{"quantum"}.}

\item{qty.out}{character Allowed values \code{"transmittance"},
\code{"absorptance"}, and \code{"absorbance"} and in the method for
\code{object_spct}, also \code{"reflectance"} (.}
}
\value{
A \code{data.frame} with the spectra joined by, possibly
  interpolated, wavelength, with rows sorted by wavelength (variable
  \code{w.length}) and data columns named according to the names of members
  in \code{x}, by default made unique and valid.
}
\description{
Join all the spectra contained in a homogeneous collection, returning a data
frame with spectral-data columns named according to the names of the spectra
in the collection. By default a full join is done within the overlapping range
of wavelengths, after interpolating the spectra to a shared set of wavelength
values, and discarding data for wavelength not shared. Alternatively, filling
the spectral data for wavelengths outside the overlapping range with with
\code{NA} when data is not available.
}
\section{Methods (by class)}{
\itemize{
\item \code{join_mspct(default)}: 

\item \code{join_mspct(generic_mspct)}: 

\item \code{join_mspct(source_mspct)}: 

\item \code{join_mspct(response_mspct)}: 

\item \code{join_mspct(filter_mspct)}: 

\item \code{join_mspct(reflector_mspct)}: 

\item \code{join_mspct(object_mspct)}: 

\item \code{join_mspct(solute_mspct)}: 

}}
\note{
Currently only \code{generic_spct}, \code{source_mspct},
  \code{response_mspct}, \code{filter_mspct}, \code{reflector_mspct},
  \code{object_mspct} and \code{solute_mspct} classes have this method
  implemented.
}
\examples{
my.mspct <- solute_mspct(list(water = water.spct, pha = phenylalanine.spct))
join_mspct(my.mspct, type = "inner")
join_mspct(my.mspct, type = "full")

}
\concept{conversion of collections of spectra}
