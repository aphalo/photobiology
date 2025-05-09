% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zmspct.classes.R
\name{subset2mspct}
\alias{subset2mspct}
\title{Convert 'long' or tidy spectral data into a collection of spectra}
\usage{
subset2mspct(
  x,
  member.class = NULL,
  idx.var = NULL,
  drop.idx = TRUE,
  ncol = 1,
  byrow = FALSE,
  ...
)
}
\arguments{
\item{x}{a generic_spct object or of a derived class, or a data frame, or a
generic_mspct object or of a derived class.}

\item{member.class}{character string.}

\item{idx.var}{character Name of column containing data to be copied
unchanged to each spct object or used for member names. If \code{NULL},
the default, the name is retrieved from x or its members when possible.}

\item{drop.idx}{logical Flag indicating whether to drop or keep idx.var in
the collection members.}

\item{ncol}{integer Number of 'virtual' columns in data.}

\item{byrow}{logical If \code{ncol > 1} how to read in the data.}

\item{...}{additional named arguments passed to the member constructor
function.}
}
\value{
A collection of spectral objects, each with attributes set if x is a
  spectral object in long form with metadata attributes. If this object
  was created by row binding with 'photobiology' 0.9.14 or later then
  all metadata for each individual spectrum will be preserved, except for
  unique comments which are merged.
}
\description{
Convert a data frame object or spectral object into a collection
  of spectra object of the matching class. For data frames converting numeric
  columns other than wavelength into individual spct objects. For collection
  of spectra objects, subset/expand long-form members into multiple members
  of the same collection.
}
\note{
A non-null value for \code{member.class} is mandatory only when
  \code{x} is a data frame.
}
\seealso{
Other Coercion methods for collections of spectra: 
\code{\link{as.calibration_mspct}()},
\code{\link{as.chroma_mspct}()},
\code{\link{as.cps_mspct}()},
\code{\link{as.filter_mspct}()},
\code{\link{as.generic_mspct}()},
\code{\link{as.object_mspct}()},
\code{\link{as.raw_mspct}()},
\code{\link{as.reflector_mspct}()},
\code{\link{as.response_mspct}()},
\code{\link{as.solute_mspct}()},
\code{\link{as.source_mspct}()},
\code{\link{split2mspct}()}
}
\concept{Coercion methods for collections of spectra}
