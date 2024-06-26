% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zmspct.classes.R
\name{is.generic_mspct}
\alias{is.generic_mspct}
\alias{is.calibration_mspct}
\alias{is.raw_mspct}
\alias{is.cps_mspct}
\alias{is.source_mspct}
\alias{is.response_mspct}
\alias{is.filter_mspct}
\alias{is.reflector_mspct}
\alias{is.object_mspct}
\alias{is.solute_mspct}
\alias{is.chroma_mspct}
\alias{is.any_mspct}
\title{Query class of spectrum objects}
\usage{
is.generic_mspct(x)

is.calibration_mspct(x)

is.raw_mspct(x)

is.cps_mspct(x)

is.source_mspct(x)

is.response_mspct(x)

is.filter_mspct(x)

is.reflector_mspct(x)

is.object_mspct(x)

is.solute_mspct(x)

is.chroma_mspct(x)

is.any_mspct(x)
}
\arguments{
\item{x}{an R object.}
}
\value{
These functions return \code{TRUE} if its argument is a of the queried type
  of spectrum and \code{FALSE} otherwise.
}
\description{
Functions to check if an object is of a given type of spectrum, or coerce it if
possible.
}
\note{
Derived types also return TRUE for a query for a base type such as
\code{generic_mspct}.
}
\examples{
my.mspct <- filter_mspct(list(polyester.spct, yellow_gel.spct))
is.any_mspct(my.mspct)
is.filter_mspct(my.mspct)
is.source_mspct(my.mspct)

}
