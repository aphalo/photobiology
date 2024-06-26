% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spct.clean.r
\name{clean}
\alias{clean}
\alias{clean.default}
\alias{clean.source_spct}
\alias{clean.filter_spct}
\alias{clean.reflector_spct}
\alias{clean.solute_spct}
\alias{clean.object_spct}
\alias{clean.response_spct}
\alias{clean.cps_spct}
\alias{clean.raw_spct}
\alias{clean.generic_spct}
\alias{clean.source_mspct}
\alias{clean.filter_mspct}
\alias{clean.reflector_mspct}
\alias{clean.object_mspct}
\alias{clean.solute_mspct}
\alias{clean.response_mspct}
\alias{clean.cps_mspct}
\alias{clean.raw_mspct}
\alias{clean.generic_mspct}
\title{Clean (=replace) off-range values in a spectrum}
\usage{
clean(x, range, range.s.data, fill, ...)

\method{clean}{default}(x, range, range.s.data, fill, ...)

\method{clean}{source_spct}(
  x,
  range = x,
  range.s.data = c(0, NA),
  fill = range.s.data,
  unit.out = getOption("photobiology.radiation.unit", default = "energy"),
  ...
)

\method{clean}{filter_spct}(
  x,
  range = x,
  range.s.data = NULL,
  fill = range.s.data,
  qty.out = getOption("photobiology.filter.qty", default = "transmittance"),
  ...
)

\method{clean}{reflector_spct}(x, range = x, range.s.data = c(0, 1), fill = range.s.data, ...)

\method{clean}{solute_spct}(x, range = x, range.s.data = c(0, NA), fill = range.s.data, ...)

\method{clean}{object_spct}(
  x,
  range = x,
  range.s.data = c(0, 1),
  fill = range.s.data,
  min.Afr = NULL,
  ...
)

\method{clean}{response_spct}(
  x,
  range = x,
  range.s.data = c(0, NA),
  fill = range.s.data,
  unit.out = getOption("photobiology.radiation.unit", default = "energy"),
  ...
)

\method{clean}{cps_spct}(x, range = x, range.s.data = c(0, NA), fill = range.s.data, ...)

\method{clean}{raw_spct}(
  x,
  range = x,
  range.s.data = c(NA_real_, NA_real_),
  fill = range.s.data,
  ...
)

\method{clean}{generic_spct}(
  x,
  range = x,
  range.s.data = c(NA_real_, NA_real_),
  fill = range.s.data,
  col.names,
  ...
)

\method{clean}{source_mspct}(
  x,
  range = NULL,
  range.s.data = c(0, NA),
  fill = range.s.data,
  unit.out = getOption("photobiology.radiation.unit", default = "energy"),
  ...,
  .parallel = FALSE,
  .paropts = NULL
)

\method{clean}{filter_mspct}(
  x,
  range = NULL,
  range.s.data = NULL,
  fill = range.s.data,
  qty.out = getOption("photobiology.filter.qty", default = "transmittance"),
  ...,
  .parallel = FALSE,
  .paropts = NULL
)

\method{clean}{reflector_mspct}(
  x,
  range = NULL,
  range.s.data = c(0, 1),
  fill = range.s.data,
  ...,
  .parallel = FALSE,
  .paropts = NULL
)

\method{clean}{object_mspct}(
  x,
  range = NULL,
  range.s.data = c(0, 1),
  fill = range.s.data,
  min.Afr = NULL,
  ...,
  .parallel = FALSE,
  .paropts = NULL
)

\method{clean}{solute_mspct}(
  x,
  range = NULL,
  range.s.data = c(0, NA),
  fill = range.s.data,
  ...,
  .parallel = FALSE,
  .paropts = NULL
)

\method{clean}{response_mspct}(
  x,
  range = NULL,
  range.s.data = c(0, NA),
  fill = range.s.data,
  unit.out = getOption("photobiology.radiation.unit", default = "energy"),
  ...,
  .parallel = FALSE,
  .paropts = NULL
)

\method{clean}{cps_mspct}(
  x,
  range = NULL,
  range.s.data = c(0, NA),
  fill = range.s.data,
  ...,
  .parallel = FALSE,
  .paropts = NULL
)

\method{clean}{raw_mspct}(
  x,
  range = NULL,
  range.s.data = c(0, NA),
  fill = range.s.data,
  ...,
  .parallel = FALSE,
  .paropts = NULL
)

\method{clean}{generic_mspct}(
  x,
  range = x,
  range.s.data = c(NA_real_, NA_real_),
  fill = range.s.data,
  col.names,
  ...,
  .parallel = FALSE,
  .paropts = NULL
)
}
\arguments{
\item{x}{an R object}

\item{range}{numeric vector of wavelengths}

\item{range.s.data}{numeric vector of length two giving the allowable
range for the spectral data.}

\item{fill}{numeric vector of length 1 or 2, giving the replacement
values to use at each extreme of the range.}

\item{...}{currently ignored}

\item{unit.out}{character string with allowed values "energy", and "photon",
or its alias "quantum"}

\item{qty.out}{character string with allowed values "energy", and "photon",
or its alias "quantum"}

\item{min.Afr}{numeric Gives the minimum value accepted for the computed
absorptance. The default \code{NULL} sets a valid value (Afr >= 0) with
a warning. If an integer value is passed to \code{digits} values are
adjusted silently.}

\item{col.names}{character The name of the variable to clean}

\item{.parallel}{if TRUE, apply function in parallel, using parallel backend
provided by foreach}

\item{.paropts}{a list of additional options passed into the foreach function
when parallel computation is enabled. This is important if (for example)
your code relies on external data or packages: use the .export and
.packages arguments to supply them so that all cluster nodes have the
correct environment set up for computing.}
}
\value{
A copy of \code{x}, possibly with some of the spectral data values
  replaced by the value passed to \code{fill}.
}
\description{
These functions implement the equivalent of replace() but for spectral
objects instead of vectors.
}
\section{Methods (by class)}{
\itemize{
\item \code{clean(default)}: Default for generic function

\item \code{clean(source_spct)}: Replace off-range values in a source spectrum

\item \code{clean(filter_spct)}: Replace off-range values in a filter spectrum

\item \code{clean(reflector_spct)}: Replace off-range values in a reflector spectrum

\item \code{clean(solute_spct)}: Replace off-range values in a solute spectrum

\item \code{clean(object_spct)}: Replace off-range values in an object spectrum

\item \code{clean(response_spct)}: Replace off-range values in a response spectrum

\item \code{clean(cps_spct)}: Replace off-range values in a counts per second spectrum

\item \code{clean(raw_spct)}: Replace off-range values in a raw counts spectrum

\item \code{clean(generic_spct)}: Replace off-range values in a generic spectrum

\item \code{clean(source_mspct)}: 

\item \code{clean(filter_mspct)}: 

\item \code{clean(reflector_mspct)}: 

\item \code{clean(object_mspct)}: 

\item \code{clean(solute_mspct)}: 

\item \code{clean(response_mspct)}: 

\item \code{clean(cps_mspct)}: 

\item \code{clean(raw_mspct)}: 

\item \code{clean(generic_mspct)}: 

}}
\note{
In the case of \code{object_spct} objects, cleaning is done first
  on the Rfr and Tfr columns and subsequently Afr estimated and if needed
  half of deviation of Afr from the expected minimum value subtracted from
  each of Rfr and Tfr.
}
