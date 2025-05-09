% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spct-conversion-qty-unit.R
\name{Afr2T}
\alias{Afr2T}
\alias{Afr2T.default}
\alias{Afr2T.numeric}
\alias{Afr2T.filter_spct}
\alias{Afr2T.object_spct}
\alias{Afr2T.filter_mspct}
\alias{Afr2T.object_mspct}
\title{Convert transmittance into absorptance.}
\usage{
Afr2T(x, action, byref, clean, ...)

\method{Afr2T}{default}(x, action = NULL, byref = FALSE, clean = FALSE, ...)

\method{Afr2T}{numeric}(x, action = NULL, byref = FALSE, clean = FALSE, Rfr = NA_real_, ...)

\method{Afr2T}{filter_spct}(x, action = "add", byref = FALSE, clean = FALSE, ...)

\method{Afr2T}{object_spct}(x, action = "add", byref = FALSE, clean = FALSE, ...)

\method{Afr2T}{filter_mspct}(
  x,
  action = "add",
  byref = FALSE,
  clean = FALSE,
  ...,
  .parallel = FALSE,
  .paropts = NULL
)

\method{Afr2T}{object_mspct}(
  x,
  action = "add",
  byref = FALSE,
  clean = FALSE,
  ...,
  .parallel = FALSE,
  .paropts = NULL
)
}
\arguments{
\item{x}{an R object}

\item{action}{character Allowed values "replace" and "add"}

\item{byref}{logical indicating if new object will be created by reference or by copy of x}

\item{clean}{logical replace off-boundary values before conversion}

\item{...}{not used in current version}

\item{Rfr}{numeric vector. Spectral reflectance o reflectance factor.
Set to zero if \code{x} is internal reflectance,}

\item{.parallel}{if TRUE, apply function in parallel, using parallel backend
provided by foreach}

\item{.paropts}{a list of additional options passed into the foreach function
when parallel computation is enabled. This is important if (for example)
your code relies on external data or packages: use the .export and
.packages arguments to supply them so that all cluster nodes have the
correct environment set up for computing.}
}
\value{
A copy of \code{x} with a column \code{Tfr} added and other columns
  possibly deleted except for \code{w.length}. If \code{action = "replace"},
  in all cases, the additional columns are removed, even if no column needs
  to be added.
}
\description{
Function that converts transmittance (fraction) into absorptance (fraction).
If reflectance (fraction) is available, it allows conversions between
internal and total absorptance.
}
\section{Methods (by class)}{
\itemize{
\item \code{Afr2T(default)}: Default method for generic function

\item \code{Afr2T(numeric)}: Default method for generic function

\item \code{Afr2T(filter_spct)}: Method for filter spectra

\item \code{Afr2T(object_spct)}: Method for object spectra

\item \code{Afr2T(filter_mspct)}: Method for collections of filter spectra

\item \code{Afr2T(object_mspct)}: Method for collections of object spectra

}}
\examples{
T2Afr(Ler_leaf.spct)

}
\seealso{
Other quantity conversion functions: 
\code{\link{A2T}()},
\code{\link{T2A}()},
\code{\link{T2Afr}()},
\code{\link{any2T}()},
\code{\link{as_quantum}()},
\code{\link{e2q}()},
\code{\link{e2qmol_multipliers}()},
\code{\link{e2quantum_multipliers}()},
\code{\link{q2e}()}
}
\concept{quantity conversion functions}
