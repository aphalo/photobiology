% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spct.fscale.r
\name{setScaled}
\alias{setScaled}
\alias{setScaled.default}
\alias{setScaled.generic_spct}
\alias{setScaled.summary_generic_spct}
\alias{setScaled.generic_mspct}
\title{Set the "scaled" attribute}
\usage{
setScaled(x, ...)

\method{setScaled}{default}(x, ...)

\method{setScaled}{generic_spct}(x, ..., scaled = FALSE)

\method{setScaled}{summary_generic_spct}(x, ..., scaled = FALSE)

\method{setScaled}{generic_mspct}(x, ..., scaled = FALSE)
}
\arguments{
\item{x}{a generic_spct object.}

\item{...}{currently ignored.}

\item{scaled}{logical with \code{FALSE} meaning that values are expressed in
absolute physical units and \code{TRUE} meaning that relative units are
used. If \code{NULL} the attribute is not modified.}
}
\value{
a new object of the same class as \code{x}.

a new object of the same class as \code{x}.

a new object of the same class as \code{x}.

a new object of the same class as \code{x}.
}
\description{
Function to write the "scaled" attribute of an existing generic_spct
object.
}
\section{Methods (by class)}{
\itemize{
\item \code{setScaled(default)}: Default for generic function

\item \code{setScaled(generic_spct)}: Specialization for generic_spct

\item \code{setScaled(summary_generic_spct)}: Specialization for summary_generic_spct

\item \code{setScaled(generic_mspct)}: Specialization for generic_mspct

}}
\note{
if x is not a \code{generic_spct} object, x is not modified.
}
\seealso{
Other rescaling functions: 
\code{\link{fscale}()},
\code{\link{fshift}()},
\code{\link{getNormalized}()},
\code{\link{getScaled}()},
\code{\link{is_normalized}()},
\code{\link{is_scaled}()},
\code{\link{normalize}()},
\code{\link{setNormalized}()}
}
\concept{rescaling functions}
