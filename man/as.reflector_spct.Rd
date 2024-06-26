% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spct.new.r
\name{as.reflector_spct}
\alias{as.reflector_spct}
\alias{as.reflector_spct.default}
\title{Coerce to a spectrum}
\usage{
as.reflector_spct(x, ...)

\method{as.reflector_spct}{default}(
  x,
  Rfr.type = c("total", "specular"),
  strict.range = getOption("photobiology.strict.range", default = FALSE),
  ...
)
}
\arguments{
\item{x}{an R object.}

\item{...}{other arguments passed to "set" functions.}

\item{Rfr.type}{a character string, either \code{"total"} or
\code{"specular"}.}

\item{strict.range}{logical Flag indicating whether off-range values result
in an error instead of a warning.}
}
\value{
A copy of \code{x} converted into a \code{reflector_spct} object.
}
\description{
Return a copy of an R object with its class set to a given type of spectrum.
}
\section{Methods (by class)}{
\itemize{
\item \code{as.reflector_spct(default)}: 

}}
\seealso{
\code{\link{setGenericSpct}}

Other constructors of spectral objects: 
\code{\link{as.calibration_spct}()},
\code{\link{as.chroma_spct}()},
\code{\link{as.cps_spct}()},
\code{\link{as.filter_spct}()},
\code{\link{as.generic_spct}()},
\code{\link{as.object_spct}()},
\code{\link{as.raw_spct}()},
\code{\link{as.response_spct}()},
\code{\link{as.solute_spct}()},
\code{\link{as.source_spct}()},
\code{\link{source_spct}()}
}
\concept{constructors of spectral objects}
