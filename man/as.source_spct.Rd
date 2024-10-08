% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spct.new.r
\name{as.source_spct}
\alias{as.source_spct}
\alias{as.source_spct.default}
\title{Coerce to a spectrum}
\usage{
as.source_spct(x, ...)

\method{as.source_spct}{default}(
  x,
  time.unit = c("second", "day", "exposure"),
  bswf.used = c("none", "unknown"),
  strict.range = getOption("photobiology.strict.range", default = FALSE),
  ...
)
}
\arguments{
\item{x}{an R object.}

\item{...}{other arguments passed to "set" functions.}

\item{time.unit}{character string indicating the time unit used for spectral
irradiance or exposure (\code{"second"}, \code{"day"} or \code{"exposure"})
or an object of class duration as defined in package lubridate.}

\item{bswf.used}{character A string indicating the BSWF used, if any, for
spectral effective irradiance or exposure (\code{"none"} or the name of the
BSWF).}

\item{strict.range}{logical Flag indicating whether off-range values result
in an error instead of a warning.}
}
\value{
A copy of \code{x} converted into a \code{source_spct} object.
}
\description{
Return a copy of an R object with its class set to a given type of spectrum.
}
\section{Methods (by class)}{
\itemize{
\item \code{as.source_spct(default)}: 

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
\code{\link{as.reflector_spct}()},
\code{\link{as.response_spct}()},
\code{\link{as.solute_spct}()},
\code{\link{source_spct}()}
}
\concept{constructors of spectral objects}
