% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spct-math.r
\name{MathFun}
\alias{MathFun}
\alias{sqrt.generic_spct}
\alias{abs.generic_spct}
\title{Miscellaneous Mathematical Functions}
\usage{
\method{sqrt}{generic_spct}(x)

\method{abs}{generic_spct}(x)
}
\arguments{
\item{x}{an object of class "generic_spct"}
}
\description{
\code{abs(x)} computes the absolute value of \code{x}, \code{sqrt(x)}
computes the (principal) square root of \code{x}. The functions are applied
to the spectral data, not the wavelengths. The quantity in the spectrum to
which the function is applied depends on the class of \code{x} and the
current value of output options.
}
\seealso{
Other math operators and functions: 
\code{\link{^.generic_spct}()},
\code{\link{convolve_each}()},
\code{\link{div-.generic_spct}},
\code{\link{log}()},
\code{\link{minus-.generic_spct}},
\code{\link{mod-.generic_spct}},
\code{\link{plus-.generic_spct}},
\code{\link{round}()},
\code{\link{sign}()},
\code{\link{slash-.generic_spct}},
\code{\link{times-.generic_spct}}
}
\concept{math operators and functions}
