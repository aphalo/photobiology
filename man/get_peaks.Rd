% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spct.deprecated.r
\name{get_peaks}
\alias{get_peaks}
\alias{get_valleys}
\title{Get peaks and valleys from a spectrum}
\usage{
get_peaks(
  x,
  y,
  global.threshold = 0,
  span = 5,
  strict = TRUE,
  x_unit = "",
  x_digits = 3,
  na.rm = FALSE
)

get_valleys(
  x,
  y,
  global.threshold = 0,
  span = 5,
  strict = TRUE,
  x_unit = "",
  x_digits = 3,
  na.rm = FALSE
)
}
\arguments{
\item{x, y}{numeric}

\item{global.threshold}{numeric A value belonging to class \code{"AsIs"} is
interpreted as an absolute minimum height or depth expressed in data units.
A bare \code{numeric} value (normally between 0.0 and 1.0), is interpreted
as relative to \code{threshold.range}. In both cases it sets a
\emph{global} height (depth) threshold below which peaks (valleys) are
ignored. A bare negative \code{numeric} value indicates the \emph{global}
height (depth) threshold below which peaks (valleys) are be ignored. If
\code{global.threshold = NULL}, no threshold is applied and all peaks
returned.}

\item{span}{odd positive integer A peak is defined as an element in a
sequence which is greater than all other elements within a moving window of
width \code{span} centred at that element. The default value is 5, meaning
that a peak is taller than its four nearest neighbours. \code{span = NULL}
extends the span to the whole length of \code{x}.}

\item{strict}{logical flag: if \code{TRUE}, an element must be strictly
greater than all other values in its window to be considered a peak.}

\item{x_unit}{character Vector of texts to be pasted at end of labels built
from x value at peaks.}

\item{x_digits}{numeric Number of significant digits in wavelength label.}

\item{na.rm}{logical indicating whether \code{NA} values should be stripped
before searching for peaks.}
}
\value{
A data frame with variables w.length and s.irrad with their values at
  the peaks or valleys plus a character variable of labels.
}
\description{
These functions "get" (or extract) peaks (maxima) and valleys (minima) in two
vectors, usually a spectral quantity and wavelength, using a user selectable
span for window width and global and local (within moving window) size
thresholds. They also generate \code{character} values for \code{x}.
}
\details{
As \code{\link[photobiology]{find_valleys}},
  \code{\link[photobiology]{peaks}} and \code{\link[photobiology]{valleys}}
  call \code{\link[photobiology]{find_peaks}} to search for peaks and
  valleys, this explanation applies to the four functions. It also applies to
  \code{\link[ggspectra]{stat_peaks}} and
  \code{\link[ggspectra:stat_peaks]{stat_valleys}}. Function
  \code{find_peaks} is a wrapper built onto function
  \code{\link[splus2R]{peaks}} from \pkg{splus2R}, adds support for peak
  height thresholds and handles \code{span = NULL} and non-finite (including
  NA) values differently than \code{splus2R::peaks}. Instead of giving an
  error when \code{na.rm = FALSE} and \code{x} contains \code{NA} values,
  \code{NA} values are replaced with the smallest finite value in \code{x}.
  \code{span = NULL} is treated as a special case and selects \code{max(x)}.
  Passing \code{strict = TRUE} ensures that non-unique global and within window
  maxima are ignored, and can result in no peaks being returned.

  Two tests make it possible to ignore irrelevant peaks. One test
  (\code{global.threshold}) is based on the absolute height of the peaks and
  can be used in all cases to ignore globally low peaks. A second test
  (\code{local.threshold}) is available when the window defined by `span`
  does not include all observations and can be used to ignore peaks that are
  not locally prominent. In this second approach the height of each peak is
  compared to a summary computed from other values within the window of width
  equal to \code{span} where it was found. In this second case, the reference
  value used within each window containing a peak is given by the argument
  passed to \code{local.reference}. Parameter \code{threshold.range}
  determines how the values passed as argument to \code{global.threshold} and
  \code{local.threshold} are scaled. The default, \code{NULL} uses the range
  of \code{x}. Thresholds for ignoring too small peaks are applied after
  peaks are searched for, and threshold values can in some cases result in no
  peaks being returned.

  The \code{local.threshold} argument is used \emph{as is} when
  \code{local.reference} is \code{"median"} or \code{"farthest"}, i.e., the
  same distance between peak and reference is used as cut-off irrespective of
  the value of the reference. In cases when the prominence of peaks is
  positively correlated with the baseline, a \code{local.threshold} that
  increases together with increasing computed within window median or
  farthest value applies apply a less stringent height requirement in regions
  with overall low height. In this case, natural logarithm or square root
  weighting can be requested with \code{local.reference} arguments
  \code{"median.log"}, \code{"farthest.log"}, \code{"median.sqrt"}, and
  \code{"farthest.sqrt"} as arguments for \code{local.reference}.

  While functions \code{\link{find_peaks}} and \code{\link{find_valleys}}
  accept as input a \code{numeric} vector and return a \code{logical} vector,
  methods \code{\link{peaks}} and \code{\link{valleys}} accept as input
  different R objects, including spectra and collections of spectra and
  return a subset of the object. These methods are implemented using calls to
  functions \code{find_peaks}, \code{find_valleys} and
  \code{\link{fit_peaks}}.
}
\note{
The use of these two functions is deprecated. They are retained for
  backwards compatibility and will be removed in the near future.
}
\seealso{
Other peaks and valleys functions: 
\code{\link{find_peaks}()},
\code{\link{find_spikes}()},
\code{\link{peaks}()},
\code{\link{replace_bad_pixs}()},
\code{\link{spikes}()},
\code{\link{valleys}()},
\code{\link{wls_at_target}()}
}
\concept{peaks and valleys functions}
