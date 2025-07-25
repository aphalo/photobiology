% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source-spct-fractions.R
\name{e_fraction}
\alias{e_fraction}
\alias{e_fraction.default}
\alias{e_fraction.source_spct}
\alias{e_fraction.source_mspct}
\title{Energy:energy fraction}
\usage{
e_fraction(
  spct,
  w.band.num,
  w.band.denom,
  scale.factor,
  wb.trim,
  use.cached.mult,
  use.hinges,
  ...
)

\method{e_fraction}{default}(
  spct,
  w.band.num,
  w.band.denom,
  scale.factor,
  wb.trim,
  use.cached.mult,
  use.hinges,
  ...
)

\method{e_fraction}{source_spct}(
  spct,
  w.band.num = NULL,
  w.band.denom = NULL,
  scale.factor = 1,
  wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
  use.cached.mult = FALSE,
  use.hinges = NULL,
  quantity = "total",
  naming = "short",
  name.tag = NULL,
  ...
)

\method{e_fraction}{source_mspct}(
  spct,
  w.band.num = NULL,
  w.band.denom = NULL,
  scale.factor = 1,
  wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
  use.cached.mult = FALSE,
  use.hinges = NULL,
  quantity = "total",
  naming = "short",
  name.tag = ifelse(naming != "none", "[e:e]", ""),
  ...,
  attr2tb = NULL,
  idx = "spct.idx",
  .parallel = FALSE,
  .paropts = NULL
)
}
\arguments{
\item{spct}{source_spct}

\item{w.band.num}{waveband object or a list of waveband objects used to
compute the numerator(s) and denominator(s) of the fraction(s).}

\item{w.band.denom}{waveband object or a list of waveband objects used to
compute the denominator(s) of the fraction(s).}

\item{scale.factor}{numeric vector of length 1, or length equal to that of
\code{w.band}. Numeric multiplier applied to returned values.}

\item{wb.trim}{logical if TRUE wavebands crossing spectral data boundaries
are trimmed, if FALSE, they are discarded}

\item{use.cached.mult}{logical Flag telling whether multiplier values should
be cached between calls.}

\item{use.hinges}{logical Flag indicating whether to insert "hinges" into the
spectral data before integration so as to reduce interpolation errors at
the boundaries of the wavebands.}

\item{...}{other arguments (possibly used by derived methods).}

\item{quantity}{character One of "total", "average" or "mean".}

\item{naming}{character one of "long", "default", "short" or "none". Used to
select the type of names to assign to returned value.}

\item{name.tag}{character Used to tag the name of the returned values.}

\item{attr2tb}{character vector, see \code{\link{add_attr2tb}} for the syntax
for \code{attr2tb} passed as is to formal parameter \code{col.names}.}

\item{idx}{character Name of the column with the names of the members of the
collection of spectra.}

\item{.parallel}{if TRUE, apply function in parallel, using parallel backend
provided by foreach.}

\item{.paropts}{a list of additional options passed into the foreach function
when parallel computation is enabled. This is important if (for example)
your code relies on external data or packages: use the .export and
.packages arguments to supply them so that all cluster nodes have the
correct environment set up for computing.}
}
\value{
In the case of methods for individual spectra, a \code{numeric}
  vector with name attribute set. The name is based on the name of the
  wavebands unless a named list of wavebands is supplied in which case the
  names of the list elements are used. "[e:e]" is appended if \code{quantity
  = "total"} and "[e(wl):e(wl)]" if \code{quantity = "mean"} or
  \code{quantity = "average"}.

  A \code{data.frame} is returned in the case of collections of spectra,
  containing one column for each fraction definition, an index column with
  the names of the spectra, and optionally additional columns with metadata
  values retrieved from the attributes of the member spectra.

  Fraction definitions are "assembled" from the arguments passed to
  \code{w.band.num} and \code{w.band.denom}. If both arguments are lists of
  waveband definitions, with an equal number of members, then the wavebands
  are paired to obtain as many fractions as the number of wavebands in each
  list. Recycling for wavebands takes place when the number of denominator
  and numerator wavebands differ.
}
\description{
This function returns the energy fraction for a given pair of wavebands of a
light source spectrum.
}
\details{
With the default \code{quantity = "total"} the fraction is based on
  two \strong{energy irradiances}, one computed for each waveband.

  \deqn{\frac{E(s, wb_\mathrm{num})}{E(s, wb_\mathrm{denom}) + E(s, wb_\mathrm{num})}}

If the argument is set to \code{quantity = "mean"} or
 \code{quantity = "average"} the ratio is based on two
 \strong{mean spectral energy irradiances}, one computed for each waveband.

  \deqn{\frac{\overline{Q_\lambda}(s, wb_\mathrm{num})}{\overline{Q_\lambda}(s, wb_\mathrm{denom}) + \overline{Q_\lambda}(s, wb_\mathrm{num})}}

Only if the wavelength expanse of the two wavebands is the same, these two
ratios are numerically identical.
}
\section{Methods (by class)}{
\itemize{
\item \code{e_fraction(default)}: Default for generic function

\item \code{e_fraction(source_spct)}: Method for \code{source_spct} objects

\item \code{e_fraction(source_mspct)}: Calculates energy:energy fraction from a
\code{source_mspct} object.

}}
\note{
Recycling for wavebands takes place when the number of denominator and
  denominator wavebands differ. The last two parameters control speed
  optimizations. The defaults should be suitable in most cases. If you will
  use repeatedly the same SWFs on many spectra measured at exactly the same
  wavelengths you may obtain some speed up by setting
  \code{use.cached.mult=TRUE}. However, be aware that you are responsible for
  ensuring that the wavelengths are the same in each call, as the only test
  done is for the length of the \code{w.length} vector.
}
\examples{
e_fraction(sun.spct, new_waveband(400,700), new_waveband(400,500))

}
\seealso{
Other photon and energy ratio functions: 
\code{\link{e_ratio}()},
\code{\link{eq_ratio}()},
\code{\link{q_fraction}()},
\code{\link{q_ratio}()},
\code{\link{qe_ratio}()}
}
\concept{photon and energy ratio functions}
