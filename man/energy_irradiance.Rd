% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/energy.irradiance.r
\name{energy_irradiance}
\alias{energy_irradiance}
\title{Calculate (energy) irradiance from spectral irradiance}
\usage{
energy_irradiance(
  w.length,
  s.irrad,
  w.band = NULL,
  unit.in = "energy",
  check.spectrum = TRUE,
  use.cached.mult = FALSE,
  use.hinges = getOption("photobiology.use.hinges", default = NULL)
)
}
\arguments{
\item{w.length}{numeric vector of wavelength (nm).}

\item{s.irrad}{numeric vector of spectral irradiances, by default as energy (W
m-2 nm-1).}

\item{w.band}{waveband.}

\item{unit.in}{a character Allowed values "photon" or "energy", default is
"energy".}

\item{check.spectrum}{logical Flag indicating whether to sanity check input
data, default is TRUE.}

\item{use.cached.mult}{logical Flag indicating whether multiplier values
should be cached between calls.}

\item{use.hinges}{logical Flag indicating whether to insert "hinges" into the
spectral data before integration so as to reduce interpolation errors at
the boundaries of the wavebands.}
}
\value{
A single numeric value with no change in scale factor: [W m-2 nm-1]
  -> [W m-2].
}
\description{
Energy irradiance for a waveband from a radiation spectrum, optionally
applying a "biological spectral weighting function" or BSWF.
}
\examples{
with(sun.data, energy_irradiance(w.length, s.e.irrad))
with(sun.data, energy_irradiance(w.length, s.e.irrad, new_waveband(400,700)))

}
\seealso{
Other low-level functions operating on numeric vectors.: 
\code{\link{as_energy}()},
\code{\link{as_quantum_mol}()},
\code{\link{calc_multipliers}()},
\code{\link{div_spectra}()},
\code{\link{energy_ratio}()},
\code{\link{insert_hinges}()},
\code{\link{integrate_xy}()},
\code{\link{interpolate_spectrum}()},
\code{\link{irradiance}()},
\code{\link{l_insert_hinges}()},
\code{\link{oper_spectra}()},
\code{\link{photon_irradiance}()},
\code{\link{photon_ratio}()},
\code{\link{photons_energy_ratio}()},
\code{\link{prod_spectra}()},
\code{\link{s_e_irrad2rgb}()},
\code{\link{split_energy_irradiance}()},
\code{\link{split_photon_irradiance}()},
\code{\link{subt_spectra}()},
\code{\link{sum_spectra}()},
\code{\link{trim_tails}()},
\code{\link{v_insert_hinges}()},
\code{\link{v_replace_hinges}()}
}
\concept{low-level functions operating on numeric vectors.}
