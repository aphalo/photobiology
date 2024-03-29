% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.energy.r
\name{as_energy}
\alias{as_energy}
\title{Convert spectral photon irradiance into spectral energy irradiance}
\usage{
as_energy(w.length, s.qmol.irrad)
}
\arguments{
\item{w.length}{numeric vector of wavelengths [\eqn{nm}]).}

\item{s.qmol.irrad}{numeric vector of spectral photon irradiance values.}
}
\value{
A numeric vector of spectral (energy) irradiances.
}
\description{
Convert a spectral photon irradiance [\eqn{mol\,s^{-1}\,m^{-2}\,nm^{-1}}{mol
s-1 m-2 nm-1}] into a spectral energy irradiance [\eqn{W\,m^{-2}\,nm^{-1}}{W
m-2 nm-1}].
}
\examples{
with(sun.spct, as_energy(w.length, s.q.irrad))

}
\seealso{
Other low-level functions operating on numeric vectors.: 
\code{\link{as_quantum_mol}()},
\code{\link{calc_multipliers}()},
\code{\link{div_spectra}()},
\code{\link{energy_irradiance}()},
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
