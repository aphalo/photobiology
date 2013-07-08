as_photon_mol <- as_quantum_mol <- function(wavelengths,spectral.e.irradiance){
  return(spectral.e.irradiance / e2qmol_multipliers(wavelengths))
}
