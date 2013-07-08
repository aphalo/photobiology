e2qmol_multipliers <- function(wavelength){
  h = 6.626e-34 # Plank's contsant (Js)
  c = 2.998e8 # speed of light in vaccum (m/s), so we need to convert wavelength to metres
  Na = 6.022e23 # Avogadro's number, photons per mol
  return((h * c / (wavelength * 1e-9)) * Na)
}
