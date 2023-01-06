#' Wavelength conversions
#'
#' Convert wavelength (nm) into wave number, frequency (Hz) or energy per photon
#' (J, or eV) and back.
#'
#' @details These functions always expect as input and return wavelengths
#'   expressed in nanometres (nm) as all other functions in the R for
#'   photobiology suite of packages. Conversions depend on Plank's constant,
#'   \emph{h}, the speed of light in vacuum, \emph{c}, and Avogadro's number,
#'   \eqn{N_A}. The values used for these constants have at least nine
#'   significant digits.
#'
#' @param w.length numeric wavelength (nm)
#' @param unit.exponent integer Exponent of the scale multiplier implicit in
#'   result, e.g., use 3 for kJ.
#' @param unit character One of "joule" or "eV".
#' @param frequency numeric Frequency in Hz, possibly
#'   with the scale factor according to \code{unit.exponent}.
#' @param photon.energy numeric Energy of one photon in joule or eV, possibly
#'   with a scale factor according to \code{unit.exponent}.
#' @param wavenumber numeric Wave number in waves per metre, possibly
#'   with a scale factor according to \code{unit.exponent}.
#'
#' @export
#'
#' @examples
#'
#' wl2wavenumber(600) # wavelength in nm -> wave number
#' wavenumber2wl(1666666.66) # wave number -> wavelength in nm
#' wl2frequency(600) # wavelength in nm -> wave frequency (Hz)
#' frequency2wl(499654096666667) # wave frequency (Hz) -> wavelength in nm
#' wl2energy(600) # wavelength in nm -> energy of one photon (J)
#' wl2energy(600, unit = "eV") # wavelength in nm -> energy of one photon (eV)
#' wl2energy(600,
#'           unit.exponent = -3,
#'           unit = "eV")  # wavelength in nm -> energy of one photon (meV)
#' energy2wl(2066.40330,
#'           unit.exponent = -3,
#'           unit = "eV")  # energy of one photon (meV) -> wavelength (nm)
#'
# In each function we first compute the combined constant value that does
# not depend on wavelength, and then use it do the conversion, which usually
# involves a long vector.
wl2wavenumber <- function(w.length,
                          unit.exponent = 0) {
  k <- 1e9 * 10^unit.exponent
  k / w.length # wave number [1/m, by default]
}

#' @rdname wl2wavenumber
#'
#' @export
#'
wavenumber2wl <- function(wavenumber,
                          unit.exponent = 0) {
  k <- 1e9 * 10^unit.exponent
  k / wavenumber # wavelength [nm]
}

#' @rdname wl2wavenumber
#'
#' @export
#'
wl2frequency <- function(w.length,
                         unit.exponent = 0) {
  k <- 299792458 * 1e9 / 10^unit.exponent # speed of light [m/s] / w.length [m]
  k / w.length # frequency, [Hz, by default]
}

#' @rdname wl2wavenumber
#'
#' @export
#'
frequency2wl <- function(frequency,
                         unit.exponent = 0) {
  k <- 299792458 * 1e9 / 10^unit.exponent # speed of light [m/s] / frequency [Hz]
  k / frequency # wavelength [nm]
}

#' @rdname wl2wavenumber
#'
#' @export
#'
wl2energy <- function(w.length,
                      unit.exponent = 0,
                      unit = "joule") {
  if (unit %in% c("SI", "joule")) {
    # Plank's constant [J/Hz] * speed of light [m/s] / w.length [m] -> joule
    k <- 6.62607015e-34 * 299792458 * 1e9 / 10^unit.exponent
  } else if (unit == "eV") {
    # Plank's constant [eV/Hz] * speed of light [m/s] / w.length [m] -> joule
    k <- 4.135667696e-15 * 299792458 * 1e9 / 10^unit.exponent
  }
  k / w.length
}

#' @rdname wl2wavenumber
#'
#' @export
#'
energy2wl <- function(photon.energy,
                      unit.exponent = 0,
                      unit = "joule") {
  if (unit %in% c("SI", "joule")) {
    # Plank's constant [J/Hz] * speed of light [m/s] / energy [J]
    k <- 6.62607015e-34 * 299792458 * 1e9 / 10^unit.exponent
  } else if (unit == "eV") {
    # Plank's constant [eV/Hz] * speed of light [m/s] / energy [eV]
    k <- 4.135667696e-15 * 299792458 * 1e9 / 10^unit.exponent
  }
  k / photon.energy # wavelength [nm]
}

