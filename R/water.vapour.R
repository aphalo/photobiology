#' Water vapour pressure
#'
#' Approximate water pressure in air as a function of temperature, and its inverse
#' the calculation of dewpoint.
#'
#' @param temperature numeric vector of air temperatures (C).
#' @param water.vp numeric vector of water vapour pressure in air (Pa).
#' @param water.mvc numeric vector of water vapour concnetration as mass per volume (g m-3).
#' @param over.ice logical Is the estimate for equilibrium with liquid water or with ice.
#' @param method character Currently "tetens" and modified "magnus" equations
#'   are supported.
#'
#' @details This is an implementation of Tetens (1930) equation for the cases of
#'   equilibrium with a water and an ice surface and of the modified Magnus
#'   equations of Alduchov and Eskridge (1996, eqs. 21 and 23).
#'
#'   The equations are approximations, and in spite of the different names, have
#'   the same form with the only difference in the values of the parameters.
#'   However, the modified Magnus equation is more accurate as Tetens equation
#'   suffers from some bias errors at extreme low temperatures (< -40 C). In
#'   contrast Magnus equations with recently fitted values for the parameters
#'   are usable for temperatures from -80 C to +50 C over water and -80 C to 0 C
#'   over ice. The switch between equations for ice or water cannot be based on
#'   air temperature, as it depends on the presence or not of a surface of
#'   liquid water.
#'
#' @note Tetens equation is still very frequently used, and is for example the
#'   one recommended by FAO for computing potential evapotranspiration. For
#'   this reason it is used as default here.
#'
#' @return A numeric vector of partial pressures in pascal (P) for
#'   \code{water_vp_sat} and \code{water_mvc2vp}, a numeric vector of dew point
#'   temperatures (C) for \code{water_dp} and numeric vector of mass per
#'   volume concentrations (g m-3) for \code{water_vp2mvc}.
#'
#' @references
#'  Tetens, O., 1930. Uber einige meteorologische Begriffe. Zeitschrift fur
#'  Geophysik, Vol. 6:297.
#'
#'  Alduchov, O. A., Eskridge, R. E., 1996. Improved Magnus Form Approximation of
#'  Saturation Vapor Pressure. Journal of Applied Meteorology, 35: 601-609 .
#'
#'  Monteith, J., Unsworth, M. (2008) Principles of Environmental Physics.
#'  Academic Press, Amsterdam.
#'
#'  [Equations describing the physical properties of moist air](http://www.conservationphysics.org/atmcalc/atmoclc2.pdf)
#'
#' @export
#'
#' @examples
#' water_vp_sat(20) # C -> Pa
#' water_vp_sat(-10) # over water!!
#' water_vp_sat(-10, over.ice = TRUE)
#' water_vp_sat(20) / 100 # C -> mbar
#'
#' water_dp(1000) # Pa -> C
#'
#' water_vp2mvc(1000, 20) # Pa -> g m-3
#'
#' water_mvc2vp(30, 40) # g m-3 -> Pa
#'
#' water_dp(water_mvc2vp(10, 30)) # g m-3 -> C
#'
water_vp_sat <- function(temperature,
                         over.ice = FALSE,
                         method = "tetens") {
  method <- tolower(method)
  if (any(temperature > 0) && over.ice) {
    warning("At temperature > 0 C, ice surface will be wet.")
  }
  if (method == "magnus") {
    if (over.ice) {
      z <- 611.21 * exp(22.587 * temperature / (273.86 + temperature))
    } else {
      z <- 610.94 * exp(17.625 * temperature / (243.04 + temperature))
    }
  } else if (method == "tetens") {
    if (over.ice) {
      z <- 610.78 * exp(21.875 * temperature / (265.5 + temperature))
    } else {
      z <- 610.78 * exp(17.269 * temperature / (237.3 + temperature))
    }
  } else {
    warning("Method '", method, "' is unknown.")
    z <- rep(NA_real_, length(temperature))
  }
  z
}

#' @rdname water_vp_sat
#'
#' @export
#'
water_dp <- function(water.vp,
                     over.ice = FALSE,
                     method = "tetens") {
  method <- tolower(method)
  if (method == "magnus") {
    if (over.ice) {
      z <- 273.86 * log(water.vp / 611.21) / (22.587 - log(water.vp / 611.21))
    } else {
      z <- 243.04 * log(water.vp / 610.94) / (17.625 - log(water.vp / 610.94))
    }
  } else if (method == "tetens") {
    if (over.ice) {
      z <- 265.5  * log(water.vp / 610.78) / (21.875 - log(water.vp / 610.78))
    } else {
      z <- 237.3  * log(water.vp / 610.78) / (17.269 - log(water.vp / 610.78))
    }
  } else {
    warning("Method '", method, "' is unknown.")
    z <- rep(NA_real_, length(water.vp))
  }
  method <- tolower(method)
  if (any(z > 0) && over.ice) {
    warning("At dew point temperature > 0 C, ice surface will be wet.")
  }
  z
}

#' @rdname water_vp_sat
#'
#' @export
#'
water_vp2mvc <- function(water.vp,
                         temperature) {
  2.166 * water.vp / (temperature + 273.16)
  # 273.16 K is the temperature of the triple point
  # while 0 C, the temperature of ice, is 273.15 K
}

#' @rdname water_vp_sat
#'
#' @export
#'
water_mvc2vp <- function(water.mvc,
                         temperature) {
  (temperature + 273.16) / 2.166 * water.mvc
}
