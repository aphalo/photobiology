#' Water vapour pressure
#'
#' Approximate water pressure in air as a function of temperature, and its
#' inverse the calculation of dewpoint.
#'
#' @param temperature numeric vector of air temperatures (C).
#' @param water.vp numeric vector of water vapour pressure in air (Pa).
#' @param water.mvc numeric vector of water vapour concnetration as mass per
#'   volume (\eqn{g m^{-3}}).
#' @param relative.humidity numeric Relative humidity as % (default) or as a
#'   fraction of 1.
#' @param over.ice logical vector Is the estimate for equilibrium with liquid
#'   water or with ice.
#' @param method character Currently "tetens", modified "magnus", "wexler" and
#'   "goff.gratch" equations  are supported.
#' @param pc logical flag for result returned as percent or not.
#' @param check.range logical Flag indicating whether to check or not that
#'   arguments for temperature are within the range of validity of the
#'   \code{method} used.
#'
#' @details Function \code{water_vp_sat()} provides implementations of several
#'   well known equations for the estimation of saturation vapor pressure in
#'   air. Functions \code{water_dp()} and \code{water_fp()} use the inverse of
#'   these equations to compute the dew point or frost point from water vapour
#'   pressure in air. The inverse functions are either analytical solutions or
#'   fitted approximations. None of these functions are solved numerically by
#'   iteration.
#'
#'   Method \code{"tetens"} implements Tetens' (1930) equation for the cases of
#'   equilibrium with a water and an ice surface. Method \code{"magnus"}
#'   implements the modified Magnus equations of Alduchov and Eskridge (1996,
#'   eqs. 21 and 23). Method \code{"wexler"} implements the equations proposed
#'   by Wexler (1976, 1977), and their inverse according to Hardy (1998). Method
#'   \code{"goff.gratch"} implements the equations of Groff and Gratch (1946)
#'   with the minor updates of Groff (1956).
#'
#'   The equations are approximations, and in spite of their different names,
#'   Tetens' and Magnus' equations have the same form with the only difference
#'   in the values of the parameters. However, the modified Magnus equation is
#'   more accurate as Tetens equation suffers from some bias errors at extreme
#'   low temperatures (< -40 C). In contrast Magnus equations with recently
#'   fitted values for the parameters are usable for temperatures from -80 C to
#'   +50 C over water and -80 C to 0 C over ice. The Groff Gratch equation is
#'   more complex and is frequently used as a reference in comparison as it is
#'   considered reliable over a broad range of temperatures. Wexler's equations
#'   are computationally simpler and fitted to relatively recent data. There is
#'   little difference at temperatures in the range -20 C to +50 C, and
#'   differences become large at extreme temperatures. Temperatures outside the
#'   range where estimations are highly reliable for each equation return
#'   \code{NA}, unless extrapolation is enabled by passing \code{FALSE} as
#'   argument to parameter \code{check.range}.
#'
#'   The switch between equations for ice or water cannot be based on
#'   air temperature, as it depends on the presence or not of a surface of
#'   liquid water. It must be set by passing an argument to parameter
#'   \code{over.ice} which defaults to \code{FALSE}.
#'
#'   Tetens equation is still very frequently used, and is for example the
#'   one recommended by FAO for computing potential evapotranspiration. For this
#'   reason it is used as default here.
#'
#' @note The inverse of the Groff Gratch equation has yet to be implemented.
#'
#' @return A numeric vector of partial pressures in pascal (Pa) for
#'   \code{water_vp_sat()} and \code{water_mvc2vp()}, a numeric vector of dew point
#'   temperatures (C) for \code{water_dp()} and numeric vector of mass per volume
#'   concentrations (\eqn{g m^{-3}}) for \code{water_vp2mvc()}.  \code{water_vp_sat()} and
#'   \code{psychrometric_constant()} both return numeric vectors of pressure per
#'   degree of temperature (\eqn{Pa C^{-1}})
#'
#' @references Tetens, O., 1930. Uber einige meteorologische Begriffe.
#'   Zeitschrift fur Geophysik, Vol. 6:297.
#'
#'   Goff, J. A., and S. Gratch (1946) Low-pressure properties of water from
#'   -160 to 212 F, in Transactions of the American Society of Heating and
#'   Ventilating Engineers, pp 95-122, presented at the 52nd annual meeting of
#'   the American Society of Heating and Ventilating Engineers, New York, 1946.
#'
#'   Wexler, A. (1976) Vapor Pressure Formulation for Water in Range 0 to 100°C.
#'   A Revision, Journal of Research ofthe National Bureau of Standards: A.
#'   Physics and Chemistry, September-December 1976, Vol. 80A, Nos.5 and 6,
#'   775-785
#'
#'   Wexler, A., (1977) Vapor Pressure Formulation for Ice, Journal of Research of the
#'   National Bureau of Standards - A. Physics and Chemistry, Vol. 81A, No. 1, 5-19
#'
#'   Alduchov, O. A., Eskridge, R. E., 1996. Improved Magnus Form Approximation
#'   of Saturation Vapor Pressure. Journal of Applied Meteorology, 35: 601-609 .
#'
#'   Hardy, Bob (1998) ITS-90 formulations for vapor pressure, frostpoint
#'   temperature, dewpoint temperature, andenhancement factors in the range -100
#'   TO +100 C. The Proceedings of the Third International Symposium on Humidity
#'   & Moisture, Teddington, London, England, April 1998.
#'   \url{https://www.decatur.de/javascript/dew/resources/its90formulas.pdf}
#'
#'   Monteith, J., Unsworth, M. (2008) Principles of Environmental Physics.
#'   Academic Press, Amsterdam.
#'
#'   Allen R G, Pereira L S, Raes D, Smith M. (1998) Crop evapotranspiration:
#'   Guidelines for computing crop water requirements. FAO Irrigation and
#'   drainage paper 56. Rome: FAO.
#'
#'   [Equations describing the physical properties of moist
#'   air](http://www.conservationphysics.org/atmcalc/atmoclc2.pdf)
#'
#'
#' @export
#'
#' @examples
#' water_vp_sat(20) # C -> Pa
#' water_vp_sat(temperature = c(0, 10, 20, 30, 40)) # C -> Pa
#' water_vp_sat(temperature = -10) # over water!!
#' water_vp_sat(temperature = -10, over.ice = TRUE)
#' water_vp_sat(temperature = 20) / 100 # C -> mbar
#'
#' water_vp_sat(temperature = 20, method = "magnus") # C -> Pa
#' water_vp_sat(temperature = 20, method = "tetens") # C -> Pa
#' water_vp_sat(temperature = 20, method = "wexler") # C -> Pa
#' water_vp_sat(temperature = 20, method = "goff.gratch") # C -> Pa
#'
#' water_vp_sat(temperature = -20, over.ice = TRUE, method = "magnus") # C -> Pa
#' water_vp_sat(temperature = -20, over.ice = TRUE, method = "tetens") # C -> Pa
#' water_vp_sat(temperature = -20, over.ice = TRUE, method = "wexler") # C -> Pa
#' water_vp_sat(temperature = -20, over.ice = TRUE, method = "goff.gratch") # C -> Pa
#'
#' water_dp(water.vp = 1000) # Pa -> C
#' water_dp(water.vp = 1000, method = "magnus") # Pa -> C
#' water_dp(water.vp = 1000, method = "wexler") # Pa -> C
#' water_dp(water.vp = 500, over.ice = TRUE) # Pa -> C
#' water_dp(water.vp = 500, method = "wexler", over.ice = TRUE) # Pa -> C
#'
#' water_fp(water.vp = 300) # Pa -> C
#' water_dp(water.vp = 300, over.ice = TRUE) # Pa -> C
#'
#' water_vp2RH(water.vp = 1500, temperature = 20) # Pa, C -> RH %
#' water_vp2RH(water.vp = 1500, temperature = c(20, 30)) # Pa, C -> RH %
#' water_vp2RH(water.vp = c(600, 1500), temperature = 20) # Pa, C -> RH %
#'
#' water_vp2mvc(water.vp = 1000, temperature = 20) # Pa -> g m-3
#'
#' water_mvc2vp(water.mvc = 30, temperature = 40) # g m-3 -> Pa
#'
#' water_dp(water.vp = water_mvc2vp(water.mvc = 10, temperature = 30)) # g m-3 -> C
#'
water_vp_sat <- function(temperature,
                         over.ice = FALSE,
                         method = "tetens",
                         check.range = TRUE) {
  method <- tolower(method)
  if (length(method) > 1L) {
    if (length(unique(method)) > 1L) {
      stop("Only one method can be used per function call.")
    } else {
      method <- method[1L]
    }
  }
  if (length(over.ice) == 1L && length(temperature > 1)) {
    over.ice <- rep_len(over.ice, length(temperature))
  }
  if (any(temperature > 0 & over.ice)) {
    warning("At temperature > 0 C, ice surface will be wet.")
  }
  if (method == "magnus") {
    if (check.range &&
        any(!is.na(temperature) & (temperature < -80 | temperature > 50))) {
      warning("Out of bounds temperature value(s) set to NA, range: -80 C to +50 C.")
      temperature <- ifelse(temperature < -80, NA_real_, temperature)
    }
    z <-
      ifelse(over.ice,
             611.21 * exp(22.587 * temperature / (273.86 + temperature)),
             610.94 * exp(17.625 * temperature / (243.04 + temperature)))
  } else if (method == "tetens") {
    if (check.range && any(!is.na(temperature) & (temperature < -40))) {
      warning("Out of bounds temperature value(s) set to NA, range: -40 C to +50 C.")
      temperature <- ifelse(temperature < -40, NA_real_, temperature)
    }
    z <-
      ifelse(over.ice,
             610.78 * exp(21.875 * temperature / (265.5 + temperature)),
             610.78 * exp(17.269 * temperature / (237.3 + temperature)))
  } else if (method == "wexler") {
    if (check.range && any(!is.na(temperature) & (temperature < -100 | temperature > 110))) {
      warning("Out of bounds temperature value(s) set to NA, range: -100 C to +100 C")
      temperature <-
        ifelse(temperature < -100 | temperature > 110, NA_real_, temperature)
    }
    temperature.K <- temperature + 273.15
    wexler.ice <-  function(temperature.K) {
      # g_ITS68 <- c(-5.8653696e3, 2.224103300e1, 1.3749042e-2, -3.4031775e-5,
      #        2.6967687e-8, 6.918651e-1)
      g <- c(-5.8666426e3, 2.232870244e1, 1.39387003e-2, -3.4262402e-5,
             2.7040955e-8, 6.7063522e-1)
      exp(sum(temperature.K^((0:4) - 1) * g[1:5]) + g[6] * log(temperature.K))
    }
    wexler.water <-  function(temperature.K) {
      # g_ITS68 <- c(-2.9912729e3, -6.0170128e3, 1.887643854e1, -2.8354721e-2,
      #        1.7838301e-5, -8.4150417e-10, 4.4412543e-13, 2.858487)
      g <- c(-2.8365744e3, -6.028076559e3, 1.954263612e1, -2.737830188e-2,
             1.6261698e-5, 7.0229056e-10, -1.8680009e-13, 2.7150305)
      exp(sum(temperature.K^((0:6) - 2) * g[1:7]) + g[8] * log(temperature.K))
    }
    z <-
      ifelse(over.ice,
        sapply(temperature.K, wexler.ice),
        sapply(temperature.K, wexler.water))
  } else if (method == "goff.gratch") {
    if (check.range && any(temperature < -50)) {
      warning("Out of bounds temperature value(s) set to NA, range: -50 C to +100 C")
      temperature <- ifelse(temperature < -50, NA_real_, temperature)
    }
    temperature.K <- temperature + 273.15
    z <-
      ifelse(over.ice,
             10^(-9.09718 * (273.16 / temperature.K - 1) -
                   3.56654 * log10(273.16 / temperature.K) +
                   0.876793 * (1 - temperature.K / 273.16) +
                   log10(6.1173) ),
             10^(-7.90298 * (373.16 / temperature.K - 1) +
                   5.02808 * log10(373.16 / temperature.K) -
                   1.3816e-7 * (10^(11.344 * (1 - temperature.K / 373.16)) - 1) +
                   8.1328e-3 * (10^(-3.49149 * (373.16 / temperature.K - 1)) - 1) +
                   log10(1013.25)))
    z <- z * 1e2 # hPa -> Pa
  } else {
    warning("Method '", method,
            "' unavailable; use 'magnus', 'tetens', 'wexler' or 'goff.gratch'.")
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
                     method = "tetens",
                     check.range = TRUE) {
  method <- tolower(method)
  if (length(method) > 1L) {
    if (length(unique(method)) > 1L) {
      stop("Only one method can be used per function call.")
    } else {
      method <- method[1L]
    }
  }
  if (length(over.ice) == 1L && length(water.vp > 1)) {
    over.ice <- rep_len(over.ice, length(water.vp))
  }
  if (any(water.vp <= 0)) {
    warning("Dew point is not defined for vapour pressure  <= 0 Pa.")
    water.vp <- ifelse(water.vp <= 0, NA_real_, water.vp)
  }
  if (method == "magnus") {
    z <-
      ifelse(over.ice,
             273.86 * log(water.vp / 611.21) / (22.587 - log(water.vp / 611.21)),
             243.04 * log(water.vp / 610.94) / (17.625 - log(water.vp / 610.94)))
    if (check.range && any((!is.na(z)) & z < -80 | z > 50)) {
      warning("Out of bounds temperature value(s) set to NA, range: -80 C to +50 C.")
      z <- ifelse((!is.na(z)) &  z < -80 | z > 50, NA_real_, z)
    }
  } else if (method == "tetens") {
    z <-
      ifelse(over.ice,
             265.5  * log(water.vp / 610.78) / (21.875 - log(water.vp / 610.78)),
             237.3  * log(water.vp / 610.78) / (17.269 - log(water.vp / 610.78)))
    if (check.range && any((!is.na(z)) & z < -30)) {
      warning("Out of bounds temperature value(s) set to NA, range: -30 C to +50 C.")
      z <- ifelse((!is.na(z)) & z < -30 | z > 50, NA_real_, z)
    }
  } else if (method == "wexler") {
    wexler.inv.ice <-  function(water.vp) {
      c <- c(2.1257969e2, -1.0264612e1, 1.4354796e-1)
      d <- c(1, -8.2871619e-2, 2.3540411e-3, -2.4363951e-5)
      sum(c * log(water.vp)^(0:2)) / sum(d * log(water.vp)^(0:3))
    }
    wexler.inv.water <-  function(water.vp) {
      c <- c(2.0798233e2, -2.0156028e1, 4.6778925e-1, -9.2288067e-6)
      d <- c(1, -1.3319669e-1, 5.6577518e-3, -7.5172865e-5)
      sum(c * log(water.vp)^(0:3)) / sum(d * log(water.vp)^(0:3))
    }
    z <-
      ifelse(over.ice,
             sapply(water.vp, wexler.inv.ice) - 273.15,
             sapply(water.vp, wexler.inv.water) - 273.15)
    if (check.range && any((!is.na(z)) & z < -100 | z > 100)) {
      warning("Out of bounds temperature value(s) set to NA, -100 C to + 100 C.")
      z <- ifelse((!is.na(z)) & z < -100 | z > 100, NA_real_, z)
    }
  } else {
    warning("Method '", method,
            "' unavailable; use 'magnus', 'tetens' or 'wexler'.")
    z <- rep(NA_real_, length(water.vp))
  }
  method <- tolower(method)
  if (any(z > 0 & over.ice)) {
    warning("At dew point temperature > 0 C, ice surface will be wet.")
  }
  z
}

#' @rdname water_vp_sat
#'
#' @export
#'
water_fp <- function(water.vp,
                     over.ice = TRUE,
                     method = "tetens",
                     check.range = TRUE) {
  water_dp(water.vp = water.vp,
           over.ice = over.ice,
           method = method,
           check.range = check.range)
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

#' @rdname water_vp_sat
#'
#' @export
#'
water_vp2RH <- function(water.vp,
                     temperature,
                     over.ice = FALSE,
                     method = "tetens",
                     pc = TRUE,
                     check.range = TRUE) {
  z <- water.vp / water_vp_sat(temperature,
                               over.ice = over.ice,
                               method = method,
                               check.range = check.range)
  if (pc) {
    z * 100
  } else {
    z
  }
}

#' @rdname water_vp_sat
#'
#' @export
#'
water_RH2vp <- function(relative.humidity,
                        temperature,
                        over.ice = FALSE,
                        method = "tetens",
                        pc = TRUE,
                        check.range = TRUE) {
  if (pc) {
    relative.humidity <- relative.humidity * 1e-2
  }
  relative.humidity * water_vp_sat(temperature,
                                   over.ice = over.ice,
                                   method = method,
                                   check.range = check.range)
}

#' @rdname water_vp_sat
#'
#' @param temperature.step numeric Delta or step used to estimate the slope
#'   as a finite difference (C).
#'
#' @export
#'
#' @examples
#'
#' water_vp_sat_slope(temperature = 20) # C -> Pa / C
#'
water_vp_sat_slope <-  function(temperature,
                                over.ice = FALSE,
                                method = "tetens",
                                check.range = TRUE,
                                temperature.step = 0.1) {
  vp_sat1 <- water_vp_sat(temperature + temperature.step / 2,
                          over.ice = over.ice,
                          method = method,
                          check.range = check.range)
  vp_sat2 <- water_vp_sat(temperature - temperature.step / 2,
                          over.ice = over.ice,
                          method = method,
                          check.range = check.range)
  (vp_sat1 - vp_sat2) / temperature.step
}

#' @rdname water_vp_sat
#'
#' @param atmospheric.pressure numeric Atmospheric pressure (Pa).
#'
#' @export
#'
#' @examples
#'
#' psychrometric_constant(atmospheric.pressure = 81.8e3) # Pa -> Pa / C
#'
psychrometric_constant <- function(atmospheric.pressure = 101325) {
  # latent heat of vaporization, 2.45 [MJ kg-1 ],
  lambda <- 2.45
  # specific heat at constant pressure, 1.013 10 -3 [MJ kg-1 °C-1 ]
  C.p <- 1.013e-3
  # ratio molecular weight of water vapour/dry air = 0.622
  epsilon <- 0.622

  # Pa -> Pa / C
  (C.p * atmospheric.pressure) / (epsilon * lambda)
}
