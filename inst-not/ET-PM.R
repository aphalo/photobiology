#' Potential evapotranspiration
#'
#' Formulation of the original Penman-Monteith equation used in the computation
#' of the forest fire index (FFI) by the Finnish Meteorological Institute.
#' Compared to the simplified formulation used by FAO this formualtion requires
#' downwelling long wave radiation as input. The present version is based on the
#' FORTRAN code by Ari Venäläinen from 1996. However, we compute saturated water
#' vapour pressure and its slope vs. temperature using Tetens' formulation as
#' implemented in functions \code{water_vp_sat()} and
#' \code{water_vp_sat_slope()}.
#'
#' @param temperature numeric vector of air temperatures (C) at 2 m height.
#' @param water.vp numeric vector of water vapour pressure in air (Pa).
#' @param wind.speed numeric Wind speed (m/s) at 2 m height.
#' @param net.irradiance numeric Net radiation balance (W/m2).
#' @param lw.emissivity numeric Emissivity for long wave radiation (/1)
#' @param Rs numeric Surface resistance (s/m).
#' @param atmospheric.pressure numeric Atmospheric pressure (Pa).
#' @param method character The name of an estimation method.
#' @param check.range logical Flag indicating whether to check or not that
#'   arguments for temperature are within range of method. Passed to
#'   function calls to \code{water_vp_sat()} and \code{water_vp_sat_slope()}.
#'
#' @return A numeric vector of reference evapotranspiration estimates expressed
#'   as mm / h.
#'
#' @export
#'
#' @family Evapotranspiration and energy balance related functions.
#'
#' @examples
#' # instantaneous
#' ET_PM(temperature = 20,
#'       water.vp = water_RH2vp(0.7, 20),
#'       wind.speed = 0,
#'       net.irradiance = 10,
#'       Rs = 25)
#'
ET_PM <- function(temperature,
                  water.vp,
                  wind.speed,
                  net.irradiance,
                  lw.emissivity = 0.98,
                  Rs = 70,
                  atmospheric.pressure = 1013e-2,
                  method = "FMI.PM",
                  check.range = TRUE) {
  if (method == "PM") {

  } else if (method == "FMI.PM") {
    wind.speed <- max(wind.speed, 0.5)
    sigma <- 5.670374419e-8 # Stefan–Boltzmann constant (W⋅m−2⋅)−
    rho <- 1.2923 # mean air density
    cp <- 1004 # specific heat
    l <- 2.5e6
    k <- 421.0825 # 0.6 * ((log(8 / 0.0002) / 0.4)^2)
    Ra <-  k / wind.speed
    e.sat <- water_vp_sat(temperature,
                          over.ice = temperature <= -5,
                          check.range = check.range) * 10e-2 # Pa -> mb
    es.slope <- water_vp_sat_slope(temperature,
                                   over.ice = temperature <= -5,
                                   check.range = check.range) * 10e-2  # Pa / K -> mb / K
    height.cor <- 4 * sigma * lw.emissivity * ((273.15 + temperature)^3)
    (es.slope * net.irradiance + rho * cp * (1 + height.cor * Ra / (rho * cp)) *
        (e.sat - water.vp) / Ra) /
      (es.slope + 0.66 * (1 + height.cor * Ra / (rho * cp))) / l * 3600 * 3 # s -> h
  } else {
    warning("Method '", method,
            "' unavailable; use 'FMI.PM'")
    return(NA_real_, length(temperature))
  }
}
