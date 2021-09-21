#' Evapotranspiration
#'
#' Compute an estimate of potential evapotranspiration from
#' meteorologial data. Evapotranspiration from vegetation includes
#' transpiraction by plants plus evaporation from the soil or other wet
#' surfaces. $ET_{ref}$ is the reference value assuming no limitation to
#' transpiration due to soil water, similar to potential evapotranspiration
#' (PET). An actual value $ET_{ref}$ can be estimated only if additional
#' information on the plants and soil is available.
#'
#' @details Currently three methods, based on the Penman-Monteith equation
#'   formulated as recommended by FAO56 (Allen et al., 1998) as well as modified
#'   in 2005 for tall and short vegetation according to ASCE-EWRI. The
#'   computations rely on data measured according WHO standards at 2 m above
#'   ground level to estimate reference evapotranspiration ($ET_{ref}$). The
#'   formulations are those for ET expressed in mm/h, but modified to use as
#'   input flux rates in W/m2 and pressures expressed in Pa.
#'
#' @param temperature numeric vector of air temperatures (C) at 2 m height.
#' @param water.vp numeric vector of water vapour pressure in air (Pa).
#' @param atmospheric.pressure numeric Atmospheric pressure (Pa).
#' @param wind.speed numeric Wind speed (m/s) at 2 m height.
#' @param net.irradiance numeric Long wave and short wave balance (W/m2).
#' @param net.radiation numeric Long wave and short wave balance (J/m2/day).
#' @param nighttime logical Used only for methods that distinguish between
#'   daytime- and nighttime canopy conductances.
#' @param soil.heat.flux numeric Soil heat flux (W/m2), positive if soil
#'   temperature is increasing.
#' @param method character The name of an estimation method.
#' @param check.range logical Flag indicating whether to check or not that
#'   arguments for temperature are within range of method. Passed to
#'   function calls to \code{water_vp_sat()} and \code{water_vp_sat_slope()}.
#'
#' @return A numeric vector of $ET_{ref}$ estimates expressed as mm/h.
#'
#' @references
#'   Allen R G, Pereira L S, Raes D, Smith M. 1998. Crop evapotranspiration:
#'   Guidelines for computing crop water requirements. Rome: FAO.
#'
#' @family Evapotranspiration and energy balance related functions.
#'
#' @examples
#' # instantaneous
#' ET_ref(temperature = 20,
#'        water.vp = water_RH2vp(0.7, 20),
#'        wind.speed = 0,
#'        net.irradiance = 10)
#'
#' ET_ref(temperature = c(5, 20, 35),
#'        water.vp = water_RH2vp(0.7, 20),
#'        wind.speed = 0,
#'        net.irradiance = 10)
#'
#' # Hot and dry air
#' ET_ref(temperature = 35,
#'        water.vp = water_RH2vp(0.1, 35),
#'        wind.speed = 5,
#'        net.irradiance = 400)
#'
#' ET_ref(temperature = 35,
#'        water.vp = water_RH2vp(0.10, 35),
#'        wind.speed = 5,
#'        net.irradiance = 400,
#'        method = "FAO.PM")
#'
#' ET_ref(temperature = 35,
#'        water.vp = water_RH2vp(0.1, 35),
#'        wind.speed = 5,
#'        net.irradiance = 400,
#'        method = "ASCE.PM.short")
#'
#' ET_ref(temperature = 35,
#'        water.vp = water_RH2vp(0.10, 35),
#'        wind.speed = 5,
#'        net.irradiance = 400,
#'        method = "ASCE.PM.tall")
#'
#' # Low temperature and high humidity
#' ET_ref(temperature = 5,
#'        water.vp = water_RH2vp(0.95, 5),
#'        wind.speed = 0.5,
#'        net.irradiance = -10,
#'        nighttime = TRUE,
#'        method = "ASCE.PM.short")
#'
#' ET_ref_day(temperature = 35,
#'            water.vp = water_RH2vp(0.10, 35),
#'            wind.speed = 5,
#'            net.radiation = 35e6) # 35 MJ / d / m2
#'
#' @export
#'
ET_ref <- function(temperature,
                   water.vp,
                   wind.speed,
                   net.irradiance,
                   nighttime = FALSE,
                   atmospheric.pressure = 1013e-2,
                   soil.heat.flux = 0,
                   method = "FAO.PM",
                   check.range = TRUE) {
  if (method == "FAO.PM") {
    vp.method <- "tetens"
    Cd <- 0.34
    Cn <- 37
    k <- 0.480
  } else if (method == "ASCE.PM.short") {
    vp.method <- "tetens"
    Cd <- ifelse(nighttime, 0.96, 0.24)
    Cn <- 37
    k <- 0.408
  } else if (method == "ASCE.PM.tall") {
    vp.method <- "tetens"
    Cd <- ifelse(nighttime, 1.7, 0.25)
    Cn <- 66
    k <- 0.408
  } else {
    warning("Method '", method,
            "' unavailable; use 'FAO.PM', 'ASCE.PM.short', or 'ASCE.PM.tall'.")
    return(NA_real_, length(temperature))
  }
  # slope of water vapour pressure curve (kPa C-1)
  delta <- water_vp_sat_slope(temperature,
                              over.ice = temperature <= -5,
                              method = vp.method,
                              check.range = check.range) * 1e-3 # Pa -> kPa
  # psychrometric constant (kPa C-1)
  gamma <- psychrometric_constant(atmospheric.pressure) * 1e-3 # Pa -> kPa
  # water vapour pressure deficit (kPa)
  vp.sat <- water_vp_sat(temperature,
                         over.ice = temperature <= -5,
                         method = vp.method,
                         check.range = check.range)
  vpd <- (vp.sat - water.vp) * 1e-3 # Pa -> kPa
  vpd <- ifelse(vpd < 0, 0, vpd)
  # net radiation
  radiation <- net.irradiance * 1e-6 * 3600 # W / m2 = J / m2 /s -> MJ / m2 / h
  # ET0
  (k * delta * (radiation - soil.heat.flux) +
      gamma * (Cn / (temperature + 273)) * wind.speed * vpd) /
    (delta + gamma * (1 + Cd * wind.speed))
}

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
#' @param Rs numeric Surface resistance (?).
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
#'       net.irradiance = 20)
#'
ET_PM <- function(temperature,
                  water.vp,
                  wind.speed,
                  net.irradiance,
                  lw.emissivity = 0.98,
                  Rs = 0,
                  atmospheric.pressure = 1013e-2,
                  method = "FMI.PM",
                  check.range = TRUE) {
  if (method == "FMI.PM") {
    wind.speed <- max(wind.speed, 0.5)
    sigma <- 5.670374419e-8 # Stefan–Boltzmann constant (W⋅m−2⋅)−
    rho <- 1.2923 # mean air density
    cp <- 1004 # specific heat
    l <- 2.5e6
    Rs <- 0
    k <- 421.0825 # 0.6 * ((log(8 / 0.0002) / 0.4)^2)
    Ra <-  k / wind.speed
    e.sat <- water_vp_sat(temperature,
                          over.ice = temperature <= -5,
                          check.range = check.range) * 10e-2 # Pa -> mbar
    es.slope <- water_vp_sat_slope(temperature,
                                   over.ice = temperature <= -5,
                                   check.range = check.range) * 10e-2  # Pa / K -> mbar / K
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

#' Net radiation flux
#'
#' Estimate net radiation balance expressed as a flux in W/m2. If
#' \code{lw.down.irradiance} is passed a value in W / m2 the difference is
#' computed directly and if not an approximate value is estimated, using
#' \code{R_rel = 0.75} which corresponds to clear sky, i.e., uncorrected for
#' cloudiness. This is the approach to estimation is that recommended by FAO for
#' hourly estimates while here we use it for instantaneous or mean flux rates.
#'
#' @param temperature numeric vector of air temperatures (C) at 2 m height.
#' @param sw.down.irradiance,lw.down.irradiance numeric Down-welling short wave
#'   and long wave radiation radiation (W/m2).
#' @param sw.albedo numeric Albedo as a fraction of one (/1).
#' @param lw.emissivity numeric Emissivity of the surface (ground or vegetation)
#'   for long wave radiation.
#' @param water.vp numeric vector of water vapour pressure in air (Pa), ignored
#'   if \code{lw.down.irradiance} is available.
#' @param R_rel numeric The ratio of actual and clear sky short wave irradiance
#'   (/1).
#'
#' @return A numeric vector of evapotranspiration estimates expressed as
#'   W / m-2.
#'
#' @export
#'
#' @family Evapotranspiration and energy balance related functions.
#'
net_irradiance <- function(temperature,
                           sw.down.irradiance,
                           lw.down.irradiance = NULL,
                           sw.albedo = 0.23,
                           lw.emissivity = 0.98,
                           water.vp = 0,
                           R_rel = 1) {
  stopifnot(all(is.na(temperature) | temperature >= -273.16))
  sigma <- 5.670374419e-8 # Stefan–Boltzmann constant (W⋅m−2 K)
  lw.up.irradiance <-
    sigma * (temperature + 273.16)^4 * lw.emissivity
  if (!is.null(lw.down.irradiance)) {
    net.lw.irradiance <- lw.down.irradiance - lw.up.irradiance
  } else {
    # we guess lw.down.irradiance
    net.lw.irradiance <-
      -lw.up.irradiance *
      (0.34 - 0.14 * sqrt(water.vp * 1e-3)) *
      (1.35 * R_rel - 0.35)
  }
  sw.down.irradiance * (1 - sw.albedo) + net.lw.irradiance
}

#' @rdname ET_ref
#'
#' @export
#'
ET_ref_day <- function(temperature,
                       water.vp,
                       wind.speed,
                       net.radiation,
                       nighttime = FALSE,
                       atmospheric.pressure = 1013e-2,
                       soil.heat.flux = 0,
                       method = "FAO.PM",
                       check.range = TRUE) {
  if (method == "FAO.PM") {
    vp.method <- "tetens"
    Cd <- 0.34
    Cn <- 900
    k <- 0.480
  } else {
    warning("Method '", method,
            "' unavailable; use 'FAO.PM'.")
    return(NA_real_, length(temperature))
  }
  # slope of water vapour pressure curve (kPa C-1)
  delta <- water_vp_sat_slope(temperature,
                              method = vp.method,
                              check.range = check.range) * 1e-3 # Pa -> kPa
  # psychrometric constant (kPa C-1)
  gamma <- psychrometric_constant(atmospheric.pressure) * 1e-3 # Pa -> kPa
  # water vapour pressure deficit (kPa)
  vp.sat <- water_vp_sat(temperature,
                         method = vp.method,
                         check.range = check.range)
  vpd <- (vp.sat - water.vp) * 1e-3 # Pa -> kPa
  vpd <- ifelse(vpd < 0, 0, vpd)
  # net radiation
  radiation <- net.radiation * 1e-6 # J / m2 / d -> MJ / m2 / d
  # ET0
  (k * delta * (radiation - soil.heat.flux) +
      gamma * (Cn / (temperature + 273)) * wind.speed * vpd) /
    (delta + gamma * (1 + Cd * wind.speed))
}



