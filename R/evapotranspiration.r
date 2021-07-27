#' Evapotranspiration
#'
#' Compute an estimate of actual or potential evapotranspiration from
#' meteorologial data. Evapotranspiration from vegetation includes
#' transpiraction by plants plus evaporation from the soil or other wet
#' surfaces. $ET_{ref}$ is the reference value assuming no limitation to
#' transpiration due to soil water, similar to potential evapotranspiration
#' (PET). An actual value $ET_{ref}$ can be estimated only if additional
#' information on the plants and soil is available.
#'
#' @details Currently only one method, based on the Penman-Monteith equation
#'   formulated as recommended by FAO (Allen et al., 1998), is implemented. It
#'   relies on data meassured according WHO standards at 2 m above ground level
#'   to estimate reference evapotranspiration ($ET_{ref}$).
#'   Corrections are included to account for the assumption that the vegetation
#'   has a height of 0.12 m. These input data can be averages over different
#'   lengths of time. Allen et al. (1998) describes different approaches to
#'   estimate input data when not available. The auxiliary function
#'   \code{net_lw_radiation()} can be used to estimate the net lw radiation
#'   uncorrected for cloudiness.
#'
#' @param temperature numeric vector of air temperatures (C) at 2 m height.
#' @param water.vp numeric vector of water vapour pressure in air (Pa).
#' @param atmospheric.pressure numeric Atmospheric pressure (Pa).
#' @param wind.speed numeric Wind speed (m/s) at 2 m height.
#' @param radiation numeric Global radiation (W/m2).
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
#' @examples
#' # instantaneous
#' ET_ref(temperature = 20,
#'        water.vp = water_RH2vp(0.7, 20),
#'        wind.speed = 0,
#'        radiation = 10)
#'
#' ET_ref(temperature = 35,
#'        water.vp = water_RH2vp(0.1, 35),
#'        wind.speed = 5,
#'        radiation = 400)
#'
#' ET_ref(temperature = 35,
#'        water.vp = water_RH2vp(0.1, 35),
#'        wind.speed = 5,
#'        radiation = 400,
#'        method = "ASCE.PM.short")
#'
#' ET_ref(temperature = 5,
#'        water.vp = water_RH2vp(0.95, 5),
#'        wind.speed = 0.5,
#'        radiation = -10,
#'        nighttime = TRUE,
#'        method = "ASCE.PM.short")
#'
#' ET_ref(temperature = 35,
#'        water.vp = water_RH2vp(0.10, 35),
#'        wind.speed = 5,
#'        radiation = 400,
#'        method = "ASCE.PM.tall")
#'
#' @export
#'
ET_ref <- function(temperature,
                   water.vp,
                   wind.speed,
                   radiation,
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
                              method = vp.method,
                              check.range = check.range) * 1e-3
  # psychrometric constant (kPa C-1)
  gamma <- psychrometric_constant(atmospheric.pressure) * 1e-3
  # water vapour pressure deficit (kPa)
  vpd <- water_vp_sat(temperature,
                      method = vp.method,
                      check.range = check.range) - water.vp
  vpd <- ifelse(vpd < 0, 0, water_vp_sat(temperature) - water.vp) * 1e-3
  # net radiation
  radiation <- radiation * 1e-6 * 3600 # W / m2 = J / m2 /s -> Mj / m2 / h
  # ET0
  (k * delta * (radiation - soil.heat.flux) +
      gamma * (Cn / (temperature + 273)) * wind.speed * vpd) /
    (delta + gamma * (1 + Cd * wind.speed))
}

#' @rdname ET_ref
#'
#' @return A numeric vector of $R_{nl}$ estimates expressed as W / m-2.
#'
#' @export
#'
net_lw_radiation <- function(temperature,
                             water.vp) {
  sigma <- 4.903e-9 / 86400

  sigma * (temperature + 273.16)^4 * (0.34 - 0.14 * sqrt(water.vp * 1e-3))
}
