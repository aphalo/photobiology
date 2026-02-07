#' Simulate light sensor response
#'
#' @param source.spct source_spct Light source spectral irradiance.
#' @param sensor.mspct response_spct or response_mspct Light sensor
#'   spectral responsiveness.
#' @param range numeric vector of length two Range of wavelengths to use
#'   (nanometres, nm)
#' @inheritParams response
#' @inheritParams normalize
#'
#' @details Compute sensor response spectrum by convolution of light source
#'   emission spectrum or illumination spectrum and the responsiveness spectrum
#'   of a sensor with one or more channels. Return the integral over
#'   wavelengths' \code{range} for each sensor channel, or sensor in
#'   \code{sensor.mspct}.
#'
#' @return A data frame
#'
#' @export
#'
#' @import photobiology
#'
#' @examples
#' simul_sensor_response(sun.spct, ccd.spct)
#' simul_sensor_response(sun.spct, ccd.spct, range = c(400, 700))
#' simul_sensor_response(sun.spct, ccd.spct, unit.out = "photon")
#'
simul_sensor_response <-
  function(source.spct,
           sensor.mspct,
           norm = "skip",
           range = NULL,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           time.unit = NULL,
           scale.factor = 1,
           attr2tb = NULL) {
    if (!is.null(range)) {
      range <- range(range)
    }
    # handle spectra in long form individually and as members of a collection
    sensor.mspct <- subset2mspct(sensor.mspct)
    sensor.mspct <-
      normalise(sensor.mspct, range = range, norm = norm)
    stopifnot("'source.spct' wrong class" =
                is.source_spct(source.spct))
    stopifnot("'sensor.mspct' wrong class" =
                is.response_mspct(sensor.mspct))
    source.spct <- trim_wl(source.spct, range = range, fill = 0)
    sensor.mspct <- trim_wl(sensor.mspct, range = range, fill = 0)

    channel.responses.mspct <- response_mspct()

    # convolution
    for (ch in names(sensor.mspct)) {
      channel.responses.mspct[[ch]] <- sensor.mspct[[ch]] * source.spct
    }
    # integration
    z <-
      response(spct = channel.responses.mspct,
               w.band = range,
               unit.out = unit.out,
               quantity = "total",
               time.unit = time.unit,
               scale.factor = scale.factor,
               attr2tb = attr2tb)

    comment(z) <-
      paste("Simulated light sensor response.",
            "\nillumination: ", what_measured(source.spct),
            "\nsensor: ", what_measured(sensor.mspct), sep = "")
    z
  }
