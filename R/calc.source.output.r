#' Calculate light source output by interpolation from lamp data
#' 
#' @description
#' Calculate interpolated values by interpolation from user-supplied spectral emission
#' data or by name for light source data included in the packages photobiologySun,
#' photobiologyLamps, or photobiologyLEDs, scaling the values.
#'
#' @usage calc_source_output(w.length.out, source.name=NULL, 
#'                                w.length.in=NULL, s.irrad.in=NULL, 
#'                                unit.in="energy",
#'                                scaled=NULL, fill=NA)
#' 
#' @param w.length.out numeric vector of wavelengths (nm) for output
#' @param source.name a character string giving the name of a lamp data set, default is NULL
#' @param w.length.in numeric vector of wavelengths (nm) for input
#' @param s.irrad.in numeric vector of spectral transmittance value (fractions or percent)
#' @param unit.in a character string "energy" or "photon"
#' @param scaled NULL, "peak", "area"; div ignored if !is.null(scaled) 
#' @param fill if NA, no extrapolation is done, and NA is returned for wavelengths outside the range of the input. If NULL then the tails are deleted. If 0 then the tails are set to zero.
#'  
#' @return a dataframe with four numeric vectors with wavelength values (w.length), scaled and interpolated spectral energy irradiance (s.e.irrad), 
#'    scaled and interpolated spectral photon irradiance values (s.q.irrad).
#' @keywords manip misc
#' @export
#' 
#' @note This is a convenience function that adds no new functionality but makes it a little easier to plot lamp spectral emission data consistently.
#' It automates interpolation, extrapolation/trimming and scaling.
#' @examples
#' with(sun.data, calc_source_output(290:1100, w.length.in=w.length, s.irrad.in=s.e.irrad))
#' calc_source_output(290:1100, w.length.in=sun.data$w.length, s.irrad.in=sun.data$s.e.irrad, fill=0.0)
#' calc_source_output(200:4000, "sun")
#' calc_source_output(500:600, "sun")
#' 
calc_source_output <- function(w.length.out,
                               source.name=NULL, 
                               w.length.in=NULL, s.irrad.in=NULL, 
                               unit.in="energy",
                               scaled=NULL, fill=NA) {  
  # we first check the different possible inputs and convert to
  # two vectors w.length.in and s.irrad.in
  
  if (is.null(w.length.in) | is.null(s.irrad.in)) {
    if (is.null(source.name)) {
      return(NA) 
    } else {
      lamp.object.name <- paste(source.name, "dt", sep=".")
      if (!exists(lamp.object.name)) {
        lamp.object.name <- paste(source.name, "data", sep=".")
        if (!exists(lamp.object.name)) {
          warning("No data for lamp with name: ", lamp.object.name)
          return(NA)
        }
      }
      lamp.object <- get(lamp.object.name)
      w.length.in <- lamp.object$w.length
      if (with(lamp.object, exists("s.e.irrad"))) {
        unit.in <- "energy"
        s.irrad.in <- lamp.object$s.e.irrad
      } else if (with(lamp.object, exists("s.q.irrad"))) {
        unit.in <- "photon"
        s.irrad.in <- lamp.object$s.q.irrad
      } else {
        return(NA)
      }
    }
  } else if (!check_spectrum(w.length.in, s.irrad.in)) {
      return(NA)
  }
  
  # we interpolate using a spline or linear interpolation
  out.fill.selector <- w.length.out < w.length.in[1] | w.length.out > w.length.in[length(w.length.in)]
  if (is.null(fill)) {
    w.length.out <- w.length.out[!out.fill.selector]
    out.fill.selector <- rep(FALSE, length(w.length.out))                  
  }
  s.irrad.out <- numeric(length(w.length.out))
  
  if (length(w.length.out) < 25) {
    # cubic spline
    s.irrad.out[!out.fill.selector] <- spline(w.length.in, s.irrad.in, xout=w.length.out[!out.fill.selector])$y
  } else {
    # linear interpolation
    s.irrad.out[!out.fill.selector] <- approx(x = w.length.in, y = s.irrad.in, xout = w.length.out[!out.fill.selector], ties = "ordered")$y
  }

  # we check unit.in and and convert the output spectrum accordingly
  
  if (unit.in == "energy") {
    out.data <- data.frame(w.length = w.length.out, s.e.irrad = s.irrad.out, s.q.irrad = as_quantum_mol(w.length.out, s.irrad.out))
    } else if (unit.in == "photon") {
    out.data <- data.frame(w.length = w.length.out, s.e.irrad = as_energy(w.length.out, s.irrad.out), s.q.irrad = s.irrad.out)
    } else {
    warning("Bad argument for unit.in: ", unit.in)
    return(NA)
  }
  
  # do scaling
  
  if (!is.null(scaled)) {
    if (scaled=="peak") {
      e.div <- max(out.data$s.e.irrad, na.rm=TRUE)
      q.div <- max(out.data$s.q.irrad, na.rm=TRUE)
    } else if (scaled=="area") {
      s.irrad.na.sub <- out.data$s.e.irrad
      s.irrad.na.sub[is.na(s.irrad.na.sub)] <- 0.0
      e.div <- integrate_irradiance(w.length.out, s.irrad.na.sub)
      s.irrad.na.sub <- out.data$s.q.irrad
      s.irrad.na.sub[is.na(s.irrad.na.sub)] <- 0.0
      q.div <- integrate_irradiance(w.length.out, s.irrad.na.sub)
    } else {
      warning("Ignoring unsupported scaled argument: ", scaled)
      e.div <- q.div <- 1.0
    }
    out.data$s.e.irrad[!out.fill.selector] <- out.data$s.e.irrad[!out.fill.selector] / e.div
    out.data$s.q.irrad[!out.fill.selector] <- out.data$s.q.irrad[!out.fill.selector] / q.div
  }
  out.data$s.e.irrad[out.fill.selector] <- fill
  out.data$s.q.irrad[out.fill.selector] <- fill
  
  return(out.data)
}