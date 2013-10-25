#' Calculate spectral values at a different set of wavelengths
#'
#' For example an spectrum [W m-2 nm-1] is converted into a spectrum [mol s-1 m-2 nm-1]
#'
#' @usage interpolate_spectrum(w.length.in, s.irrad, w.length.out, fill.value=NA)
#' 
#' @param w.length.in numeric array of wavelengths (nm)
#' @param s.irrad a numeric array of spectral values
#' @param w.length.out numeric array of wavelengths (nm)
#' @param fill.value a value to be assigned to out of range wavelengths
#'
#' @return a numeric array of interpolated spectral values
#' 
#' @export
#' @keywords manip misc
#' @examples
#' data(sun.data)
#' my.w.length <- 300:700
#' my.s.e.irrad <- with(sun.data, interpolate_spectrum(w.length, s.e.irrad, my.w.length))
#' plot(my.s.e.irrad ~ my.w.length)
#' lines(s.e.irrad ~ w.length, data=sun.data)

interpolate_spectrum <- function(w.length.in, s.irrad, w.length.out, fill.value=NA) {
  if (is.null(fill.value) && (w.length.out[1] < w.length.in[1] || 
                                w.length.out[length(w.length.out)] > w.length.in[length(w.length.in)])) {
    stop("Extrapolation attempted with fill.value=NULL")
  }
  selector <- w.length.out >= w.length.in[1] & w.length.out <= w.length.in[length(w.length.in)]
  s.irrad.out <- numeric(length(w.length.out))
  if (!is.null(fill.value)){
    s.irrad.out[!selector] <- fill.value
  }
  s.irrad.out[selector] <- spline(w.length.in, s.irrad, xout=w.length.out[selector])$y
  return(s.irrad.out)
}