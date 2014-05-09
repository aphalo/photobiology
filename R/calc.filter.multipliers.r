#' Calculate multipliers by interpolation from filter data
#' 
#' @description
#' Calculate multipliers by interpolation from filter data, from user-supplied spectral transmittance 
#' data or by name for data included in the package
#'
#' @usage calc_filter_multipliers(w.length.out, filter.name="clear", 
#'                                w.length.in=NULL, transmittance.in=NULL, 
#'                                pc=FALSE, div=1.0)
#' 
#' @param w.length.out numeric vector of wavelengths (nm) for output
#' @param filter.name a character string giving the name of a filter data set, default is a 'clear' filter (100 \% transmittance)
#' @param w.length.in numeric vector of wavelengths (nm) for input
#' @param transmittance.in numeric vector of spectral transmittance value (as percent)
#' @param pc logical value indicating whether transmittances are expressed as percentages or fractions (default is to return a fraction)
#' @param div numeric value default is 100.0 if pc=TRUE, but if pc=FALSE, default is 1.0, but can be changed in the call
#'  
#' @return a numeric vector of fractional transmittances 
#' @keywords manip misc
#' @export
#' @examples
#' require(photobiologyFilters)
#' data(polythene.new.data)
#' with(polythene.new.data, calc_filter_multipliers(400:500, w.length.in=w.length, transmittance.in=transmittance))
#' calc_filter_multipliers(400:500, "polythene.new")
#' calc_filter_multipliers(400:500, "polythene.new", pc=TRUE)
#' calc_filter_multipliers(400:500)
#' 
calc_filter_multipliers <- function(w.length.out,
                                    filter.name="clear", 
                                    w.length.in=NULL, transmittance.in=NULL, 
                                    pc=FALSE, div = 1.0) {
  if (!pc) div <- div * 100.0 # filter transmittance data is assumed to be stored as percentages
  if (is.null(w.length.in) | is.null(transmittance.in)) {
    filter.object <- paste(filter.name, "data", sep=".")
    return(with(get(filter.object), spline(w.length, transmittance, xout=w.length.out)$y) / div)
  } else {
    if (!check_spectrum(w.length.in, transmittance.in)) {
      return(NA)
    }
    return(spline(w.length.in,transmittance.in,xout=w.length.out)$y / div)
  }
}