#' Calculate multipliers by interpolation from filter data
#'
#' @description
#' Calculate multipliers by interpolation from filter data, from user-supplied spectral transmittance
#' data or by name for data included in the package
#'
#' @usage calc_filter_multipliers(w.length.out, filter=clear.spct,
#'                                w.length.in=NULL, transmittance.in=NULL,
#'                                pc.in=TRUE, pc.out=FALSE, div=1.0)
#'
#' @param w.length.out numeric vector of wavelengths (nm) for output
#' @param filter a character string giving the name of a filter data set, or an object of class "filter.spct"
#'   default is 'clear.spct' a clear filter (T = 1.0)
#' @param w.length.in numeric vector of wavelengths (nm) for input
#' @param transmittance.in numeric vector of spectral transmittance value (as percent)
#' @param pc.in logical value indicating whether transmittances are expressed as percentages or fractions (default is to receive a percent)
#' @param pc.out logical value indicating whether transmittances are expressed as percentages or fractions (default is to return a fraction)
#' @param div numeric value default is 100.0 if pc=TRUE, but if pc=FALSE, default is 1.0, but can be changed in the call
#'
#' @return a numeric vector of fractional transmittances
#' @keywords manip misc
#' @export
#' @examples
#' require(photobiologyFilters)
#' data(polythene.new.spct)
#' with(polythene.new.spct, calc_filter_multipliers(400:500, w.length.in=w.length, transmittance.in=Tfr))
#' calc_filter_multipliers(400:500, "polythene.new")
#' calc_filter_multipliers(400:500, polythene.new.spct)
#' calc_filter_multipliers(400:500, "polythene.new", pc.out=TRUE)
#' calc_filter_multipliers(400:500)
#'
calc_filter_multipliers <- function(w.length.out,
                                    filter=clear.spct,
                                    w.length.in=NULL, transmittance.in=NULL,
                                    pc.in=TRUE,
                                    pc.out=FALSE, div = 1.0) {
  if (is.null(w.length.in) | is.null(transmittance.in)) {
    if (is.character(filter)) {
      filter.name <- filter
      filter.object.name <- paste(filter.name, "spct", sep=".")
      if (!exists(filter.object.name)) {
        filter.object.name <- paste(filter.name, "data", sep=".")
        if (!exists(filter.object.name)) {
          warning("No data for filter with name: ", filter.object.name)
          return(NA)
        }
      }
      filter.object <- get(filter.object.name)
    } else {
      filter.object <- filter
    }
    if (!is(filter.object, "filter.spct")) {
      warning('Object is not of class "filter.spct".')
      return(NA)
    }
    if (pc.out) div <- div / 100.0 # we use fractional Tfr from filter object
    return(with(filter.object, spline(w.length, Tfr, xout=w.length.out)$y) / div)
  } else {
    if (pc.in!=pc.out){
      if (pc.in ) { # & !pc.out
        div <- div * 100.0
      } else {  # !pc,in & pc.out
        div <- div / 100.0
      }
    }
    if (!check_spectrum(w.length.in, transmittance.in)) {
      return(NA)
    }
    return(spline(w.length.in, transmittance.in, xout=w.length.out)$y / div)
  }
}
