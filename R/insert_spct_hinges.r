##' Insert new wavelengths into a spectrum, interpolating the corresponding spectral data values.
##'
##' Inserting wavelengths values immediately before and after a discontinuity in the SWF,
##' greatly reduces the errors caused by interpolating the weighted irradiance during
##' integration of the effective spectral irradiance. This is specially true when data
##' has a large wavelength step size.
##'
##' @usage insert_spct_hinges(spct, hinges=NULL)
##'
##' @param spct an object of class "generic.spct"
##' @param hinges a numeric array giving the wavelengths (nm) at which the s.irrad should be inserted by
##' interpolation, no interpolation is indicated by an empty array (numeric(0))
##'
##' @return a data.frame with variables \code{w.length} and \code{s.irrad}
##' @keywords manip misc
##' @export
##' @examples
##' data(sun.spct)
##' insert_spct_hinges(sun.spct, c(399.99,400.00,699.99,700.00))
##' insert_spct_hinges(sun.spct, c(199.99,200.00,399.50,399.99,400.00,699.99,700.00,799.99,1000.00))
insert_spct_hinges <- function(spct, hinges=NULL) {
  hinges <- hinges[hinges > min(spct) & hinges < max(spct)]
  if (length(hinges) > 0) {
    names.spct <- names(spct)
    names.data <- names.spct[names.spct != "w.length"]
    class.spct <- class(spct)
    comment.spct <- comment(spct)
    hinges <- unique(sort(hinges))
    first.iter <- TRUE
    for (data.col in names.data) {
      temp.data <- insert_hinges(spct[["w.length"]], spct[[eval(data.col)]], hinges)
      if (first.iter) {
        new.spct <- data.table(w.length = temp.data[["w.length"]])
        first.iter <- FALSE
      }
      new.spct[ , eval(data.col) := temp.data[["s.irrad"]] ]
    }
    setattr(new.spct, "comment", comment.spct)
    setattr(new.spct, "class", class.spct)
    return(new.spct)
  } else {
    return(spct)
  }
}
