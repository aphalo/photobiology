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
  if (is.null(hinges)) {
    return(spct)
  }
  hinges <- hinges[hinges > min(spct) & hinges < max(spct)]
  if (length(hinges) > 0) {
    name <- substitute(spct)
    names.spct <- names(spct)
    names.data <- names.spct != "w.length"
    idx.wl <- which(!names.data)
    idx.data <- which(names.data)
    class.spct <- class(spct)
    comment.spct <- comment(spct)
    if  (is(x, "source.spct")) {
      time.unit.spct <- attr(spct, "time.unit", exact=TRUE)
    }
    hinges <- unique(sort(hinges))
    first.iter <- TRUE
#    setDF(spct)
    for (data.col in idx.data) {
      temp.data <- insert_hinges(spct[["w.length"]], spct[[data.col]], hinges)
      if (first.iter) {
        new.spct <- data.table(w.length = temp.data[[1]])
        first.iter <- FALSE
      }
      new.spct[ , names.spct[data.col] := temp.data[[2]] ]
    }
    setattr(new.spct, "comment", comment.spct)
    if(class.spct[1] == "source.spct") {
      setSourceSpct(new.spct)
      if (!is.null(time.unit.spct)) {
        setattr(new.spct, "unit.out", time.unit.spct)
      }
    } else if (class.spct[1] == "filter.spct") {
      setFilterSpct(new.spct)
    } else if (class.spct[1] == "reflector.spct") {
      setReflectorSpct(new.spct)
    } else if (class.spct[1] == "response.spct") {
      setResponseSpct(new.spct)
    } else if (class.spct[1] == "chroma.spct") {
      setChromaSpct(new.spct)
    } else if (class.spct[1] == "generic.spct") {
      setGenericSpct(new.spct)
    }
    setattr(new.spct, "comment", comment.spct)
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, new.spct, parent.frame(), inherits = TRUE)
    }
    invisible(return(new.spct))
  } else {
    invisible(return(spct))
  }
}
