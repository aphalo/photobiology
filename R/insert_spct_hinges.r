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
##' @return a generic.spct or a derived type with variables \code{w.length} and other numeric variables.
##' @keywords manip misc
##' @export
##' @examples
##' data(sun.spct)
##' insert_spct_hinges(sun.spct, c(399.99,400.00,699.99,700.00))
##' insert_spct_hinges(sun.spct,
##'                    c(199.99,200.00,399.50,399.99,400.00,699.99,
##'                          700.00,799.99,1000.00))
insert_spct_hinges <- function(spct, hinges=NULL) {
  if (is.null(hinges)) {
    return(spct)
  }
  hinges <- hinges[hinges > min(spct) & hinges < max(spct)]
  hinges <- unique(sort(hinges))
  old.w.length <- spct[["w.length"]]
  hinges <- setdiff(hinges, old.w.length)
  if (length(hinges) > 0) {
    new.w.length <- unique(sort(c(hinges, old.w.length)))
    name <- substitute(spct)
    names.spct <- names(spct)
    names.data <- names.spct != "w.length"
    idx.wl <- which(!names.data)
    idx.data <- which(names.data)
    class.spct <- class(spct)
    comment.spct <- comment(spct)
    if  (is(spct, "source.spct")) {
      time.unit.spct <- attr(spct, "time.unit", exact=TRUE)
    }
    new.spct <- data.table(w.length = new.w.length)
    first.iter <- TRUE
    for (data.col in idx.data) {
      temp.data <- spct[[data.col]]
      if (is.numeric(temp.data)) {
        new.spct[ , names.spct[data.col] := put_hinges(old.w.length, temp.data, hinges)]
      } else {
        new.spct[ , names.spct[data.col] := NA]
      }
    }
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
    return(new.spct)
  } else {
    return(spct)
  }
}
