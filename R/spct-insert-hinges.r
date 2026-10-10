#' Insert new wavelength values into a spectrum
#'
#' Insert new wavelength values into a spectrum interpolating the corresponding
#' spectral data values.
#'
#' @param spct an object of class "generic_spct"
#' @param hinges numeric vector of wavelengths (nm) at which the
#'   s.irrad should be inserted by interpolation, no interpolation is indicated
#'   by an empty vector (numeric(0))
#' @param byref logical indicating if new object will be created by reference or
#'   by copy of spct
#'
#' @return a generic_spct or a derived type with variables \code{w.length} and
#'   other numeric variables.
#'
#' @note Inserting wavelengths values "hinges" immediately before and after a
#'   discontinuity in the SWF, greatly reduces the errors caused by
#'   interpolating the weighted irradiance during integration of the effective
#'   spectral irradiance. This is specially true when data has a large
#'   wavelength step size.
#'
#' @export
#'
#' @examples
#'
#' insert_spct_hinges(sun.spct, c(399.99,400.00,699.99,700.00))
#' insert_spct_hinges(sun.spct,
#'                    c(199.99,200.00,399.50,399.99,400.00,699.99,
#'                          700.00,799.99,1000.00))
insert_spct_hinges <- function(spct, hinges=NULL, byref = FALSE) {
  if (!is.generic_spct(spct)) {
    warning("Only objects derived from 'generic_spct' are supported")
    return(spct)
  }
  if (is.null(hinges) || length(hinges) == 0) {
    return(spct)
  } else {
    name <- substitute(spct)
    old.w.length <- spct[["w.length"]]

    colnames.spct <- names(sapply(spct, is.numeric))
    if (length(colnames.spct) != ncol(spct)) {
      warning("Dropping non-numeric columns: ", setdiff(colnames(spct), colnames.spct))
    }

    idx.data <- which(colnames.spct != "w.length")

    first.iter <- TRUE
    for (data.col in idx.data) {
      temp.data <- spct[[data.col]]
      if (first.iter) {
        new.spct <- insert_hinges(old.w.length, temp.data, hinges)
        names(new.spct) <- c("w.length", colnames.spct[data.col])
        first.iter <- FALSE
      } else {
        new.spct[ , colnames.spct[data.col] ] <-
          v_insert_hinges(old.w.length, temp.data, hinges)
      }
    }
    new.spct <- copy_attributes(spct, new.spct, copy.class = TRUE)
    if (byref && is.name(name)) {
      name <- as.character(name)
      assign(name, new.spct, parent.frame(), inherits = TRUE)
    }
    return(new.spct)
  }
}
