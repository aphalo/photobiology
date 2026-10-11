# convolution ---------------------------------------------------------------

#' Convolve function for collections of spectra
#'
#' Convolve function for collections of spectra which applies an operation on
#' all the individual members of the collection(s) of spectra.
#'
#' @param e1 an object of class \code{generic_mspct} or \code{generic_scpt} or
#'   \code{numeric}
#' @param e2 an object of class \code{generic_mspct} or \code{generic_scpt} or
#'   \code{numeric}
#' @param oper function, usually but not necessarily an operator with two
#'   arguments.
#' @param sep character Used when pasting the names of members of \code{e1} and
#'   \code{e2} to form the names of members of the returned collection of
#'   spectra.
#' @param ... additional arguments passed to \code{oper} if present.
#'
#' @note At least one of e1 and e2 must be a \code{generic_mspct} object or
#'   derived.
#'
#' @export
#'
#' @family math operators and functions
#'
convolve_each <- function(e1, e2, oper = `*`, sep = "_", ...) {
  e3 <- list()
  if (is.any_mspct(e1) && !is.any_mspct(e2)) {
    for (spct.name in names(e1)) {
      e3[[spct.name]] <- oper(e1[[spct.name]], e2, ...)
    }
    z <- generic_mspct(e3, class = shared_member_class(e3),
                       ncol = ncol(e1),
                       byrow = attr(e1, "mspct.byrow", exact = TRUE))
  } else if (!is.any_mspct(e1) && is.any_mspct(e2)) {
    for (spct.name in names(e2)) {
      e3[[spct.name]] <- oper(e1, e2[[spct.name]], ...)
    }
    z <- generic_mspct(e3, class = shared_member_class(e3),
                       ncol = ncol(e2),
                       byrow = attr(e2, "mspct.byrow", exact = TRUE))
  } else if (is.any_mspct(e1) && is.any_mspct(e2)) {
    for (spct.name1 in names(e1)) {
      for (spct.name2 in names(e2)) {
        combined.name <- paste(spct.name1, spct.name2, sep = sep)
        e3[[combined.name]] <- oper(e1[[spct.name1]], e2[[spct.name2]], ...)
      }
    }
    z <- generic_mspct(e3, class = shared_member_class(e3),
                       ncol = nrow(e2),
                       byrow = FALSE)
    attr(z, "mspct.dimnames")  <- list(names(e1), names(e2))
  } else {
    stop("At least one of 'e1' and 'e2' should be a collection of spectra.")
  }
  z
}
