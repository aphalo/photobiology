#' Generic function
#'
#' Remove tags from an R object containing the expected data members.
#'
#' @param x an R object
#' @param ... not used in current version
#' @export untag
untag <- function(x, ...) UseMethod("untag")

#' Default for generic function
#'
#' Remove tags from an R object containing the expected data members.
#'
#' @param x an R object
#' @param ... not used in current version
#' @export untag.default
untag.default <- function(x, ...) {
  return(x)
}

#' Specialization for generic.spct
#'
#' Remove tags from a generic.spct object.
#'
#' @param x a generic.spct object
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @param ... not used in current version
#'
#' @return if \code{x} contains tag data they are removed and the "spct.tags" atrribute is set to \code{NA},
#' while if \code{x} has no tags, it is not modified. In either case, the byref argument is respected:
#' in all cases if \code{byref=FALSE} a copy of \code{x} is returned.
#'
#' @export untag.generic.spct
#'
untag.generic.spct <- function(x,
                               byref=TRUE, ...) {
  if (!byref) {
    x <- copy(x)
    name <- NA
  } else {
    name <- substitute(x)
  }
  if (!is.tagged(x)) {
    return(x)
  }
  x[ , wl.color := NULL]
  x[ , wb.f := NULL ]
  tag.data <- NA
  setattr(x, "spct.tags", tag.data)
  # to work by reference we need to assign the new DT to the old one
  if (byref & is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}