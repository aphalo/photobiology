# names of all multi spectral classes -------------------------------------------

#' Function that returns a vector containing the names of multi-spectra classes.
#'
#' @export
#'
#' @return A \code{character} vector of class names.
#'
multi_spct_classes <- function() {
  c("generic_mspct", "cps_mspct",
    "filter_mspct", "reflector_mspct",
    "source_mspct", "object_mspct",
    "response_mspct", "chroma_mspct")
}

#' @title Constructors of multi_spct Objects
#'
#' @description Converts a list of spectral objects into a "multi spectrum"
#'   object by setting the class attibute of the list of spectra to the
#'   corresponding multi-spct class, check that components of the list belong to
#'   the expected class.
#'
#' @param l list of generic_spct or derived classes
#' @param class character The class expected for the elements of l
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#' @param ... additional arguments
#'
#' @export
#' @exportClass generic_mspct
#'
#' @family collections of spectra classes family
#'
generic_mspct <- function(l, class = "generic_spct", ncol = 1, byrow = FALSE,  ...) {
  stopifnot(is.list(l))

  multi_class <- paste(sub("_spct", "_mspct", class))
  if (class(l)[1] == multi_class) {
    warning("Class is already set to '", multi_class, "'")
  } else {
    stopifnot(class %in% spct_classes())
    for (spct in l) {
      stopifnot(class %in% class_spct(spct))
    }
    setattr(l, "class", c(multi_class,
                          ifelse(multi_class != "generic_mspct", "generic_mspct", NULL),
                          "list"))
  }
  setattr(l, "mspct.version", 1)

  setattr(l, "ncol", ncol)
  setattr(l, "byrow", byrow)
  l
}

#' @describeIn generic_mspct Specialization for collections of \code{cps_spct} objects.
#'
#' @export
#' @exportClass cps_mspct
#'
cps_mspct <- function(l, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "cps_spct", ncol = ncol, byrow = byrow, ...)
}

#' @describeIn generic_mspct Specialization for collections of \code{source_spct} objects.
#'
#' @export
#' @exportClass source_mspct
#'
source_mspct <- function(l, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "source_spct", ncol = ncol, byrow = byrow, ...)
}

#' @describeIn generic_mspct Specialization for collections of \code{filter_spct} objects.
#'
#' @export
#' @exportClass filter_mspct
#'
filter_mspct <- function(l, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "filter_spct", ncol = ncol, byrow = byrow, ...)
}

#' @describeIn generic_mspct Specialization for collections of \code{reflector_spct} objects.
#'
#' @export
#' @exportClass reflector_mspct
#'
reflector_mspct <- function(l, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "reflector_spct", ncol = ncol, byrow = byrow, ...)
}

#' @describeIn generic_mspct Specialization for collections of \code{object_spct} objects.
#'
#' @export
#' @exportClass object_mspct
#'
object_mspct <- function(l, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "object_spct", ncol = ncol, byrow = byrow, ...)
}

#' @describeIn generic_mspct Specialization for collections of \code{response_spct} objects.
#'
#' @export
#' @exportClass response_mspct
#'
response_mspct <- function(l, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "response_spct", ncol = ncol, byrow = byrow, ...)
}

#' @describeIn generic_mspct Specialization for collections of \code{chroma_spct} objects.
#'
#' @export
#' @exportClass chroma_mspct
#'
chroma_mspct <- function(l, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "chroma_spct", ncol = ncol, byrow = byrow, ...)
}
