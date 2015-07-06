# names of all multi spectral classes -------------------------------------------

#' Function that returns a vector containing the names of multi-spectra classes.
#'
#' @export
#'
#' @return A \code{character} vector of class names.
#'
multi_spct_classes <- function() {
  c("generic_multi_spct", "cps_multi_spct",
    "filter_multi_spct", "reflector_multi_spct",
    "source_multi_spct", "object_multi_spct",
    "response_multi_spct", "chroma_multi_spct")
}

#' @title Convert a list of spectral objects into a "multi spectrum" object
#'
#' @description Sets the class attibute of a list of spectra to the
#'   corresponding multi-spct class, also setting the target class for summaries.
#'
#' @param l list of generic_spct or derived classes
#' @param class the class expected for the elements of l
#' @param ncol integer Number of 'virtual' columns in data
#'
#' @export
#' @exportClass generic_multi_spct
#' @exportClass cps_multi_spct
#' @exportClass source_multi_spct
#' @exportClass response_multi_spct
#' @exportClass filter_multi_spct
#' @exportClass reflector_multi_spct
#' @exportClass object_multi_spct
#' @exportClass response_multi_spct
#' @exportClass chroma_multi_spct
#'
#' @family collections of spectra classes family
#'
generic_multi_spct <- function(l, class = "generic_spct", ncol = 1, byrow = FALSE,  ...) {
  stopifnot(is.list(l))

  multi_class <- paste(sub("_spct", "_multi_spct", class))
  if (class(l)[1] == multi_class) {
    warning("Class is already set to '", multi_class, "'")
  } else {
    stopifnot(class %in% spct_classes())
    for (spct in l) {
      stopifnot(class %in% class_spct(spct))
    }
    setattr(l, "class", c(multi_class,
                          ifelse(multi_class != "generic_multi_spct", "generic_multi_spct", NULL),
                          "list"))
  }
  setattr(l, "ncol", ncol)
  setattr(l, "byrow", byrow)
  l
}

#' @describeIn generic_multi_spct Specialization for collections of \code{cps_spct} objects.
#'
#' @export
#'
cps_multi_spct <- function(l, ncol = 1, byrow = FALSE, idfactor = TRUE, ...) {
  generic_multi_spct(l, class = "cps_spct", ncol = ncol, byrow = byrow, ...)
}

#' @describeIn generic_multi_spct Specialization for collections of \code{source_spct} objects.
#'
#' @export
#'
source_multi_spct <- function(l, ncol = 1, byrow = FALSE, idfactor = TRUE, ...) {
  generic_multi_spct(l, class = "source_spct", ncol = ncol, byrow = byrow, ...)
}

#' @describeIn generic_multi_spct Specialization for collections of \code{filter_spct} objects.
#'
#' @export
#'
filter_multi_spct <- function(l, ncol = 1, byrow = FALSE, idfactor = TRUE, ...) {
  generic_multi_spct(l, class = "filter_spct", ncol = ncol, byrow = byrow, ...)
}

#' @describeIn generic_multi_spct Specialization for collections of \code{reflector_spct} objects.
#'
#' @export
#'
reflector_multi_spct <- function(l, ncol = 1, byrow = FALSE, idfactor = TRUE, ...) {
  generic_multi_spct(l, class = "reflector_spct", ncol = ncol, byrow = byrow, ...)
}

#' @describeIn generic_multi_spct Specialization for collections of \code{object_spct} objects.
#'
#' @export
#'
object_multi_spct <- function(l, ncol = 1, byrow = FALSE, idfactor = TRUE, ...) {
  generic_multi_spct(l, class = "object_spct", ncol = ncol, byrow = byrow, ...)
}

#' @describeIn generic_multi_spct Specialization for collections of \code{response_spct} objects.
#'
#' @export
#'
response_multi_spct <- function(l, ncol = 1, byrow = FALSE, idfactor = TRUE, ...) {
  generic_multi_spct(l, class = "response_spct", ncol = ncol, byrow = byrow, ...)
}

#' @describeIn generic_multi_spct Specialization for collections of \code{chroma_spct} objects.
#'
#' @export
#'
chroma_multi_spct <- function(l, ncol = 1, byrow = FALSE, idfactor = TRUE, ...) {
  generic_multi_spct(l, class = "chroma_spct", ncol = ncol, byrow = byrow, ...)
}
