# names of all multi spectral classes -------------------------------------------

#' Function that returns a vector containing the names of multi-spectra classes.
#'
#' @export
#'
#' @return A \code{character} vector of class names.
#'
mspct_classes <- function() {
  c("generic_mspct", "cps_mspct",
    "filter_mspct", "reflector_mspct",
    "source_mspct", "object_mspct",
    "response_mspct", "chroma_mspct")
}

# remove mspct class attributes --------------------------------------------

#' Remove "generic_mspct" and derived class attributes.
#'
#' Removes from an spectrum object the class attibutes "generic_mspct" and any
#' derived class attribute such as "source_mspct". \strong{This operation is done
#' by reference!}
#'
#' @param x an R object.
#' @export
#'
#' @note If \code{x} is an object of any of the multi spectral classes defined
#'   in this package, this function changes by reference the multi spectrum
#'   object into the underlying lis object. Otherwise, it just leaves \code{x}
#'   unchanged. The modified \code{x} is also returned invisibly.
#'
#' @return A character vector containing the removed class attribute values.
#'   This is different to the behaviour of function \code{unlist} in base R!
#'
#' @family set and unset 'multi spectral' class functions
#'
rmDerivedMspct <- function(x) {
  name <- substitute(x)
  mspctclasses <- mspct_classes()
  allclasses <- class(x)
  setattr(x, "class", setdiff(allclasses, mspctclasses))
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(setdiff(allclasses, class(x)))
}

# Constructors ------------------------------------------------------------

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

# is functions for mmspct classes --------------------------------------------

#' Query class of spectrum objects
#'
#' Functions to check if an object is of a given type of spectrum, or coerce it if
#' possible.
#'
#' @param x an R object.
#'
#' @return These functions return \code{TRUE} if its argument is a of the queried type
#'   of spectrum and \code{FALSE} otherwise.
#'
#' @note Derived types also return TRUE for a query for a base type such as
#' \code{generic_mspct}.
#'
#' @export
#' @rdname is.generic_mspct
#'
is.generic_mmspct <- function(x) inherits(x, "generic_mmspct")

#' @rdname is.generic_mspct
#' @export
#'
is.cps_mspct <- function(x) inherits(x, "cps_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.source_mspct <- function(x) inherits(x, "source_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.response_mspct <- function(x) inherits(x, "response_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.filter_mspct <- function(x) inherits(x, "filter_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.reflector_mspct <- function(x) inherits(x, "reflector_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.object_mspct <- function(x) inherits(x, "object_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.chroma_mspct <- function(x) inherits(x, "chroma_mspct")

#' @rdname is.generic_mspct
#'
#' @export
#'
is.any_mspct <- function(x) {
  inherits(x, mspct_classes())
}

# as functions for mspct classes --------------------------------------------

#' Return a copy of an R object as an spectrum object
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object
#'
#' @return These functions return a copy of \code{x} converted into a given
#'   class of spectral object, if \code{x} is a valid argument to the
#'   correcponding set function.
#'
#' @export
#'
#' @family creation of spectral objects functions
#' @rdname as.generic_mspct
#'
as.generic_mspct <- function(x) {
  y <- copy(x)
  rmDerivedMspct(y)
  z <- lapply(y, setGenericSpct)
  generic_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @export
#'
as.cps_mspct <- function(x) {
  y <- copy(x)
  rmDerivedMspct(y)
  z <- lapply(y, setCpsSpct)
  cps_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @param time.unit character A string, "second", "day" or "exposure"
#' @param bswf.used character
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning
#'
#' @export
#'
as.source_mspct <- function(x,
                           time.unit=c("second", "day", "exposure"),
                           bswf.used=c("none", "unknown"),
                           strict.range = FALSE) {
  y <- copy(x)
  rmDerivedMspct(y)
  z <- lapply(y, setSourceSpct, time.unit = time.unit, strict.range = strict.range, bswf.used = bswf.used)
  source_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @export
#'
as.response_mspct <- function(x, time.unit = "second") {
  y <- copy(x)
  rmDerivedMspct(y)
  z <- lapply(y, setResponseSpct, time.unit = time.unit)
  reponse_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @param Tfr.type a character string, either "total" or "internal"
#'
#' @export
#'
as.filter_mspct <- function(x, Tfr.type=c("total", "internal"), strict.range = TRUE) {
  y <- copy(x)
  rmDerivedMspct(y)
  z <- lapply(y, setFilterSpct, Tfr.type = Tfr.type, strict.range = strict.range)
  filter_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @param Rfr.type a character string, either "total" or "specular"
#'
#' @export
#'
as.reflector_mspct <- function(x, Rfr.type = c("total", "specular"), strict.range = TRUE) {
  y <- copy(x)
  rmDerivedMspct(y)
  z <- lapply(y, setReflectorSpct, Rfr.type = Rfr.type, strict.range = strict.range)
  reflector_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @export
#'
as.object_mspct <- function(x,
                           Tfr.type=c("total", "internal"),
                           Rfr.type=c("total", "specular"),
                           strict.range = TRUE) {
  y <- copy(x)
  rmDerivedMspct(y)
  z <- lapply(y, setObjectSpct, Tfr.type = Tfr.type, Rfr.type = Rfr.type,
              strict.range = strict.range)
  object_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @export
#'
as.chroma_mspct <- function(x) {
  y <- copy(x)
  rmDerivedMspct(y)
  z <- lapply(y, setChromaSpct)
  chroma_mspct(z)
}


