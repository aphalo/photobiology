
# attributes used by spectral classes -------------------------------------

private.attributes <- c("spct.version",
                        "spct.tags")

common.attributes <- c("comment",
                       "what.measured",
                       "when.measured",
                       "where.measured",
                       "normalized",
                       "scaled",
                       "multiple.wl")

source.attributes <- c("bswf.used", "time.unit")
filter.attributes <- c("Tfr.type")
reflector.attributes <- c("Rfr.type")
object.attributes <- union(filter.attributes, reflector.attributes)
instr.attributes <- c("instr.descriptor", "instr.settings")
calibration.attributes <- NULL
all.attributes <- unique(c(private.attributes,
                           common.attributes,
                           source.attributes,
                           filter.attributes,
                           reflector.attributes,
                           instr.attributes))

# copy_attributes ---------------------------------------------------------

#' Copy attributes from one R object to another
#'
#' Copy attributes from \code{x} to \code{y}. Methods defined for spectral
#' and waveband objects of classes from package 'photobiology'.
#'
#' @param x,y R objects
#' @param which character
#' @param ... not used
#'
#' @return A copy of \code{y} with additional attributes set.
#'
#' @export
#'
copy_attributes <- function(x, y, which, ...) UseMethod("copy_attributes")

#' @describeIn copy_attributes Default for generic function
#'
#' @export
#'
copy_attributes.default <- function(x, y,
                                    which = NULL,
                                    ...) {
  warning("'copy_attributes' is not defined for objects of class ", class(x)[1])
  y
}

#' @describeIn copy_attributes
#'
#' @param copy.class logical If TRUE class attributes are also copied.
#'
#' @export
#'
copy_attributes.generic_spct <- function(x, y,
                                         which = NULL,
                                         copy.class = FALSE,
                                         ...) {
  if (copy.class) {
    class(y) <- class(x)
    check_spct(y)
  }
  stopifnot(is.generic_spct(y))
  if (length(which) == 0) {
    which <- c("comment",
               "instr.desc",
               "instr.settings",
               "when.measured",
               "where.measured",
               "what.measured",
               "spct.tags",
               "normalized",
               "scaled",
               "multiple.wl",
               "spct.version")
    which.add <- c(switch, class(y),
                   generic_spct,
                   raw_spct = "linearized",
                   cps_spct = "linearized",
                   source_spct = c("time.unit", "bswf.used"),
                   response_spct = c("time.unit", "bswf.used"),
                   # need to be copied in case class of object_spct
                   # is changed temporarily
                   filter_spct = c("Tfr.type", "Rfr.type", "Afr.type"),
                   reflector_spct = c("Tfr.type", "Rfr.type"),
                   object_spct = c("Tfr.type", "Rfr.type", "Afr.type"),
                   chroma_spct = character()
    )
    which <- c(which, which.add)
  }
  attr.x <- attributes(x)
  which.x <- intersect(names(attr.x), which)
  # this is likely to be very slow
  for (w in which.x) {
        attr(y, w) <- attr(x, w, exact = TRUE)
  }
  y
}

#' @describeIn copy_attributes
#'
#' @export
#'
copy_attributes.waveband <- function(x, y, which = NULL, ...) {
  stopifnot(is.waveband(y))
  if (length(which) == 0L) {
    which <- "comment"
  }
  attr.x <- attributes(x)
  which.x <- intersect(names(attr.x), which)
  for (w in which.x) {
    attr(y, w) <- attr(x, w, exact = TRUE)
  }
  y
}

# get_attributes -----------------------------------------------------------

#' Get the metadata attributes
#'
#' Method returning attributes of an object of class generic_spct or derived,
#' or of class waveband. Only attributes defined and/or set by package
#' 'photobiology' for objects of the corresponding class are returned.
#'
#' @param x a generic_spct object.
#' @param which character vector Names of attributes to retrieve.
#' @param allowed character vector Names of attributes accepted by \code{which}.
#' @param ... currently ignored
#'
#' @return Named \code{list} of attribute values.
#'
#' @export
#'
#' @family measurement metadata functions
#'
get_attributes <-
  function(x, which, ...) UseMethod("get_attributes")

#' @describeIn get_attributes generic_spct
#' @export
#'
get_attributes.generic_spct <-
  function(x,
           which = NULL,
           allowed = all.attributes,
           ...) {
    if (length(which) == 0L) {
      which <- allowed
    }
    spct.attr <- attributes(x)
    spct.attr[names(spct.attr) %in% intersect(which, allowed)]
  }

#' @describeIn get_attributes source_spct
#' @export
#'
get_attributes.source_spct <- function(x,
                                       which = NULL,
                                       ...) {
  get_attributes.generic_spct(x, which = which,
                              allowed = c(source.attributes,
                                          common.attributes),
                              ...)
}

#' @describeIn get_attributes filter_spct
#' @export
#'
get_attributes.filter_spct <- function(x,
                                       which = NULL,
                                       ...) {
  get_attributes.generic_spct(x, which = which,
                              allowed = c(filter.attributes,
                                          common.attributes),
                              ...)
}


#' @describeIn get_attributes reflector_spct
#' @export
#'
get_attributes.reflector_spct <- function(x,
                                          which = NULL,
                                          ...) {
  get_attributes.generic_spct(x, which = which,
                              allowed = c(reflector.attributes,
                                          common.attributes),
                              ...)
}


#' @describeIn get_attributes object_spct
#' @export
#'
get_attributes.object_spct <- function(x,
                                       which = NULL,
                                       ...) {
  get_attributes.generic_spct(x, which = which,
                              allowed = c(object.attributes,
                                          common.attributes),
                              ...)
}

#' @describeIn get_attributes waveband
#' @export
#'
get_attributes.waveband <- function(x,
                                    which = NULL,
                                    ...) {
  if (length(which) == 0L || which == "comment")
  list(comment = attr(x, "comment", exact = TRUE))
}

