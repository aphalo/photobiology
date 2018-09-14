
# attributes used by spectral classes -------------------------------------

all_spct_attr.ls <-
  list(
    private = c("spct.version",
                "spct.tags"),
    generic_spct = c("comment",
                     "instr.desc",
                     "instr.settings",
                     "when.measured",
                     "where.measured",
                     "what.measured",
                     "spct.tags",
                     "normalized",
                     "scaled",
                     "multiple.wl",
                     "spct.version"),
    raw_spct = c("time.unit", "linearized"),
    cps_spct = c("time.unit", "linearized"),
    source_spct = c("time.unit", "bswf.used"),
    response_spct = c("time.unit", "bswf.used"),
    object_spct = c("Tfr.type", "Rfr.type", "Afr.type"),
    filter_spct = c("Tfr.type", "Rfr.type", "Afr.type"),
    reflector_spct = c("Tfr.type", "Rfr.type", "Afr.type"),
    calibration_spct = character(),
    chroma_spct = character())

private.attributes <- all_spct_attr.ls[["private"]]

common.attributes <- all_spct_attr.ls[["generic_spct"]]

all.attributes <- unique(unlist(all_spct_attr.ls, use.names = FALSE))

# copy_attributes ---------------------------------------------------------

#' Copy attributes
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
    which <- c(all_spct_attr.l[["private"]],
               all_spct_attr.ls[["generic_spct"]],
               all_spct_attr.ls[[class(y)[1]]])
  }
  attr.x <- attributes(x)
  which.x <- intersect(names(attr.x), which)
  # this is likely to be slow
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

# merge_attributes ---------------------------------------------------------

#' Merge and copy attributes
#'
#' Merge attributes from \code{x} and \code{y} and copy them to \code{z}.
#' Methods defined for spectral objects of classes from package 'photobiology'.
#'
#' @param x,y,z R objects
#' @param which character
#' @param ... not used
#'
#' @return A copy of \code{z} with additional attributes set.
#'
#' @export
#'
merge_attributes <- function(x, y, z, which, ...) UseMethod("merge_attributes")

#' @describeIn merge_attributes Default for generic function
#'
#' @export
#'
merge_attributes.default <- function(x, y, z,
                                    which = NULL,
                                    ...) {
  warning("'merge_attributes' is not defined for objects of class ", class(x)[1])
  z
}

#' @describeIn merge_attributes
#'
#' @param copy.class logical If TRUE class attributes are also copied.
#'
#' @export
#'
merge_attributes.generic_spct <- function(x, y, z,
                                         which = NULL,
                                         copy.class = FALSE,
                                         ...) {
  if (copy.class) {
    stopifnot(class(x) == class(y))
    class(z) <- class(x)
    check_spct(z)
  }
  stopifnot(is.generic_spct(y) && is.generic_spct(z))
  if (length(which) == 0) {
    which <- c(all_spct_attr.l[["private"]],
               all_spct_attr.ls[["generic_spct"]],
               all_spct_attr.ls[[class(y)[1]]])
  }
  # this is likely to be slow
  for (w in which) {
    att.x <- attr(x, w, exact = TRUE)
    att.y <- attr(y, w, exact = TRUE)
    if (length(att.x) == 0L && length(att.y) == 0L) {
      attr(z, w) <- NULL
    } else if (length(att.x) == 0L) {
      attr(z, w) <- att.y
    } else if (length(att.y) == 0L) {
      attr(z, w) <- att.x
    } else if (is.na(att.x) || is.na(att.y) ||
               class(att.x) != class(att.y) ||
               length(att.x) != length(att.y) ||
               xor(is.atomic(att.x), is.atomic(att.y)) ||
               (is.atomic(att.x) && any(att.x != att.y)) ||
               isFALSE(all.equal(att.x, att.y))) {
      attr(z, w) <- ifelse(w %in% c("comment", "time.unit"), NA_character_, NA)
    } else {
      attr(z, w) <- att.x
    }
  }
  setMultipleWl(z)
  z
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
                              allowed = c(all_spct_attr.ls["generic_spct"],
                                          all_spct_attr.ls["source_spct"]),
                              ...)
}

#' @describeIn get_attributes filter_spct
#' @export
#'
get_attributes.filter_spct <- function(x,
                                       which = NULL,
                                       ...) {
  get_attributes.generic_spct(x, which = which,
                              allowed = c(all_spct_attr.ls["generic_spct"],
                                          all_spct_attr.ls["filter_spct"]),
                              ...)
}


#' @describeIn get_attributes reflector_spct
#' @export
#'
get_attributes.reflector_spct <- function(x,
                                          which = NULL,
                                          ...) {
  get_attributes.generic_spct(x, which = which,
                              allowed = c(all_spct_attr.ls["generic_spct"],
                                          all_spct_attr.ls["reflector_spct"]),
                              ...)
}


#' @describeIn get_attributes object_spct
#' @export
#'
get_attributes.object_spct <- function(x,
                                       which = NULL,
                                       ...) {
  get_attributes.generic_spct(x, which = which,
                              allowed = c(all_spct_attr.ls["generic_spct"],
                                          all_spct_attr.ls["object_spct"]),
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

