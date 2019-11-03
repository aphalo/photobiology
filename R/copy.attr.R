
# attributes used by spectral classes -------------------------------------

all_spct_attr.ls <-
  list(
    private = c("spct.version",
                "spct.tags",
                "na.action"),
    generic_spct = c("comment",
                     "instr.desc",
                     "instr.settings",
                     "when.measured",
                     "where.measured",
                     "what.measured",
                     "how.measured",
#                     "spct.tags",
#                     "spct.version",
                     "normalized",
                     "scaled",
                     "multiple.wl"),
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

# select_spct_attributes --------------------------------------------------

#' Merge user supplied attribute names with default ones
#'
#' Allow users to add and subract from default attributes in addition to
#' providing a given set of attributes.
#'
#' @param attributes,attributes.default character vector or a list of character
#'   vectors.
#'
#' @details Vectors of character strings passed as argument to \code{attributes}
#'   are parsed so that if the first member string is \code{"+"}, the remaining
#'   members are added to those in \code{attributes.default}; if it is
#'   \code{"-"} the remaining members are removed from in
#'   \code{attributes.default}; and if it is \code{"="} the remaining members
#'   replace those in in \code{attributes.default}. If the first member is none
#'   of these three strings, the behaviour is the same as when the first string
#'   is \code{"="}. If \code{attributes} is \code{NULL} all the attributes in
#'   \code{attributes.default} are used and if it is \code{""} no attribute
#'   names are returned, \code{""} has precedence over other member values. The
#'   order of the names of annotations has no meaning: the vector is interpreted
#'   as a set except for the three possible "operators" at position 1.
#'
#' @return A character vector of attribute names.
#'
#' @seealso \code{\link{get_attributes}}
#'
#' @export
#'
#' @family measurement metadata functions
#'
select_spct_attributes <- function(attributes,
                                   attributes.default = spct_attributes()) {
  if (length(attributes) == 0L) { # handle character(0) and NULL without delay
    return(attributes.default)
  } else if (is.list(attributes)) {
    attributes.ls <- attributes
  } else if (is.character(attributes)) {
    attributes.ls <- list(attributes)
  }
  attributes <- NULL

  for (attributes in attributes.ls) {
    stopifnot(is.character(attributes))
    # we can receive character(0) from preceeding iteration
    if (length(attributes) == 0L || attributes[1] == "*") {
      z <- attributes.default
    } else if ("" %in% attributes) {
      # no annotations
      z <- ""
    } else if (attributes[1] == "-") {
      # remove exact matches
      z <- setdiff(attributes.default, attributes[-1])
    } else if (attributes[1] == "+") {
      attributes <- attributes[-1]
      # merge default with addition
      z <- union(attributes.default, attributes)
    } else if (attributes[1] == "=") {
      # replace
      z <- attributes[-1]
      # handle character(0), using "" is a kludge but produces intuitive behaviour
      if (length(z) == 0L) {
        z <- ""
      }
    } else {
      z <- attributes
    }
    attributes.default <- z
  }

  unique(z) # remove duplicates
}

#' @rdname select_spct_attributes
#'
#' @param .class character Name of spectral class.
#' @export
#'
spct_attributes <- function(.class = "all", attributes = "*") {
  if ("all" %in% .class) {
    attributes.default <- unlist(all_spct_attr.ls, use.names = FALSE)
  } else {
    attributes.default <- unlist(all_spct_attr.ls[union("generic_spct", .class)], use.names = FALSE)
  }
  select_spct_attributes(attributes = attributes,
                         attributes.default = attributes.default)
}

# copy_attributes ---------------------------------------------------------

#' Copy attributes
#'
#' Copy attributes from \code{x} to \code{y}. Methods defined for spectral
#' and waveband objects of classes from package 'photobiology'.
#'
#' @param x,y R objects
#' @param which character Names of attributes to copy, if NULL all those
#'    relevant according to the class of \code{x} is used as defaul,
#' @param which.not character Names of attributes not to be copied. The
#'    names passed here are removed from the list for \code{which}, which
#'    is most useful when we want to modify the default.
#' @param copy.class logical If TRUE class attributes are also copied.
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
#' @export
#'
copy_attributes.generic_spct <- function(x, y,
                                         which = NULL,
                                         which.not = NULL,
                                         copy.class = FALSE,
                                         ...) {
  if (copy.class) {
    class(y) <- class(x)
    check_spct(y)
  }
  stopifnot(is.generic_spct(y))
  if (length(which) == 0) {
    which <- c(all_spct_attr.ls[["private"]],
               all_spct_attr.ls[["generic_spct"]],
               all_spct_attr.ls[[class(y)[1]]])
  }
  which <- setdiff(which, which.not)
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
copy_attributes.generic_mspct <- function(x, y,
                                         which = NULL,
                                         which.not = NULL,
                                         copy.class = FALSE,
                                         ...) {
  stopifnot(length(x) == length(y))
  for (i in seq_along(x)) {
    y[[i]] <- copy_attributes(x[[i]], y[[i]],
                              which = which,
                              which.not = which.not,
                              copy.class = copy.class,
                              ...)
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
    stopifnot(class_spct(x)[1] == class_spct(y)[1])
    class(z) <- class(x)
    check_spct(z)
  }
  stopifnot(is.generic_spct(y) && is.generic_spct(z))
  if (length(which) == 0) {
    which <- c(all_spct_attr.ls[["private"]],
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
    } else if (any(is.na(att.x)) || any(is.na(att.y)) ||
               class(att.x)[1] != class(att.y)[1] ||
               length(att.x) != length(att.y) ||
               xor(is.atomic(att.x), is.atomic(att.y))) {
      attr(z, w) <- ifelse(w %in% c("comment", "time.unit"), NA_character_, NA)
    } else {
      ## Add generic test of equality to warning
      if (getOption("photobiology.verbose", default = FALSE)) {
        warning("Keeping attribute ", w, "'s value from lhs operand.")
      }
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
#' 'photobiology' for objects of the corresponding class are returned. Parameter
#' \code{which} can be used to subset the list of attributes.
#'
#' @param x a generic_spct object.
#' @param which character vector Names of attributes to retrieve.
#' @param allowed character vector Names of attributes accepted by \code{which}.
#' @param ... currently ignored
#'
#' @return Named \code{list} of attribute values.
#'
#' @details Vectors of character strings passed as argument to \code{which} are
#'   parsed so that if the first member string is \code{"-"} the remaining
#'   members are removed from the \code{allowed}; and if it is \code{"="} the
#'   remaining members are used if in \code{allowed}. If the first member is
#'   none of these three strings, the behaviour is the same as if the first
#'   string is \code{"="}. If \code{which} is \code{NULL} all the attributes in
#'   \code{allowed} are used. The string \code{""} means no attributes, and has
#'   precedence over any other values in the character vector. The order of the
#'   names of annotations has no meaning: the vector is interpreted as a set
#'   except for the three possible "operators" at position 1.
#'
#' @seealso \code{\link{select_spct_attributes}}
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
      which <- "*"
    }
    target.attributes <- select_spct_attributes(attributes = which,
                                                attributes.default = allowed)
    spct.attr <- attributes(x)
    spct.attr[names(spct.attr) %in% target.attributes]
  }

#' @describeIn get_attributes source_spct
#' @export
#'
get_attributes.source_spct <- function(x,
                                       which = NULL,
                                       ...) {
  get_attributes.generic_spct(x, which = which,
                              allowed = c(all_spct_attr.ls[["generic_spct"]],
                                          all_spct_attr.ls[["source_spct"]]),
                              ...)
}

#' @describeIn get_attributes filter_spct
#' @export
#'
get_attributes.filter_spct <- function(x,
                                       which = NULL,
                                       ...) {
  get_attributes.generic_spct(x, which = which,
                              allowed = c(all_spct_attr.ls[["generic_spct"]],
                                          all_spct_attr.ls[["filter_spct"]]),
                              ...)
}


#' @describeIn get_attributes reflector_spct
#' @export
#'
get_attributes.reflector_spct <- function(x,
                                          which = NULL,
                                          ...) {
  get_attributes.generic_spct(x, which = which,
                              allowed = c(all_spct_attr.ls[["generic_spct"]],
                                          all_spct_attr.ls[["reflector_spct"]]),
                              ...)
}


#' @describeIn get_attributes object_spct
#' @export
#'
get_attributes.object_spct <- function(x,
                                       which = NULL,
                                       ...) {
  get_attributes.generic_spct(x, which = which,
                              allowed = c(all_spct_attr.ls[["generic_spct"]],
                                          all_spct_attr.ls[["object_spct"]]),
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

# attributes2tb -----------------------------------------------------------

#' Copy attributes into a tibble
#'
#' Method returning attributes of an object of class generic_spct or derived,
#' or of class waveband. Only attributes defined and/or set by package
#' 'photobiology' for objects of the corresponding class are returned.
#'
#' @param x a generic_spct object.
#' @param which character vector Names of attributes to retrieve.
#' @param ... currently ignored
#'
#' @return A tibble with the values stored in the attributes whose names were
#'   selected through the argument to \code{which} if present in \code{x}.
#'
#' @export
#'
#' @family measurement metadata functions
#'
spct_attr2tb <-
  function(x,
           which = c("-", "names", "row.names", "spct.tags", "spct.version", "comment"),
           ...) {
    spct.attr <- get_attributes(x = x, which = which, ...)
    as_tibble(spct.attr)
  }

