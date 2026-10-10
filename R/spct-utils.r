# spct and mspct utility methods ----------------------------------------------

#' Extract all members from a collection
#'
#' Extract all members from a collection into separate objects in the parent
#' frame of the call.
#'
#' @param x An R object
#' @param ... additional named arguments passed down to \code{f}.
#'
#' @return Utility used for its side effects, invisibly returns a character
#'   vector with the names of the objects created.
#'
#' @export
#'
#' @examples
#'
#' my.mscpt <- source_mspct(list(sun1.spct = sun.spct, sun2.spct = sun.spct))
#' uncollect2spct(my.mscpt)
#' ls(pattern = "*.spct")
#'
#' @family experimental utility functions
#'
uncollect2spct <- function(x, ...) UseMethod("uncollect2spct")

#' @rdname uncollect2spct
#'
#' @export
#'
uncollect2spct.default <- function(x, ...) {
  warning("'uncollect2spct()' is not defined for objects of class '", class(x)[1], "'.")
  invisible(character())
}

#' @rdname uncollect2spct
#'
#' @param name.tag character. A string used as tag for the names of the objects.
#'   If of length zero, names of members are used as named of objects. Otherwise
#'   the tag is appended, unless already present in the member name.
#' @param ignore.case	logical. If FALSE, the pattern matching used for \code{name.tag} is
#'   case sensitive and if TRUE, case is ignored during matching.
#' @param check.names logical. If TRUE then the names of the objects created are
#'   checked to ensure that they are syntactically valid variable names and
#'   unique. If necessary they are adjusted (by make.names) so that they are,
#'   and if FALSE names are used as is.
#' @param check.overwrite logical. If TRUE trigger an error if an exisitng
#'   object would be overwritten, and if FALSE silently overwrite objects.
#'
#' @export
#'
uncollect2spct.generic_mspct <- function(x,
                                    name.tag = ".spct",
                                    ignore.case = FALSE,
                                    check.names = TRUE,
                                    check.overwrite = TRUE,
                                    ...) {
  stopifnot(rlang::is_named(x))
  # make object names
  if (length(name.tag) == 0L) {
    obj.names <- names(x)
  } else {
    obj.names <- character()
    for (member.name in names(x)) {
      if (!grepl(pattern = paste("*", name.tag, "$", sep = ""),
                 x = member.name,
                 ignore.case = ignore.case)) {
        obj.name <- paste(member.name, name.tag, sep = "")
      } else {
        obj.name <- member.name
      }
      obj.names <- c(obj.names, obj.name)
    }
  }
  # by checking a vector of object names we protect from repeated names
  if (check.names) {
    obj.names <- make.names(obj.names, unique = TRUE)
  }
  names(obj.names) <- names(x)

  # check for exisiting objects
  if (check.overwrite) {
    existing.names <-
      ls(envir = parent.frame(2), pattern = paste("*", name.tag, "$", sep = ""), sorted = FALSE)
    if (length(intersect(obj.names, existing.names)) > 0L) {
      stop("Overwriting would take place, if intended, change default.")
    }
  }

  # create objects
  for (member.name in names(x)) {
    assign(obj.names[member.name], x[[member.name]], envir = parent.frame(2))
  }
  invisible(unname(obj.names))
}

#' Form a new collection
#'
#' Form a collection of spectra from separate objects in the parent
#' frame of the call.
#'
#' @details This is a convenience function that simplifies the creation of
#'   collections from existing objects of class \code{generic_spct} or a derived
#'   class. A list of objects con be passed as argument, or a search pattern. If
#'   a list is passed, no search is done. If \code{collection.class} is
#'   \code{NULL}, then all objects of class \code{generic_spct} or of a class
#'   derived from it are added to the collection. If objects of only one derived
#'   class are to be collected this class or that of the matching collection
#'   should be passed as argument to \code{collection.class}. Objects of other R
#'   classes are silently discarded, which simplifies the specification of
#'   search patterns. By default, i.e., if \code{collection.class} is
#'   \code{NULL}, if all the objects collected belong to the same class then the
#'   corresponding collection class will be returned, otherwise a
#'   \code{generic_mspct} object with heterogeneous members will be returned. To
#'   force the return of a \code{generic_mspct} even when the collected spectra
#'   all belong to the same class, pass \code{generic_mspct} as argument to
#'   \code{collection.class}. If the argument to \code{collection.class} is a
#'   vector containing two of more class names, only the matching spectra will
#'   be collected, and a \code{generic_mspct} will be returned. The returned
#'   object is created with the constructor for the class, and validated.
#'
#' @param .list list of R objects
#' @param pattern character an optional regular expression, ignored if
#'   \code{.list} is not \code{NULL}.
#' @param collection.class character vector
#' @param ... additional named arguments passed down to the collection
#'   constructor.
#'
#' @return By default a collection of spectra.
#'
#' @export
#'
#' @examples
#' collect2mspct() # returns empty generic_mspct object
#'
#' sun1.spct <- sun.spct
#' sun2.spct <- sun.spct
#' kk.spct <- 10:30 # ignored
#' collect2mspct()
#' collect2mspct(collection.class = "generic_mspct")
#'
#' pet1.spct <- polyester.spct
#' collect2mspct()
#' collect2mspct(collection.class = "source_mspct")
#' collect2mspct(collection.class = "filter_mspct")
#' collect2mspct(collection.class = "response_mspct")
#'
#' @family experimental utility functions
#'
collect2mspct <- function(.list = NULL,
                          pattern = "*\\.spct$",
                          collection.class = NULL,
                          ...)  {
  collection.class <- gsub("_spct$", "_mspct", collection.class) # NULL -> character(0)!!
  if (length(collection.class) > 0L &&
      !all(collection.class %in% mspct_classes())) {
    warning("Discarding unrecognized class(es) in 'collection.class'")
    collection.class <- intersect(collection.class, spct_classes())
  }
  if (length(.list) == 0L) {
    names <- ls(pattern = pattern, envir = parent.frame())
    .list <- mget(x = names, envir = parent.frame())
  }
  if (length(.list) >= 1L) {
    if (length(collection.class) == 0L) {
      # we keep objects of class generic_spct or derived
      .list <- .list[sapply(.list, is.any_spct)]
      members.class <- shared_member_class(.list)[1] # derived class may be used
      collection.class <- gsub("_spct$", "_mspct", members.class)
    } else {
      # we keep objects of the classes in collection.class
      members.class <- gsub("_mspct$", "_spct", collection.class)
      if ("generic_spct" %in% members.class) {
        # we keep objects of class generic_spct or derived
        .list <-
          .list[sapply(X = .list, FUN = is.any_spct)]
        collection.class <- "generic_mspct" # forced
      } else {
        # we keep objects of classes listed in collection.class
        .list <-
          .list[sapply(.list, function(x) {class(x)[1]}) %in% members.class]
        if (length(collection.class) > 1L) {
          collection.class <- "generic_mspct"
        }
      }
    }
    # call the constructor for the class
    do.call(what = collection.class, args = list(l = .list))
  } else if (length(collection.class) == 1L) {
    do.call(what = collection.class, args = list())
  } else {
    generic_mspct()
  }
}

# drop user columns -------------------------------------------------------

#' Drop user columns
#'
#' Remove from spectral object additional columns that are user defined.
#'
#' @param x An R object
#' @param keep.also character Additionlal columns to preserve.
#' @param ... needed to allow derivation.
#'
#' @return A copy of \code{x} possibly with some columns removed.
#'
#' @export
#'
#' @family experimental utility functions
#'
drop_user_cols <- function(x, keep.also, ...) UseMethod("drop_user_cols")

#' @rdname drop_user_cols
#'
#' @export
#'
drop_user_cols.default <- function(x, keep.also = NULL, ...) {
  warning("'drop_user_cols()' is not defined for objects of class '", class(x)[1], "'.")
  x
}

#' @rdname drop_user_cols
#'
#' @export
#'
drop_user_cols.generic_spct <- function(x, keep.also, ...) {
  stopifnot(length(keep.also) >= 1L)
  default.cols <- c("w.length")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @rdname drop_user_cols
#'
#' @export
#'
drop_user_cols.source_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "s.e.irrad", "s.q.irrad")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @rdname drop_user_cols
#'
#' @export
#'
drop_user_cols.response_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "s.e.response", "s.q.response")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @rdname drop_user_cols
#'
#' @export
#'
drop_user_cols.object_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "Tfr", "Rfr", "Afr")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @rdname drop_user_cols
#'
#' @export
#'
drop_user_cols.filter_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "Tfr", "Afr", "A")
  if ("Rfr" %in% colnames(x)) {
    warning("Deleting 'object_spct' columns from 'filter_spct'.")
  }
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @rdname drop_user_cols
#'
#' @export
#'
drop_user_cols.reflector_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "Rfr")
  if (any(c("Tfr", "Afr", "A") %in% colnames(x))) {
    warning("Deleting 'object_spct' columns from 'reflector_spct'.")
  }
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @rdname drop_user_cols
#'
#' @export
#'
drop_user_cols.solute_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "K.mole", "K.mass")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @rdname drop_user_cols
#'
#' @export
#'
drop_user_cols.chroma_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "x", "y", "z")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @rdname drop_user_cols
#'
#' @export
#'
drop_user_cols.calibration_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "irrad.mult")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @rdname drop_user_cols
#'
#' @export
#'
drop_user_cols.cps_spct <- function(x, keep.also = NULL, ...) {
  cps.cols <- grep("^cps", colnames(x), value = TRUE)
  default.cols <- c("w.length", cps.cols)
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @rdname drop_user_cols
#'
#' @export
#'
drop_user_cols.raw_spct <- function(x, keep.also = NULL, ...) {
  counts.cols <- grep("^counts", colnames(x), value = TRUE)
  default.cols <- c("w.length", counts.cols)
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @rdname drop_user_cols
#'
#' @export
#'
drop_user_cols.generic_mspct <-
  function(x, keep.also = NULL, ...) {
    msmsply(x, drop_user_cols, keep.also = keep.also, ...)
}

