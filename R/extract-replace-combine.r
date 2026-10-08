
# Subset ------------------------------------------------------------------

# subset.data.frame() should work as expected with all spectral classes as it
# calls the Extract methods defined below on the object passed to x!
#
# However the methods defined bellow fail to retain attributes when j = TRUE,
# which is what subset() passes.
# The extract methods behave as R data.frame does, this may need to be changed
# but meanwhile we include here our own definition of subset to retain the
# expected behaviour for subset().

#' Subsetting spectra
#'
#' Return subsets of spectra stored in class \code{generic_spct} or derived from
#' it.
#'
#' @param x object to be subsetted.
#' @param subset logical expression indicating elements or rows to keep: missing
#'   values are taken as false.
#' @param drop passed on to \code{[} indexing operator.
#' @param select expression, indicating columns to select from a spectrum.
#' @param ...	further arguments to be passed to or from other methods.
#'
#' @return An object similar to \code{x} containing just the selected rows and
#'   columns. Depending on the columns remaining after subsetting the class of
#'   the object will be simplified to the most derived parent class.
#'
#' @export
#'
#' @method subset generic_spct
#'
#' @name Subset
#' @rdname subset
#'
#' @note This method is copied from \code{base::subset.data.frame()} but ensures
#'   that all metadata stored in attributes of spectral objects are copied to
#'   the returned value.
#'
#' @examples
#'
#' subset(sun.spct, w.length > 400)
#'
subset.generic_spct <- function(x, subset, select, drop = FALSE, ...) {
  r <- if (missing(subset))
    rep_len(TRUE, nrow(x))
  else {
    e <- substitute(subset)
    r <- eval(e, x, parent.frame())
    if (!is.logical(r))
      stop("'subset' must be logical")
    r & !is.na(r)
  }
  vars <- if (missing(select))
    rep_len(TRUE, ncol(x))
  else {
    nl <- as.list(seq_along(x))
    names(nl) <- names(x)
    eval(substitute(select), nl, parent.frame())
  }
  z <- x[r, vars, drop = drop]
  z <- copy_attributes(x, z)
  id.factor <- getIdFactor(x)
  if (!is.na(id.factor)) {
    # drop unused levels
    z[[id.factor]] <- factor(z[[id.factor]])
    multiple_wl(z) <- length(levels(z[[id.factor]]))
    # keep attributes matching remaining spectra
    z <- subset_attributes(z, to.keep = levels(z[[id.factor]]))
  }
  z
}

# Extract ------------------------------------------------------------------

# $ operator for extraction does not need any wrapping as it always extracts
# single columns returning objects of the underlying classes (e.g. numeric)
# rather than spectral objects.
#
# [ needs special handling as it can be used to extract rows, or groups of
# columns which are returned as spectral objects. Such returned objects
# can easily become invalid, for example, lack a w.length variable.

#' Extract or replace parts of a spectrum
#'
#' Just like extraction and replacement with indexes in base R, but preserving
#' the special attributes used in spectral classes and checking for validity of
#' remaining spectral data.
#'
#' @param x	spectral object from which to extract element(s) or in which to replace element(s)
#' @param i index for rows,
#' @param j index for columns, specifying elements to extract or replace. Indices are
#'   numeric or character vectors or empty (missing) or NULL. Please, see
#'   \code{\link[base]{Extract}} for more details.
#' @param drop logical. If TRUE the result is coerced to the lowest possible
#'   dimension. The default is FALSE unless the result is a single column.
#'
#' @details These methods are just wrappers on the method for data.frame objects
#'   which copy the additional attributes used by these classes, and validate
#'   the extracted object as a spectral object. When drop is TRUE and the
#'   returned object has only one column, then a vector is returned. If the
#'   extracted columns are more than one but do not include \code{w.length}, a
#'   data frame is returned instead of a spectral object.
#'
#' @return An object of the same class as \code{x} but containing only the
#'   subset of rows and columns that are selected. See details for special
#'   cases.
#'
#' @note If any argument is passed to \code{j}, even \code{TRUE}, some metadata
#'   attributes are removed from the returned object. This is how the
#'   extraction operator works with \code{data.frames} in R. For the time
#'   being we retain this behaviour for spectra, but it may change in the
#'   future.
#'
#' @method [ generic_spct
#'
#' @examples
#' sun.spct[sun.spct[["w.length"]] > 400, ]
#' subset(sun.spct, w.length > 400)
#'
#' tmp.spct <- sun.spct
#' tmp.spct[tmp.spct[["s.e.irrad"]] < 1e-5 , "s.e.irrad"] <- 0
#' e2q(tmp.spct[ , c("w.length", "s.e.irrad")]) # restore data consistency!
#'
#' @rdname extract
#' @name Extract
#'
#' @seealso \code{\link[base]{subset}} and \code{\link{trim_spct}}
#'
"[.generic_spct" <-
  function(x, i, j, drop = NULL) {
    if (is.null(drop)) {
      xx <- `[.data.frame`(x, i, j)
    } else {
      xx <- `[.data.frame`(x, i, j, drop = drop)
    }
    if (is.data.frame(xx)) {
      if ("w.length" %in% names(xx)) {
        # still a generic_spct object
        xx <- copy_attributes(x, xx)
        if (!(getMultipleWl(x) == 1L || nrow(xx) == nrow(x))) {
          # subsetting of rows can decrease the number of spectra
          id.factor <- getIdFactor(x)
          if (!is.na(id.factor)) {
            # drop unused levels only if needed for performance
            if (length(unique(xx[[id.factor]])) != length(levels(x[[id.factor]]))) {
              xx[[id.factor]] <- factor(xx[[id.factor]])
              multiple_wl(xx) <- length(levels(xx[[id.factor]]))
              # keep attributes matching remaining spectra
              xx <- subset_attributes(xx, to.keep = levels(xx[[id.factor]]))
            }
            # disable check of wavelengths as known good
            xx <- check_spct(xx, strict.range = NULL, multiple.wl = NULL)
          } else {
            xx <- check_spct(xx, strict.range = NULL)
          }
        } else {
          # disable check of wavelengths as known good
          xx <- check_spct(xx, strict.range = NULL, multiple.wl = NULL)
        }
      } else {
        # no longer a valid spectrum or spectra object
        rmDerivedSpct(xx)
      }
    }
    xx
  }

#' @export
#' @rdname extract
#'
"[.raw_spct" <-
  function(x, i, j, drop = NULL) {
    if (is.null(drop)) {
      xx <- `[.data.frame`(x, i, j)
    } else {
      xx <- `[.data.frame`(x, i, j, drop = drop)
    }
    if (is.data.frame(xx)) {
      if ("w.length" %in% names(xx)) {
        if (!(getMultipleWl(x) == 1L || nrow(xx) == nrow(x))) {
          # subsetting of rows can decrease the number of spectra
          multiple.wl <- findMultipleWl(xx, same.wls = FALSE)
          xx <- setMultipleWl(xx, multiple.wl)
        }
        if (ncol(xx) != ncol(x)) {
          xx <- copy_attributes(x, xx)
        }
        xx <- check_spct(xx)
      } else {
        rmDerivedSpct(xx)
      }
    }
    xx
  }

#' @export
#' @rdname extract
#'
"[.cps_spct" <-
  function(x, i, j, drop = NULL) {
    if (is.null(drop)) {
      xx <- `[.data.frame`(x, i, j)
    } else {
      xx <- `[.data.frame`(x, i, j, drop = drop)
    }
    if (is.data.frame(xx)) {
      if ("w.length" %in% names(xx)) {
        if (!(getMultipleWl(x) == 1L || nrow(xx) == nrow(x))) {
          # subsetting of rows can decrease the number of spectra
          multiple.wl <- findMultipleWl(xx, same.wls = FALSE)
          xx <- setMultipleWl(xx, multiple.wl)
        }
        if (ncol(xx) != ncol(x)) {
          xx <- copy_attributes(x, xx)
        }
        xx <- check_spct(xx)
      } else {
        rmDerivedSpct(xx)
      }
    }
    xx
  }

#' @export
#' @rdname extract
#'
"[.source_spct" <- `[.generic_spct`

#' @export
#' @rdname extract
#'
"[.response_spct" <-`[.generic_spct`

#' @export
#' @rdname extract
#'
"[.filter_spct" <-`[.generic_spct`

#' @export
#' @rdname extract
#'
"[.reflector_spct" <- `[.generic_spct`

#' @export
#' @rdname extract
#'
"[.solute_spct" <- `[.generic_spct`

#' @export
#' @rdname extract
#'
"[.object_spct" <- `[.generic_spct`

#' @export
#' @rdname extract
#'
"[.chroma_spct" <- `[.generic_spct`

# replace -----------------------------------------------------------------

# We need to wrap the replace functions adding a call to our check method
# to make sure that the object is still a valid spectrum after the
# replacement.

#' @param value	A suitable replacement value: it will be repeated a whole number
#'   of times if necessary and it may be coerced: see the Coercion section. If
#'   NULL, deletes the column if a single column is selected.
#'
#' @export
#' @method [<- generic_spct
#' @rdname extract
#'
"[<-.generic_spct" <- function(x, i, j, value) {
  check_spct(`[<-.data.frame`(x, i, j, value), byref = FALSE)
}

#' @param name A literal character string or a name (possibly backtick quoted).
#'   For extraction, this is normally (see under 'Environments') partially
#'   matched to the names of the object.
#'
#' @export
#' @method $<- generic_spct
#' @rdname extract
#'
"$<-.generic_spct" <- function(x, name, value) {
  check_spct(`$<-.data.frame`(x, name, value), byref = FALSE)
}

# Extract ------------------------------------------------------------------

# $ operator for extraction does not need any wrapping as it always extracts
# single objects of the underlying classes (e.g. generic_spct)
# rather than collections of spectral objects.
#
# [ needs special handling as it can be used to extract members, or groups of
# members which must be returned as collections of spectral objects.
#
# In the case of replacement, collections of objects can easily become invalid,
# if the replacement or added member belongs to a class other than the expected
# one(s) for the collection.

#' Extract or replace members of a collection of spectra
#'
#' Just like extraction and replacement with indexes for base R lists, but
#' preserving the special attributes used in spectral classes.
#'
#' @param x	Collection of spectra object from which to extract member(s) or in
#'   which to replace member(s)
#' @param i Index specifying elements to extract or replace. Indices are numeric
#'   or character vectors. Please, see \code{\link[base]{Extract}} for
#'   more details.
#' @param drop If TRUE the result is coerced to the lowest possible dimension
#'   (see the examples). This only works for extracting elements, not for the
#'   replacement.
#'
#' @details This method is a wrapper on base R's extract method for lists that
#'   sets additional attributes used by these classes.
#'
#' @return An object of the same class as \code{x} but containing only the
#'   subset of members that are selected.
#'
#' @method [ generic_mspct
#' @export
#'
#' @rdname extract_mspct
#' @name Extract_mspct
#'
"[.generic_mspct" <-
  function(x, i, drop = NULL) {
    old.byrow <- attr(x, "mspct.byrow", exact = TRUE)
    if (is.null(old.byrow)) {
      old.byrow <- FALSE
    }
    old.class <- rmDerivedMspct(x)
    x <- `[`(x, i)
    class(x) <- c(old.class, class(x))
    attr(x, "mspct.dim") <- c(length(x), 1L)
    attr(x, "mspct.byrow") <- old.byrow
    attr(x, "mspct.version") <- 3
    x
  }

# Not exported
# Check if class_spct is compatible with class_mspct
#
is.member_class <- function(l, x) {
  class(l)[1] == "generic_mspct" && is.generic_spct(x) ||
    sub("_mspct", "", class(l)[1], fixed = TRUE) == sub("_spct", "", class(x)[1], fixed = TRUE)
}

#' @param value	A suitable replacement value: it will be repeated a whole number
#'   of times if necessary and it may be coerced: see the Coercion section. If
#'   NULL, deletes the column if a single column is selected.
#'
#' @export
#' @method [<- generic_mspct
#' @rdname extract_mspct
#'
"[<-.generic_mspct" <- function(x, i, value) {
  # could be improved to accept derived classes as valid for replacement.
  stopifnot(class(x) == class(value))
  # could not find a better way of avoiding infinite recursion as '[<-' is
  # a primitive with no explicit default method.
  old.byrow <- attr(x, "mspct.byrow", exact = TRUE)
  if (is.null(old.byrow)) {
    old.byrow <- FALSE
  }
  old.mspct.dim <- attr(x, "mspct.dim")
  old.class <- rmDerivedMspct(x)
  x[i] <- value
  class(x) <- c(old.class, class(x))
  attr(x, "mspct.dim") <- old.mspct.dim
  attr(x, "mspct.byrow") <- old.byrow
  attr(x, "mspct.version") <- 3
  x
}

#' @param name A literal character string or a name (possibly backtick quoted).
#'   For extraction, this is normally (see under 'Environments') partially
#'   matched to the names of the object.
#'
#' @export
#' @method $<- generic_mspct
#' @rdname extract_mspct
#'
"$<-.generic_mspct" <- function(x, name, value) {
  x[[name]] <- value
}

#' @export
#' @method [[<- generic_mspct
#' @rdname extract_mspct
#'
"[[<-.generic_mspct" <- function(x, name, value) {
  stopifnot(is.member_class(x, value) || is.null(value))
  # could not find a better way of avoiding infinite recursion as '[[<-' is
  # a primitive with no explicit default method.
  if (is.character(name) && !(name %in% names(x)) ) {
    if (ncol(x) == 1) {
      dimension <- c(nrow(x) + 1, 1)
    } else {
      stop("Appending to a matrix-like collection not supported.")
    }
  } else if (is.numeric(name) && (name > length(x)) ) {
    stop("Appending to a collection using numeric indexing not supported.")
  } else if (is.null(value)) {
    if (ncol(x) != 1) {
      stop("Deleting members from a matrix-like collection not supported.")
    } else {
      dimension <- attr(x, "mspct.dim", exact = TRUE)
      dimension[1] <- dimension[1] - 1L
    }
  } else {
    dimension <- attr(x, "mspct.dim", exact = TRUE)
  }
  old.byrow <- attr(x, "mspct.byrow", exact = TRUE)
  if (is.null(old.byrow)) {
    old.byrow <- FALSE
  }
  old.class <- rmDerivedMspct(x)
  x[[name]] <- value
  class(x) <- c(old.class, class(x))
  attr(x, "mspct.dim") <- dimension
  attr(x, "mspct.byrow") <- old.byrow
  attr(x, "mspct.version") <- 3
  x
}

# Combine -----------------------------------------------------------------

#' Combine collections of spectra
#'
#' Combine two or more generic_mspct objects into a single object.
#'
#' @param ... one or more generic_mspct objects to combine.
#' @param recursive logical ignored as nesting of collections of spectra is
#' not supported.
#' @param ncol numeric Virtual number of columns
#' @param byrow logical When object has two dimensions, how to map member
#' objects to columns and rows.
#'
#' @return A collection of spectra object belonging to the most derived class
#' shared among the combined objects.
#'
#' @name c
#'
#' @export
#' @method c generic_mspct
#'
c.generic_mspct <- function(..., recursive = FALSE, ncol = 1, byrow = FALSE) {
  l <- list(...)
  shared.class <- shared_member_class(l, target.set = mspct_classes())
  stopifnot(length(shared.class) > 0)
  shared.class <- shared.class[1]
  ul <- unlist(l, recursive = FALSE)
  do.call(shared.class, list(l = ul, ncol = ncol, byrow = byrow))
}
