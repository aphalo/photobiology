
# rbind -------------------------------------------------------------------

#' Makes one spectral object from a list or collection of many
#'
#' A wrapper on \code{dplyr::rbind_fill} that preserves class and other
#' attributes of spectral objects.
#'
#' @param l A \code{source_mspct}, \code{filter_mspct}, \code{reflector_mspct},
#'   \code{response_mspct}, \code{chroma_mspct}, \code{cps_mspct},
#'   \code{generic_mspct} object or a list containing \code{source_spct},
#'   \code{filter_spct}, \code{reflector_spct}, \code{response_spct},
#'   \code{chroma_spct}, \code{cps_spct}, or \code{generic_spct} objects.
#'
#' @param use.names logical If \code{TRUE} items will be bound by matching
#'   column names. By default \code{TRUE} for \code{rbindspct}. Columns with
#'   duplicate names are bound in the order of occurrence, similar to base. When
#'   TRUE, at least one item of the input list has to have non-null column
#'   names.
#'
#' @param fill logical If \code{TRUE} fills missing columns with NAs. By default
#'   \code{TRUE}. When \code{TRUE}, \code{use.names} has also to be \code{TRUE},
#'   and all items of the input list have to have non-null column names.
#'
#' @param idfactor logical or character Generates an index column of
#'   \code{factor} type. Default (\code{TRUE}) is to for both lists and
#'   \code{_mspct} objects. If \code{idfactor=TRUE} then the column is auto
#'   named \code{spct.idx}. Alternatively the column name can be directly
#'   provided to \code{idfactor} as a character string.
#'
#' @details Each item of \code{l} should be a spectrum, including \code{NULL}
#'   (skipped) or an empty object (0 rows). \code{rbindspc} is most useful when
#'   there are a variable number of (potentially many) objects to stack.
#'   \code{rbindspct} always returns at least a \code{generic_spct} as long as
#'   all elements in l are spectra.
#'
#' @note Note that any additional 'user added' attributes that might exist on
#'   individual items of the input list will not be preserved in the result.
#'   The attributes used by the \code{photobiology} package are preserved, and
#'   if they are not consistent accross the bound spectral objetcs, a warning is
#'   issued.
#'
#' @return An spectral object of a type common to all bound items containing a
#'   concatenation of all the items passed in. If the argument 'idfactor' is
#'   TRUE, then a factor 'spct.idx' will be added to the returned spectral
#'   object.
#'
#' @export
#'
#' @note \code{dplyr::rbind_fill} is called internally and the result returned is
#'   the highest class in the inheritance hierachy which is common to all
#'   elements in the list. If not all members of the list belong to one of the
#'   \code{_spct} classes, an error is triggered. The function sets all data in
#'   \code{source_spct} and \code{response_spct} objects supplied as arguments
#'   into energy-based quantities, and all data in \code{filter_spct} objects
#'   into transmittance before the row binding is done.
#'
#' @examples
#' spct <- rbindspct(list(sun.spct, sun.spct))
#' spct
#' class(spct)
#'
#' # adds factor 'spct.idx' with letters as levels
#' spct <- rbindspct(list(sun.spct, sun.spct), idfactor = TRUE)
#' head(spct)
#' class(spct)
#'
#' # adds factor 'spct.idx' with the names given to the spectra in the list
#' # supplied as formal argument 'l' as levels
#' spct <- rbindspct(list(one = sun.spct, two = sun.spct), idfactor = TRUE)
#' head(spct)
#' class(spct)
#'
#' # adds factor 'ID' with the names given to the spectra in the list
#' # supplied as formal argument 'l' as levels
#' spct <- rbindspct(list(one = sun.spct, two = sun.spct),
#'                   idfactor = "ID")
#' head(spct)
#' class(spct)
#'
rbindspct <- function(l, use.names = TRUE, fill = TRUE, idfactor = TRUE) {
  if ( (is.null(idfactor) && (!is.null(names(l)))) ||
       (is.logical(idfactor) && idfactor ) ) {
    idfactor <- "spct.idx"
  }
  add.idfactor <- is.character(idfactor)

  if (is.null(l) || length(l) < 1) {
    return(l)
  }
  if (!is.list(l) || is.any_spct(l) || is.waveband(l)) {
    stop("Argument 'l' should be a list of spectra or an _mspct object")
    return(NULL)
  }
  # we find the lowest common class
  # and in the same loop we make sure that all spectral data use consistent units
  l.class <- spct_classes()
  photon.based.input <- any(sapply(l, FUN = is_photon_based))
  absorbance.based.input <- any(sapply(l, FUN = is_absorbance_based))
  scaled.input <- sapply(l, FUN = is_scaled)
  normalized.input <- sapply(l, FUN = is_normalized)
  effective.input <- sapply(l, FUN = is_effective)
  if (any(scaled.input) && !all(scaled.input)) {
    warning("Spectra being row-bound have been differently re-scaled")
  }
  if (any(normalized.input) && length(unique(normalized.input)) > 1L) {
    warning("Spectra being row-bound have been differently normalized")
  }
  for (i in 1:length(l)) {
    class_spct <- class(l[[i]])
    l.class <- intersect(l.class, class_spct)
    if (photon.based.input && ("source_spct" %in% class_spct ||
        "response_spct" %in% class_spct )) {
      l[[i]] <- q2e(l[[i]], action = "replace", byref = FALSE)
    }
    if (absorbance.based.input && "filter_spct" %in% class_spct) {
      l[[i]] <- A2T(l[[i]], action = "replace", byref = FALSE)
    }
  }
  if (length(l.class) < 1L) {
    warning("Argument 'l' contains objects which are not spectra")
    return(NA)
  }
  # safer than l.class[1] which depends on order of class names in l.class
  if (length(l.class) > 1L) {
    l.class <- setdiff(l.class, "generic_spct")
  }
  #  print(l.class)

  # Here we do the actual binding
  if (length(l) < 2) {
    ans <- l[[1]]
  } else {
    ans <- plyr::rbind.fill(l)
    ans <- dplyr::as_data_frame(ans)
  }
  if (is.null(ans)) {
    return(NULL)
  }

  names.spct <- names(l)
  if (is.null(names.spct) || anyNA(names.spct) || length(names.spct) < length(l)) {
    names.spct <- paste("spct", 1:length(l), sep = "_")
  }
  if (add.idfactor) {
    ans[[idfactor]] <- factor(rep.int(names.spct, sapply(l, FUN = nrow)),
                                levels = names.spct)
  }

  comment.ans <- "rbindspct: concatenated comments"
  comments.found <- FALSE

  for (i in 1:length(l)) {
    temp <- comment(l[[i]])
    comments.found <- comments.found || !is.null(temp)
    if (add.idfactor) {
      temp <- paste("\n", idfactor , "= ", names.spct[i], ":\n", comment(l[[i]]), sep = "")
    } else {
      temp <- paste("\n spectrum = ", names.spct[i], ":\n", comment(l[[i]]), sep = "")
    }
    comment.ans <- paste(comment.ans, temp)
  }
  if (!comments.found) {
    comment.ans <- NULL
  }

  add.bswf <- FALSE

  if (l.class == "source_spct") {
    time.unit <- sapply(l, FUN = getTimeUnit)
    if (length(unique(time.unit)) > 1L) {
      warning("Inconsistent time units among source spectra in rbindspct")
      return(NA)
    }
    if (any(effective.input)) {
      bswfs.input <- sapply(l, FUN = getBSWFUsed)
      if (length(unique(bswfs.input)) > 1L) {
        add.bswf <- TRUE
        bswf.used <- "multiple"
        ans[["BSWF"]] <- factor(rep.int(bswfs.input, sapply(l, FUN = nrow)), levels = bswfs.input)
      } else {
        add.bswf <- FALSE
        bswf.used <- bswfs.input[1]
      }
    } else {
      add.bswf <- FALSE
      bswf.used <- rep("none", length(l))
    }
    setSourceSpct(ans, time.unit = time.unit[1], bswf.used = bswf.used, multiple.wl = length(l))
    if (photon.based.input) {
      e2q(ans, action = "add", byref = TRUE)
    }
  } else if (l.class == "filter_spct") {
    Tfr.type <- sapply(l, FUN = getTfrType)
    if (length(unique(Tfr.type)) > 1L) {
      warning("Inconsistent 'Tfr.type' among filter spectra in rbindspct")
      return(NA)
    }
    setFilterSpct(ans, Tfr.type = Tfr.type[1], multiple.wl = length(l))
    if (absorbance.based.input) {
      T2A(ans, action = "add", byref = TRUE)
    }
  } else if (l.class == "reflector_spct") {
    Rfr.type <- sapply(l, FUN = getRfrType)
    if (length(unique(Rfr.type)) > 1L) {
      warning("Inconsistent 'Rfr.type' among reflector spectra in rbindspct")
      return(NA)
    }
    setReflectorSpct(ans, Rfr.type = Rfr.type[1], multiple.wl = length(l))
  } else if (l.class == "response_spct") {
    time.unit <- sapply(l, FUN = getTimeUnit)
    if (length(unique(time.unit)) > 1L) {
      warning("Inconsistent time units among respose spectra in rbindspct")
      return(NA)
    }
    setResponseSpct(ans, time.unit = time.unit[1], multiple.wl = length(l))
    if (photon.based.input) {
      e2q(ans, action = "add", byref = TRUE)
    }
  } else if (l.class == "chroma_spct") {
    setChromaSpct(ans, multiple.wl = length(l))
  } else if (l.class == "cps_spct") {
    setCpsSpct(ans, multiple.wl = length(l))
  } else if (l.class == "generic_spct") {
    setGenericSpct(ans, multiple.wl = length(l))
  }
  if (any(scaled.input)) {
    attr(ans, "scaled") <- TRUE
  }
  if (any(normalized.input)) {
    attr(ans, "normalized") <- TRUE
  }
  if (!is.null(comment.ans)) {
    comment(ans) <- comment.ans
  }
  return(ans)
}

# Subset ------------------------------------------------------------------

# subset.data.frame works as expected with all spectral classes as it
# calls the Extract methods defined below on the object passed!

# Extract ------------------------------------------------------------------

# $ operator for extraction does not need any wrapping as it always extracts
# single columns returning objects of the underlying classes (e.g. numeric)
# rather than spectral objects.
#
# [ needs special handling as it can be used to extract rows, or groups of
# columns which are returned as spectral objects. Such returned objects
# can easily become invalid, for example, lack a w.length variable.

#' Extract or Replace Parts of a Spectrum
#'
#' Just like extraction and replacement with indexes in base R, but preserving
#' the special attributes used in spectral classes and checking for validity of
#' remaining spectral data.
#'
#' @param x	spectral object from which to extract element(s) or in which to replace element(s)
#' @param i index for rows,
#' @param j index for columns, specifying elements to extract or replace. Indices are
#'   numeric or character vectors or empty (missing) or NULL. Please, see
#'   \code{\link[base]{Extract.data.frame}} for more details.
#' @param drop logical. If TRUE the result is coerced to the lowest possible
#'   dimension. The default is FALSE unless the result is a single column.
#'
#' @details These methods are just wrappers on the method for data.frame objects
#'   which copy the additional attributes used by these classes, and validate
#'   the extracted object as a spectral object. When drop is TRUE and the
#'   returned object has only one column, then a vector is returned. If the
#'   extrated columns are more than one but do not include \code{w.length}, a
#'   data frame is returned instead of a spectral object.
#'
#' @return An object of the same class as \code{x} but containing only the
#'   subset of rows and columns that are selected. See details for special
#'   cases.
#'
#' @method "[" generic_spct
#'
#' @note Currently only extract methods are implemented. Replacement methods
#'   may be implemented in the future.
#'
#' @examples
#' sun.spct[sun.spct$w.length > 400, ]
#' subset(sun.spct, w.length > 400)
#'
#' tmp.spct <- sun.spct
#' tmp.spct[tmp.spct$s.e.irrad < 1e-5 , "s.e.irrad"] <- 0
#' e2q(tmp.spct[ , c("w.length", "s.e.irrad")]) # restore data consistency!
#'
#' @rdname extract
#' @name Extract
#'
#' @seealso \code{\link[base]{subset.data.frame}} and \code{\link{trim_spct}}
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
        setGenericSpct(xx, multiple.wl = getMultipleWl(x))
        setNormalized(xx, getNormalized(x))
        setScaled(xx, getScaled(x))
        comment <- comment(x)
        if (!is.null(comment)) {
          comment(xx) <- comment
        }
      } else {
        xx <- dplyr::as_data_frame(xx)
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
        setCpsSpct(xx, multiple.wl = getMultipleWl(x))
        setNormalized(xx, getNormalized(x))
        setScaled(xx, getScaled(x))
        comment <- comment(x)
        if (!is.null(comment)) {
          comment(xx) <- comment
        }
      } else {
        xx <- dplyr::as_data_frame(xx)
      }
    }
    xx
  }

#' @export
#' @rdname extract
#'
"[.source_spct" <-
  function(x, i, j, drop = NULL) {
    if (is.null(drop)) {
      xx <- `[.data.frame`(x, i, j)
    } else {
      xx <- `[.data.frame`(x, i, j, drop = drop)
    }
    if (is.data.frame(xx)) {
      if ("w.length" %in% names(xx)) {
        time.unit <- getTimeUnit(x)
        bswf.used <- getBSWFUsed(x)
        setSourceSpct(x = xx, time.unit = time.unit, bswf.used = bswf.used,
                      multiple.wl = getMultipleWl(x), strict.range = NULL)
        setNormalized(xx, getNormalized(x))
        setScaled(xx, getScaled(x))
        comment <- comment(x)
        if (!is.null(comment)) {
          comment(xx) <- comment
        }

      } else {
        xx <- dplyr::as_data_frame(xx)
      }
    }
    xx
  }

#' @export
#' @rdname extract
#'
"[.response_spct" <-
  function(x, i, j, drop = NULL) {
    if (is.null(drop)) {
      xx <- `[.data.frame`(x, i, j)
    } else {
      xx <- `[.data.frame`(x, i, j, drop = drop)
    }
    if (is.data.frame(xx)) {
      if ("w.length" %in% names(xx)) {
        time.unit <- getTimeUnit(x)
        setResponseSpct(x = xx, time.unit = time.unit,
                        multiple.wl = getMultipleWl(x))
        setNormalized(xx, getNormalized(x))
        setScaled(xx, getScaled(x))
        comment <- comment(x)
        if (!is.null(comment)) {
          comment(xx) <- comment
        }
      } else {
        xx <- dplyr::as_data_frame(xx)
      }
    }
    xx
  }

#' @export
#' @rdname extract
#'
"[.filter_spct" <-
  function(x, i, j, drop = NULL) {
    if (is.null(drop)) {
      xx <- `[.data.frame`(x, i, j)
    } else {
      xx <- `[.data.frame`(x, i, j, drop = drop)
    }
    if (is.data.frame(xx)) {
      if ("w.length" %in% names(xx)) {
        Tfr.type <- getTfrType(x)
        setFilterSpct(x = xx, Tfr.type = Tfr.type,
                      multiple.wl = getMultipleWl(x), strict.range = NULL)
        setNormalized(xx, getNormalized(x))
        setScaled(xx, getScaled(x))
        comment <- comment(x)
        if (!is.null(comment)) {
          comment(xx) <- comment
        }
      } else {
        xx <- dplyr::as_data_frame(xx)
      }
    }
    xx
  }

#' @export
#' @rdname extract
#'
"[.reflector_spct" <-
  function(x, i, j, drop = NULL) {
    if (is.null(drop)) {
      xx <- `[.data.frame`(x, i, j)
    } else {
      xx <- `[.data.frame`(x, i, j, drop = drop)
    }
    if (is.data.frame(xx)) {
      if ("w.length" %in% names(xx)) {
        Rfr.type <- getRfrType(x)
        setReflectorSpct(x = xx, Rfr.type = Rfr.type,
                         multiple.wl = getMultipleWl(x), strict.range = NULL)
        setNormalized(xx, getNormalized(x))
        setScaled(xx, getScaled(x))
        comment <- comment(x)
        if (!is.null(comment)) {
          comment(xx) <- comment
        }
      } else {
        xx <- dplyr::as_data_frame(xx)
      }
    }
    xx
  }

#' @export
#' @rdname extract
#'
"[.object_spct" <-
  function(x, i, j, drop = NULL) {
    if (is.null(drop)) {
      xx <- `[.data.frame`(x, i, j)
    } else {
      xx <- `[.data.frame`(x, i, j, drop = drop)
    }
    if (is.data.frame(xx)) {
      if ("w.length" %in% names(xx)) {
        Tfr.type <- getTfrType(x)
        Rfr.type <- getRfrType(x)
        setObjectSpct(x = xx, Tfr.type = Tfr.type, Rfr.type = Rfr.type,
                      multiple.wl = getMultipleWl(x), strict.range = NULL)
        setNormalized(xx, getNormalized(x))
        setScaled(xx, getScaled(x))
        comment <- comment(x)
        if (!is.null(comment)) {
          comment(xx) <- comment
        } else {
          xx <- dplyr::as_data_frame(xx)
        }
      }
    }
    xx
  }

#' @export
#' @rdname extract
#'
"[.chroma_spct" <-
  function(x, i, j, drop = NULL) {
    if (is.null(drop)) {
      xx <- `[.data.frame`(x, i, j)
    } else {
      xx <- `[.data.frame`(x, i, j, drop = drop)
    }
    if (is.data.frame(xx)) {
      if ("w.length" %in% names(xx)) {
        setChromaSpct(xx, multiple.wl = getMultipleWl(x))
        setNormalized(xx, getNormalized(x))
        setScaled(xx, getScaled(x))
        comment <- comment(x)
        if (!is.null(comment)) {
          comment(xx) <- comment
        }
      } else {
        xx <- dplyr::as_data_frame(xx)
      }
    }
    xx
  }


# replace -----------------------------------------------------------------

# We need to wrap the replace functions adding a call to our check method
# to make sure that the object is still a valid spectrum after the
# replacement.

#' @param value	A suitable replacement value: it will be repeated a whole number
#'   of times if necessary and it may be coerced: see the Coercion section. If
#'   NULL, deletes the column if a single column is selected.
#'
#' @export
#' @method "[<-" generic_spct
#' @rdname extract
#'
"[<-.generic_spct" <- function(x, i, j, value) {
  check(`[<-.data.frame`(x, i, j, value), byref = FALSE)
}

#' @export
#' @method "$<-" generic_spct
#' @rdname extract
#'
"$<-.generic_spct" <- function(x, name, value) {
  check(`$<-.data.frame`(x, name, value), byref = FALSE)
}
