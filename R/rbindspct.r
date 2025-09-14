
# rbind -------------------------------------------------------------------

#' Row-bind spectra
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
#'   \code{factor} type. Default is (\code{idfactor=TRUE}) for both lists and
#'   \code{_mspct} objects. If \code{idfactor=TRUE} then the column is auto
#'   named \code{spct.idx}. Alternatively the column name can be directly
#'   provided to \code{idfactor} as a character string.
#'
#' @param attrs.source integer Index into the members of the list from which
#'   attributes should be copied. If \code{NULL}, all attributes are collected
#'   into named lists, except that unique comments are pasted.
#'
#' @param attrs.simplify logical Flag indicating that when all values of an
#'   attribute are equal for all members, the named list will be replaced by
#'   a single copy of the value.
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
#'   if they are not consistent across the bound spectral objects, a warning is
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
#'   the highest class in the inheritance hierarchy which is common to all
#'   elements in the list. If not all members of the list belong to one of the
#'   \code{_spct} classes, an error is triggered. The function sets all data in
#'   \code{source_spct} and \code{response_spct} objects supplied as arguments
#'   into energy-based quantities, and all data in \code{filter_spct} objects
#'   into transmittance before the row binding is done. If any member spectrum
#'   is tagged, it is untagged before row binding.
#'
#' @examples
#' # default, adds factor 'spct.idx' with letters as levels
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
rbindspct <- function(l,
                      use.names = TRUE,
                      fill = TRUE,
                      idfactor = TRUE,
                      attrs.source = NULL,
                      attrs.simplify = FALSE) {
  if (is.null(l) || !is.list(l) || length(l) < 1) {
    # _mspct classes are derived from "list"
    warning("Argument 'l' should be a non-empty list or ",
            "a collection of spectra.")
    return(generic_spct())
  }

  if ((is.null(idfactor) &&
         (!is.null(names(l)))) || is.logical(idfactor) && idfactor) {
    idfactor <- "spct.idx"
  } else if (is.logical(idfactor) && !idfactor) {
    idfactor <- NULL
  }

  # inefficient but simpler to implement, and ensures proper naming
  # make sure each member spct object contains a single spectrum
  if (any(sapply(l, getMultipleWl) > 1L)) {
    l <- subset2mspct(l)
  }

  # we skip spectra with no rows
  selector <- unname(sapply(l, nrow)) > 0

  if (use.names && !rlang::is_named(l)) {
    names(l) <- paste("spct", seq_along(l), sep = "_")
  }
  add.idfactor <- is.character(idfactor)

  # We find the most derived common class for spectra
  l.class <- shared_member_class(l)
  if (length(l.class) < 1L) {
    stop("Argument 'l' should contain spectra.")
  } else {
    l.class <- l.class[1L]
  }
  if (!any(selector)) {
    return(do.call(what = l.class, args = list()))
  }
  if (length(l[selector]) == 1L) {
    z <- l[selector][[1L]]
    if (add.idfactor) {
      z[[idfactor]] <- factor(rep(names(l[selector]), times = nrow(z)))
      setIdFactor(z, idfactor)
    }
    return(z)
  }

  # list may have members which already have multiple spectra in long form
  mltpl.wl <- sum(sapply(l, FUN = getMultipleWl))

  # we check that all spectral data contain consistent quantities
  if (l.class %in% c("source_spct", "response_spct")) {
    photon.based <- sapply(l, FUN = is_photon_based)
    energy.based <- sapply(l, FUN = is_energy_based)
    qe.consistent.based <-
      all(photon.based) && !any(energy.based) ||
      all(energy.based) && !any(photon.based) ||
      all(energy.based) && all(photon.based)
  } else {
    qe.consistent.based <- NA
  }

  if (l.class == "filter_spct") {
    absorbance.based <- sapply(l, FUN = is_absorbance_based)
    transmittance.based <- sapply(l, FUN = is_transmittance_based)
    absorptance.based <- sapply(l, FUN = is_absorptance_based)
    TA.consistent.based <- all(absorbance.based) ||
      all(absorptance.based) ||
      all(transmittance.based)
  } else {
    TA.consistent.based <- NA
  }

  # check for transformed data
  scaled.input <- sapply(l, FUN = is_scaled)
  normalized.input <- sapply(l, FUN = is_normalized)
  effective.input <- sapply(l, FUN = is_effective)

  if (any(scaled.input) && !all(scaled.input)) {
    warning("Spectra being row-bound have been differently re-scaled")
  }
  if (any(normalized.input) && length(unique(normalized.input)) > 1L) {
    warning("Spectra being row-bound have been differently normalized")
  }

  for (i in seq_along(l)) {
    class_spct <- class(l[[i]])[1]
    l.class <- intersect(l.class, class_spct)
    if (is_tagged(l[[i]])) {
      l[[i]] <- untag(l[[i]])
    }
    if (!is.na(qe.consistent.based) && !qe.consistent.based) {
      l[[i]] <- q2e(l[[i]], action = "replace", byref = FALSE)
    }
    if (!is.na(TA.consistent.based) && !TA.consistent.based) {
      l[[i]] <- A2T(l[[i]], action = "replace", byref = FALSE)
    }
  }

  # check class is same for all spectra
  #  print(l.class)
  if (length(l.class) != 1L) {
    stop("All spectra in 'l' should belong to the same spectral class.")
  }

  # Here we do the actual binding
  if (length(l) == 1) {
    ans <- l[[1]]
  } else {
    ans <- plyr::rbind.fill(l)
    ans <- tibble::as_tibble(ans)
  }
  if (is.null(ans)) {
    return(generic_spct())
  }

  names.spct <- names(l)
  if (is.null(names.spct) || anyNA(names.spct) || length(names.spct) < length(l)) {
    names.spct <- paste("spct", seq_along(l), sep = "_")
  } else {
    if (anyDuplicated(names.spct)) {
      warning("Duplicated member names have been de-ambiguated before binding spectra.")
      names.spct <- make.unique(names.spct, sep = "_")
      names(l) <- names.spct
    }
  }
  if (add.idfactor) {
    ans[[idfactor]] <- factor(rep(names.spct, times = sapply(l, FUN = nrow)),
                              levels = names.spct)
  }

  comment.ans <- "rbindspct: concatenated comments"
  comments.found <- FALSE

  if (length(attrs.source)) {
    idxs <- intersect(seq_along(l), attrs.source)
  } else {
    idxs <- seq_along(l)
  }

  # get methods and functions return NA if attr is not set
  if (length(idxs) == 1L) {
    comment.ans <- comment(l[[idxs]])
    instr.desc <- getInstrDesc(l[[idxs]])
    instr.settings <- getInstrSettings(l[[idxs]])
    when.measured <- getWhenMeasured(l[[idxs]])
    where.measured <- getWhereMeasured(l[[idxs]])
    what.measured <- getWhatMeasured(l[[idxs]])
    how.measured <- getHowMeasured(l[[idxs]])
    normalized <- getNormalized(l[[idxs]])
    normalization <- getNormalization(l[[idxs]])
  } else {
    # we avoid duplicating the attributes when possible
    comments <- lapply(l[idxs], comment)
    comment.ans <- paste(unique(comments))

    instr.desc <- lapply(l[idxs], getInstrDesc)
    if (attrs.simplify && length(unique(instr.desc)) == 1) {
      instr.desc <- instr.desc[[1]]
    } else {
      names(instr.desc) <- names.spct[idxs]
    }

    instr.settings <- lapply(l[idxs], getInstrSettings)
    if (attrs.simplify && length(unique(instr.settings)) == 1) {
      instr.settings <- instr.settings[[1]]
    } else {
      names(instr.settings) <- names.spct[idxs]
    }

    when.measured <- lapply(l[idxs], getWhenMeasured)
    names(when.measured) <- names.spct[idxs]

    where.measured <-
      dplyr::bind_rows(lapply(l[idxs], getWhereMeasured), .id = idfactor)
    if (attrs.simplify &&
          (all(is.na(where.measured$lon)) ||
             length(unique(where.measured$lon)) == 1) &&
          (all(is.na(where.measured$lat)) ||
             length(unique(where.measured$lat)) == 1) &&
          (all(is.na(where.measured$address)) ||
             length(unique(where.measured$address)) == 1)) {
      where.measured <- where.measured[1, ]
    # do not remove columns by numeric index!!
    }

    what.measured <- lapply(l[idxs], getWhatMeasured)
    if (attrs.simplify && length(unique(what.measured)) == 1) {
      what.measured <- what.measured[[1]]
    } else {
      names(what.measured) <- names.spct[idxs]
    }

    how.measured <- lapply(l[idxs], getHowMeasured)
    if (attrs.simplify && length(unique(how.measured)) == 1) {
      how.measured <- how.measured[[1]]
    } else {
      names(how.measured) <- names.spct[idxs]
    }

    normalized <- lapply(l[idxs], getNormalized)
    names(normalized) <- names.spct[idxs]

    normalization <- lapply(l[idxs], getNormalization)
    names(normalization) <- names.spct[idxs]

  }

  if (l.class == "source_spct") {
    time.unit <- sapply(l, FUN = getTimeUnit)
    names(time.unit) <- NULL
    time.unit <- unique(time.unit)
    if (length(time.unit) > 1L) {
      warning("Inconsistent time units among source spectra ",
              "passed to rbindspct")
      return(source_spct())
    }
    if (any(effective.input)) {
      bswfs.input <- sapply(l, FUN = getBSWFUsed)
      if (length(unique(bswfs.input)) > 1L) {
        bswf.used <- "multiple"
        ans[["BSWF"]] <-
          factor(rep(bswfs.input, times = sapply(l, FUN = nrow)),
                 levels = bswfs.input)
      } else {
        bswf.used <- bswfs.input[1]
      }
    } else {
      bswf.used <- "none"
    }
    setSourceSpct(ans,
                  time.unit = time.unit[1],
                  bswf.used = bswf.used,
                  multiple.wl = mltpl.wl)
    if (!qe.consistent.based) {
      e2q(ans, action = "add", byref = TRUE)
    }
  } else if (l.class == "filter_spct") {
    Tfr.type <- sapply(l, FUN = getTfrType)
    names(Tfr.type) <- NULL
    Tfr.type <- unique(Tfr.type)
    if (length(Tfr.type) > 1L) {
      warning("Inconsistent 'Tfr.type' among filter spectra ",
              "passed to rbindspct")
      return(filter_spct())
    }
    filter.descriptor <-
      sapply(l, FUN = getFilterProperties, return.null = TRUE)
    # TODO merge it if possible
    # and then set
    setFilterSpct(ans, Tfr.type = Tfr.type[1], multiple.wl = mltpl.wl)
    if (!TA.consistent.based) {
      T2A(ans, action = "add", byref = TRUE)
    }
  } else if (l.class == "reflector_spct") {
    Rfr.type <- sapply(l, FUN = getRfrType)
    names(Rfr.type) <- NULL
    Rfr.type <- unique(Rfr.type)
    if (length(Rfr.type) > 1L) {
      warning("Inconsistent 'Rfr.type' among reflector spectra in rbindspct")
      return(reflector_spct())
    }
    setReflectorSpct(ans, Rfr.type = Rfr.type[1], multiple.wl = mltpl.wl)
  } else if (l.class == "object_spct") {
    Tfr.type <- sapply(l, FUN = getTfrType)
    names(Tfr.type) <- NULL
    Tfr.type <- unique(Tfr.type)
    Rfr.type <- sapply(l, FUN = getRfrType)
    names(Rfr.type) <- NULL
    Rfr.type <- unique(Rfr.type)
    if (length(Tfr.type) > 1L) {
      warning("Inconsistent 'Tfr.type' among filter spectra ",
              "passed to rbindspct")
      return(filter_spct())
    }
    if (length(Rfr.type) > 1L) {
      warning("Inconsistent 'Rfr.type' among reflector spectra ",
              "passed to rbindspct")
      return(reflector_spct())
    }
    setObjectSpct(ans, Tfr.type = Tfr.type[1], Rfr.type = Rfr.type[1],
                  multiple.wl = mltpl.wl)
  } else if (l.class == "solute_spct") {
    K.type <- sapply(l, FUN = getKType)
    names(K.type) <- NULL
    K.type <- unique(K.type)
    if (length(K.type) > 1L) {
      warning("Inconsistent 'K.type' among solute spectra in rbindspct")
      return(reflector_spct())
    }
    setSoluteSpct(ans, K.type = K.type, multiple.wl = mltpl.wl)
  } else if (l.class == "response_spct") {
    time.unit <- sapply(l, FUN = getTimeUnit)
    names(time.unit) <- NULL
    time.unit <- unique(time.unit)
    if (length(time.unit) > 1L) {
      warning("Inconsistent time units among response spectra in rbindspct")
      return(response_spct())
    }
    setResponseSpct(ans, time.unit = time.unit[1], multiple.wl = mltpl.wl)
    if (!qe.consistent.based) {
      e2q(ans, action = "add", byref = TRUE)
    }
  } else if (l.class == "chroma_spct") {
    setChromaSpct(ans, multiple.wl = mltpl.wl)
  } else if (l.class == "cps_spct") {
    setCpsSpct(ans, multiple.wl = mltpl.wl)
  } else if (l.class == "raw_spct") {
    setRawSpct(ans, multiple.wl = mltpl.wl)
  } else if (l.class == "generic_spct") {
    setGenericSpct(ans, multiple.wl = mltpl.wl)
  }
  if (any(scaled.input)) {
    attr(ans, "scaled") <- TRUE
  }
  if (!is.null(comment.ans)) {
    comment(ans) <- comment.ans
  }
  if (is.character(idfactor)) {
    setIdFactor(ans, idfactor)
  }
  setWhenMeasured(ans, when.measured)
  setWhereMeasured(ans, where.measured)
  setWhatMeasured(ans, what.measured)
  setHowMeasured(ans, how.measured)
  attr(ans, "normalized") <- normalized
  if (any(normalized.input)) {
    attr(ans, "normalization") <- normalization
  }
  if (!all(is.na(instr.desc))) {
    setInstrDesc(ans, instr.desc)
  }
  if (!all(is.na(instr.settings))) {
    setInstrSettings(ans, instr.settings)
  }
  ans
}

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
    attr(x, "mspct.version") <- 2
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
  attr(x, "mspct.version") <- 2
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
  attr(x, "mspct.version") <- 2
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
