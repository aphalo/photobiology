
# rbind -------------------------------------------------------------------

#' Makes one spectral object from a list of many
#'
#' Same as \code{rbindlist} from package data.table but preserves class of
#' spectral objects. Has different defaults for use names and fill.
#'
#' @param l A list containing \code{source_spct}, \code{filter_spct},
#'   \code{reflector_spct}, \code{response_spct}, \code{chroma_spct},
#'   \code{cps_spct}, \code{generic_spct}, \code{data.table}, \code{data.frame}
#'   or \code{list} objects. At least one of the inputs should have column names
#'   set. \code{\dots} is the same but you pass the objects by name separately.
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
#'   \code{factor} type. Default (\code{FALSE}) is not to. If
#'   \code{idfactor=TRUE} then the column is auto named \code{spct.idx}.
#'   Alternatively the column name can be directly provided to \code{idfactor}.
#'
#' @details Each item of \code{l} can be a spectrum, \code{data.table},
#' \code{data.frame} or \code{list}, including \code{NULL} (skipped) or an empty
#' object (0 rows). \code{rbindspc} is most useful when there are a variable
#' number of (potentially many) objects to stack. \code{rbind} (not implemented
#' yet for spectra) however is most useful to stack two or three objects which
#' you know in advance. \code{rbindspct} always returns at least a
#' \code{generic_spct} as long as all elements in l are spectra, otherwise a
#' \code{data.frame} is returned even when stacking a \code{list} with a
#' \code{data.frame}, for example. The difference between \code{rbindspct(l)}
#' and \code{rbindlist(l)} from package \pkg{data.table} is in their
#' \emph{default value for formal argument} \code{use.names}, and in that
#' \code{rbindlist} will NOT return a spct object even when the list l contains
#' only spct objects. In other words it drops derived classes, so its use should
#' be avoided for spectral objects, and \code{rbindspct(l)} should be always
#' used when working with spectral objects.
#'
#' @note Note that any additional 'user added' attributes that might exist on
#'   individual items of the input list would not be preserved in the result.
#'   The attributes used by the \code{photobiology} package are preserved, and
#'   if they are not consistent accross the bound spectral objetcs, a warning is
#'   issued.
#'
#' @return An spectral object of a type common to all bound items or a
#'   \code{data.table} containing a concatenation of all the items passed in. If
#'   the argument 'add.factor' is true, then a factor 'spct.idx' will be added
#'   to the returned spectral object.
#'
#' @export
#'
#' @seealso  \code{\link{data.table}}
#'
#' @note data.table::rbindlist is called internally and the result returned is
#'   the highest class in the inheritance hierachy which is common to all
#'   elements in the list. If not all members of the list belong to one of the
#'   \code{.spct} classes, an error is triggered. The function sets all data in
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
rbindspct <- function(l, use.names = TRUE, fill = TRUE, idfactor = NULL) {
  # original rbindlist from data.table strips attributes and sets class to data.table
  if (is.null(l) || length(l) < 1) {
    return(l)
  }
  if (!is.list(l) || is.any_spct(l) || is.waveband(l)) {
    stop("Argument 'l' should be a list of spectra")
    return(NULL)
  }
  # we find the lowest common class
  # and in the same loop we make sure that all spectral data uses consistent units
  l.class <- c( "source_spct", "filter_spct", "reflector_spct", "response_spct", "chroma_spct",
                "generic_spct")
  photon.based.input <- any(sapply(l, FUN=is_photon_based))
  absorbance.based.input <- any(sapply(l, FUN=is_absorbance_based))
  scaled.input <- sapply(l, FUN = is_scaled)
  normalized.input <- sapply(l, FUN = is_normalized)
  effective.input <- sapply(l, FUN = is_effective)
  if (any(scaled.input) && !all(scaled.input)) {
    warning("Only some of the spectra being row-bound have been previously scaled")
  }
  if (any(normalized.input) && length(unique(normalized.input)) > 1L) {
    warning("Only some of the spectra being row-bound have been previously normalized")
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
  for (spct in l) {
    l.class <- intersect(l.class, class(spct))
  }
  if (length(l.class) < 1L) {
    warning("Argument 'l' contains objects which are not spectra")
    return(NA)
  }
  l.class <- l.class[1]
  #  print(l.class)

    names.spct <- names(l)
  if (is.null(names.spct) || anyNA(names.spct)) {
    names.spct <- LETTERS[1:length(l)]
  }

  # Here we do the actual binding
  if (length(l) < 2) {
    ans <- l[[1]]
  } else {
    ans <- data.table::rbindlist(l, use.names, fill)
  }
  if (is.null(ans)) {
    return(NULL)
  }

  add.factor <- !is.null(idfactor)
  if (add.factor) {
    if (is.character(idfactor)) {
      factor.name <- idfactor
    } else {
      factor.name <- "spct.idx"
    }
    ans[ , (factor.name) := factor(rep.int(names.spct, sapply(l, FUN = nrow)), levels = names.spct)]
  }

  comment.ans <- "rbindspct: concatenated comments"
  comments.found <- FALSE

  for (i in 1:length(l)) {
    temp <- comment(l[[i]])
    comments.found <- comments.found || !is.null(temp)
    if (add.factor) {
      temp <- paste("\n", factor.name , "= ", names.spct[i], ":\n", comment(l[[i]]), sep="")
    } else {
      temp <- paste("\n spectrum = ", names.spct[i], ":\n", comment(l[[i]]), sep="")
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
        bswf.used <- "unknown"
        ans[ , BSWF := factor(rep.int(bswfs.input, sapply(l, FUN = nrow)), levels = bswfs.input)]
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
    setChromSpct(ans, multiple.wl = length(l))
  } else if (l.class == "cps_spct") {
    setCpsSpct(ans, multiple.wl = length(l))
  } else if (l.class == "generic_spct") {
    setGenericSpct(ans, multiple.wl = length(l))
  }
  if (any(scaled.input)) {
    setattr(ans, "scaled", "TRUE")
  }
  if (any(normalized.input)) {
    setattr(ans, "normalized", "TRUE")
  }
  if (add.factor && !add.bswf) {
    keys <- c(factor.name, "w.length")
    setkeyv(ans, keys)
  } else if (!add.factor && add.bswf) {
    keys <- c("BSWF", "w.length")
    setkeyv(ans, keys)
  } else if (add.factor && add.bswf) {
    keys <- c("BSWF", factor.name, "w.length")
    setkeyv(ans, keys)
  } # else we keep the default "w.length"
  if (!is.null(comment.ans)) setattr(ans, "comment", comment.ans)
  return(ans)
}


# subset ------------------------------------------------------------------


#' Subsetting methods for spectra
#'
#' Just like \code{subset} in base R, but preserves the special attributes used
#' in spectral classes.
#'
#' @param x	generic_spct to subset
#' @param subset logical expression indicating elements or rows to keep
#' @param select expression indicating columns to select from x (IGNORED)
#' @param idx integer vector of indexes of rows to keep
#' @param ...	further arguments to be passed to or from other methods
#'
#' @details The subset argument works on the rows and will be evaluated in the
#'   generic_spct so columns can be referred to (by name) as variables in the
#'   expression The generic_spct that is returned will maintain the original
#'   attributes and keys as long as they are not select-ed out.
#'
#' @return An object of the same class as \code{x} but containing only the
#'   subset of rows and columns that are selected.
#'
#' @method subset generic_spct
#'
#' @note Current implementation restores object class after subsetting, as
#' subsetting using subset.data.table() strips the object of attributes
#' including derived classes.
#'
#' @examples
#' subset(sun.spct, w.length > 400)
#'
#' @seealso \code{\link{subset}} and \code{\link{trim_spct}}
#'
subset.generic_spct <- function(x, subset, select, idx = NULL, ...) {
  comment <- comment(x)
  if (!is.null(idx)) {
    z <- x[idx]
  } else {
    z <- x[eval(substitute(subset))]
  }
  setGenericSpct(z)
  if (!is.null(comment)) {
    setattr(z, "comment", comment)
  }
  untag(z)
  return(z)
}

# @describeIn subset.generic_spct Subset for counts per second spectra.
#'
#' @export
#' @rdname subset.generic_spct
#'
subset.cps_spct <- function(x, subset, select = NULL, idx = NULL, ...) {
  comment <- comment(x)
  if (!is.null(idx)) {
    z <- x[idx]
  } else {
    z <- x[eval(substitute(subset))]
  }
  setCpsSPct(z)
  if (!is.null(comment)) {
    setattr(z, "comment", comment)
  }
  untag(z)
  return(z)
}

# @describeIn subset.generic_spct Subset for light source spectra.
#'
#' @export
#' @rdname subset.generic_spct
#'
subset.source_spct <- function(x, subset, select = NULL, idx = NULL, ...) {
  time.unit <- getTimeUnit(x)
  bswf.used <- getBSWFUsed(x)
  comment <- comment(x)
  if (!is.null(idx)) {
    z <- x[idx]
  } else {
    z <- x[eval(substitute(subset))]
  }
  setSourceSpct(x = z, time.unit = time.unit, bswf.used = bswf.used, multiple.wl = Inf)
  if (!is.null(comment)) {
    setattr(z, "comment", comment)
  }
  untag(z)
  return(z)
}

# @describeIn subset.generic_spct Subset for light filter spectra.
#'
#' @export
#' @rdname subset.generic_spct
#'
subset.filter_spct <- function(x, subset, select = NULL, idx = NULL, ...) {
  Tfr.type <- getTfrType(x)
  comment <- comment(x)
  if (!is.null(idx)) {
    z <- x[idx]
  } else {
    z <- x[eval(substitute(subset))]
  }
  setFilterSpct(x = z, Tfr.type = Tfr.type, multiple.wl = Inf)
  if (!is.null(comment)) {
    setattr(z, "comment", comment)
  }
  untag(z)
  return(z)
}

# @describeIn subset.generic_spct Subset for light reflector spectra.
#'
#' @export
#' @rdname subset.generic_spct
#'
subset.reflector_spct <- function(x, subset, select = NULL, idx = NULL, ...) {
  Rfr.type <- getRfrType(x)
  comment <- comment(x)
  if (!is.null(idx)) {
    z <- x[idx]
  } else {
    z <- x[eval(substitute(subset))]
  }
  setReflectorSpct(x = z, Rfr.type = Rfr.type, multiple.wl = Inf)
  if (!is.null(comment)) {
    setattr(z, "comment", comment)
  }
  untag(z)
  return(z)
}

# @describeIn subset.generic_spct Subset for light response spectra.
#'
#' @export
#' @rdname subset.generic_spct
#'
subset.response_spct <- function(x, subset, select = NULL, idx = NULL, ...) {
  time.unit <- getTimeUnit(x)
  comment <- comment(x)
  if (!is.null(idx)) {
    z <- x[idx]
  } else {
    z <- x[eval(substitute(subset))]
  }
  setResponseSpct(x = z, time.unit = time.unit, multiple.wl = Inf)
  if (!is.null(comment)) {
    setattr(z, "comment", comment)
  }
  untag(z)
  return(z)
}

# @describeIn subset.generic_spct Subset for light response spectra.
#'
#' @export
#' @rdname subset.generic_spct
#'
subset.object_spct <- function(x, subset, select = NULL, idx = NULL, ...) {
  Tfr.type <- getTfrType(x)
  Rfr.type <- getRfrType(x)
  comment <- comment(x)
  if (!is.null(idx)) {
    z <- x[idx]
  } else {
    z <- x[eval(substitute(subset))]
  }
  setObjectSpct(x = z, Tfr.type = Tfr.type, Rfr.type = Rfr.type, multiple.wl = Inf)
  if (!is.null(comment)) {
    setattr(z, "comment", comment)
  }
  untag(z)
  return(z)
}


