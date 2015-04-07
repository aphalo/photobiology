
# rbind -------------------------------------------------------------------

#' Makes one spectral object from a list of many
#'
#' Same as \code{rbindlist} from package data.table but preserves class of spectral objects. Has different defaults
#' for use names and fill.
#'
#' @usage rbindspct(l, use.names = TRUE, fill=TRUE, idfactor = NULL)
#'
#' @param l A list containing \code{source.spct}, \code{filter.spct}, \code{reflector.spct}, \code{response.spct},
#' \code{chroma.spct}, \code{generic.spct}, \code{data.table}, \code{data.frame} or \code{list} objects.
#' At least one of the inputs should have column names set. \code{\dots} is the same but you pass the objects by name separately.
#'
#' @param use.names If \code{TRUE} items will be bound by matching column names. By default \code{TRUE} for
#' \code{rbindspct}. Columns with duplicate names are bound in the order of occurrence, similar to base.
#' When TRUE, at least one item of the input list has to have non-null column names.
#'
#' @param fill If \code{TRUE} fills missing columns with NAs. By default \code{TRUE}. When \code{TRUE},
#' \code{use.names} has also to be \code{TRUE}, and all items of the input list have to have non-null column names.
#'
#' @param idfactor Generates an index column of \code{factor} type. Default (\code{FALSE}) is not to.
#' If \code{idfactor=TRUE} then the column is auto named \code{spct.idx}. Alternatively the column name can be
#' directly provided to \code{idfactor}.
#'
#' @details
#' Each item of \code{l} can be a spectrum, \code{data.table}, \code{data.frame} or \code{list}, including \code{NULL} (skipped)
#' or an empty object (0 rows). \code{rbindspc} is most useful when there are a variable number of (potentially many)
#' objects to stack. \code{rbind} (not implemented yet for spectra) however is most useful to
#' stack two or three objects which you know in advance. \code{rbindspct} always
#' returns at least a \code{generic.spct} as long as all elements in l are spectra, otherwise a \code{data.frame}
#' is returned even when stacking a \code{list} with a \code{data.frame}, for example.
#  The difference between \code{rbindspct(l)} and \code{rbindlist(l)} from package \code{data.table} is in their
#' \emph{default value for formal argument} \code{use.names}, and in that \code{rbindlist} will NOT return an spct
#' object even when the list l contains only spct objects. In other words it drops derived classes, so its use
#' should be avoided for spectral objects, and \code{rbindspct(l)} should be always used when working with
#' spectral objects.
#'
#' Note that any additional 'user added' attributes that might exist on individual items of the input list would not
#' be preserved in the result. The attributes used by the \code{photobiology} package are preserved, and if they are
#' not consistent accross the bound spectral objetcs, a warning is issued.
#'
#' @return An spectral object of a type common to all bound items or a \code{data.table} containing a concatenation of
#' all the items passed in. If the argument 'add.factor' is true, then a factor 'spct.idx' will be added to the
#' returned spectral object.
#'
#' @export
#'
#' @seealso  \code{\link{data.table}}
#'
#' @note data.table::rbindlist is called internally and the result returned is the highest class in the inheritance
#' hierachy which is common to all elements in the list. If not all members of the list belong to one of the
#' \code{.spct} classes, an error is triggered. The function sets all data in \code{source.spct} and \code{response.spct}
#' objects supplied as arguments into energy-based quantities, and all data in \code{filter.spct} objects into
#' transmittance before the row binding is done.
#'
#' @examples
#'
#' # examples for spectra
#'
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
#' @keywords data
#'

rbindspct <- function(l, use.names = TRUE, fill = TRUE, idfactor = NULL) {
  # original rbindlist from data.table strips attributes and sets class to data.table
  if (is.null(l) || length(l) < 1) {
    return(l)
  }
  if (!is.list(l) || is.any.spct(l) || is.waveband(l)) {
    stop("Argument 'l' should be a list of spectra")
    return(NULL)
  }
  # we find the lowest common class
  # and in the same loop we make sure that all spectral data uses consistent units
  l.class <- c( "source.spct", "filter.spct", "reflector.spct", "response.spct", "chroma.spct",
                "generic.spct")
  photon.based.input <- any(sapply(l, FUN=is.photon.based))
  absorbance.based.input <- any(sapply(l, FUN=is.absorbance.based))
  rescaled.input <- sapply(l, FUN = is.rescaled)
  normalized.input <- sapply(l, FUN = is.normalized)
  effective.input <- sapply(l, FUN = is.effective)
  if (any(rescaled.input) && !all(rescaled.input)) {
    warning("Only some of the spectra being row-bound have been previously rescaled")
  }
  if (any(normalized.input) && length(unique(normalized.input)) > 1L) {
    warning("Only some of the spectra being row-bound have been previously normalized")
  }
  for (i in 1:length(l)) {
    class.spct <- class(l[[i]])
    l.class <- intersect(l.class, class.spct)
    if (photon.based.input && ("source.spct" %in% class.spct ||
        "response.spct" %in% class.spct )) {
      l[[i]] <- q2e(l[[i]], action = "replace", byref = FALSE)
    }
    if (absorbance.based.input && "filter.spct" %in% class.spct) {
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

  if (l.class == "source.spct") {
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
    setSourceSpct(ans, time.unit = time.unit[1], bswf.used = bswf.used)
    if (photon.based.input) {
      e2q(ans, action = "add", byref = TRUE)
    }
  } else if (l.class == "filter.spct") {
    Tfr.type <- sapply(l, FUN = getTfrType)
    if (length(unique(Tfr.type)) > 1L) {
      warning("Inconsistent 'Tfr.type' among filter spectra in rbindspct")
      return(NA)
    }
    setFilterSpct(ans, Tfr.type = Tfr.type[1])
    if (absorbance.based.input) {
      T2A(ans, action = "add", byref = TRUE)
    }
  } else if (l.class == "reflector.spct") {
    Rfr.type <- sapply(l, FUN = getRfrType)
    if (length(unique(Rfr.type)) > 1L) {
      warning("Inconsistent 'Rfr.type' among reflector spectra in rbindspct")
      return(NA)
    }
    setReflectorSpct(ans, Rfr.type = Rfr.type[1])
  } else if (l.class == "response.spct") {
    time.unit <- sapply(l, FUN = getTimeUnit)
    if (length(unique(time.unit)) > 1L) {
      warning("Inconsistent time units among respose spectra in rbindspct")
      return(NA)
    }
    setResponseSpct(ans, time.unit = time.unit[1])
    if (photon.based.input) {
      e2q(ans, action = "add", byref = TRUE)
    }
  } else if (l.class == "chroma.spct") {
    setChromSpct(ans)
  } else if (l.class == "generic.spct") {
    setGenericSpct(ans)
  }
  if (any(rescaled.input)) {
    setattr(ans, "rescaled", "TRUE")
  }
  if (any(normalized.input)) {
    setattr(ans, "normalized", "TRUE")
  }
  if (add.factor && !add.bswf) {
    keys <- c(factor.name, "w.length")
    setkeyv(ans, keys)
    if (!is.null(comment.ans)) setattr(ans, "comment", comment.ans)
  } else if (!add.factor && add.bswf) {
    keys <- c("BSWF", "w.length")
    setkeyv(ans, keys)
    if (!is.null(comment.ans)) setattr(ans, "comment", comment.ans)
  } else if (add.factor && add.bswf) {
    keys <- c("BSWF", factor.name, "w.length")
    setkeyv(ans, keys)
    if (!is.null(comment.ans)) setattr(ans, "comment", comment.ans)
  } # else we keep the default "w.length"
  return(ans)
}


# subset ------------------------------------------------------------------


#' Subsetting generic.spct
#'
#' Returns subsets of a generic.spct
#'
#' @usage subset(x, subset, select, ...)
#'
#' @param x	generic.spct to subset
#' @param subset logical expression indicating elements or rows to keep
#' @param select expression indicating columns to select from x (IGNORED)
#' @param ...	further arguments to be passed to or from other methods
#'
#' @details The subset argument works on the rows and will be evaluated in the
#' generic.spct so columns can be referred to (by name) as variables in the expression
#' The generic.spct that is returned will maintain the original attributes and keys as long as they are not select-ed out.
#'
#' @return A generic.spct containing the subset of rows and columns that are selected.
#'
#' @method subset generic.spct
#'
#' @seealso \code{\link{subset}}
#'

subset.generic.spct <- function(x, subset, select, ...) {
  comment <- comment(x)
  x.out <- x[eval(substitute(subset))]
  if (!is.null(comment)) {
    setattr(x.out, "comment", comment)
  }
  return(x.out)
}

#' Subsetting source.spct
#'
#' Returns subsets of a source.spct
#'
#' @usage subset(x, subset, select, ...)
#'
#' @param x	source.spct to subset
#' @param subset logical expression indicating elements or rows to keep
#' @param select expression indicating columns to select from x (IGNORED)
#' @param ...	further arguments to be passed to or from other methods
#'
#' @details The subset argument works on the rows and will be evaluated in the
#' source.spct so columns can be referred to (by name) as variables in the expression
#' The source.spct that is returned will maintain the original attributes and keys as long as they are not select-ed out.
#'
#' @return A source.spct containing the subset of rows and columns that are selected.
#'
#' @export
#'
#' @examples
#'
#' subset(sun.spct, w.length > 400)
#'
#' @seealso \code{\link{subset}}
#'

subset.source.spct <- function(x, subset, select = NULL, ...) {
  time.unit <- getTimeUnit(x)
  bswf.used <- getBSWFUsed(x)
  comment <- comment(x)
  x.out <- x[eval(substitute(subset))]
  setSourceSpct(x = x.out, time.unit = time.unit, bswf.used = bswf.used)
  if (!is.null(comment)) {
    setattr(x.out, "comment", comment)
  }
  return(x.out)
}

#' Subsetting filter.spct
#'
#' Returns subsets of a filter.spct
#'
#' @usage subset(x, subset, select, ...)
#'
#' @param x	filter.spct to subset
#' @param subset logical expression indicating elements or rows to keep
#' @param select expression indicating columns to select from x (IGNORED)
#' @param ...	further arguments to be passed to or from other methods
#'
#' @details The subset argument works on the rows and will be evaluated in the
#' filter.spct so columns can be referred to (by name) as variables in the expression
#' The filter.spct that is returned will maintain the original attributes and keys as long as they are not select-ed out.
#'
#' @return A filter.spct containing the subset of rows and columns that are selected.
#'
#' @export
#'
#' @seealso \code{\link{subset}}
#'

subset.filter.spct <- function(x, subset, select = NULL, ...) {
  Tfr.type <- getTfrType(x)
  comment <- comment(x)
  x.out <- x[eval(substitute(subset))]
  setFilterSpct(x = x.out, Tfr.type = Tfr.type)
  if (!is.null(comment)) {
    setattr(x.out, "comment", comment)
  }
  return(x.out)
}

#' Subsetting reflector.spct
#'
#' Returns subsets of a reflector.spct
#'
#' @usage subset(x, subset, select, ...)
#'
#' @param x	reflector.spct to subset
#' @param subset logical expression indicating elements or rows to keep
#' @param select expression indicating columns to select from x (IGNORED)
#' @param ...	further arguments to be passed to or from other methods
#'
#' @details The subset argument works on the rows and will be evaluated in the
#' reflector.spct so columns can be referred to (by name) as variables in the expression
#' The reflector.spct that is returned will maintain the original attributes and keys as long as they are not select-ed out.
#'
#' @return A reflector.spct containing the subset of rows and columns that are selected.
#'
#' @export
#'
#' @seealso \code{\link{subset}}
#'

subset.reflector.spct <- function(x, subset, select = NULL, ...) {
  Rfr.type <- getRfrType(x)
  comment <- comment(x)
  x.out <- x[eval(substitute(subset))]
  setReflectorSpct(x = x.out, Rfr.type = Rfr.type)
  if (!is.null(comment)) {
    setattr(x.out, "comment", comment)
  }
  return(x.out)
}

#' Subsetting response.spct
#'
#' Returns subsets of a response.spct
#'
#' @usage subset(x, subset, select, ...)
#'
#' @param x	response.spct to subset
#' @param subset logical expression indicating elements or rows to keep
#' @param select expression indicating columns to select from x (IGNORED)
#' @param ...	further arguments to be passed to or from other methods
#'
#' @details The subset argument works on the rows and will be evaluated in the
#' response.spct so columns can be referred to (by name) as variables in the expression
#' The response.spct that is returned will maintain the original attributes and keys as long as they are not select-ed out.
#'
#' @return A response.spct containing the subset of rows and columns that are selected.
#'
#' @export
#'
#' @seealso \code{\link{subset}}
#'

subset.response.spct <- function(x, subset, select = NULL, ...) {
  time.unit <- getTimeUnit(x)
  comment <- comment(x)
  x.out <- x[eval(substitute(subset))]
  setResponseSpct(x = x.out, time.unit = time.unit)
  if (!is.null(comment)) {
    setattr(x.out, "comment", comment)
  }
  return(x.out)
}

#' Subsetting object.spct
#'
#' Returns subsets of a object.spct
#'
#' @usage subset(x, subset, select, ...)
#'
#' @param x	object.spct to subset
#' @param subset logical expression indicating elements or rows to keep
#' @param select expression indicating columns to select from x (IGNORED)
#' @param ...	further arguments to be passed to or from other methods
#'
#' @details The subset argument works on the rows and will be evaluated in the
#' object.spct so columns can be referred to (by name) as variables in the expression
#' The object.spct that is returned will maintain the original attributes and keys as long as they are not select-ed out.
#'
#' @return A object.spct containing the subset of rows and columns that are selected.
#'
#' @export
#'
#' @seealso \code{\link{subset}}
#'

subset.object.spct <- function(x, subset, select = NULL, ...) {
  Tfr.type <- getTfrType(x)
  Rfr.type <- getRfrType(x)
  comment <- comment(x)
  x.out <- x[eval(substitute(subset))]
  setObjectSpct(x = x.out, Tfr.type = Tfr.type, Rfr.type = Rfr.type)
  if (!is.null(comment)) {
    setattr(x.out, "comment", comment)
  }
  return(x.out)
}
