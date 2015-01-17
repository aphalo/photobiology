#' Makes one spectral object from a list of many
#'
#' Same as \code{rbindlist} from package data.table but preserves class of spectral objects. Has different defaults
#' for use names and fill.
#'
#' @usage rbindspct(l, use.names = TRUE, fill=TRUE, add.factor = FALSE)
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
#' @param add.factor logical indicating if a factor should be added to distinguish data from each spectrum
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
#' If column \code{i} of the different input items do not all have the same type; e.g, a \code{generic.spct} may be
#' bound with a \code{list} or a column is \code{factor} while others are \code{character} types, they are coerced
#' to the highest type (SEXPTYPE).
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
#' @note data.table::rbindlist is called internally and the result returned as is unless all elements in the list
#' belong to one of the \code{.spct} classes.
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
#' spct <- rbindspct(list(sun.spct, sun.spct), add.factor = TRUE)
#' head(spct)
#' class(spct)
#'
#' # adds factor 'spct.idx' with the names given to the spectra in the list
#' # supplied as formal argument 'l' as levels
#' spct <- rbindspct(list(one = sun.spct, two = sun.spct), add.factor = TRUE)
#' head(spct)
#' class(spct)
#'
#' @keywords data
#'

rbindspct <- function(l, use.names = TRUE, fill = TRUE, add.factor = FALSE) {
  # original rbindlist from data.table strips attributes and sets class to data.table
  if (is.null(l) || length(l) < 1) {
    return(l)
  }
  if (!is.list(l) || is.any.spct(l) || is.waveband(l)) {
    stop("Argument 'l' should be a list of spectra")
    return(NULL)
  }
  l.class <- c( "source.spct", "filter.spct", "reflector.spct", "response.spct", "chroma.spct",
                "generic.spct")
  for (spct in l) {
    l.class <- intersect(l.class, class(spct))
  }
  if (length(l.class) < 1) {
    warning("Argument 'l' contains objects which are not spectra")
    return(NA)
  }
  l.class <- l.class[1]
  #  print(l.class)
  if (add.factor) {
    names.spct <- names(l)
    if (is.null(names.spct) || anyNA(names.spct)) {
      names.spct <- LETTERS[1:length(l)]
    }
#    comment.spct <- ""
    for (i in 1:length(l)) { # transversing the list with spct
      l[[i]][ , spct.idx := names.spct[i] ]
#      comment.spct <- paste(comment.spct, "\nspct.idx: ", names.spct[i], "\n", comment(l[[i]]))
    }
  }
  if (length(l) < 2) {
    ans <- l[[1]]
  } else {
    ans <- data.table::rbindlist(l, use.names, fill)
  }
  if (is.null(ans)) {
    return(NULL)
  }
  if (l.class == "source.spct") {
    time.unit <- character(length(l))
    i <- 0L
    for (spct in l) {
      i <- i + 1
      time.unit[i] <- ifelse(is.null(attr(spct, "time.unit")), "unknown", attr(spct, "time.unit"))
      if (i > 1) {
        if (!(time.unit[i-1] == time.unit[i])) {
          warning("Inconsistent time units among source spectra in rbindspct")
          return(NA)
        }
      }
    }
    setSourceSpct(ans, time.unit = time.unit[1])
  } else if (l.class == "filter.spct") {
    Tfr.type <- character(length(l))
    i <- 0L
    for (spct in l) {
      i <- i + 1
      Tfr.type[i] <- ifelse(is.null(attr(spct, "Tfr.type")), "unknown", attr(spct, "Tfr.type"))
      if (i > 1) {
        if (!(Tfr.type[i-1] == Tfr.type[i])) {
          warning("Inconsistent 'Tfr.type' among filter spectra in rbindspct")
          return(NA)
        }
      }
    }
    setFilterSpct(ans, Tfr.type = Tfr.type[1])
  } else if (l.class == "reflector.spct") {
    setReflectorSpct(ans)
  } else if (l.class == "response.spct") {
    setResponseSpct(ans)
  } else if (l.class == "chroma.spct") {
    setChromSpct(ans)
  } else if (l.class == "generic.spct") {
    setGenericSpct(ans)
  }
  if (add.factor) {
    ans[ , spct.idx := factor(spct.idx)]
    setkey(ans, spct.idx, w.length)
#    settr(ans, "comment", comment.spct)
  }
  return(ans)
}
