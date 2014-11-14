#' Makes one spectral object from a list of many
#'
#' Same as \code{rbindlist} from package data.table but preserves class of spectral objects. Has different defaults
#' for use names and fill.
#'
#' @usage rbindspct(l, use.names=TRUE, fill=TRUE)
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
#' @details
#' Each item of \code{l} can be a spectrum, \code{data.table}, \code{data.frame} or \code{list}, including \code{NULL} (skipped)
#' or an empty object (0 rows). \code{rbindspc} is most useful when there are a variable number of (potentially many)
#' objects to stack. \code{rbind} (not implemented yet for spectra) however is most useful to
#' stack two or three objects which you know in advance. \code{\dots} should contain at least one \code{data.table}
#' for \code{rbind(...)} to call the fast method and return a \code{data.table}, whereas \code{rbindspct} always
#' returns at least a \code{generic.spct} as long as all elements in in are spectra, atherwise a \code{data.frame}
#' is returned even when stacking
#' a plain \code{list} with a \code{data.frame}, for example.
#  The only difference between \code{rbindspct(l)} and \code{rbindlist(l)} from package \code{data.frame} is in their
#' \emph{default arguments} \code{use.names}, and in that \code{rbindlist} will always return a \code{data.frame} even
#' when the list l contains only spectra.
#' If column \code{i} of input items do not all have the same type; e.g, a \code{data.table} may be bound with a
#' \code{list} or a column is \code{factor} while others are \code{character} types, they are coerced to the highest
#' type (SEXPTYPE).
#'
#' Note that any additional attributes that might exist on individual items of the input list would not be preserved
#' in the result.
#'
#' @return An spectral object of a type common to all bound items or a \code{data.table} containing a concatenation of
#' all the items passed in.
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
#' class(rbindspct(list(sun.spct[1:100], sun.spct[300:400])))
#'
#' @keywords data
#'

rbindspct <- function(l, use.names = TRUE, fill = TRUE) {
  # original rbindlist from data.table strips attributes and sets class to data.table
  l.class <- c( "source.spct", "filter.spct", "reflector.spct", "response.spct", "chroma.spct",
                "generic.spct", "data.table", "data.frame")
  for (spct in l) {
    l.class <- intersect(l.class, class(spct))
  }
  l.class <- l.class[1]
#  print(l.class)
  ans <- data.table::rbindlist(l, use.names, fill)
  if (is.null(ans)) {
    NULL
  } else if (l.class == "source.spct") {
    time.unit <- character(length(l))
    i <- 0L
    for (scpt in l) {
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
    for (scpt in l) {
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
  invisible(ans)
}
