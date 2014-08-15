#' Makes one spectral object from a list of many
#'
#' Same as \code{rbindlist} from package data.table but preserves class of spectral objects. Has different defaults
#' for use names and fill.
#'
#' @usage rbindspct(l, use.names=fill, fill=TRUE)
#'
#' @usage rbindlist(l, use.names=fill, fill=FALSE)
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
#' \code{use.names} has to be \code{TRUE}, and all items of the input list has to have non-null column names.
#'
#' @details
#' Each item of \code{l} can be a spectrum, \code{data.table}, \code{data.frame} or \code{list}, including \code{NULL} (skipped)
#' or an empty object (0 rows). \code{rbindspc} is most useful when there are a variable number of (potentially many)
#' objects to stack. \code{rbind} (not implemented yet for spectra) however is most useful to
#' stack two or three objects which you know in advance. \code{\dots} should contain at least one \code{data.table}
#' for \code{rbind(...)} to call the fast method and return a \code{data.table}, whereas \code{rbindlist(l)} always
#' returns a \code{data.table}, and \code{rbindspct} always returns at least a \code{generic.spct}, even when stacking
#' a plain \code{list} with a \code{data.frame}, for example.
#  With these changes, the only difference between \code{rbindspct(l)} and \code{rbindlist(l)} is their
#' \emph{default argument} \code{use.names}.
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
#' @export rbindspct rbindlist
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
#' rbindspct(list(sun.spct[1:100], sun.spct[300:400]))
#'
#' rbindlist(list(sun.spct[1:100], sun.spct[300:400]))
#'
#' # examples from package data.table
#'
#' # default case
#' DT1 = data.table(A=1:3,B=letters[1:3])
#' DT2 = data.table(A=4:5,B=letters[4:5])
#' l = list(DT1,DT2)
#' rbindlist(l)
#'
#' # bind correctly by names
#' DT1 = data.table(A=1:3,B=letters[1:3])
#' DT2 = data.table(B=letters[4:5],A=4:5)
#' l = list(DT1,DT2)
#' rbindlist(l, use.names=TRUE)
#'
#' # fill missing columns, and match by col names
#' DT1 = data.table(A=1:3,B=letters[1:3])
#' DT2 = data.table(B=letters[4:5],C=factor(1:2))
#' l = list(DT1,DT2)
#' rbindlist(l, use.names=TRUE, fill=TRUE)
#'
#'
#' @keywords data
#'
#' @aliases rbindspc rbindlist
#'

rbindspct <- function(l, use.names = fill, fill = TRUE) {
  # rbindlist strips attributes and sets class to data.table
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
    setSourceSpct(ans)
  } else if (l.class == "filter.spct") {
    setFilterSpct(ans)
  } else if (l.class == "reflector.spct") {
    setReflectorSpct(ans)
  } else if (l.class == "response.spct") {
    setResponseSpct(ans)
  } else if (l.class == "chroma.spct") {
    setChromSpct(ans)
  } else if (l.class == "generic.spct") {
    setGenericSpct(ans)
  }
  return(ans)
}

rbindlist <- function(l, use.names = fill, fill = FALSE) {
  rbindspct(l, use.names, fill)
}
