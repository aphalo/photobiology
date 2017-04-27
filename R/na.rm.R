#' Handle Missing Values in Objects
#'
#' These methods are useful for dealing with NAs in e.g., \code{source_spct},
#' \code{response_spct}, \code{filter_spct} and \code{reflector_spct}.
#'
#' @param object an R object
#' @param na.action character One of "omit" or "exclude"
#' @param ... further arguments other special methods could require
#'
#' @details If \code{na.omit} removes cases, the row numbers of the cases form
#'   the \code{"na.action"} attribute of the result, of class \code{"omit"}.
#'
#'   \code{na.exclude} differs from \code{na.omit} only in the class of the
#'   "na.action" attribute of the result, which is \code{"exclude"}.
#'
#' @export
#' @importFrom stats na.omit na.exclude
#'
#' @note \code{na.fail} and \code{na.pass} do not require a specialisation
#'   for spectral objetcs. R's definitions work as expected with no need to
#'   override them.
#'
#' @seealso \code{\link[stats]{na.fail}} and \code{\link[stats]{na.action}}
#'
#' @examples
#' my_sun.spct <- sun.spct
#' my_sun.spct[3, "s.e.irrad"] <- NA
#' my_sun.spct[5, "s.q.irrad"] <- NA
#' na.omit(my_sun.spct)
#' na.action(na.omit(my_sun.spct))
#'
na.omit.source_spct <- function(object, na.action = "omit", ...) {
  stopifnot(na.action %in% c("omit", "exclude"))
  data_cols <- which(names(object) %in% c("s.e.irrad", "s.q.irrad"))
  rows_to_omit <- integer()
  for (col in data_cols) {
    rows_to_omit <- union(rows_to_omit, which(is.na(object[[col]])))
  }
  rows_to_keep <- setdiff(1:nrow(object), rows_to_omit)
  z <- dplyr::slice(object, rows_to_keep)
  z <- as.source_spct(z)
  z <- copy_attributes(object, z)
  class(rows_to_omit) <- na.action
  attr(z, "na.action") <- rows_to_omit
  z
}

#' @rdname na.omit.source_spct
#'
#' @export
#'
na.omit.response_spct <- function(object, na.action = "omit", ...) {
  stopifnot(na.action %in% c("omit", "exclude"))
  data_cols <- which(names(object) %in% c("s.e.response", "s.q.response"))
  rows_to_omit <- integer()
  for (col in data_cols) {
    rows_to_omit <- union(rows_to_omit, which(is.na(object[[col]])))
  }
  rows_to_keep <- setdiff(1:nrow(object), rows_to_omit)
  z <- dplyr::slice(object, rows_to_keep)
  z <- as.response_spct(z)
  z <- copy_attributes(object, z)
  class(rows_to_omit) <- na.action
  attr(z, "na.action") <- rows_to_omit
  z
}

#' @rdname na.omit.source_spct
#'
#' @export
#'
na.omit.filter_spct <- function(object, na.action = "omit", ...) {
  stopifnot(na.action %in% c("omit", "exclude"))
  data_cols <- which(names(object) %in% c("Tfr", "A"))
  rows_to_omit <- integer()
  for (col in data_cols) {
    rows_to_omit <- union(rows_to_omit, which(is.na(object[[col]])))
  }
  rows_to_keep <- setdiff(1:nrow(object), rows_to_omit)
  z <- dplyr::slice(object, rows_to_keep)
  z <- as.filter_spct(z)
  z <- copy_attributes(object, z)
  class(rows_to_omit) <- na.action
  attr(z, "na.action") <- rows_to_omit
  z
}

#' @rdname na.omit.source_spct
#'
#' @export
#'
na.omit.reflector_spct <- function(object, na.action = "omit", ...) {
  stopifnot(na.action %in% c("omit", "exclude"))
  data_col <- which(names(object) == "Rfr")
  rows_to_omit <- which(is.na(object[[data_col]]))
  rows_to_keep <- setdiff(1:nrow(object), rows_to_omit)
  z <- dplyr::slice(object, rows_to_keep)
  z <- as.reflector_spct(z)
  z <- copy_attributes(object, z)
  class(rows_to_omit) <- na.action
  attr(z, "na.action") <- rows_to_omit
  z
}

#' @rdname na.omit.source_spct
#'
#' @export
#'
na.omit.cps_spct <- function(object, na.action = "omit", ...) {
  stopifnot(na.action %in% c("omit", "exclude"))
  data_col <- which(names(object) == "cps")
  if (length(data_col) == 0) {
    warning("No spectral data column 'cps' found.")
    return(object)
  }
  rows_to_omit <- which(is.na(object[[data_col]]))
  rows_to_keep <- setdiff(1:nrow(object), rows_to_omit)
  z <- dplyr::slice(object, rows_to_keep)
  z <- as.cps_spct(z)
  z <- copy_attributes(object, z)
  class(rows_to_omit) <- na.action
  attr(z, "na.action") <- rows_to_omit
  z
}

#' @rdname na.omit.source_spct
#'
#' @export
#'
na.omit.raw_spct <- function(object,
                             na.action = "omit",
                             ...) {
  stopifnot(na.action %in% c("omit", "exclude"))
  data_col <- which(names(object) == "counts")
  rows_to_omit <- which(is.na(object[[data_col]]))
  if (length(data_col) == 0) {
    warning("No spectral data column 'counts' found.")
    return(object)
  }
  rows_to_keep <- setdiff(1:nrow(object), rows_to_omit)
  z <- dplyr::slice(object, rows_to_keep)
  z <- as.raw_spct(z)
  z <- copy_attributes(object, z)
  rows_to_omit <- sort(rows_to_omit)
  class(rows_to_omit) <- na.action
  attr(z, "na.action") <- rows_to_omit
  z
}

#' @rdname na.omit.source_spct
#'
#' @export
#'
na.omit.chroma_spct <- function(object, na.action = "omit", ...) {
  stopifnot(na.action %in% c("omit", "exclude"))
  data_cols <- which(names(object) %in% c("x", "y", "z"))
  rows_to_omit <- integer()
  for (col in data_cols) {
    rows_to_omit <- union(rows_to_omit, which(is.na(object[[col]])))
  }
  rows_to_keep <- setdiff(1:nrow(object), rows_to_omit)
  z <- dplyr::slice(object, rows_to_keep)
  z <- as.chroma_spct(z)
  z <- copy_attributes(object, z)
  class(rows_to_omit) <- na.action
  attr(z, "na.action") <- rows_to_omit
  z
}

#' @rdname na.omit.source_spct
#'
#' @export
#'
na.exclude.generic_spct <- function(object, na.action = "exclude", ...) {
  na.omit(object, na.action = na.action, ...)
}
