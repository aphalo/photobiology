#' Handle Missing Values in Objects
#'
#' These methods are useful for dealing with NAs in e.g., \code{source_spct},
#' \code{filter_spct} and \code{reflector_spct}
#'
#' @param object an R object
#' @param ... further arguments other special methods could require
#'
#' @details If na.omit removes cases, the row numbers of the cases form the
#' "na.action" attribute of the result, of class "omit".
#'
#' @export
#' @importFrom stats na.omit
#'
#' @examples
#' my_sun.spct <- sun.spct
#' my_sun.spct[3, "s.e.irrad"] <- NA
#' my_sun.spct[5, "s.q.irrad"] <- NA
#' na.omit(my_sun.spct)
#' na.action(na.omit(my_sun.spct))
#'
na.omit.source_spct <- function(object, ...) {
  data_cols <- which(names(object) %in% c("s.e.irrad", "s.q.irrad"))
  rows_to_omit <- integer()
  for (col in data_cols) {
    rows_to_omit <- union(rows_to_omit, which(is.na(object[[col]])))
  }
  rows_to_keep <- setdiff(1:nrow(object), rows_to_omit)
  z <- dplyr::slice(object, rows_to_keep)
  z <- as.source_spct(z)
  z <- copy_attributes(object, z)
  class(rows_to_omit) <- "omit"
  attr(z, "na.action") <- rows_to_omit
  z
}

#' @rdname na.omit.source_spct
#'
#' @export
#'
na.omit.response_spct <- function(object, ...) {
  data_cols <- which(names(object) %in% c("s.e.response", "s.q.response"))
  rows_to_omit <- integer()
  for (col in data_cols) {
    rows_to_omit <- union(rows_to_omit, which(is.na(object[[col]])))
  }
  rows_to_keep <- setdiff(1:nrow(object), rows_to_omit)
  z <- dplyr::slice(object, rows_to_keep)
  z <- as.response_spct(z)
  z <- copy_attributes(object, z)
  class(rows_to_omit) <- "omit"
  attr(z, "na.action") <- rows_to_omit
  z
}

#' @rdname na.omit.source_spct
#'
#' @export
#'
na.omit.filter_spct <- function(object, ...) {
  data_cols <- which(names(object) %in% c("Tfr", "A"))
  rows_to_omit <- integer()
  for (col in data_cols) {
    rows_to_omit <- union(rows_to_omit, which(is.na(object[[col]])))
  }
  rows_to_keep <- setdiff(1:nrow(object), rows_to_omit)
  z <- dplyr::slice(object, rows_to_keep)
  z <- as.filter_spct(z)
  z <- copy_attributes(object, z)
  class(rows_to_omit) <- "omit"
  attr(z, "na.action") <- rows_to_omit
  z
}

#' @rdname na.omit.source_spct
#'
#' @export
#'
na.omit.reflector_spct <- function(object, ...) {
  data_col <- which(names(object) == "Rfr")
  rows_to_omit <- which(is.na(object[[data_col]]))
  rows_to_keep <- setdiff(1:nrow(object), rows_to_omit)
  z <- dplyr::slice(object, rows_to_keep)
  z <- as.reflector_spct(z)
  z <- copy_attributes(object, z)
  class(rows_to_omit) <- "omit"
  attr(z, "na.action") <- rows_to_omit
  z
}

