#' Handle Missing Values in Objects
#'
#' These methods are useful for dealing with NAs in e.g., \code{source_spct},
#' \code{response_spct}, \code{filter_spct} and \code{reflector_spct}.
#'
#' @param object an R object
#' @param na.action character One of "omit", "exclude" or "replace".
#' @param fill numeric Value used to replace NAs unless NULL, in which case
#'    interpolation is attempted.
#' @param target.colnames character Vector of names for the target columns
#'    to operate upon, if present in \code{object}.
#' @param ... further arguments other special methods could require
#'
#' @details If \code{na.omit} removes cases, the row numbers of the cases form
#'   the \code{"na.action"} attribute of the result, of class \code{"omit"}.
#'
#'   \code{na.exclude} differs from \code{na.omit} only in the class of the
#'   "na.action" attribute of the result, which is \code{"exclude"}.
#'
#' @note \code{na.fail} and \code{na.pass} do not require a specialisation
#'   for spectral objects. R's definitions work as expected with no need to
#'   override them. We do not define a method \code{na.replace}, just pass
#'   \code{"replace"} as argument. The current implementation replaces by
#'   interpolation only individual NAs which are flanked on both sides by
#'   valid data. Runs of multiple NAs con only replaced by a constant value
#'   passed through parameter \code{fill}.
#'
#' @seealso \code{\link[stats]{na.fail}} and \code{\link[stats]{na.action}}
#'
#' @export
#'
#' @importFrom stats na.omit na.exclude
#'
#' @name na.omit
#'
#' @examples
#' my_sun.spct <- sun.spct
#' my_sun.spct[3, "s.e.irrad"] <- NA
#' my_sun.spct[5, "s.q.irrad"] <- NA
#'
#' head(my_sun.spct)
#'
#' # rows omitted
#' zo <- na.omit(my_sun.spct)
#' head(zo)
#' na.action(zo)
#'
#' # rows excluded
#' ze <- na.exclude(my_sun.spct)
#' head(ze)
#' na.action(ze)
#'
#' # data in both rows replaced
#' zr <- na.omit(my_sun.spct, na.action = "replace")
#' head(zr)
#' na.action(zr)
#'
# generic used also as worker
na.omit.generic_spct <- function(object,
                                 na.action = "omit",
                                 fill = NULL,
                                 target.colnames,
                                 ...) {
  stopifnot(na.action %in% c("omit", "exclude", "replace"))
  data_cols <- which(colnames(object) %in% target.colnames)
  if (length(data_cols) == 0) {
    warning("No columns matching :", target.colnames, " found.")
    return(object)
  }
  rows_to_omit <- integer()
  for (col in data_cols) {
    rows_to_omit <- union(rows_to_omit, which(is.na(object[[col]])))
  }
  rows_to_keep <- setdiff(seq_len(nrow(object)), rows_to_omit)
  if (na.action == "replace") {
    z <- object
    if (!is.null(fill)) {
      z[rows_to_omit, data_cols] <- as.numeric(fill[1])
    } else {
      for (col in data_cols) {
        # replace existing NA values with interpolated values
        z[[col]] <- v_replace_hinges(z[["w.length"]], z[[col]], z[["w.length"]][rows_to_omit])
      }
    }
  } else {
    # removes rows with NAs in data
 #   z <- dplyr::slice(.data = object, rows_to_keep)
    z <- object[rows_to_keep, ]
#    z <- copy_attributes(object, z, copy.class = TRUE)
  }
  z <- copy_attributes(object, z)
  class(rows_to_omit) <- na.action
  attr(z, "na.action") <- rows_to_omit
  z
}

# omit

#' @rdname na.omit
#'
#' @export
#'
na.omit.source_spct <- function(object, na.action = "omit", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = c("s.e.irrad", "s.q.irrad"))
}

#' @rdname na.omit
#'
#' @export
#'
na.omit.response_spct <- function(object, na.action = "omit", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = c("s.e.response", "s.q.response"))
 }

#' @rdname na.omit
#'
#' @export
#'
na.omit.filter_spct <- function(object, na.action = "omit", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = c("Tfr", "A"))
}

#' @rdname na.omit
#'
#' @export
#'
na.omit.reflector_spct <- function(object, na.action = "omit", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = "Rfr")
}

#' @rdname na.omit
#'
#' @export
#'
na.omit.object_spct <- function(object, na.action = "omit", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = c("Trf", "Rfr"))
}

#' @rdname na.omit
#'
#' @export
#'
na.omit.solute_spct <- function(object, na.action = "omit", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = c("K.mole", "K.mass"))
}

#' @rdname na.omit
#'
#' @export
#'
na.omit.cps_spct <- function(object, na.action = "omit", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = "cps")
}

#' @rdname na.omit
#'
#' @export
#'
na.omit.raw_spct <- function(object, na.action = "omit", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = "counts")
}

#' @rdname na.omit
#'
#' @export
#'
na.omit.chroma_spct <- function(object, na.action = "omit", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = c("x", "y", "z"))
}

#' @rdname na.omit
#'
#' @export
#'
na.omit.generic_mspct <- function(object, na.action = "omit", fill = NULL, ...) {
  msmsply(object, na.omit, na.action = na.action, fill = fill, ...)
}

# exclude

#' @rdname na.omit
#'
#' @export
#'
na.exclude.generic_spct <- function(object, na.action = "exclude", fill = NULL, target.colnames, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = target.colnames)
}

#' @rdname na.omit
#'
#' @export
#'
na.exclude.source_spct <- function(object, na.action = "exclude", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = c("s.e.irrad", "s.q.irrad"))
}

#' @rdname na.omit
#'
#' @export
#'
na.exclude.response_spct <- function(object, na.action = "exclude", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = c("s.e.response", "s.q.response"))
}

#' @rdname na.omit
#'
#' @export
#'
na.exclude.filter_spct <- function(object, na.action = "exclude", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = c("Tfr", "A"))
}

#' @rdname na.omit
#'
#' @export
#'
na.exclude.reflector_spct <- function(object, na.action = "exclude", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = "Rfr")
}

#' @rdname na.omit
#'
#' @export
#'
na.exclude.object_spct <- function(object, na.action = "exclude", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = c("Trf", "Rfr"))
}

#' @rdname na.omit
#'
#' @export
#'
na.exclude.solute_spct <- function(object, na.action = "exclude", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = c("K.mole", "K.mass"))
}

#' @rdname na.omit
#'
#' @export
#'
na.exclude.cps_spct <- function(object, na.action = "exclude", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = "cps")
}

#' @rdname na.omit
#'
#' @export
#'
na.exclude.raw_spct <- function(object, na.action = "exclude", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = "counts")
}

#' @rdname na.omit
#'
#' @export
#'
na.exclude.chroma_spct <- function(object, na.action = "exclude", fill = NULL, ...) {
  na.omit.generic_spct(object = object,
                       na.action = na.action,
                       fill = fill,
                       target.colnames = c("x", "y", "z"))
}

#' @rdname na.omit
#'
#' @export
#'
na.exclude.generic_mspct <- function(object, na.action = "exclude", fill = NULL, ...) {
  msmsply(object, na.exclude, na.action = na.action, fill = fill, ...)
}

