#' Upgrade spectral objects
#'
#' Update the spectral class names of objects to those used in photobiology (>=
#' 0.6.0).
#'
#' @param object generic.spct
#' @param ... not used
#'
#' @note The object is modified by reference. The class names with ending
#'   ".spct" replaced by their new equivalents ending in "_spct".
#'
#' @return The modified object (invisibly).
#'
#' @export upgrade.generic.spct
#' @exportClass generic.spct
#'
upgrade.generic.spct <-
  function(object, ...) {
    name <- substitute(object)
    class(object) <- gsub(".spct", "_spct", class(object), fixed = TRUE)
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, object, parent.frame(), inherits = TRUE)
    }
    invisible(object)
  }

