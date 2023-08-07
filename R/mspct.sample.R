#' Random sample of spectra
#'
#' A method to extract a random sample of spectra from a collection of spectra
#' or a spectrum object containing multiple spectra in long form.
#'
#' @param x generic_spct or generic_mspct An object possibly containing multiple
#'   spectra.
#' @param size integer The number of spectra to extract, if available.
#' @param replace logical Sample with or without replacement.
#' @param keep.order logical Return the spectra ordered as in \code{x} or in
#' random order.
#'
#' @return If \code{x} is an spectrum object, such as a
#'   \code{"filter_spct"} object, the returned object is of the same class but
#'   in most cases containing fewer spectra in long form than \code{x}.
#'    If \code{x} is a collection of spectrum objecta, such as a
#'   \code{"filter_mspct"} object, the returned object is of the same class but
#'   in most cases containing fewer member spectra than \code{x}.
#'
#' @seealso See \code{\link[base]{sample}} for the method used for
#'   the sampling.
#'
#' @export
#'
sample_spct <- function(x,
                        size = 1,
                        replace = FALSE,
                        recursive = FALSE,
                        keep.order = TRUE) {
  stopifnot(inherits(x, what = "generic_spct"))
  num.spectra <- getMultipleWl(x)
  if (num.spectra <= size) {
    # nothing to do
    return(x)
  }
  # brute force method ensures handling of metadata attributes in subset2mspct()
  mspct <- subset2mspct(x)
  z <- sample_mspct(x = mspct,
                    size = size,
                    replace = replace,
                    recursive = recursive,
                    keep.order = keep.order)
  rbindspct(z)
}

#' @rdname sample_spct
#'
#' @param recursive logical Sample members of x directly, or first split
#'    members containing multiple spectra.
#' @param simplify logical When \code{n = 1} return a collection with a single
#'    member or a spectrum object.
#'
#' @export
#'
sample_mspct <- function(x,
                         size = 1,
                         replace = FALSE,
                         recursive = FALSE,
                         keep.order = TRUE,
                         simplify = FALSE) {
  stopifnot(inherits(x, what = "generic_mspct"))
  if (length(x) <= size) {
    # nothing to do
    return(x)
  }
  if (recursive) {
    # separate multiple spectra within individual members
    x <- subset2mspct(x)
  }
  selector.idx <- sample(x = length(x), size = size, replace = replace)
  if (keep.order) {
    selector.idx <- sort(selector.idx)
  }
  selector.names <- names(x)[selector.idx]
  if (simplify && size == 1) {
    x[[selector.names]]
  } else {
    x[selector.names]
  }
}
