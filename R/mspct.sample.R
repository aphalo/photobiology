#' Random sample from collection of spectra
#'
#' A method to extract a random sample of spectra from a collection of spectra
#' or a spectrum object containing multiple spectra in long form.
#'
#'
#'
#' @param x generic_spct An object possibly containing multiple spectra.
#' @param size integer The number of spectra to extract, if available.
#'
#' @return If \code{x} is an object, such as a
#'   "filter_spct" object, the returned object is of same class but possibly
#'   containing fewer spectra in long form than \code{x}.
#'
#' @seealso See \code{\link[base]{sample}} for the method used for
#'   the sampling.
#'
#' @export
#'
sample_spct <- function(x, size = 1) {
  stopifnot(inherits(x, what = "generic_spct"))
  num.spectra <- getMultipleWl(x)
  if (num.spectra <= size) {
    # nothing to do
    return(x)
  }
  # brute force method allows handling of metadata attributes in subset2mspct()
  mspct <- subset2mspct(x)
  z <- sample_mspct(x = mspct, size = size, recursive = FALSE)
  rbindspct(z)
}

#' @rdname sample_spct
#'
#' @param recursive logical Sample members of x directly, or first split
#'    members containing multiple spectra.
#'
#' @export
#'
sample_mspct <- function(x, size = 1, recursive = FALSE) {
  stopifnot(inherits(x, what = "generic_mspct"))
  if (length(x) <= size) {
    # nothing to do
    return(x)
  }
  if (recursive) {
    # separate multiple spectra within individual members
    x <- subset2mspct(x)
  }
  selector.idx <- sample(x = length(x), size = size)
  selector.names <- names(x)[selector.idx]
  x[selector.names]
}
