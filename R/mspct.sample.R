#' Random sample of spectra
#'
#' A method to extract a random sample of members from a list, a collection of
#' spectra or a spectrum object containing multiple spectra in long form.
#'
#' @param x An R object possibly containing multiple spectra or other
#'   components.
#' @param size integer The number of spectra to extract, if available.
#' @param replace logical Sample with or without replacement.
#' @param recursive logical If code{x} is a collection, expand or not member
#'   spectra containing multiple spectra in long form into individual members
#'   before sampling.
#' @param keep.order logical Return the spectra ordered as in \code{x} or in
#'   random order.
#' @param simplify logical If \code{size = 1}, and code{x} is a collection
#'   return the spectrum object instead of a colelction with it as only member.
#' @param ... currently ignored.
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
#' @examples
#' a.list <- as.list(letters)
#' names(a.list) <- LETTERS
#' set.seed(12345678)
#' pull_sample(a.list, size = 8)
#' pull_sample(a.list, size = 8, keep.order = FALSE)
#' pull_sample(a.list, size = 8, replace = TRUE)
#' pull_sample(a.list, size = 8, replace = TRUE, keep.order = FALSE)
#' pull_sample(a.list, size = 1)
#' pull_sample(a.list, size = 1, simplify = TRUE)
#'
#' @export
#'
pull_sample <- function(x, size, ...) {
  UseMethod("pull_sample")
}

#' @describeIn pull_sample Default for generic function
#'
#' @export
#'
pull_sample.default <- function(x, size, ...) {
  warning("'pull_sample' is not defined for objects of class ", class(x)[1])
  return(generic_mspct())
}

#' @describeIn pull_sample Specialization for generic_spct
#'
#' @export
#'
pull_sample.list <- function(x,
                             size = 1,
                             replace = FALSE,
                             keep.order = TRUE,
                             simplify = FALSE,
                             ...) {
  if (length(x) <= size) {
    # nothing to do
    return(x)
  }
  selector.idx <- sample(x = length(x), size = size, replace = replace)
  if (keep.order) {
    selector.idx <- sort(selector.idx)
  }
  if (simplify && size == 1) {
    z <- x[[selector.idx]]
  } else {
    z <- x[selector.idx]
    if (replace && length(names(x))) {
      names(z) <- make.unique(names(x)[selector.idx], sep = ".copy")
    }
  }
  z
}

#' @describeIn pull_sample Specialization for generic_spct
#'
#' @export
#'
pull_sample.generic_spct <- function(x,
                                     size = 1,
                                     replace = FALSE,
                                     recursive = FALSE,
                                     keep.order = TRUE,
                                     ...) {
  num.spectra <- getMultipleWl(x)
  if (num.spectra <= size) {
    # nothing to do
    return(x)
  }
  # brute force method ensures handling of metadata attributes
  mspct <- subset2mspct(x)
  z <- pull_sample(x = mspct,
                   size = size,
                   replace = replace,
                   recursive = recursive,
                   keep.order = keep.order)
  rbindspct(z)
}

#' @describeIn pull_sample Specialization for generic_mspct
#'
#' @export
#'
pull_sample.generic_mspct <- function(x,
                                     size = 1,
                                     replace = FALSE,
                                     recursive = FALSE,
                                     keep.order = TRUE,
                                     simplify = FALSE,
                                     ...) {
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
  if (simplify && size == 1) {
    z <- x[[selector.idx]]
  } else {
    z <- x[selector.idx]
    if (replace && length(names(x))) {
      names(z) <- make.unique(names(x)[selector.idx], sep = ".copy")
    }
  }
  z
}
