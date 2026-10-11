#' Random sample of spectra
#'
#' A method to extract a random sample of members from a list, a collection of
#' spectra or a spectrum object containing multiple spectra in long form.
#'
#' @param x An R object possibly containing multiple spectra or other
#'   components.
#' @param size integer The number of spectra to extract, if available.
#' @param replace logical Sample with or without replacement.
#' @param recursive logical If \code{x} is a collection, expand or not member
#'   spectra containing multiple spectra in long form into individual members
#'   before sampling.
#' @param keep.order logical Return the spectra ordered as in \code{x} or in
#'   random order.
#' @param method character Use \code{"random"} for random sampling or
#'   \code{"equal.steps"} for systematic sampling.
#' @param simplify logical If \code{size = 1}, and \code{x} is a collection
#'   return the spectrum object instead of a collection with it as only member.
#' @param ... currently ignored.
#'
#' @details
#' This function calls \code{\link[base]{sample}()} to generate a random set of
#' indexes, or a regular sequence, and uses it to extract members from lists or
#' collections. Method \code{"equal.steps"} is most useful for time series
#' of spectra. With \code{recursive = TRUE} lists and list-like collections
#' are first flattened, and with \code{recursive = FALSE} sampling is applied
#' to the topmost level only.
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
#' pull_sample(a.list, size = 7, method = "equal.steps")
#' pull_sample(a.list, size = 8, keep.order = FALSE)
#' pull_sample(a.list, size = 8, replace = TRUE)
#' pull_sample(a.list, size = 8, replace = TRUE, keep.order = FALSE)
#' pull_sample(a.list, size = 1)
#' pull_sample(a.list, size = 1, simplify = TRUE)
#'
#' set.seed(12345678)
#' pull_sample(sun_evening.spct, 2)
#' set.seed(12345678)
#' pull_sample(sun_evening.mspct, 2)
#'
#' @export
#'
pull_sample <- function(x, size, ...) {
  UseMethod("pull_sample")
}

#' @rdname pull_sample
#'
#' @export
#'
pull_sample.default <- function(x, size, ...) {
  warning("'pull_sample' is not defined for objects of class ", class(x)[1])
  generic_mspct()
}

#' @rdname pull_sample
#'
#' @export
#'
pull_sample.list <- function(x,
                             size = 1,
                             replace = FALSE,
                             keep.order = TRUE,
                             method = "random",
                             simplify = FALSE,
                             ...) {
  size <- as.integer(size)
  if (length(x) <= size) {
    # nothing to do
    return(x)
  }
  if (method == "random") {
    selector.idx <- sample(x = length(x), size = size, replace = replace)
    if (keep.order) {
      selector.idx <- sort(selector.idx)
    }
  } else if (method == "equal.steps") {
    step <- length(x) %/% size
    selector.idx <- seq(from = 1, by = step, length.out = size)
  } else {
    stop("Bad method: \"", method, "\" instead of \"random\" or \"equal.steps\"")
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

#' @rdname pull_sample
#'
#' @export
#'
pull_sample.generic_spct <- function(x,
                                     size = 1,
                                     replace = FALSE,
                                     keep.order = TRUE,
                                     method = "random",
                                     ...) {
  size <- as.integer(size)
  num.spectra <- getMultipleWl(x)
  if (num.spectra <= size) {
    # nothing to do
    return(x)
  }
  if (method == "random") {
    selector.idx <- sample(x = num.spectra, size = size, replace = replace)
    if (keep.order) {
      selector.idx <- sort(selector.idx)
    }
  } else if (method == "equal.steps") {
    step <- length(x) %/% size
    selector.idx <- seq(from = 1, by = step, length.out = size)
  } else {
    stop("Bad method: \"", method, "\" instead of \"random\" or \"equal.steps\"")
  }
  id.factor <- x[[getIdFactor(x)]]
  pulled.ids <- as.character(unique(id.factor)[selector.idx])

  x[id.factor %in% pulled.ids, ]
}

#' @rdname pull_sample
#'
#' @export
#'
pull_sample.generic_mspct <- function(x,
                                      size = 1,
                                      replace = FALSE,
                                      keep.order = TRUE,
                                      method = "random",
                                      recursive = FALSE,
                                      simplify = FALSE,
                                      ...) {
  if (recursive) {
    # separate multiple spectra within individual members
    x <- subset2mspct(x)
  }
  if (length(x) <= size) {
    # nothing to do
    return(x)
  }
  if (method == "random") {
    selector.idx <- sample(x = length(x), size = size, replace = replace)
    if (keep.order) {
      selector.idx <- sort(selector.idx)
    }
  } else if (method == "equal.steps") {
    step <- length(x) %/% size
    selector.idx <- seq(from = 1, by = step, length.out = size)
  } else {
    stop("Bad method: \"", method, "\" instead of \"random\" or \"equal.steps\"")
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
