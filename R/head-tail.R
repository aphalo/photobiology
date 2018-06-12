#' Return the First and Last Part of an Object
#'
#' Returns the first and last "parts" (rows or members) of a spectrum,
#' dataframe, vector, function, table or ftable. In other words, the combined
#' output from methods \code{head} and \code{tail}.
#'
#' @param x an R object.
#' @param n a single integer. If positive, half the size for the resulting
#'   object: number of elements for a vector (including lists), rows for a
#'   matrix or data frame or lines for a function, to be retained from each
#'   end. If negative, all but the n last and n first number of elements of x.
#' @param ... arguments to be passed to or from other methods.
#'
#' @return An object (usually) like x but generally smaller. For ftable objects
#'  \code{x}, a transformed \code{format(x)}.
#'
#' @details The value returned by head_tail() is equivalent to row binding the
#'  the values returned by head() and tail(), although not implemented in this
#'  way. The same specializations as defined in package 'utils' for head() and
#'  tail() have been implemented.
#'
#' @seealso \code{\link[utils]{head}}, and compare the examples and the values
#'  returned to the examples below.
#'
#' @note For some types of input, like functions, the output may be confusing,
#'  however, we have opted for consistency with existing functions. The code
#'  is in part a revision of that of head() and tail() from package 'utils'. I
#'  have been missing this method especially when checking spectral data, as
#'  both ends are of interest.
#'
#' @export
#'
#' @examples
#'
#' head_tail(letters)
#' head_tail(letters, n = -6L)
#' head_tail(freeny.x, n = 10L)
#' head_tail(freeny.y)
#'
#' head_tail(stats::ftable(Titanic))
#'
#' @keywords
#'
head_tail <- function(x, n, ...) UseMethod("head_tail")

#' @describeIn head_tail
#'
#' @export
#'
head_tail.default <- function(x, n = 3L, ...) {
  stopifnot(length(n) == 1L)
  if (n > 0) {
    if ((2 * n) < length(x)) {
      selector <- unique(c(1:n, (length(x) - n + 1):length(x)))
    } else {
      selector <- TRUE
    }
  } else if (n <= 0) {
    if ((2 * -n) < length(x)) {
      selector <- (1 - n):(length(x) + n)
    } else {
      selector <- FALSE
    }
  }
  x[selector]
}

#' @describeIn head_tail
#'
#' @export
#'
head_tail.data.frame <- function(x, n = 3L, ...) {
  stopifnot(length(n) == 1L)
  if (n > 0) {
    if ((2 * n) < nrow(x)) {
      selector <- unique(c(1:n, (nrow(x) - n + 1):nrow(x)))
    } else {
      selector <- TRUE
    }
  } else if (n <= 0) {
    if ((2 * -n) < nrow(x)) {
      selector <- (1 - n):(nrow(x) + n)
    } else {
      selector <- FALSE
    }
  }
  x[selector, , drop = FALSE]
}

#' @describeIn head_tail
#'
#' @export
#'
head_tail.matrix <- head_tail.data.frame

#' @describeIn head_tail
#'
#' @export
#'
#' @note \code{head_tail()} methods for function, table and ftable classes, are
#'   wrappers for head() method.
#'
head_tail.function <- function (x, n = 6L, ...) {
  lines <- as.matrix(deparse(x))
  dimnames(lines) <- list(seq_along(lines), "")
  noquote(head_tail(lines, n = n))
}

#' @describeIn head_tail
#'
#' @export
#'
head_tail.table <- function (x, n = 6L, ...) {
  (if (length(dim(x)) == 2L)
    head_tail.matrix
   else head_tail.default)(x, n = n)
}

#' @describeIn head_tail
#'
#' @export
#'
head_tail.ftable <-  function (x, n = 6L, ...) {
  r <- format(x)
  dimnames(r) <- list(rep.int("", nrow(r)), rep.int("", ncol(r)))
  noquote(head_tail.matrix(r, n = n + nrow(r) - nrow(x), ...))
}
