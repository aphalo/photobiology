#' Defunct functions and methods
#'
#' Functions listed here have been removed or deleted, and temporarily replaced
#' by stubs that report this when they are called.
#'
#' @name defunct
#'
#' @param ... ignored
#'
#' @export
#'
#' @note Function \code{f_mscpt} has been renamed \code{msdply()}.
#'
f_mspct <- function(...) {
  stop("Function 'f_mspct()' is defunct, please use 'msdply()' instead.")
}

#' @rdname defunct
#'
#' @export
#'
#' @note Function \code{mutate_mscpt} has been renamed \code{msmsply()}.
#'
mutate_mspct <- function(...) {
  stop("Function 'mutate_mspct()' is defunct, please use 'msdply()' instead.")
}
