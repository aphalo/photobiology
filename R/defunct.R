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
#' @note Function \code{f_mspct()} has been renamed \code{msdply()}.
#'
f_mspct <- function(...) {
  stop("Function 'photobiology::f_mspct()' is defunct, please use 'msdply()' instead.")
}

#' @rdname defunct
#'
#' @export
#'
#' @note Function \code{mutate_mspct()} has been renamed \code{msmsply()}.
#'
mutate_mspct <- function(...) {
  stop("Function 'photobiology::mutate_mspct()' is defunct, please use 'msdply()' instead.")
}

#' @rdname defunct
#'
#' @export
#'
#' @note Function \code{calc_filter_multipliers()} has been removed.
#'
calc_filter_multipliers <- function(...) {
  stop("Function 'photobiology::calc_filter_multipliers()' is defunct, please use
'interpolate_wl()' instead.")
}

#' @rdname defunct
#'
#' @export
#'
#' @note Function \code{calc_filter_multipliers()} has been removed.
#'
T2T <- function(...) {
  stop("Function 'photobiology::T2T()' is defunct, please use
'convertTfrType()' instead.")
}

#' @rdname defunct
#'
#' @export
#'
#' @note Method \code{getAfrType()} has been removed.
#'
getAfrType <- function(...) {
  stop("Function 'photobiology::getAfrType()' is defunct, the \"Afr.type\"
attribute is redundant and no longer used.")
}

#' @rdname defunct
#'
#' @export
#'
#' @note Method \code{setAfrType()} has been removed.
#'
setAfrType <- function(...) {
  stop("Function 'photobiology::setAfrType()' is defunct, the \"Afr.type\"
attribute is redundant and no longer used.")
}

