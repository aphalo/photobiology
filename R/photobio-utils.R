#' Compute range and format it
#'
#' Compute the range of an R object, and format it as string suitable for
#' printing.
#'
#' @param x an R object
#' @param na.rm logical, indicating if \code{\link[base]{NA}}'s should be
#'   omitted.
#' @param digits,nsmall numeric, passed to same name parameters of
#'   \code{format()}.
#' @param collapse character, passed to same name parameter of \code{paste()}.
#'
#' @seealso \code{\link[base]{range}}, \code{\link[base]{format}} and
#'   \code{\link[base]{paste}}.
#'
#' @export
#'
#' @examples
#' formatted_range(c(1, 3.5, -0.01))
#'
formatted_range <- function(x,
                            na.rm = TRUE,
                            digits = 3,
                            nsmall = 2,
                            collapse = "..") {
  paste(format(range(x, na.rm = na.rm),
               digits = digits,
               nsmall = nsmall,
               trim = TRUE),
        collapse = collapse)
}
