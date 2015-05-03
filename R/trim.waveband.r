#' Trim tails of the wavebands based on wavelength limits.
#'
#' Trimming is needed for example when the spectral data does not cover the
#' whole waveband.
#'
#' @param w.band an object of class "waveband" or a list of such objects
#' @param range a numeric vector of length two, or any other object for which
#'   function range() will return two
#' @param low.limit shortest wavelength to be kept (defaults to shortest
#'   w.length value)
#' @param high.limit longest wavelength to be kept (defaults to longest w.length
#'   value)
#' @param trim logical (default is FALSE, which just discards off range
#'   wavebands)
#' @param use.hinges logical, if TRUE (the default)
#'
#' @return a waveband objects or a list of waveband objects trimmed or filtered
#'
#' @keywords manip misc
#' @family trim functions
#' @export
#' @examples
#' trim_waveband(waveband(c(200,1000)), c(400,700))
#'
trim_waveband <- function(w.band,
                          range = NULL,
                          low.limit = NULL, high.limit = NULL,
                          trim = getOption("photobiology.waveband.trim", default = TRUE),
                          use.hinges = TRUE)
{
  if (!is.null(w.band) && is.waveband(w.band)) {
    w.band <- list(w.band)
  }
  verbose <- TRUE
  if (!is.null(range)) {
    low.limit <- ifelse(!is.null(low.limit), max(min(range), low.limit), min(range))
    high.limit <- ifelse(!is.null(high.limit), min(max(range), high.limit), max(range))
  }
  if (is.null(low.limit)) {
    low.limit <- min(w.band)
  }
  if (is.null(high.limit)) {
    high.limit <- max(w.band)
  }
  w.band.out <- list()
  i <- 0
  for (wb in w.band) {
    if (min(wb) >= high.limit || max(wb) <= low.limit) {
      next
    }
    if (min(wb) >= low.limit && max(wb) <= high.limit) {
      i <- i + 1L
      w.band.out[i] <- list(wb)
    } else if (trim) {
      trimmed.wb <- wb
      trimmed.high <- trimmed.low <- FALSE
      if (min(wb) < low.limit) {
        trimmed.wb$low <- low.limit
        trimmed.wb$hinges <- unique(sort(c(low.limit - 1e-4, low.limit, wb$hinges[wb$hinges>=low.limit])))
        trimmed.low <- TRUE
      }
      if (max(wb) > high.limit) {
        trimmed.wb$high <- high.limit
        trimmed.wb$hinges <- unique(sort(c(wb$hinges[wb$hinges<=high.limit], high.limit - 1e-4, high.limit)))
        trimmed.high <- TRUE
      }
      if (trimmed.low || trimmed.high) {
        trimmed.tag <-  paste("tr", ifelse(trimmed.low, ".lo", ""),
                              ifelse(trimmed.high, ".hi", ""), sep="")
        trimmed.wb$label <- paste(wb$label, trimmed.tag, sep=".")
        trimmed.wb$name <- paste(wb$name, trimmed.tag, sep=".")
        i <- i + 1L
        w.band.out[[i]] <- trimmed.wb
      }
    }
  }
  if (length(w.band.out) < 1L) {
    return(NULL)
  } else {
    return(w.band.out)
  }
}
