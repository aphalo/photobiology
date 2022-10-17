#' Convert spectrum from wide to long form
#'
#' @details Only objects of classes raw_spct, cps_spct, and object_spct normally contain
#' multiple columns of spectral data. These are supported as well as
#' generic_spct. Is the wide spectra contain multiple spectra in long form,
#' the original \code{idfactor} is preserved.
#'
#' Spectra that are already in long form, if passed as argument, are returned
#' unchanged.
#'
#' Because the classes defined for spectra have a well defined format, and
#' known column names we can define a rather simple function for this
#' operation.
#'
#' @param spct An object with spectral data.
#' @param fixed.cols character Names of variables that should be copied
#'   unchanged for each spectrum.
#' @param idfactor character The name of the factor to be added to the long-form
#'   object and used to store the original name of the columns as an index
#'   to the different spectra.
#' @param rm.spct.class logical If true the returned object is a data frame.
#' @param ... Currently ignored.
#'
#' @return An object of the same class as \code{spct} or a \code{data.frame}
#'   with derived classes removed.
#'
#' @examples
#'
#' spct_wide2long(white_led.raw_spct)
#' spct_wide2long(white_led.cps_spct)
#' spct_wide2long(Ler_leaf.spct)
#'
#' @export
#'
spct_wide2long <- function(spct,
                           fixed.cols = "w.length",
                           idfactor = "spct.idx",
                           rm.spct.class = FALSE,
                           ...)
{
  orig.idfactor <- getIdFactor(spct)
  if (!is.na(orig.idfactor)) {
    if (idfactor == orig.idfactor) {
      stop("'idfactor' argument matches column name already in use")
    }
    fixed.cols <- c(fixed.cols, orig.idfactor)
  } else if (getMultipleWl(spct) > 1L) {
    if ("spct.idx" %in% colnames(spct)) {
      # created with old version
      orig.idfactor <- "spct.idx"
    } else {
      stop("'idfactor' not found for multiple spectra in long form")
    }
  }

  varying.cols <- setdiff(colnames(spct), fixed.cols)
  if (is.raw_spct(spct)) {
    varying.cols <- grep("^counts[0-9]*", sort(varying.cols), value = TRUE)
    values.to <- "counts"
  } else if (is.cps_spct(spct)) {
    varying.cols <- grep("^cps[0-9]*", sort(varying.cols), value = TRUE)
    values.to <- "cps"
  } else if (is.object_spct(spct)) {
    varying.cols <- intersect(c("Rfr", "Afr", "Tfr"), varying.cols)
    values.to <- "value"
  } else if (is.calibration_spct(spct)) {
    varying.cols <- grep("^irrad.mult", sort(varying.cols), value = TRUE)
    values.to <- "irrad.mult"
  } else if (is.generic_spct(spct)) {
    values.to <- "value"
  }
  num.cols <- length(varying.cols)
  if (num.cols > 1L) {
    orig.multiple.wl <- getMultipleWl(spct)
    old.class <- rmDerivedSpct(spct)
    spct <- spct[ , c(fixed.cols, varying.cols)]
    values.col <- length(fixed.cols) + 1L
    long.spct <- data.frame()
    for (i in varying.cols) {
      spct.slice <- spct[ , c(fixed.cols, i)]
      names(spct.slice)[values.col] <- values.to
      spct.slice[[idfactor]] <- i
      long.spct <- rbind(long.spct, spct.slice)
    }
    long.spct[[idfactor]] <- factor(long.spct[[idfactor]], levels = varying.cols)
    if (!rm.spct.class) {
      if (old.class[1] == "raw_spct") {
        setRawSpct(long.spct,
                   multiple.wl = num.cols * orig.multiple.wl,
                   idfactor = ifelse(is.na(orig.idfactor), idfactor, orig.idfactor))
      } else if (old.class[1] == "cps_spct") {
        setCpsSpct(long.spct,
                   multiple.wl = num.cols * orig.multiple.wl,
                   idfactor = ifelse(is.na(orig.idfactor), idfactor, orig.idfactor))
      } else if (old.class[1] == "calibration_spct") {
        setCalibrationSpct(long.spct,
                   multiple.wl = num.cols * orig.multiple.wl,
                   idfactor = ifelse(is.na(orig.idfactor), idfactor, orig.idfactor))
      } else {
        setGenericSpct(long.spct,
                       multiple.wl = num.cols * orig.multiple.wl,
                       idfactor = ifelse(is.na(orig.idfactor), idfactor, orig.idfactor))
      }
    }
    long.spct
  } else {
    spct
  }
}
