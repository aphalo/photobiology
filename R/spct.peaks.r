#' Find local maxima or global maximum (peaks)
#'
#' These functions find peaks (local maxima) and valleys (local minima) in a
#' numeric vector, using a user selectable span or window. Global and local
#' size thresholds based on different criteria make it possible restrict
#' the returned peaks to those more prominent. A \code{logical} vector is
#' returned.
#'
#' @param x numeric vector.
#' @param global.threshold numeric A value belonging to class \code{"AsIs"} is
#'   interpreted as an absolute minimum height or depth expressed in data units.
#'   A bare \code{numeric} value (normally between 0.0 and 1.0), is interpreted
#'   as relative to \code{threshold.range}. In both cases it sets a
#'   \emph{global} height (depth) threshold below which peaks (valleys) are
#'   ignored. A bare negative \code{numeric} value indicates the \emph{global}
#'   height (depth) threshold below which peaks (valleys) are be ignored. If
#'   \code{global.threshold = NULL}, no threshold is applied and all peaks
#'   returned.
#' @param local.threshold numeric A value belonging to class \code{"AsIs"} is
#'   interpreted as an absolute minimum height (depth) expressed in data units
#'   relative to a within-window computed reference value. A bare \code{numeric}
#'   value (normally between 0.0 and 1.0), is interpreted as expressed in units
#'   relative to \code{threshold.range}. In both cases \code{local.threshold}
#'   sets a \emph{local} height (depth) threshold below which peaks (valleys)
#'   are ignored. If \code{local.threshold = NULL} or if \code{span} spans the
#'   whole of \code{x}, no threshold is applied.
#' @param local.reference character One of \code{"median"}, \code{"median.log"},
#'   \code{"median.sqrt"}, \code{"farthest"}, \code{"farthest.log"} or
#'   \code{"farthest.sqrt"}. The reference used to assess the height of the
#'   peak, either the minimum/maximum value within the window or the median of
#'   all values in the window.
#' @param threshold.range numeric vector If of length 2 or a longer vector
#'   \code{range(threshold.range)} is used to scale both thresholds. With
#'   \code{NULL}, the default, \code{range(x)} is used, and with a vector of
#'   length one \code{range(threshold.range, x)} is used, i.e., the range
#'   is expanded.
#' @param span odd positive integer A peak is defined as an element in a
#'   sequence which is greater than all other elements within a moving window of
#'   width \code{span} centred at that element. The default value is 5, meaning
#'   that a peak is taller than its four nearest neighbours. \code{span = NULL}
#'   extends the span to the whole length of \code{x}.
#' @param strict logical flag: if \code{TRUE}, an element must be strictly
#'   greater than all other values in its window to be considered a peak.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for peaks.
#'
#' @return A vector of logical values of the same length as \code{x}. Values
#'   that are TRUE correspond to local peaks in vector \code{x} and can be used
#'   to extract the rows corresponding to peaks from a data frame.
#'
#' @details  As \code{\link[photobiology]{find_valleys}},
#'   \code{\link[photobiology]{peaks}} and \code{\link[photobiology]{valleys}}
#'   call \code{\link[photobiology]{find_peaks}} to search for peaks and
#'   valleys, this explanation applies to the four functions. It also applies to
#'   \code{\link[ggspectra]{stat_peaks}} and
#'   \code{\link[ggspectra:stat_peaks]{stat_valleys}}. Function
#'   \code{find_peaks} is a wrapper built onto function
#'   \code{\link[splus2R]{peaks}} from \pkg{splus2R}, adds support for peak
#'   height thresholds and handles \code{span = NULL} and non-finite (including
#'   NA) values differently than \code{splus2R::peaks}. Instead of giving an
#'   error when \code{na.rm = FALSE} and \code{x} contains \code{NA} values,
#'   \code{NA} values are replaced with the smallest finite value in \code{x}.
#'   \code{span = NULL} is treated as a special case and selects \code{max(x)}.
#'   Passing \code{strict = TRUE} ensures that non-unique global and within window
#'   maxima are ignored, and can result in no peaks being returned.
#'
#'   Two tests make it possible to ignore irrelevant peaks. One test
#'   (\code{global.threshold}) is based on the absolute height of the peaks and
#'   can be used in all cases to ignore globally low peaks. A second test
#'   (\code{local.threshold}) is available when the window defined by `span`
#'   does not include all observations and can be used to ignore peaks that are
#'   not locally prominent. In this second approach the height of each peak is
#'   compared to a summary computed from other values within the window of width
#'   equal to \code{span} where it was found. In this second case, the reference
#'   value used within each window containing a peak is given by the argument
#'   passed to \code{local.reference}. Parameter \code{threshold.range}
#'   determines how the values passed as argument to \code{global.threshold} and
#'   \code{local.threshold} are scaled. The default, \code{NULL} uses the range
#'   of \code{x}. Thresholds for ignoring too small peaks are applied after
#'   peaks are searched for, and threshold values can in some cases result in no
#'   peaks being returned.
#'
#'   The \code{local.threshold} argument is used \emph{as is} when
#'   \code{local.reference} is \code{"median"} or \code{"farthest"}, i.e., the
#'   same distance between peak and reference is used as cut-off irrespective of
#'   the value of the reference. In cases when the prominence of peaks is
#'   positively correlated with the baseline, a \code{local.threshold} that
#'   increases together with increasing computed within window median or
#'   farthest value applies apply a less stringent height requirement in regions
#'   with overall low height. In this case, natural logarithm or square root
#'   weighting can be requested with \code{local.reference} arguments
#'   \code{"median.log"}, \code{"farthest.log"}, \code{"median.sqrt"}, and
#'   \code{"farthest.sqrt"} as arguments for \code{local.reference}.
#'
#'   While functions \code{\link{find_peaks}} and \code{\link{find_valleys}}
#'   accept as input a \code{numeric} vector and return a \code{logical} vector,
#'   methods \code{\link{peaks}} and \code{\link{valleys}} accept as input
#'   different R objects, including spectra and collections of spectra and
#'   return a subset of the object. These methods are implemented using calls to
#'   functions \code{find_peaks}, \code{find_valleys} and
#'   \code{\link{fit_peaks}}.
#'
#' @note The default for parameter \code{strict} is \code{FALSE} in functions
#'   \code{\link{find_peaks}} and \code{\link{find_valleys}}, while the default
#'   in \code{\link[splus2R]{peaks}} is \code{strict = TRUE}.
#'
#' @seealso \code{\link[splus2R]{peaks}}.
#'
#' @family peaks and valleys functions
#'
#' @export
#'
#' @examples
#' with(sun.data, which(find_peaks(s.e.irrad, span = NULL)))
#' with(sun.data, which(find_peaks(s.e.irrad, span = 51)))
#' with(sun.data, w.length[find_peaks(s.e.irrad, span = 51)])
#' with(sun.data, sum(find_peaks(s.e.irrad, span = NULL, strict = TRUE)))
#'
#' with(sun.data, which(find_valleys(s.e.irrad, span = NULL)))
#' with(sun.data, which(find_valleys(s.e.irrad, span = 51)))
#'
find_peaks <-
  function(x,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           span = 3,
           strict = FALSE,
           na.rm = FALSE) {
    # validate parameters
    if (length(unique(na.omit(x))) < 2L) {
      return(logical(length(x)))
    }
    if (is.null(span) || span >= length(x)) {
      # ignore local.threshold argument when not applicable
      local.threshold <- NULL
    }
    if (!is.null(span) && (is.na(span) || !is.numeric(span) || span < 1)) {
      stop("'span' must be NULL or a positive odd integer, not: ", format(span))
    }
    if (!is.null(local.threshold)) {
      if (is.na(local.threshold)) {
        return(rep(NA, length(x)))
      }
      if (!is.numeric(local.threshold) || length(local.threshold) > 1L ||
          local.threshold < 0 || local.threshold > 1) {
        stop("'local.threshold' must be NULL or a single number in [0..1], not: ",
             local.threshold)
      }
    }

    if (!is.null(global.threshold)) {
      if (is.na(global.threshold)) {
        return(rep(NA, length(x)))
      }
      if (!is.numeric(global.threshold) || length(global.threshold) > 1L) {
        stop("'global.threshold' must be NULL or a single number, not: ",
             local.threshold)
      } else if (!inherits(global.threshold, "AsIs") &&
                 (global.threshold < -1 || global.threshold > 1)) {
        stop("'global.threshold' when not \"AsIs\" must be a number in [-1..1]",
             "not: ", global.threshold)
      }
      if (inherits(global.threshold, "AsIs") && !is.finite(global.threshold)) {
        # accept all peaks/valleys
        global.threshold <- NULL
      }
    }

    # Replace NA, NaN and Inf with smallest finite value
    if (na.rm) {
      x <- ifelse(!is.finite(x), min(x, na.rm = TRUE), x)
    }

    # compute threshold range only if needed
    if (!is.null(global.threshold) || !is.null(local.threshold)) {
      if (is.null(threshold.range)) {
        threshold.range <- range(x, na.rm = TRUE)
      } else if (length(threshold.range) == 1L) {
        threshold.range <- range(threshold.range, x, na.rm = TRUE)
      } else {
        threshold.range <- range(threshold.range, na.rm = TRUE)
      }
      if (length(threshold.range) != 2L ||
          (threshold.range[2] - threshold.range[1]) < 1e-16) {
        warning("Skipping! Bad 'threshold.range': [",
                paste(threshold.range, collapse = ", "), "]")
        return(logical(length(x)))
      }
    }
    # compute threshold range only if needed
    if (!is.null(global.threshold) || !is.null(local.threshold)) {
      if (is.null(threshold.range)) {
        threshold.range <- range(x, na.rm = TRUE)
      } else if (length(threshold.range) == 1L) {
        threshold.range <- range(threshold.range, x, na.rm = TRUE)
      } else {
        threshold.range <- range(threshold.range, na.rm = TRUE)
      }
      if (length(threshold.range) != 2L ||
          (threshold.range[2] - threshold.range[1]) < 1e-16) {
        warning("Skipping! Bad 'threshold.range': [",
                paste(threshold.range, collapse = ", "), "]")
        return(logical(length(x)))
      }
    }
    # compute global multiplier and base only if needed
    if (!is.null(global.threshold)) {
      if (inherits(global.threshold, "AsIs")) {
        global.multiplier <- 1
        global.base <- 0
      } else {
        global.multiplier <- abs(threshold.range[2] - threshold.range[1])
        global.base <- threshold.range[1]
      }
    }
    # compute local multiplier only if needed
    if (!is.null(local.threshold)) {
      if (inherits(local.threshold, "AsIs")) {
        local.multiplier <- 1
      } else {
        local.multiplier <- threshold.range[2] - threshold.range[1]
      }
    }

    # search for global or local maxima
    if (is.null(span) || span >= length(x)) {
      # find maximum
      pks <- (x == max(x, na.rm = na.rm))
      if (strict && sum(pks) > 1L) {
        message("Non-unique extreme value discarded as 'strict = TRUE'")
        pks <- logical(length(x)) # all FALSE
      }
    } else {
      # search for local maxima
      pks <-
        splus2R::peaks(x = x, span = span, strict = strict, endbehavior = 0)
    }

    # apply global height threshold test
    if (!is.null(global.threshold)) {
      if (global.threshold >= 0 || inherits(global.threshold, "AsIs")) {
        pks <- pks &
          x > global.base + global.threshold * global.multiplier
      } else {
        pks <- pks &
          x <= global.base + abs(global.threshold) * global.multiplier
      }
    }

    # apply local.threshold height test to found peaks
    if (!is.null(local.threshold)) {
      # apply local height threshold test
      local.base <-
        switch(local.reference,
               median =,
               median.sqrt =,
               median.log = stats::runmed(x, k = span, endrule = "median"),
               farthest =,
               farthest.sqrt =,
               farthest.log = caTools::runmin(x, k = span, endrule = "min"),
               stop("Bad 'local.reference': ", local.reference))
      if (grepl("\\.log$", local.reference)) {
        pks <- pks & log(x - local.base) >
          (local.threshold * log(local.multiplier))
      } else if (grepl("\\.sqrt$", local.reference)) {
        pks <- pks & sqrt(x - local.base) >
          (local.threshold * sqrt(local.multiplier))
      } else {
        pks <- pks & (x - local.base) > (local.threshold * local.multiplier)
      }
    }
    pks
  }

#' @rdname find_peaks
#'
#' @export
#'
find_valleys <-
  function(x,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           span = 3,
           strict = FALSE,
           na.rm = FALSE) {
    x <- -x
    if (inherits(global.threshold, "AsIs")) {
      global.threshold <- I(-global.threshold)
    }
    if (inherits(local.threshold, "AsIs")) {
      local.threshold <- I(-local.threshold)
    }
    if (!is.null(threshold.range)) {
      threshold.range <- rev(-threshold.range)
    }
    find_peaks(x = x,
               global.threshold = global.threshold,
               local.threshold = local.threshold,
               local.reference = local.reference,
               threshold.range = threshold.range,
               span = span,
               strict = strict,
               na.rm = na.rm)
  }

# fit peaks ---------------------------------------------------------------

#' Refine position and value of extremes by fitting
#'
#' Functions implementing fitting of peaks in a class-agnostic way. The fitting
#' refines the location of peaks and value of peaks based on the location of
#' maxima and minima supplied. This function is to be used together with
#' \code{find_peaks()} or \code{find_valleys()}.
#'
#' @param x generic_spct or data.frame object.
#' @param peaks.idx,valleys.idx logical or integer Indexes into \code{x}
#'   selecting global or local extremes.
#' @param span odd integer The span used when refining the location of maxima or
#'   minima of \code{x}.
#' @param x.col.name,y.col.name character Name of the column of \code{x} on
#'   which to operate.
#' @param method character The method to use for the fit.
#' @param max.span odd integer The maximum number of data points used when
#'   when refining the location of maxima and minima.
#' @param maximum logical A flag indicating whether to search for maxima or
#'   minima.
#' @param keep.cols logical Keep unrecognized columns in data frames
#'
#' @note These functions are not meant for everyday use. Use option
#'   \code{refine.wl = TRUE} of methods \code{peaks()} and \code{valleys()}
#'   instead.
#'
#' @details The only method currently implemented is \code{"spline"} based on
#'   a call to \code{\link[stats]{splinefun}} in a window of width \code{span}
#'   centred on each peak pointed at by \code{peaks.idx}. A spline fitted to
#'   a narrow window will usually locate the position of the peak in the
#'   column named by the argument passed to \code{x.col.name} better than
#'   estimating the true height of the peak in the column named by the argument
#'   passed to \code{y.col.name}.
#'
#' @return An R object of the same class as \code{x} containing the fitted
#'   values for the peaks, and optionally the unmodified values at the rows
#'   matching \code{peaks.idx} or \code{valleys.idx} for other retained columns.
#'
#' @examples
#'
#' peaks <- find_peaks(sun.spct[["s.e.irrad"]], span = 31)
#' fit_peaks(sun.spct, peaks, span = 31,
#'           y.col.name = "s.e.irrad", method = "spline")
#'
#' @export
#'
fit_peaks <- function(x,
                      peaks.idx,
                      span,
                      x.col.name = NULL,
                      y.col.name,
                      method,
                      max.span = 5L,
                      maximum = TRUE,
                      keep.cols = NULL) {
  if (is.null(span)) {
    span <- max.span
  }
  if (is.logical(peaks.idx)) {
    peaks.idx <- which(peaks.idx)
  }
  if (is.null(x.col.name) && is.any_spct(x)) {
    x.col.name <- "w.length"
  }
  if (is.any_spct(x)) {
    x <- untag(x)
  }
  if (is.null(keep.cols) && !is.any_spct(x)) {
    # treat data frames as special case, needed by ggspectra::stat_peaks().
    z <- x[peaks.idx, ]
  } else {
    # delete numeric columns which would be invalidated
    cols2rm <- names(x)[sapply(X = x, FUN = is.numeric)]
    cols2rm <- setdiff(cols2rm, c(x.col.name, y.col.name, keep.cols))
    z <- x[peaks.idx , setdiff(colnames(x), cols2rm)]
  }
  if (method == "spline") {
    f <- stats::splinefun(x[[x.col.name]], x[[y.col.name]])
  } else {
    stop("'method' ", method, " is not implemented")
  }
  w.length <- numeric()
  var <- numeric()
  # interval should not be wider than span used to locate maxima
  half.interval <- min(span %/% 2L, max.span)
  for (p in peaks.idx) {
    # we need to avoid off-range indexes!
    interval.p <- c(x[[x.col.name]][max(p - half.interval, 0L)],
                    x[[x.col.name]][min(p + half.interval, nrow(x))])
    temp <- stats::optimize(f,
                            interval = interval.p,
                            maximum = maximum)
    if (maximum) {
      w.length <- c(w.length, temp[["maximum"]])
    } else {
      w.length <- c(w.length, temp[["minimum"]])
    }
    var <- c(var, temp[["objective"]])
  }
  # replace columns of same name
  z[[x.col.name]] <- w.length
  z[[y.col.name]] <- var
  if (is.any_spct(x)) {
     z <- copy_attributes(x, z, copy.class = FALSE)
  }
  z
}

#' @rdname fit_peaks
#'
#' @export
#'
fit_valleys <- function(x,
                        valleys.idx,
                        span,
                        x.col.name = NULL,
                        y.col.name,
                        method,
                        max.span = 5L,
                        maximum = FALSE,
                        keep.cols = NULL) {
  fit_peaks(x = x,
            peaks.idx = valleys.idx,
            span = span,
            x.col.name = x.col.name,
            y.col.name = y.col.name,
            method = method,
            max.span = max.span,
            maximum = maximum,
            keep.cols = keep.cols)
}

# peaks -------------------------------------------------------------------

#' Peaks or local maxima
#'
#' Function that returns a subset of an R object with observations corresponding
#' to local maxima.
#'
#' @inheritParams find_peaks
#' @param var.name,x.var.name,y.var.name character Name of column where to look
#'   for peaks.
#' @param refine.wl logical Flag indicating if peak location should be refined
#'   by fitting a function.
#' @param method character String with the name of a method. Currently only
#'   spline interpolation is implemented.
#' @param ... ignored
#'
#' @inherit find_peaks details note seealso
#'
#' @return A subset of \code{x} with rows corresponding to local maxima.
#'
#' @export
#'
#' @examples
#' # default span = 5
#' peaks(sun.spct)
#' # global maximum
#' peaks(sun.spct, span = NULL)
#' peaks(sun.spct, span = NULL)$w.length
#' # fitted peak wavelength
#' peaks(sun.spct, span = NULL, refine.wl = TRUE)
#' peaks(sun.spct, span = NULL, refine.wl = TRUE)$w.length
#' # a wider window
#' peaks(sun.spct, span = 51)
#' # global threshold relative to the range of s.e.irrad values
#' peaks(sun.spct, global.threshold = 0.7)
#' peaks(sun.spct, global.threshold = -0.3)
#' # global threshold in actual s.e.irrad values
#' peaks(sun.spct, global.threshold = 0.7, threshold.range = c(0, 1))
#' # local threshold  relative to the range of s.e.irrad values
#' peaks(sun.spct, local.threshold = 0.1)
#' # local threshold in actual s.e.irrad values
#' peaks(sun.spct, local.threshold = 0.1, threshold.range = c(0, 1))
#' # local threshold  relative to the range of s.e.irrad values, using window
#' # median instead of window minimum
#' peaks(sun.spct, local.threshold = 0.05, local.reference = "median")
#' # minimum, the default.
#' peaks(sun.spct, local.threshold = 0.05, local.reference = "farthest")
#'
#' @family peaks and valleys functions
#'
peaks <- function(x,
                  span,
                  global.threshold,
                  local.threshold,
                  local.reference,
                  threshold.range,
                  strict,
                  na.rm,
                  ...) UseMethod("peaks")

#' @rdname peaks
#'
#' @export
#'
peaks.default <-
  function(x,
           span = NA,
           global.threshold = NA,
           local.threshold = NA,
           local.reference = NA,
           threshold.range = NA,
           strict = NA,
           na.rm = FALSE,
           ...) {
  warning("Method 'peaks' not implemented for objects of class ", class(x)[1])
  x[NA]
}

#' @rdname peaks
#'
#' @export
#'
peaks.numeric <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           ...) {
  x[find_peaks(x = x,
               span = span,
               global.threshold = global.threshold,
               local.threshold = local.threshold,
               local.reference = local.reference,
               threshold.range = threshold.range,
               strict = strict,
               na.rm = na.rm)]
}

#' @rdname peaks
#'
#' @export
#'
peaks.data.frame <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           x.var.name = NULL,
           y.var.name = NULL,
           var.name = y.var.name,
           refine.wl = FALSE,
           method = "spline",
           ...) {
    if (is.null(var.name) || (refine.wl && is.null(x.var.name))) {
      warning("Variable (column) names required.")
      return(x[NA, ])
    }

    # vector is not necessarily sorted
    if (!is.null(x.var.name)) {
      if (is.unsorted(x[[x.var.name]])) {
        message("In 'peaks.data.frame()' x variable is not sorted")
      } else {
        check_wl_stepsize(x = x[[x.var.name]],
                          span = span,
                          na.rm = TRUE,
                          min.stepsize = 3)
      }
    }

    peaks.idx <-
      which(find_peaks(x[[var.name]],
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm))
    if (refine.wl && length(peaks.idx > 0L)) {
      fit_peaks(x = x,
                peaks.idx = peaks.idx,
                span = span,
                x.col.name = x.var.name,
                y.col.name = var.name,
                method = method)
    } else {
      x[peaks.idx,  , drop = FALSE]
    }
  }

#' @rdname peaks
#'
#' @export
#'
peaks.generic_spct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           var.name = NULL,
           refine.wl = FALSE,
           method = "spline",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- peaks(x = mspct,
                     span = span,
                     global.threshold = global.threshold,
                     local.threshold = local.threshold,
                     local.reference = local.reference,
                     threshold.range = threshold.range,
                     strict = strict,
                     na.rm = na.rm,
                     var.name = var.name,
                     refine.wl = refine.wl,
                     method = method,
                     ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

    if (is.null(var.name)) {
      # find target variable
      var.name <- names(x)
      var.name <- subset(var.name, sapply(x, is.numeric))
      var.name <- setdiff(var.name, "w.length")
      if (length(var.name) > 1L) {
        warning("Multiple numeric data columns found, explicit argument to 'var.name' required.")
        return(x[NA, ])
      }
    }
    peaks.idx <-
      which(find_peaks(x[[var.name]],
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm))
    if (refine.wl && length(peaks.idx > 0L)) {
      fit_peaks(x = x,
                peaks.idx = peaks.idx,
                span = span,
                y.col.name = var.name,
                method = method)
    } else {
      x[peaks.idx,  , drop = FALSE]
    }
  }

#' @rdname peaks
#'
#' @param unit.out character One of "energy" or "photon"
#'
#' @export
#'
#' @examples
#' peaks(sun.spct)
#'
peaks.source_spct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- peaks(x = mspct,
                     span = span,
                     global.threshold = global.threshold,
                     local.threshold = local.threshold,
                     local.reference = local.reference,
                     threshold.range = threshold.range,
                     strict = strict,
                     na.rm = na.rm,
                     unit.out = unit.out,
                     refine.wl = refine.wl,
                     method = method,
                     ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

    if (unit.out == "energy") {
      z <- q2e(x, action = "replace", byref = FALSE)
      col.name <- "s.e.irrad"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, action = "replace", byref = FALSE)
      col.name <- "s.q.irrad"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    peaks.idx <-
      which(find_peaks(z[[col.name]],
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm))
    if (refine.wl && length(peaks.idx > 0L)) {
      fit_peaks(x = z,
                peaks.idx = peaks.idx,
                span = span,
                y.col.name = col.name,
                method = method)
    } else {
      z[peaks.idx,  , drop = FALSE]
    }
  }

#' @rdname peaks
#'
#' @export
#'
peaks.response_spct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- peaks(x = mspct,
                     span = span,
                     global.threshold = global.threshold,
                     local.threshold = local.threshold,
                     local.reference = local.reference,
                     threshold.range = threshold.range,
                     strict = strict,
                     na.rm = na.rm,
                     unit.out = unit.out,
                     refine.wl = refine.wl,
                     method = method,
                     ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.response"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.response"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    peaks.idx <-
      which(find_peaks(z[[col.name]],
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm))
    if (refine.wl && length(peaks.idx > 0L)) {
      fit_peaks(x = z,
                peaks.idx = peaks.idx,
                span = span,
                y.col.name = col.name,
                method = method)
    } else {
      z[peaks.idx,  , drop = FALSE]
    }
  }

#' @rdname peaks
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
peaks.filter_spct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty",
                                  default = "transmittance"),
           refine.wl = FALSE,
           method = "spline",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- peaks(x = mspct,
                     span = span,
                     global.threshold = global.threshold,
                     local.threshold = local.threshold,
                     local.reference = local.reference,
                     threshold.range = threshold.range,
                     strict = strict,
                     na.rm = na.rm,
                     filter.qty = filter.qty,
                     refine.wl = refine.wl,
                     method = method,
                     ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

    if (filter.qty == "transmittance") {
      z <- A2T(x, "replace", FALSE)
      col.name <- "Tfr"
    } else if (filter.qty == "absorbance") {
      z <- T2A(x, "replace", FALSE)
      col.name <- "A"
    } else {
      stop("Unrecognized 'filter.qty': ", filter.qty)
    }
    peaks.idx <-
      which(find_peaks(z[[col.name]],
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm))
    if (refine.wl && length(peaks.idx > 0L)) {
      fit_peaks(x = z,
                peaks.idx = peaks.idx,
                span = span,
                y.col.name = col.name,
                method = method)
    } else {
      z[peaks.idx,  , drop = FALSE]
    }
  }

#' @rdname peaks
#'
#' @export
#'
peaks.reflector_spct <- function(x,
                                 span = 5,
                                 global.threshold = NULL,
                                 local.threshold = NULL,
                                 local.reference = "median",
                                 threshold.range = NULL,
                                 strict = FALSE,
                                 na.rm = FALSE,
                                 refine.wl = FALSE,
                                 method = "spline",
                                 ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(x) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    # call method on the collection
    mspct <- peaks(x = mspct,
                   span = span,
                   global.threshold = global.threshold,
                   local.threshold = local.threshold,
                   local.reference = local.reference,
                   threshold.range = threshold.range,
                   strict = strict,
                   na.rm = na.rm,
                   refine.wl = refine.wl,
                   method = method,
                   ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

  col.name <- "Rfr"
  peaks.idx <-
    which(find_peaks(x[[col.name]],
                     span = span,
                     global.threshold = global.threshold,
                     local.threshold = local.threshold,
                     local.reference = local.reference,
                     threshold.range = threshold.range,
                     strict = strict,
                     na.rm = na.rm))
  if (refine.wl && length(peaks.idx > 0L)) {
    fit_peaks(x = x,
              peaks.idx = peaks.idx,
              span = span,
              y.col.name = col.name,
              method = method)
  } else {
    x[peaks.idx,  , drop = FALSE]
  }
}

#' @rdname peaks
#'
#' @export
#'
peaks.solute_spct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           refine.wl = FALSE,
           method = "spline",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- peaks(x = mspct,
                     span = span,
                     global.threshold = global.threshold,
                     local.threshold = local.threshold,
                     local.reference = local.reference,
                     threshold.range = threshold.range,
                     strict = strict,
                     na.rm = na.rm,
                     refine.wl = refine.wl,
                     method = method,
                     ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

    cols <- intersect(c("K.mole", "K.mass"), names(x))
    if (length(cols) == 1) {
      col.name <- cols
      z <- x
    } else {
      stop("Invalid number of columns found:", length(cols))
    }
    peaks.idx <-
      which(find_peaks(z[[col.name]],
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm))
    if (refine.wl && length(peaks.idx > 0L)) {
      fit_peaks(x = z,
                peaks.idx = peaks.idx,
                span = span,
                y.col.name = col.name,
                method = method)
    } else {
      z[peaks.idx,  , drop = FALSE]
    }
  }

#' @rdname peaks
#'
#' @export
#'
peaks.cps_spct <- function(x,
                           span = 5,
                           global.threshold = NULL,
                           local.threshold = NULL,
                           local.reference = "median",
                           threshold.range = NULL,
                           strict = FALSE,
                           na.rm = FALSE,
                           var.name = "cps",
                           refine.wl = FALSE,
                           method = "spline",
                           ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(x) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    # call method on the collection
    mspct <- peaks(x = mspct,
                   span = span,
                   global.threshold = global.threshold,
                   local.threshold = local.threshold,
                   local.reference = local.reference,
                   threshold.range = threshold.range,
                   strict = strict,
                   na.rm = na.rm,
                   var.name = var.name,
                   refine.wl = refine.wl,
                   method = method,
                   ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

  peaks.idx <-
    which(find_peaks(x[[var.name]],
                     span = span,
                     global.threshold = global.threshold,
                     local.threshold = local.threshold,
                     local.reference = local.reference,
                     threshold.range = threshold.range,
                     strict = strict,
                     na.rm = na.rm))
  if (refine.wl && length(peaks.idx > 0L)) {
    fit_peaks(x = x,
              peaks.idx = peaks.idx,
              span = span,
              y.col.name = var.name,
              method = method)
  } else {
    x[peaks.idx,  , drop = FALSE]
  }
}

#' @rdname peaks
#'
#' @export
#'
peaks.raw_spct <- function(x, span = 5,
                           global.threshold = NULL,
                           local.threshold = NULL,
                           local.reference = "median",
                           threshold.range = NULL,
                           strict = FALSE,
                           na.rm = FALSE,
                           var.name = "counts",
                           refine.wl = FALSE,
                           method = "spline",
                           ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(x) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    # call method on the collection
    mspct <- peaks(x = mspct,
                   span = span,
                   global.threshold = global.threshold,
                   local.threshold = local.threshold,
                   local.reference = local.reference,
                   threshold.range = threshold.range,
                   strict = strict,
                   na.rm = na.rm,
                   var.name = var.name,
                   refine.wl = refine.wl,
                   method = method,
                   ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

  peaks.idx <-
    which(find_peaks(x[[var.name]],
                     span = span,
                     global.threshold = global.threshold,
                     local.threshold = local.threshold,
                     local.reference = local.reference,
                     threshold.range = threshold.range,
                     strict = strict,
                     na.rm = na.rm))
  if (refine.wl && length(peaks.idx > 0L)) {
    fit_peaks(x = x,
              peaks.idx = peaks.idx,
              span = span,
              y.col.name = var.name,
              method = method)
  } else {
    x[peaks.idx,  , drop = FALSE]
  }
}


# _mspct ------------------------------------------------------------------


#' @rdname peaks
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
peaks.generic_mspct <- function(x,
                                span = 5,
                                global.threshold = NULL,
                                local.threshold = NULL,
                                local.reference = "median",
                                threshold.range = NULL,
                                strict = FALSE,
                                na.rm = FALSE,
                                var.name = NULL,
                                refine.wl = FALSE,
                                method = "spline",
                                ...,
                                .parallel = FALSE,
                                .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(x,
          .fun = peaks,
          span = span,
          global.threshold = global.threshold,
          local.threshold = local.threshold,
          local.reference = local.reference,
          threshold.range = threshold.range,
          strict = strict,
          na.rm = na.rm,
          var.name = var.name,
          refine.wl = refine.wl,
          method = method,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
  }

#' @rdname peaks
#'
#' @export
#'
peaks.source_mspct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = peaks,
            span = span,
            global.threshold = global.threshold,
            local.threshold = local.threshold,
            local.reference = local.reference,
            threshold.range = threshold.range,
            strict = strict,
            unit.out = unit.out,
            na.rm = na.rm,
            refine.wl = refine.wl,
            method = method,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @rdname peaks
#'
#' @export
#'
peaks.response_mspct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = peaks,
            span = span,
            global.threshold = global.threshold,
            local.threshold = local.threshold,
            local.reference = local.reference,
            threshold.range = threshold.range,
            strict = strict,
            unit.out = unit.out,
            na.rm = na.rm,
            refine.wl = refine.wl,
            method = method,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @rdname peaks
#'
#' @export
#'
peaks.filter_mspct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty",
                                  default = "transmittance"),
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = peaks,
            span = span,
            global.threshold = global.threshold,
            local.threshold = local.threshold,
            local.reference = local.reference,
            threshold.range = threshold.range,
            strict = strict,
            filter.qty = filter.qty,
            na.rm = na.rm,
            refine.wl = refine.wl,
            method = method,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }


#' @rdname peaks
#'
#' @export
#'
peaks.reflector_mspct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = peaks,
            span = span,
            global.threshold = global.threshold,
            local.threshold = local.threshold,
            local.reference = local.reference,
            threshold.range = threshold.range,
            strict = strict,
            na.rm = na.rm,
            refine.wl = refine.wl,
            method = method,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }


#' @rdname peaks
#'
#' @export
#'
peaks.solute_mspct <- peaks.reflector_mspct

#' @rdname peaks
#'
#' @export
#'
peaks.cps_mspct <- function(x,
                            span = 5,
                            global.threshold = NULL,
                            local.threshold = NULL,
                            local.reference = "median",
                            threshold.range = NULL,
                            strict = FALSE,
                            na.rm = FALSE,
                            var.name = "cps",
                            refine.wl = FALSE,
                            method = "spline",
                            ...,
                            .parallel = FALSE,
                            .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(x,
          .fun = peaks,
          span = span,
          global.threshold = global.threshold,
          local.threshold = local.threshold,
          local.reference = local.reference,
          threshold.range = threshold.range,
          strict = strict,
          na.rm = na.rm,
          var.name = var.name,
          refine.wl = refine.wl,
          method = method,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @rdname peaks
#'
#' @export
#'
peaks.raw_mspct <- function(x,
                            span = 5,
                            global.threshold = NULL,
                            local.threshold = NULL,
                            local.reference = "median",
                            threshold.range = NULL,
                            strict = FALSE,
                            na.rm = FALSE,
                            var.name = "counts",
                            refine.wl = FALSE,
                            method = "spline",
                            ...,
                            .parallel = FALSE,
                            .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(x,
          .fun = peaks,
          span = span,
          global.threshold = global.threshold,
          local.threshold = local.threshold,
          local.reference = local.reference,
          threshold.range = threshold.range,
          strict = strict,
          na.rm = na.rm,
          var.name = var.name,
          refine.wl = refine.wl,
          method = method,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

# valleys -------------------------------------------------------------------

#' Valleys or local minima
#'
#' Function that returns a subset of an R object with observations corresponding
#' to local maxima.
#'
#' @inheritParams find_peaks
#' @param var.name,x.var.name,y.var.name character Name of column where to look
#'   for valleys.
#' @param refine.wl logical Flag indicating if valley location should be refined by
#'   fitting a function.
#' @param method character String with the name of a method. Currently only
#'   spline interpolation is implemented.
#' @param ... ignored
#'
#' @inherit find_peaks details note seealso
#'
#' @return A subset of \code{x} with rows corresponding to local minima or
#'   global minimum.
#'
#' @examples
#' # default span = 5
#' valleys(sun.spct)
#' # global minimum
#' valleys(sun.spct, span = NULL)
#' valleys(sun.spct, span = NULL, strict = FALSE)
#' # a wider window
#' valleys(sun.spct, span = 51)
#' # global threshold relative to the range of s.e.irrad values
#' valleys(sun.spct, global.threshold = -0.2)
#' # global threshold in actual s.e.irrad values
#' valleys(sun.spct, global.threshold = -0.2, threshold.range = c(0, 1))
#' # local threshold  relative to the range of s.e.irrad values
#' valleys(sun.spct, local.threshold = 0.1)
#' # local threshold in actual s.e.irrad values
#' valleys(sun.spct, local.threshold = 0.1, threshold.range = c(0, 1))
#' # local threshold  relative to the range of s.e.irrad values, using window
#' # median instead of window minimum
#' valleys(sun.spct, local.threshold = 0.1, local.reference = "median")
#' # minimum, the default.
#' valleys(sun.spct, local.threshold = 0.1, local.reference = "farthest")
#'
#' @export
#'
#' @family peaks and valleys functions
#'
valleys <- function(x,
                    span,
                    global.threshold,
                    local.threshold,
                    local.reference,
                    threshold.range,
                    strict,
                    ...) UseMethod("valleys")

#' @rdname valleys
#'
#' @export
valleys.default <- function(x,
                            span,
                            global.threshold = NA,
                            local.threshold = NA,
                            local.reference = NA,
                            threshold.range = NA,
                            strict,
                            ...) {
  warning("Method 'valleys' not implemented for objects of class ", class(x)[1])
  x[NA]
}

#' @rdname valleys
#'
#' @export
valleys.numeric <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           ...) {
    x[find_valleys(x = x,
                   span = span,
                   global.threshold = global.threshold,
                   local.threshold = local.threshold,
                   local.reference = local.reference,
                   threshold.range = threshold.range,
                   strict = strict,
                   na.rm = na.rm)]
  }

#' @rdname valleys
#'
#' @export
#'
valleys.data.frame <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           x.var.name = NULL,
           y.var.name = NULL,
           var.name = y.var.name,
           refine.wl = FALSE,
           method = "spline",
           ...) {

    if (is.null(var.name) || (refine.wl && is.null(x.var.name))) {
      warning("Variable (column) names required.")
      return(x[NA, ])
    }
    valleys.idx <-
      which(find_valleys(x = x[[var.name]],
                         span = span,
                         global.threshold = global.threshold,
                         local.threshold = local.threshold,
                         local.reference = local.reference,
                         threshold.range = threshold.range,
                         strict = strict))
    if (refine.wl && length(valleys.idx > 0L)) {
      fit_valleys(x = x,
                valleys.idx = valleys.idx,
                span = span,
                x.col.name = x.var.name,
                y.col.name = y.var.name,
                method = method)
    } else {
      x[valleys.idx,  , drop = FALSE]
    }
  }

#' @rdname valleys
#'
#' @export
#'
valleys.generic_spct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           var.name = NULL,
           refine.wl = FALSE,
           method = "spline",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- valleys(x = mspct,
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm,
                       var.name = var.name,
                       refine.wl = refine.wl,
                       method = method,
                       ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

    if (is.null(var.name)) {
      # find target variable
      var.name <- names(x)
      var.name <- subset(var.name, sapply(x, is.numeric))
      var.name <- setdiff(var.name, "w.length")
      if (length(var.name) > 1L) {
        warning("Multiple numeric data columns found, explicit argument to 'var.name' required.")
        return(x[NA, ])
      }
    }
    valleys.idx <-
      which(find_valleys(x = x[[var.name]],
                         span = span,
                         global.threshold = global.threshold,
                         local.threshold = local.threshold,
                         local.reference = local.reference,
                         threshold.range = threshold.range,
                         strict = strict,
                         na.rm = na.rm))
    if (refine.wl && length(valleys.idx > 0L)) {
      fit_valleys(x = x,
                  valleys.idx = valleys.idx,
                  span = span,
                  y.col.name = var.name,
                  method = method)
    } else {
      x[valleys.idx,  , drop = FALSE]
    }
  }

#' @rdname valleys
#'
#' @param unit.out character One of "energy" or "photon"
#'
#' @export
#'
#' @examples
#' valleys(sun.spct)
#'
valleys.source_spct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "median",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- valleys(x = mspct,
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm,
                       unit.out = unit.out,
                       refine.wl = refine.wl,
                       method = method,
                       ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.irrad"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.irrad"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    valleys.idx <-
      which(find_valleys(x = z[[col.name]],
                         span = span,
                         global.threshold = global.threshold,
                         local.threshold = local.threshold,
                         local.reference = local.reference,
                         threshold.range = threshold.range,
                         strict = strict,
                         na.rm = na.rm))
    if (refine.wl && length(valleys.idx > 0L)) {
      fit_valleys(x = z,
                  valleys.idx = valleys.idx,
                  span = span,
                  y.col.name = col.name,
                  method = method)
    } else {
      z[valleys.idx,  , drop = FALSE]
    }
  }

#' @rdname valleys
#'
#' @export
#'
valleys.response_spct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- valleys(x = mspct,
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm,
                       unit.out = unit.out,
                       refine.wl = refine.wl,
                       method = method,
                       ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.response"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.response"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    valleys.idx <-
      which(find_valleys(z[[col.name]],
                         span = span,
                         global.threshold = global.threshold,
                         local.threshold = local.threshold,
                         local.reference = local.reference,
                         threshold.range = threshold.range,
                         strict = strict,
                         na.rm = na.rm))
    if (refine.wl && length(valleys.idx > 0L)) {
      fit_valleys(x = z,
                  valleys.idx = valleys.idx,
                  span = span,
                  y.col.name = col.name,
                  method = method)
    } else {
      z[valleys.idx,  , drop = FALSE]
    }
  }

#' @rdname valleys
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
valleys.filter_spct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty", default = "transmittance"),
           refine.wl = FALSE,
           method = "spline",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- valleys(x = mspct,
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm,
                       filter.qty = filter.qty,
                       refine.wl = refine.wl,
                       method = method,
                       ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

    if (filter.qty == "transmittance") {
      z <- A2T(x, "replace", FALSE)
      col.name <- "Tfr"
    } else if (filter.qty == "absorbance") {
      z <- T2A(x, "replace", FALSE)
      col.name <- "A"
    } else {
      stop("Unrecognized 'filter.qty': ", filter.qty)
    }
    valleys.idx <-
      which(find_valleys(z[[col.name]],
                         span = span,
                         global.threshold = global.threshold,
                         local.threshold = local.threshold,
                         local.reference = local.reference,
                         threshold.range = threshold.range,
                         strict = strict,
                         na.rm = na.rm))
    if (refine.wl && length(valleys.idx > 0L)) {
      fit_valleys(x = z,
                  valleys.idx = valleys.idx,
                  span = span,
                  y.col.name = col.name,
                  method = method)
    } else {
      z[valleys.idx,  , drop = FALSE]
    }
  }

#' @rdname valleys
#'
#' @export
#'
valleys.reflector_spct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           refine.wl = FALSE,
           method = "spline",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- valleys(x = mspct,
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm,
                       refine.wl = refine.wl,
                       method = method,
                       ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

    col.name <- "Rfr"
    valleys.idx <-
      which(find_valleys(x[[col.name]],
                         span = span,
                         global.threshold = global.threshold,
                         local.threshold = local.threshold,
                         local.reference = local.reference,
                         threshold.range = threshold.range,
                         strict = strict,
                         na.rm = na.rm))
    if (refine.wl && length(valleys.idx > 0L)) {
      fit_valleys(x = x,
                  valleys.idx = valleys.idx,
                  span = span,
                  y.col.name = col.name,
                  method = method)
    } else {
      x[valleys.idx,  , drop = FALSE]
    }
  }

#' @rdname valleys
#'
#' @export
#'
valleys.solute_spct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           refine.wl = FALSE,
           method = "spline",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- valleys(x = mspct,
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm,
                       refine.wl = refine.wl,
                       method = method,
                       ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

    cols <- intersect(c("K.mole", "K.mass"), names(x))
    if (length(cols) == 1) {
      col.name <- cols
      z <- x
    } else {
      stop("Invalid number of columns found:", length(cols))
    }
    valleys.idx <-
      which(find_valleys(z[[col.name]],
                         span = span,
                         global.threshold = global.threshold,
                         local.threshold = local.threshold,
                         local.reference = local.reference,
                         threshold.range = threshold.range,
                         strict = strict,
                         na.rm = na.rm))
    if (refine.wl && length(valleys.idx > 0L)) {
      fit_valleys(x = z,
                  valleys.idx = valleys.idx,
                  span = span,
                  y.col.name = col.name,
                  method = method)
    } else {
      z[valleys.idx,  , drop = FALSE]
    }
  }

#' @rdname valleys
#'
#' @export
#'
valleys.cps_spct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           var.name = "cps",
           refine.wl = FALSE,
           method = "spline",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- valleys(x = mspct,
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm,
                       var.name = var.name,
                       refine.wl = refine.wl,
                       method = method,
                       ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

    valleys.idx <-
      which(find_valleys(x[[var.name]],
                         span = span,
                         global.threshold = global.threshold,
                         local.threshold = local.threshold,
                         local.reference = local.reference,
                         threshold.range = threshold.range,
                         strict = strict,
                         na.rm = na.rm))
    if (refine.wl && length(valleys.idx > 0L)) {
      fit_valleys(x = x,
                  valleys.idx = valleys.idx,
                  span = span,
                  y.col.name = var.name,
                  method = method)
    } else {
      x[valleys.idx,  , drop = FALSE]
    }
  }

#' @rdname valleys
#'
#' @export
#'
valleys.raw_spct <- function(x, span = 5,
                             global.threshold = NULL,
                             local.threshold = NULL,
                             local.reference = "minimum",
                             threshold.range = NULL,
                             strict = FALSE,
                             na.rm = FALSE,
                             var.name = "counts",
                             refine.wl = FALSE,
                             method = "spline",
                             ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(x) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    # call method on the collection
    mspct <- valleys(x = mspct,
                     span = span,
                     global.threshold = global.threshold,
                     local.threshold = local.threshold,
                     local.reference = local.reference,
                     threshold.range = threshold.range,
                     strict = strict,
                     na.rm = na.rm,
                     var.name = var.name,
                     refine.wl = refine.wl,
                     method = method,
                     ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  check_wl_stepsize(x = x, span = span, na.rm = TRUE, min.stepsize = 3)

  valleys.idx <-
    which(find_valleys(x[[var.name]],
                       span = span,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       strict = strict,
                       na.rm = na.rm))
  if (refine.wl && length(valleys.idx > 0L)) {
    fit_valleys(x = x,
                valleys.idx = valleys.idx,
                span = span,
                y.col.name = var.name,
                method = method)
  } else {
    x[valleys.idx,  , drop = FALSE]
  }
}


# _mspct ------------------------------------------------------------------


#' @rdname valleys
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
valleys.generic_mspct <- function(x,
                                  span = 5,
                                  global.threshold = NULL,
                                  local.threshold = NULL,
                                  local.reference = "minimum",
                                  threshold.range = NULL,
                                  strict = FALSE,
                                  na.rm = FALSE,
                                  var.name = NULL,
                                  refine.wl = FALSE,
                                  method = "spline",
                                  ...,
                                  .parallel = FALSE,
                                  .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(x,
          .fun = valleys,
          span = span,
          global.threshold = global.threshold,
          local.threshold = local.threshold,
          local.reference = local.reference,
          threshold.range = threshold.range,
          strict = strict,
          na.rm = na.rm,
          var.name = var.name,
          refine.wl = refine.wl,
          method = method,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @rdname valleys
#'
#' @export
#'
valleys.source_mspct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = valleys,
            span = span,
            global.threshold = global.threshold,
            local.threshold = local.threshold,
            local.reference = local.reference,
            threshold.range = threshold.range,
            strict = strict,
            unit.out = unit.out,
            na.rm = na.rm,
            refine.wl = refine.wl,
            method = method,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @rdname valleys
#'
#' @export
#'
valleys.response_mspct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = valleys,
            span = span,
            global.threshold = global.threshold,
            local.threshold = local.threshold,
            local.reference = local.reference,
            threshold.range = threshold.range,
            strict = strict,
            unit.out = unit.out,
            na.rm = na.rm,
            refine.wl = refine.wl,
            method = method,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @rdname valleys
#'
#' @export
#'
valleys.filter_mspct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty",
                                  default = "transmittance"),
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = valleys,
            span = span,
            global.threshold = global.threshold,
            local.threshold = local.threshold,
            local.reference = local.reference,
            threshold.range = threshold.range,
            strict = strict,
            filter.qty = filter.qty,
            na.rm = na.rm,
            refine.wl = refine.wl,
            method = method,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }


#' @rdname valleys
#'
#' @export
#'
valleys.reflector_mspct <-
  function(x,
           span = 5,
           global.threshold = NULL,
           local.threshold = NULL,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = FALSE,
           na.rm = FALSE,
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = valleys,
            span = span,
            global.threshold = global.threshold,
            local.threshold = local.threshold,
            local.reference = local.reference,
            threshold.range = threshold.range,
            strict = strict,
            na.rm = na.rm,
            refine.wl = refine.wl,
            method = method,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }


#' @rdname valleys
#'
#' @export
#'
valleys.solute_mspct <- valleys.reflector_mspct

#' @rdname valleys
#'
#' @export
#'
valleys.cps_mspct <- function(x,
                              span = 5,
                              global.threshold = NULL,
                              local.threshold = NULL,
                              local.reference = "minimum",
                              threshold.range = NULL,
                              strict = FALSE,
                              na.rm = FALSE,
                              var.name = "cps",
                              refine.wl = FALSE,
                              method = "spline",
                              ...,
                              .parallel = FALSE,
                              .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(x,
          .fun = valleys,
          span = span,
          global.threshold = global.threshold,
          local.threshold = local.threshold,
          local.reference = local.reference,
          threshold.range = threshold.range,
          strict = strict,
          na.rm = na.rm,
          var.name = var.name,
          refine.wl = refine.wl,
          method = method,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @rdname valleys
#'
#' @export
#'
valleys.raw_mspct <- function(x,
                              span = 5,
                              global.threshold = NULL,
                              local.threshold = NULL,
                              local.reference = "minimum",
                              threshold.range = NULL,
                              strict = FALSE,
                              na.rm = FALSE,
                              var.name = "counts",
                              refine.wl = FALSE,
                              method = "spline",
                              ...,
                              .parallel = FALSE,
                              .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(x,
          .fun = valleys,
          span = span,
          global.threshold = global.threshold,
          local.threshold = local.threshold,
          local.reference = local.reference,
          threshold.range = threshold.range,
          strict = strict,
          na.rm = na.rm,
          var.name = var.name,
          refine.wl = refine.wl,
          method = method,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

# find wavelengths for a target y ----------------------------------------------

#' Find wavelength values in a spectrum
#'
#' Find wavelength values corresponding to a target y value in any spectrum. The
#' name of the column of the spectral data to be used to match the target needs
#' to be passed as argument unless the spectrum contains a single numerical
#' variable in addition to "w.length".
#'
#' @param x an R object
#' @param target numeric or character. A numeric value indicates the spectral
#'   quantity value for which wavelengths are to be searched. A character
#'   representing a number is converted to a number. A character value
#'   representing a number followed by a function name, will be also accepted
#'   and decoded, such that \code{"0.1max"} is interpreted as targetting one
#'   tenthof the maximum value in a column. The character strings "half.maximum"
#'   and "HM" are synonyms for "0.5max" while "half.range" and "HR" are synonyms
#'   for "0.5range". These synonyms are converted to the cannonical form before
#'   saving them to the returned value.
#' @param col.name.x character The name of the column in which to the
#'   independent variable is stored. Defaults to "w.length" for objects of class
#'   \code{"generic_spct"} or derived.
#' @param col.name character The name of the column in which to search for the
#'   target value.
#' @param .fun function A binary comparison function or operator.
#' @param interpolate logical Indicating whether the nearest wavelength value in
#'   \code{x} should be returned or a value calculated by linear interpolation
#'   between wavelength values stradling the target.
#' @param idfactor logical or character Generates an index column of factor
#'   type. If \code{idfactor = TRUE} then the column is auto named target.idx.
#'   Alternatively the column name can be directly passed as argument to
#'   \code{idfactor} as a character string.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for the target.
#'
#' @note This function is used internally by method \code{wls_at_target()}, and
#'   these methods should be preferred in user code and scripts.
#'
#' @return A spectrum object of the same class as \code{x} with fewer rows,
#'   possibly even no rows. If \code{FALSE} is passed to \code{interpolate} a
#'   subset of \code{x} is returned, otherwise a new object of the same class
#'   containing interpolated wavelenths for the \code{target} value is returned.
#'
#' @examples
#' find_wls(white_led.source_spct)
#' find_wls(white_led.source_spct, target = "0.5max")
#' find_wls(white_led.source_spct, target = 0.4)
#' find_wls(white_led.source_spct, target = 0.4, interpolate = TRUE)
#' find_wls(white_led.source_spct, target = c(0.3, 0.4))
#' find_wls(white_led.source_spct, target = c(0.3, 0.4), idfactor = "target")
#' find_wls(white_led.source_spct, target = c(0.3, 0.4), idfactor = TRUE)
#' find_wls(white_led.source_spct, target = "0.5max")
#' find_wls(white_led.source_spct, target = "0.05max")
#' find_wls(white_led.source_spct, target = "0.5range")
#'
#' led.df <- as.data.frame(white_led.source_spct)
#' find_wls(led.df)
#' find_wls(led.df, col.name = "s.e.irrad", col.name.x = "w.length")
#' find_wls(led.df, col.name = "s.e.irrad", col.name.x = "w.length",
#'          target = 0.4)
#' find_wls(led.df, col.name = "s.e.irrad", col.name.x = "w.length",
#'          target = c(0.3, 0.4))
#' find_wls(led.df, col.name = "s.e.irrad", col.name.x = "w.length",
#'          target = 0.4, idfactor = "target")
#'
#' @export
#'
find_wls <- function(x,
                     target = NULL,
                     col.name.x = NULL,
                     col.name = NULL,
                     .fun = `<=`,
                     interpolate = FALSE,
                     idfactor = length(target) > 1,
                     na.rm = FALSE) {
  stopifnot(is.data.frame(x))
  x.class <- class(x)[1]
  target <- na.omit(target)
  if (!length(target)) {
    return(x[NULL, ])
  }

  if (is.null(col.name.x)) {
    if (is.any_spct(x)) {
      col.name.x <- "w.length"
    } else {
      warning("Object is not a \"generic spectrum\" explicit argument to 'col.name.x' required.")
      return(x[NULL, ])
    }
  }
  if (is.null(col.name)) {
    # find target variable
    col.name <- names(x)
    col.name <- subset(col.name, sapply(x, is.numeric))
    col.name <- setdiff(col.name, col.name.x)
    if (length(col.name) > 1L) {
      warning("Multiple numeric data columns found, explicit argument to 'col.name' required.")
      return(x[NULL, ])
    }
  }
  # This will make a mess of idxs
  # if (na.rm) {
  #   x <- na.omit(x)
  # }
  # .fun may not be vectorized over targets so we need to iterate
  collector.ls <- list()
  target.names <- character()
  targets <- target
  for (target in targets) {
    # keeping this inside the loop allows target's argument can be a list
    if (is.character(target)) {
      target <- gsub("half.maximum|HM", "0.5max", target)
      target <- gsub("half.range|HR", "0.5range", target)
      target <- gsub("[ ]*", "", target)
      ref.fun.name <- gsub("^[0-9.]*", "", target)
      num.in.target <- as.numeric(gsub(ref.fun.name, "", target))
      if (ref.fun.name == "") { # target is absolute
        target.num <- num.in.target
      } else { # target is relative
        ref.fun <- match.fun(ref.fun.name)
        if (na.rm) {
          ref.num <- ref.fun(na.omit(x[[col.name]]))
        } else {
          ref.num <- ref.fun(x[[col.name]])
        }
        if (length(ref.num) == 1L) { # summary value
          target.num <- ref.num * num.in.target
        } else if (length(ref.num) == 2L) { # two sorted values as range
          target.num <- ref.num[1] + (ref.num[2] - ref.num[1]) * num.in.target
        } else {
          warning("Target function '", ref.fun.name, "' returned value of length > 2")
          target.num <- NA_real_
        }
      }
    } else if (is.numeric(target)) {
      target.num <- target
    }
    # test all rows for the condition
    true.rows <- .fun(x[[col.name]], target.num)
    # use run length to find transition points
    runs <- rle(true.rows)
    if (length(runs[["lengths"]]) < 2) {
      next()
      #    return(do.call(x.class, args = list()))
    }
    # accumulate run lengths to get index positions
    opening.idx <- cumsum(runs[["lengths"]][-length(runs[["lengths"]])])
    closing.idx <- opening.idx + 1L
    if (max(closing.idx, na.rm = na.rm) > nrow(x)) {
      closing.idx[length(closing.idx)] <- nrow(x)
    }
    if (interpolate) {
      # do vectorized interpolation to fetch true intersects
      delta.wl <- x[[col.name.x]][closing.idx] - x[[col.name.x]][opening.idx]
      delta.col <- x[[col.name]][closing.idx] - x[[col.name]][opening.idx]
      delta.col.target <- target.num - x[[col.name]][opening.idx]
      wl.increment <- delta.wl * abs(delta.col.target / delta.col)
      wls <- x[[col.name.x]][opening.idx] + wl.increment

      # return as a "short" spectrum containing only matching wls and target values
      z <- tibble::tibble(wls, target = target.num)
      names(z) <- c(col.name.x, col.name)
      if (x.class %in% spct_classes()) {
        z <- do.call(paste("as", x.class, sep = "."), args = list(x = z))
        # we need to copy our private attributes as we are building a new object
        z <- copy_attributes(x, z)
      }
    } else {
      # extract nearest wl value for target
      idxs <- ifelse(abs((x[[col.name]][closing.idx] - target.num) /
                           (x[[col.name]][closing.idx] - x[[col.name]][opening.idx])) > 0.5,
                     opening.idx,
                     closing.idx)
      # if the target value is close to a peak or valley, we may pick the same idx on both sides of it.
      z <- x[unique(idxs), ]
    }
    current.name <- if (is.character(target)) target else as.character(target)
    target.names <- c(target.names, rep(current.name, nrow(z)))
    collector.ls <- c(collector.ls, list(z))
  }
  if (!length(collector.ls)) {
    # we will still bind it to respect idfactor and ensure invariant columns
    collector.ls[[1L]] <- x[NULL, ]
  }
  if (is.any_spct(x)) {
    zz <- rbindspct(collector.ls,
                    fill = FALSE,
                    idfactor = FALSE,
                    attrs.source = 1L)
  } else { # data.frame or tibble
    zz <- dplyr::bind_rows(collector.ls, .id = NULL)
  }
  if (is.logical(idfactor) && idfactor) {
    idfactor <- "target.idx"
  }
  if (is.character(idfactor)) {
    zz[[idfactor]] <- factor(target.names)
  }
  # order by wavelength for spectra
  zz[order(zz[[col.name.x]]), ]
}

# find wavelengths for a target y ----------------------------------------------

#' Find wavelengths values corresponding to a target spectral value
#'
#' Find wavelength values corresponding to a target spectral value in a spectrum.
#' The name of the column of the spectral data to be used is inferred from the
#' class of \code{x} and the argument passed to \code{unit.out} or
#' \code{filter.qty} or their defaults that depend on R options set.
#'
#' @param x data.frame or spectrum object.
#' @param target numeric or character vector. A numeric value indicates the spectral
#'   quantity value for which wavelengths are to be searched. A character
#'   string representing a number is converted to numeric. A character value
#'   representing a number followed by a function name, will be also accepted
#'   and decoded, such that \code{"0.1max"} is interpreted as targeting one
#'   tenth of the maximum value in the column. The character
#'   strings "half.maximum" and "HM" are synonyms for "0.5max" while
#'   "half.range" and "HR" are synonyms for "0.5range".
#' @param interpolate logical Indicating whether the nearest wavelength value
#'   in \code{x} should be returned or a value calculated by linear
#'   interpolation between wavelength values straddling the target.
#' @param x.var.name,y.var.name,col.name character The name of the columns in
#'   which to search for the target value. Use of \code{col.name} is deprecated,
#'   and is a synonym for \code{y.var.name}.
#' @param idfactor logical or character Generates an index column of factor
#'   type. If \code{idfactor = TRUE} then the column is auto named target.idx.
#'   Alternatively the column name can be directly passed as argument to
#'   \code{idfactor} as a character string.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for the target.
#' @param ... currently ignored.
#'
#' @return A data.frame, a spectrum object or a collection of spectra object of
#'   the same class as \code{x} with fewer rows, possibly even no rows. If
#'   \code{FALSE} is passed to \code{interpolate} a subset of \code{x} is
#'   returned, otherwise a new object of the same class containing interpolated
#'   wavelengths for the \code{target} value is returned. As `target` accepts
#'   a vector or list as argument, a factor can be added to the output with
#'   the corresponding target value.
#'
#' @note When interpolation is used, only column \code{w.length} and the column
#'   against which the target value was compared are included in the returned
#'   object, otherwise, all columns in \code{x} are returned. We implement
#'   support for \code{data.frame} to simplify the coding of 'ggplot2' stats
#'   using this function.
#'
#' @examples
#' wls_at_target(sun.spct, target = 0.1)
#' wls_at_target(sun.spct, target = 2e-6, unit.out = "photon")
#' wls_at_target(polyester.spct, target = "HM")
#' wls_at_target(polyester.spct, target = "HM", interpolate = TRUE)
#' wls_at_target(polyester.spct, target = "HM", idfactor = "target")
#' wls_at_target(polyester.spct, target = "HM", filter.qty = "absorbance")
#'
#' @export
#'
#' @family peaks and valleys functions
#'
wls_at_target <- function(x,
                          target,
                          interpolate,
                          idfactor,
                          na.rm,
                          ...) UseMethod("wls_at_target")

#' @rdname wls_at_target
#'
#' @export
#'
wls_at_target.default <-
  function(x,
           target = NULL,
           interpolate = FALSE,
           idfactor = length(target) > 1,
           na.rm = FALSE,
           ...) {
    warning("Method 'wls_at_target' not implemented for objects of class ", class(x)[1])
    x[NULL]
  }

#' @rdname wls_at_target
#'
#' @export
#'
wls_at_target.data.frame <-
  function(x,
           target = "0.5max",
           interpolate = FALSE,
           idfactor = length(target) > 1,
           na.rm = FALSE,
           x.var.name = NULL,
           y.var.name = NULL,
           ...) {
    if (is.null(y.var.name) || is.null(x.var.name)) {
      warning("Variable (column) names required.")
      return(x[NA, ])
    }
    find_wls(x,
             target = target,
             col.name.x = x.var.name,
             col.name = y.var.name,
             .fun = `<=`,
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm)
  }

#' @rdname wls_at_target
#'
#' @export
#'
wls_at_target.generic_spct <-
  function(x,
           target = "0.5max",
           interpolate = FALSE,
           idfactor = length(target) > 1,
           na.rm = FALSE,
           col.name = NULL,
           y.var.name = col.name,
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- wls_at_target(x = mspct,
                             target = target,
                             interpolate = interpolate,
                             idfactor = idfactor,
                             na.rm = na.rm,
                             col.name = col.name,
                             y.var.name = y.var.name,
                             ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    find_wls(x,
             target = target,
             col.name = col.name,
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm,
             ...)
  }

#' @rdname wls_at_target
#'
#' @param unit.out character One of "energy" or "photon"
#'
#' @export
#'
wls_at_target.source_spct <-
  function(x,
           target = "0.5max",
           interpolate = FALSE,
           idfactor = length(target) > 1,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- wls_at_target(x = mspct,
                             target = target,
                             interpolate = interpolate,
                             idfactor = idfactor,
                             na.rm = na.rm,
                             unit.out = unit.out,
                             ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.irrad"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.irrad"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    find_wls(z,
             target = target,
             col.name = col.name,
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm)
  }

#' @rdname wls_at_target
#'
#' @export
#'
wls_at_target.response_spct <-
  function(x,
           target = "0.5max",
           interpolate = FALSE,
           idfactor = length(target) > 1,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- wls_at_target(x = mspct,
                             target = target,
                             interpolate = interpolate,
                             idfactor = idfactor,
                             na.rm = na.rm,
                             unit.out = unit.out,
                             ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.response"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.response"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    find_wls(z,
             target = target,
             col.name = col.name,
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm)
  }

#' @rdname wls_at_target
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
wls_at_target.filter_spct <-
  function(x,
           target = "0.5max",
           interpolate = FALSE,
           idfactor = length(target) > 1,
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty",
                                  default = "transmittance"),
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- wls_at_target(x = mspct,
                             target = target,
                             interpolate = interpolate,
                             idfactor = idfactor,
                             na.rm = na.rm,
                             filter.qty = filter.qty,
                             ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    if (filter.qty == "transmittance") {
      z <- A2T(x, "replace", FALSE)
      col.name <- "Tfr"
    } else if (filter.qty == "absorbance") {
      z <- T2A(x, "replace", FALSE)
      col.name <- "A"
    } else {
      stop("Unrecognized 'filter.qty': ", filter.qty)
    }
    find_wls(z,
             target = target,
             col.name = col.name,
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm)
  }

#' @rdname wls_at_target
#'
#' @export
#'
wls_at_target.reflector_spct <-
  function(x,
           target = "0.5max",
           interpolate = FALSE,
           idfactor = length(target) > 1,
           na.rm = FALSE,
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- wls_at_target(x = mspct,
                             target = target,
                             interpolate = interpolate,
                             idfactor = idfactor,
                             na.rm = na.rm,
                             ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    find_wls(x,
             target = target,
             col.name = "Rfr",
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm)
  }

#' @rdname wls_at_target
#'
#' @export
#'
wls_at_target.solute_spct <-
  function(x,
           target = "0.5max",
           interpolate = FALSE,
           idfactor = length(target) > 1,
           na.rm = FALSE,
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- wls_at_target(x = mspct,
                             target = target,
                             interpolate = interpolate,
                             idfactor = idfactor,
                             na.rm = na.rm,
                             ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    cols <- intersect(c("K.mole", "K.mass"), names(x))
    if (length(cols) == 1) {
      col.name <- cols
      z <- x
    } else {
      stop("Invalid number of columns found:", length(cols))
    }
    find_wls(z,
             target = target,
             col.name = col.name,
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm)
  }

#' @rdname wls_at_target
#'
#' @export
#'
wls_at_target.cps_spct <-
  function(x,
           target = "0.5max",
           interpolate = FALSE,
           idfactor = length(target) > 1,
           na.rm = FALSE,
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- wls_at_target(x = mspct,
                             target = target,
                             interpolate = interpolate,
                             idfactor = idfactor,
                             na.rm = na.rm,
                             ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    find_wls(x,
             target = target,
             col.name = "cps",
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm)
  }


#' @rdname wls_at_target
#'
#' @export
#'
wls_at_target.raw_spct <-
  function(x,
           target = "0.5max",
           interpolate = FALSE,
           idfactor = length(target) > 1,
           na.rm = FALSE,
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- wls_at_target(x = mspct,
                             target = target,
                             interpolate = interpolate,
                             idfactor = idfactor,
                             na.rm = na.rm,
                             ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    find_wls(x,
             target = target,
             col.name = "counts",
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm)
  }


# _mspct ------------------------------------------------------------------


#' @rdname wls_at_target
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
wls_at_target.generic_mspct <- function(x,
                                        target = "0.5max",
                                        interpolate = FALSE,
                                        idfactor = length(target) > 1,
                                        na.rm = FALSE,
                                        ...,
                                        .parallel = FALSE,
                                        .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(x,
          .fun = wls_at_target,
          target = target,
          interpolate = interpolate,
          idfactor = idfactor,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}
