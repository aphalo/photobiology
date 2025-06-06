#' Find local maxima or global maximum (peaks)
#'
#' This function finds peaks (local maxima) in a numeric vector, using a user
#' selectable span and global and local size thresholds, returning a
#' \code{logical} vector.
#'
#' @param x numeric vector. Hint: to find valleys, change the sign of the
#'   argument with the unary operator \code{-}.
#' @param global.threshold numeric A value between 0.0 and 1.0,
#'   relative to \code{threshold.range} indicating the \emph{global} height
#'   (depth) threshold below which peaks (valleys) will be ignored, or a
#'   negative value, between 0.0 and -1.0 indicating the \emph{global} height
#'   (depth) threshold above which peaks (valleys) will be ignored. If
#'   \code{threshold.range = 0} or the value passed as argument belongs to class
#'   \code{"AsIs"} the value is interpreted as an absolute value expressed in
#'   data units.
#' @param local.threshold numeric A value between 0.0 and 1.0, relative to
#'   \code{threshold.range}, indicating the
#'   \emph{within-window} height (depth) threshold below which peaks (valleys)
#'   will be ignored.  If \code{threshold.range = 0} or the value passed
#'   as argument belongs to class \code{"AsIs"} the value is interpreted as an
#'   absolute value expressed in data units.
#' @param local.reference character One of \code{"minimum"} (eqv.
#'   \code{"maximum"}) or \code{"median"}. The reference used to assess the
#'   height of the peak, either the minimum (maximum) value within the window or
#'   the median of all values in the window.
#' @param threshold.range numeric vector of length 2 or a longer vector or list
#'   on which a call to \code{range()} returns a numeric vector of length 2. If
#'   \code{NULL}, the default, \code{range(x)} is used.
#' @param span odd integer A peak is defined as an element in a sequence which
#'   is greater than all other elements within a moving window of width
#'   \code{span} centred at that element. The default value is 5, meaning that a
#'   peak is taller than its four nearest neighbours. \code{span = NULL} extends
#'   the span to the whole length of \code{x}.
#' @param strict logical flag: if TRUE, an element must be strictly greater than
#'   all other values in its window to be considered a peak. Default: TRUE.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for peaks.
#'
#' @return A vector of logical values of the same length as \code{x}. Values
#'   that are TRUE correspond to local peaks in vector \code{x} and can be used
#'   to extract the rows corresponding to peaks from a data frame.
#'
#' @details Function \code{find_peaks} is a wrapper built onto function
#'   \code{\link[splus2R]{peaks}} from \pkg{splus2R}, adds support for peak
#'   height thresholds and handles \code{span = NULL} and non-finite (including
#'   NA) values differently than \code{splus2R::peaks}. Instead of giving an
#'   error when \code{na.rm = FALSE} and \code{x} contains \code{NA} values,
#'   \code{NA} values are replaced with the smallest finite value in \code{x}.
#'   \code{span = NULL} is treated as a special case and selects \code{max(x)}.
#'   Two tests are optional, one based on the absolute height of the peaks
#'   (\code{global.threshold}) and another based on the height of the peaks
#'   compared to other values within the window of width equal to \code{span}
#'   (\code{local.threshold}). The reference value used within each window
#'   containing a peak is given by \code{local.reference}. Parameter
#'   \code{threshold.range} determines how the values passed as argument to
#'   \code{global.threshold} and \code{local.threshold} are scaled. The default,
#'   \code{NULL} uses the range of \code{x}. Thresholds for ignoring too small
#'   peaks are applied after peaks are searched for, and negative threshold
#'   values can in some cases result in no peaks being returned.
#'
#'   While function \code{find_peaks} accepts as input a \code{numeric} vector
#'   and returns a \code{logical} vector, methods \code{peaks} and
#'   \code{valleys} accept as input different R objects, including spectra and
#'   collections of spectra and return a subset of the object. These methods
#'   are implemented using calls to functions \code{find_peaks} and
#'   \code{\link{fit_peaks}}.
#'
#' @note The default for parameter \code{strict} is \code{FALSE} in functions
#'   \code{peaks()} and \code{find_peaks()}, as in \code{stat_peaks()} and in
#'   \code{stat_valleys()}, while the default in \code{\link[splus2R]{peaks}}
#'   is \code{strict = TRUE}.
#'
#' @seealso \code{\link[splus2R]{peaks}}.
#'
#' @family peaks and valleys functions
#'
#' @export
#'
#' @examples
#' with(sun.data, w.length[find_peaks(s.e.irrad)])
#'
find_peaks <-
  function(x,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           span = 3,
           strict = FALSE,
           na.rm = FALSE) {
    # keep track
    threshold.delta.computed <- FALSE
    # find peaks
    if(is.null(span) || span >= length(x)) {
      pks <- x == max(x, na.rm = na.rm)
      if (strict && sum(pks) != 1L) {
        pks <- logical(length(x)) # all FALSE
      }
    } else {
      pks <- splus2R::peaks(x = x, span = span, strict = strict)
    }

    x <- ifelse(!is.finite(x), min(x, na.rm = TRUE), x)
    # discard peaks that are low or not locally prominent
    if (abs(global.threshold) >= 1e-5 || abs(local.threshold) >= 1e-5) {
      if (is.null(threshold.range)) {
        threshold.range <- range(x)
      } else {
        threshold.range <- range(threshold.range)
      }
      if (all(abs(threshold.range < 1e-5))) {
        if (!inherits(global.threshold, "AsIs")) {
          global.threshold <- I(global.threshold)
        }
      } else {
        threshold.delta <- threshold.range[2] - threshold.range[1]
        threshold.delta.computed <- TRUE
        if (length(unique(threshold.range)) != 2L ||
            !all(is.finite(threshold.range))) {
          threshold.range <- signif(threshold.range, digits = 3)
          stop("Bad 'threshold.range' value: (",
               threshold.range[1], ", ", threshold.range[2], ")")
        }
      }
      # apply global height threshold test to found peaks
      if (abs(global.threshold) >= 1e-5) {
        # this can cater for the case when max_x < 0, as with logs
        if (inherits(global.threshold, "AsIs")) {
          pks <- ifelse(x > global.threshold, pks , FALSE)
        } else {
          scaled.global.threshold <- threshold.delta * abs(global.threshold)
          if (global.threshold > 0.0) {
            pks <- ifelse(x - threshold.range[1] > scaled.global.threshold,
                          pks ,
                          FALSE)
          } else {
            pks <- ifelse(x - threshold.range[1] <= scaled.global.threshold,
                          pks ,
                          FALSE)
          }
        }
      }
      # apply local.threshold height test to found peaks
      if (local.threshold > 0) {
        # we always search for maxima, even in the case of valleys
        if (local.reference == "maximum") {
          local.reference <- "minimum"
        }

        if (all(threshold.range == 0) && !inherits(local.threshold, "AsIs")) {
          local.threshold <- I(local.threshold)
        } else if (!threshold.delta.computed) {
          threshold.delta <- threshold.range[2] - threshold.range[1]
          if (length(unique(threshold.range)) != 2L ||
              !all(is.finite(threshold.range))) {
            threshold.range <- signif(threshold.range, digits = 3)
            stop("Bad 'threshold.range' value: (",
                 threshold.range[1], ", ",
                 threshold.range[2], ")")
          }
        }
        # apply local height threshold test to found peaks
        if (abs(local.threshold) >= 1e-5) {
          smooth_x <-
            switch(local.reference,
                   minimum = caTools::runmin(x, k = span, endrule = "min"),
                   median = stats::runmed(x, k = span, endrule = "median"))

          if (inherits(local.threshold, "AsIs")) {
            pks <- ifelse(x - smooth_x > local.threshold, pks , FALSE)
          } else {
            scaled.local.threshold <- threshold.delta * abs(local.threshold)
            pks <- ifelse(x - smooth_x > scaled.local.threshold, pks , FALSE)
          }
        }
      }
    }
    pks
  }

#' Get peaks and valleys from a spectrum
#'
#' These functions "get" (or extract) peaks (maxima) and valleys (minima) in two
#' vectors, usually a spectral quantity and wavelength, using a user selectable
#' span for window width and global and local (within moving window) size
#' thresholds. They also generate \code{character} values for \code{x}.
#'
#' @param x,y numeric
#' @param x_unit character Vector of texts to be pasted at end of labels built
#'   from x value at peaks.
#' @param x_digits numeric Number of significant digits in wavelength label.
#' @inheritParams find_peaks
#'
#' @inherit find_peaks details seealso
#'
#' @note The use of these two functions is deprecated. They are retained for
#'   backwards compatibility and will be removed in the near future.
#'
#' @return A data frame with variables w.length and s.irrad with their values at
#'   the peaks or valleys plus a character variable of labels.
#'
#' @export
#'
#' @family peaks and valleys functions
#'
get_peaks <- function(x,
                      y,
                      global.threshold = 0,
                      local.threshold = 0,
                      local.reference = "minimum",
                      threshold.range = NULL,
                      span = 5,
                      strict = TRUE,
                      x_unit = "",
                      x_digits = 3,
                      na.rm = FALSE) {
  warning("Functions 'get_peaks()' and 'get_valeys()' have been deprecated: please use 'peaks()' and 'valleys()', or 'find_peaks()', instead.")
  stopifnot(length(x) == length(y))
  selector <- find_peaks(x = y,
                         global.threshold = global.threshold,
                         local.threshold = local.threshold,
                         local.reference = local.reference,
                         threshold.range = threshold.range,
                         span = span,
                         strict = strict,
                         na.rm = na.rm)
  if (sum(selector) < 1) {
    return(data.frame(
      x = numeric(0),
      y = numeric(0),
      label = character(0)
    ))
  } else {
    peaks.x <- x[selector]
    peaks.y <- y[selector]
    return(data.frame(
      x = peaks.x,
      y = peaks.y,
      label = paste(as.character(signif(
        x = peaks.x, digits = x_digits
      )), x_unit, sep = "")
    ))
  }
}

#' @rdname get_peaks
#' @export
#'
get_valleys <- function(x, y,
                        global.threshold = 0,
                        local.threshold = 0,
                        local.reference = "minimum",
                        threshold.range = NULL,
                        span = 5,
                        strict = TRUE,
                        x_unit = "",
                        x_digits = 3,
                        na.rm = FALSE) {
  xy.data <- get_peaks(x = x, y = -y,
                       global.threshold = global.threshold,
                       local.threshold = local.threshold,
                       local.reference = local.reference,
                       threshold.range = threshold.range,
                       span = span,
                       strict = strict,
                       x_unit = x_unit,
                       x_digits = x_digits,
                       na.rm = na.rm)
  xy.data[["y"]] <- -xy.data[["y"]]
  return(xy.data)
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
#' peaks(sun.spct, local.threshold = 0.05, local.reference = "minimum")
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

#' @describeIn peaks Default returning always NA.
#' @export
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

#' @describeIn peaks Default function usable on numeric vectors.
#' @export
peaks.numeric <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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

#' @describeIn peaks  Method for "data.frame" objects.
#'
#' @export
#'
peaks.data.frame <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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

#' @describeIn peaks  Method for "generic_spct" objects.
#'
#' @export
#'
peaks.generic_spct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span)

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

#' @describeIn peaks  Method for "source_spct" objects.
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
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span)

    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.irrad"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
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

#' @describeIn peaks  Method for "response_spct" objects.
#'
#' @export
#'
peaks.response_spct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span)

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

#' @describeIn peaks  Method for "filter_spct" objects.
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
peaks.filter_spct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span)

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

#' @describeIn peaks  Method for "reflector_spct" objects.
#'
#' @export
#'
peaks.reflector_spct <- function(x,
                                 span = 5,
                                 global.threshold = 0,
                                 local.threshold = 0,
                                 local.reference = "minimum",
                                 threshold.range = NULL,
                                 strict = TRUE,
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
    return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
  }

  check_wl_stepsize(x = x, span = span)

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

#' @describeIn peaks  Method for "solute_spct" objects.
#'
#' @export
#'
peaks.solute_spct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span)

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

#' @describeIn peaks  Method for "cps_spct" objects.
#'
#' @export
#'
peaks.cps_spct <- function(x,
                           span = 5,
                           global.threshold = 0,
                           local.threshold = 0,
                           local.reference = "minimum",
                           threshold.range = NULL,
                           strict = TRUE,
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
    return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
  }

  check_wl_stepsize(x = x, span = span)

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

#' @describeIn peaks  Method for "raw_spct" objects.
#'
#' @export
#'
peaks.raw_spct <- function(x, span = 5,
                           global.threshold = 0,
                           local.threshold = 0,
                           local.reference = "minimum",
                           threshold.range = NULL,
                           strict = TRUE,
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
    return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
  }

  check_wl_stepsize(x = x, span = span)

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


#' @describeIn peaks  Method for "generic_mspct" objects.
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
                                global.threshold = 0,
                                local.threshold = 0,
                                local.reference = "minimum",
                                threshold.range = NULL,
                                strict = TRUE,
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

#' @describeIn peaks  Method for "source_mspct" objects.
#'
#' @export
#'
peaks.source_mspct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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

#' @describeIn peaks  Method for "cps_mspct" objects.
#'
#' @export
#'
peaks.response_mspct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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

#' @describeIn peaks  Method for "filter_mspct" objects.
#'
#' @export
#'
peaks.filter_mspct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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


#' @describeIn peaks  Method for "reflector_mspct" objects.
#'
#' @export
#'
peaks.reflector_mspct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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


#' @describeIn peaks  Method for "solute_mspct" objects.
#'
#' @export
#'
peaks.solute_mspct <- peaks.reflector_mspct

#' @describeIn peaks  Method for "cps_mspct" objects.
#'
#' @export
#'
peaks.cps_mspct <- function(x,
                            span = 5,
                            global.threshold = 0,
                            local.threshold = 0,
                            local.reference = "minimum",
                            threshold.range = NULL,
                            strict = TRUE,
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

#' @describeIn peaks  Method for "raw_mspct" objects.
#'
#' @export
#'
peaks.raw_mspct <- function(x,
                            span = 5,
                            global.threshold = 0,
                            local.threshold = 0,
                            local.reference = "minimum",
                            threshold.range = NULL,
                            strict = TRUE,
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
#' valleys(sun.spct, local.threshold = 0.1, local.reference = "minimum")
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

#' @describeIn valleys Default returning always NA.
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

#' @describeIn valleys Default function usable on numeric vectors.
#' @export
valleys.numeric <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
           na.rm = FALSE,
           ...) {
    x[find_peaks(x = -x,
                 span = span,
                 global.threshold = global.threshold,
                 local.threshold = local.threshold,
                 local.reference = local.reference,
                 threshold.range = threshold.range,
                 strict = strict,
                 na.rm = na.rm)]
  }

#' @describeIn valleys  Method for "data.frame" objects.
#'
#' @export
#'
valleys.data.frame <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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
      which(find_peaks(-x[[var.name]],
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

#' @describeIn valleys  Method for "generic_spct" objects.
#'
#' @export
#'
valleys.generic_spct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span)

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
      which(find_peaks(-x[[var.name]],
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

#' @describeIn valleys  Method for "source_spct" objects.
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
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span)

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
      which(find_peaks(-z[[col.name]],
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

#' @describeIn valleys  Method for "response_spct" objects.
#'
#' @export
#'
valleys.response_spct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span)

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
      which(find_peaks(-z[[col.name]],
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

#' @describeIn valleys  Method for "filter_spct" objects.
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
valleys.filter_spct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span)

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
      which(find_peaks(-z[[col.name]],
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

#' @describeIn valleys  Method for "reflector_spct".
#'
#' @export
#'
valleys.reflector_spct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span)

    col.name <- "Rfr"
    valleys.idx <-
      which(find_peaks(-x[[col.name]],
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

#' @describeIn valleys  Method for "solute_spct" objects.
#'
#' @export
#'
valleys.solute_spct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span)

    cols <- intersect(c("K.mole", "K.mass"), names(x))
    if (length(cols) == 1) {
      col.name <- cols
      z <- x
    } else {
      stop("Invalid number of columns found:", length(cols))
    }
    valleys.idx <-
      which(find_peaks(-z[[col.name]],
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

#' @describeIn valleys  Method for "cps_spct" objects.
#'
#' @export
#'
valleys.cps_spct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    check_wl_stepsize(x = x, span = span)

    valleys.idx <-
      which(find_peaks(-x[[var.name]],
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

#' @describeIn valleys  Method for "raw_spct" objects.
#'
#' @export
#'
valleys.raw_spct <- function(x, span = 5,
                             global.threshold = 0,
                             local.threshold = 0,
                             local.reference = "minimum",
                             threshold.range = NULL,
                             strict = TRUE,
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
    return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
  }

  check_wl_stepsize(x = x, span = span)

  valleys.idx <-
    which(find_peaks(-x[[var.name]],
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


#' @describeIn valleys  Method for "generic_mspct" objects.
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
                                  global.threshold = 0,
                                  local.threshold = 0,
                                  local.reference = "minimum",
                                  threshold.range = NULL,
                                  strict = TRUE,
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

#' @describeIn valleys  Method for "source_mspct" objects.
#'
#' @export
#'
valleys.source_mspct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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

#' @describeIn valleys  Method for "cps_mspct" objects.
#'
#' @export
#'
valleys.response_mspct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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

#' @describeIn valleys  Method for "filter_mspct" objects.
#'
#' @export
#'
valleys.filter_mspct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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


#' @describeIn valleys  Method for "reflector_mspct" objects.
#'
#' @export
#'
valleys.reflector_mspct <-
  function(x,
           span = 5,
           global.threshold = 0,
           local.threshold = 0,
           local.reference = "minimum",
           threshold.range = NULL,
           strict = TRUE,
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


#' @describeIn valleys  Method for "solute_mspct" objects.
#'
#' @export
#'
valleys.solute_mspct <- valleys.reflector_mspct

#' @describeIn valleys  Method for "cps_mspct" objects.
#'
#' @export
#'
valleys.cps_mspct <- function(x,
                              span = 5,
                              global.threshold = 0,
                              local.threshold = 0,
                              local.reference = "minimum",
                              threshold.range = NULL,
                              strict = TRUE,
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

#' @describeIn valleys  Method for "raw_mspct" objects.
#'
#' @export
#'
valleys.raw_mspct <- function(x,
                              span = 5,
                              global.threshold = 0,
                              local.threshold = 0,
                              local.reference = "minimum",
                              threshold.range = NULL,
                              strict = TRUE,
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

#' @describeIn wls_at_target Default returning always an empty object of the
#'   same class as \code{x}.
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

#' @describeIn wls_at_target  Method for "data.frame" objects.
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

#' @describeIn wls_at_target Method for "generic_spct" objects.
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    find_wls(x,
             target = target,
             col.name = col.name,
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm,
             ...)
  }

#' @describeIn wls_at_target Method for "source_spct" objects.
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
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

#' @describeIn wls_at_target Method for "response_spct" objects.
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
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

#' @describeIn wls_at_target Method for "filter_spct" objects.
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
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

#' @describeIn wls_at_target Method for "reflector_spct" objects.
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    find_wls(x,
             target = target,
             col.name = "Rfr",
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm)
  }

#' @describeIn wls_at_target Method for "solute_spct" objects.
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
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

#' @describeIn wls_at_target Method for "cps_spct" objects.
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    find_wls(x,
             target = target,
             col.name = "cps",
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm)
  }


#' @describeIn wls_at_target Method for "raw_spct" objects.
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
      return(rbindspct(mspct, idfactor = FALSE, attrs.simplify = TRUE))
    }

    find_wls(x,
             target = target,
             col.name = "counts",
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm)
  }


# _mspct ------------------------------------------------------------------


#' @describeIn wls_at_target  Method for "generic_mspct" objects.
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


# utils -------------------------------------------------------------------

#' Check wavelength stepsize
#'
#' @inheritParams peaks
#'
#' @details As the search for peaks uses a window based on a fixed number of
#'   observations at neighbouring wavelengths, if the wavelength step between
#'   observations varies drastically, the window expressed in nanometres of
#'   wavelength becomes very irregular. With the default \code{span = 5} the
#'   search for peaks in most cases still works.
#'
#'   The typical case are spectra returned by \code{thin_wl()}, which retain
#'   the original local maxima, and a reasonably narrow wavelength maximum
#'   step size when using default arguments.
#'
#' @return logical \code{TRUE} if check is passed and otherwise \code{FALSE}
#' with a warning.
#'
#' @keywords internal
#'
check_wl_stepsize <-
  function(x,
           span) {
    if (!is.null(span) && span > 5) {
      step.size.range <- wl_stepsize(x)
      step.size.range[1] <- max(1, step.size.range[1]) # ignore wl steps < 1 nm
      if ((step.size.range[2] / step.size.range[1]) > 2) {
        warning("Peaks cannot be reliably searched for when the wavelength step varies")
        return(invisible(FALSE))
      }
    }
    invisible(TRUE)
  }


