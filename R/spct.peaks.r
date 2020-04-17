#' Find peaks in a spectrum
#'
#' This function finds all peaks (local maxima) in a spectrum, using a user
#' provided size threshold relative to the tallest peak (global maximum) bellow
#' which found peaks are ignored---i.e., not included in the returned value. This
#' is a wrapper built on top of function \code{peaks()} from package 'splus2R'.
#'
#' @param x numeric vector
#' @param ignore_threshold numeric Value between 0.0 and 1.0 indicating the
#'   relative size compared to tallest peak threshold below which peaks will be
#'   ignored. Negative values set a threshold so that the tallest peaks are
#'   ignored, instead of the shortest.
#' @param span integer A peak is defined as an element in a sequence which is
#'   greater than all other elements within a window of width \code{span}
#'   centered at that element. Use \code{NULL} for the global peak.
#' @param strict logical If \code{TRUE}, an element must be strictly greater
#'   than all other values in its window to be considered a peak.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for peaks.
#'
#' @return A logical vector of the same length as \code{x}. Values that are
#'   \code{TRUE} correspond to local peaks in the data.
#'
#' @export
#' @examples
#' with(sun.data, w.length[find_peaks(s.e.irrad)])
#'
#' @note This function is a wrapper built on function
#'   \code{\link[splus2R]{peaks}} from \pkg{splus2R} and handles non-finite
#'   (including \code{NA}) values differently than \code{splus2R::peaks},
#'   instead of giving an error they are replaced with the smallest finite value
#'   in \code{x}.
#'
#' @seealso \code{\link[splus2R]{peaks}}
#'
#' @family peaks and valleys functions
#'
find_peaks <-
  function(x,
           ignore_threshold = 0,
           span = 3,
           strict = TRUE,
           na.rm = FALSE) {
    if (na.rm) {
      x <- na.omit(x)
    }
    # find peaks
    if(is.null(span)) {
      pks <- x == max(x)
      if (strict && sum(pks) != 1L) {
        pks <- logical(length(x))
      }
    } else {
      pks <- splus2R::peaks(x = x, span = span, strict = strict)
    }
    # apply threshold to found peaks
    if (abs(ignore_threshold) < 1e-5) {
      pks
    } else {
      range_x <- range(x, finite = TRUE)
      min_x <- range_x[1]
      max_x <- range_x[2]
      x <- ifelse(!is.finite(x), min_x, x)
      # this can cater for the case when max_x < 0, as with logs
      delta <- max_x - min_x
      top_flag <- ignore_threshold > 0.0
      scaled_threshold <- delta * abs(ignore_threshold)
      if (top_flag) {
        ifelse(x - min_x > scaled_threshold, pks , FALSE)
      } else {
        ifelse(max_x - x > scaled_threshold, pks , FALSE)
      }
    }
  }

#' Get peaks and valleys in a spectrum
#'
#' These functions find peaks (local maxima) or valleys (local minima) in a
#' spectrum, using a user selectable size threshold relative to the tallest peak
#' (global maximum). This a wrapper built on top of function peaks from package
#' splus2R.
#'
#' @param x numeric
#' @param y numeric
#' @param ignore_threshold numeric Value between 0.0 and 1.0 indicating the
#'   relative size compared to tallest peak threshold below which peaks will be
#'   ignored. Negative values set a threshold so that the tallest peaks are
#'   ignored, instead of the shortest.
#' @param span integer A peak is defined as an element in a sequence which is
#'   greater than all other elements within a window of width \code{span}
#'   centered at that element. Use \code{NULL} for the global peak.
#' @param strict logical If \code{TRUE}, an element must be strictly greater
#'   than all other values in its window to be considered a peak.
#' @param x_unit character Vector of texts to be pasted at end of labels built
#'   from x value at peaks.
#' @param x_digits numeric Number of significant digits in wavelength label.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for peaks.
#'
#' @return A data frame with variables w.length and s.irrad with their values at
#'   the peaks or valleys plus a character variable of labels.
#'
#' @export
#' @examples
#' with(sun.spct, get_peaks(w.length, s.e.irrad))
#' with(sun.spct, get_valleys(w.length, s.e.irrad))
#'
#' @family peaks and valleys functions
#'
get_peaks <- function(x,
                      y,
                      ignore_threshold = 0,
                      span = 5,
                      strict = TRUE,
                      x_unit = "",
                      x_digits = 3,
                      na.rm = FALSE) {
  stopifnot(length(x) == length(y))
  selector <- find_peaks(x = y,
                         ignore_threshold = ignore_threshold,
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
                        ignore_threshold = 0,
                        span = 5,
                        strict = TRUE,
                        x_unit = "",
                        x_digits = 3,
                        na.rm = FALSE) {
  xy.data <- get_peaks(x = x, y = -y,
                       ignore_threshold = -ignore_threshold,
                       span = span,
                       strict = strict,
                       x_unit = x_unit,
                       x_digits = x_digits,
                       na.rm = na.rm)
  xy.data$y <- -xy.data$y
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
#'   \code{refine.wl = TRUE} of methods \code{peaks()} and \code{valleys()} instead.
#'
#' @return An R object of the same class as \code{x} containing the fitted
#'   values for the peaks, and optionally the values for at \code{peaks.idx} or
#'   \code{valleys.idx} for other retained columns.
#'
#' @examples
#'
#' peaks <- find_peaks(sun.spct$s.e.irrad, span = 31)
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
#' @param x an R object
#' @param ignore_threshold numeric Value between 0.0 and 1.0 indicating the
#'   relative size compared to tallest peak threshold below which peaks will be
#'   ignored. Negative values set a threshold so that the tallest peaks are
#'   ignored, instead of the shortest.
#' @param span integer A peak is defined as an element in a sequence which is
#'   greater than all other elements within a window of width \code{span}
#'   centered at that element. Use \code{NULL} for the global peak.
#' @param strict logical If \code{TRUE}, an element must be strictly greater
#'   than all other values in its window to be considered a peak.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for peaks.
#' @param var.name,x.var.name,y.var.name character Name of column where to look
#'   for peaks.
#' @param refine.wl logical Flag indicating if peak location should be refined by
#'   fitting a function.
#' @param method character String with the name of a method. Currently only
#'   spline interpolation is implemented.
#' @param ... ignored
#'
#' @return A subset of \code{x} with rows corresponding to local maxima.
#'
#' @note Thresholds for ignoring peaks are applied after peaks are searched for,
#' and negative threshold values can in some cases result in no peaks being
#' returned.
#'
#' @export
#'
#' @examples
#' peaks(sun.spct, span = 51)
#' peaks(sun.spct, span = NULL)
#' peaks(sun.spct, span = 51, refine.wl = TRUE)
#'
#' @family peaks and valleys functions
#'
peaks <- function(x, span, ignore_threshold, strict, na.rm, ...) UseMethod("peaks")

#' @describeIn peaks Default returning always NA.
#' @export
peaks.default <-
  function(x, span = NA, ignore_threshold = NA, strict = NA, na.rm = FALSE, ...) {
  warning("Method 'peaks' not implemented for objects of class ", class(x)[1])
  x[NA]
}

#' @describeIn peaks Default function usable on numeric vectors.
#' @export
peaks.numeric <-
  function(x, span = 5, ignore_threshold = NA, strict = TRUE, na.rm = FALSE, ...) {
  x[find_peaks(x = x, span = span, strict = strict, na.rm = na.rm)]
}

#' @describeIn peaks  Method for "data.frame" objects.
#'
#' @export
#'
peaks.data.frame <-
  function(x,
           span = 5,
           ignore_threshold = 0,
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
                       ignore_threshold = ignore_threshold,
                       strict = strict))
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           var.name = NULL,
           refine.wl = FALSE,
           method = "spline",
           ...) {
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
                       span = span, ignore_threshold = ignore_threshold,
                       strict = strict))
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...) {
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
                       span = span, ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0.0,
           strict = TRUE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...) {
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
                            span = span, ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty",
                                  default = "transmittance"),
           refine.wl = FALSE,
           method = "spline",
           ...) {
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
                            span = span, ignore_threshold = ignore_threshold,
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
                                 ignore_threshold = 0,
                                 strict = TRUE,
                                 na.rm = FALSE,
                                 refine.wl = FALSE,
                                 method = "spline",
                                 ...) {
  col.name <- "Rfr"
  peaks.idx <-
    which(find_peaks(x[[col.name]],
                     span = span, ignore_threshold = ignore_threshold,
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

#' @describeIn peaks  Method for "cps_spct" objects.
#'
#' @export
#'
peaks.cps_spct <- function(x, span = 5,
                           ignore_threshold = 0,
                           strict = TRUE,
                           na.rm = FALSE,
                           var.name = "cps",
                           refine.wl = FALSE,
                           method = "spline",
                           ...) {
  peaks.idx <-
    which(find_peaks(x[[var.name]],
                     span = span, ignore_threshold = ignore_threshold,
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
                           ignore_threshold = 0,
                           strict = TRUE,
                           na.rm = FALSE,
                           var.name = "counts",
                           refine.wl = FALSE,
                           method = "spline",
                           ...) {
  peaks.idx <-
    which(find_peaks(x[[var.name]],
                     span = span, ignore_threshold = ignore_threshold,
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
                                ignore_threshold = 0,
                                strict = TRUE,
                                na.rm = FALSE,
                                var.name = NULL,
                                refine.wl = FALSE,
                                method = "spline",
                                ...,
                                .parallel = FALSE,
                                .paropts = NULL) {
  msmsply(x,
          .fun = peaks,
          span = span,
          ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = peaks,
            span = span,
            ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = peaks,
            span = span,
            ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty",
                                  default = "transmittance"),
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = peaks,
            span = span,
            ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = peaks,
            span = span,
            ignore_threshold = ignore_threshold,
            strict = strict,
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
peaks.cps_mspct <- function(x,
                            span = 5,
                            ignore_threshold = 0,
                            strict = TRUE,
                            na.rm = FALSE,
                            var.name = "cps",
                            refine.wl = FALSE,
                            method = "spline",
                            ...,
                            .parallel = FALSE,
                            .paropts = NULL) {
  msmsply(x,
          .fun = peaks,
          span = span,
          ignore_threshold = ignore_threshold,
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
                            ignore_threshold = 0,
                            strict = TRUE,
                            na.rm = FALSE,
                            var.name = "counts",
                            refine.wl = FALSE,
                            method = "spline",
                            ...,
                            .parallel = FALSE,
                            .paropts = NULL) {
  msmsply(x,
          .fun = peaks,
          span = span,
          ignore_threshold = ignore_threshold,
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
#' @param x an R object
#' @param ignore_threshold numeric Value between 0.0 and 1.0 indicating the
#'   relative size compared to tallest peak threshold below which peaks will be
#'   ignored. Negative values set a threshold so that the tallest peaks are
#'   ignored, instead of the shortest.
#' @param span integer A valley is defined as an element in a sequence which is smaller
#'   than all other elements within a window of width \code{span} centered at that
#'   element. Use \code{NULL} for the global peak.
#' @param strict logical If \code{TRUE}, an element must be strictly greater
#'   than all other values in its window to be considered a peak.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for peaks.
#' @param var.name,x.var.name,y.var.name character Name of column where to look
#'   for valleys.
#' @param refine.wl logical Flag indicating if valley location should be refined by
#'   fitting a function.
#' @param method character String with the name of a method. Currently only
#'   spline interpolation is implemented.
#' @param ... ignored
#'
#' @return A subset of \code{x} with rows corresponding to local minima.
#'
#' @examples
#' valleys(sun.spct, span = 50)
#'
#' @export
#'
#' @family peaks and valleys functions
#'
valleys <- function(x, span, ignore_threshold, strict, ...) UseMethod("valleys")

#' @describeIn valleys Default function usable on numeric vectors.
#' @export
valleys.default <- function(x, span, ignore_threshold, strict, ...) {
  x[NA]
}

#' @describeIn valleys Default returning always NA.
#' @export
valleys.default <-
  function(x, span = NA, ignore_threshold = NA, strict = NA, na.rm = FALSE, ...) {
    warning("Method 'valleys' not implemented for objects of class ", class(x)[1])
    x[NA]
  }

#' @describeIn valleys Default function usable on numeric vectors.
#' @export
valleys.numeric <-
  function(x, span = 5, ignore_threshold, strict = TRUE, na.rm = FALSE, ...) {
    x[find_peaks(x = -x, span = span, strict = strict, na.rm = na.rm)]
  }

#' @describeIn valleys  Method for "data.frame" objects.
#'
#' @export
#'
valleys.data.frame <-
  function(x,
           span = 5,
           ignore_threshold = 0,
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
                       ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           var.name = NULL,
           refine.wl = FALSE,
           method = "spline",
           ...) {
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
                       ignore_threshold = ignore_threshold,
                       strict = strict))
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...) {
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
                       ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...) {
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
                       ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty", default = "transmittance"),
           refine.wl = FALSE,
           method = "spline",
           ...) {
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
                       ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           refine.wl = FALSE,
           method = "spline",
           ...) {
    col.name <- "Rfr"
    valleys.idx <-
      which(find_peaks(-x[[col.name]],
                       span = span,
                       ignore_threshold = ignore_threshold,
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

#' @describeIn valleys  Method for "cps_spct" objects.
#'
#' @export
#'
valleys.cps_spct <-
  function(x,
           span = 5,
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           refine.wl = FALSE,
           method = "spline",
           ...) {
    col.name <- "cps"
    valleys.idx <-
      which(find_peaks(-x[[col.name]],
                       span = span,
                       ignore_threshold = ignore_threshold,
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

#' @describeIn valleys  Method for "raw_spct" objects.
#'
#' @export
#'
valleys.raw_spct <- function(x, span = 5,
                             ignore_threshold = 0,
                             strict = TRUE,
                             na.rm = FALSE,
                             var.name = "counts",
                             refine.wl = FALSE,
                             method = "spline",
                             ...) {
  valleys.idx <-
    which(find_peaks(-x[[var.name]],
                     span = span, ignore_threshold = ignore_threshold,
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
                                  ignore_threshold = 0,
                                  strict = TRUE,
                                  na.rm = FALSE,
                                  var.name = NULL,
                                  refine.wl = FALSE,
                                  method = "spline",
                                  ...,
                                  .parallel = FALSE,
                                  .paropts = NULL) {
  msmsply(x,
          .fun = valleys,
          span = span,
          ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = valleys,
            span = span,
            ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = valleys,
            span = span,
            ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty",
                                  default = "transmittance"),
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = valleys,
            span = span,
            ignore_threshold = ignore_threshold,
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
           ignore_threshold = 0,
           strict = TRUE,
           na.rm = FALSE,
           refine.wl = FALSE,
           method = "spline",
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = valleys,
            span = span,
            ignore_threshold = ignore_threshold,
            strict = strict,
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
valleys.cps_mspct <- function(x,
                              span = 5,
                              ignore_threshold = 0,
                              strict = TRUE,
                              na.rm = FALSE,
                              var.name = "cps",
                              refine.wl = FALSE,
                              method = "spline",
                              ...,
                              .parallel = FALSE,
                              .paropts = NULL) {
  msmsply(x,
          .fun = valleys,
          span = span,
          ignore_threshold = ignore_threshold,
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
                              ignore_threshold = 0,
                              strict = TRUE,
                              na.rm = FALSE,
                              var.name = "counts",
                              refine.wl = FALSE,
                              method = "spline",
                              ...,
                              .parallel = FALSE,
                              .paropts = NULL) {
  msmsply(x,
          .fun = valleys,
          span = span,
          ignore_threshold = ignore_threshold,
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
#' @param target numeric value indicating the spectral quantity value for which
#'   wavelengths are to be searched and interpolated if need. The character
#'   strings "half.maximum" and "half.range" are also accepted as arguments.
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
#'   type. If \code{idfactor = TRUE} then the column is auto named spct.idx.
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
#' find_wls(white_led.source_spct, target = "half.maximum")
#' find_wls(white_led.source_spct, target = 0.4)
#' find_wls(white_led.source_spct, target = 0.4, interpolate = TRUE)
#' find_wls(white_led.source_spct, target = c(0.3, 0.4))
#' find_wls(white_led.source_spct, target = c(0.3, 0.4), idfactor = "target")
#' find_wls(white_led.source_spct, target = c(0.3, 0.4), idfactor = TRUE)
#' find_wls(white_led.source_spct, target = c("HM", "HR"))
#' find_wls(white_led.source_spct, target = c("HM", "HR"), interpolate = TRUE)
#'
#' led.df <- as.data.frame(white_led.source_spct)
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
                     idfactor = FALSE,
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
  if (na.rm) {
    x <- na.omit(x)
  }
  # .fun may not be vectorized over targets so we need to iterate
  collector.ls <- list()
  targets <- target
  for (target in targets) {
    if (is.character(target)) {
      if (target %in% c("half.maximum", "HM")) {
        target <- max(x[[col.name]]) / 2
      } else if (target %in% c("half.range", "HR")) {
        target <- mean(range(x[[col.name]]))
      } else {
        warning("Unrecognized character string: '", target, "' passed to 'target'", sep = "")
        target <- NA_real_
      }
      if (is.na(target)) {
        next()
      }
    }
    # test all rows for the condition
    true.rows <- .fun(x[[col.name]], target)
    # use run length to find transition points
    runs <- rle(true.rows)
    if (length(runs$lengths) < 2) {
      next()
      #    return(do.call(x.class, args = list()))
    }
    # accumulate run lengths to get index positions
    opening.idx <- cumsum(runs$lengths[-length(runs$lengths)])
    closing.idx <- opening.idx + 1L
    if (max(closing.idx) > nrow(x)) {
      closing.idx[length(closing.idx)] <- nrow(x)
    }
    if (interpolate) {
      # do vectorized interpolation to fetch true intersects
      delta.wl <- x[[col.name.x]][closing.idx] - x[[col.name.x]][opening.idx]
      delta.col <- x[[col.name]][closing.idx] - x[[col.name]][opening.idx]
      delta.col.target <- target - x[[col.name]][opening.idx]
      wl.increment <- delta.wl * abs(delta.col.target / delta.col)
      wls <- x[[col.name.x]][opening.idx] + wl.increment

      # return as a "short" spectrum containing only matching wls and target values
      z <- tibble::tibble(wls, target)
      names(z) <- c(col.name.x, col.name)
      if (x.class %in% spct_classes()) {
        z <- do.call(paste("as", x.class, sep = "."), args = list(x = z))
        # we need to copy our private attributes as we are building a new object
        z <- copy_attributes(x, z)
      }
    } else {
      # extract nearest wl value for target
      idxs <- ifelse(abs((x[[col.name]][closing.idx] - target) /
                           (x[[col.name]][closing.idx] - x[[col.name]][opening.idx])) > 0.5,
                     opening.idx,
                     closing.idx)
      # if the target value is close to a peak or valley, we may pick the same idx on both sides of it.
      z <- x[unique(idxs), ]
    }
    collector.ls[[as.character(target)]] <- z
  }
  if (!length(collector.ls)) {
    # we will still bind it to respect idfactor and ensure invariant columns
    collector.ls[[1L]] <- x[NULL, ]
  }
  if (is.any_spct(x)) {
    rbindspct(collector.ls,
              fill = FALSE,
              idfactor = idfactor,
              attrs.source = 1L)
  } else { # data.frame or tibble
    if (is.logical(idfactor)) {
      # convert to equivalent dplyr argument
      if (!idfactor) {
        # avoid setting of row names, for consistency
        names(collector.ls) <- NULL
      }
      idfactor <- if (idfactor) "spct.idx" else NULL
    }
    dplyr::bind_rows(collector.ls, .id = idfactor)
  }
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
#' @param target numeric value indicating the spectral quantity value for which
#'   wavelengths are to be searched and interpolated if need. The character
#'   string "half.maximum" is also accepted as argument.
#' @param interpolate logical Indicating whether the nearest wavelength value
#'   in \code{x} should be returned or a value calculated by linear
#'   interpolation between wavelength values straddling the target.
#' @param x.var.name,y.var.name,col.name character The name of the columns in
#'   which to search for the target value. Use of \code{col.name} is deprecated,
#'   and is a synonym for \code{y.var.name}.
#' @param idfactor logical or character Generates an index column of factor
#'   type. If \code{idfactor = TRUE} then the column is auto named spct.idx.
#'   Alternatively the column name can be directly passed as argument to
#'   \code{idfactor} as a character string.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for the target.
#' @param ... currently ignored.
#'
#' @return A data.frame or a spectrum object of the same class as \code{x} with
#'   fewer rows, possibly even no rows. If \code{FALSE} is passed to
#'   \code{interpolate} a subset of \code{x} is returned, otherwise a new object
#'   of the same class containing interpolated wavelengths for the \code{target}
#'   value is returned.
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
                          target = NULL,
                          interpolate = FALSE,
                          idfactor = FALSE,
                          na.rm = FALSE,
                          ...) UseMethod("wls_at_target")

#' @describeIn wls_at_target Default returning always an empty object of the
#'   same class as \code{x}.
#' @export
#'
wls_at_target.default <-
  function(x,
           target = NULL,
           interpolate = FALSE,
           idfactor = FALSE,
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
           target = "half.maximum",
           interpolate = FALSE,
           idfactor = FALSE,
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
           target = "half.maximum",
           interpolate = FALSE,
           idfactor = FALSE,
           na.rm = FALSE,
           col.name = NULL,
           y.var.name = col.name,
           ...) {
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
           target = "half.maximum",
           interpolate = FALSE,
           idfactor = FALSE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {
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
           target = "half.maximum",
           interpolate = FALSE,
           idfactor = FALSE,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...) {
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
           target = "half.maximum",
           interpolate = FALSE,
           idfactor = FALSE,
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty",
                                  default = "transmittance"),
           ...) {
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
           target = "half.maximum",
           interpolate = FALSE,
           idfactor = FALSE,
           na.rm = FALSE,
           ...) {
    find_wls(x,
             target = target,
             col.name = "Rfr",
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
           target = "half.maximum",
           interpolate = FALSE,
           idfactor = FALSE,
           na.rm = FALSE,
           ...) {
    find_wls(x,
             target = target,
             col.name = "cps",
             interpolate = interpolate,
             idfactor = idfactor,
             na.rm = na.rm)
  }

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
                                        target = "half.maximum",
                                        interpolate = FALSE,
                                        idfactor = FALSE,
                                        na.rm = FALSE,
                                        ...,
                                        .parallel = FALSE,
                                        .paropts = NULL) {
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

