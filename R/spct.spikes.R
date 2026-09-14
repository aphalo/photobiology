#' Find spikes in vector
#'
#' Find spikes in a numeric vector using the algorithm of
#' Whitaker and Hayes (2018). Spikes are values in spectra that are unusually
#' high or low compared to neighbours. They are usually individual values or very
#' short runs of similar "unusual" values. Spikes caused by cosmic radiation are
#' a frequent problem in Raman spectra. Another source of spikes are "hot
#' pixels" in CCD and diode arrays. Other kinds of accidental "outliers" can
#' be also detected.
#'
#' @section Spike detection:
#'   Spikes are detected based on a modified \eqn{Z} score calculated
#'   from the differenced spectrum. The \eqn{Z} threshold used should be
#'   adjusted to the characteristics of the input and desired sensitivity. The
#'   lower the threshold the more stringent the test becomes, with shorter
#'   spikes being detected.
#'
#'   \strong{The algorithms assume a consistent step size for the underlying
#'   independent variable, e.g., wavelength or time, and should not be applied
#'   if the data do not fulfil this assumption, at least approximately. As
#'   \code{find_spkikes()} operates on a single vector, checking this remains
#'   the responsibility of calling functions or methods such as
#'   \code{\link{spikes}()} and \code{\link{despike}()}.}
#'
#'   The algorithm uses running differences to detect abrupt changes in value,
#'   compared to an estimate of the baseline variation of the differences,
#'   approximating a baseline \eqn{Z} from MAD and a baseline value from the
#'   median differences. Currently, a single estimate of MAD is used but running
#'   medians, when possible, as baseline. This comparison detects running
#'   differences that are unusually large, in most cases signalling a transition
#'   between values near the baseline and far from it, in both directions.
#'
#'   Transitions into- and out of spikes are distinguished based on the median
#'   of the non-differenced values, as a descriptor of the data baseline. As for
#'   the median of the differences, a running median is used when possible.
#'
#'   This function thus detects the start and end of each spike, and
#'   distinguishes upward and downward spikes.
#'
#'   \code{k} is the width in number of observations of the window used for
#'   running median smoothing to extract the baseline. A value several times the
#'   width of the broader spike but narrow enough to track broader peaks needs
#'   to be manually set in most cases.
#'
#'   With \code{na.rm = TRUE}, \code{NA} values are omitted before searching for
#'   spikes and set to \code{0L} in the returned vector.
#'
#'   If all spikes are guaranteed to be one observation-wide and either going up
#'   or down from the baseline, it is possible to detect them based purely on
#'   the \code{z.threshold} by passing \code{height.threshold = NA} and either
#'   \code{spike.direction = "up"} or \code{spike.direction = "down"}, which
#'   ensures very fast computation.
#'
#'   Parameters of the algorithm need to be adjusted depending on the data, so
#'   inspection of returned values is needed together with adjustment by trial
#'   and error of suitable values for \code{z.threshold},
#'   \code{height.threshold}, and \code{k}.
#'
#'   Parameter \code{max.spike.width} searches for too wide spikes in the
#'   output of the algorithms described above and ignores them. This is
#'   possibly redundant, but maintained for partial backwards compatibility.
#'
#' @param x numeric vector containing the data.
#' @param x.is.delta logical Flag indicating whether \code{x} contains
#'   differences or original values.
#' @param height.threshold numeric The minimum height of spikes expressed
#'   relative to the median amplitude of the baseline local variation of
#'   \code{x}.
#' @param z.threshold numeric Modified local \eqn{Z} values larger than
#'   \code{z.threshold} are detected as boundaries of spikes.
#' @param k integer width of median window used for smoothing; must be odd
#' @param spike.direction character Controls the direction of spikes to be
#'   detected. Accepted arguments are \code{"up"}, \code{"down"},
#'   \code{"both"}.
#' @param return.numeric logical If \code{TRUE} a numeric vector is returned
#'   and otherwise a logical one.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for spikes.
#' @param max.spike.width integer The width of the widest spike to be detected,
#'   \code{NA} puts no limit.
#'
#' @return An integer vector of the same length as \code{x}. Values that are
#'   \code{0}, \code{+1} or \code{-1} corresponding to no-spike, upwards-spike,
#'   and downwards-spike in the data. Conversion to logical with
#'   \code{as.logical()} results in a vector with \code{TRUE} for spikes and
#'   \code{FALSE} otherwise.
#'
#' @references
#' Whitaker, D. A.; Hayes, K. (2018) A simple algorithm for despiking Raman
#' spectra. Chemometrics and Intelligent Laboratory Systems, 179, 82-84.
#' \doi{10.1016/j.chemolab.2018.06.009}.
#'
#' @export
#'
#' @family peaks and valleys functions
#'
find_spikes <-
  function(x,
           x.is.delta = FALSE,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           return.numeric = FALSE,
           na.rm = FALSE,
           max.spike.width = NA) {
    if (is.null(height.threshold)) {
      height.threshold <- 10
    } else if (!is.na(height.threshold) && height.threshold < 2) {
      warning("'height.threshold < 2' set to 2")
      height.threshold <- 2
    }
    if (is.null(k)) {
      k <- 20
    } else if (k %% 2 == 0) {
      k <- k + 1
    }
    x.len.original <- length(x)
    if (na.rm) {
      na.idx <- which(is.na(x))
      x <- na.omit(x)
    }
    if (x.is.delta) {
      d.var <- x
      x <- stats::diffinv(x)
    } else {
      d.var <- diff(x)
      x <- x - x[1]
    }
    # running median is used to detect spikes relative to the local baseline
    if (k > length(x) / 2) {
      x.median <- stats::median(x)
      d.var.median <- stats::median(d.var)
    } else {
      x.median <- stats::runmed(x,
                                k = k,
                                na.action = "na.omit",
                                endrule = "constant")
      d.var.median <- stats::runmed(d.var,
                                    k = k,
                                    na.action = "na.omit",
                                    endrule = "constant")
    }
    z <- (d.var - d.var.median) / stats::mad(d.var) * 0.6745
    outcomes.up <- c(FALSE, z > z.threshold)
    outcomes.down <- c(FALSE, z < -z.threshold)

    if (is.na(height.threshold)) {
      spikes.up <- outcomes.up
      spikes.down <- outcomes.down
    } else {
      scaled.threshold <- stats::median(abs(d.var.median)) * height.threshold
      if (spike.direction %in% c("up", "both")) {
        outcomes.head.up <-
          outcomes.up & x > x.median + scaled.threshold
        temp <-
          outcomes.down &
          # near the baseline
          x <= x.median + scaled.threshold &
          x >= x.median - scaled.threshold
        outcomes.tail.up <- logical(length(temp))
        outcomes.tail.up[which(temp) - 1L] <- TRUE

        # fill gaps
        spk.starts <- which(outcomes.head.up)
        spk.ends <- which(outcomes.tail.up)

        if (length(spk.starts) > 1 && length(spk.ends) > 1) {
          # check if data ends or starts in a spike
          if (spk.ends[1] < spk.starts[1]) {
            if (spk.ends[1] > 1) {
              spk.starts <- c(1, spk.starts)
            } else {
              spk.ends <- spk.ends[-1L]
            }
          }
          if (spk.starts[length(spk.starts)] > spk.ends[length(spk.ends)]) {
            if (spk.ends[length(spk.ends)] < length(x)) {
              spk.ends <- c(spk.ends, length(x))
            } else {
              spk.starts <- spk.starts[-length(spk.starts)]
            }
          }
          outcomes.middle.up <- logical(length(x))
          i <- j <- 0
          i.max <- length(spk.starts)
          j.max <- length(spk.ends)
          while (i < i.max && j < j.max) {
            i <- i + 1
            j <- j + 1
            # skip narrow spikes
            while (i < i.max && spk.starts[i + 1] < spk.ends[j]) i <- i + 1
            while (j < j.max && spk.ends[j + 1] < spk.starts[i]) j <- j + 1
            # fill in the middle of wide spikes
            if (spk.starts[i] + 1 < spk.ends[j]) {
              outcomes.middle.up[(spk.starts[i] + 1):(spk.ends[j] - 1)] <- TRUE
            }
          }
          spikes.up <- outcomes.head.up | outcomes.tail.up | outcomes.middle.up
        } else {
          spikes.up <- outcomes.up
        }
        spikes.up <-
          spikes.up & x > x.median + scaled.threshold
      }

      if (!is.null(max.spike.width) && !is.na(max.spike.width)) {
        widths <- rle(spikes.up)
        widths$values[which(widths$lengths > max.spike.width)] <- FALSE
        spikes.up <- inverse.rle(widths)
      }

      if (spike.direction %in% c("down", "both")) {
        outcomes.head.down <-
          outcomes.up & x < x.median - scaled.threshold
        temp <-
          outcomes.up &
          # near the baseline
          x <= x.median + scaled.threshold &
          x >= x.median - scaled.threshold
        outcomes.tail.down <- logical(length(temp))
        outcomes.tail.down[which(temp) - 1L] <- TRUE

        # fill gaps
        spk.starts <- which(outcomes.head.down)
        spk.ends <- which(outcomes.tail.down)

        if (length(spk.starts) > 1 && length(spk.ends) > 1) {
          # check if data ends or starts in a spike
          if (spk.ends[1] < spk.starts[1]) {
            if (spk.ends[1] > 1) {
              spk.starts <- c(1, spk.starts)
            } else {
              spk.ends <- spk.ends[-1L]
            }
          }
          if (spk.starts[length(spk.starts)] > spk.ends[length(spk.ends)]) {
            if (spk.ends[length(spk.ends)] < length(x)) {
              spk.ends <- c(spk.ends, length(x))
            } else {
              spk.starts <- spk.starts[-length(spk.starts)]
            }
          }
          outcomes.middle.down <- logical(length(x))
          i <- j <- 0
          i.max <- length(spk.starts)
          j.max <- length(spk.ends)
          while (i < i.max && j < j.max) {
            i <- i + 1
            j <- j + 1
            # skip narrow spikes
            while (i < i.max && spk.starts[i + 1] < spk.ends[j]) i <- i + 1
            while (j < j.max && spk.ends[j + 1] < spk.starts[i]) j <- j + 1
            # fill in the middle of wide spikes
            if (spk.starts[i] + 1 < spk.ends[j]) {
              outcomes.middle.down[(spk.starts[i] + 1):(spk.ends[j] - 1)] <- TRUE
            }
          }
          spikes.down <- outcomes.head.down | outcomes.tail.down | outcomes.middle.down
        } else {
          spikes.down <- outcomes.down
        }
        spikes.down <-
          spikes.down & x < x.median - scaled.threshold

        temp <-
          outcomes.up &
          # near the baseline
          x <= x.median + scaled.threshold&
          x >= x.median - scaled.threshold
        outcomes.tail.down <- logical(length(temp))
        outcomes.tail.down[which(temp) - 1L] <- TRUE
        spikes.down <- outcomes.down | outcomes.tail.down
        spikes.down <-
          spikes.down & x < x.median - scaled.threshold
      }
    }

    if (!is.null(max.spike.width) && !is.na(max.spike.width)) {
      widths <- rle(spikes.down)
      widths$values[which(widths$lengths > max.spike.width)] <- FALSE
      spikes.down <- inverse.rle(widths)
    }

    outcomes <-
      switch(spike.direction,
             "up" = spikes.up * 1L,
             "down" = spikes.down * -1L,
             "both" = spikes.up + spikes.down * -1L,
             "skip" = integer(length(x)),
             {
               warning("'spike.direction' must be \"up\", \"down\", \"both\", or \"skip\", not \"",
                       spike.direction, "\"")
               integer(length(x))
             }
      )

    if (na.rm) {
      # restore length of logical vector
      for (i in na.idx) {
        outcomes <- append(outcomes, FALSE, after = i - 1L)
      }
    }
    # check assertion
    stopifnot(length(outcomes) == x.len.original)
    if (return.numeric) {
      outcomes
    } else {
      as.logical(outcomes)
    }
  }

#' Replace bad pixels in a spectrum
#'
#' This function replaces data for bad pixels by a local estimate, by either
#' simple interpolation or using the algorithm of Whitaker and Hayes (2018).
#'
#' @section Replacement values:
#' Simple interpolation enabled by \code{method = "adj.mean"} replaces values of
#' isolated bad pixels by the mean of their two closest neighbours. The running
#' mean approach enabled by \code{method = "run.mean"} allows the replacement of
#' short runs of bad pixels by the running mean of neighboring pixels within a
#' window of user-specified width. The first approach works well for spectra
#' from array spectrometers to correct for hot and dead pixels in an instrument.
#' The second approach is most suitable for Raman spectra in which spikes
#' triggered by radiation are wider than a single pixel but usually not more
#' than five pixels wide.
#'
#' Simple interpolation can replace spikes at any position in \code{x}, using
#' a single neighbour as replacement at the extremes of \code{x} instead of the
#' mean of two neighbours. The
#' running mean approach does not replace those pixels whose distance to the
#' first or last member of \code{x} is less than half the window used for the
#' running mean, issuing a warning.
#'
#' When \code{na.rm = TRUE}, \code{NA} values are considered "bad pixels" and
#' replaced as such rather than discarded with no replacement. This is the
#' default behaviour.
#'
#' @param x numeric vector containing spectral data.
#' @param bad.pix.idx logical vector or integer. Index into bad pixels in
#'   \code{x}.
#' @param window.width integer. The full width of the window used for the
#'   running mean.
#' @param method character The name of the method: \code{"run.mean"} is running
#'  mean as described in Whitaker and Hayes (2018); \code{"adj.mean"} is mean
#'  of adjacent neighbors (isolated bad pixels only).
#' @param na.rm logical Treat \code{NA} values as additional bad pixels and
#'  replace them.
#'
#' @note In the current implementation \code{NA} values are not removed, and
#'   if they are in the neighborhood of bad pixels, they will result in the
#'   generation of additional \code{NA}s during their replacement. On the other
#'   hand if the \code{NA}s locations are listed in \code{bad.pix.idx} they
#'   will be replaced as any other bad pixel.
#'
#' @return A logical vector of the same length as \code{x}. Values that are TRUE
#'   correspond to local spikes in the data.
#'
#' @references
#' Whitaker, D. A.; Hayes, K. (2018) A simple algorithm for despiking Raman
#' spectra. Chemometrics and Intelligent Laboratory Systems, 179, 82-84.
#'
#' @examples
#' # in a vector
#' replace_bad_pixs(c(1, 2, NA, 4, 5))
#'
#' # in a vector
#' replace_bad_pixs(c(1, 2, 100, 4, 5),
#'                  method = "adj.mean",
#'                  bad.pix.idx = c(FALSE, FALSE, TRUE, FALSE, FALSE))
#'
#' replace_bad_pixs(c(1, 2, 100, 4, 5),
#'                  method = "adj.mean",
#'                  bad.pix.idx = 3)
#'
#' # in a vector
#' replace_bad_pixs(c(0, 1, 2, 100, 4, 5, 6),
#'                  method = "run.mean",
#'                  bad.pix.idx = 4)
#'
#' # in a vector
#' replace_bad_pixs(c(1, 1, NA, 1, 1),
#'                  method = "run.mean",
#'                  bad.pix.idx = 3)
#'
#' # in a vector
#' replace_bad_pixs(c(1, 1, NA, 1, 1),
#'                  method = "run.mean",
#'                  bad.pix.idx = 1, na.rm = FALSE)
#'
#' # In spectrum
#' # before replacement
#' white_led.raw_spct$counts_3[120:125]
#'
#' # replacing bad pixels at index positions 123 and 1994
#' with(white_led.raw_spct,
#'      replace_bad_pixs(counts_3, bad.pix.idx = c(123, 1994)))[120:125]
#'
#' @export
#'
#' @family peaks and valleys functions
#'
replace_bad_pixs <-
  function(x,
           bad.pix.idx = FALSE,
           window.width =  min(11, length(x) - 1),
           method = "run.mean",
           na.rm = TRUE) {
    if (is.logical(bad.pix.idx)) {
      if (length(bad.pix.idx) == length(x)) {
         bad.pix.idx <- which(bad.pix.idx)
      } else if (length(bad.pix.idx) == 1L) {
        if (bad.pix.idx) {
          return(rep(NA_real_, length(x)))
        } else {
          bad.pix.idx <- integer(0)
        }
      } else {
        stop("Logical 'bad.pix.idx' has wrong length.")
      }
    }
    if (na.rm) {
      bad.pix.idx <- union(bad.pix.idx, which(is.na(x)))
    }
    if (length(bad.pix.idx) == 0L) {
      # nothing to do
      return(x)
    }
    if (length(window.width) == 0L) {
      # force computation of a wide enough window
      window.width <- 0L
    }
    n <- length(x)
    z <- x
    if (method == "run.mean") {
      if (length(bad.pix.idx) > 1L) {
        max.spike.width <- max(rle(diff(bad.pix.idx))[["lengths"]] + 1L)
      } else {
        max.spike.width <- 1L
      }
      needed.window.width <- 2L * max.spike.width + 1L
      if (window.width < needed.window.width) {
        if (window.width > 0L) {
          warning("Increasing 'window.width' from ", window.width,
                  " to ", needed.window.width)
        }
        window.width <- needed.window.width
      }
      half.window.width <- window.width %/% 2 # half window
      # running mean method of Whitaker and Hayes (2018)
      # fast but biased.
      for(i in bad.pix.idx) {
        window.idx <- seq(max(1 , i - half.window.width),
                          min(n, i + half.window.width))
        window.idx <- setdiff(window.idx, bad.pix.idx)
        if (any(window.idx < 1 | window.idx > n)) {
          warning("Bad pixel at position ", i,
                  "not replaced! Too near edge.")
        }
        z[i] = mean(x[window.idx])
      }
    } else if (method == "adj.mean") {
      # simple mean of neighbors, for isolated bad pixels.
      z[bad.pix.idx] <- NA_integer_
      if (1L %in% bad.pix.idx) {
        z[1L] <- z[2L]
        bad.pix.idx <- setdiff(bad.pix.idx, 1L)
      }
      if (n %in% bad.pix.idx) {
        z[n] <- z[n - 1L]
        bad.pix.idx <- setdiff(bad.pix.idx, n)
      }
      z[bad.pix.idx] <- (z[bad.pix.idx - 1] + z[bad.pix.idx + 1]) / 2
    }
    z
  }

# despike -------------------------------------------------------------------

#' Remove spikes from spectrum
#'
#' Function that returns an R object with observations corresponding to spikes
#' replaced by values computed from neighboring pixels. Spikes are values in
#' spectra that are unusually high compared to neighbors. They are usually
#' individual values or very short runs of similar "unusual" values.
#'
#' @inheritSection find_spikes Spike detection
#'
#' @inheritSection replace_bad_pixs Replacement values
#'
#' @inheritParams find_spikes
#' @inheritParams replace_bad_pixs
#' @param var.name,y.var.name character Names of columns where to look
#'   for spikes to remove.
#' @param ... passed in recursive calls.
#'
#' @return A copy of the object passed as argument to \code{x} with values
#'   detected as spikes replaced by a local average of neighbours
#'   outside the spike.
#'
#' @seealso See \code{\link{find_spikes}()} for locating spikes in a vector,
#'   \code{\link{spikes}()} for extracting/detecting spikes in spectra and
#'   and \code{\link{replace_bad_pixs}()} for replacing by interpolation
#'   missing or bad values in a vector.
#'
#' @export
#'
#' @examples
#'
#' white_led.raw_spct[120:125, ]
#'
#' # find and replace spike at 245.93 nm
#' despike(white_led.raw_spct,
#'         z.threshold = 5,
#'         window.width = 7)[120:125, ]
#'
#' # A high z.threshold value detects more extreme spikes
#' despike(white_led.raw_spct,
#'         z.threshold = 50,
#'         window.width = 7)[120:125, ]
#'
#' @family despike and valleys functions
#'
despike <- function(x,
                    height.threshold,
                    z.threshold,
                    k,
                    spike.direction,
                    window.width,
                    method,
                    na.rm,
                    max.spike.width,
                    ...) UseMethod("despike")

#' @rdname despike
#'
#' @export
despike.default <-
  function(x,
           height.threshold,
           z.threshold = NA,
           k = NA,
           spike.direction = NA,
           window.width = NA,
           method = "run.mean",
           na.rm = FALSE,
           max.spike.width = NA,
           ...) {
    warning("Method 'despike' not implemented for objects of class ",
            class(x)[1])
    x[NA]
  }

#' @rdname despike
#'
#' @export
despike.numeric <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           max.spike.width = NA,
           ...) {
   spike.idxs <- find_spikes(x = x,
                             height.threshold = height.threshold,
                             z.threshold = z.threshold,
                             k = k,
                             spike.direction = spike.direction,
                             na.rm = na.rm,
                             max.spike.width = max.spike.width) |>
     as.logical()

   replace_bad_pixs(x,
                    bad.pix.idx = spike.idxs,
                    window.width = window.width,
                    method = method,
                    na.rm = na.rm,
                    ...)
  }

#' @rdname despike
#'
#' @export
#'
despike.data.frame <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           max.spike.width = NA,
           ...,
           y.var.name = NULL,
           var.name = y.var.name) {
    if (is.null(var.name)) {
      warning("Variable (column) names required.")
      return(x[NA, ])
    }
    for (col.name in var.name) {
      if (!is.numeric(x[[col.name]])) {
        next()
      }
      x[[col.name]] <- despike(x[[col.name]],
                               height.threshold = height.threshold,
                               z.threshold = z.threshold,
                               k = k,
                               spike.direction = spike.direction,
                               window.width = window.width,
                               method = method,
                               na.rm = na.rm,
                               max.spike.width = max.spike.width,
                               ...)
    }
    x
  }

#' @rdname despike
#'
#' @export
#'
despike.generic_spct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           max.spike.width = NA,
           ...,
           y.var.name = NULL,
           var.name = y.var.name) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- despike(x = mspct,
                       height.threshold = height.threshold,
                       z.threshold = z.threshold,
                       k = k,
                       spike.direction = spike.direction,
                       window.width = window.width,
                       method = method,
                       na.rm = na.rm,
                       max.spike.width = max.spike.width,
                       y.var.name = y.var.name,
                       var.name = var.name,
                       ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    if (!check_wl_stepsize(x, span = 15, min.stepsize = 3)) {
      warning("Despike skipped!")
      return(x)
    }

    if (is.null(var.name)) {
      # find target variable
      var.name <- names(x)
      var.name <- subset(var.name, sapply(x, is.numeric))
      var.name <- setdiff(var.name, "w.length")
    }
    if (length(var.name) == 0L) {
      warning("No data columns found, skipping.")
    }

    for (col.name in var.name) {
      if (!is.numeric(x[[col.name]])) {
        next()
      }
      x[[col.name]] <- despike(x[[col.name]],
                               height.threshold = height.threshold,
                               z.threshold = z.threshold,
                               k = k,
                               spike.direction = spike.direction,
                               window.width = window.width,
                               method = method,
                               na.rm = na.rm,
                               max.spike.width = max.spike.width,
                               ...
      )
    }
    x
  }

#' @rdname despike
#'
#' @param unit.out character One of "energy" or "photon"
#'
#' @export
#'
despike.source_spct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           max.spike.width = NA,
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
      mspct <- despike(x = mspct,
                       height.threshold = height.threshold,
                       z.threshold = z.threshold,
                       k = k,
                       spike.direction = spike.direction,
                       window.width = window.width,
                       method = method,
                       na.rm = na.rm,
                       max.spike.width = max.spike.width,
                       unit.out = unit.out,
                       ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    if (unit.out == "energy") {
      z <- q2e(x, action = "replace", byref = FALSE)
      col.name <- "s.e.irrad"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, action = "replace", byref = FALSE)
      col.name <- "s.q.irrad"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }

    if (!check_wl_stepsize(z, span = 15, min.stepsize = 3)) {
      warning("Despike skipped!")
      return(z)
    }

    z[[col.name]] <- despike(z[[col.name]],
                             height.threshold = height.threshold,
                             z.threshold = z.threshold,
                             k = k,
                             spike.direction = spike.direction,
                             window.width = window.width,
                             method = method,
                             na.rm = na.rm,
                             max.spike.width = max.spike.width,
                             ...)
    z
  }

#' @rdname despike
#'
#' @export
#'
despike.response_spct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           max.spike.width = NA,
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
      mspct <- despike(x = mspct,
                       height.threshold = height.threshold,
                       z.threshold = z.threshold,
                       k = k,
                       spike.direction = spike.direction,
                       window.width = window.width,
                       method = method,
                       na.rm = na.rm,
                       max.spike.width = max.spike.width,
                       unit.out = unit.out,
                       ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    if (unit.out == "energy") {
      z <- q2e(x, action = "replace", byref = FALSE)
      col.name <- "s.e.response"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, action = "replace", byref = FALSE)
      col.name <- "s.q.response"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }

    if (!check_wl_stepsize(z, span = 15, min.stepsize = 3)) {
      warning("Despike skipped!")
      return(z)
    }

    z[[col.name]] <- despike(z[[col.name]],
                             height.threshold = height.threshold,
                             z.threshold = z.threshold,
                             k = k,
                             spike.direction = spike.direction,
                             window.width = window.width,
                             method = method,
                             na.rm = na.rm,
                             max.spike.width = max.spike.width,
                             ...)
    z
  }

#' @rdname despike
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
despike.filter_spct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           max.spike.width = NA,
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
      mspct <- despike(x = mspct,
                       height.threshold = height.threshold,
                       z.threshold = z.threshold,
                       k = k,
                       spike.direction = spike.direction,
                       window.width = window.width,
                       method = method,
                       na.rm = na.rm,
                       max.spike.width = max.spike.width,
                       filter.qty = filter.qty,
                       ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    if (filter.qty == "transmittance") {
      z <- A2T(x, action = "replace", byref = FALSE)
      col.name <- "Tfr"
    } else if (filter.qty == "absorbance") {
      z <- T2A(x, action = "replace", byref = FALSE)
      col.name <- "A"
    }  else if (filter.qty == "absorptance") {
      z <- T2Afr(x, action = "replace", byref = FALSE)
      col.name <- "Afr"
    } else {
      stop("Unrecognized 'filter.qty': ", filter.qty)
    }

    if (!check_wl_stepsize(z, span = 15, min.stepsize = 3)) {
      warning("Despike skipped!")
      return(z)
    }

    z[[col.name]] <- despike(z[[col.name]],
                             height.threshold = height.threshold,
                             z.threshold = z.threshold,
                             k = k,
                             spike.direction = spike.direction,
                             window.width = window.width,
                             method = method,
                             na.rm = na.rm,
                             max.spike.width = max.spike.width,
                             ...)
    z
  }

#' @rdname despike
#'
#' @export
#'
despike.reflector_spct <- function(x,
                                   height.threshold = 10,
                                   z.threshold = 5,
                                   k = 20,
                                   spike.direction = "both",
                                   window.width = 11,
                                   method = "run.mean",
                                   na.rm = FALSE,
                                   max.spike.width = NA,
                                   ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(x) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    # call method on the collection
    mspct <- despike(x = mspct,
                     height.threshold = height.threshold,
                     z.threshold = z.threshold,
                     k = k,
                     spike.direction = spike.direction,
                     window.width = window.width,
                     method = method,
                     na.rm = na.rm,
                     max.spike.width = max.spike.width,
                     ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  if (!check_wl_stepsize(x, span = 15, min.stepsize = 3)) {
    warning("Despike skipped!")
  }

  col.name <- "Rfr"
  x[[col.name]] <- despike(x[[col.name]],
                           height.threshold = height.threshold,
                           z.threshold = z.threshold,
                           k = k,
                           spike.direction = spike.direction,
                           method = method,
                           na.rm = na.rm,
                           max.spike.width = max.spike.width,
                           ...
  )
  x
}

#' @rdname despike
#'
#' @export
#'
despike.solute_spct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           max.spike.width = NA,
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- despike(x = mspct,
                       height.threshold = height.threshold,
                       z.threshold = z.threshold,
                       k = k,
                       spike.direction = spike.direction,
                       window.width = window.width,
                       method = method,
                       na.rm = na.rm,
                       max.spike.width = max.spike.width,
                       ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    if (!check_wl_stepsize(x, span = 15, min.stepsize = 3)) {
      warning("Despike skipped!")
    }

    cols <- intersect(c("K.mole", "K.mass"), names(x))
    if (length(cols) == 1) {
      col.name <- cols
      z <- x
    } else {
      stop("Invalid number of columns found:", length(cols))
    }
    z[[col.name]] <- despike(z[[col.name]],
                             height.threshold = height.threshold,
                             z.threshold = z.threshold,
                             k = k,
                             spike.direction = spike.direction,
                             window.width = window.width,
                             method = method,
                             na.rm = na.rm,
                             max.spike.width = max.spike.width,
                             ...)
    z
  }

#' @rdname despike
#'
#' @export
#'
despike.cps_spct <- function(x,
                             height.threshold = 10,
                             z.threshold = 5,
                             k = 20,
                             spike.direction = "both",
                             window.width = 11,
                             method = "run.mean",
                             na.rm = FALSE,
                             max.spike.width = NA,
                             ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(x) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    # call method on the collection
    mspct <- despike(x = mspct,
                     height.threshold = height.threshold,
                     z.threshold = z.threshold,
                     k = k,
                     spike.direction = spike.direction,
                     window.width = window.width,
                     method = method,
                     na.rm = na.rm,
                     max.spike.width = max.spike.width,
                     ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  if (!check_wl_stepsize(x, span = 15, min.stepsize = 3)) {
    warning("Despike skipped!")
  }

  var.name <- grep("cps", colnames(x), value = TRUE)
  for (col.name in var.name) {
    x[[col.name]] <- despike(x[[col.name]],
                             height.threshold = height.threshold,
                             z.threshold = z.threshold,
                             k = k,
                             spike.direction = spike.direction,
                             window.width = window.width,
                             method = method,
                             na.rm = na.rm,
                             max.spike.width = max.spike.width,
                             ...
    )
  }
  x
}

#' @rdname despike
#'
#' @export
#'
despike.raw_spct <- function(x,
                             height.threshold = 10,
                             z.threshold = 5,
                             k = 20,
                             spike.direction = "both",
                             window.width = 11,
                             method = "run.mean",
                             na.rm = FALSE,
                             max.spike.width = NA,
                             ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(x) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    # call method on the collection
    mspct <- despike(x = mspct,
                     height.threshold = height.threshold,
                     z.threshold = z.threshold,
                     k = k,
                     spike.direction = spike.direction,
                     window.width = window.width,
                     method = method,
                     na.rm = na.rm,
                     max.spike.width = max.spike.width,
                     ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  if (!check_wl_stepsize(x, span = 15, min.stepsize = 3)) {
    warning("Despike skipped!")
  }

  var.name <- grep("counts", colnames(x), value = TRUE)
  for (col.name in var.name) {
    x[[col.name]] <- despike(x[[col.name]],
                             height.threshold = height.threshold,
                             z.threshold = z.threshold,
                             k = k,
                             spike.direction = spike.direction,
                             window.width = window.width,
                             method = method,
                             na.rm = na.rm,
                             max.spike.width = max.spike.width,
                             ...
    )
  }
  x
}


# _mspct methods ----------------------------------------------------------

#' @rdname despike
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
despike.generic_mspct <- function(x,
                                  height.threshold = 10,
                                  z.threshold = 5,
                                  k = 20,
                                  spike.direction = "both",
                                  window.width = 11,
                                  method = "run.mean",
                                  na.rm = FALSE,
                                  max.spike.width = NA,
                                  ...,
                                  y.var.name = NULL,
                                  var.name = y.var.name,
                                  .parallel = FALSE,
                                  .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(x,
          .fun = despike,
          height.threshold = height.threshold,
          z.threshold = z.threshold,
          k = k,
          spike.direction = spike.direction,
          window.width = window.width,
          method = method,
          na.rm = na.rm,
          var.name = var.name,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @rdname despike
#'
#' @export
#'
despike.source_mspct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           max.spike.width = NA,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = despike,
            height.threshold = height.threshold,
            z.threshold = z.threshold,
            k = k,
            spike.direction = spike.direction,
            window.width = window.width,
            method = method,
            na.rm = na.rm,
            unit.out = unit.out,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @rdname despike
#'
#' @export
#'
despike.response_mspct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           max.spike.width = NA,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = despike,
            height.threshold = height.threshold,
            z.threshold = z.threshold,
            k = k,
            spike.direction = spike.direction,
            window.width = window.width,
            method = method,
            na.rm = na.rm,
            unit.out = unit.out,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @rdname despike
#'
#' @export
#'
despike.filter_mspct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           max.spike.width = NA,
           filter.qty = getOption("photobiology.filter.qty",
                                  default = "transmittance"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = despike,
            height.threshold = height.threshold,
            z.threshold = z.threshold,
            k = k,
            spike.direction = spike.direction,
            window.width = window.width,
            method = method,
            filter.qty = filter.qty,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }


#' @rdname despike
#'
#' @export
#'
despike.reflector_mspct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           max.spike.width = NA,
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = despike,
            height.threshold = height.threshold,
            z.threshold = z.threshold,
            k = k,
            spike.direction = spike.direction,
            window.width = window.width,
            method = method,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @rdname despike
#'
#' @export
#'
despike.solute_mspct <- despike.reflector_mspct

#' @rdname despike
#'
#' @export
#'
despike.cps_mspct <- function(x,
                              height.threshold = 10,
                              z.threshold = 5,
                              k = 20,
                              spike.direction = "both",
                              window.width = 11,
                              method = "run.mean",
                              na.rm = FALSE,
                              max.spike.width = NA,
                              ...,
                              .parallel = FALSE,
                              .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(x,
          .fun = despike,
          height.threshold = height.threshold,
          z.threshold = z.threshold,
          k = k,
          spike.direction = spike.direction,
          window.width = window.width,
          method = method,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @rdname despike
#'
#' @export
#'
despike.raw_mspct <- function(x,
                              height.threshold = 10,
                              z.threshold = 5,
                              k = 20,
                              spike.direction = "both",
                              window.width = 11,
                              method = "run.mean",
                              na.rm = FALSE,
                              max.spike.width = NA,
                              ...,
                              .parallel = FALSE,
                              .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(x,
          .fun = despike,
          height.threshold = height.threshold,
          z.threshold = z.threshold,
          k = k,
          spike.direction = spike.direction,
          window.width = window.width,
          method = method,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

# spikes -------------------------------------------------------------------

#' Spikes
#'
#' Function that returns a subset of an R object with observations corresponding
#' to spikes. Spikes are values in spectra that are unusually high compared to
#' neighbors. They are usually individual values or very short runs of similar
#' "unusual" values.
#'
#' @inheritSection find_spikes Spike detection
#'
#' @inheritParams find_spikes
#' @param var.name,y.var.name character Name of column where to look
#'   for spikes.
#' @param ... ignored
#'
#' @return A subset of the object passed as argument to \code{x} with rows
#'   corresponding to spikes.
#'
#' @seealso See \code{\link{find_spikes}()} for locating spikes in a vector,
#'   \code{\link{despike}()} for replacement of spikes by interpolation in
#'   spectra and and \code{\link{replace_bad_pixs}()} for replacing by
#'   interpolation missing or bad values in a vector.
#'
#' @export
#'
#' @examples
#' spikes(sun.spct)
#'
#' @family peaks and valleys functions
#'
spikes <- function(x,
                   height.threshold,
                   z.threshold,
                   k,
                   spike.direction,
                   na.rm,
                   max.spike.width,
                   ...) UseMethod("spikes")

#' @rdname spikes
#'
#' @export
spikes.default <-
  function(x,
           height.threshold = NA,
           z.threshold = NA,
           k = NA,
           spike.direction = NA,
           na.rm = FALSE,
           max.spike.width = NA,
           ...) {
    warning("Method 'spikes' not implemented for objects of class ",
            class(x)[1])
    x[NA]
  }

#' @rdname spikes
#'
#' @export
spikes.numeric <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           na.rm = FALSE,
           max.spike.width = NA,
           ...) {
    x[find_spikes(x = x,
                  height.threshold = height.threshold,
                  z.threshold = z.threshold,
                  k = k,
                  spike.direction = spike.direction,
                  na.rm = na.rm,
                  max.spike.width = max.spike.width)]
  }

#' @rdname spikes
#'
#' @export
#'
spikes.data.frame <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           na.rm = FALSE,
           max.spike.width = NA,
           ...,
           y.var.name = NULL,
           var.name = y.var.name) {
    if (is.null(var.name)) {
      warning("Variable (column) names required.")
      return(x[NA, ])
    }
    spikes.idx <-
      which(find_spikes(x[[var.name]],
                        height.threshold = height.threshold,
                        z.threshold = z.threshold,
                        k = k,
                        spike.direction = spike.direction,
                        na.rm = na.rm,
                        max.spike.width = max.spike.width))
    x[spikes.idx,  , drop = FALSE]
  }

#' @rdname spikes
#'
#' @export
#'
spikes.generic_spct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           na.rm = FALSE,
           max.spike.width = NA,
           var.name = NULL,
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- spikes(x = mspct,
                      height.threshold = height.threshold,
                      z.threshold = z.threshold,
                      k = k,
                      spike.direction = spike.direction,
                      na.rm = na.rm,
                      max.spike.width = max.spike.width,
                      var.name = var.name,
                      ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    if (!check_wl_stepsize(x, span = 15, min.stepsize = 3)) {
      warning("Detection of spikes in unreliable!")
    }

    if (is.null(var.name)) {
      # find target variable
      var.name <- names(x)
      var.name <- subset(var.name, sapply(x, is.numeric))
      var.name <- setdiff(var.name, "w.length")
      if (length(var.name) > 1L) {
        warning("Multiple numeric data columns found, explicit argument to",
                "'var.name' required.")
        return(x[NA, ])
      }
    }
    spikes.idx <-
      which(find_spikes(x[[var.name]],
                        height.threshold = height.threshold,
                        z.threshold = z.threshold,
                        k = k,
                        spike.direction = spike.direction,
                        na.rm = na.rm,
                        max.spike.width = max.spike.width))
    x[spikes.idx,  , drop = FALSE]
  }

#' @rdname spikes
#'
#' @param unit.out character One of "energy" or "photon"
#'
#' @export
#'
spikes.source_spct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           na.rm = FALSE,
           max.spike.width = NA,
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
      mspct <- spikes(x = mspct,
                      height.threshold = height.threshold,
                      z.threshold = z.threshold,
                      k = k,
                      spike.direction = spike.direction,
                      na.rm = na.rm,
                      max.spike.width = max.spike.width,
                      unit.out = unit.out,
                      ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    if (!check_wl_stepsize(x, span = 15, min.stepsize = 3)) {
      warning("Detection of spikes in unreliable!")
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
    spikes.idx <-
      which(find_spikes(z[[col.name]],
                        height.threshold = height.threshold,
                        z.threshold = z.threshold,
                        k = k,
                        spike.direction = spike.direction,
                        na.rm = na.rm,
                        max.spike.width = max.spike.width))
    z[spikes.idx,  , drop = FALSE]
  }

#' @rdname spikes
#'
#' @export
#'
spikes.response_spct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           na.rm = FALSE,
           max.spike.width = NA,
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
      mspct <- spikes(x = mspct,
                      height.threshold = height.threshold,
                      z.threshold = z.threshold,
                      k = k,
                      spike.direction = spike.direction,
                      na.rm = na.rm,
                      max.spike.width = max.spike.width,
                      unit.out = unit.out,
                      ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    if (!check_wl_stepsize(x, span = 15, min.stepsize = 3)) {
      warning("Detection of spikes in unreliable!")
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
    spikes.idx <-
      which(find_spikes(z[[col.name]],
                        height.threshold = 10,
                        z.threshold = 5,
                        k = 20,
                        spike.direction = "both",
                        na.rm = na.rm,
                        max.spike.width = max.spike.width))
    z[spikes.idx,  , drop = FALSE]
  }

#' @rdname spikes
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
spikes.filter_spct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           na.rm = FALSE,
           max.spike.width = NA,
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
      mspct <- spikes(x = mspct,
                      height.threshold = height.threshold,
                      z.threshold = z.threshold,
                      k = k,
                      spike.direction = spike.direction,
                      na.rm = na.rm,
                      max.spike.width = max.spike.width,
                      filter.qty = filter.qty,
                      ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    if (!check_wl_stepsize(x, span = 15, min.stepsize = 3)) {
      warning("Detection of spikes in unreliable!")
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
    spikes.idx <-
      which(find_spikes(z[[col.name]],
                        height.threshold = height.threshold,
                        z.threshold = z.threshold,
                        k = k,
                        spike.direction = spike.direction,
                        na.rm = na.rm,
                        max.spike.width = max.spike.width))
    z[spikes.idx,  , drop = FALSE]
  }

#' @rdname spikes
#'
#' @export
#'
spikes.reflector_spct <- function(x,
                                  height.threshold = 10,
                                  z.threshold = 5,
                                  k = 20,
                                  spike.direction = "both",
                                  na.rm = FALSE,
                                  max.spike.width = NA,
                                  ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(x) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    # call method on the collection
    mspct <- spikes(x = mspct,
                    height.threshold = height.threshold,
                    z.threshold = z.threshold,
                    k = k,
                    spike.direction = spike.direction,
                    na.rm = na.rm,
                    max.spike.width = max.spike.width,
                    ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  if (!check_wl_stepsize(x, span = 15, min.stepsize = 3)) {
    warning("Detection of spikes in unreliable!")
  }

  col.name <- "Rfr"
  spikes.idx <-
    which(find_spikes(x[[col.name]],
                      height.threshold = height.threshold,
                      z.threshold = z.threshold,
                      k = k,
                      spike.direction = spike.direction,
                      na.rm = na.rm))
  x[spikes.idx,  , drop = FALSE]
}

#' @rdname spikes
#'
#' @export
#'
spikes.solute_spct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           na.rm = FALSE,
           max.spike.width = NA,
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(x) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = x,
                            idx.var = getIdFactor(x),
                            drop.idx = FALSE)
      # call method on the collection
      mspct <- spikes(x = mspct,
                      height.threshold = height.threshold,
                      z.threshold = z.threshold,
                      k = k,
                      spike.direction = spike.direction,
                      na.rm = na.rm,
                      max.spike.width = max.spike.width,
                      ...)
      return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
    }

    if (!check_wl_stepsize(x, span = 15, min.stepsize = 3)) {
      warning("Detection of spikes in unreliable!")
    }

    cols <- intersect(c("K.mole", "K.mass"), names(x))
    if (length(cols) == 1) {
      col.name <- cols
      z <- x
    } else {
      stop("Invalid number of columns found:", length(cols))
    }
    spikes.idx <-
      which(find_spikes(z[[col.name]],
                        height.threshold = height.threshold,
                        z.threshold = z.threshold,
                        k = k,
                        spike.direction = spike.direction,
                        na.rm = na.rm,
                        max.spike.width = max.spike.width))
    z[spikes.idx,  , drop = FALSE]
  }

#' @rdname spikes
#'
#' @export
#'
spikes.cps_spct <- function(x,
                            height.threshold = 10,
                            z.threshold = 5,
                            k = 20,
                            spike.direction = "both",
                            na.rm = FALSE,
                            max.spike.width = NA,
                            var.name = "cps",
                            ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(x) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    # call method on the collection
    mspct <- spikes(x = mspct,
                    height.threshold = height.threshold,
                    z.threshold = z.threshold,
                    k = k,
                    spike.direction = spike.direction,
                    na.rm = na.rm,
                    max.spike.width = max.spike.width,
                    var.name = var.name,
                    ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  if (!check_wl_stepsize(x, span = 15, min.stepsize = 3)) {
    warning("Detection of spikes in unreliable!")
  }

  spikes.idx <-
    which(find_spikes(x[[var.name]],
                      height.threshold = height.threshold,
                      z.threshold = z.threshold,
                      k = k,
                      spike.direction = spike.direction,
                      na.rm = na.rm,
                      max.spike.width = max.spike.width))
  x[spikes.idx,  , drop = FALSE]
}

#' @rdname spikes
#'
#' @export
#'
spikes.raw_spct <- function(x,
                            height.threshold = 10,
                            z.threshold = 5,
                            k = 20,
                            spike.direction = "both",
                            na.rm = FALSE,
                            max.spike.width = NA,
                            var.name = "counts",
                            ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(x) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = x,
                          idx.var = getIdFactor(x),
                          drop.idx = FALSE)
    # call method on the collection
    mspct <- spikes(x = mspct,
                    height.threshold = height.threshold,
                    z.threshold = z.threshold,
                    k = k,
                    spike.direction = spike.direction,
                    na.rm = na.rm,
                    max.spike.width = max.spike.width,
                    var.name = var.name,
                    ...)
    return(rbindspct(mspct, idfactor = getIdFactor(x), attrs.simplify = TRUE))
  }

  check_wl_stepsize(x, span = 15, min.stepsize = 3)

  spikes.idx <-
    which(find_spikes(x[[var.name]],
                      height.threshold = height.threshold,
                      z.threshold = z.threshold,
                      k = k,
                      spike.direction = spike.direction,
                      na.rm = na.rm,
                      max.spike.width = max.spike.width))
  x[spikes.idx,  , drop = FALSE]
}

#' @rdname spikes
#'
#' @export
#'
spikes.generic_mspct <- function(x,
                                 height.threshold = 10,
                                 z.threshold = 5,
                                 k = 20,
                                 spike.direction = "both",
                                 na.rm = FALSE,
                                 max.spike.width = NA,
                                 ...,
                                 var.name = NULL,
                                 .parallel = FALSE,
                                 .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(x,
          .fun = spikes,
          height.threshold = height.threshold,
          z.threshold = z.threshold,
          k = k,
          spike.direction = spike.direction,
          na.rm = na.rm,
          var.name = var.name,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @rdname spikes
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
spikes.source_mspct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           na.rm = FALSE,
           max.spike.width = NA,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = spikes,
            height.threshold = height.threshold,
            z.threshold = z.threshold,
            k = k,
            spike.direction = spike.direction,
            unit.out = unit.out,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @rdname spikes
#'
#' @export
#'
spikes.response_mspct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           na.rm = FALSE,
           max.spike.width = NA,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = spikes,
            height.threshold = height.threshold,
            z.threshold = z.threshold,
            k = k,
            spike.direction = spike.direction,
            unit.out = unit.out,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @rdname spikes
#'
#' @export
#'
spikes.filter_mspct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           na.rm = FALSE,
           max.spike.width = NA,
           filter.qty = getOption("photobiology.filter.qty",
                                  default = "transmittance"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = spikes,
            height.threshold = height.threshold,
            z.threshold = z.threshold,
            k = k,
            spike.direction = spike.direction,
            filter.qty = filter.qty,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }


#' @rdname spikes
#'
#' @export
#'
spikes.reflector_mspct <-
  function(x,
           height.threshold = 10,
           z.threshold = 5,
           k = 20,
           spike.direction = "both",
           na.rm = FALSE,
           max.spike.width = NA,
           ...,
           .parallel = FALSE,
           .paropts = NULL) {

    x <- subset2mspct(x) # expand long form spectra within collection

    msmsply(x,
            .fun = spikes,
            height.threshold = height.threshold,
            z.threshold = z.threshold,
            k = k,
            spike.direction = spike.direction,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @rdname spikes
#'
#' @export
#'
spikes.solute_mspct <- spikes.reflector_mspct


#' @rdname spikes
#'
#' @export
#'
spikes.cps_mspct <- function(x,
                             height.threshold = 10,
                             z.threshold = 5,
                             k = 20,
                             spike.direction = "both",
                             na.rm = FALSE,
                             max.spike.width = NA,
                             ...,
                             var.name = "cps",
                             .parallel = FALSE,
                             .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(x,
          .fun = spikes,
          height.threshold = height.threshold,
          z.threshold = z.threshold,
          k = k,
          spike.direction = spike.direction,
          na.rm = na.rm,
          var.name = var.name,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @rdname spikes
#'
#' @export
#'
spikes.raw_mspct <- function(x,
                             height.threshold = 10,
                             z.threshold = 5,
                             k = 20,
                             spike.direction = "both",
                             na.rm = FALSE,
                             max.spike.width = NA,
                             ...,
                             var.name = "counts",
                             .parallel = FALSE,
                             .paropts = NULL) {

  x <- subset2mspct(x) # expand long form spectra within collection

  msmsply(x,
          .fun = spikes,
          height.threshold = height.threshold,
          z.threshold = z.threshold,
          k = k,
          spike.direction = spike.direction,
          na.rm = na.rm,
          var.name = var.name,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

