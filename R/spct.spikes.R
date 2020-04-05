#' Find spikes
#'
#' This function finds spikes in a numeric vector using the algorithm of
#' Whitaker and Hayes (2018). Spikes are values in spectra that are unusually
#' high or low compared to neighbors. They are usually individual values or very
#' short runs of similar "unusual" values. Spikes caused by cosmic radiation are
#' a frequent problem in Raman spectra. Another source of spikes are "hot
#' pixels" in CCD and diode arrays. Other kinds of accidental "outlayers" will
#' be also detected.
#'
#' @details Spikes are detected based on a modified Z score calculated from the
#'   differenced spectrum. The Z threshold used should be adjusted to the
#'   characteristics of the input and desired sensitivity. The lower the
#'   threshold the more stringent the test becomes, resulting in most cases in
#'   more spikes being detected. A modified version of the algorithm is used if
#'   a value different from \code{NULL} is passed as argument to
#'   \code{max.spike.width}. In such a case, an additional step filters out
#'   broader spikes (or falsely detected steep slopes) from the returned values.
#'
#' @param x numeric vector containing spectral data.
#' @param x.is.delta logical Flag indicating if x contains already differences.
#' @param z.threshold numeric Modified Z values larger than \code{z.threshold}
#'   are considered to be spikes.
#' @param max.spike.width integer Wider regions with high Z values are not detected as
#'   spikes.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for spikes.
#'
#' @return A logical vector of the same length as \code{x}. Values that are TRUE
#'   correspond to local spikes in the data.
#'
#' @references
#' Whitaker, D. A.; Hayes, K. (2018) A simple algorithm for despiking Raman
#' spectra. Chemometrics and Intelligent Laboratory Systems, 179, 82-84.
#'
#' @export
#' @examples
#'
#' with(white_led.raw_spct,
#'      which(find_spikes(counts_3, z.threshold = 30)))
#'
#' @family peaks and valleys functions
#'
find_spikes <-
  function(x,
           x.is.delta = FALSE,
           z.threshold = 9,
           max.spike.width = 8,
           na.rm = FALSE) {
    if (na.rm) {
      na.idx <- which(is.na(x))
      x <- na.omit(x)
    }
    if (x.is.delta) {
      d.var <- x
    } else {
      d.var <- diff(x)
    }
    z <- (d.var - stats::median(d.var)) / stats::mad(d.var) * 0.6745
    outcomes <- abs(z) > z.threshold
    if (!x.is.delta) {
      # ensure same length as input
      outcomes <- c(FALSE, outcomes)
    }
    if (!is.null(max.spike.width) && max.spike.width > 0) {
      # ignore broad peaks using run length encoding
      runs <- rle(outcomes)
      runs$values <- ifelse(runs$lengths > max.spike.width, FALSE, runs$values)
      outcomes <- inverse.rle(runs)
    }
    if (na.rm) {
      # restore length of logical vector
      for (i in na.idx) {
        outcomes <- append(outcomes, FALSE, after = i - 1L)
      }
    }
    # check assertion
    stopifnot(length(outcomes) == length(x))
    outcomes
  }

#' Replace bad pixels in a spectrum
#'
#' This function replaces data for bad pixels by a local estimate, by either
#' simple interpolation or using the algorithm of Whitaker and Hayes (2018).
#'
#' @details
#' Simple interpolation replaces values of isolated bad pixels by the mean of
#' their two closest neighbors. The running mean approach allows the replacement
#' of short runs of bad pixels by the running mean of neighboring pixels within
#' a window of user-specified width. The first approach works well for spectra
#' from array spectrometers to correct for hot and dead pixels in an instrument.
#' The second approach is most suitable for Raman spectra in which spikes
#' triggered by radiation are wider than a single pixel but usually not more
#' than five pixels wide.
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
#'   generation of additional \code{NA}s during their replacement.
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
#' replace_bad_pixs(c(1, 1, 45, 1, 1), bad.pix.idx = 3)
#'
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
           window.width = 11,
           method = "run.mean",
           na.rm = TRUE) {
    if (is.logical(bad.pix.idx)) {
      if (length(bad.pix.idx) == length(x)) {
         bad.pix.idx <- which(bad.pix.idx)
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
      bad.pix.idx <- unique(c(1L, bad.pix.idx, n))
      max.spike.width <- max(rle(diff(bad.pix.idx))[["lengths"]]) + 1L
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
        z[i] = mean(x[window.idx])
      }
    } else if (method == "adj.mean") {
      # simple mean of neighbors, for isolated bad pixels.
      x[bad.pix.idx] <- NA_integer_
      if (1L %in% bad.pix.idx) {
        x[1L] <- x[2L]
        bad.pix.idx <- setdiff(bad.pix.idx, 1L)
      }
      if (n %in% bad.pix.idx) {
        x[n] <- x[n - 1L]
        bad.pix.idx <- setdiff(bad.pix.idx, n)
      }
      x[bad.pix.idx] <- (x[bad.pix.idx - 1] + x[bad.pix.idx + 1]) / 2
    }
    z
  }

# despike -------------------------------------------------------------------

#' Remove spikes from spectrum
#'
#' Function that returns an R object with observations corresponding to spikes
#' replaced by values computed from neighboring pixels. Spikes are values in
#' spectra that are unusually high compared to neighbors. They are usually
#' individual values or very short runs of similar "unusual" values. Spikes
#' caused by cosmic radiation are a frequent problem in Raman spectra. Another
#' source of spikes are "hot pixels" in CCD and diode array detectors.
#'
#' @param x an R object
#' @param z.threshold numeric Modified Z values larger than \code{z.threshold}
#'   are considered to correspond to spikes.
#' @param max.spike.width integer Wider regions with high Z values are not detected as
#'   spikes.
#' @param window.width integer. The full width of the window used for the
#'   running mean used as replacement.
#' @param method character The name of the method: \code{"run.mean"} is running
#'  mean as described in Whitaker and Hayes (2018); \code{"adj.mean"} is mean
#'  of adjacent neighbors (isolated bad pixels only).
#' @param na.rm logical indicating whether \code{NA} values should be treated
#'   as spikes and replaced.
#' @param var.name,y.var.name character Names of columns where to look
#'   for spikes to remove.
#' @param ... Arguments passed by name to \code{find_spikes()}.
#'
#' @return \code{x} with rows corresponding to spikes replaced by a local
#'   average of adjacent neighbors outside the spike.
#'
#' @note Current algorithm misidentifies steep smooth slopes as spikes, so
#'   manual inspection is needed together with adjustment by trial and error
#'   of a suitable argument value for \code{z.threshold}.
#'
#' @seealso See the documentation for \code{\link{find_spikes}} and
#'   \code{\link{replace_bad_pixs}} for details of the algorithm and
#'   implementation.
#'
#' @export
#'
#' @examples
#'
#' white_led.raw_spct[120:125, ]
#'
#' # find and replace spike at 245.93 nm
#' despike(white_led.raw_spct,
#'         z.threshold = 10,
#'         window.width = 25)[120:125, ]
#'
#' @family despike and valleys functions
#'
despike <- function(x,
                    z.threshold,
                    max.spike.width,
                    window.width,
                    method,
                    na.rm,
                    ...) UseMethod("despike")

#' @describeIn despike Default returning always NA.
#' @export
despike.default <-
  function(x,
           z.threshold = NA,
           max.spike.width = NA,
           window.width = NA,
           method = "run.mean",
           na.rm = FALSE,
           ...) {
    warning("Method 'despike' not implemented for objects of class ",
            class(x)[1])
    x[NA]
  }

#' @describeIn despike Default function usable on numeric vectors.
#' @export
despike.numeric <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           ...) {
   spike.idxs <- find_spikes(x = x,
                             z.threshold = z.threshold,
                             max.spike.width = max.spike.width,
                             na.rm = na.rm)
   replace_bad_pixs(x,
                    bad.pix.idx = spike.idxs,
                    window.width = window.width,
                    method = method,
                    na.rm = na.rm,
                    ...)
  }

#' @describeIn despike  Method for "data.frame" objects.
#'
#' @export
#'
despike.data.frame <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           ...,
           y.var.name = NULL,
           var.name = y.var.name) {
    if (is.null(var.name)) {
      warning("Variable (column) names required.")
      return(x[NA, ])
    }
    for (colname in var.name) {
      if (!is.numeric(x[[colname]])) {
        next()
      }
      x[[colname]] <- despike(x[[colname]],
                              z.threshold = z.threshold,
                              max.spike.width = max.spike.width,
                              window.width = window.width,
                              method = method,
                              na.rm = na.rm,
                              ...
      )
    }
    x
  }

#' @describeIn despike  Method for "generic_spct" objects.
#'
#' @export
#'
despike.generic_spct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           ...,
           y.var.name = NULL,
           var.name = y.var.name) {
    if (is.null(var.name)) {
      # find target variable
      var.name <- names(x)
      var.name <- subset(var.name, sapply(x, is.numeric))
      var.name <- setdiff(var.name, "w.length")
    }
    if (length(var.name) == 0L) {
      warning("No data columns found, skipping.")
    }
    for (colname in var.name) {
      if (!is.numeric(x[[colname]])) {
        next()
      }
      x[[colname]] <- despike(x[[colname]],
                              z.threshold = z.threshold,
                              max.spike.width = max.spike.width,
                              window.width = window.width,
                              method = method,
                              na.rm = na.rm,
                              ...
      )
    }
    x
  }

#' @describeIn despike  Method for "source_spct" objects.
#'
#' @param unit.out character One of "energy" or "photon"
#'
#' @export
#'
despike.source_spct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...) {
    if (unit.out == "energy") {
      z <- q2e(x, action = "replace", byref = FALSE)
      colname <- "s.e.irrad"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, action = "replace", byref = FALSE)
      colname <- "s.q.irrad"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    z[[colname]] <- despike(z[[colname]],
                            z.threshold = z.threshold,
                            max.spike.width = max.spike.width,
                            window.width = window.width,
                            method = method,
                            na.rm = na.rm,
                            ...)
    z
  }

#' @describeIn despike  Method for "response_spct" objects.
#'
#' @export
#'
despike.response_spct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...) {
    if (unit.out == "energy") {
      z <- q2e(x, action = "replace", byref = FALSE)
      colname <- "s.e.response"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, action = "replace", byref = FALSE)
      colname <- "s.q.response"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    z[[colname]] <- despike(z[[colname]],
                            z.threshold = z.threshold,
                            max.spike.width = max.spike.width,
                            window.width = window.width,
                            method = method,
                            na.rm = na.rm,
                            ...)
    z
  }

#' @describeIn despike  Method for "filter_spct" objects.
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
despike.filter_spct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty",
                                  default = "transmittance"),
           ...) {
    if (filter.qty == "transmittance") {
      z <- A2T(x, action = "replace", byref = FALSE)
      colname <- "Tfr"
    } else if (filter.qty == "absorbance") {
      z <- T2A(x, action = "replace", byref = FALSE)
      colname <- "A"
    }  else if (filter.qty == "absorptance") {
      z <- T2Afr(x, action = "replace", byref = FALSE)
      colname <- "Afr"
    } else {
      stop("Unrecognized 'filter.qty': ", filter.qty)
    }
    z[[colname]] <- despike(z[[colname]],
                            z.threshold = z.threshold,
                            max.spike.width = max.spike.width,
                            window.width = window.width,
                            method = method,
                            na.rm = na.rm,
                            ...)
    z
  }

#' @describeIn despike  Method for "reflector_spct" objects.
#'
#' @export
#'
despike.reflector_spct <- function(x,
                                   z.threshold = 9,
                                   max.spike.width = 8,
                                   window.width = 11,
                                   method = "run.mean",
                                   na.rm = FALSE,
                                   ...) {
  colname <- "Rfr"
  x[[colname]] <- despike(x[[colname]],
                          z.threshold = z.threshold,
                          max.spike.width = max.spike.width,
                          window.width = window.width,
                          method = method,
                          na.rm = na.rm,
                          ...
  )
  x
}

#' @describeIn despike  Method for "cps_spct" objects.
#'
#' @export
#'
despike.cps_spct <- function(x,
                             z.threshold = 9,
                             max.spike.width = 8,
                             window.width = 11,
                             method = "run.mean",
                             na.rm = FALSE,
                             ...) {
  var.name <- grep("cps", colnames(x), value = TRUE)
  for (colname in var.name) {
    x[[colname]] <- despike(x[[colname]],
                            z.threshold = z.threshold,
                            max.spike.width = max.spike.width,
                            window.width = window.width,
                            method = method,
                            na.rm = na.rm,
                            ...
    )
  }
  x
}

#' @describeIn despike  Method for "raw_spct" objects.
#'
#' @export
#'
despike.raw_spct <- function(x,
                             z.threshold = 9,
                             max.spike.width = 8,
                             window.width = 11,
                             method = "run.mean",
                             na.rm = FALSE,
                             ...) {
  var.name <- grep("counts", colnames(x), value = TRUE)
  for (colname in var.name) {
    x[[colname]] <- despike(x[[colname]],
                            z.threshold = z.threshold,
                            max.spike.width = max.spike.width,
                            window.width = window.width,
                            method = method,
                            na.rm = na.rm,
                            ...
    )
  }
  x
}

#' @describeIn despike  Method for "generic_mspct" objects.
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
                                  z.threshold = 9,
                                  max.spike.width = 8,
                                  window.width = 11,
                                  method = "run.mean",
                                  na.rm = FALSE,
                                  ...,
                                  y.var.name = NULL,
                                  var.name = y.var.name,
                                  .parallel = FALSE,
                                  .paropts = NULL) {
  msmsply(x,
          .fun = despike,
          z.threshold = z.threshold,
          max.spike.width = max.spike.width,
          window.width = window.width,
          method = method,
          na.rm = na.rm,
          var.name = var.name,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn despike  Method for "source_mspct" objects.
#'
#' @export
#'
despike.source_mspct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = despike,
            z.threshold = z.threshold,
            max.spike.width = max.spike.width,
            window.width = window.width,
            method = method,
            na.rm = na.rm,
            unit.out = unit.out,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @describeIn despike  Method for "cps_mspct" objects.
#'
#' @export
#'
despike.response_mspct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = despike,
            z.threshold = z.threshold,
            max.spike.width = max.spike.width,
            window.width = window.width,
            method = method,
            na.rm = na.rm,
            unit.out = unit.out,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @describeIn despike  Method for "filter_mspct" objects.
#'
#' @export
#'
despike.filter_mspct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty",
                                  default = "transmittance"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = despike,
            z.threshold = z.threshold,
            max.spike.width = max.spike.width,
            window.width = window.width,
            method = method,
            filter.qty = filter.qty,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }


#' @describeIn despike  Method for "reflector_mspct" objects.
#'
#' @export
#'
despike.reflector_mspct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           window.width = 11,
           method = "run.mean",
           na.rm = FALSE,
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = despike,
            z.threshold = z.threshold,
            max.spike.width = max.spike.width,
            window.width = window.width,
            method = method,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }


#' @describeIn despike  Method for "cps_mspct" objects.
#'
#' @export
#'
despike.cps_mspct <- function(x,
                              z.threshold = 9,
                              max.spike.width = 8,
                              window.width = 11,
                              method = "run.mean",
                              na.rm = FALSE,
                              ...,
                              .parallel = FALSE,
                              .paropts = NULL) {
  msmsply(x,
          .fun = despike,
          z.threshold = z.threshold,
          max.spike.width = max.spike.width,
          window.width = window.width,
          method = method,
          na.rm = na.rm,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn despike  Method for "raw_mspct" objects.
#'
#' @export
#'
despike.raw_mspct <- function(x,
                              z.threshold = 9,
                              max.spike.width = 8,
                              window.width = 11,
                              method = "run.mean",
                              na.rm = FALSE,
                              ...,
                              .parallel = FALSE,
                              .paropts = NULL) {
  msmsply(x,
          .fun = despike,
          z.threshold = z.threshold,
          max.spike.width = max.spike.width,
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
#' "unusual" values. Spikes caused by cosmic radiation are a frequent problem in
#' Raman spectra. Another source of spikes are "hot pixels" in CCD and diode
#' arrays.
#'
#' @param x an R object
#' @param z.threshold numeric Modified Z values larger than \code{z.threshold}
#'   are considered to correspond to spikes.
#' @param max.spike.width integer Wider regions with high Z values are not
#'   detected as spikes.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for spikes.
#' @param var.name,y.var.name character Name of column where to look
#'   for spikes.
#' @param ... ignored
#'
#' @return A subset of \code{x} with rows corresponding to spikes.
#'
#' @seealso See the documentation for \code{\link{find_spikes}} for details of
#'   the algorithm and implementation.
#'
#' @export
#'
#' @examples
#' spikes(sun.spct)
#'
#' @family peaks and valleys functions
#'
spikes <- function(x, z.threshold, max.spike.width, na.rm, ...) UseMethod("spikes")

#' @describeIn spikes Default returning always NA.
#' @export
spikes.default <-
  function(x,
           z.threshold = NA,
           max.spike.width = 8,
           na.rm = FALSE,
           ...) {
    warning("Method 'spikes' not implemented for objects of class ",
            class(x)[1])
    x[NA]
  }

#' @describeIn spikes Default function usable on numeric vectors.
#' @export
spikes.numeric <-
  function(x,
           z.threshold = NA,
           max.spike.width = 8,
           na.rm = FALSE,
           ...) {
    x[find_spikes(x = x,
                  z.threshold = z.threshold,
                  max.spike.width = max.spike.width,
                  na.rm = na.rm)]
  }

#' @describeIn spikes  Method for "data.frame" objects.
#'
#' @export
#'
spikes.data.frame <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           na.rm = FALSE,
           ...,
           y.var.name = NULL,
           var.name = y.var.name) {
    if (is.null(var.name)) {
      warning("Variable (column) names required.")
      return(x[NA, ])
    }
    spikes.idx <-
      which(find_spikes(x[[var.name]],
                        z.threshold = z.threshold,
                        max.spike.width = max.spike.width,
                        na.rm = na.rm))
    x[spikes.idx,  , drop = FALSE]
  }

#' @describeIn spikes  Method for "generic_spct" objects.
#'
#' @export
#'
spikes.generic_spct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           na.rm = FALSE,
           ...,
           var.name = NULL) {
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
                        z.threshold = z.threshold,
                        max.spike.width = max.spike.width,
                        na.rm = na.rm))
    x[spikes.idx,  , drop = FALSE]
  }

#' @describeIn spikes  Method for "source_spct" objects.
#'
#' @param unit.out character One of "energy" or "photon"
#'
#' @export
#'
spikes.source_spct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
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
    spikes.idx <-
      which(find_spikes(z[[col.name]],
                        z.threshold = z.threshold,
                        max.spike.width = max.spike.width,
                        na.rm = na.rm))
    z[spikes.idx,  , drop = FALSE]
  }

#' @describeIn spikes  Method for "response_spct" objects.
#'
#' @export
#'
spikes.response_spct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
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
    spikes.idx <-
      which(find_spikes(z[[col.name]],
                        z.threshold = z.threshold,
                        max.spike.width = max.spike.width,
                        na.rm = na.rm))
    z[spikes.idx,  , drop = FALSE]
  }

#' @describeIn spikes  Method for "filter_spct" objects.
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
spikes.filter_spct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
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
    spikes.idx <-
      which(find_spikes(z[[col.name]],
                        z.threshold = z.threshold,
                        max.spike.width = max.spike.width,
                        na.rm = na.rm))
    z[spikes.idx,  , drop = FALSE]
  }

#' @describeIn spikes  Method for "reflector_spct" objects.
#'
#' @export
#'
spikes.reflector_spct <- function(x,
                                  z.threshold = 9,
                                  max.spike.width = 8,
                                  na.rm = FALSE,
                                  ...) {
  col.name <- "Rfr"
  spikes.idx <-
    which(find_spikes(x[[col.name]],
                      z.threshold = z.threshold,
                      max.spike.width = max.spike.width,
                      na.rm = na.rm))
  x[spikes.idx,  , drop = FALSE]
}

#' @describeIn spikes  Method for "cps_spct" objects.
#'
#' @export
#'
spikes.cps_spct <- function(x,
                            z.threshold = 9,
                            max.spike.width = 8,
                            na.rm = FALSE,
                            ...,
                            var.name = "cps") {
  spikes.idx <-
    which(find_spikes(x[[var.name]],
                      z.threshold = z.threshold,
                      max.spike.width = max.spike.width,
                      na.rm = na.rm))
  x[spikes.idx,  , drop = FALSE]
}

#' @describeIn spikes  Method for "raw_spct" objects.
#'
#' @export
#'
spikes.raw_spct <- function(x,
                           z.threshold = 9,
                           max.spike.width = 8,
                           na.rm = FALSE,
                            ...,
                           var.name = "counts") {
  spikes.idx <-
    which(find_spikes(x[[var.name]],
                      z.threshold = z.threshold,
                      max.spike.width = max.spike.width,
                      na.rm = na.rm))
  x[spikes.idx,  , drop = FALSE]
}

#' @describeIn spikes  Method for "generic_mspct" objects.
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
spikes.generic_mspct <- function(x,
                                 z.threshold = 9,
                                 max.spike.width = 8,
                                 na.rm = FALSE,
                                 ...,
                                 var.name = NULL,
                                 .parallel = FALSE,
                                 .paropts = NULL) {
  msmsply(x,
          .fun = spikes,
          z.threshold = z.threshold,
          max.spike.width = max.spike.width,
          na.rm = na.rm,
          var.name = var.name,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn spikes  Method for "source_mspct" objects.
#'
#' @export
#'
spikes.source_mspct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = spikes,
            z.threshold = z.threshold,
            max.spike.width = max.spike.width,
            unit.out = unit.out,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @describeIn spikes  Method for "cps_mspct" objects.
#'
#' @export
#'
spikes.response_mspct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           na.rm = FALSE,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = spikes,
            z.threshold = z.threshold,
            max.spike.width = max.spike.width,
            unit.out = unit.out,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }

#' @describeIn spikes  Method for "filter_mspct" objects.
#'
#' @export
#'
spikes.filter_mspct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           na.rm = FALSE,
           filter.qty = getOption("photobiology.filter.qty",
                                  default = "transmittance"),
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = spikes,
            z.threshold = z.threshold,
            max.spike.width = max.spike.width,
            filter.qty = filter.qty,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }


#' @describeIn spikes  Method for "reflector_mspct" objects.
#'
#' @export
#'
spikes.reflector_mspct <-
  function(x,
           z.threshold = 9,
           max.spike.width = 8,
           na.rm = FALSE,
           ...,
           .parallel = FALSE,
           .paropts = NULL) {
    msmsply(x,
            .fun = spikes,
            z.threshold = z.threshold,
            max.spike.width = max.spike.width,
            na.rm = na.rm,
            ...,
            .parallel = .parallel,
            .paropts = .paropts)
  }


#' @describeIn spikes  Method for "cps_mspct" objects.
#'
#' @export
#'
spikes.cps_mspct <- function(x,
                             z.threshold = 9,
                             max.spike.width = 8,
                             na.rm = FALSE,
                             ...,
                             var.name = "cps",
                             .parallel = FALSE,
                             .paropts = NULL) {
  msmsply(x,
          .fun = spikes,
          z.threshold = z.threshold,
          max.spike.width = max.spike.width,
          na.rm = na.rm,
          var.name = var.name,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn spikes  Method for "raw_mspct" objects.
#'
#' @export
#'
spikes.raw_mspct <- function(x,
                             z.threshold = 9,
                             max.spike.width = 8,
                             na.rm = FALSE,
                             ...,
                             var.name = "counts",
                             .parallel = FALSE,
                             .paropts = NULL) {
  msmsply(x,
          .fun = spikes,
          z.threshold = z.threshold,
          max.spike.width = max.spike.width,
          na.rm = na.rm,
          var.name = var.name,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

