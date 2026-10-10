# thinning of wavelengths --------------------------------------------------

#' Thin the density of wavelength values
#'
#' Increase the wavelength step in stored spectral data in featureless regions
#' to save storage space.
#'
#' @details The algorithm used for spectra is "naive" in an effort to keep it
#'   efficient. It works by iteratively attempting to delete every other
#'   observation along wavelengths, based on the criteria for maximum wavelength
#'   step and maximum relative step in the spectral variable between adjacent
#'   data values.
#'
#' @param x An R object
#' @param ... additional named arguments passed down to \code{f}.
#'
#' @return An object of the same class as \code{x} but with a reduced density of
#'   wavelength values in those regions were slope is shallow and featureless.
#'
#' @export
#'
#' @examples
#'
#' nrow(yellow_gel.spct)
#' wl_stepsize(yellow_gel.spct)
#' thinned.spct <- thin_wl(yellow_gel.spct)
#' nrow(thinned.spct)
#' wl_stepsize(thinned.spct)
#'
#' @family experimental utility functions
#'
thin_wl <- function(x, ...) UseMethod("thin_wl")

#' @rdname thin_wl
#'
#' @export
#'
thin_wl.default <- function(x, ...) {
  warning("'thin_wl()' is not defined for objects of class '", class(x)[1], "'.")
  x
}

#' @rdname thin_wl
#'
#' @param max.wl.step numeric. Largest allowed wavelength difference between
#'   adjacent spectral values in nanometres (nm).
#' @param max.slope.delta numeric in 0 to 1. Largest allowed change in relative
#'   slope of the spectral quantity per nm between adjacent pairs of values.
#' @param span integer A peak (or valley) is defined as an element in a sequence
#'   which is greater (or smaller) than all other elements within a window of
#'   width span centred at that element. Use NULL for the global peak.
#' @param col.names character. Name of the column of \code{x} containing the
#'   spectral data to check against \code{max.slope.delta}. Currently only one
#'   column supported.
#'
#' @note The value of \code{max.slope.delta} is expressed as relative change in
#'   the slope of spectral variable per nanometre. This means that values
#'   between 0.0005 and 0.005 tend to work reasonably well. The best value will
#'   depend on the wavelength step of the input and noise in data. A moderate
#'   smoothing before thinning can sometimes help in the case of noisy data.
#'
#'   The amount of thinning is almost always less than the value of criteria
#'   passed as argument as it is based on existing wavelength values. For
#'   example if we start with a spectrum with a uniform wavelength step of 1 nm,
#'   possible steps in the thinned spectrum are 2, 4, 8, 16, 32, etc. nm. The
#'   algorithm, does work with any step sizes, regular or variable in the input.
#'   Thinning is most effective for spectra with large "featureless" regions as
#'   the algorithm attempts not to discard information, contrary to smoothing or
#'   interpolation.
#'
#'   Local peaks and valleys are always preserved, using by default a span of 21
#'   to search for them. See \code{\link{find_peaks}}.
#'
#' @export
#'
thin_wl.generic_spct <- function(x,
                                 max.wl.step = 5.0,
                                 max.slope.delta = 0.0005,
                                 span = 15,
                                 col.names,
                                 ...) {
  # compute stopping criterion
  wl.stepsize <- wl_stepsize(x)
  max.wl.thinning <- trunc(max.wl.step / min(wl.stepsize))
  if (max.wl.thinning < 2) {
    warning("No thinning of wavelengths possible!")
    return(x)
  }
  # make code simpler by setting range to 0..1
  # we force use of parent class method as signature varies
  x.norm <- normalize.generic_spct(x,
                                   range = NULL,
                                   norm = "max",
                                   col.names = col.names,
                                   na.rm = FALSE)
  # collect peaks and valleys to ensure that they are not removed
  peaks <- find_peaks(x.norm[[col.names]], span = span, strict = FALSE)
  valleys <- find_peaks(-x.norm[[col.names]], span = span, strict = FALSE)
  extremes <- peaks | valleys

  # iterative loop
  thin.factor <- 1L
  wls.all <- x.norm[["w.length"]]
  wl.thinned <- wls.all
  var.thinned <- x.norm[[col.names]]
  wls.to.keep <- c(wls.all[c(1, length(wls.all))],
                   wls.all[extremes])
  while (thin.factor < max.wl.thinning / 2) {
    # compute differences as proxi for slope
    diff.wl <- diff(wl.thinned)
    diff.var <- diff(var.thinned)
    local.slope <- diff.var / diff.wl
    # select wavelengths to keep based on the local change in slope
    selector <-
      (abs(diff(c(0, local.slope))) > max.slope.delta) | (diff.wl > max.wl.step)
    wls.to.keep <- c(wls.to.keep,
                     wl.thinned[-1][selector])
    # keep every other value for next iteration
    wl.thinned <- wl.thinned[c(TRUE, FALSE)]
    var.thinned <- var.thinned[c(TRUE, FALSE)]
    thin.factor <- thin.factor * 2L
  }
  # add all wavelengths to ensure max wavelength step
  wls.to.keep <- unique(c(wls.to.keep, wl.thinned))
  # extract from x the rows to keep, retaining all columns keeping original
  # ordering
  x[x[["w.length"]] %in% wls.to.keep, ]
}

#' @rdname thin_wl
#'
#' @param unit.out character Allowed values "energy", and "photon", or its alias
#'   "quantum".
#'
#' @export
#'
thin_wl.source_spct <- function(x,
                                max.wl.step = 5.0,
                                max.slope.delta = 0.0005,
                                span = 15,
                                unit.out = getOption("photobiology.radiation.unit", default = "energy"),
                                ...) {
  if (unit.out == "energy") {
    thin_wl.generic_spct(x = q2e(x, action = "replace"),
                         max.wl.step = max.wl.step,
                         max.slope.delta = max.slope.delta,
                         span = span,
                         col.names = "s.e.irrad",
                         ...)
  } else if (unit.out %in% c("photon", "quantum")) {
    thin_wl.generic_spct(x = e2q(x, action = "replace"),
                         max.wl.step = max.wl.step,
                         max.slope.delta = max.slope.delta,
                         span = span,
                         col.names = "s.q.irrad",
                         ...)
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' @rdname thin_wl
#'
#' @export
#'
thin_wl.response_spct <- function(x,
                                  max.wl.step = 5.0,
                                  max.slope.delta = 0.0005,
                                  span = 15,
                                  unit.out = getOption("photobiology.radiation.unit", default = "energy"),
                                  ...) {
  if (unit.out == "energy") {
    thin_wl.generic_spct(x = q2e(x, action = "replace"),
                         max.wl.step = max.wl.step,
                         max.slope.delta = max.slope.delta,
                         span = span,
                         col.names = "s.e.response",
                         ...)
  } else if (unit.out %in% c("photon", "quantum")) {
    thin_wl.generic_spct(x = e2q(x, action = "replace"),
                         max.wl.step = max.wl.step,
                         max.slope.delta = max.slope.delta,
                         span = span,
                         col.names = "s.q.response",
                         ...)
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' @rdname thin_wl
#'
#' @param qty.out character Allowed values "transmittance", and "absorbance".
#'
#' @export
#'
thin_wl.filter_spct <- function(x,
                                max.wl.step = 5.0,
                                max.slope.delta = 0.0005,
                                span = 15,
                                qty.out = getOption("photobiology.filter.qty",
                                                    default = "transmittance"),
                                ...) {
  if (qty.out == "transmittance") {
    thin_wl.generic_spct(x = any2T(x, action = "replace"),
                         max.wl.step = max.wl.step,
                         max.slope.delta = max.slope.delta,
                         span = span,
                         col.names = "Tfr",
                         ...)
  } else if (qty.out == "absorbance") {
    thin_wl.generic_spct(x = any2A(x, action = "replace"),
                         max.wl.step = max.wl.step,
                         max.slope.delta = max.slope.delta,
                         span = span,
                         col.names = "A",
                         ...)
  } else if (qty.out == "absorptance") {
    thin_wl.generic_spct(x = any2Afr(x, action = "replace"),
                         max.wl.step = max.wl.step,
                         max.slope.delta = max.slope.delta,
                         span = span,
                         col.names = "Afr",
                         ...)
  }else {
    stop("'unit.out ", qty.out, " is unknown")
  }
}

#' @rdname thin_wl
#'
#' @export
#'
thin_wl.reflector_spct <- function(x,
                                   max.wl.step = 5.0,
                                   max.slope.delta = 0.0005,
                                   span = 15,
                                   ...) {
  thin_wl.generic_spct(x = x,
                       max.wl.step = max.wl.step,
                       max.slope.delta = max.slope.delta,
                       span = span,
                       col.names = "Rfr",
                       ...)
}

#' @rdname thin_wl
#'
#' @export
#'
thin_wl.solute_spct <- function(x,
                                max.wl.step = 5.0,
                                max.slope.delta = 0.0005,
                                span = 15,
                                ...) {
  cols <- intersect(c("K.mole", "K.mass"), names(x))
  if (length(cols) == 1) {
    col.name <- cols
  } else {
    stop("Invalid number of columns found:", length(cols))
  }
  thin_wl.generic_spct(x = x,
                       max.wl.step = max.wl.step,
                       max.slope.delta = max.slope.delta,
                       span = span,
                       col.names = col.name,
                       ...)
}

#' @rdname thin_wl
#'
#' @export
#'
thin_wl.raw_spct <- thin_wl.generic_spct

#' @rdname thin_wl
#'
#' @export
#'
thin_wl.cps_spct <- thin_wl.generic_spct

#' @rdname thin_wl
#'
#' @export
#'
thin_wl.object_spct <- thin_wl.generic_spct

#' @rdname thin_wl
#'
#' @export
#'
thin_wl.chroma_spct <- thin_wl.default

#' @rdname thin_wl
#'
#' @export
#'
thin_wl.calibration_spct <- thin_wl.default

#' @rdname thin_wl
#'
#' @export
#'
thin_wl.generic_mspct <- function(x,
                                  max.wl.step = 5.0,
                                  max.slope.delta = 0.0005,
                                  span = 15,
                                  ...) {
  msmsply(x,
          thin_wl,
          max.wl.step = max.wl.step,
          max.slope.delta = max.slope.delta,
          span = span,
          ...)
}

#' @rdname thin_wl
#'
#' @export
#'
thin_wl.chroma_mspct <- thin_wl.default

#' @rdname thin_wl
#'
#' @export
#'
thin_wl.calibration_mspct <- thin_wl.default
