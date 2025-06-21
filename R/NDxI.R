#' Calculate a normalized difference.
#'
#' This method returns a normalized difference index value for an arbitrary
#' pair of wavebands. There are many such indexes in use, such as NDVI
#' (normalized difference vegetation index), NDWI (normalized difference water
#' index), NDMI (normalized difference moisture index), etc., the only
#' difference among then is in the wavebands used.
#'
#' @param spct an R object
#' @param w.band.plus,w.band.minus waveband objects The wavebands determine the
#'   regions of the spectrum used in the calculations.
#' @param f function used for integration taking spct as first argument and a
#'   list of wavebands as second argument.
#' @param ... additional arguments passed to f
#'
#' @return A named numeric value for the index, or a tibble depending on whether
#'   a spectrum or a collection of spectra is passed as first argument. If
#'   the wavelength range of \code{spct} does not fully overlap with both
#'   wavebands \code{NA} is silently returned.
#'
#' @export
#'
#' @details \code{f} is most frequently \code{\link{reflectance}}, but also
#'   \code{\link{transmittance}}, or even \code{\link{absorbance}},
#'   \code{\link{response}}, \code{\link{irradiance}} or a user-defined function
#'   can be used if there is a good reason for it. In every case \code{spct}
#'   should be of the class expected by \code{f}. When using two wavebands of
#'   different widths do consider passing to \code{f} a suitable \code{quantity}
#'   argument, for example to compare averages rather than integrals. Wavebands
#'   can describe weighting functions if desired.
#'
#'   \deqn{\mathrm{NDxI} = \frac{f(s, wb_\mathrm{plus}) - f(s, wb_\mathrm{minus})}{f(s, wb_\mathrm{plus}) + f(s, wb_\mathrm{minus})}}
#'
#' @note Some NDxI indexes are directly based on satellite instrument data, such
#'   as those in the Landsat satellites. To simulate such indexes using spectral
#'   reflectande as input, constructors of \code{waveband} definitions from
#'   package 'photobiologyWavebands' can be useful.
#'
#' @seealso \code{\link{Rfr_normdiff}}
#'
#' @export
#'
normalized_diff_ind <-
  function(spct, w.band.plus, w.band.minus, f, ...) {
    UseMethod("normalized_diff_ind")
  }

#' @rdname normalized_diff_ind
#'
#' @note \code{normalised_diff_ind()} is a synonym for
#'   \code{normalized_diff_ind()}.
#'
#' @export
#'
normalised_diff_ind <- normalized_diff_ind

#' @rdname normalized_diff_ind
#'
#' @note \code{NDxI()} is a shorthand for \code{normalized_diff_ind()}.
#'
#' @export
#'
NDxI <- normalized_diff_ind

#' @describeIn normalized_diff_ind default
#'
#' @export
#'
normalized_diff_ind.default <-
  function(spct, w.band.plus, w.band.minus, f, ...) {
    warning("'normalized_diff_ind' is not defined for objects of class ",
            class(spct)[1])
    NA_real_
  }

#' @describeIn normalized_diff_ind
#'
#' @export
#'
normalized_diff_ind.generic_spct <- function(spct,
                                             w.band.plus,
                                             w.band.minus,
                                             f,
                                             ...) {

  # we look for multiple spectra in long form
  if (getMultipleWl(spct) > 1) {
    # convert to a collection of spectra
    mspct <- subset2mspct(x = spct,
                          idx.var = getIdFactor(spct),
                          drop.idx = FALSE)
    # call method on the collection
    return(normalized_diff_ind(spct = mspct,
                               w.band.plus =  w.band.plus,
                               w.band.minus = w.band.minus,
                               f = f,
                               ...))
  }

  # check that spectral data fully covers both wavebands
  min.wl.bands <- min(wl_min(w.band.plus), wl_min(w.band.minus))
  max.wl.bands <- max(wl_max(w.band.plus), wl_max(w.band.minus))
  if (wl_min(spct) > min.wl.bands || wl_max(spct) < max.wl.bands) {
    NA_real_
  } else {
    x <- as.numeric(f(spct, list(w.band.plus, w.band.minus), ...))
    z <- (x[1] - x[2]) / (x[1] + x[2])
    name <- paste("NDI ", as.character(substitute(f)), " [",
                  sub("range.", "", labels(w.band.plus)[["label"]]), "] - [",
                  sub("range.", "", labels(w.band.minus)[["label"]]), "]",
                  sep = "")
    names(z) <- name
    z
  }
}

#' @describeIn normalized_diff_ind
#'
#' @export
#'
normalized_diff_ind.generic_mspct <- function(spct,
                                              w.band.plus,
                                              w.band.minus,
                                              f,
                                              ...) {

  spct <- subset2mspct(spct) # expand long form spectra within collection

  msdply(mspct = spct,
         .fun = normalized_diff_ind,
         w.band.plus = w.band.plus,
         w.band.minus = w.band.minus,
         f = f,
         ...)
}
