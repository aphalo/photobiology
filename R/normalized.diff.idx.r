#' Calculate a normalized index.
#'
#' This method returns a normalized difference index value for an arbitrary
#' pair of wavebands. There are many such indexes in use, such as NDVI
#' (normalized difference vegetation index), NDWI (normalized difference water
#' index), NDMI (normalized difference moisture index), etc., the only
#' difference among then is in the wavebands used.
#'
#' @param spct an R object
#' @param plus.w.band waveband objects The waveband determine the
#'   region of the spectrum used in the calculations
#' @param minus.w.band waveband objects The waveband determine the
#'   region of the spectrum used in the calculations
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
#'   argument. Wavebands can describe weighting functions if desired.
#'
#' @note Some NDxI indexes are directly based on satellite instrument data, such
#'   as those in the Landsat satellites. To simulate such indexes using spectral
#'   reflectande as input, \code{waveband} definitions provided by package
#'   'photobiologyWavebands' can be used.
#'
#' @export
#'
normalized_diff_ind <-
  function(spct, plus.w.band, minus.w.band, f, ...) {
    UseMethod("normalized_diff_ind")
  }

#' @rdname normalized_diff_ind
#'
#' @note \code{normalised_diff_ind()} is a synonym for \code{normalized_diff_ind()}.
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
  function(spct, plus.w.band, minus.w.band, f, ...) {
    warning("'normalized_diff_ind' is not defined for objects of class ",
            class(spct)[1])
    return(spct)
}

#' @describeIn normalized_diff_ind
#'
#' @export
#'
normalized_diff_ind.generic_spct <- function(spct, plus.w.band, minus.w.band, f, ...) {
  # check that spectral data fully covers both wavebands
  min.wl.bands <- min(wl_min(plus.w.band), wl_min(minus.w.band))
  max.wl.bands <- max(wl_max(plus.w.band), wl_max(minus.w.band))
  if (wl_min(spct) > min.wl.bands || wl_max(spct) < max.wl.bands) {
    NA_real_
  } else {
    x <- as.numeric(f(spct, list(plus.w.band, minus.w.band), ...))
    z <- (x[1] - x[2]) / (x[1] + x[2])
    name <- paste("NDI ", as.character(substitute(f)), " [",
                  sub("range.", "", labels(plus.w.band)$label), "] - [",
                  sub("range.", "", labels(minus.w.band)$label), "]",
                  sep = "")
    names(z) <- name
    z
  }
}

#' @describeIn normalized_diff_ind
#'
#' @export
#'
normalized_diff_ind.generic_mspct <- function(spct, plus.w.band, minus.w.band, f, ...) {
  msdply(mspct = spct,
         plus.w.band = plus.w.band,
         minus.w.band = minus.w.band,
         f = f,
         ...)
}


