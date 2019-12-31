# comparison of spectra ----------------------------------------------

#' Coarse-grained comparison of two spectra
#'
#' Compare two spectra using a specified summary function pre-applied to
#' wavelength intervals.
#'
#' @param x A collection of two spectral objects of the same type.
#' @param .summary.fun function. The summary function to use.
#' @param .comparison.fun function. The comparison function to use.
#' @param w.band waveband object or a numeric stepsize in nanometres.
#' @param ... additional named arguments passed down to \code{.summary.fun}.
#'
#' @return A data frame containing the summary values per waveband for each
#'   spectrum and the result of applying the comparison function to these
#'   summaries.
#'
#' @export
#'
#' @examples
#'
#' compare_spct(source_mspct(list(sun1 = sun.spct, sun2 = sun.spct * 2)))
#' compare_spct(source_mspct(list(sun1 = sun.spct, sun2 = sun.spct * 2)),
#'              w.band = NULL)
#' compare_spct(source_mspct(list(sun1 = sun.spct, sun2 = sun.spct * 2)),
#'              w.band = list(waveband(c(640, 650)), waveband(c(720, 740))))
#'
#' compare_spct(filter_mspct(list(pet = polyester.spct,
#'                                yllw = yellow_gel.spct)),
#'              w.band = 50,
#'              .comparison.fun = `<`)
#'
compare_spct <- function(x,
                         w.band = 10,
                         .summary.fun = NULL,
                         .comparison.fun = `/`,
                         ...) {
  # summary function default depends on class of x
  if (is.null(.summary.fun)) {
    .summary.fun <-
      switch(class(x)[1],
             source_mspct = irrad,
             filter_mspct = transmittance,
             reflector_mspct = reflectance,
             object_mspct = absorptance
      )
    # temporary kludge
    f.name <-
      switch(class(x)[1],
             source_mspct = "irrad",
             filter_mspct = "transmittance",
             reflector_mspct = "reflectance",
             object_mspct = "absorptance"
      )
  } else {
     f.name <- substitute(.summary.fun)
     if (is.symbol(f.name)) {
       f.name <- as.character(f.name)
     } else {
       f.name <- "summary.fun"
     }
  }
  # we accept a collection with two members (easy to change to accept more)
  stopifnot(is.any_mspct(x), length(x) == 2)
  # keep overlapping wavelength range
  x <- trim2overlap(x)
  # w.band can take different arguments
  wl.range <- wl_range(x[[1]])
  if (is.null(w.band)) {
    # we summarize the whole spectrum at once
    w.band <- list(waveband(wl.range))
  } else if (is.numeric(w.band)) {
    if (length(w.band) == 1L) {
      # w.band gives the stepsize
      w.band <- split_bands(seq(from = wl.range[1], to = wl.range[2], by = w.band))
    } else {
      # w.band gives the boundaries of the bands
      w.band <- split_bands(w.band)
    }
  } else if (is.waveband(w.band)) {
    w.band <- list(w.band)
  }
  # we make sure we have a list of wavebands
  stopifnot(all(sapply(w.band, is.waveband)))
  # compute summaries
  wl.mid <- sapply(w.band, wl_midpoint)
  wl.min <- sapply(w.band, wl_min)
  wl.max <- sapply(w.band, wl_max)
  summaries.tb <- .summary.fun(x, w.band = w.band)
  z <- cbind(wl.mid, wl.min, wl.max, as.data.frame(t(summaries.tb[-1])))
  names(z) <- c("w.length", "wl.min", "wl.max", paste(names(x), f.name, sep = "."))
  z[["comparison.result"]] <- .comparison.fun(z[[5]], z[[4]])
  z
}
