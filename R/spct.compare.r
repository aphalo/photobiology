# comparison of spectra ----------------------------------------------

#' Coarse-grained comparison of two spectra
#'
#' Compare two spectra using a specified summary function pre-applied to
#' wavelength intervals.
#'
#' @details Summaries are computed for each of the wavebands in \code{w.band} by
#'   applying function \code{.summary.fun} separately to each spectrum, after
#'   trimming them to the overlapping wavelength region. Next the matching
#'   summaries are compared by means of \code{.comparison.fun}. Both the
#'   summaries and the result of the comparison are returned. Columns containing
#'   summary values are named by concatenating the name each member spectrum
#'   with the name of the argument passed to \code{.summary.fun}.
#'
#'   Tagging is useful for plotting using wavelength based colours, or when
#'   names for wavebands are used as annotations. When tagging is requested, the
#'   spectrum is passed to method \code{\link{tag}} with \code{use.hinges} and
#'   \code{short.names} as additional arguments.
#'
#' @param x A collection of two spectral objects of the same type.
#' @param .summary.fun function. The summary function to use. It must be a
#'   method accepting object \code{x} as first argument.
#' @param ... additional named arguments passed down to \code{.summary.fun}.
#' @param .comparison.fun function. The comparison function to use.
#' @param w.band waveband object or a numeric stepsize in nanometres.
#' @param returned.value character One of "data.frame", "spectrum",
#'   "tagged.spectrum".
#' @param use.hinges logical Flag indicating whether to insert "hinges" into the
#'   returned spectrum when tagging it.
#' @param short.names logical Flag indicating whether to use short or long names
#'   for wavebands when tagging.
#'
#' @return A \code{generic_spct}, tagged or not with the wavebdans, or a
#'   \code{data.frame} object containing the summary values per waveband for
#'   each spectrum and the result of applying the comparison function to these
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
#' head(
#'   compare_spct(source_mspct(list(sun1 = sun.spct, sun2 = sun.spct * 2)),
#'                returned.value = "data.frame")
#' )
#' compare_spct(source_mspct(list(sun1 = sun.spct, sun2 = sun.spct * 2)),
#'              returned.value = "tagged.spectrum")
#' compare_spct(source_mspct(list(sun1 = sun.spct, sun2 = sun.spct * 2)),
#'              returned.value = "tagged.spectrum",
#'              use.hinges = TRUE)
#'
compare_spct <- function(x,
                         w.band = 10,
                         .summary.fun = NULL,
                         ...,
                         .comparison.fun = `/`,
                         returned.value = "spectrum",
                         use.hinges = FALSE,
                         short.names = TRUE) {
  # we accept a collection with two members (easy to change to accept more)
  stopifnot(is.any_mspct(x), length(x) == 2)
  # We unset the scaled attribute to make sure default functions work
  if (is_scaled(x[[1]]) || is_scaled(x[[2]])) {
    warning("At least one of the spectra has been previously scaled.")
    if (is_scaled(x[[1]])) {
      setScaled(x[[1]], scaled = FALSE)
    }
    if (is_scaled(x[[2]])) {
      setScaled(x[[2]], scaled = FALSE)
    }
  }

  # We unset the normalized attribute to make sure default functions work
  if (is_normalized(x[[1]]) || is_normalized(x[[2]])) {
    warning("At least one of the spectra has been previously normalized")
    if (is_normalized(x[[1]])) {
      setNormalized(x[[1]], norm = FALSE)
    }
    if (is_normalized(x[[2]])) {
      setNormalized(x[[2]], norm = FALSE)
    }
  }

  # summary function default depends on class of x
  if (is.null(.summary.fun)) {
    .summary.fun <-
      switch(class(x)[1],
             source_mspct = irrad,
             response_mspct = response,
             filter_mspct = transmittance,
             reflector_mspct = reflectance,
             object_mspct = absorptance,
             NULL
      )
    # temporary kludge
    f.name <-
      switch(class(x)[1],
             source_mspct = "irrad",
             response_mspct = "response",
             filter_mspct = "transmittance",
             reflector_mspct = "reflectance",
             object_mspct = "absorptance",
             stop("No default '.summary.fun' available for class ", class(x)[1])
      )
  } else {
     f.name <- substitute(.summary.fun)
     if (is.symbol(f.name)) {
       f.name <- as.character(f.name)
     } else {
       f.name <- "summary.fun"
     }
  }
  # Skip checks for intermediate results
  prev_state <- disable_check_spct()
  on.exit(set_check_spct(prev_state), add = TRUE)
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
      w.band <- split_bands(seq(from = wl.range[1],
                                to = wl.range[2],
                                by = w.band))
    } else {
      # w.band gives the boundaries of the bands
      w.band <- split_bands(w.band)
    }
  } else if (is.waveband(w.band)) {
    w.band <- list(w.band)
  } else if (is.list(w.band)) {
  # we make sure we have a list of wavebands
    stopifnot(all(sapply(w.band, is.waveband)))
  }
  # compute summaries
  wl.mid <- sapply(w.band, wl_midpoint)
  wl.min <- sapply(w.band, wl_min)
  wl.max <- sapply(w.band, wl_max)
  summaries.tb <- .summary.fun(x, w.band = w.band, ...)
  z <- cbind(wl.mid, wl.min, wl.max, as.data.frame(t(summaries.tb[-1])))
  names(z) <- c("w.length", "wl.min", "wl.max",
                paste(names(x), f.name, sep = "."))
  z[["comparison.result"]] <- .comparison.fun(z[[5]], z[[4]])
  if (returned.value %in% c("spectrum", "tagged.spectrum")) {
    z <- as.generic_spct(z)
    if (returned.value == "tagged.spectrum") {
      z <- tag(z[ , -c(2, 3)],
               w.band = w.band,
               use.hinges = use.hinges,
               short.names = short.names)
    }
    z <- check_spct(z, force = TRUE)
  } else if (!returned.value == "data.frame") {
    warning("Returning a data frame as argument \"", returned.value,
            "\" passed to 'returned.value' is unknown.")
  }
  z
}
