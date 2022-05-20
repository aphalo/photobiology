#' Join all spectra in a collection
#'
#' Join all the spectra contained in a homogeneous collection, returning a data
#' frame with spectral-data columns named according to the names of the spectra
#' in the collection. By default a full join is done within the overlapping range
#' of wavelengths, after interpolating the spectra to a shared set of wavelength
#' values, and discarding data for wavelength not shared. Alternatively, filling
#' the spectral data for wavelengths outside the overlapping range with with
#' \code{NA} when data is not available.
#'
#' @param x generic_mspct object, or an object of a class derived from
#'   \code{generic_mspct}.
#' @param unit.out character Allowed values \code{"energy"}, and \code{"photon"}, or its alias
#'   \code{"quantum"}.
#' @param qty.out character Allowed values \code{"transmittance"},
#'   \code{"absorptance"}, and \code{"absorbance"} and in the method for
#'   \code{object_spct}, also \code{"reflectance"} (.
#' @param type character Type of join: \code{"inner"} (default) or
#'   \code{"full"}. See details for more information.
#' @param validate.names logical A flag to enable (default) or disable
#'   validation of column names with \code{\link[base]{make.names}}.
#' @param ... ignored (possibly used by derived methods).
#'
#' @return A \code{data.frame} with the spectra joined by, possibly
#'   interpolated, wavelength, with rows sorted by wavelength (variable
#'   \code{w.length}) and data columns named according to the names of members
#'   in \code{x}, possibly made unique and valid.
#'
#' @note Currently only \code{generic_spct}, \code{source_mspct},
#'   \code{response_mspct}, \code{filter_mspct}, \code{reflector_mspct},
#'   \code{object_mspct} and \code{solute_mspct} classes have this method
#'   implemented.
#'
#' @export
#'
#' @examples
#'
#' join_mspct(solute_mspct(list(water = water.spct, pha = phenylalanine.spct)),
#'            type = "inner")
#'
#' @family conversion of collections of spectra
#'
join_mspct <- function(x, type, ...) UseMethod("join_mspct")

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.default <- function(x, type = "full", ...) {
  stop("'join_mspct()' is only implemented for some collections of spectra")
}

#' @describeIn join_mspct
#'
#' @param col.name character, name of the column in the spectra to be preserved,
#'   in addition to "w.length".
#'
#' @export
#'
join_mspct.generic_mspct <- function(x,
                                     type = "full",
                                     col.name,
                                     validate.names = TRUE,
                                     ...) {

  col.exists <- function(spct, col.name) {
    any(grepl(pattern = col.name, x = names(spct)))
  }

  if (length(x) == 0L) {
    return(data.frame("w.length" = numeric()))
  }

  if (!all(sapply(x, col.exists, col.name = col.name))) {
    stop("Variable '", col.name, "' not present in all spectra")
  }

  if (validate.names) {
    names(x) <- make.names(names(x), unique = TRUE)
  }

  spct.names <- names(x)

  if (length(x) == 1L) {
    z <- as.data.frame(x[[1]])[ , c("w.length", col.name)]
    colnames(z)[2] <- spct.names
    return(z)
  } else {
    # we need to check consistency of walengths
    wl.range <- c(photobiology::wl_range(x))
    wl.stepsize <- photobiology::wl_stepsize(x)
    wl.ranges.consistent <-
      length(unique(wl.range[["min.wl"]])) == 1 &&
      length(unique(wl.range[["max.wl"]])) == 1
    wl.stepsizes.consistent <- length(unique(wl.stepsize[["min.step.wl"]])) == 1 &&
      length(unique(wl.stepsize[["max.step.wl"]])) == 1

    if (!wl.ranges.consistent || !wl.stepsizes.consistent) {
      # overlapping range
      wl.range.inner <- c(max(wl.range[["min.wl"]]),
                          min(wl.range[["max.wl"]]))
      # full range
      wl.range.full <- c(min(wl.range[["min.wl"]]),
                         max(wl.range[["max.wl"]]))
      if (any(wl.range.inner != wl.range.full)) {
        if (type == "inner") {
          wl.range.out <- wl.range.inner
          message("Trimming non-overlapping wavelengths")
        } else if (type == "full") {
          wl.range.out <- wl.range.inner
          message("Filling non-overlapping wavelengths with NA")
        }
      }
      wl.stepsize.out <- stats::median(wl.stepsize[["min.step.wl"]]) / 2
      # we try to find a nearby "nice" stepsize
      wl.stepsize.out <- ifelse(wl.stepsize.out >= 1,
                                trunc(wl.stepsize.out),
                                ifelse(wl.stepsize.out >= 0.25,
                                       trunc(wl.stepsize.out * 4) / 4,
                                       round(wl.stepsize.out, digits = 2)))
      wl.out <- seq(from = wl.range.out[1],
                    to = wl.range.out[2],
                    by = wl.stepsize.out)
      x <- photobiology::interpolate_mspct(x,
                                           w.length.out = wl.out,
                                           fill = NA_real_)
      message("Spectra interpolated and/or trimmed as wavelengths differed.")
    } else {
      wl.out <- x[[1]][["w.length"]]
    }
  }
  rmDerivedMspct(x) # convert to list
  z <- list(w.length = wl.out)
  for (i in spct.names) {
    z[[i]] <- x[[i]][[col.name]]
  }
  as.data.frame(z)
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.source_mspct <- function(x,
                                    type = "full",
                                    unit.out = "energy",
                                    validate.names = TRUE,
                                    ...) {
  if (length(x)) {
    if (unit.out == "energy") {
      x <- q2e(x, action = "replace")
      col.name <- "s.e.irrad"
    } else if (unit.out %in% c("photon", "quantum")) {
      x <- e2q(x, action = "replace")
      col.name <- "s.q.irrad"
    } else {
      stop("Unit out '", unit.out, "' unknown")
    }
  }
  class(x) <- class(x)[-1L] # convert to generic_spct
  join_mspct(x,
             type = type,
             col.name = col.name,
             validate.names = validate.names,
             ...)
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.response_mspct <- function(x,
                                      type = "full",
                                      unit.out = "energy",
                                      validate.names = TRUE,
                                      ...) {
  if (length(x)) {
    if (unit.out == "energy") {
      x <- q2e(x, action = "replace")
      col.name <- "s.e.response"
    } else if (unit.out %in% c("photon", "quantum")) {
      x <- e2q(x, action = "replace")
      col.name <- "s.q.response"
    } else {
      stop("Unit out '", unit.out, "' unknown")
    }
  }
  class(x) <- class(x)[-1L] # convert to generic_spct
  join_mspct(x,
             type = type,
             col.name = col.name,
             validate.names = validate.names,
             ...)
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.filter_mspct <- function(x,
                                    type = "full",
                                    qty.out = "transmittance",
                                    validate.names = TRUE,
                                    ...) {
  if (length(x)) {
    if (qty.out == "transmittance") {
      x <- any2T(x, action = "replace")
      col.name <- "Tfr"
    } else if (qty.out == "absorbance") {
      x <- any2A(x, action = "replace")
      col.name <- "A"
    } else if (qty.out == "absorptance") {
      x <- any2Afr(x, action = "replace")
      col.name <- "Afr"
    } else {
      stop("Unit out '", qty.out, "' unknown")
    }
  }
  class(x) <- class(x)[-1L] # convert to generic_spct
  join_mspct(x,
             type = type,
             col.name = col.name,
             validate.names = validate.names,
             ...)
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.reflector_mspct <- function(x,
                                       type = "full",
                                       validate.names = TRUE,
                                       ...) {
  class(x) <- class(x)[-1L] # convert to generic_spct
  join_mspct(x,
             type = type,
             col.name = "Rfr",
             validate.names = validate.names,
             ...)
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.object_mspct <- function(x,
                                    type = "full",
                                    qty.out,
                                    validate.names = TRUE,
                                    ...) {
  switch(qty.out,
         "transmittance" = join_mspct(as.filter_mspct(x), type = type, qty.out = qty.out, ...),
         "absorbance" = join_mspct(as.filter_mspct(x), type = type, qty.out = qty.out, ...),
         "absorbtance" = join_mspct(as.filter_mspct(x), type = type, qty.out = qty.out, ...),
         "reflectance" = join_mspct(as.reflector_mspct(x), type = type, ...),
         stop("'qty.out = ", qty.out, " not implemented.")
  )
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.solute_mspct <- function(x,
                                    type = "full",
                                    validate.names = TRUE,
                                    ...) {
  # guess column name from 1st spectrum
  if (length(x)) {
    col.name <- intersect(c("K.mole", "K.mass"), names(x[[1]]))
  } else {
    col.name <- NA_character_
  }

  class(x) <- class(x)[-1L] # convert to generic_spct
  join_mspct(x,
             type = type,
             col.name = col.name,
             validate.names = validate.names,
             ...)
}

