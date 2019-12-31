# spct and mspct utility methods ----------------------------------------------

#' Extract all members from a collection
#'
#' Extract all members from a collection into separate objects in the parent
#' frame of the call.
#'
#' @param x An R object
#' @param ... additional named arguments passed down to \code{f}.
#'
#' @return Utility used for its side effects, invisibly returns a character
#'   vector with the names of the objects created.
#'
#' @export
#'
#' @examples
#'
#' my.mscpt <- source_mspct(list(sun1 = sun.spct, sun2 = sun.spct))
#' uncollect(my.mscpt)
#' ls(pattern = "*.spct")
#'
#' @family experimental utility functions
#'
uncollect <- function(x, ...) UseMethod("uncollect")

#' @describeIn uncollect Default for generic function
#'
#' @export
#'
uncollect.default <- function(x, ...) {
  warning("'uncollect()' is not defined for objects of class '", class(x)[1], "'.")
  invisible(character())
}

#' @describeIn uncollect
#'
#' @param name.tag character. A string used as tag for the names of the objects.
#'   If of length zero, names of members are used as named of objects. Otherwise
#'   the tag is appended, unless already present in the member name.
#' @param ignore.case	logical. If FALSE, the pattern matching used for \code{name.tag} is
#'   case sensitive and if TRUE, case is ignored during matching.
#' @param check.names logical. If TRUE then the names of the objects created are
#'   checked to ensure that they are syntactically valid variable names and
#'   unique. If necessary they are adjusted (by make.names) so that they are,
#'   and if FALSE names are used as is.
#' @param check.overwrite logical. If TRUE trigger an error if an exisitng
#'   object would be overwritten, and if FALSE silently overwrite objects.
#'
#' @export
#'
uncollect.generic_mspct <- function(x,
                                    name.tag = ".spct",
                                    ignore.case = FALSE,
                                    check.names = TRUE,
                                    check.overwrite = TRUE,
                                    ...) {
  stopifnot(rlang::is_named(x))
  # make object names
  if (length(name.tag) == 0L) {
    obj.names <- names(x)
  } else {
    obj.names <- character()
    for (member.name in names(x)) {
      if (!grepl(pattern = paste("*", name.tag, "$", sep = ""),
                 x = member.name,
                 ignore.case = ignore.case)) {
        obj.name <- paste(member.name, name.tag, sep = "")
      } else {
        obj.name <- member.name
      }
      obj.names <- c(obj.names, obj.name)
    }
  }
  # by checking a vector of object names we protect from repeated names
  if (check.names) {
    obj.names <- make.names(obj.names, unique = TRUE)
  }
  names(obj.names) <- names(x)

  # check for exisiting objects
  if (check.overwrite) {
    existing.names <-
      ls(envir = parent.frame(2), pattern = paste("*", name.tag, "$", sep = ""), sorted = FALSE)
    if (length(intersect(obj.names, existing.names)) > 0L) {
      stop("Overwriting would take place, if intended, change default.")
    }
  }

  # create objects
  for (member.name in names(x)) {
    assign(obj.names[member.name], x[[member.name]], envir = parent.frame(2))
  }
  invisible(unname(obj.names))
}

# thining of wavelengths --------------------------------------------------

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
#' thinned.spct <- wl_thin(yellow_gel.spct)
#' nrow(thinned.spct)
#' wl_stepsize(thinned.spct)
#'
#' @family experimental utility functions
#'
wl_thin <- function(x, ...) UseMethod("wl_thin")

#' @describeIn wl_thin Default for generic function
#'
#' @export
#'
wl_thin.default <- function(x, ...) {
  warning("'wl_thin()' is not defined for objects of class '", class(x)[1], "'.")
  x
}

#' @describeIn wl_thin
#'
#' @param max.wl.step numeric. Largest allowed wavelength difference between
#'    adjacent spectral values in nanometres (nm).
#' @param max.slope.delta numeric in 0 to 1. Largest allowed change in relative
#'    slope of the spectral quantity per nm betweem adjacent pairs of values.
#' @param col.names character. Name of the column of \code{x} containing the
#'   spectral data to check against \code{max.slope.delta}. Currently only one
#'   column supported.
#'
#' @note The value of \code{max.slope.delta} is expressed as relative change in
#'   the slope of spectral variable per nanometre. This means that values
#'   between 0.0005 and 0.005 tend to work reasonably well. The best value
#'   will depend on the wavelength step of the input and noise in data. A
#'   moderate smoothing before thinning can sometimes help in the case of
#'   noisy data.
#'   The amount of thinning is almost always less than the value of criteria
#'   passed as argument as it is based on existing wavelength values. For
#'   example if we start with a spectrum with a uniform wavelength step of 1 nm,
#'   possible steps in the thinned spectrum are 2, 4, 8, 16, 32, etc. nm. The
#'   algorithm, does work with any step sizes, regular or variable in the input.
#'   Thinning is most effective for spectra with large "featureless" regions as
#'   the algorithm attempts not to discard information, contrary to smoothing or
#'   interpolation.
#'
#' @export
#'
wl_thin.generic_spct <- function(x,
                                 max.wl.step = 10.0,
                                 max.slope.delta = 0.001,
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
  peaks <- find_peaks(x.norm[[col.names]], span = 21, strict = FALSE)
  valleys <- find_peaks(-x.norm[[col.names]], span = 21, strict = FALSE)
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

#' @describeIn wl_thin
#'
#' @param unit.out character Allowed values "energy", and "photon", or its alias
#'   "quantum".
#'
#' @export
#'
wl_thin.source_spct <- function(x,
                                max.wl.step = 10.0,
                                max.slope.delta = 0.001,
                                unit.out = getOption("photobiology.radiation.unit", default = "energy"),
                                ...) {
  if (unit.out == "energy") {
    wl_thin.generic_spct(x = q2e(x, action = "replace"),
                         max.wl.step = max.wl.step,
                         max.slope.delta = max.slope.delta,
                         col.names = "s.e.irrad",
                         ...)
  } else if (unit.out %in% c("photon", "quantum")) {
    wl_thin.generic_spct(x = e2q(x, action = "replace"),
                         max.wl.step = max.wl.step,
                         max.slope.delta = max.slope.delta,
                         col.names = "s.q.irrad",
                         ...)
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' @describeIn wl_thin
#'
#' @export
#'
wl_thin.response_spct <- function(x,
                                  max.wl.step = 10.0,
                                  max.slope.delta = 0.001,
                                  unit.out = getOption("photobiology.radiation.unit", default = "energy"),
                                  ...) {
  if (unit.out == "energy") {
    wl_thin.generic_spct(x = q2e(x, action = "replace"),
                         max.wl.step = max.wl.step,
                         max.slope.delta = max.slope.delta,
                         col.names = "s.e.response",
                         ...)
  } else if (unit.out %in% c("photon", "quantum")) {
    wl_thin.generic_spct(x = e2q(x, action = "replace"),
                         max.wl.step = max.wl.step,
                         max.slope.delta = max.slope.delta,
                         col.names = "s.q.response",
                         ...)
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' @describeIn wl_thin
#'
#' @param qty.out character Allowed values "transmittance", and "absorbance".
#'
#' @export
#'
wl_thin.filter_spct <- function(x,
                                max.wl.step = 10.0,
                                max.slope.delta = 0.001,
                                qty.out = getOption("photobiology.filter.qty",
                                                    default = "transmittance"),
                                ...) {
  if (qty.out == "transmittance") {
    wl_thin.generic_spct(x = A2T(x, action = "replace"),
                         max.wl.step = max.wl.step,
                         max.slope.delta = max.slope.delta,
                         col.names = "Tfr",
                         ...)
  } else if (qty.out == "absorbance") {
    wl_thin.generic_spct(x = T2A(x, action = "replace"),
                         max.wl.step = max.wl.step,
                         max.slope.delta = max.slope.delta,
                         col.names = "A",
                         ...)
  } else {
    stop("'unit.out ", qty.out, " is unknown")
  }
}

#' @describeIn wl_thin
#'
#' @export
#'
wl_thin.reflector_spct <- function(x,
                                   max.wl.step = 10.0,
                                   max.slope.delta = 0.001,
                                   ...) {
  wl_thin.generic_spct(x = x,
                       max.wl.step = max.wl.step,
                       max.slope.delta = max.slope.delta,
                       col.names = "Rfr",
                       ...)
}

#' @describeIn wl_thin
#'
#' @export
#'
wl_thin.raw_spct <- wl_thin.generic_spct

#' @describeIn wl_thin
#'
#' @export
#'
wl_thin.cps_spct <- wl_thin.generic_spct

#' @describeIn wl_thin
#'
#' @export
#'
wl_thin.object_spct <- wl_thin.generic_spct

#' @describeIn wl_thin
#'
#' @export
#'
wl_thin.chroma_spct <- wl_thin.default

#' @describeIn wl_thin
#'
#' @export
#'
wl_thin.calibration_spct <- wl_thin.default

#' @describeIn wl_thin
#'
#' @export
#'
wl_thin.generic_mspct <- function(x,
                                  max.wl.step = 10.0,
                                  max.slope.delta = 0.001,
                                  ...) {
  msmsply(x,
          wl_thin,
          max.wl.step = max.wl.step,
          max.slope.delta = max.slope.delta,
          ...)
}

#' @describeIn wl_thin
#'
#' @export
#'
wl_thin.chroma_mspct <- wl_thin.default

#' @describeIn wl_thin
#'
#' @export
#'
wl_thin.calibration_mspct <- wl_thin.default
