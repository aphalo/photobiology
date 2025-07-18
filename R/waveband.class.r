# labels ------------------------------------------------------------------

#' Find labels from "waveband" object
#'
#' A method specialization that extracts the name and label of objects of class
#' \code{waveband}.
#'
#' @param object an object of class "waveband"
#' @param ... not used in current version
#'
#' @export
#'
#' @name labels
#'
#' @family waveband attributes
#'
labels.waveband <- function(object, ...) {
  return(list(label = object[["label"]], name = object[["name"]]))
}

#' @describeIn labels
#'
#' @export
#'
#' @examples
#' labels(sun.spct)
#'
labels.generic_spct <- function(object, ...) {
  return(names(object))
}

# range -------------------------------------------------------------------

#' @rdname range
#'
#' @param x generic_spct, generic_mspct or waveband object.
#'
#' @export
#'
wl_range <- function(x, na.rm = FALSE) {
  stopifnot(is.any_spct(x) || is.any_mspct(x) || is.waveband(x))
  range(x, na.rm = na.rm)
}

#' Wavelength range
#'
#' A method specialization that returns the wavelength range [\eqn{nm}] from
#' objects of classes \code{waveband} or of class \code{generic_spct} or
#' derived.
#'
#' @param ... a single R object
#' @param na.rm ignored
#' @export
#'
#' @name range
#'
#' @return a length-two vector for individual objects or numeric vectors
#'   or a data frame for collections of spectra.
#'
#' @family wavelength summaries
#'
#' @examples
#' range(sun.spct)
#' wl_range(sun.spct)
#'
range.waveband <- function(..., na.rm = FALSE) {
  x <- c(...)
  return(c(x[["low"]], x[["high"]])) # we are using double precision
}

#' @describeIn range
#'
#' @export
#'
#' @examples
#' range(sun.spct)
#'
range.generic_spct <- function(..., na.rm = FALSE) {
  wl <- list(...)[[1]][["w.length"]]
  # guaranteed to be sorted
  wl[c(1, length(wl))]
  #  range(x[["w.length"]], na.rm = na.rm)
}

#' @describeIn range
#'
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#'
#' @export
#'
range.generic_mspct <- function(..., na.rm = FALSE, idx = "spct.idx") {
  mspct <- list(...)[[1]]
  if (is.null(idx)) {
    idx <- !is.null(names(mspct))
  }
  msdply(mspct = mspct, .fun = range, na.rm = na.rm, idx = idx)
}

# min ---------------------------------------------------------------------

#' @rdname min
#'
#' @param x generic_spct, generic_mspct or waveband object.
#'
#' @export
#'
wl_min <- function(x, na.rm = FALSE) {
  stopifnot(is.any_spct(x) || is.any_mspct(x) || is.waveband(x))
  min(x, na.rm = na.rm)
}

#' Wavelength minimum
#'
#' A method specialization that returns the wavelength minimum [\eqn{nm}] from
#' objects of classes \code{waveband} or of class \code{generic_spct} or
#' derived.
#'
#' @param ... not used in current version
#' @param na.rm ignored
#' @export
#'
#' @return a length-one vector for individual objects or numeric vectors
#'   or a data frame for collections of spectra.
#'
#' @name min
#'
#' @family wavelength summaries
#'
#' @examples
#' min(sun.spct)
#' wl_min(sun.spct)
#'
min.waveband <- function(..., na.rm = FALSE) {
  x <- c(...)
    return(x[["low"]])
}

#' @describeIn min
#'
#' @export
#'
min.generic_spct <- function(..., na.rm = FALSE) {
  wl <- list(...)[[1]][["w.length"]]
  # guaranteed to be sorted
  wl[1]
}

#' @describeIn min
#'
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#'
#' @export
#'
min.generic_mspct <- function(..., na.rm = FALSE, idx = "spct.idx") {
  mspct <- list(...)[[1]]
  if (is.null(idx)) {
    idx <- !is.null(names(mspct))
  }
  msdply(mspct = mspct, .fun = min, na.rm = na.rm, idx = idx)
}

# max ---------------------------------------------------------------------

#' @rdname max
#'
#' @param x generic_spct, generic_mspct or waveband object.
#'
#' @export
#'
wl_max <- function(x, na.rm = FALSE) {
  stopifnot(is.any_spct(x) || is.any_mspct(x) || is.waveband(x))
  max(x, na.rm = na.rm)
}

#' Wavelength maximum
#'
#' A method specialization that returns the wavelength maximum [\eqn{nm}] from
#' objects of classes \code{waveband} or of class \code{generic_spct} or
#' derived.
#'
#' @param ... numeric, \code{waveband} or \code{generic_spct} arguments.
#' @param na.rm ignored
#' @export
#'
#' @return a length-one vector for individual objects or numeric vectors
#'   or a data frame for collections of spectra.
#'
#' @name max
#'
#' @examples
#' max(sun.spct)
#' wl_max(sun.spct)
#'
max.waveband <- function(..., na.rm = FALSE) {
  x <- c(...)
  return(x[["high"]])
}

#' @describeIn max
#'
#' @export
#'
max.generic_spct <- function(..., na.rm=FALSE) {
  wl <- list(...)[[1]][["w.length"]]
  # guaranteed to be sorted
  wl[length(wl)]
}

#' @describeIn max
#'
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#'
#' @export
#'
max.generic_mspct <- function(..., na.rm = FALSE, idx = "spct.idx") {
  mspct <- list(...)[[1]]
  if (is.null(idx)) {
    idx <- !is.null(names(mspct))
  }
  msdply(mspct = mspct, .fun = max, ..., na.rm = na.rm, idx = idx)
}


# midpoint ------------------------------------------------------------------

#' @rdname midpoint
#'
#' @export
#'
wl_midpoint <- function(x, ...) {
  stopifnot(is.any_spct(x) || is.any_mspct(x) || is.waveband(x))
  midpoint(x, ...)
}

#' Midpoint
#'
#' A method that returns the wavelength [\eqn{nm}] (or value) at the center of
#' the wavelength range of objects of classes \code{waveband} or of class
#' \code{generic_spct} or derived (or the midpoint from a \code{numeric}
#' vector).
#'
#' @param x an R object
#' @param ... not used in current version
#' @export
#'
#' @return A numeric value equal to \code{max(x) - min(x)) / 2}. In the case of
#'   spectral objects a wavelength [\eqn{nm}]. For any other R object, according
#'   to available definitions of \code{\link{min}} and \code{\link{max}}.
#'
#' @family wavelength summaries
#'
#' @examples
#' midpoint(10:20)
#' midpoint(sun.spct)
#' wl_midpoint(sun.spct)
#'
midpoint <- function(x, ...) UseMethod("midpoint")

#' @describeIn midpoint Default method for generic function
#'
#' @export
#'
#' @family wavelength summaries
#'
midpoint.default <- function(x, ...) {
  warning("'midpoint()' not implemented for class '", class(x), "'.")
  NA_real_
}

#' @describeIn midpoint Default method for generic function
#'
#' @export
#'
#' @family wavelength summaries
#'
midpoint.numeric <- function(x, ...) {
  if (length(x) > 0) {
    min(x) + (max(x) - min(x)) / 2
  } else {
    NA_real_
  }
}

#' @describeIn midpoint Wavelength at center of a "waveband".
#'
#' @export
#'
midpoint.waveband <- function(x, ...) {
  return(x[["low"]] + (x[["high"]] - x[["low"]]) / 2)
}

#' @describeIn midpoint Method for "generic_spct".
#'
#' @export
#'
#' @examples
#' midpoint(sun.spct)
#'
midpoint.generic_spct <- function(x, ...) {
  wl <- x[["w.length"]]
  wl[1] + (wl[length(wl)] - wl[1]) / 2
}

#' @describeIn midpoint Method for "generic_mspct" objects.
#'
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#'
#' @export
#'
midpoint.generic_mspct <- function(x, ..., idx = "spct.idx") {
  if (is.null(idx)) {
    idx <- !is.null(names(x))
  }
  msdply(mspct = x, .fun = midpoint, ..., idx = idx)
}

# expanse ------------------------------------------------------------------

#' @rdname expanse
#'
#' @export
#'
spread <- function(x, ...) {
  message("Use of method photobiology::spread() is deprecated. It has been ",
          "renamed into expanse() to avoid a name clash with 'tidyr::spread()'.")
  expanse(x, ...)
}

#' @rdname expanse
#'
#' @export
#'
wl_expanse <- function(x, ...) {
  stopifnot(is.any_spct(x) || is.any_mspct(x) || is.waveband(x))
  expanse(x, ...)
}

#' Expanse
#'
#' A method that returns the expanse (\eqn{max(x) - min(x)}) for R objects. In
#' particular the wavelength [\eqn{nm}] expanse of the wavelength range of
#' objects of classes \code{waveband} or of class \code{generic_spct} or derived
#' (or the expanse of values in a \code{numeric} vector).
#'
#' @param x an R object
#' @param ... not used in current version
#'
#' @return A numeric value equal to \code{max(x) - min(x)}. In the case of
#'   spectral objects wavelength difference [\eqn{nm}]. For any other R object,
#'   according to available specialised methods of \code{\link{min}} and
#'   \code{\link{max}}.
#'
#' @export
#'
#' @examples
#' expanse(10:20)
#' expanse(sun.spct)
#' wl_expanse(sun.spct)
#'
expanse <- function(x, ...) UseMethod("expanse")

#' @describeIn expanse Default method for generic function
#'
#' @export
#'
expanse.default <- function(x, ...) {
  warning("'expanse()' not defined for class '", paste(class(x), collapse = " "), "'")
  NA
}

#' @describeIn expanse Method for "numeric"
#'
#' @export
#'
expanse.numeric <- function(x, ...) {
  if (length(x) > 0) {
    return(max(x) - min(x))
  } else {
    return(NA_real_)
  }
}

#' @describeIn expanse Method for "waveband"
#'
#' @export
#'
expanse.waveband <- function(x, ...) {
  return(x[["high"]] - x[["low"]])
}

#' @describeIn expanse  Method for "generic_spct"
#'
#' @export
#'
#' @examples
#' expanse(sun.spct)
#'
expanse.generic_spct <- function(x, ...) {
  wl <- x[["w.length"]]
  wl[length(wl)] - wl[1]
}

#' @describeIn expanse  Method for "generic_mspct" objects.
#'
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#'
#' @export
#'
expanse.generic_mspct <- function(x, ..., idx = "spct.idx") {
  if (is.null(idx)) {
    idx <- !is.null(names(x))
  }
  msdply(mspct = x, .fun = expanse, ..., idx = idx)
}

# normalization -----------------------------------------------------------

#' Normalization of an R object
#'
#' Normalization wavelength [\eqn{nm}] and other normalization metadata of an R
#' object, retrieved from the object's attributes.
#'
#' @details In the case of wavebands for spectral weighting functions
#'   (\code{waveband} objects), the normalization wavelength is returned. For
#'   spectral objects (\code{generic_spct} and derived ), the normalization
#'   descriptor, a list object, is returned. This list contains in addition to
#'   the normalization wavelength, the multiplier used and type of normalization
#'   applied. These metadata makes it possible to "undo" the normalization and
#'   to "update" the normalization after a transformation, such as conversion to
#'   a related physical quantity, of the spectral data.
#'
#' @param x an R object
#' @export
#'
#' @return A single numeric value of wavelength [\eqn{nm}] or a list with
#'   with members.
#'
#' @note Older versions of the package stored only a subset of the metadata or
#'   only a flag to indicate that normalization had been applied. For such
#'   objects some or even all fields in the returned list are set to \code{NA}.
#'
#' @examples
#'
#' is_normalized(sun.spct)
#' normalization(sun.spct)
#' sun_norm.spct <- normalize(sun.spct)
#' is_normalized(sun_norm.spct)
#' normalization(sun_norm.spct)
#'
#' my_wband <- waveband(c(400,700))
#' is_normalized(my_wband)
#' normalization(my_wband)
#'
#' @family waveband attributes
#'
normalization <- function(x) UseMethod("normalization")

#' @describeIn normalization Default methods.
#'
#' @export
#'
normalization.default <- function(x) {
  warning("'normalization()' not implemented for class '", class(x)[1], "'.")
  return(NA_real_)
}

#' @describeIn normalization Normalization of a \code{\link{waveband}} object.
#'
#' @export
#'
normalization.waveband <- function(x) {
  return(ifelse(is.null(x[["norm"]]), NA_real_, x[["norm"]]))
}

#' @describeIn normalization Normalization of a \code{\link{generic_spct}} object.
#'
#' @export
#'
normalization.generic_spct <- function(x) {
  getNormalization(x)
}

#' @describeIn normalization Normalization of a \code{\link{summary.generic_spct}} object.
#'
#' @export
#'
normalization.summary_generic_spct <- function(x) {
  getNormalization(x)
}

#' @describeIn normalization Normalization of a \code{\link{generic_mspct}} object.
#'
#' @export
#'
normalization.generic_mspct <- function(x) {
  getNormalization(x)
}

# is_effective -----------------------------------------------------------

#' Is an R object "effective"
#'
#' A generic function for querying if a biological spectral weighting function
#' (BSWF) has been applied to an object or is included in its definition.
#'
#' @param x an R object
#'
#' @return A \code{logical}.
#'
#' @export
#'
#' @family waveband attributes
#'
is_effective <- function(x) UseMethod("is_effective")

#' @describeIn is_effective Default method.
#'
#' @export
#'
is_effective.default <- function(x) {
  warning("'is_effective()' not implemented for class '", class(x), "'.")
  NA_integer_
}

#' @describeIn is_effective Is a \code{waveband} object defining a method for
#'   calculating effective irradiance.
#'
#' @export
#'
is_effective.waveband <- function(x) {
  x[["weight"]] != "none"
}

#' @describeIn is_effective Does a \code{source_spct} object contain effective
#'   spectral irradiance values.
#'
#' @export
#'
is_effective.generic_spct <- function(x) {
  FALSE
}

#' @describeIn is_effective Does a \code{source_spct} object contain effective
#'   spectral irradiance values.
#'
#' @export
#'
is_effective.source_spct <- function(x) {
  bswf.used <- getBSWFUsed(x)
  !is.null(bswf.used) && (bswf.used != "none")
}

#' @describeIn is_effective Method for "summary_generic_spct".
#'
#' @export
#' @examples
#' is_effective(summary(sun.spct))
#'
is_effective.summary_generic_spct <- function(x) {
  FALSE
}

#' @describeIn is_effective Method for "summary_source_spct".
#'
#' @export
#'
is_effective.summary_source_spct <- function(x) {
  bswf.used <- getBSWFUsed(x)
  !is.null(bswf.used) && (bswf.used != "none")
}

# w.length summaries ------------------------------------------------------

#' @rdname stepsize
#'
#' @export
#'
wl_stepsize <- function(x, ...) {
  stopifnot(is.any_spct(x) || is.any_mspct(x))
  stepsize(x, ...)
}

#' Stepsize
#'
#' Method returning the range of step sizes in an object; i.e., the Range of
#' differences between successive sorted values.  In particular the wavelength
#' step sizes [\eqn{nm}] of objects of class \code{generic_spct} or derived (or
#' the step sizes of values in a \code{numeric} vector).
#'
#' @param x an R object
#' @param ... not used in current version
#'
#' @return A numeric vector of length 2 with min and maximum stepsize values.
#' @export
#' @family wavelength summaries
#'
#' @examples
#' stepsize(sun.spct)
#' wl_stepsize(sun.spct)
#'
stepsize <- function(x, ...) UseMethod("stepsize")

#' @describeIn stepsize Default function usable on numeric vectors.
#' @export
stepsize.default <- function(x, ...) {
  warning("'stepsize()' not implemented for class '", class(x), "'.")
  c(NA_real_, NA_real_)
}

#' @describeIn stepsize Method for numeric vectors.
#' @export
stepsize.numeric <- function(x, ...) {
  stopifnot(!is.unsorted(x))
  if (length(x) > 1) {
    range(diff(x))
  } else {
    c(NA_real_, NA_real_)
  }
}

#' @describeIn stepsize  Method for "generic_spct" objects.
#'
#' @export
#'
#' @examples
#' stepsize(sun.spct)
#'
stepsize.generic_spct <- function(x, ...) {
  num.spectra <- getMultipleWl(x)
  if (num.spectra > 1) {
    wl <- unique(x[["w.length"]])
  } else {
    wl <- x[["w.length"]]
  }
  if (length(wl) > 1) {
    range(diff(wl))
  } else {
    c(NA_real_, NA_real_)
  }
}

#' @describeIn stepsize  Method for "generic_mspct" objects.
#'
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#'
#' @export
#'
stepsize.generic_mspct <- function(x, ..., idx = "spct.idx") {
  if (is.null(idx)) {
    idx <- !is.null(names(x))
  }
  msdply(mspct = x, .fun = stepsize, ..., idx = idx)
}
