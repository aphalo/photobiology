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
#' my.mscpt <- source_mspct(list(sun1.spct = sun.spct, sun2.spct = sun.spct))
#' uncollect2spct(my.mscpt)
#' ls(pattern = "*.spct")
#'
#' @family experimental utility functions
#'
uncollect2spct <- function(x, ...) UseMethod("uncollect2spct")

#' @describeIn uncollect2spct Default for generic function
#'
#' @export
#'
uncollect2spct.default <- function(x, ...) {
  warning("'uncollect2spct()' is not defined for objects of class '", class(x)[1], "'.")
  invisible(character())
}

#' @describeIn uncollect2spct
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
uncollect2spct.generic_mspct <- function(x,
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

#' Form a new collection
#'
#' Form a collection of spectra from separate objects in the parent
#' frame of the call.
#'
#' @details This is a convenience function that simplifies the creation of
#'   collections from existing objects of class \code{generic_spct} or a derived
#'   class. A list of objects con be passed as argument, or a search pattern. If
#'   a list is passed, no search is done. If \code{collection.class} is
#'   \code{NULL}, then all objects of class \code{generic_spct} or of a class
#'   derived from it are added to the collection. If objects of only one derived
#'   class are to be collected this class or that of the matching collection
#'   should be passed as argument to \code{collection.class}. Objects of other R
#'   classes are silently discarded, which simplifies the specification of
#'   search patterns. By default, i.e., if \code{collection.class} is
#'   \code{NULL}, if all the objects collected belong to the same class then the
#'   corresponding collection class will be returned, otherwise a
#'   \code{generic_mspct} object with heterogeneous members will be returned. To
#'   force the return of a \code{generic_mspct} even when the collected spectra
#'   all belong to the same class, pass \code{generic_mspct} as argument to
#'   \code{collection.class}. If the argument to \code{collection.class} is a
#'   vector containing two of more class names, only the matching spectra will
#'   be collected, and a \code{generic_mspct} will be returned. The returned
#'   object is created with the constructor for the class, and validated.
#'
#' @param .list list of R objects
#' @param pattern character an optional regular expression, ignored if
#'   \code{.list} is not \code{NULL}.
#' @param collection.class character vector
#' @param ... additional named arguments passed down to the collection
#'   constructor.
#'
#' @return By default a collection of spectra.
#'
#' @export
#'
#' @examples
#' collect2mspct() # returns empty generic_mspct object
#'
#' sun1.spct <- sun.spct
#' sun2.spct <- sun.spct
#' kk.spct <- 10:30 # ignored
#' collect2mspct()
#' collect2mspct(collection.class = "generic_mspct")
#'
#' pet1.spct <- polyester.spct
#' collect2mspct()
#' collect2mspct(collection.class = "source_mspct")
#' collect2mspct(collection.class = "filter_mspct")
#' collect2mspct(collection.class = "response_mspct")
#'
#' @family experimental utility functions
#'
collect2mspct <- function(.list = NULL,
                          pattern = "*\\.spct$",
                          collection.class = NULL,
                          ...)  {
  collection.class <- gsub("_spct$", "_mspct", collection.class) # NULL -> character(0)!!
  if (length(collection.class) > 0L &&
      !all(collection.class %in% mspct_classes())) {
    warning("Discarding unrecognized class(es) in 'collection.class'")
    collection.class <- intersect(collection.class, spct_classes())
  }
  if (length(.list) == 0L) {
    names <- ls(pattern = pattern, envir = parent.frame())
    .list <- mget(x = names, envir = parent.frame())
  }
  if (length(.list) >= 1L) {
    if (length(collection.class) == 0L) {
      # we keep objects of class generic_spct or derived
      .list <- .list[sapply(.list, is.any_spct)]
      members.class <- shared_member_class(.list)[1] # derived class may be used
      collection.class <- gsub("_spct$", "_mspct", members.class)
    } else {
      # we keep objects of the classes in collection.class
      members.class <- gsub("_mspct$", "_spct", collection.class)
      if ("generic_spct" %in% members.class) {
        # we keep objects of class generic_spct or derived
        .list <-
          .list[sapply(X = .list, FUN = is.any_spct)]
        collection.class <- "generic_mspct" # forced
      } else {
        # we keep objects of classes listed in collection.class
        .list <-
          .list[sapply(.list, function(x) {class(x)[1]}) %in% members.class]
        if (length(collection.class) > 1L) {
          collection.class <- "generic_mspct"
        }
      }
    }
    # call the constructor for the class
    do.call(what = collection.class, args = list(l = .list))
  } else if (length(collection.class) == 1L) {
    do.call(what = collection.class, args = list())
  } else {
    generic_mspct()
  }
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
#' thinned.spct <- thin_wl(yellow_gel.spct)
#' nrow(thinned.spct)
#' wl_stepsize(thinned.spct)
#'
#' @family experimental utility functions
#'
thin_wl <- function(x, ...) UseMethod("thin_wl")

#' @describeIn thin_wl Default for generic function
#'
#' @export
#'
thin_wl.default <- function(x, ...) {
  warning("'thin_wl()' is not defined for objects of class '", class(x)[1], "'.")
  x
}

#' @describeIn thin_wl
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
                                 max.wl.step = 10.0,
                                 max.slope.delta = 0.001,
                                 span = 21,
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

#' @describeIn thin_wl
#'
#' @param unit.out character Allowed values "energy", and "photon", or its alias
#'   "quantum".
#'
#' @export
#'
thin_wl.source_spct <- function(x,
                                max.wl.step = 10.0,
                                max.slope.delta = 0.001,
                                span = 21,
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

#' @describeIn thin_wl
#'
#' @export
#'
thin_wl.response_spct <- function(x,
                                  max.wl.step = 10.0,
                                  max.slope.delta = 0.001,
                                  span = 21,
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

#' @describeIn thin_wl
#'
#' @param qty.out character Allowed values "transmittance", and "absorbance".
#'
#' @export
#'
thin_wl.filter_spct <- function(x,
                                max.wl.step = 10.0,
                                max.slope.delta = 0.001,
                                span = 21,
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

#' @describeIn thin_wl
#'
#' @export
#'
thin_wl.reflector_spct <- function(x,
                                   max.wl.step = 10.0,
                                   max.slope.delta = 0.001,
                                   span = 21,
                                   ...) {
  thin_wl.generic_spct(x = x,
                       max.wl.step = max.wl.step,
                       max.slope.delta = max.slope.delta,
                       span = span,
                       col.names = "Rfr",
                       ...)
}

#' @describeIn thin_wl
#'
#' @export
#'
thin_wl.solute_spct <- function(x,
                                max.wl.step = 10.0,
                                max.slope.delta = 0.001,
                                span = 21,
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

#' @describeIn thin_wl
#'
#' @export
#'
thin_wl.raw_spct <- thin_wl.generic_spct

#' @describeIn thin_wl
#'
#' @export
#'
thin_wl.cps_spct <- thin_wl.generic_spct

#' @describeIn thin_wl
#'
#' @export
#'
thin_wl.object_spct <- thin_wl.generic_spct

#' @describeIn thin_wl
#'
#' @export
#'
thin_wl.chroma_spct <- thin_wl.default

#' @describeIn thin_wl
#'
#' @export
#'
thin_wl.calibration_spct <- thin_wl.default

#' @describeIn thin_wl
#'
#' @export
#'
thin_wl.generic_mspct <- function(x,
                                  max.wl.step = 10.0,
                                  max.slope.delta = 0.001,
                                  span = 21,
                                  ...) {
  msmsply(x,
          thin_wl,
          max.wl.step = max.wl.step,
          max.slope.delta = max.slope.delta,
          span = span,
          ...)
}

#' @describeIn thin_wl
#'
#' @export
#'
thin_wl.chroma_mspct <- thin_wl.default

#' @describeIn thin_wl
#'
#' @export
#'
thin_wl.calibration_mspct <- thin_wl.default

# drop user columns -------------------------------------------------------

#' Drop user columns
#'
#' Remove from spectral object additional columns that are user defined.
#'
#' @param x An R object
#' @param keep.also character Additionlal columns to preserve.
#' @param ... needed to allow derivation.
#'
#' @return A copy of \code{x} possibly with some columns removed.
#'
#' @export
#'
#' @family experimental utility functions
#'
drop_user_cols <- function(x, keep.also, ...) UseMethod("drop_user_cols")

#' @describeIn drop_user_cols
#'
#' @export
#'
drop_user_cols.default <- function(x, keep.also = NULL, ...) {
  warning("'drop_user_cols()' is not defined for objects of class '", class(x)[1], "'.")
  x
}

#' @describeIn drop_user_cols
#'
#' @export
#'
drop_user_cols.generic_spct <- function(x, keep.also, ...) {
  stopifnot(length(keep.also) >= 1L)
  default.cols <- c("w.length")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @describeIn drop_user_cols
#'
#' @export
#'
drop_user_cols.source_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "s.e.irrad", "s.q.irrad")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @describeIn drop_user_cols
#'
#' @export
#'
drop_user_cols.response_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "s.e.response", "s.q.response")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @describeIn drop_user_cols
#'
#' @export
#'
drop_user_cols.object_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "Tfr", "Rfr", "Afr")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @describeIn drop_user_cols
#'
#' @export
#'
drop_user_cols.filter_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "Tfr", "Afr", "A")
  if ("Rfr" %in% colnames(x)) {
    warning("Deleting 'object_spct' columns from 'filter_spct'.")
  }
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @describeIn drop_user_cols
#'
#' @export
#'
drop_user_cols.reflector_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "Rfr")
  if (any(c("Tfr", "Afr", "A") %in% colnames(x))) {
    warning("Deleting 'object_spct' columns from 'reflector_spct'.")
  }
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @describeIn drop_user_cols
#'
#' @export
#'
drop_user_cols.solute_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "K.mole", "K.mass")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @describeIn drop_user_cols
#'
#' @export
#'
drop_user_cols.chroma_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "x", "y", "z")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @describeIn drop_user_cols
#'
#' @export
#'
drop_user_cols.calibration_spct <- function(x, keep.also = NULL, ...) {
  default.cols <- c("w.length", "irrad.mult")
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @describeIn drop_user_cols
#'
#' @export
#'
drop_user_cols.cps_spct <- function(x, keep.also = NULL, ...) {
  cps.cols <- grep("^cps", colnames(x), value = TRUE)
  default.cols <- c("w.length", cps.cols)
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @describeIn drop_user_cols
#'
#' @export
#'
drop_user_cols.raw_spct <- function(x, keep.also = NULL, ...) {
  counts.cols <- grep("^counts", colnames(x), value = TRUE)
  default.cols <- c("w.length", counts.cols)
  cols.to.keep <- unique(c(default.cols, keep.also))
  selector <- colnames(x) %in% cols.to.keep
  x[ , selector]
}

#' @describeIn drop_user_cols
#'
#' @export
#'
drop_user_cols.generic_mspct <-
  function(x, keep.also = NULL, ...) {
    msmsply(x, drop_user_cols, keep.also = keep.also, ...)
}

