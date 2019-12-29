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
#' Increase the wavelength step in stored spectral data in featureless
#' regions to save storage space.
#'
#' @param x An R object
#' @param ... additional named arguments passed down to \code{f}.
#'
#' @return An object of the same class as \code{x} but with a reduced density
#' of wavelength values in those regions were slope is shallow and featureless.
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
wl_thin <- function(x, ...) UseMethod("wl_thin")

#' @describeIn wl_thin Default for generic function
#'
#' @export
#'
wl_thin.default <- function(x, ...) {
  warning("'wl_thin()' is not defined for objects of class '", class(x)[1], "'.")
  invisible(character())
}

#' @describeIn wl_thin
#'
#' @param max.wl.step numeric. Largest allowed wavelength difference between
#'    adjacent spectral values.
#' @param delta.var.max numeric in 0 to 1. Largest allowed relative difference
#'    in spectral quantity betweem adjacent spectral values.
#' @param variable character. Name of the variable conatining the spectral
#'    data.
#'
#' @note Wavelengths are always stored expressed in nanometres in variable
#'    \code{"w.length"}.
#'
#' @export
#'
wl_thin.generic_spct <- function(x,
                                 max.wl.step = 2.0,
                                 delta.var.max = 0.01,
                                 variable,
                                 ...) {
  wl.stepsize <- wl_stepsize(x)
  # make code simpler as range is fixed
  x.norm <- normalize(x)
  # collect peaks and valleys to ensure that they are not removed
  peaks <- find_peaks(x.norm[[variable]], span = 21, strict = FALSE)
  valleys <- find_peaks(-x.norm[[variable]], span = 21, strict = FALSE)
  extremes <- peaks | valleys

  # compute differences as proxi for slope
  # using all values
  wl.step1 <- x.norm[["w.length"]]
  var.step1 <- x.norm[[variable]]
  diff.wl.step1 <- diff(wl.step1)
  diff.var.step1 <- diff(var.step1)
  local.slope.step1 <- diff.var.step1 / diff.wl.step1
  # using every other value
  wl.step2 <- wl.step1[c(TRUE, FALSE)]
  var.step2 <- var.step1[c(TRUE, FALSE)]
  diff.wl.step2 <- diff(wl.step2)
  diff.var.step2 <- diff(var.step2)
  local.slope.step2 <- diff.var.step2 / diff.wl.step2
  # using every fourth value
  wl.step4 <- wl.step2[c(TRUE, FALSE)]
  var.step4 <- var.step2[c(TRUE, FALSE)]
  diff.wl.step4 <- diff(wl.step4)
  diff.var.step4 <- diff(var.step4)
  local.slope.step4 <- diff.var.step4 / diff.wl.step4
  # using every eights value
  wl.step8 <- wl.step4[c(TRUE, FALSE)]
  # check var and wl deltas
  # if TRUE we cannot discard intermediate skipped values
  step1.test <- (abs(local.slope.step1) > delta.var.max) | (diff.wl.step1 > max.wl.step)
  step2.test <- (abs(local.slope.step2) > delta.var.max) | (diff.wl.step2 > max.wl.step)
  step4.test <- (abs(local.slope.step4) > delta.var.max) | (diff.wl.step4 > max.wl.step)
  # decide wavelengths to drop
  wls.to.keep <- unique(c(wl.step1[c(1, length(wl.step1))],
                          wl.step1[extremes],
                          wl.step1[-1][step1.test],
                          wl.step2[-1][step2.test],
                          wl.step4[-1][step4.test],
                          wl.step8[TRUE]))
  x[x[["w.length"]] %in% wls.to.keep, ]
}

