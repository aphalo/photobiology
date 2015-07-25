
# names of all spectral classes -------------------------------------------

#' Function that returns a vector containing the names of spectra classes.
#'
#' @export
#'
#' @return A \code{character} vector of class names.
#'
spct_classes <- function() {
  c("cps_spct",
    "filter_spct", "reflector_spct",
    "source_spct", "object_spct",
    "response_spct", "chroma_spct", "generic_spct")
}

# check -------------------------------------------------------------------

#' Check validity of spectral objects
#'
#' Check that an R object contains the expected data members.
#'
#' @param x An R object
#' @param byref logical indicating if new object will be created by reference or
#'   by copy of \code{x}
#' @param strict.range logical indicating whether off-range values result in an
#'   error instead of a warning
#' @param ... additional param possible derived methods
#' @export
#'
#' @family data validity check functions
#'
check <- function(x, byref, strict.range, ...) UseMethod("check")

#' @describeIn check Default for generic function.
#' @export
check.default <- function(x, byref=FALSE, strict.range=TRUE, ...) {
  return(x)
}

#' @describeIn check Specialization for generic_spct.
#'
#' @param multiple.wl numeric Maximum number of repeated w.length entries with same value.
#'
#' @export
check.generic_spct <- function(x, byref=TRUE, strict.range=TRUE, multiple.wl = 1L, ...) {
  if (exists("w.length", x, mode = "numeric", inherits=FALSE)) {
    if (is.unsorted(x[["w.lengtgh"]], na.rm = TRUE, strictly = TRUE)) {
      stop("'w.length' must be sorted in ascending order and have unique values")
    }
  } else if (exists("wl", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "wl", "w.length")
  } else if (exists("wavelength", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "wavelength", "w.length")
  } else if (exists("Wavelength", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "Wavelength", "w.length")
  } else {
    warning("No wavelength data found in generic_spct")
    x[ , w.length := NA]
  }
  wl.min <- min(x[["w.length"]], na.rm = TRUE)
#  wl.max <- max(x$w.length, na.rm = TRUE)
  if (wl.min == Inf) {
    warning("No valid 'w.length' values, probably a spectrum of length zero")
  } else if (wl.min < 99.999 || wl.min > 5e3) {
    stop("Off-range minimum w.length value ", wl.min, " instead of within 100 nm and 5000 nm")
  }
  wl.reps <- with(x, length(w.length) / length(unique(w.length)) )
  if (wl.reps > 1) {
    if ((wl.reps - trunc(wl.reps)) < 1e-5) {
      wl.reps <- trunc(wl.reps)
      if (with(x, length(w.length)) >= 2 && wl.reps > multiple.wl) {
        warning("'w.length' values are not unique: ", wl.reps, " copies of each.")
      }
    } else {
      warning("'w.length' values are not all unique with ",
              round((wl.reps - trunc(wl.reps)) * with(x, length(w.length)), 0) ,
              " duplicate values.")
    }
  }
  return(x)
}

#' @describeIn check Specialization for cps_spct.
#' @export
check.cps_spct <- function(x, byref=TRUE, strict.range = TRUE, ...) {

  range_check <- function(x, strict.range) {
    cps.min <- min(x$cps, na.rm = TRUE)
    if (!is.null(strict.range) & (cps.min < 0)) {
      message.text <- paste0("Off-range cps values:", signif(cps.min, 2))
      if (strict.range) {
        stop(message.text)
      } else {
        warning(message.text)
      }
    }
  }
  if (exists("cps", x, mode = "numeric", inherits=FALSE)) {
    return(x)
  } else if (exists("counts.per.second", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "counts.per.second", "cps")
    warning("Found variable 'counts.per.second', renamed to 'cps'")
    return(x)
  } else {
    warning("No counts per second data found in cps_spct")
    x[["cps"]] = NA
    return(x)
  }
}

#' @describeIn check Specialization for filter_spct.
#' @export
check.filter_spct <- function(x, byref=TRUE, strict.range = TRUE, multiple.wl = 1L, ...) {

  range_check <- function(x, strict.range) {
    Tfr.min <- min(x[["Tfr"]], na.rm = TRUE)
    Tfr.max <- max(x[["Tfr"]], na.rm = TRUE)
    if (!is.null(strict.range) & (Tfr.min < 0 || Tfr.max > 1)) {
      message.text <- paste("Off-range transmittance values [", signif(Tfr.min, 2),
                            "...", signif(Tfr.max, 2), "] instead of  [0..1]", sep="")
      if (strict.range) {
        stop(message.text)
      } else {
        warning(message.text)
      }
    }
  }

  if (is.null(getTfrType(x))) {
    setTfrType(x, "total")
    warning("Missing Tfr.type attribute replaced by 'total'")
  }
  # check and replace 'other' quantity names
  if (exists("transmittance", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "transmittance", "Tpc")
    warning("Found varaible 'transmittance', I am assuming it is expressed as percent")
  }
  if (exists("absorbance", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "absorbance", "A")
    warning("Found varaible 'absorbance', I am assuming it is in log10-based absorbance units")
  }
  # look for percentages and change them into fractions of one
  if (exists("Tfr", x, mode = "numeric", inherits=FALSE)) {
    range_check(x, strict.range=strict.range)
    return(x)
  } else if (exists("Tpc", x, mode = "numeric", inherits=FALSE)) {
    x[["Tfr"]] <- x[["Tpc"]] / 100
    x[["Tpc"]] <-  NULL
    range_check(x, strict.range=strict.range)
    return(x)
  } else if (exists("A", x, mode = "numeric", inherits=FALSE)) {
#    x[ , Tfr := A2T(A)]
    if (min(x$A, na.rm = TRUE) < 0) {
      warning("Off-range min absorbance value: ", signif(min(x$A, na.rm = TRUE), 2), " instead of 0")
    }
    return(x)
  } else {
    warning("No transmittance or absorbance data found in filter_spct")
    x[["Tfr"]] <- NA
    return(x)
  }
}

#' @describeIn check Specialization for reflector_spct.
#' @export
check.reflector_spct <- function(x, byref=TRUE, strict.range = TRUE, ...) {

  range_check <- function(x, strict.range) {
    Rfr.min <- min(x$Rfr, na.rm = TRUE)
    Rfr.max <- max(x$Rfr, na.rm = TRUE)
    if (!is.null(strict.range) & (Rfr.min < 0 ||  Rfr.max > 1)) {
      message.text <- paste0("Off-range reflectance values [", signif(Rfr.min, 2), "...",
                             signif(Rfr.max, 2), "] instead of  [0..1]", sep="")
      if (strict.range) {
        stop(message.text)
      } else {
        warning(message.text)
      }
    }
  }

  if (is.null(getRfrType(x))) {
    setRfrType(x, "total")
    warning("Missing Rfr.type attribute replaced by 'total'")
  }
  if (exists("reflectance", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "reflectance", "Rpc")
    warning("Found variable 'reflectance', I am assuming it is expressed as percent")
  }
  if (exists("Rfr", x, mode = "numeric", inherits=FALSE)) {
    range_check(x, strict.range=strict.range)
    return(x)
  } else if (exists("Rpc", x, mode = "numeric", inherits=FALSE)) {
    x[["Rfr"]] <- x[["Rpc"]] / 100
    x[["Rpc"]] <- NULL
    range_check(x, strict.range=strict.range)
    return(x)
  } else {
    warning("No reflectance data found in reflector_spct")
    x[["Rfr"]] <- NA
    return(x)
  }
}

#' @describeIn check Specialization for object_spct.
#' @export
check.object_spct <- function(x, byref=TRUE, strict.range = TRUE, multiple.wl = 1L, ...) {

  range_check <- function(x, strict.range) {
    Rfr.min <- min(x[["Rfr"]], na.rm = TRUE)
    Rfr.max <- max(x[["Rfr"]], na.rm = TRUE)
    if (!is.na(Rfr.min) && !is.na(Rfr.max)) {
      if (!is.null(strict.range) & (Rfr.min < 0 ||  Rfr.max > 1)) {
        message.text <- paste0("Off-range reflectance values [", signif(Rfr.min, 2), "...",
                               signif(Rfr.max, 2), "] instead of  [0..1]", sep="")
        if (strict.range) {
          stop(message.text)
        } else {
          warning(message.text)
        }
      }
    }
    Tfr.min <- min(x[["Tfr"]], na.rm = TRUE)
    Tfr.max <- max(x[["Tfr"]], na.rm = TRUE)
    if (!is.na(Tfr.min) && !is.na(Tfr.max)) {
      if (!is.null(strict.range) & (Tfr.min < 0 ||  Tfr.max > 1)) {
      message.text <- paste0("Off-range Transmittance values [", signif(Tfr.min, 2), "...",
                             signif(Tfr.max, 2), "] instead of  [0..1]", sep="")
      if (strict.range) {
        stop(message.text)
      } else {
        warning(message.text)
      }
    }
    }
  }

  if (is.null(getTfrType(x))) {
    setTfrType(x, "total")
    warning("Missing Tfr.type attribute replaced by 'total'")
  }
  if (is.null(getRfrType(x))) {
    setRfrType(x, "total")
    warning("Missing Rfr.type attribute replaced by 'total'")
  }
  if (exists("reflectance", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "reflectance", "Rpc")
    warning("Found variable 'reflectance', I am assuming it is expressed as percent")
  }
  if (exists("Rfr", x, mode = "numeric", inherits=FALSE)) {
  } else if (exists("Rpc", x, mode = "numeric", inherits=FALSE)) {
    x[["Rfr"]] <- x[["Rpc"]] / 100
    x[["Rpc"]] <- NULL
  } else {
    warning("No reflectance data found in object_spct")
    x[["Rfr"]] <- NA
  }

  if (exists("transmittance", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "transmittance", "Tpc")
    warning("Found varaible 'transmittance', I am assuming it expressed as percent")
  }
  if (exists("Tfr", x, mode = "numeric", inherits=FALSE)) {
    range_check(x, strict.range = strict.range)
    return(x)
  } else if (exists("Tpc", x, mode = "numeric", inherits=FALSE)) {
    x[["Tfr"]] <- x[["Tpc"]] / 100
    x[["Tpc"]] <- NULL
    range_check(x, strict.range = strict.range)
    return(x)
  } else {
    warning("No transmittance data found in object_spct")
    x[["Tfr"]] <- NA
    return(x)
  }
  range_check(x, strict.range = strict.range)
  return(x)
}

#' @describeIn check Specialization for response_spct.
#' @export
check.response_spct <- function(x, byref=TRUE, strict.range=TRUE, multiple.wl = 1L, ...) {

  x <- checkTimeUnit(x)

  if (exists("s.e.response", x, mode = "numeric", inherits=FALSE)) {
    return(x)
  } else if (exists("s.q.response", x, mode = "numeric", inherits=FALSE)) {
    return(x)
  } else if (exists("response", x, mode = "numeric", inherits=FALSE)) {
    x[["s.e.response"]] <- x[["response"]]
    x[["response"]] <- NULL
    warning("Found variable 'response', I am assuming it is expressed on an energy basis")
    return(x)
  } else if (exists("signal", x, mode = "numeric", inherits=FALSE)) {
    x[["s.e.response"]] <- x[["signal"]]
    x[["signal"]] <- NULL
    warning("Found variable 'signal', I am assuming it is expressed on an energy basis")
    return(x)
  } else {
    warning("No response data found in response_spct")
    x[["s.e.response"]] <- NA
    return(x)
  }
}

#' @describeIn check Specialization for source_spct.
#' @export
check.source_spct <- function(x, byref=TRUE, strict.range=FALSE, multiple.wl = 1L, ...) {

  range_check <- function(x, strict.range) {
    if (is.null(strict.range)) {
      return()
    }
    if (exists("s.e.irrad", x, inherits = FALSE)) {
      s.e.min <- min(x$s.e.irrad, na.rm = TRUE)
      if (s.e.min < 0) {
        message.text <- paste("Negative spectral energy irradiance values; minimun s.e.irrad =", signif(s.e.min, 2))
        if (strict.range) {
          stop(message.text)
        } else {
          warning(message.text)
        }
      }
    }
    if (exists("s.q.irrad", x, inherits = FALSE)) {
      s.q.min <- min(x$s.q.irrad, na.rm = TRUE)
      if (s.q.min < 0) {
        message.text <- paste("Negative spectral photon irradiance values; minimun s.q.irrad =", signif(s.q.min, 2))
        if (strict.range) {
          stop(message.text)
        } else {
          warning(message.text)
        }
      }
    }
  }

  x <- checkTimeUnit(x)

  if (is.null(is_effective(x))) {
    setBSWFUsed(x, "none")
    warning("Missing atrribute 'bswf.used' set to 'none'")
  }
  if (exists("s.e.irrad", x, mode = "numeric", inherits=FALSE)) {
    NULL
  } else if (exists("s.q.irrad", x, mode = "numeric", inherits=FALSE)) {
    NULL
  } else if (exists("irradiance", x, mode = "numeric", inherits=FALSE)) {
    x[["s.e.irradiance"]] <- x[["irradiance"]]
    x[["irradiance"]] <- NULL
    warning("Found variable 'irradiance', I am assuming it is expressed on an energy basis")
  } else {
    warning("No spectral irradiance data found in source_spct")
    x[["s.e.irrad"]] <- NA
    return(x)
  }
  range_check(x, strict.range = strict.range)
  return(x)
}

#' @describeIn check Specialization for chroma_spct.
#' @export
check.chroma_spct <- function(x, byref=TRUE, strict.range=TRUE, multiple.wl = 1L, ...) {
  names_x <- names(x)
  idxs <- grep("[XYZ]", names_x)
  names2lc <- names_x[idxs]
  setnames(x, names2lc, tolower(names2lc))
  if (!exists("x", x, mode="numeric", inherits=FALSE)) {
    warning("Chromaticity coordinate 'x' data missing")
    x[["x"]] <- NA
  }
  if (!exists("y", x, mode="numeric", inherits=FALSE)) {
    warning("Chromaticity coordinate 'y' data missing")
    x[["y"]] <- NA
  }
  if (!exists("z", x, mode="numeric", inherits=FALSE)) {
    warning("Chromaticity coordinate 'z' data missing")
    x[["z"]] <- NA
  }
  return(x)
}


# set class ---------------------------------------------------------------

#' Remove "generic_spct" and derived class attributes.
#'
#' Removes from an spectrum object the class attibutes "generic_spct" and any
#' derived class attribute such as "source_spct". \strong{This operation is done
#' by reference!}
#'
#' @param x an R object.
#' @export
#'
#' @note If \code{x} is an object of any of the spectral classes defined
#'   in this package, this function changes by reference the spectrum
#'   object into the underlying data.table object. Otherwise, it just leaves \code{x}
#'   unchanged.
#'
#' @return A character vector containing the removed class attribute values.
#'   This is different to the behaviour of function \code{unlist} in base R!
#'
#' @family set and unset spectral class functions
#'
rmDerivedSpct <- function(x) {
  name <- substitute(x)
  spctclasses <- spct_classes()
  allclasses <- class(x)
  setattr(x, "class", setdiff(allclasses, spctclasses))
  setattr(x, "spct.version", NULL)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(setdiff(allclasses, class(x)))
}

#' Convert an R object into a spectrum object.
#'
#' Sets the class attibute of a data.frame or data.table or a different spectral
#' object to "generic_spct". If the object is a data.frame is is made a
#' data.table in the process.
#'
#' @param x data.frame, data.table, list or generic_spct and derived classes
#' @param multiple.wl numeric Maximum number of repeated w.length entries with same value.
#'
#' @export
#' @exportClass generic_spct
#' @family set and unset spectral class functions
#'
setGenericSpct <- function(x, multiple.wl = 1L) {
  name <- substitute(x)
  rmDerivedSpct(x)
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  if (!is.generic_spct(x)){
    setattr(x, "class", c("generic_spct", class(x)))
    setattr(x, "spct.tags", NA)
  }
  x <- check(x, multiple.wl = multiple.wl)
  setattr(x, "spct.version", 2)
#  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' @describeIn setGenericSpct Set class of a an object to "cps_spct".
#'
#' @export
#' @exportClass cps_spct
#'
setCpsSpct <- function(x, strict.range = TRUE, multiple.wl = 1L) {
  name <- substitute(x)
  rmDerivedSpct(x)
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.cps_spct(x)) {
    setattr(x, "class", c("cps_spct", class(x)))
  }
  x <- check(x, strict.range = strict.range)
#  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' @describeIn setGenericSpct Set class of an object to "filter_spct".
#'
#' @param Tfr.type character A string, either "total" or "internal".
#' @param strict.range logical Flag indicating whether off-range values result in an
#'   error instead of a warning.
#' @export
#' @exportClass filter_spct
#'
setFilterSpct <- function(x, Tfr.type=c("total", "internal"),
                          strict.range = TRUE, multiple.wl = 1L) {
  name <- substitute(x)
  if ((is.object_spct(x) || is.filter_spct(x)) && getTfrType(x) != "unknown") {
    if (length(Tfr.type) > 1) {
      Tfr.type <- getTfrType(x)
    } else {
      warning("Replacing existing attribute 'Tfr.type' ", getTfrType(x))
    }
  }
  rmDerivedSpct(x)
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.filter_spct(x)) {
    setattr(x, "class", c("filter_spct", class(x)))
  }
  setTfrType(x, Tfr.type[1])
  x <- check(x, strict.range = strict.range)
#  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' @describeIn setGenericSpct Set class of a an object to "reflector_spct".
#'
#' @param Rfr.type character A string, either "total" or "specular".
#' @export
#' @exportClass reflector_spct
#'
setReflectorSpct <- function(x, Rfr.type=c("total", "specular"),
                             strict.range = TRUE, multiple.wl = 1L) {
  name <- substitute(x)
  if ((is.object_spct(x) || is.reflector_spct(c)) && getRfrType(x) != "unknown") {
    if (length(Rfr.type) > 1) {
      Rfr.type <- getRfrType(x)
    } else {
      warning("Replacing existing attribute 'Rfr.type' ", getRfrType(x))
    }
  }
  rmDerivedSpct(x)
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.reflector_spct(x)) {
    setattr(x, "class", c("reflector_spct", class(x)))
  }
  setRfrType(x, Rfr.type[1])
  x <- check(x, strict.range = strict.range)
#  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' @describeIn setGenericSpct Set class of an object to "object_spct".
#'
#' @export
#' @exportClass object_spct
#'
setObjectSpct <- function(x,
                          Tfr.type=c("total", "internal"),
                          Rfr.type=c("total", "specular"),
                          strict.range = TRUE, multiple.wl = 1L) {
  name <- substitute(x)
  if ((is.filter_spct(x) || is.object_spct(x)) && getTfrType(x) != "unknown") {
    if (length(Tfr.type) > 1) {
      Tfr.type <- getTfrType(x)
    } else {
      warning("Replacing existing attribute 'Tfr.type' ", getTfrType(x))
    }
  }
  if ((is.reflector_spct(x) || is.object_spct(x)) && getRfrType(x) != "unknown") {
    if (length(Rfr.type) > 1) {
      Rfr.type <- getRfrType(x)
    } else {
      warning("Replacing existing attribute 'Rfr.type' ", getRfrType(x))
    }
  }
  rmDerivedSpct(x)
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.object_spct(x)) {
    setattr(x, "class", c("object_spct", class(x)))
  }
  setTfrType(x, Tfr.type)
  setRfrType(x, Rfr.type)
  x <- check(x, strict.range = strict.range)
#  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' @describeIn setGenericSpct Set class of an object to "response_spct".
#'
#' @param time.unit character A string "second", "day" or "exposure".
#' @export
#' @exportClass response_spct
#'
setResponseSpct <- function(x, time.unit="second", multiple.wl = 1L) {
  name <- substitute(x)
  rmDerivedSpct(x)
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.response_spct(x)) {
    setattr(x, "class", c("response_spct", class(x)))
  }
  setTimeUnit(x, time.unit)
  x <- check(x)
#  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' @describeIn setGenericSpct Set class of an object to "source_spct".
#'
#' @param bswf.used character A string, either "none" or the name of a BSWF.
#' @export
#' @exportClass source_spct
#'
setSourceSpct <- function(x, time.unit="second", bswf.used=c("none", "unknown"),
                          strict.range = FALSE, multiple.wl = 1L) {
  name <- substitute(x)
  rmDerivedSpct(x)
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.source_spct(x)) {
    setattr(x, "class", c("source_spct", class(x)))
  }
  setTimeUnit(x, time.unit)
  setBSWFUsed(x, bswf.used = bswf.used)
  x <- check(x, strict.range = strict.range)
#  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' @describeIn setGenericSpct Set class of an object to "chroma_spct".
#'
#' @export
#' @exportClass chroma_spct
#'
setChromaSpct <- function(x, multiple.wl = 1L) {
  name <- substitute(x)
  rmDerivedSpct(x)
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.chroma_spct(x)) {
    setattr(x, "class", c("chroma_spct", class(x)))
  }
  x <- check(x)
#  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}


# is functions for spct classes --------------------------------------------

#' Query class of spectrum objects
#'
#' Functions to check if an object is of a given type of spectrum, or coerce it if
#' possible.
#'
#' @param x an R object.
#'
#' @return These functions return \code{TRUE} if its argument is a of the queried type
#'   of spectrum and \code{FALSE} otherwise.
#'
#' @note Derived types also return TRUE for a query for a base type such as
#' \code{generic_spct}.
#'
#' @examples
#' is.source_spct(sun.spct)
#' is.filter_spct(sun.spct)
#' is.generic_spct(sun.spct)
#' is.any_spct(sun.spct)
#'
#' @export is.generic_spct
#' @rdname is.generic_spct
#'
is.generic_spct <- function(x) inherits(x, "generic_spct")

#' @rdname is.generic_spct
#' @export
#'
is.cps_spct <- function(x) inherits(x, "cps_spct")

#' @rdname is.generic_spct
#' @export
#'
is.source_spct <- function(x) inherits(x, "source_spct")

#' @rdname is.generic_spct
#' @export
#'
is.response_spct <- function(x) inherits(x, "response_spct")

#' @rdname is.generic_spct
#' @export
#'
is.filter_spct <- function(x) inherits(x, "filter_spct")

#' @rdname is.generic_spct
#' @export
#'
is.reflector_spct <- function(x) inherits(x, "reflector_spct")

#' @rdname is.generic_spct
#' @export
#'
is.object_spct <- function(x) inherits(x, "object_spct")

#' @rdname is.generic_spct
#' @export
#'
is.chroma_spct <- function(x) inherits(x, "chroma_spct")

#' @rdname is.generic_spct
#'
#' @export
#'
is.any_spct <- function(x) {
  inherits(x, spct_classes())
}

#' Query which is the class of an spectrum
#'
#' Functions to check if an object is a generic spectrum, or coerce it if
#' possible.
#'
#' @param x any R object
#'
#' @return class_spct returns a vector containing all matching xxxx.spct
#'   classes.
#'
#' @export
#'
class_spct <- function(x) {
#  intersect(spct_classes(), class(x)) # alters order!
  class(x)[class(x) %in% spct_classes()] # maintains order
}

#' Query if it is an spectrum is tagged
#'
#' Functions to check if an spct object contains tags.
#'
#' @param x any R object
#'
#' @return is_tagged returns TRUE if its argument is a an spectrum
#' that contains tags and FALSE if it is an untagged spectrun, but
#' returns NA for any other R object.
#'
#' @export
#'
#' @family tagging and related functions
#'
is_tagged <- function(x) {
  if (!is.any_spct(x)) {
    return(NA)
  } else {
    tags <- attr(x, "spct.tags", exact=TRUE)
    return(!is.null(tags) && length(tags) > 0 && !is.na(tags[[1]]))
  }
}

# is_photon_based ---------------------------------------------------------

#' Query if a spectrum contains photon- or energy-based data.
#'
#' Functions to check if \code{source_spct} and \code{response_spct} objects
#' contains photon-based or energy-based data.
#'
#' @param x any R object
#'
#' @return \code{is_photon_based} returns \code{TRUE} if its argument is a a
#'   \code{source_spct} or a \code{response_spct} object that contains photon
#'   base data and \code{FALSE} if such an object does not contain such data,
#'   but returns \code{NA} for any other R object, including those belonging
#'   other \code{generic_spct}-derived classes.
#'
#' @export
#' @family query units functions
#'
#' @rdname is_photon_based
#'
is_photon_based <- function(x) {
  if (is.source_spct(x)) {
    return("s.q.irrad" %in% names(x))
  } else if (is.response_spct(x)) {
    return("s.q.response" %in% names(x))
  } else {
    return(NA)
  }
}

# is_energy_based ---------------------------------------------------------

#' @rdname is_photon_based
#'
#' @return \code{is_energy_based} returns \code{TRUE} if its argument is a a \code{source_spct} or
#' a \code{response_spct} object that contains energy base data and \code{FALSE} if such an
#' object does not contain such data, but returns \code{NA} for any other R object,
#' including those belonging other \code{generic_spct}-derived classes
#'
#' @export
#'
is_energy_based <- function(x) {
  if (is.source_spct(x)) {
    return("s.e.irrad" %in% names(x))
  } else if (is.response_spct(x)) {
    return("s.e.response" %in% names(x))
  } else {
    return(NA)
  }
}

# is_absorbance_based ---------------------------------------------------------

#' Query if a spectrum contains absorbance or transmittance data
#'
#' Functions to check if an filter spectrum contains spectral absorbance data or
#' spectral transmittance data.
#'
#' @param x an R object
#'
#' @return \code{is_absorbance_based} returns TRUE if its argument is a \code{filter_spct}
#' object that contains spectral absorbance data and FALSE if it does not contain
#' such data, but returns NA for any other R object, including those belonging
#' other \code{generic_spct}-derived classes.
#'
#' @export
#' @family query units functions
#'
#' @rdname is_absorbance_based
#'
is_absorbance_based <- function(x) {
  if (is.filter_spct(x)) {
    return("A" %in% names(x))
  } else {
    return(NA)
  }
}

# is_transmittance_based ---------------------------------------------------------

#' @rdname is_absorbance_based
#'
#' @return \code{is_transmittance_based} returns TRUE if its argument is a a \code{filter_spct}
#' object that contains spectral transmittance data and FALSE if it does not contain
#' such data, but returns NA for any other R object, including those belonging
#' other \code{generic_spct}-derived classes.
#'
#' @export
#'
is_transmittance_based <- function(x) {
  if (is.filter_spct(x)) {
    return("Tfr" %in% names(x))
  } else {
    return(NA)
  }
}

# as functions for spct classes --------------------------------------------

#' Return a copy of an R object as an spectrum object
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object
#'
#' @return These functions return a copy of \code{x} converted into a given
#'   class of spectral object, if \code{x} is a valid argument to the
#'   correcponding set function.
#'
#' @export
#'
#' @family creation of spectral objects functions
#' @rdname as.generic_spct
#'
as.generic_spct <- function(x) {
  y <- copy(x)
  setGenericSpct(y)
}

#' @rdname as.generic_spct
#'
#' @export
#'
as.cps_spct <- function(x) {
  y <- copy(x)
  setCpsSpct(y)
}

#' @rdname as.generic_spct
#'
#' @param time.unit character A string, "second", "day" or "exposure"
#' @param bswf.used character
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning
#'
#' @export
#'
as.source_spct <- function(x,
                           time.unit=c("second", "day", "exposure"),
                           bswf.used=c("none", "unknown"),
                           strict.range = FALSE) {
  y <- copy(x)
  setSourceSpct(y, time.unit, strict.range = strict.range, bswf.used = bswf.used)
}

#' @rdname as.generic_spct
#'
#' @export
#'
as.response_spct <- function(x, time.unit = "second") {
  y <- copy(x)
  setResponseSpct(y, time.unit = time.unit)
}

#' @rdname as.generic_spct
#'
#' @param Tfr.type a character string, either "total" or "internal"
#'
#' @export
#'
as.filter_spct <- function(x, Tfr.type=c("total", "internal"), strict.range = TRUE) {
  y <- copy(x)
  setFilterSpct(y, Tfr.type, strict.range = strict.range)
}

#' @rdname as.generic_spct
#'
#' @param Rfr.type a character string, either "total" or "specular"
#'
#' @export
#'
as.reflector_spct <- function(x, Rfr.type = c("total", "specular"), strict.range = TRUE) {
  y <- copy(x)
  setReflectorSpct(y, Rfr.type = Rfr.type, strict.range = strict.range)
}

#' @rdname as.generic_spct
#'
#' @export
#'
as.object_spct <- function(x,
                           Tfr.type=c("total", "internal"),
                           Rfr.type=c("total", "specular"),
                           strict.range = TRUE) {
  y <- copy(x)
  setObjectSpct(y, Tfr.type = Tfr.type, Rfr.type = Rfr.type,
                strict.range = strict.range)
}

#' @rdname as.generic_spct
#'
#' @export
#'
as.chroma_spct <- function(x) {
  y <- copy(x)
  setChromaSpct(y)
}


# time.unit attribute -----------------------------------------------------

#' Set the "time.unit" attribute of an existing source_spct object
#'
#' Funtion to set by reference the "time.unit" attribute
#'
#' @param x a source_spct object
#' @param time.unit a character string, either "second", "hour", "day",
#'   "exposure" or "none", or a lubridate::duration
#' @param override.ok logical Flag that can be used to silence warning when
#'   overwritting an existing attribute value (used internally)
#'
#' @return x
#'
#' @note if x is not a source_spct or response_spct object, x is not modified.
#'   The behaviour of this function is 'unusual' in that the default for
#'   parameter \code{time.unit} is used only if \code{x} does not already have
#'   this attribute set. \code{time.unit = "hour"} is currently not fully
#'   supported.
#'
#' @export
#' @family time attribute functions
#'
setTimeUnit <- function(x,
                        time.unit = c("second", "hour", "day", "exposure", "none"),
                        override.ok = FALSE) {
  old.time.unit <- getTimeUnit(x)
  override.ok <- override.ok ||
    is.character(old.time.unit) &&
    old.time.unit %in% c("unknown", "none", time.unit)
  if (!override.ok) {
    warning("Overrriding existing 'time.unit' '", old.time.unit,
            "' with '", time.unit, "' may invalidate data!")
  }

  if (length(time.unit) > 1) {
    if (getTimeUnit(x) != "unknown") {
      time.unit <- getTimeUnit(x)
    } else {
      time.unit <- time.unit[[1]]
    }
  }
  if (is.source_spct(x) || is.response_spct(x)) {
    if (is.character(time.unit)) {
      if (!(time.unit %in% c("second", "hour", "day", "none", "exposure", "unknown"))) {
        warning("Unrecognized 'time.unit' argument ", time.unit, " set to 'unknown'.")
        time.unit <- "unknown"
      }
      else if (lubridate::is.duration(time.unit)) {
        if (time.unit <= duration(0, "seconds")) {
          stop("When 'time.unit' is a duration, it must be > 0")
        }
      }
    }
    setattr(x, "time.unit", time.unit)
  }
  invisible(x)
}

#' Get the "time.unit" attribute of an existing source_spct object
#'
#' Funtion to read the "time.unit" attribute
#'
#' @param x a source_spct object
#' @param force.duration logical If TRUE a lubridate::duration is returned even
#'   if the object attribute is a character string, if no conversion is possible
#'   NA is returned.
#'
#' @return character string or a lubridate::duration
#'
#' @note if x is not a \code{source_spct} or a \code{response_spct} object, NA
#' is returned
#'
#' @export
#' @family time attribute functions
#'
getTimeUnit <- function(x, force.duration = FALSE) {
  if (is.source_spct(x) || is.response_spct(x)) {
    time.unit <- attr(x, "time.unit", exact = TRUE)
    if (is.null(time.unit)) {
      # need to handle objects created with old versions
      time.unit <- "unknown"
    }
    # needed because of bad handling of defaults in constructor
    if (force.duration && is.character(time.unit)) {
      time.unit <- time.unit[[1]]
      time.unit <- char2duration(time.unit)
    }
    return(time.unit)
  } else {
    return(NA)
  }
}

#' Convert the "time.unit" attribute of an existing source_spct object
#'
#' Funtion to set the "time.unit" attribute and simultaneously rescaling the
#' spectral data to be expressed in the new time unit. The change is done
#' by reference ('in place')
#'
#' @param x a source_spct object
#' @param time.unit a character string, either "second", "hour", "day",
#'   "exposure" or "none", or a lubridate::duration
#' @param byref logical indicating if new object will be created by reference or
#'   by copy of \code{x}
#'
#' @return x possibly with the \code{time.unit} attribute modified
#'
#' @note if x is not a \code{source_spct} or a \code{response_spct} object, or
#'   time.unit is NULL x is returned unchanged, if the existing or new time.unit
#'   cannot be converted to a duration, then the returned spectrum will contain
#'   NAs.
#'
#' @export
#' @family time attribute functions
#'
convertTimeUnit <- function(x, time.unit = NULL, byref = TRUE) {
#  if (byref) x.out <- x else x.out <- copy(x) # needs fixing!
  if (is.null(time.unit) || !is.any_spct(x)) {
    return(invisible(x))
  }
  x.out <- checkTimeUnit(x)

  new.time.unit <- char2duration(time.unit)
  old.time.unit <- getTimeUnit(x.out, force.duration = TRUE)

  factor <- as.numeric(new.time.unit) / as.numeric(old.time.unit)

  columns <- intersect(names(x.out), c("s.e.irrad", "s.q.irrad", "s.e.response", "s.q.response") )

  for (col in columns) {
    x.out[[col]] <- x.out[[col]] * factor
  }

  if (length(columns) > 0) {
    setTimeUnit(x.out, time.unit, override.ok = TRUE)
  } else {
    warning("Nothing to convert to new time unit")
  }

  invisible(x.out)
}


#' Check the "time.unit" attribute of an existing source_spct object
#'
#' Funtion to read the "time.unit" attribute
#'
#' @param x a source_spct object
#'
#' @return x possibly with the \code{time.unit} attribute modified
#'
#' @note if x is not a \code{source_spct} or a \code{response_spct} object, NA
#' is returned
#'
#' @export
#' @family time attribute functions
#'
checkTimeUnit <- function(x) {
  if (is.source_spct(x) || is.response_spct(x)) {
    time.unit <- getTimeUnit(x)
    if (is.null(time.unit)) {
      setTimeUnit(x, "second")
      warning("Missing attribute 'time.unit' set to 'second'")
    }

    if (is.character(time.unit) &&
        !(time.unit %in% c("second", "minute", "hour", "day", "exposure", "none", "unknown"))) {
      stop("'time.unit' ",  time.unit, " is unknown")
    } else if (lubridate::is.duration(time.unit)) {
      if (time.unit <= lubridate::duration(0, "seconds")) {
        stop("When 'time.unit' is a duration, it must be > 0")
      }
    } else if (!is.character(time.unit)) {
      stop("'time.unit' must be of class character or lubridate::duration")
    }
  }
  invisible(x)
}

# private
char2duration <- function(time.unit) {
  if (is.character(time.unit)) {
    time.duration <- switch(time.unit,
                        second  = lubridate::duration(1, "seconds"),
                        minute  = lubridate::duration(1, "minutes"),
                        hour    = lubridate::duration(1, "hours"),
                        day     = lubridate::duration(1, "days"),
                        exposure = lubridate::duration(NA),
                        none    = lubridate::duration(NA),
                        unknown = lubridate::duration(NA)
    )
  } else if (lubridate::is.duration(time.unit)) {
    time.duration <- time.unit
  }
  return(time.duration)
}


# bswf attribute -----------------------------------------------------

#' Set the "bswf.used" attribute
#'
#' Funtion to set by reference the "time.unit" attribute of an existing
#' source_spct object
#'
#' @param x a source_spct object
#' @param bswf.used a character string, either "none" or the name of a BSWF
#'
#' @return x
#'
#' @note if x is not a source_spct, x is not modified. The behaviour of this
#'   function is 'unusual' in that the default for parameter \code{bswf.used} is
#'   used only if \code{x} does not already have this attribute set.
#'   \code{time.unit = "hour"} is currently not fully supported.
#'
#' @export
#' @family BSWF attribute functions
#'
setBSWFUsed <- function(x, bswf.used=c("none", "unknown")) {
  if (is.null(bswf.used) || length(bswf.used) < 1) {
    bswf.used <- "none"
  }
  if (length(bswf.used) > 1) {
    if (is_effective(x)) {
      bswf.used <- getBSWFUsed(x)
    } else {
      bswf.used <- bswf.used[[1]]
    }
  }
  if (is.source_spct(x)) {
    if  (!(is.character(bswf.used))) {
      warning("Only character strings are valid vlues for 'bswf.used' argument")
      bswf.used <- "unknown"
    }
    setattr(x, "bswf.used", bswf.used)
  }
  return(x)
}

#' Get the "bswf.used" attribute
#'
#' Funtion to read the "time.unit" attribute of an existing source_spct object
#'
#' @param x a source_spct object
#'
#' @return character string
#'
#' @note if x is not a \code{source_spct} object, NA is retruned
#'
#' @export
#' @family BSWF attribute functions
#'
getBSWFUsed <- function(x) {
  if (is.source_spct(x)) {
    bswf.used <- attr(x, "bswf.used", exact = TRUE)
    if (is.null(bswf.used) || length(bswf.used) < 1) {
      # need to handle objects created with old versions
      bswf.used <- "none"
    }
    return(bswf.used[[1]])
  } else {
    return(NA)
  }
}

# is_effective.source_spct defined in file "waveband.class.r" to avoid the need
# of using colate to get the documentation in the correct order.

# Tfr.type attribute ------------------------------------------------------

#' Set the "Tfr.type" attribute
#'
#' Funtion to set by reference the "Tfr.type" attribute of an existing
#' filter_spct or object_spct object
#'
#' @param x a filter_spct or an object_spct object
#' @param Tfr.type a character string, either "total" or "internal"
#'
#' @return x
#'
#' @note if x is not a filter_spct or an object_spct object, x is not modified
#'   The behaviour of this function is 'unusual' in that the default for
#'   parameter \code{Tfr.type} is used only if \code{x} does not already have
#'   this attribute set.
#'
#' @export
#' @family Tfr attribute functions
#'
setTfrType <- function(x, Tfr.type=c("total", "internal")) {
  if (length(Tfr.type) > 1) {
    if (getTfrType(x) != "unknown") {
      Tfr.type <- getTfrType(x)
    } else {
      Tfr.type <- Tfr.type[[1]]
    }
  }
  if (is.filter_spct(x) || is.object_spct(x)) {
    if  (!(Tfr.type %in% c("total", "internal", "unknown"))) {
      warning("Invalid 'Tfr.type' argument, only 'total' and 'internal' supported.")
      return(x)
    }
    setattr(x, "Tfr.type", Tfr.type)
  }
  return(x)
}

#' Get the "Tfr.type" attribute
#'
#' Funtion to read the "Tfr.type" attribute of an existing filter_spct or
#' object_spct object.
#'
#' @param x a filter_spct or object_spct object
#'
#' @return character string
#'
#' @note If x is not a \code{filter_spct} or an \code{object_spct} object,
#'   \code{NA} is returned.
#'
#' @export
#' @family Tfr attribute functions
#'
getTfrType <- function(x) {
  if (is.filter_spct(x) || is.object_spct(x)) {
    Tfr.type <- attr(x, "Tfr.type", exact = TRUE)
    if (is.null(Tfr.type) || is.na(Tfr.type)) {
      # need to handle objects created with old versions
      Tfr.type <- "unknown"
    }
    return(Tfr.type[[1]])
  } else {
    return(NA)
  }
}

# Rfr.type attribute ------------------------------------------------------

#' Set the "Rfr.type" attribute
#'
#' Funtion to set by reference the "Rfr.type" attribute  of an existing
#' reflector_spct or object_spct object.
#'
#' @param x a reflector_spct or an object_spct object
#' @param Rfr.type a character string, either "total" or "specular"
#'
#' @return x
#'
#' @note if x is not a reflector_spct or object_spct object, x is not modified.
#'   The behaviour of this function is 'unusual' in that the default for
#'   parameter Rfr.type is used only if \code{x} does not already have this
#'   attribute set.
#'
#' @export
#' @family Rfr attribute functions
#'
setRfrType <- function(x, Rfr.type=c("total", "specular")) {
  if (length(Rfr.type) > 1) {
    if (getRfrType(x) != "unknown") {
      Rfr.type <- getRfrType(x)
    } else {
      Rfr.type <- Rfr.type[[1]]
    }
  }
  if (is.reflector_spct(x) || is.object_spct(x)) {
    if  (!(Rfr.type %in% c("total", "specular", "unknown"))) {
      warning("Invalid 'Rfr.type' argument, only 'total' and 'internal' supported.")
      return(x)
    }
    setattr(x, "Rfr.type", Rfr.type)
  }
  return(x)
}

#' Get the "Rfr.type" attribute
#'
#' Funtion to read the "Rfr.type" attribute of an existing reflector_spct
#' object.
#'
#' @param x a source_spct object
#'
#' @return character string
#'
#' @note if x is not a \code{filter_spct} object, \code{NA} is returned
#'
#' @export
#' @family Rfr attribute functions
#'
getRfrType <- function(x) {
  if (is.reflector_spct(x) || is.object_spct(x)) {
    Rfr.type <- attr(x, "Rfr.type", exact = TRUE)
    if (is.null(Rfr.type) || is.na(Rfr.type)) {
      # need to handle objects created with old versions
      Rfr.type <- "unknown"
    }
    return(Rfr.type[[1]])
  } else {
    return(NA)
  }
}

#' Get the "spct.version" attribute
#'
#' Funtion to read the "spct.version" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#'
#' @return numeric value
#'
#' @note if x is not a \code{generic_spct} object, \code{NA} is returned,
#'   and if it the attribute is missing, zero is returned with a warning.
#'
#' @export
#'
getSpctVersion <- function(x) {
  if (is.any_spct(x)) {
    version <- attr(x, "spct.version", exact = TRUE)
    if (is.null(version)) {
      # need to handle objects created with old versions
      version <- 0L
    }
  } else {
    version <- NA
  }
  version
}

#' Check that the "spct.version" attribute is set
#'
#' Funtion to check the "spct.version" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#'
#' @return numeric value
#'
#' @note if x is not a \code{generic_spct} object, \code{NA} is returned,
#'   and if it the attribute is missing, zero is returned with a warning.
#'
#' @keywords internal
#'
checkSpctVersion <- function(x) {
  version <- getObjectVersion(x)
  stopifnot(!is.na(version))
  if (version < 1L) {
    warning("The object '", as.character(substitute(x)),
            "' was created in a version (< 0.7.0) or has become corrupted")
  }
}
