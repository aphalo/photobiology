
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
check.default <- function(x, byref=FALSE, strict.range = FALSE, ...) {
  return(x)
}

#' @describeIn check Specialization for generic_spct.
#'
#' @param multiple.wl numeric Maximum number of repeated w.length entries with same value.
#'
#' @export
check.generic_spct <- function(x, byref=TRUE, strict.range = FALSE, multiple.wl = getMultipleWl(x), ...) {
  # fix old class attributes
  class.x <- class_spct(x)
  if (!("tbl_df") %in% class(x)) {
    x <- dplyr::as_data_frame(x)
  }
  class(x) <- union(class.x, class(x))
  # check variables
  if (exists("wl", x, mode = "numeric", inherits = FALSE)) {
    x <- dplyr::rename(x, w.length = wl)
  } else if (exists("wavelength", x, mode = "numeric", inherits = FALSE)) {
    x <- dplyr::rename(x, w.length = wavelength)
  } else if (exists("Wavelength", x, mode = "numeric", inherits = FALSE)) {
    x <- dplyr::rename(x, w.length = Wavelength)
  }

  if (exists("w.length", x, mode = "numeric", inherits = FALSE)) {
    if (nrow(x) && is.unsorted(x[["w.lengtgh"]], na.rm = TRUE, strictly = TRUE)) {
      stop("'w.length' must be sorted in ascending order and have unique values")
    }
  } else {
      stop("No wavelength data found in generic_spct")
  }

  if (nrow(x)) {
    wl.min <- min(x[["w.length"]], na.rm = TRUE)
    #  wl.max <- max(x$w.length, na.rm = TRUE)
    if (wl.min == Inf) {
      warning("No valid 'w.length' values found")
    } else if ((wl.min < 99.999 || wl.min > 5e3)) {
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
  }
  x
}

#' @describeIn check Specialization for cps_spct.
#' @export
check.cps_spct <- function(x, byref=TRUE, strict.range = FALSE, multiple.wl = getMultipleWl(x), ...) {

  range_check <- function(x) {
    cps.min <- min(x$cps, na.rm = TRUE)
    if (!is.null(strict.range) & (cps.min < -1e-8)) {
      message.text <- paste0("Off-range cps values:", signif(cps.min, 2))
      if (strict.range) {
        stop(message.text)
      } else {
        warning(message.text)
      }
    }
  }

  x <- check.generic_spct(x, multiple.wl = multiple.wl)

  if (exists("cps", x, mode = "numeric", inherits=FALSE)) {
    range_check(x)
    return(x)
  } else if (exists("counts.per.second", x, mode = "numeric", inherits=FALSE)) {
    x <- dplyr::rename(x, cps = counts.per.second)
    warning("Found variable 'counts.per.second', renamed it to 'cps'")
    range_check(x)
    return(x)
  } else {
    warning("No counts per second data found in cps_spct")
    x[["cps"]] = NA_real_
    return(x)
  }
}

#' @describeIn check Specialization for filter_spct.
#' @export
check.filter_spct <- function(x, byref=TRUE, strict.range = FALSE, multiple.wl = getMultipleWl(x), ...) {

  range_check <- function(x, strict.range) {
    if (!all(is.na(x[["Tfr"]]))) {
      Tfr.min <- min(x[["Tfr"]], na.rm = TRUE)
      Tfr.max <- max(x[["Tfr"]], na.rm = TRUE)
      if (!is.null(strict.range) & (Tfr.min < -1e-8 || Tfr.max > 1)) {
        message.text <- paste("Off-range transmittance values [", signif(Tfr.min, 2),
                              "...", signif(Tfr.max, 2), "] instead of  [0..1]", sep="")
        if (strict.range) {
          stop(message.text)
        } else {
          warning(message.text)
        }
      }
    }
  }

  x <- check.generic_spct(x, multiple.wl = multiple.wl)

  range_check_A <- function(x, strict.range) {
    if (!all(is.na(x[["A"]]))) {
      A.min <- min(x[["A"]], na.rm = TRUE)
      A.max <- max(x[["A"]], na.rm = TRUE)
      if (!is.null(strict.range) & (A.min < -1e-8 || A.max > 20)) {
        message.text <- paste("Off-range absorbance values [", signif(A.min, 2),
                              "...", signif(A.max, 2), "] instead of  [0..1]", sep="")
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
  if (is.null(getTfrType(x))) {
    setTfrType(x, "total")
    warning("Missing Tfr.type attribute replaced by 'total'")
  }
  # check and replace 'other' quantity names
  if (exists("transmittance", x, mode = "numeric", inherits = FALSE)) {
    x <- dplyr::rename(x, Tpc = transmittance)
    warning("Found variable 'transmittance', I am assuming it is expressed as percent")
  }
  if (exists("absorbance", x, mode = "numeric", inherits = FALSE)) {
    x <- dplyr::rename(x, A = absorbance)
    warning("Found variable 'absorbance', I am assuming it is in log10-based absorbance units")
  } else if (exists("Absorbance", x, mode = "numeric", inherits = FALSE)) {
    x <- dplyr::rename(x, A = Absorbance)
    warning("Found variable 'Absorbance', I am assuming it is in log10-based absorbance units")
  }
  # look for percentages and change them into fractions of one
  if (exists("Tfr", x, mode = "numeric", inherits=FALSE)) {
    range_check(x, strict.range=strict.range)
  } else if (exists("Tpc", x, mode = "numeric", inherits=FALSE)) {
    x[["Tfr"]] <- x[["Tpc"]] / 100
    x[["Tpc"]] <-  NULL
    range_check(x, strict.range=strict.range)
  } else if (exists("A", x, mode = "numeric", inherits=FALSE)) {
    range_check_A(x, strict.range=strict.range)
  } else {
    warning("No transmittance or absorbance data found in filter_spct")
    x[["Tfr"]] <- NA_real_
  }
  x
}

#' @describeIn check Specialization for reflector_spct.
#' @export
check.reflector_spct <- function(x, byref = TRUE, strict.range = FALSE, multiple.wl = getMultipleWl(x), ...) {

  range_check <- function(x, strict.range) {
    if (!all(is.na(x$Rfr))) {
      Rfr.min <- min(x$Rfr, na.rm = TRUE)
      Rfr.max <- max(x$Rfr, na.rm = TRUE)
      if (!is.null(strict.range) & (Rfr.min < -1e-8 ||  Rfr.max > 1)) {
        message.text <- paste0("Off-range reflectance values [", signif(Rfr.min, 2), "...",
                               signif(Rfr.max, 2), "] instead of  [0..1]", sep="")
        if (strict.range) {
          stop(message.text)
        } else {
          warning(message.text)
        }
      }
    }
  }

  x <- check.generic_spct(x, multiple.wl = multiple.wl)

  if (is.null(getRfrType(x))) {
    setRfrType(x, "total")
    warning("Missing Rfr.type attribute replaced by 'total'")
  }
  if (exists("reflectance", x, mode = "numeric", inherits=FALSE)) {
    x <- dplyr::rename(x, Rpc = reflectance)
    warning("Found variable 'reflectance', I am assuming it is expressed as percent")
  }
  if (exists("Rfr", x, mode = "numeric", inherits=FALSE)) {
    range_check(x, strict.range=strict.range)
  } else if (exists("Rpc", x, mode = "numeric", inherits=FALSE)) {
    x[["Rfr"]] <- x[["Rpc"]] / 100
    x[["Rpc"]] <- NULL
    range_check(x, strict.range=strict.range)
  } else {
    warning("No reflectance data found in reflector_spct")
    x[["Rfr"]] <- NA_real_
  }
  x
}

#' @describeIn check Specialization for object_spct.
#' @export
check.object_spct <- function(x, byref=TRUE, strict.range = FALSE, multiple.wl = getMultipleWl(x), ...) {

  range_check_Tfr <- function(x, strict.range) {
    if (!all(is.na(x[["Tfr"]]))) {
      Tfr.min <- min(x[["Tfr"]], na.rm = TRUE)
      Tfr.max <- max(x[["Tfr"]], na.rm = TRUE)
      if (!is.null(strict.range) & (Tfr.min < -1e-8 || Tfr.max > 1)) {
        message.text <- paste("Off-range transmittance values [", signif(Tfr.min, 2),
                              "...", signif(Tfr.max, 2), "] instead of  [0..1]", sep="")
        if (strict.range) {
          stop(message.text)
        } else {
          warning(message.text)
        }
      }
    }
  }

  x <- check.generic_spct(x, multiple.wl = multiple.wl)

  range_check_Rfr <- function(x, strict.range) {
    if (!all(is.na(x$Rfr))) {
      Rfr.min <- min(x[["Rfr"]], na.rm = TRUE)
      Rfr.max <- max(x[["Rfr"]], na.rm = TRUE)
      if (!is.na(Rfr.min) && !is.na(Rfr.max)) {
        if (!is.null(strict.range) & (Rfr.min < -1e-8 ||  Rfr.max > 1)) {
          message.text <- paste0("Off-range reflectance values [", signif(Rfr.min, 2), "...",
                                 signif(Rfr.max, 2), "] instead of  [0..1]", sep="")
          if (strict.range) {
            stop(message.text)
          } else {
            warning(message.text)
          }
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
    x <- dplyr::rename(x, Rpc = reflectance)
    warning("Found variable 'reflectance', I am assuming it is expressed as percent")
  }
  if (exists("Rfr", x, mode = "numeric", inherits=FALSE)) {
    range_check_Rfr(x, strict.range=strict.range)
  } else if (exists("Rpc", x, mode = "numeric", inherits=FALSE)) {
    x[["Rfr"]] <- x[["Rpc"]] / 100
    x[["Rpc"]] <- NULL
    range_check_Rfr(x, strict.range=strict.range)
  } else {
    warning("No reflectance data found in object_spct")
    x[["Rfr"]] <- NA
  }

  if (exists("transmittance", x, mode = "numeric", inherits=FALSE)) {
    x <- dplyr::rename(x, Tpc = transmittance)
    warning("Found variable 'transmittance', I am assuming it expressed as percent")
  }
  if (exists("Tfr", x, mode = "numeric", inherits=FALSE)) {
    range_check_Tfr(x, strict.range=strict.range)
  } else if (exists("Tpc", x, mode = "numeric", inherits=FALSE)) {
    x[["Tfr"]] <- x[["Tpc"]] / 100
    x[["Tpc"]] <- NULL
    range_check_Tfr(x, strict.range=strict.range)
  } else {
    warning("No transmittance or absorptance data found in object_spct")
    x[["Tfr"]] <- NA_real_
  }

  x
}

#' @describeIn check Specialization for response_spct.
#' @export
check.response_spct <- function(x, byref=TRUE, strict.range = FALSE, multiple.wl = getMultipleWl(x), ...) {

  x <- check.generic_spct(x, multiple.wl = multiple.wl)

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
    x[["s.e.response"]] <- NA_real_
    return(x)
  }
}

#' @describeIn check Specialization for source_spct.
#' @export
check.source_spct <- function(x, byref=TRUE, strict.range = FALSE, multiple.wl = getMultipleWl(x), ...) {

  range_check <- function(x, strict.range) {
    if (is.null(strict.range)) {
      return()
    }
    if (exists("s.e.irrad", x, inherits = FALSE)) {
      s.e.min <- min(x$s.e.irrad, na.rm = TRUE)
      if (s.e.min < -1e-8) {
        message.text <- paste("Negative spectral energy irradiance values; minimun s.e.irrad =",
                              signif(s.e.min, 2))
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

  x <- check.generic_spct(x, multiple.wl = multiple.wl)
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
    x[["s.e.irrad"]] <- NA_real_
    return(x)
  }
  range_check(x, strict.range = strict.range)
  return(x)
}

#' @describeIn check Specialization for chroma_spct.
#' @export
check.chroma_spct <- function(x, byref=TRUE, strict.range = FALSE, multiple.wl = getMultipleWl(x), ...) {
  names_x <- names(x)

  x <- check.generic_spct(x, multiple.wl = multiple.wl)

  idxs <- grep("[XYZ]", names_x)
  names(x)[idxs] <- tolower(names_x[idxs])
  if (!exists("x", x, mode="numeric", inherits=FALSE)) {
    warning("Chromaticity coordinate 'x' data missing")
    x[["x"]] <- NA_real_
  }
  if (!exists("y", x, mode="numeric", inherits=FALSE)) {
    warning("Chromaticity coordinate 'y' data missing")
    x[["y"]] <- NA_real_
  }
  if (!exists("z", x, mode="numeric", inherits=FALSE)) {
    warning("Chromaticity coordinate 'z' data missing")
    x[["z"]] <- NA_real_
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
#'   object into the underlying data.frame object. Otherwise, it just leaves \code{x}
#'   unchanged.
#'
#' @return A character vector containing the removed class attribute values.
#'   This is different to the behaviour of function \code{unlist} in base R!
#'
#' @family set and unset spectral class functions
#'
rmDerivedSpct <- function(x) {
  name <- substitute(x)
  allclasses <- class(x)
  class(x) <- setdiff(allclasses, spct_classes())
  attr(x, "spct.version") <- NULL
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(setdiff(allclasses, class(x)))
}

#' Convert an R object into a spectrum object.
#'
#' Sets the class attibute of a data.frame or an object of a derived
#' class to "generic_spct".
#'
#' @param x data.frame, list or generic_spct and derived classes
#' @param multiple.wl numeric Maximum number of repeated w.length entries with same value.
#'
#' @export
#' @exportClass generic_spct
#' @family set and unset spectral class functions
#'
setGenericSpct <- function(x, multiple.wl = 1L) {
  name <- substitute(x)
  rmDerivedSpct(x)
  if (!is.data.frame(x) || inherits(x, "data.table")) {
    x <- dplyr::as_data_frame(x)
  }
  if (!is.generic_spct(x)){
    class(x) <- c("generic_spct", class(x))
    attr(x, "spct.tags") <- NA
    setMultipleWl(x, multiple.wl = multiple.wl)
  }
  x <- check(x)
  attr(x, "spct.version") <- 2
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
setCpsSpct <- function(x, strict.range = FALSE, multiple.wl = 1L) {
  name <- substitute(x)
  rmDerivedSpct(x)
  if (!is.data.frame(x) || inherits(x, "data.table")) {
    x <- dplyr::as_data_frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.cps_spct(x)) {
    class(x) <- c("cps_spct", class(x))
  }
  x <- check(x, strict.range = strict.range)
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
                          strict.range = FALSE, multiple.wl = 1L) {
  name <- substitute(x)
  if ((is.object_spct(x) || is.filter_spct(x)) && getTfrType(x) != "unknown") {
    if (length(Tfr.type) > 1) {
      Tfr.type <- getTfrType(x)
    } else if (Tfr.type != getTfrType(x)) {
      warning("Changing attribute 'Tfr.type' from ", getTfrType(x),
              " into ", Tfr.type)
    }
  }
  rmDerivedSpct(x)
  if (!is.data.frame(x) || inherits(x, "data.table")) {
    x <- dplyr::as_data_frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.filter_spct(x)) {
    class(x) <- c("filter_spct", class(x))
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
                             strict.range = FALSE, multiple.wl = 1L) {
  name <- substitute(x)
  if ((is.object_spct(x) || is.reflector_spct(c)) && getRfrType(x) != "unknown") {
    if (length(Rfr.type) > 1) {
      Rfr.type <- getRfrType(x)
    } else if (Rfr.type != getRfrType(x)) {
      warning("Changing attribute 'Rfr.type' from ", getRfrType(x),
              " into ", Rfr.type)
    }
  }
  rmDerivedSpct(x)
  if (!is.data.frame(x) || inherits(x, "data.table")) {
    x <- dplyr::as_data_frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.reflector_spct(x)) {
    class(x) <- c("reflector_spct", class(x))
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
                          strict.range = FALSE, multiple.wl = 1L) {
  name <- substitute(x)
  if ((is.filter_spct(x) || is.object_spct(x)) && getTfrType(x) != "unknown") {
    if (length(Tfr.type) > 1) {
      Tfr.type <- getTfrType(x)
    } else if (Tfr.type != getTfrType(x)) {
      warning("Changing attribute 'Tfr.type' from ", getTfrType(x),
              " into ", Tfr.type)
    }
  }
  if ((is.reflector_spct(x) || is.object_spct(x)) && getRfrType(x) != "unknown") {
    if (length(Rfr.type) > 1) {
      Rfr.type <- getRfrType(x)
    } else if (Rfr.type != getRfrType(x)) {
      warning("Changing attribute 'Rfr.type' from ", getRfrType(x),
              " into ", Rfr.type)
    }
  }
  rmDerivedSpct(x)
  if (!is.data.frame(x) || inherits(x, "data.table")) {
    x <- dplyr::as_data_frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.object_spct(x)) {
    class(x) <- c("object_spct", class(x))
  }
  setTfrType(x, Tfr.type)
  setRfrType(x, Rfr.type)
  x <- check(x, strict.range = strict.range)
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
  if (!is.data.frame(x) || inherits(x, "data.table")) {
    x <- dplyr::as_data_frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.response_spct(x)) {
    class(x) <- c("response_spct", class(x))
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
  if (!is.data.frame(x) || inherits(x, "data.table")) {
    x <- dplyr::as_data_frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.source_spct(x)) {
    class(x) <- c("source_spct", class(x))
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
  if (!is.data.frame(x) || inherits(x, "data.table")) {
    x <- dplyr::as_data_frame(x)
  }
  if (!is.generic_spct(x)) {
    setGenericSpct(x, multiple.wl = multiple.wl)
  }
  if (!is.chroma_spct(x)) {
    class(x) <- c("chroma_spct", class(x))
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
  y <- x
  setGenericSpct(y)
}

#' @rdname as.generic_spct
#'
#' @export
#'
as.cps_spct <- function(x) {
  y <- x
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
  y <- x
  setSourceSpct(y, time.unit, strict.range = strict.range, bswf.used = bswf.used)
}

#' @rdname as.generic_spct
#'
#' @export
#'
as.response_spct <- function(x, time.unit = "second") {
  y <- x
  setResponseSpct(y, time.unit = time.unit)
}

#' @rdname as.generic_spct
#'
#' @param Tfr.type a character string, either "total" or "internal"
#'
#' @export
#'
as.filter_spct <- function(x, Tfr.type=c("total", "internal"), strict.range = FALSE) {
  y <- x
  setFilterSpct(y, Tfr.type, strict.range = strict.range)
}

#' @rdname as.generic_spct
#'
#' @param Rfr.type a character string, either "total" or "specular"
#'
#' @export
#'
as.reflector_spct <- function(x, Rfr.type = c("total", "specular"), strict.range = FALSE) {
  y <- x
  setReflectorSpct(y, Rfr.type = Rfr.type, strict.range = strict.range)
}

#' @rdname as.generic_spct
#'
#' @export
#'
as.object_spct <- function(x,
                           Tfr.type=c("total", "internal"),
                           Rfr.type=c("total", "specular"),
                           strict.range = FALSE) {
  y <- x
  setObjectSpct(y, Tfr.type = Tfr.type, Rfr.type = Rfr.type,
                strict.range = strict.range)
}

#' @rdname as.generic_spct
#'
#' @export
#'
as.chroma_spct <- function(x) {
  y <- x
  setChromaSpct(y)
}


# time.unit attribute -----------------------------------------------------

#' Set the "time.unit" attribute of an existing source_spct object
#'
#' Function to set by reference the "time.unit" attribute
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
  if (!(class(x)[1] %in% c("source_spct", "response_spct"))) {
    return(invisible(x))
  }
  name <- substitute(x)
  if (length(time.unit) > 1) {
    if (getTimeUnit(x) != "unknown") {
      time.unit <- getTimeUnit(x)
    } else {
      time.unit <- time.unit[[1]]
    }
  } else {
    old.time.unit <- getTimeUnit(x)
    override.ok <- override.ok ||
      is.character(old.time.unit) && old.time.unit %in% c("unknown", "none", time.unit)
    if (!override.ok && old.time.unit != time.unit[1]) {
      warning("Overrriding existing 'time.unit' '", old.time.unit,
              "' with '", time.unit, "' may invalidate data!")
    }
  }
  if (is.character(time.unit)) {
    if (!(time.unit %in% c("second", "hour", "day", "none", "exposure", "unknown"))) {
      warning("Unrecognized 'time.unit' argument ", time.unit, " set to 'unknown'.")
      time.unit <- "unknown"
    }
  } else if (lubridate::is.duration(time.unit)) {
    if (time.unit <= duration(0, "seconds")) {
      stop("When 'time.unit' is a duration, it must be > 0")
    }
  }
  attr(x, "time.unit") <- time.unit
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' Get the "time.unit" attribute of an existing source_spct object
#'
#' Function to read the "time.unit" attribute
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
#' Function to set the "time.unit" attribute and simultaneously rescaling the
#' spectral data to be expressed in the new time unit. The change is done
#' by reference ('in place')
#'
#' @param x a source_spct object
#' @param time.unit a character string, either "second", "hour", "day",
#'   "exposure" or "none", or a lubridate::duration
#' @param byref logical indicating if new object will be created by reference or
#'   by copy of \code{x} (currently ignored)
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
convertTimeUnit <- function(x, time.unit = NULL, byref = FALSE) {
  #  if (!byref) x.out <- x else x.out <- copy(x) # needs fixing!
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
#' Function to read the "time.unit" attribute
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
#' Function to set by reference the "time.unit" attribute of an existing
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
    name <- substitute(x)
    if  (!(is.character(bswf.used))) {
      warning("Only character strings are valid vlues for 'bswf.used' argument")
      bswf.used <- "unknown"
    }
    attr(x, "bswf.used") <- bswf.used
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' Get the "bswf.used" attribute
#'
#' Function to read the "time.unit" attribute of an existing source_spct object
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
#' Function to set by reference the "Tfr.type" attribute of an existing
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
  name <- substitute(x)
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
    attr(x, "Tfr.type") <- Tfr.type
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' Get the "Tfr.type" attribute
#'
#' Function to read the "Tfr.type" attribute of an existing filter_spct or
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
#' Function to set by reference the "Rfr.type" attribute  of an existing
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
  name <- substitute(x)
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
    attr(x, "Rfr.type") <- Rfr.type
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' Get the "Rfr.type" attribute
#'
#' Function to read the "Rfr.type" attribute of an existing reflector_spct
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
#' Function to read the "spct.version" attribute of an existing generic_spct
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
  if (is.any_spct(x) || is.old_spct(x)) {
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
#' Function to check the "spct.version" attribute of an existing generic_spct
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
  version <- getSpctVersion(x)
  stopifnot(!is.na(version))
  if (version < 1L) {
    warning("The object '", as.character(substitute(x)),
            "' was created in a version (< 0.7.0) or has become corrupted")
  }
}


# multiple wl -------------------------------------------------------------

#' Set the "multiple.wl" attribute
#'
#' Function to set by reference the "multiple.wl" attribute  of an existing
#' generic_spct or an object of a class derived from generic_spct.
#'
#' @param x a generic_spct object
#' @param multiple.wl numeric >= 1 If \code{multiple.wl} is \code{NULL}, the
#'   default, the attribute is not modified if it is already present and valid,
#'   and set to 1 otherwise.
#'
#' @return x
#'
#' @note if x is not a generic_spct or an object of a class derived from
#'   generic_spct, x is not modified. If \code{multiple.wl}
#'
#' @export
#' @family multiple.wl attribute functions
#'
setMultipleWl <- function(x, multiple.wl = NULL) {
  stopifnot(is.any_spct(x))
  name <- substitute(x)
  if (is.null(multiple.wl)) {
    multiple.wl <- getMultipleWl(x)
  } else {
    multiple.wl <- trunc(multiple.wl)
    stopifnot(multiple.wl > 0)
  }
  attr(x, "multiple.wl") <- multiple.wl
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' Get the "multiple.wl" attribute
#'
#' Function to read the "Tfr.type" attribute of an existing filter_spct or
#' object_spct object.
#'
#' @param x a generic_spct object
#'
#' @return integer
#'
#' @note If x is not a \code{generic_spct} or an object of a derived class
#'   \code{NA} is returned.
#'
#' @export
#' @family multiple.wl attribute functions
#'
getMultipleWl <- function(x) {
  if (is.any_spct(x)) {
    multiple.wl <- attr(x, "multiple.wl", exact = TRUE)
    if (is.null(multiple.wl) || is.na(multiple.wl) || !is.numeric(multiple.wl)) {
      # need to handle objects created with old versions
      multiple.wl <- 1
    }
    return(multiple.wl)
  } else {
    return(NA)
  }
}

