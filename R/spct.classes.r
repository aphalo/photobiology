
# names of all spectral classes -------------------------------------------

#' Function that returns a vector containing the names of spectra classes.
#'
#' @usage spct.classes()
#'
#' @export
#'
#' @return a character vector of class names
#'
spct.classes <- function() {
  c("generic.spct", "private.spct",
    "filter.spct", "reflector.spct",
    "source.spct", "object.spct",
    "response.spct", "chroma.spct")
}


# conditional setkey ------------------------------------------------------

#' Stolen from data.table except that test added so that if the same key is
#' already set setkeyv is not called.
#'
#' @usage setkey_spct(x, ..., verbose = getOption("datatable.verbose"), physical = TRUE)
#'
#' @param x spct object
#' @param ... columns
#' @param verbose logical
#' @param physical logical
#'
#' @note \code{\link[data.table]{setkey}}
#'
#' @keywords internal
#'
setkey_spct <- function (x, ..., verbose = getOption("datatable.verbose"), physical = TRUE)
{
  if (is.character(x))
    stop("x may no longer be the character name of the data.table. The possibility was undocumented and has been removed.")
  cols = as.character(substitute(list(...))[-1])
  if (!length(cols))
    cols = colnames(x)
  else if (identical(cols, "NULL"))
    cols = NULL
  if (is.any.spct(x) && !is.null(key(x)) && identical(cols, key(x)))
    invisible(x)
  setkeyv(x, cols, verbose = verbose, physical = physical)
}

# check -------------------------------------------------------------------

#' Generic check function.
#'
#' Check that an R object contains the expected data members.
#'
#' @usage check(x, byref, strict.range)
#'
#' @param x An R object
#' @param byref logical indicating if new object will be created by reference or by copy of \code{x}
#' @param strict.range logical indicating whether off-range values result in an error instead of a warning
#' @export
#'
#' @family data validity check functions
#'
check <- function(x, byref, strict.range) UseMethod("check")

#' @describeIn check Default for generic function.
#' @export
check.default <- function(x, byref=FALSE, strict.range=TRUE) {
  return(x)
}

#' check Specialization for private.spct.
#' @export
#' @keywords internal
check.private.spct <- function(x, byref=TRUE, strict.range=TRUE) {
  if (exists("w.length", x, mode = "numeric", inherits=FALSE) &&
        exists("numbers", x, mode = "numeric", inherits=FALSE)) {
    return(x)
  }
}

#' @describeIn check Specialization for generic.spct.
#' @export
check.generic.spct <- function(x, byref=TRUE, strict.range=TRUE) {
  if (exists("w.length", x, mode = "numeric", inherits=FALSE)) {
    NULL
  } else if (exists("wl", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "wl", "w.length")
  } else if (exists("wavelength", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "wavelength", "w.length")
  } else if (exists("Wavelength", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "Wavelength", "w.length")
  } else {
    warning("No wavelength data found in generic.spct")
    x[ , w.length := NA]
  }
  wl.min <- min(x$w.length, na.rm = TRUE)
#  wl.max <- max(x$w.length, na.rm = TRUE)
  if (wl.min == Inf) {
    warning("No valid 'w.length' values, probably a spectrum of length zero")
  } else if (wl.min < 99.999 || wl.min > 5e3) {
    stop("Off-range minimum w.length value ", wl.min, " instead of within 100 nm and 5000 nm")
  }
  return(x)
}

#' @describeIn check Specialization for filter.spct.
#' @export
check.filter.spct <- function(x, byref=TRUE, strict.range = TRUE) {

  range_check <- function(x, strict.range) {
    Tfr.min <- min(x$Tfr, na.rm = TRUE)
    Tfr.max <- max(x$Tfr, na.rm = TRUE)
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
    warning("Found varaible 'transmittance', I am assuming it expressed as percent")
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
    x[ , Tfr := Tpc / 100]
    x[ , Tpc := NULL]
    range_check(x, strict.range=strict.range)
    return(x)
  } else if (exists("A", x, mode = "numeric", inherits=FALSE)) {
#    x[ , Tfr := A2T(A)]
    if (min(x$A, na.rm = TRUE) < 0) {
      warning("Off-range min absorbance value: ", signif(min(x$A, na.rm = TRUE), 2), " instead of 0")
    }
    return(x)
  } else {
    warning("No transmittance or absorbance data found in filter.spct")
    x[ , Tfr := NA]
    return(x)
  }
}

#' @describeIn check Specialization for reflector.spct.
#' @export
check.reflector.spct <- function(x, byref=TRUE, strict.range = TRUE) {

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
    x[ , Rfr := Rpc / 100]
    x[ , Rpc := NULL]
    range_check(x, strict.range=strict.range)
    return(x)
  } else {
    warning("No reflectance data found in reflector.spct")
    x[ , Rfr := NA]
    return(x)
  }
}

#' @describeIn check Specialization for object.spct.
#' @export
check.object.spct <- function(x, byref=TRUE, strict.range = TRUE) {

  range_check <- function(x, strict.range) {
    Rfr.min <- min(x$Rfr, na.rm = TRUE)
    Rfr.max <- max(x$Rfr, na.rm = TRUE)
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
    Tfr.min <- min(x$Tfr, na.rm = TRUE)
    Tfr.max <- max(x$Tfr, na.rm = TRUE)
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
    x[ , Rfr := Rpc / 100]
    x[ , Rpc := NULL]
  } else {
    warning("No reflectance data found in object.spct")
    x[ , Rfr := NA]
  }

  if (exists("transmittance", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "transmittance", "Tpc")
    warning("Found varaible 'transmittance', I am assuming it expressed as percent")
  }
  if (exists("Tfr", x, mode = "numeric", inherits=FALSE)) {
    range_check(x, strict.range=strict.range)
    return(x)
  } else if (exists("Tpc", x, mode = "numeric", inherits=FALSE)) {
    x[ , Tfr := Tpc / 100]
    x[ , Tpc := NULL]
    range_check(x, strict.range=strict.range)
    return(x)
  } else {
    warning("No transmittance data found in object.spct")
    x[ , Tfr := NA]
    return(x)
  }
  range_check(x, strict.range=strict.range)
  return(x)
}

#' @describeIn check Specialization for response.spct.
#' @export
check.response.spct <- function(x, byref=TRUE, strict.range=TRUE) {
  if (exists("s.e.response", x, mode = "numeric", inherits=FALSE)) {
    return(x)
  } else if (exists("s.q.response", x, mode = "numeric", inherits=FALSE)) {
    return(x)
  } else if (exists("response", x, mode = "numeric", inherits=FALSE)) {
    x[ , s.e.response := response]
    x[ , response := NULL]
    warning("Found variable 'response', I am assuming it is expressed on an energy basis")
    return(x)
  } else if (exists("signal", x, mode = "numeric", inherits=FALSE)) {
    x[ , s.e.response := signal]
    x[ , signal := NULL]
    warning("Found variable 'signal', I am assuming it is expressed on an energy basis")
    return(x)
  } else {
    warning("No response data found in response.spct")
    x[ , s.e.response := NA]
    return(x)
  }
}

#' @describeIn check Specialization for source.spct.
#' @export
check.source.spct <- function(x, byref=TRUE, strict.range=FALSE) {

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

  if (is.null(getTimeUnit(x))) {
    setTimeUnit(x, "second")
    warning("Missing attribute 'time.unit' set to 'second'")
  }
  if (is.null(is.effective(x))) {
    setBSWFUsed(x, "none")
    warning("Missing atrribute 'bswf.used' set to 'none'")
  }
  if (exists("s.e.irrad", x, mode = "numeric", inherits=FALSE)) {
    NULL
  } else if (exists("s.q.irrad", x, mode = "numeric", inherits=FALSE)) {
    NULL
  } else if (exists("irradiance", x, mode = "numeric", inherits=FALSE)) {
    x[ , s.e.irradiance := irradiance]
    x[ , irradiance := NULL]
    warning("Found variable 'irradiance', I am assuming it is expressed on an energy basis")
  } else {
    warning("No spectral irradiance data found in source.spct")
    x[ , s.e.irrad := NA]
    return(x)
  }
  range_check(x, strict.range = strict.range)
  return(x)
}

#' @describeIn check Specialization for chroma.spct.
#' @export
check.chroma.spct <- function(x, byref=TRUE, strict.range=TRUE) {
  names_x <- names(x)
  idxs <- grep("[XYZ]", names_x)
  names2lc <- names_x[idxs]
  setnames(x, names2lc, tolower(names2lc))
  if (exists("x", x, mode="numeric", inherits=FALSE) &&
        exists("y", x, mode="numeric", inherits=FALSE) &&
        exists("z", x, mode="numeric", inherits=FALSE) ) {
    return(x)
  } else {
    warning("No spectral chromaticity coordinates data found in chroma.spct")
    return(x[ , c(x, y, z) := NA])
  }
}


# set class ---------------------------------------------------------------

#' Remove generic.spct and derived xxx.spct class attributes from a spectrum object.
#'
#' Removes from the class attibute of a xxxx.spct atributes.
#'
#' @param x An R object.
#' @export
#'
rmDerivedSpct <- function(x) {
  name <- substitute(x)
  spctclasses <- spct.classes()
  allclasses <- class(x)
  setattr(x, "class", setdiff(allclasses, spctclasses))
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' Set class of a data.frame or data.table object to generic.spct.
#'
#' Sets the class attibute of a data.frame or data.table object to "generic.spct" an object to store spectra.
#' If the object is a data.frame is is made a data.table in the process.
#'
#' @param x a data.frame or data.table
#' @export setGenericSpct setGenSpct
#' @exportClass generic.spct
#' @aliases setGenericSpct setGenSpct
#' @family setSpct functions
setGenSpct <- function(x) {
  name <- substitute(x)
  rmDerivedSpct(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.generic.spct(x)){
    setattr(x, "class", c("generic.spct", class(x)))
    setattr(x, "spct.tags", NA)
  }
  x <- check(x)
  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

setGenericSpct <- setGenSpct

#' set class of a data.frame or data.table object to "private.spct".
#'
#' Sets the class attibute of a data.frame or data.table object to "generic.spct" an object to store spectra.
#' If the object is a data.frame is is made a data.table in the process.
#'
#' @param x a data.frame or data.table
#' @family setSpct functions
#' @keywords internal
setPrivateSpct <- function(x) {
  name <- substitute(x)
  rmDerivedSpct(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.generic.spct(x)) {
    setGenSpct(x)
  }
  if (!is.private.spct(x)) {
    setattr(x, "class", c("private.spct", class(x)))
  }
  x <- check(x)
  setkey_spct(x, w.length)
  x <- copy(x[,list(w.length,numbers)])
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}


#' set class of a data.frame or data.table or generic.spct object to "filter.spct".
#'
#' Sets the class attibute of a data.frame or data.table object to "filter.spct" an object to store spectra.
#' If the object is a data.frame is is also made a data.table in the process.
#'
#' @param x a data.frame or data.table
#' @param Tfr.type a character string, either "total" or "internal"
#' @param strict.range logical indicating whether off-range values result in an error instead of a warning
#' @export
#' @exportClass filter.spct
#' @family setSpct functions
#'
setFilterSpct <- function(x, Tfr.type=c("total", "internal"), strict.range = TRUE) {
  name <- substitute(x)
  if ((is.object.spct(x) || is.filter.spct(x)) && getTfrType(x) != "unknown") {
    if (length(Tfr.type) > 1) {
      Tfr.type <- getTfrType(x)
    } else {
      warning("Replacing existing attribute 'Tfr.type' ", getTfrType(x))
    }
  }
  rmDerivedSpct(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.generic.spct(x)) {
    setGenSpct(x)
  }
  if (!is.filter.spct(x)) {
    setattr(x, "class", c("filter.spct", class(x)))
  }
  setTfrType(x, Tfr.type[1])
  x <- check(x, strict.range=strict.range)
  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "reflector.spct".
#'
#' Sets the class attibute of a data.frame or data.table object to "reflector.spct" an object to store spectra.
#' If the object is a data.frame is is also made a data.table in the process.
#'
#' @param x a data.frame or data.table
#' @param Rfr.type a character string, either "total" or "specular"
#' @param strict.range logical indicating whether off-range values result in an error instead of a warning
#' @export
#' @exportClass reflector.spct
#' @family setSpct functions
#'
setReflectorSpct <- function(x, Rfr.type=c("total", "specular"), strict.range = TRUE) {
  name <- substitute(x)
  if ((is.object.spct(x) || is.reflector.spct(c)) && getRfrType(x) != "unknown") {
    if (length(Rfr.type) > 1) {
      Rfr.type <- getRfrType(x)
    } else {
      warning("Replacing existing attribute 'Rfr.type' ", getRfrType(x))
    }
  }
  rmDerivedSpct(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.generic.spct(x)) {
    setGenSpct(x)
  }
  if (!is.reflector.spct(x)) {
    setattr(x, "class", c("reflector.spct", class(x)))
  }
  setRfrType(x, Rfr.type[1])
  x <- check(x, strict.range=strict.range)
  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "object.spct".
#'
#' Sets the class attibute of a data.frame or data.table object to "reflector.spct" an object to store spectra.
#' If the object is a data.frame is is also made a data.table in the process.
#'
#' @param x a data.frame or data.table
#' @param Tfr.type a character string, either "total" or "internal"
#' @param Rfr.type a character string, either "total" or "specular"
#' @param strict.range logical indicating whether off-range values result in an error instead of a warning
#' @export
#' @exportClass object.spct
#' @family setSpct functions
#'
setObjectSpct <- function(x, Tfr.type=c("total", "internal"),
                             Rfr.type=c("total", "specular"), strict.range = TRUE) {
  name <- substitute(x)
  if ((is.filter.spct(x) || is.object.spct(x)) && getTfrType(x) != "unknown") {
    if (length(Tfr.type) > 1) {
      Tfr.type <- getTfrType(x)
    } else {
      warning("Replacing existing attribute 'Tfr.type' ", getTfrType(x))
    }
  }
  if ((is.reflector.spct(x) || is.object.spct(x)) && getRfrType(x) != "unknown") {
    if (length(Rfr.type) > 1) {
      Rfr.type <- getRfrType(x)
    } else {
      warning("Replacing existing attribute 'Rfr.type' ", getRfrType(x))
    }
  }
  rmDerivedSpct(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.generic.spct(x)) {
    setGenSpct(x)
  }
  if (!is.object.spct(x)) {
    setattr(x, "class", c("object.spct", class(x)))
  }
  setTfrType(x, Tfr.type)
  setRfrType(x, Rfr.type)
  x <- check(x, strict.range=strict.range)
  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "response.spct".
#'
#' Sets the class attibute of a data.frame or data.table object to "response.spct" object,
#' used to store response spectra
#'
#' @param x a data.frame or data.table.
#' @param time.unit character string "second" or "day"
#' @export
#' @exportClass response.spct
#' @family setSpct functions
#'
setResponseSpct <- function(x, time.unit="none") {
  name <- substitute(x)
  rmDerivedSpct(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.generic.spct(x)) {
    setGenSpct(x)
  }
  if (!is.response.spct(x)) {
    setattr(x, "class", c("response.spct", class(x)))
  }
  setTimeUnit(x, time.unit)
  x <- check(x)
  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "source.spct".
#'
#' Sets the class attibute of a data.frame or data.table object to "source.spct" an object to store spectra.
#' If the object is a data.frame is is also made a data.table in the process.
#'
#' @param x a data.frame or data.table
#' @param time.unit character string "second" or "day"
#' @param bswf.used a character string, either "none" or the name of a BSWF
#' @param strict.range logical indicating whether off-range values result in an error instead of a warning
#' @export
#' @exportClass source.spct
#' @family setSpct functions
#'
setSourceSpct <- function(x, time.unit="second", bswf.used=c("none", "unknown"),
                          strict.range = FALSE) {
  name <- substitute(x)
  rmDerivedSpct(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.generic.spct(x)) {
    setGenSpct(x)
  }
  if (!is.source.spct(x)) {
    setattr(x, "class", c("source.spct", class(x)))
  }
  setTimeUnit(x, time.unit)
  setBSWFUsed(x, bswf.used = bswf.used)
  x <- check(x, strict.range = strict.range)
  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "chroma.spct".
#'
#' Sets the class attibute of a data.frame or data.table object to "chroma.spct" an object to store spectra.
#' If the object is a data.frame is is also made a data.table in the process.
#'
#' @param x a data.frame or data.table
#' @export
#' @exportClass chroma.spct
#' @family setSpct functions
#'
setChromaSpct <- function(x) {
  name <- substitute(x)
  rmDerivedSpct(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.generic.spct(x)) {
    setGenSpct(x)
  }
  if (!is.chroma.spct(x)) {
    setattr(x, "class", c("chroma.spct", class(x)))
  }
  x <- check(x)
  setkey_spct(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}


# is functions for spct classes --------------------------------------------

#' Query class of spectrum objects.
#'
#' Functions to check if an object is a given type of spectrum, or coerce it if possible.
#'
#' @usage is.generic.spct(x)
#'
#' @param x An R object.
#'
#' @return These functions return TRUE if its argument is a of the queried type of spectrum and FALSE otherwise.
#'
#' @method is generic.spct
#'
#' @export is.generic.spct
#'
is.generic.spct <- function(x) inherits(x, "generic.spct")

#' Query if an R object is a \code{private.spct}
#'
#' Functions to check if an object is a private spectrum, or coerce it if possible.
#'
#' @return is.private.spct returns TRUE if its argument is a private spectrum  and FALSE otherwise.
#' @export
#' @keywords internal
#'
is.private.spct <- function(x) inherits(x, "private.spct")

#' @describeIn is.generic.spct Query if an R object is a \code{source.spct}
#' @export
#'
is.source.spct <- function(x) inherits(x, "source.spct")

#' @describeIn is.generic.spct Query if an R object is a \code{filter.spct}
#' @export
#'
is.filter.spct <- function(x) inherits(x, "filter.spct")

#' @describeIn is.generic.spct Query if an R object is a \code{generic.spct}
#' @export
#'
is.reflector.spct <- function(x) inherits(x, "reflector.spct")

#' @describeIn is.generic.spct Query if an R object is an \code{object.spct}
#' @export
#'
is.object.spct <- function(x) inherits(x, "object.spct")

#' @describeIn is.generic.spct Query if an R object is a \code{response.spct}
#' @export
#'
is.response.spct <- function(x) inherits(x, "response.spct")

#' @describeIn is.generic.spct Query if an R object is a \code{chrome.spct}
#' @export
#'
is.chroma.spct <- function(x) inherits(x, "chroma.spct")

#' @describeIn is.generic.spct Query if an R object is an spectrum of any of the above classes.
#' @export
#'
is.any.spct <- function(x) {
  inherits(x, spct.classes())
}

#' Query which is the class of an spectrum
#'
#' Functions to check if an object is a generic spectrum, or coerce it if possible.
#'
#' @usage class.spct(x)
#'
#' @param x any R object
#'
#' @return class.spct returns a vector containing all matching xxxx.spct classes.
#'
#' @export
#'
class.spct <- function(x) {
#  intersect(spct.classes(), class(x)) # alters order!
  class(x)[class(x) %in% spct.classes()] # maintains order
}

#' Query if it is an spectrum is tagged
#'
#' Functions to check if an spct object contains tags.
#'
#' @usage is.tagged(x)
#'
#' @param x any R object
#'
#' @return is.tagged returns TRUE if its argument is a an spectrum
#' that contains tags and FALSE if it is an untagged spectrun, but
#' returns NA for any other R object.
#'
#' @export
#'
is.tagged <- function(x) {
  if (!is.any.spct(x)) {
    return(NA)
  } else {
    tags <- attr(x, "spct.tags", exact=TRUE)
    return(!is.null(tags) && length(tags) > 0 && !is.na(tags[[1]]))
  }
}

# is.photon.based ---------------------------------------------------------

#' Query if it is an spectrum contains photon based data
#'
#' Functions to check if an spct contains photon based data
#'
#' @usage is.photon.based(x)
#'
#' @param x any R object
#'
#' @return is.photon.based returns TRUE if its argument is a a source.spct or
#' a response.spct object that contains photon base data and FALSE if such an
#' object does not contain such data, but returns NA for any other R object,
#' including those belinging other \code{generic.spct}-derived classes.
#'
#' @export
#' @family query units functions
#'
is.photon.based <- function(x) {
  if (is.source.spct(x)) {
    return("s.q.irrad" %in% names(x))
  } else if (is.response.spct(x)) {
    return("s.q.response" %in% names(x))
  } else {
    return(NA)
  }
}

# is.energy.based ---------------------------------------------------------

#' Query if it is an spectrum contains energy based data
#'
#' Functions to check if an spct contains energy based data
#'
#' @usage is.energy.based(x)
#'
#' @param x any R object
#'
#' @return is.energy.based returns TRUE if its argument is a a source.spct or
#' a response.spct object that contains energy base data and FALSE if such an
#' object does not contain such data, but returns NA for any other R object,
#' including those belinging other \code{generic.spct}-derived classes
#'
#' @export
#' @family query units functions
#'
is.energy.based <- function(x) {
  if (is.source.spct(x)) {
    return("s.e.irrad" %in% names(x))
  } else if (is.response.spct(x)) {
    return("s.e.response" %in% names(x))
  } else {
    return(NA)
  }
}

# is.absorbance.based ---------------------------------------------------------

#' Query if it is an spectrum contains absorbance data
#'
#' Functions to check if an spct contains absorbance data
#'
#' @usage is.absorbance.based(x)
#'
#' @param x any R object
#'
#' @return \code{is.absorbance.based} returns TRUE if its argument is a a \code{filter.spct}
#' object that contains spectral absorbance data and FALSE if it does not contain
#' such data, but returns NA for any other R object, including those belonging
#' other \code{generic.spct}-derived classes.
#'
#' @export
#' @family query units functions
#'
is.absorbance.based <- function(x) {
  if (is.filter.spct(x)) {
    return("A" %in% names(x))
  } else {
    return(NA)
  }
}

# is.transmittance.based ---------------------------------------------------------

#' Query if it is an spectrum contains transmittance data
#'
#' Functions to check if an spct contains transmittance data
#'
#' @usage is.transmittance.based(x)
#'
#' @param x any R object
#'
#' @return \code{is.transmittance.based} returns TRUE if its argument is a a \code{filter.spct}
#' object that contains spectral transmittance data and FALSE if it does not contain
#' such data, but returns NA for any other R object, including those belonging
#' other \code{generic.spct}-derived classes.
#'
#' @export
#' @family query units functions
#'
is.transmittance.based <- function(x) {
  if (is.filter.spct(x)) {
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
#' @usage as.generic.spct(x)
#'
#' @param x An R object
#'
#' @return These functions return a copy of \code{x} converted into a given
#'   class of spectral object, if \code{x} is a valid argument to the
#'   correcponding set function.
#'
#' @export
#'
#' @family creation of spectral objects functions
#'
as.generic.spct <- function(x) {
  y <- copy(x)
  setGenericSpct(y)
}

#' Return a copy of an R object with its class set to private.spct
#'
#' Function that returns a converted copy of a spectrum object.
#'
#' @usage as.private.spct(x)
#'
#' @param x any R object
#'
#' @return as.private.spct returns a "private.spct" if possible.
#'
#' @keywords internal
#'
as.private.spct <- function(x) {
  y <- copy(x)
  setPrivateSpct(y)
}

#' @describeIn as.generic.spct Copy and convert an R object into a \code{source.spct}
#'
#' @usage as.source.spct(x, time.unit=c("second", "day"), bswf.used=c("none",
#'   "unknown"), strict.range = FALSE)
#'
#' @param time.unit character string, "second" or "day"
#' @param strict.range logical indicating whether off-range values result in an
#'   error instead of a warning
#'
#' @export
#'
as.source.spct <- function(x,
                           time.unit=c("second", "day"),
                           bswf.used=c("none", "unknown"),
                           strict.range = FALSE) {
  y <- copy(x)
  setSourceSpct(y, time.unit, strict.range = strict.range, bswf.used = bswf.used)
}

#' @describeIn as.generic.spct Copy and convert an R object into a
#'   \code{filter.spct}
#'
#' @usage as.filter.spct(x, Tfr.type=c("total", "internal"), strict.range =
#'   TRUE)
#'
#' @param Tfr.type a character string, either "total" or "internal"
#' @param strict.range logical indicating whether off-range values result in an
#'   error instead of a warning
#'
#' @export
#'
as.filter.spct <- function(x, Tfr.type=c("total", "internal"), strict.range = TRUE) {
  y <- copy(x)
  setFilterSpct(y, Tfr.type, strict.range = strict.range)
}

#' @describeIn as.generic.spct Copy and convert an R object into a
#'   \code{reflector.spct}
#'
#' @usage as.reflector.spct(x, Rfr.type = c("total", "specular"), strict.range =
#'   TRUE)
#'
#' @param Rfr.type a character string, either "total" or "specular"
#' @param strict.range logical indicating whether off-range values result in an
#'   error instead of a warning
#'
#' @export
#'
as.reflector.spct <- function(x, Rfr.type = c("total", "specular"), strict.range = TRUE) {
  y <- copy(x)
  setReflectorSpct(y, Rfr.type = Rfr.type, strict.range = strict.range)
}

#' @describeIn as.generic.spct Copy and convert an R object into a
#'   \code{object.spct}
#'
#' @usage as.object.spct(x, Tfr.type=c("total", "internal"), Rfr.type=c("total",
#'   "specular"), strict.range = TRUE)
#'
#' @export
#'
as.object.spct <- function(x,
                           Tfr.type=c("total", "internal"),
                           Rfr.type=c("total", "specular"),
                           strict.range = TRUE) {
  y <- copy(x)
  setObjectSpct(y, Tfr.type = Tfr.type, Rfr.type = Rfr.type,
                strict.range = strict.range)
}

#' @describeIn as.generic.spct Copy and convert an R object into a
#'   \code{response.spct}
#'
#' @usage as.response.spct(x, time.unit = "none")
#'
#' @param time.unit character string "second" or "day"
#'
#' @export
#'
as.response.spct <- function(x, time.unit = "none") {
  y <- copy(x)
  setResponseSpct(y, time.unit = time.unit)
}

#' @describeIn as.generic.spct Copy and convert an R object into a
#'   \code{source.spct}
#'
#' @usage as.chroma.spct(x)
#'
#' @export
#'
as.chroma.spct <- function(x) {
  y <- copy(x)
  setChromaSpct(y)
}


# time.unit attribute -----------------------------------------------------

#' Set the "time.unit" attribute of an existing source.spct object
#'
#' Funtion to set by reference the "time.unit" attribute
#'
#' @usage setTimeUnit(x, time.unit=c("second", "hour", "day", "none"))
#'
#' @param x a source.spct object
#' @param time.unit a character string, either "second", "hour", "day", or "none"
#'
#' @return x
#'
#' @note if x is not a source.spct or response.spct object, x is not modified.
#' The behaviour of this function is 'unusual' in that the default for parameter
#' \code{time.unit} is used only if \code{x} does not already have this attribute set.
#' \code{time.unit = "hour"} is currently not fully supported.
#'
#' @export
#' @family time attribute functions
setTimeUnit <- function(x, time.unit=c("second", "hour", "day", "none")) {
  if (length(time.unit) > 1) {
    if (getTimeUnit(x) != "unknown") {
      time.unit <- getTimeUnit(x)
    } else {
      time.unit <- time.unit[[1]]
    }
  }
  if (is.source.spct(x) || is.response.spct(x)) {
    if  (!(time.unit %in% c("second", "hour", "day", "none", "unknown"))) {
      warning("Invalid 'time.unit' argument, only 'second', 'hour', 'day', and 'none' supported.")
      time.unit <- "unknown"
    }
    setattr(x, "time.unit", time.unit)
  }
  return(x)
}

#' Get the "time.unit" attribute of an existing source.spct object
#'
#' Funtion to read the "time.unit" attribute
#'
#' @usage getTimeUnit(x)
#'
#' @param x a source.spct object
#'
#' @return character string
#'
#' @note if x is not a \code{filter.spct} or a \code{response.spct} object, NA
#' is retruned
#'
#' @export
#' @family time attribute functions
#'
getTimeUnit <- function(x) {
  if (is.source.spct(x) || is.response.spct(x)) {
    time.unit <- attr(x, "time.unit", exact = TRUE)
    if (is.null(time.unit)) {
      # need to handle objects created with old versions
      time.unit <- "unknown"
    }
    return(time.unit[[1]])
  } else {
    return(NA)
  }
}


# bswf attribute -----------------------------------------------------

#' Set the "bswf.used" attribute of an existing source.spct object
#'
#' Funtion to set by reference the "time.unit" attribute
#'
#' @usage setBSWFUsed(x, bswf.used=c("none", "unknown"))
#'
#' @param x a source.spct object
#' @param bswf.used a character string, either "none" or the name of a BSWF
#'
#' @return x
#'
#' @note if x is not a source.spct, x is not modified.
#' The behaviour of this function is 'unusual' in that the default for parameter
#' \code{bswf.used} is used only if \code{x} does not already have this attribute set.
#' \code{time.unit = "hour"} is currently not fully supported.
#'
#' @export
#' @family BSWF attribute functions
#'
setBSWFUsed <- function(x, bswf.used=c("none", "unknown")) {
  if (is.null(bswf.used) || length(bswf.used) < 1) {
    bswf.used <- "none"
  }
  if (length(bswf.used) > 1) {
    if (is.effective(x)) {
      bswf.used <- getBSWFUsed(x)
    } else {
      bswf.used <- bswf.used[[1]]
    }
  }
  if (is.source.spct(x)) {
    if  (!(is.character(bswf.used))) {
      warning("Only character strings are valid vlues for 'bswf.used' argument")
      bswf.used <- "unknown"
    }
    setattr(x, "bswf.used", bswf.used)
  }
  return(x)
}

#' Get the "bswf.used" attribute of an existing source.spct object
#'
#' Funtion to read the "time.unit" attribute
#'
#' @usage getBSWFUsed(x)
#'
#' @param x a source.spct object
#'
#' @return character string
#'
#' @note if x is not a \code{source.spct} object, NA
#' is retruned
#'
#' @export
#' @family BSWF attribute functions
#'
getBSWFUsed <- function(x) {
  if (is.source.spct(x)) {
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

# is.effective.source.spct defined in file "waveband.class.r" to avoid the need
# of using colate to get the documentation in the correct order.

# Tfr.type attribute ------------------------------------------------------

#' Set the "Tfr.type" attribute of an existing filter.spct or object.spct object
#'
#' Funtion to set by reference the "Tfr.type" attribute
#'
#' @usage setTfrType(x, Tfr.type=c("total", "internal"))
#'
#' @param x a filter.spct or an object.spct object
#' @param Tfr.type a character string, either "total" or "internal"
#'
#' @return x
#'
#' @note if x is not a filter.spct or an object.spct object, x is not modified
#' The behaviour of this function is 'unusual' in that the default for parameter
#' \code{Tfr.type} is used only if \code{x} does not already have this attribute set.
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
  if (is.filter.spct(x) || is.object.spct(x)) {
    if  (!(Tfr.type %in% c("total", "internal", "unknown"))) {
      warning("Invalid 'Tfr.type' argument, only 'total' and 'internal' supported.")
      return(x)
    }
    setattr(x, "Tfr.type", Tfr.type)
  }
  return(x)
}

#' Get the "Tfr.type" attribute of an existing filter.spct or object.spct object
#'
#' Funtion to read the "Tfr.type" attribute
#'
#' @usage setTfrType(x, Tfr.type=c("total", "internal"))
#'
#' @param x a filter.spct or object.spct object
#'
#' @return character string
#'
#' @note if x is not a \code{filter.spct} or an \code{object.spct} object, \code{NA} is returned
#'
#' @export
#' @family Tfr attribute functions
#'
getTfrType <- function(x) {
  if (is.filter.spct(x) || is.object.spct(x)) {
    Tfr.type <- attr(x, "Tfr.type", exact = TRUE)
    if (is.null(Tfr.type)) {
      # need to handle objects created with old versions
      Tfr.type <- "unknown"
    }
    return(Tfr.type[[1]])
  } else {
    return(NA)
  }
}

# Rfr.type attribute ------------------------------------------------------

#' Set the "Rfr.type" attribute of an existing reflector.spct or object.spct object
#'
#' Funtion to set by reference the "Rfr.type" attribute
#'
#' @usage setRfrType(x, Rfr.type=c("total", "specular"))
#'
#' @param x a reflector.spct or an object.spct object
#' @param Rfr.type a character string, either "total" or "specular"
#'
#' @return x
#'
#' @note if x is not a reflector.spct or object.spct object, x is not modified.
#' The behaviour of this function is 'unusual' in that the default for parameter
#' Rfr.type is used only if \code{x} does not already have this attribute set.
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
  if (is.reflector.spct(x) || is.object.spct(x)) {
    if  (!(Rfr.type %in% c("total", "specular", "unknown"))) {
      warning("Invalid 'Rfr.type' argument, only 'total' and 'internal' supported.")
      return(x)
    }
    setattr(x, "Rfr.type", Rfr.type)
  }
  return(x)
}

#' Get the "Rfr.type" attribute of an existing reflector.spct object
#'
#' Funtion to read the "Rfr.type" attribute
#'
#' @usage getRfrType(x)
#'
#' @param x a source.spct object
#'
#' @return character string
#'
#' @note if x is not a \code{filter.spct} object, \code{NA} is returned
#'
#' @export
#' @family Rfr attribute functions
#'
getRfrType <- function(x) {
  if (is.reflector.spct(x) || is.object.spct(x)) {
    Rfr.type <- attr(x, "Rfr.type", exact = TRUE)
    if (is.null(Rfr.type)) {
      # need to handle objects created with old versions
      Rfr.type <- "unknown"
    }
    return(Rfr.type[[1]])
  } else {
    return(NA)
  }
}

