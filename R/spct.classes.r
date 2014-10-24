
# names of all spectral classes -------------------------------------------

#' Function that returns a vector containing the names of spectra classes
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
    "source.spct",
    "response.spct", "chroma.spct")
}

# check -------------------------------------------------------------------

#' Generic function
#'
#' Check that an R object contains the expected data members.
#'
#' @param x an R object
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export check
check <- function(x, byref) UseMethod("check")

#' Default for generic function
#'
#' Check that an R object contains the expected data members.
#'
#' @param x an R object
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export check.default
check.default <- function(x, byref=FALSE) {
  return(x)
}

#' Specialization for private.spct
#'
#' Check that an R object contains the expected data members.
#'
#' @param x a private.spct object
#' @param byref logical indicating if new object will be created by reference or by copy of x
check.private.spct <- function(x, byref=TRUE) {
  if (exists("w.length", x, mode = "numeric", inherits=FALSE) &&
        exists("numbers", x, mode = "numeric", inherits=FALSE)) {
    invisible(x)
  }
}

#' Specialization for generic.spct
#'
#' Check that a generic.spct object contains the expected data members.
#'
#' @param x a generic.spct object
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export check.generic.spct
check.generic.spct <- function(x, byref=TRUE) {
  if (exists("w.length", x, mode = "numeric", inherits=FALSE)) {
    invisible(x)
  } else if (exists("wl", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "wl", "w.length")
    invisible(x)
  } else if (exists("wavelength", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "wavelength", "w.length")
    invisible(x)
  } else {
    warning("No wavelength data found in generic.spct")
    x[ , w.length := NA]
    invisible(x)
  }
}

#' Specialization for filter.spct
#'
#' Check that an R object contains the expected data members.
#'
#' @param x an R object
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export check.filter.spct
check.filter.spct <- function(x, byref=TRUE) {
  if (exists("Tfr", x, mode = "numeric", inherits=FALSE)) {
    invisible(x)
  } else if (exists("Tpc", x, mode = "numeric", inherits=FALSE)) {
    x[ , Tfr := Tpc / 100]
    invisible(x)
  } else if (exists("A", x, mode = "numeric", inherits=FALSE)) {
    x[ , Tfr := A2T(A)]
    invisible(x)
  } else {
    warning("No transmittance or absorbance data found in filter.spct")
    x[ , Tfr := NA]
    invisible(x)
  }
}

#' Specialization for reflector.spct
#'
#' Check that a reflector.spct object contains the expected data members.
#'
#' @param x a reflector.spct object
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export check.reflector.spct
check.reflector.spct <- function(x, byref=TRUE) {
  if (exists("Rfr", x, mode = "numeric", inherits=FALSE)) {
    invisible(x)
  } else if (exists("Rpc", x, mode = "numeric", inherits=FALSE)) {
    x[ , Rfr := Rpc / 100]
    invisible(x)
  } else {
    warning("No reflectance data found in filter.spct")
    x[ , Rfr := NA]
    invisible(x)
  }
}

#' Specialization for response.spct
#'
#' Check that a response.spct object contains the expected data members.
#'
#' @param x a response.spct object
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export check.response.spct
check.response.spct <- function(x, byref=TRUE) {
  if (exists("s.e.response", x, mode = "numeric", inherits=FALSE)) {
    invisible(x)
  } else if (exists("response", x, mode = "numeric", inherits=FALSE)) {
    x[ , s.e.response := response]
    x[ , response := NULL]
    invisible(x)
  } else if (exists("signal", x, mode = "numeric", inherits=FALSE)) {
    x[ , s.e.response := signal]
    x[ , signal := NULL]
    invisible(x)
  } else if (exists("s.q.response", x, mode = "numeric", inherits=FALSE)) {
    invisible(q2e(x, action="add", byref=byref))
  } else {
    warning("No response data found in filter.spct")
    x[ , s.e.response := NA]
    invisible(x)
  }
}

#' Specialization for source.spct
#'
#' Check that a source.spct object contains the expected data members.
#'
#' @param x a source.spct object
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export check.source.spct
check.source.spct <- function(x, byref=TRUE) {
  if (exists("s.e.irrad", x, mode = "numeric", inherits=FALSE)) {
    invisible(x)
  } else if (exists("s.q.irrad", x, mode = "numeric", inherits=FALSE)) {
    invisible(q2e(x, action="add", byref=byref))
  } else {
    warning("No spectral irradiance data found in source.spct")
    x[ , s.e.irrad := NA]
    invisible(x)
  }
}

#' Specialization for chroma.spct
#'
#' Check that a chroma.spct object contains the expected data members.
#'
#' @param x a source.spct object
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export check.source.spct
check.chroma.spct <- function(x, byref=TRUE) {
  names_x <- names(x)
  idxs <- grep("[XYZ]", names_x)
  names2lc <- names_x[idxs]
  setnames(x, names2lc, tolower(names2lc))
  if (exists("x", x, mode="numeric", inherits=FALSE) &&
        exists("y", x, mode="numeric", inherits=FALSE) &&
        exists("z", x, mode="numeric", inherits=FALSE) ) {
    invisible(x)
  } else {
    warning("No spectral chromaticity coordinates data found in chroma.spct")
    invisible(x[ , c(x, y, z) := NA])
  }
}


# set class ---------------------------------------------------------------


#' set class of a data.frame or data.table object to "generic.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "generic.spct" an object to store spectra.
#' If the object is a data.frame is is made a data.table in the process.
#'
#' @param x a data.frame or data.table
#' @export setGenericSpct setGenSpct
#' @exportClass generic.spct
#' @aliases setGenericSpct setGenSpct
#'
setGenSpct <- function(x) {
  name <- substitute(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.generic.spct(x)){
    setattr(x, "class", c("generic.spct", class(x)))
  }
  x <- check(x)
  setkey(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

setGenericSpct <- setGenSpct

#' set class of a data.frame or data.table object to "private.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "generic.spct" an object to store spectra.
#' If the object is a data.frame is is made a data.table in the process.
#'
#' @param x a data.frame or data.table
#'
setPrivateSpct <- function(x) {
  name <- substitute(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.private.spct(x)) {
    setattr(x, "class", c("private.spct", class(x)))
  }
  x <- check(x)
  setkey(x, w.length)
  x <- copy(x[,list(w.length,numbers)])
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}


#' set class of a data.frame or data.table or generic.spct object to "filter.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "filter.spct" an object to store spectra.
#' If the object is a data.frame is is also made a data.table in the process.
#'
#' @param x a data.frame or data.table
#' @export
#' @exportClass filter.spct
#'
setFilterSpct <- function(x) {
  name <- substitute(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.generic.spct(x)) {
    setGenSpct(x)
  }
  if (!is.filter.spct(x)) {
    setattr(x, "class", c("filter.spct", class(x)))
  }
  x <- check(x)
  setkey(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "reflector.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "reflector.spct" an object to store spectra.
#' If the object is a data.frame is is also made a data.table in the process.
#'
#' @param x a data.frame or data.table
#' @export
#' @exportClass filter.spct
#'
setReflectorSpct <- function(x) {
  name <- substitute(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.generic.spct(x)) {
    setGenSpct(x)
  }
  if (!is.reflector.spct(x)) {
    setattr(x, "class", c("reflector.spct", class(x)))
  }
  x <- check(x)
  setkey(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "response.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "response.spct" object,
#' used to store response spectra
#'
#' @param x a data.frame or data.table.
#' @export
#' @exportClass filter.spct
#'
setResponseSpct <- function(x) {
  name <- substitute(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.generic.spct(x)) {
    setGenSpct(x)
  }
  if (!is.response.spct(x)) {
    setattr(x, "class", c("response.spct", class(x)))
  }
  x <- check(x)
  setkey(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "source.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "source.spct" an object to store spectra.
#' If the object is a data.frame is is also made a data.table in the process.
#'
#' @param x a data.frame or data.table
#' @param time.unit character string "second" or "day"
#' @export
#' @exportClass source.spct
#'
setSourceSpct <- function(x, time.unit="second") {
  name <- substitute(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is.generic.spct(x)) {
    setGenSpct(x)
  }
  if (!is.source.spct(x)) {
    setattr(x, "class", c("source.spct", class(x)))
  }
  x <- check(x)
  setattr(x, "time.unit", time.unit)
  setkey(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "chroma.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "chroma.spct" an object to store spectra.
#' If the object is a data.frame is is also made a data.table in the process.
#'
#' @param x a data.frame or data.table
#' @export
#' @exportClass chroma.spct
#'
setChromaSpct <- function(x) {
  name <- substitute(x)
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
  setkey(x, w.length)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}


# is function for spct classes --------------------------------------------

#' Query class of a generic spectrum
#'
#' Functions to check if an object is a generic spectrum, or coerce it if possible.
#'
#' @usage is.generic.spct(x)
#'
#' @param x any R object
#'
#' @return is.generic.spct returns TRUE if its argument is a generic spectrum (that is, has "generic.spct" amongst its classes) and FALSE otherwise.
#'
#' @export
#'
is.generic.spct <- function(x) inherits(x, "generic.spct")

#' Query class of a private spectrum
#'
#' Functions to check if an object is a private spectrum, or coerce it if possible.
#'
#' @usage is.private.spct(x)
#'
#' @param x any R object
#'
#' @return is.private.spct returns TRUE if its argument is a private spectrum (that is, has "private.spct" amongst its classes) and FALSE otherwise.
#'
#' @export
#'
is.private.spct <- function(x) inherits(x, "private.spct")

#' Query class of a source spectrum
#'
#' Functions to check if an object is a source spectrum, or coerce it if possible.
#'
#' @usage is.source.spct(x)
#'
#' @param x any R object
#'
#' @return is.source.spct returns TRUE if its argument is a source spectrum (that is, has "source.spct" amongst its classes) and FALSE otherwise.
#'
#' @export
#'
is.source.spct <- function(x) inherits(x, "source.spct")

#' Query class of a filter spectrum
#'
#' Functions to check if an object is a filter spectrum, or coerce it if possible.
#'
#' @usage is.filter.spct(x)
#'
#' @param x any R object
#'
#' @return is.filter.spct returns TRUE if its argument is a filter spectrum (that is, has "filter.spct" amongst its classes) and FALSE otherwise.
#'
#' @export
#'
is.filter.spct <- function(x) inherits(x, "filter.spct")

#' Query class of a reflector spectrum
#'
#' Functions to check if an object is a reflector spectrum, or coerce it if possible.
#'
#' @usage is.reflector.spct(x)
#'
#' @param x any R object
#'
#' @return is.reflector.spct returns TRUE if its argument is a reflector spectrum (that is, has "reflector.spct" amongst its classes) and FALSE otherwise.
#'
#' @export
#'
is.reflector.spct <- function(x) inherits(x, "reflector.spct")

#' Query class of a response spectrum
#'
#' Functions to check if an object is a response spectrum, or coerce it if possible.
#'
#' @usage is.response.spct(x)
#'
#' @param x any R object
#'
#' @return is.response.spct returns TRUE if its argument is a response spectrum (that is, has "response.spct" amongst its classes) and FALSE otherwise.
#'
#' @export
#'
is.response.spct <- function(x) inherits(x, "response.spct")

#' Query class of a chroma spectrum
#'
#' Functions to check if an object is a chroma spectrum, or coerce it if possible.
#'
#' @usage is.chroma.spct(x)
#'
#' @param x any R object
#'
#' @return is.chroma.spct returns TRUE if its argument is a chroma spectrum (that is, has "chroma.spct" amongst its classes) and FALSE otherwise.
#'
#' @export
#'
is.chroma.spct <- function(x) inherits(x, "chroma.spct")

#' Query if it is an spectrum
#'
#' Functions to check if an object is a generic spectrum, or coerce it if possible.
#'
#' @usage is.any.spct(x)
#'
#' @param x any R object
#'
#' @return is.any.spct returns TRUE if its argument is a an spectrum and FALSE otherwise.
#'
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
  spct.cls <- spct.classes()
  spct.cls[inherits(x, spct.cls, TRUE)]
}

