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

#' Specialization for generic.spct
#'
#' Check that a generic.spct object contains the expected data members.
#'
#' @param x a generic.spct object
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export check.generic.spct
check.generic.spct <- function(x, byref=TRUE) {
  if (exists("w.length", x, mode = "numeric", inherits=FALSE)) {
    invisible(return(x))
  } else if (exists("wl", x, mode = "numeric", inherits=FALSE)) {
    setnames(x, "wl", "w.length")
    invisible(return(x))
  } else {
    warning("No wavelength data found in generic.spct")
    x[ , w.length := NA]
    invisible(return(x))
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
    invisible(return(x))
  } else if (exists("Tpc", x, mode = "numeric", inherits=FALSE)) {
    x[ , Tfr := Tpc / 100]
    invisible(return(x))
  } else if (exists("A", x, mode = "numeric", inherits=FALSE)) {
    x[ , Tfr := A2T(A)]
    invisible(return(x))
  } else {
    warning("No transmittance or absorbance data found in filter.spct")
    x[ , Tfr := NA]
    invisible(return(x))
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
    invisible(return(x))
  } else if (exists("Rpc", x, mode = "numeric", inherits=FALSE)) {
    x[ , Rfr := Rpc / 100]
    invisible(return(x))
  } else {
    warning("No reflectance data found in filter.spct")
    x[ , Rfr := NA]
    invisible(return(x))
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
  if (exists("response", x, mode = "numeric", inherits=FALSE)) {
    invisible(return(x))
  } else if (exists("signal", x, mode = "numeric", inherits=FALSE)) {
    x[ , response := signal]
    x[ , signal := NULL]
    invisible(return(x))
  } else {
    warning("No response data found in filter.spct")
    x[ , response := NA]
    invisible(return(x))
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
    invisible(return(x))
  } else if (exists("s.q.irrad", x, mode = "numeric", inherits=FALSE)) {
    invisible(return(q2e(x, action="add", byref=byref)))
  } else {
    warning("No spectral irradiance data found in source.spct")
    x[ , s.e.irrad := NA]
    invisible(return(x))
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
    invisible(return(x))
  } else {
    warning("No spectral chromaticity coordinates data found in chroma.spct")
    invisible(return(x[ , c(x, y, z) := NA]))
  }
}


# set class ---------------------------------------------------------------


#' set class of a data.frame or data.table object to "generic.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "generic.spct" an object to store spectra
#' if the object is a data.frame is is mane a data.table
#'
#' @param x a data.frame
#' @export
#' @exportClass generic.spct
#'
setGenSpct <- function(x) {
  name = substitute(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  setattr(x, "class", c("generic.spct", class(x)))
  x <- check(x)
  setkey(x, w.length)
  if (is.name(name)) {
    name = as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "filter.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "filter.spct" an object to store spectra
#' if the object is a data.frame is is mane a data.table
#'
#' @param x a data.frame
#' @export
#' @exportClass filter.spct
#'
setFilterSpct <- function(x) {
  name = substitute(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is(x, "generic.spct")) {
    setGenSpct(x)
  }
  if (!is(x, "filter.spct")) {
    setattr(x, "class", c("filter.spct", class(x)))
  }
  x <- check(x)
  setkey(x, w.length)
  if (is.name(name)) {
    name = as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "reflector.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "reflector.spct" an object to store spectra
#' if the object is a data.frame is is mane a data.table
#'
#' @param x a data.frame
#' @export
#' @exportClass filter.spct
#'
setReflectorSpct <- function(x) {
  name = substitute(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is(x, "generic.spct")) {
    setGenSpct(x)
  }
  if (!is(x, "reflector.spct")) {
    setattr(x, "class", c("reflector.spct", class(x)))
  }
  x <- check(x)
  setkey(x, w.length)
  if (is.name(name)) {
    name = as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "response.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "response.spct" object,
#' used to store response spectra
#'
#' @param x a data.frame
#' @export
#' @exportClass filter.spct
#'
setResponseSpct <- function(x) {
  name = substitute(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is(x, "generic.spct")) {
    setGenSpct(x)
  }
  if (!is(x, "response.spct")) {
    setattr(x, "class", c("response.spct", class(x)))
  }
  x <- check(x)
  setkey(x, w.length)
  if (is.name(name)) {
    name = as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "source.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "source.spct" an object to store spectra
#' if the object is a data.frame is is mane a data.table
#'
#' @param x a data.frame
#' @export
#' @exportClass source.spct
#'
setSourceSpct <- function(x) {
  name = substitute(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is(x, "generic.spct")) {
    setGenSpct(x)
  }
  if (!is(x, "sourcer_spct")) {
    setattr(x, "class", c("source.spct", class(x)))
  }
  x <- check(x)
  setkey(x, w.length)
  if (is.name(name)) {
    name = as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' set class of a data.frame or data.table or generic.spct object to "chroma.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "chroma.spct" an object to store spectra
#' if the object is a data.frame is is mane a data.table
#'
#' @param x a data.frame
#' @export
#' @exportClass chroma.spct
#'
setChromaSpct <- function(x) {
  name = substitute(x)
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is(x, "generic.spct")) {
    setGenSpct(x)
  }
  if (!is(x, "chroma.spct")) {
    setattr(x, "class", c("chroma.spct", class(x)))
  }
  x <- check(x)
  setkey(x, w.length)
  if (is.name(name)) {
    name = as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}


# multiplication ----------------------------------------------------------

#' "*" operator for generic spectra
#'
#' Multiplication operator for generic spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 an object of class "generic.spct"
#' @name times-.generic.spct
#' @export
#'
'*.generic.spct' <- function(e1, e2) {
  if (is(e1, "chroma.spct")) {
    if (is(e2, "source.spct")) {
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$s.e.irrad, bin.oper=`*`, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$s.e.irrad, bin.oper=`*`, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$s.e.irrad, bin.oper=`*`, trim="intersection")
      out.spct <- data.frame(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]])
      setChromaSpct(out.spct)
      return(out.spct)
    } else if (is(e2, "filter.spct")) {
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$Tfr, bin.oper=`*`, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$Tfr, bin.oper=`*`, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$Tfr, bin.oper=`*`, trim="intersection")
      out.spct <- data.frame(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]])
      setChromaSpct(out.spct)
      return(out.spct)
    } else if (is(e2, "reflector.spct")) {
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$Rfr, bin.oper=`*`, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$Rfr, bin.oper=`*`, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$Rfr, bin.oper=`*`, trim="intersection")
      out.spct <- data.frame(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]])
      setChromaSpct(out.spct)
      return(out.spct)
    } else if (is.numeric(e2)) {
      e3 <- copy(e1)
      if (length(e2) == 3 && names(e2) == c("x", "y", "z")) {
        e3[ , `:=`(x = x * e2["x"] , y = y * e2["y"], z = z  * e2["z"])]
        return(e3)
      } else {
        e3[ , `:=`(x = x * e2 , y = y * e2, z = z  * e2)]
        return(e3)
      }
    } else {
      return(NA)
    }
  } else if (is(e1, "filter.spct")) {
    if (is(e2, "filter.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=`*`, trim="intersection")
      setnames(z, 2, "Tfr")
      setFilterSpct(z)
      return(z)
    } else if(is(e2, "source.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.e.irrad, bin.oper=`*`, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , Tfr := Tfr * e2]
      if (exists("Tpc", z, inherits=FALSE)) {
        z[ , Tpc := NULL]
      }
      if (exists("A", z, inherits=FALSE)) {
        z[ , A := NULL]
      }
      return(z)
    } else {
      return(NA)
    }
  } else if(is(e1, "reflector.spct")) {
    if (is(e2, "reflector.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$Rfr, bin.oper=`*`, trim="intersection")
      setnames(z, 2, "Rfr")
      setReflectorSpct(z)
      return(z)
    } else if(is(e2, "source.spct")) {
      if (!exists("s.e.irrad", e2, inherits=FALSE)) {
        q2e(e2, "replace")
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$s.e.irrad, bin.oper=`*`, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , Rfr := Rfr * e2]
      if (exists("Rpc", z, inherits=FALSE)) {
        z[ , Rpc := NULL]
      }
      return(z)
    } else {
      return(NA)
    }
  } else if (is(e1, "response.spct")) {
    if (is(e2, "filter.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$response, e2$Tfr, bin.oper=`*`, trim="intersection")
      setnames(z, 2, "response")
      setResponseSpct(z)
      return(z)
    } else if (is(e2, "reflector.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$response, e2$Rfr, bin.oper=`*`, trim="intersection")
      setnames(z, 2, "response")
      setResponseSpct(z)
      return(z)
    }else if(is(e2, "source.spct")) {
      if (!exists("s.e.irrad", e2, inherits=FALSE)) {
        q2e(e2, "replace")
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$response, e2$s.e.irrad, bin.oper=`*`, trim="intersection")
      setnames(z, 2, "response")
      setResponseSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , response := response * e2]
      return(z)
    } else {
      return(NA)
    }
  } else if (is(e1, "source.spct")) {
    if (is(e2, "source.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=`*`, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      return(z)
    } else if (is(e2, "filter.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$Tfr, bin.oper=`*`, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      return(z)
    } else if (is(e2, "reflector.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$Rfr, bin.oper=`*`, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      return(z)
    } else if (is(e2, "response.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$response, bin.oper=`*`, trim="intersection")
      setnames(z, 2, "response")
      setResponseSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , s.e.irrad := s.e.irrad * e2]
      if (exists("s.q.irrad", z, inherits=FALSE)) {
        z[ , s.q.irrad := NULL]
      }
      return(z)
    } else if (is(e2, "waveband")) {
      z <- copy(e1)
      z$s.e.irrad <- z$s.e.irrad * calc_multipliers(z$w.length, e2, unit.out="energy", unit.in="energy")
      if (exists("s.q.irrad", z, inherits=FALSE)) {
        z[ , s.q.irrad := NULL]
      }
      return(z)
    } else {
      return(NA)
    }
  }
}

# division ----------------------------------------------------------------


#' "/" operator for generic spectra
#'
#' Division operator for generic spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 an object of class "generic.spct"
#' @name slash-.generic.spct
#' @export
#'
'/.generic.spct' <- function(e1, e2) {
  if (is(e1, "chroma.spct")) {
    if (is(e2, "source.spct")) {
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$s.e.irrad, bin.oper=`/`, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$s.e.irrad, bin.oper=`/`, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$s.e.irrad, bin.oper=`/`, trim="intersection")
      out.spct <- data.frame(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]])
      setChromaSpct(out.spct)
      return(out.spct)
    } else if (is(e2, "filter.spct")) {
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$Tfr, bin.oper=`/`, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$Tfr, bin.oper=`/`, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$Tfr, bin.oper=`/`, trim="intersection")
      out.spct <- data.frame(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]])
      setChromaSpct(out.spct)
      return(out.spct)
    } else if (is(e2, "reflector.spct")) {
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$Rfr, bin.oper=`/`, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$Rfr, bin.oper=`/`, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$Rfr, bin.oper=`/`, trim="intersection")
      out.spct <- data.frame(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]])
      setChromaSpct(out.spct)
      return(out.spct)
    } else if (is.numeric(e2)) {
      e3 <- copy(e1)
      if (length(e2) == 3 && names(e2) == c("x", "y", "z")) {
        e3[ , `:=`(x = x / e2["x"] , y = y / e2["y"], z = z  / e2["z"])]
        return(e3)
      } else {
        e3[ , `:=`(x = x / e2 , y = y / e2, z = z  / e2)]
        return(e3)
      }
    } else {
      return(NA)
    }
  } else if (is(e1, "filter.spct")) {
    if (is(e2, "filter.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=`/`, trim="intersection")
      setnames(z, 2, "Tfr")
      setFilterSpct(z)
      return(z)
    } else if(is(e2, "source.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.e.irrad, bin.oper=`/`, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , Tfr := Tfr / e2]
      if (exists("Tpc", z, inherits=FALSE)) {
        z[ , Tpc := NULL]
      }
      if (exists("A", z, inherits=FALSE)) {
        z[ , A := NULL]
      }
      return(z)
    } else {
      return(NA)
    }
  } else if(is(e1, "reflector.spct")) {
    if (is(e2, "reflector.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$Rfr, bin.oper=`/`, trim="intersection")
      setnames(z, 2, "Rfr")
      setReflectorSpct(z)
      return(z)
    } else if(is(e2, "source.spct")) {
      if (!exists("s.e.irrad", e2, inherits=FALSE)) {
        q2e(e2, "replace")
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$s.e.irrad, bin.oper=`/`, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , Rfr := Rfr / e2]
      if (exists("Rpc", z, inherits=FALSE)) {
        z[ , Rpc := NULL]
      }
      return(z)
    } else {
      return(NA)
    }
  } else if (is(e1, "response.spct")) {
    if (is(e2, "filter.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$response, e2$Tfr, bin.oper=`/`, trim="intersection")
      setnames(z, 2, "response")
      setResponseSpct(z)
      return(z)
    } else if (is(e2, "reflector.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$response, e2$Rfr, bin.oper=`/`, trim="intersection")
      setnames(z, 2, "response")
      setResponseSpct(z)
      return(z)
    }else if(is(e2, "source.spct")) {
      if (!exists("s.e.irrad", e2, inherits=FALSE)) {
        q2e(e2, "replace")
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$response, e2$s.e.irrad, bin.oper=`/`, trim="intersection")
      setnames(z, 2, "response")
      setResponseSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , response := response / e2]
      return(z)
    } else {
      return(NA)
    }
  } else if (is(e1, "source.spct")) {
    if (is(e2, "source.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=`/`, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      return(z)
    } else if (is(e2, "filter.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$Tfr, bin.oper=`/`, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      return(z)
    } else if (is(e2, "reflector.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$Rfr, bin.oper=`/`, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      return(z)
    } else if (is(e2, "response.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$response, bin.oper=`/`, trim="intersection")
      setnames(z, 2, "response")
      setResponseSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , s.e.irrad := s.e.irrad / e2]
      if (exists("s.q.irrad", z, inherits=FALSE)) {
        z[ , s.q.irrad := NULL]
      }
      return(z)
    } else if (is(e2, "waveband")) {
      z <- copy(e1)
      z$s.e.irrad <- z$s.e.irrad / calc_multipliers(z$w.length, e2, unit.out="energy", unit.in="energy")
      if (exists("s.q.irrad", z, inherits=FALSE)) {
        z[ , s.q.irrad := NULL]
      }
      return(z)
    } else {
      return(NA)
    }
  }
}


# Sum ---------------------------------------------------------------

#' "+" operator for generic spectra
#'
#' Division operator for generic spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 an object of class "generic.spct"
#' @name plus-.generic.spct
#' @export
#'
'+.generic.spct' <- function(e1, e2) {
  if (is(e1, "chroma.spct")) {
    if (is(e2, "filter.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=`+`, trim="intersection")
      setnames(z, 2, "Tfr")
      setFilterSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , Tfr := Tfr + e2]
      if (exists("Tpc", z, inherits=FALSE)) {
        z[ , Tpc := NULL]
      }
      if (exists("A", z, inherits=FALSE)) {
        z[ , A := NULL]
      }
      return(z)
    } else {
      return(NA)
    }
  } else if (is(e1, "filterr.spct")) {
    if (is(e2, "filter.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=`+`, trim="intersection")
      setnames(z, 2, "Tfr")
      setFilterSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , Tfr := Tfr + e2]
      if (exists("Tpc", z, inherits=FALSE)) {
        z[ , Tpc := NULL]
      }
      if (exists("A", z, inherits=FALSE)) {
        z[ , A := NULL]
      }
      return(z)
    } else {
      return(NA)
    }
  } else if (is(e1, "reflector.spct")) {
    if (is(e2, "reflector.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$Rfr, bin.oper=`+`, trim="intersection")
      setnames(z, 2, "Rfr")
      setReflectorSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , Rfr := Rfr + e2]
      if (exists("Rpc", z, inherits=FALSE)) {
        z[ , Rpc := NULL]
      }
      return(z)
    } else {
      return(NA)
    }
  } else if (is(e1, "response.spct")) {
    if(is(e2, "response.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$response, e2$response, bin.oper=`+`, trim="intersection")
      setnames(z, 2, "response")
      setResponseSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , response := response + e2]
      return(z)
    } else {
      return(NA)
    }
  } else if (is(e1, "source.spct")) {
    if (!exists("s.e.irrad", e1, inherits=FALSE)) {
      q2e(e2, "replace")
    }

    if (is(e2, "source.spct")) {
      if (!exists("s.e.irrad", e2, inherits=FALSE)) {
        q2e(e2, "replace")
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=`+`, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , s.e.irrad := s.e.irrad + e2]
      if (exists("s.q.irrad", z, inherits=FALSE)) {
        z[ , s.q.irrad := NULL]
      }
      return(z)
    }  else {
      return(NA)
    }
  }
}


# Minus -------------------------------------------------------------------

#' "-" operator for generic spectra
#'
#' Substraction operator for generic spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 an object of class "generic.spct"
#' @name minus-.generic.spct
#' @export
#'
'-.generic.spct' <- function(e1, e2) {
  if (is(e1, "chroma.spct")) {
    if (is(e2, "filter.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=`-`, trim="intersection")
      setnames(z, 2, "Tfr")
      setFilterSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , Tfr := Tfr - e2]
      if (exists("Tpc", z, inherits=FALSE)) {
        z[ , Tpc := NULL]
      }
      if (exists("A", z, inherits=FALSE)) {
        z[ , A := NULL]
      }
      return(z)
    } else {
      return(NA)
    }
  } else if (is(e1, "filterr.spct")) {
    if (is(e2, "filter.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=`-`, trim="intersection")
      setnames(z, 2, "Tfr")
      setFilterSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , Tfr := Tfr - e2]
      if (exists("Tpc", z, inherits=FALSE)) {
        z[ , Tpc := NULL]
      }
      if (exists("A", z, inherits=FALSE)) {
        z[ , A := NULL]
      }
      return(z)
    } else {
      return(NA)
    }
  } else if (is(e1, "reflector.spct")) {
    if (is(e2, "reflector.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$Rfr, bin.oper=`-`, trim="intersection")
      setnames(z, 2, "Rfr")
      setReflectorSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , Rfr := Rfr - e2]
      if (exists("Rpc", z, inherits=FALSE)) {
        z[ , Rpc := NULL]
      }
      return(z)
    } else {
      return(NA)
    }
  } else if (is(e1, "response.spct")) {
    if(is(e2, "response.spct")) {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$response, e2$response, bin.oper=`-`, trim="intersection")
      setnames(z, 2, "response")
      setResponseSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , response := response - e2]
      return(z)
    } else {
      return(NA)
    }
  } else if (is(e1, "source.spct")) {
    if (!exists("s.e.irrad", e1, inherits=FALSE)) {
      q2e(e2, "replace")
    }
    if (is(e2, "source.spct")) {
      if (!exists("s.e.irrad", e2, inherits=FALSE)) {
        q2e(e2, "replace")
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=`-`, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      return(z)
    } else if (is.numeric(e2)) {
      z <- copy(e1)
      z[ , s.e.irrad := s.e.irrad - e2]
      if (exists("s.q.irrad", z, inherits=FALSE)) {
        z[ , s.q.irrad := NULL]
      }
      return(z)
    }  else {
      return(NA)
    }
  }
}

# other operators  ---------------------------------------------------------------------



#' "^" operator for spectra
#'
#' Power operator for spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 a numeric vector. possibly of length one.
#' @export
#'
'^.generic.spct' <- function(e1, e2) {
  if(is(e1, "filter.spct") && is.numeric(e2)) {
    z <- copy(e1)
    z[ , Tfr := Tfr^e2]
    if (exists("Tpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    if (exists("A", z, inherits=FALSE)) {
      z[ , A := NULL]
    }
    return(z)
  } else if(is(e1, "reflector.spct") && is.numeric(e2)) {
    z <- copy(e1)
    z$Rfr <- z$Rfr^e2
    if (exists("Rpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    return(z)
  } else if(is(e1, "source.spct") && is.numeric(e2)) {
    z <- copy(e1)
    z[ , s.e.irrad := s.e.irrad^e2]
    if (exists("s.q.irrad", z, inherits=FALSE)) {
      z[ , s.q.irrad := NULL]
    }
    return(z)
  } else {
    return(NA)
  }
}

#' "log" function for spectra
#'
#' Logarirthm function for spectra.
#'
#' @param x an object of class "generic.spct"
#' @param base a positive number: the base with respect to which logarithms are computed. Defaults to e=exp(1).
#' @export
#'
'log.generic.spct' <- function(x, base = exp(1)) {
  if(is(x, "filter.spct")) {
    z <- copy(x)
    z$Tfr <- log(z$Tfr, base)
    if (exists("Tpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    if (exists("A", z, inherits=FALSE)) {
      z[ , A := NULL]
    }
    return(z)
  } else if(is(x, "reflector.spct")) {
    z <- copy(x)
    z$Rfr <- log(z$Rfr, base)
    if (exists("Rpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    z$Rpc <- z$Rfr * 100
    return(z)
  } else if(is(x, "source.spct")) {
    z <- copy(x)
    z$s.e.irrad <- log(z$s.e.irrad, base)
    if (exists("s.q.irrad", z, inherits=FALSE)) {
      z[ , s.q.irrad := NULL]
    }
    return(z)
  } else {
    return(NA)
  }
}

#' "log10" function for spectra
#'
#' Base 10 logairthm function for spectra.
#'
#' @param x an object of class "generic.spct"
#' @export
#'
'log10.generic.spct' <- function(x) {
  if(is(x, "filter.spct")) {
    z <- copy(x)
    z[ , Tfr := log10(Tfr)]
    if (exists("Tpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    if (exists("A", z, inherits=FALSE)) {
      z[ , A := NULL]
    }
    return(z)
  } else if(is(x, "reflector.spct")) {
    z <- copy(x)
    z[ , Rfr := log10(Rfr)]
    if (exists("Rpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    return(z)
  } else if(is(x, "source.spct")) {
    z <- copy(x)
    z[, s.e.irrad := log10(s.e.irrad)]
    if (exists("s.q.irrad", z, inherits=FALSE)) {
      z[ , s.q.irrad := NULL]
    }
    return(z)
  } else {
    return(NA)
  }
}

#' "sqrt" function for spectra
#'
#' Square root function for spectra.
#'
#' @param x an object of class "generic.spct"
#' @export
#'
'sqrt.generic.spct' <- function(x) {
  if(is(x, "filter.spct")) {
    z <- copy(x)
    z[ , Tfr := sqrt(Tfr)]
    if (exists("Tpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    if (exists("A", z, inherits=FALSE)) {
      z[ , A := NULL]
    }
    return(z)
  } else if(is(x, "reflector.spct")) {
    z <- copy(x)
    z[ , Rfr := sqrt(Rfr)]
    if (exists("Rpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    return(z)
  } else if(is(x, "source.spct")) {
    z <- copy(x)
    z[ , s.e.irrad := sqrt(s.e.irrad)]
    if (exists("s.q.irrad", z, inherits=FALSE)) {
      z[ , s.q.irrad := NULL]
    }
    return(z)
  } else {
    return(NA)
  }
}

#' "exp" function for spectra
#'
#' Exponential function for spectra.
#'
#' @param x an object of class "generic.spct"
#' @export
#'
'exp.generic.spct' <- function(x) {
  if(is(x, "filter.spct")) {
    z <- copy(x)
    z[ , Tfr := exp(Tfr)]
    if (exists("Tpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    if (exists("A", z, inherits=FALSE)) {
      z[ , A := NULL]
    }
    return(z)
  } else   if(is(x, "reflector.spct")) {
    z <- copy(x)
    z[ , Rfr <- exp(Rfr)]
    if (exists("Rpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    return(z)
  } else if(is(x, "source.spct")) {
    z <- copy(x)
    z[ , s.e.irrad := exp(s.e.irrad)]
    if (exists("s.q.irrad", z, inherits=FALSE)) {
      z[ , s.q.irrad := NULL]
    }
    return(z)
  } else {
    return(NA)
  }
}

#' "range" function for spectra
#'
#' Range function for spectra, returning wavelength range.
#'
#' @param ... not used in current version
#' @param na.rm a logical indicating whether missing values should be removed.
#' @export
#'
range.generic.spct <- function(..., na.rm=FALSE) {
  x <- c(...)
  return(range(x$w.length, na.rm=na.rm))
}

#' "max" function for spectra
#'
#' Maximun function for spectra, returning wavelength maximum.
#'
#' @param ... not used in current version
#' @param na.rm a logical indicating whether missing values should be removed.
#' @export
#'
max.generic.spct <- function(..., na.rm=FALSE) {
  x <- c(...)
  return(max(x$w.length, na.rm=na.rm))
}

#' "min" function for spectra
#'
#' Minimun function for spectra, returning wavelength minimum.
#'
#' @param ... not used in current version
#' @param na.rm a logical indicating whether missing values should be removed.
#' @export
#'
min.generic.spct <- function(..., na.rm=FALSE) {
  x <- c(...)
  return(min(x$w.length, na.rm=na.rm))
}

#' Labels of a "generic.spct" object.
#'
#' A function to obtain the labels of a spectrum. Currently returns 'names'.
#'
#' @param object an object of generic.spct
#' @param ... not used in current version
#'
#' @export
#'
labels.generic.spct <- function(object, ...) {
  return(names(object))
}


# transmittance and absorbance --------------------------------------------


# A2T ---------------------------------------------------------------------


#' Generic function
#'
#' Function that coverts absorbance into transmittance (fraction).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export A2T
A2T <- function(x, action, byref) UseMethod("A2T")

#' Default for generic function
#'
#' Function that coverts absorbance into transmittance (fraction).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export A2T.default
A2T.default <- function(x, action=NULL, byref=FALSE) {
  return(10^-x)
}

#' "generic.spct" function
#'
#' Function that coverts absorbance into transmittance (fraction).
#'
#' @param x a "filter.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export A2T.filter.spct
#'
A2T.filter.spct <- function(x, action="add", byref=FALSE) {
  if (byref) {
    z <- x
  } else {
    z <- copy(x)
  }
  if (exists("A", z, inherits=FALSE)) {
    z[ , Tfr := 10^-A]
    z[ , Tpc := Tfr * 100]
  } else {
    z[ , Tfr := NA]
    z[ , Tpc := NA]
  }
  if (action=="replace") {
    z[ , A := NULL]
  }
  return(z)
}


# T2A ---------------------------------------------------------------------


#' Generic function
#'
#' Function that coverts transmittance into absorbance (fraction).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export T2A
T2A <- function(x, action, byref) UseMethod("T2A")

#' Default for generic function
#'
#' Function that coverts transmittance into absorbance (fraction).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export T2A.default
T2A.default <- function(x, action=NULL, byref=FALSE) {
  return(-log10(x))
}

#' "filter.spct" function
#'
#' Function that coverts transmittance into absorbance (fraction).
#'
#' @param x a "filter.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export T2A.filter.spct
#'
T2A.filter.spct <- function(x, action="add", byref=FALSE) {
  if (is(x, "filter.spct")) {
    if (byref) {
      z <- x
    } else {
      z <- copy(x)
    }
    if (exists("Tfr", z, inherits=FALSE)) {
      z[ , A := -log10(Tfr)]
    } else {
      z[ , A := NA]
    }
    if (action=="replace") {
      z[ , Tfr := NULL]
      z[ , Tpc := NULL]
    }
    return(z)
  } else {
    return(NA)
  }
}


# energy - photon and photon - energy conversions -------------------------

# energy to photon ---------------------------------------------------------------------


#' Generic function
#'
#' Function that coverts spectral energy irradiance into spectral photon irradiance (molar).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export e2q
e2q <- function(x, action, byref) UseMethod("e2q")

#' Default for generic function
#'
#' Function that coverts spectral energy irradiance into spectral photon irradiance (molar).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export e2q.default
e2q.default <- function(x, action="add", byref=FALSE) {
  return(NA)
}

#' "source.spct" function
#'
#' Function that coverts spectral energy irradiance into spectral photon irradiance (molar).
#'
#' @param x a "source.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export e2q.source.spct
#'
e2q.source.spct <- function(x, action="add", byref=FALSE) {
  if (is(x, "source.spct")) {
    if (byref) {
      z <- x
    } else {
      z <- copy(x)
    }
    if (exists("s.q.irrad", z, inherits=FALSE)) {
      return(z)
    }
    if (exists("s.e.irrad", z, inherits=FALSE)) {
      z[ , s.q.irrad := s.e.irrad * e2qmol_multipliers(w.length)]
    } else {
      z[ , s.q.irrad := rep(NA, length(w.length))]
    }
    if (action=="replace") {
      z[ , s.e.irrad := NULL]
    }
    return(z)
  } else {
    return(NA)
  }
}

# photon to energy ---------------------------------------------------------------------


#' Generic function
#'
#' Function that coverts spectral photon irradiance (molar) into spectral energy irradiance.
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export q2e
q2e <- function(x, action, byref) UseMethod("q2e")

#' Default for generic function
#'
#' Function that coverts spectral photon irradiance (molar) into spectral energy irradiance.
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export q2e.default
q2e.default <- function(x, action="add", byref=FALSE) {
  return(NA)
}

#' "generic.spct" function
#'
#' Function that coverts spectral photon irradiance (molar) into spectral energy irradiance.
#'
#' @param x a "source.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export q2e.source.spct
#'
q2e.source.spct <- function(x, action="add", byref=FALSE) {
  if (is(x, "source.spct")) {
    if (byref) {
      z <- x
    } else {
      z <- copy(x)
    }
    if (exists("s.e.irrad", z, inherits=FALSE)) {
      return(z)
    }
    if (exists("s.e.irrad", z, inherits=FALSE)) {
      z[ , s.e.irrad := s.q.irrad / q2emol_multipliers(w.length)]
    } else {
      z[ , s.e.irrad := rep(NA, length(w.length))]
    }
    if (action=="replace") {
      z[ , s.q.irrad := NULL]
    }
    return(z)
  } else {
    return(NA)
  }
}
