


oper.generic.spct <- function(e1, e2, oper) {
  class1 <- class.spct(e1)[1]
  class2 <- class.spct(e2)[1]
  if (class1 == "source.spct") {
    q2e(e1, action = "add")
    if (is.numeric(e2)) {
      invisible(source.spct(w.length=e1$w.length, s.e.irrad=oper(e1$s.e.irrad, e2)))
    } else if (class2 == "source.spct") {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      invisible(z)
    } else if (class2 == "filter.spct") {
      if (!identical(oper, `*`)) invisible(NA)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$Tfr, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      invisible(z)
    } else if (class2 == "reflector.spct") {
      if (!identical(oper, `*`)) invisible(NA)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      invisible(z)
    } else if (class2 == "response.spct") {
      q2e(e2, action = "add")
      if (!identical(oper, `*`)) invisible(NA)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.response, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.e.response")
      setResponseSpct(z)
      invisible(z)
    } else if (class2 == "chroma.spct") {
      if (!identical(oper, `*`)) invisible(NA)
      x <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$x, bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$y, bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$z, bin.oper=oper, trim="intersection")
      out.spct <- data.table(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]])
      setChromaSpct(out.spct)
      invisible(out.spct)
    } else { # this traps also e2 == "generic.spct"
      invisible(NA)
    }
  } else if (class1 == "filter.spct") {
    A2T(e1)
    if (is.numeric(e2)) {
      invisible(filter.spct(w.length=e1$w.length, Tfr=oper(e1$Tfr, e2)))
    } else if (class2 == "source.spct") {
      if (!identical(oper, `*`)) invisible(NA)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      invisible(z)
    } else if (class2 == "filter.spct") {
      A2T(e2)
      if (!identical(oper, `*`)) invisible(NA)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=oper, trim="intersection")
      setnames(z, 2, "Tfr")
      setSourceSpct(z)
      invisible(z)
    } else { # this traps optically illegal operations
      invisible(NA)
    }
  } else if (class1 == "reflector.spct") {
    A2T(e1)
    if (is.numeric(e2)) {
      invisible(reflector.spct(w.length=e1$w.length, Rfr=oper(e1$Rfr, e2)))
    } else if (class2 == "source.spct") {
      if (!identical(oper, `*`)) invisible(NA)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z)
      invisible(z)
    } else { # this traps optically illegal operations
      invisible(NA)
    }
  } else if (class1 == "response.spct") {
    q2e(e1)
    if (is.numeric(e2)) {
      invisible(response.spct(w.length=e1$w.length, s.e.response=oper(e1$s.e.response, e2)))
    } else if (class2 == "source.spct") {
      if (!identical(oper, `*`)) invisible(NA)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.response, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.e.response")
      setResponseSpct(z)
      invisible(z)
    } else { # this traps optically illegal operations
      invisible(NA)
    }
  } else if (class1 == "chroma.spct") {
    if (is.numeric(e2)) {
      e3 <- copy(e1)
      if (length(e2) == 3 && names(e2) == c("x", "y", "z")) {
        e3[ , `:=`(x = x * e2["x"] , y = y * e2["y"], z = z  * e2["z"])]
        invisible(e3)
      } else {
        e3[ , `:=`(x = x * e2 , y = y * e2, z = z  * e2)]
        invisible(e3)
      }
    } else if (class2 == "source.spct") {
        x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$s.e.irrad, bin.oper=oper, trim="intersection")
        y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$s.e.irrad, bin.oper=oper, trim="intersection")
        z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$s.e.irrad, bin.oper=oper, trim="intersection")
        invisible(chroma.spct(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]]))
      } else {
        invisible(NA)
      }
  } else {
    invisible(NA)
  }
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
  invisible(oper.generic.spct(e1, e2, `*`))
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
  invisible(oper.generic.spct(e1, e2, `/`))
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
  invisible(oper.generic.spct(e1, e2, `+`))
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
    invisible(oper.generic.spct(e1, e2, `-`))
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
  invisible(oper.generic.spct(e1, e2, `^`))
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
    A2T(x)
    z <- copy(x)
    z$Tfr <- log(z$Tfr, base)
    if (exists("Tpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    if (exists("A", z, inherits=FALSE)) {
      z[ , A := NULL]
    }
    invisible(z)
  } else if(is(x, "reflector.spct")) {
    z <- copy(x)
    z$Rfr <- log(z$Rfr, base)
    if (exists("Rpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    z$Rpc <- z$Rfr * 100
    invisible(z)
  } else if(is(x, "source.spct")) {
    q2e(x)
    z <- copy(x)
    z$s.e.irrad <- log(z$s.e.irrad, base)
    if (exists("s.q.irrad", z, inherits=FALSE)) {
      z[ , s.q.irrad := NULL]
    }
    invisible(z)
  } else if(is(x, "response.spct")) {
    z <- copy(x)
    z$s.e.response <- log(z$s.e.response, base)
    if (exists("s.q.irrad", z, inherits=FALSE)) {
      z[ , s.q.irrad := NULL]
    }
    invisible(z)
  }else {
    invisible(NA)
  }
}

#' "log10" function for spectra
#'
#' Base 10 logarithm function for spectra.
#'
#' @param x an object of class "generic.spct"
#' @export
#'
'log10.generic.spct' <- function(x) {
  invisible(log.generic.spct(x, base = 10))
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
    A2T(x)
    z <- copy(x)
    z$Tfr <- sqrt(z$Tfr)
    if (exists("Tpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    if (exists("A", z, inherits=FALSE)) {
      z[ , A := NULL]
    }
    invisible(z)
  } else if (is(x, "reflector.spct")) {
    z <- copy(x)
    z$Rfr <- sqrt(z$Rfr)
    if (exists("Rpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
      z$Rpc <- z$Rfr * 100
    }
    invisible(z)
  } else if(is(x, "source.spct")) {
    q2e(x)
    z <- copy(x)
    z$s.e.irrad <- sqrt(z$s.e.irrad)
    if (exists("s.q.irrad", z, inherits=FALSE)) {
      z[ , s.q.irrad := NULL]
    }
    invisible(z)
  } else if(is(x, "response.spct")) {
    z <- copy(x)
    z$s.e.response <- sqrt(z$s.e.response)
    if (exists("s.q.irrad", z, inherits=FALSE)) {
      z[ , s.q.irrad := NULL]
    }
    invisible(z)
  }else {
    invisible(NA)
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
    A2T(x)
    z <- copy(x)
    z$Tfr <- exp(z$Tfr)
    if (exists("Tpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
    }
    if (exists("A", z, inherits=FALSE)) {
      z[ , A := NULL]
    }
    invisible(z)
  } else if (is(x, "reflector.spct")) {
    z <- copy(x)
    z$Rfr <- exp(z$Rfr)
    if (exists("Rpc", z, inherits=FALSE)) {
      z[ , Tpc := NULL]
      z$Rpc <- z$Rfr * 100
    }
    invisible(z)
  } else if(is(x, "source.spct")) {
    q2e(x)
    z <- copy(x)
    z$s.e.irrad <- exp(z$s.e.irrad)
    if (exists("s.q.irrad", z, inherits=FALSE)) {
      z[ , s.q.irrad := NULL]
    }
    invisible(z)
  } else if(is(x, "response.spct")) {
    z <- copy(x)
    z$s.e.response <- exp(z$s.e.response)
    if (exists("s.q.irrad", z, inherits=FALSE)) {
      z[ , s.q.irrad := NULL]
    }
    invisible(z)
  }else {
    invisible(NA)
  }
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
    name <- substitute(x)
  } else {
    x <- copy(x)
  }
  if (exists("Tfr", x, inherits=FALSE)) {
    NULL
  } else if (exists("A", x, inherits=FALSE)) {
    x[ , Tfr := 10^-A]
  } else {
    x[ , Tfr := NA]
  }
  if (action=="replace" && exists("A", x, inherits=FALSE)) {
    x[ , A := NULL]
  }
  if (byref && is.name(name)) { # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
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
  if (byref) {
    name <- substitute(x)
  } else {
    x <- copy(x)
  }
  if (exists("A", x, inherits=FALSE)) {
    NULL
  } else if (exists("Tfr", x, inherits=FALSE)) {
    x[ , A := -log10(Tfr)]
  } else {
    x[ , A := NA]
  }
  if (action=="replace" && exists("Tfr", x, inherits=FALSE)) {
    x[ , Tfr := NULL]
  }
  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
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
  if (byref) {
    name <- substitute(x)
  } else {
    x <- copy(x)
  }
  if (exists("s.q.irrad", x, inherits=FALSE)) {
    NULL
  } else if (exists("s.e.irrad", x, inherits=FALSE)) {
    x[ , s.q.irrad := s.e.irrad * e2qmol_multipliers(w.length)]
  } else {
    x[ , s.q.irrad := NA]
  }
  if (action=="replace" && exists("s.e.irrad", x, inherits=FALSE)) {
    x[ , s.e.irrad := NULL]
  }
  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' "response.spct" function
#'
#' Function that coverts response to spectral energy irradiance into response to spectral photon irradiance (molar).
#'
#' @param x a "response.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export e2q.response.spct
#'
e2q.response.spct <- function(x, action="add", byref=FALSE) {
  if (byref) {
    name <- substitute(x)
  } else {
    x <- copy(x)
  }
  if (exists("s.q.response", x, inherits=FALSE)) {
    NULL
  } else if (exists("s.e.response", x, inherits=FALSE)) {
    x[ , s.q.response := s.e.response / e2qmol_multipliers(w.length)]
  } else {
    x[ , s.q.response := NA]
  }
  if (action=="replace" && exists("s.e.response", x, inherits=FALSE)) {
    x[ , s.e.response := NULL]
  }
  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
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

#' "source.spct" function
#'
#' Function that coverts spectral photon irradiance (molar) into spectral energy irradiance.
#'
#' @param x a "source.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export q2e.source.spct
#'
q2e.source.spct <- function(x, action="add", byref=FALSE) {
  if (byref) {
    name <- substitute(x)
  } else {
    x <- copy(x)
  }
  if (exists("s.e.irrad", x, inherits=FALSE)) {
    NULL
  } else if (exists("s.q.irrad", x, inherits=FALSE)) {
    x[ , s.e.irrad := s.q.irrad / e2qmol_multipliers(w.length)]
  } else {
    x[ , s.e.irrad := NA]
  }
  if (action=="replace" && exists("s.q.irrad", x, inherits=FALSE)) {
    x[ , s.q.irrad := NULL]
  }
  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' "response.spct" function
#'
#' Function that coverts response to spectral photon irradiance (molar) into response to spectral energy irradiance.
#'
#' @param x a "response.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export q2e.response.spct
#'
q2e.response.spct <- function(x, action="add", byref=FALSE) {
  if (byref) {
    name <- substitute(x)
  } else {
    x <- copy(x)
  }
  if (exists("s.e.response", x, inherits=FALSE)) {
    NULL
  } else if (exists("s.q.response", x, inherits=FALSE)) {
    x[ , s.e.response := s.q.response * e2qmol_multipliers(w.length)]
  } else {
    x[ , s.e.response := NA]
  }
  if (action=="replace" && exists("s.q.response", x, inherits=FALSE)) {
    x[ , s.q.irrad := NULL]
  }
  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}
