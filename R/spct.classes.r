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
  if (!is.data.table(x)) {
      setDT(x)
  }
  setattr(x, "class", c("generic.spct", class(x)))
}

#' set class of a data.frame or data.table or generic.spct object to "filter.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "generic.spct" an object to store spectra
#' if the object is a data.frame is is mane a data.table
#'
#' @param x a data.frame
#' @export
#' @exportClass filter.spct
#'
setFilterSpct <- function(x) {
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is(x, "generic.spct")) {
    setGenSpct(x)
  }
  if (!is(x, "filter.spct")) {
    setattr(x, "class", c("filter.spct", class(x)))
  }
}

#' set class of a data.frame or data.table or generic.spct object to "source.spct"
#'
#' Sets the class attibute of a data.frame or data.table object to "generic.spct" an object to store spectra
#' if the object is a data.frame is is mane a data.table
#'
#' @param x a data.frame
#' @export
#' @exportClass source.spct
#'
setSourceSpct <- function(x) {
  if (!is.data.table(x)) {
    setDT(x)
  }
  if (!is(x, "generic.spct")) {
    setGenSpct(x)
  }
  if (!is(x, "sourcer_spct")) {
    setattr(x, "class", c("source.spct", class(x)))
  }
}

#' "*" operator for spectra
#'
#' Multiplication operator for spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 an object of class "generic.spct"
#' @name times-.generic.spct
#' @export
#'
'*.generic.spct' <- function(e1, e2) {
  if (is(e1, "filter.spct") && is(e2, "filter.spct")) {
    z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=`*`, trim="intersection")
    setDT(z)
    setnames(z, 2, "Tfr")
    z[ , Tpc := Tfr * 100]
    setFilterSpct(z)
    return(z)
  } else if(is(e1, "filter.spct") && is(e2, "source.spct")) {
    z1 <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.e.irrad, bin.oper=`*`, trim="intersection")
    z2 <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.q.irrad, bin.oper=`*`, trim="intersection")
    z <- data.table(w.length = z1$w.length, s.e.irrad = z1$s.irrad, s.q.irrad = z2$s.irrad)
    setSourceSpct(z)
    return(z)
  } else if(is(e2, "filter.spct") && is(e1, "source.spct")) {
    z1 <- oper_spectra(e2$w.length, e1$w.length, e2$Tfr, e1$s.e.irrad, bin.oper=`*`, trim="intersection")
    z2 <- oper_spectra(e2$w.length, e1$w.length, e2$Tfr, e1$s.q.irrad, bin.oper=`*`, trim="intersection")
    z <- data.table(w.length = z1$w.length, s.e.irrad = z1$s.irrad, s.q.irrad = z2$s.irrad)
    setSourceSpct(z)
    return(z)
  } else if(is(e1, "source.spct") && is(e2, "source.spct")) {
    z1 <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=`*`, trim="intersection")
    z2 <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.q.irrad, bin.oper=`*`, trim="intersection")
    z <- data.table(w.length = z1$w.length, s.e.irrad = z1$s.irrad, s.q.irrad = z2$s.irrad)
    setSourceSpct(z)
    return(z)
  } else if(is(e1, "filter.spct") && is.numeric(e2)) {
    z <- copy(e1)
    z$Tfr <- z$Tfr * e2
    z$Tpc <- z$Tfr * 100
    return(z)
  } else if(is(e1, "source.spct") && is.numeric(e2)) {
    z <- copy(e1)
    z$s.e.irrad <- z$s.e.irrad * e2
    z$s.q.irrad <- z$s.q.irrad * e2
    return(z)
  }else {
    return(NA)
  }
}

#' "/" operator for spectra
#'
#' Division operators for spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 an object of class "generic.spct"
#' @name slash-.generic.spct
#' @export
#'
'/.generic.spct' <- function(e1, e2) {
  if (is(e1, "filter.spct") && is(e2, "filter.spct")) {
    z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=`/`, trim="intersection")
    setDT(z)
    setnames(z, 2, "Tfr")
    z[ , Tpc := Tfr * 100]
    setFilterSpct(z)
    return(z)
  } else if(is(e1, "filter.spct") && is(e2, "source.spct")) {
    z1 <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.e.irrad, bin.oper=`/`, trim="intersection")
    z2 <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.q.irrad, bin.oper=`/`, trim="intersection")
    z <- data.table(w.length = z1$w.length, s.e.irrad = z1$s.irrad, s.q.irrad = z2$s.irrad)
    setSourceSpct(z)
    return(z)
  } else if(is(e2, "filter.spct") && is(e1, "source.spct")) {
    z1 <- oper_spectra(e2$w.length, e1$w.length, e2$Tfr, e1$s.e.irrad, bin.oper=`/`, trim="intersection")
    z2 <- oper_spectra(e2$w.length, e1$w.length, e2$Tfr, e1$s.q.irrad, bin.oper=`/`, trim="intersection")
    z <- data.table(w.length = z1$w.length, s.e.irrad = z1$s.irrad, s.q.irrad = z2$s.irrad)
    setSourceSpct(z)
    return(z)
  } else if(is(e1, "source.spct") && is(e2, "source.spct")) {
    z1 <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=`/`, trim="intersection")
    z2 <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.q.irrad, bin.oper=`/`, trim="intersection")
    z <- data.table(w.length = z1$w.length, s.e.irrad = z1$s.irrad, s.q.irrad = z2$s.irrad)
    setSourceSpct(z)
    return(z)
  } else if(is(e1, "filter.spct") && is.numeric(e2)) {
    z <- copy(e1)
    z$Tfr <- z$Tfr / e2
    z$Tpc <- z$Tfr * 100
    return(z)
  } else if(is(e1, "source.spct") && is.numeric(e2)) {
    z <- copy(e1)
    z$s.e.irrad <- z$s.e.irrad / e2
    z$s.q.irrad <- z$s.q.irrad / e2
    return(z)
  } else {
    return(NA)
  }
}

#' "+" operator for spectra
#'
#' Summation operator for spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 an object of class "generic.spct"
#' @name plus-.generic.spct
#' @export
#'
'+.generic.spct' <- function(e1, e2) {
  if (is(e1, "filter.spct") && is(e2, "filter.spct")) {
    z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=`+`, trim="intersection")
    setDT(z)
    setnames(z, 2, "Tfr")
    z[ , Tpc := Tfr * 100]
    setFilterSpct(z)
    return(z)
  } else if(is(e1, "filter.spct") && is(e2, "source.spct")) {
    z1 <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.e.irrad, bin.oper=`+`, trim="intersection")
    z2 <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.q.irrad, bin.oper=`+`, trim="intersection")
    z <- data.table(w.length = z1$w.length, s.e.irrad = z1$s.irrad, s.q.irrad = z2$s.irrad)
    setSourceSpct(z)
    return(z)
  } else if(is(e2, "filter.spct") && is(e1, "source.spct")) {
    z1 <- oper_spectra(e2$w.length, e1$w.length, e2$Tfr, e1$s.e.irrad, bin.oper=`+`, trim="intersection")
    z2 <- oper_spectra(e2$w.length, e1$w.length, e2$Tfr, e1$s.q.irrad, bin.oper=`+`, trim="intersection")
    z <- data.table(w.length = z1$w.length, s.e.irrad = z1$s.irrad, s.q.irrad = z2$s.irrad)
    setSourceSpct(z)
    return(z)
  } else if(is(e1, "source.spct") && is(e2, "source.spct")) {
    z1 <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=`+`, trim="intersection")
    z2 <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.q.irrad, bin.oper=`+`, trim="intersection")
    z <- data.table(w.length = z1$w.length, s.e.irrad = z1$s.irrad, s.q.irrad = z2$s.irrad)
    setSourceSpct(z)
    return(z)
  } else if(is(e1, "filter.spct") && is.numeric(e2)) {
    z <- copy(e1)
    z$Tfr <- z$Tfr + e2
    z$Tpc <- z$Tfr * 100
    return(z)
  } else if(is(e1, "source.spct") && is.numeric(e2)) {
    z <- copy(e1)
    z$s.e.irrad <- z$s.e.irrad + e2
    z$s.q.irrad <- NA
    return(z)
  } else {
    return(NA)
  }
}

#' "-" operator for spectra
#'
#' Unary negation and binary substraction operator for spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 an object of class "generic.spct"
#' @name minus-.generic.spct
#' @export
#'
'-.generic.spct' <- function(e1, e2=NULL) {
  # unary negation operator
  if (is.null(e2)) {
    if (is(e1, "filter.spct")) {
      z <- copy(e1)
      z[ , Tfr := -Tfr]
      z[ , Tpc := -Tpc]
      return(z)
    } else if(is(e1, "source.spct")) {
      z <- copy(e1)
      z[ , s.e.irrad := -s.e.irrad]
      z[ , s.q.irrad := -s.q.irrad]
      return(z)
    } else {
      return(NA)
    }
  }
  # binary substraction operator
  if (is(e1, "filter.spct") && is(e2, "filter.spct")) {
    z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=`-`, trim="intersection")
    setDT(z)
    setnames(z, 2, "Tfr")
    z[ , Tpc := Tfr * 100]
    setFilterSpct(z)
    return(z)
  } else if(is(e1, "filter.spct") && is(e2, "source.spct")) {
    z1 <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.e.irrad, bin.oper=`-`, trim="intersection")
    z2 <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.q.irrad, bin.oper=`-`, trim="intersection")
    z <- data.table(w.length = z1$w.length, s.e.irrad = z1$s.irrad, s.q.irrad = z2$s.irrad)
    setSourceSpct(z)
    return(z)
  } else if(is(e2, "filter.spct") && is(e1, "source.spct")) {
    z1 <- oper_spectra(e2$w.length, e1$w.length, e2$Tfr, e1$s.e.irrad, bin.oper=`-`, trim="intersection")
    z2 <- oper_spectra(e2$w.length, e1$w.length, e2$Tfr, e1$s.q.irrad, bin.oper=`-`, trim="intersection")
    z <- data.table(w.length = z1$w.length, s.e.irrad = z1$s.irrad, s.q.irrad = z2$s.irrad)
    setSourceSpct(z)
    return(z)
  } else if(is(e1, "source.spct") && is(e2, "source.spct")) {
    z1 <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=`-`, trim="intersection")
    z2 <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.q.irrad, bin.oper=`-`, trim="intersection")
    z <- data.table(w.length = z1$w.length, s.e.irrad = z1$s.irrad, s.q.irrad = z2$s.irrad)
    setSourceSpct(z)
    return(z)
  } else if(is(e1, "filter.spct") && is.numeric(e2)) {
    z <- copy(e1)
    z$Tfr <- z$Tfr - e2
    z$Tpc <- z$Tfr * 100
    return(z)
  } else if(is(e1, "source.spct") && is.numeric(e2)) {
    z <- copy(e1)
    z$s.e.irrad <- z$s.e.irrad - e2
    z$s.q.irrad <- NA
    return(z)
  } else {
    return(NA)
  }
}

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
    z$Tfr <- z$Tfr ^ e2
    z$Tpc <- z$Tfr * 100
    return(z)
  } else if(is(e1, "source.spct") && is.numeric(e2)) {
    z <- copy(e1)
    z$s.e.irrad <- z$s.e.irrad ^ e2
    z$s.q.irrad <- z$s.q.irrad ^ e2
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
    z$Tpc <- z$Tfr * 100
    return(z)
  } else if(is(x, "source.spct")) {
    z <- copy(x)
    z$s.e.irrad <- log(z$s.e.irrad, base)
    z$s.q.irrad <- log(z$s.q.irrad, base)
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
    z$Tfr <- log10(z$Tfr)
    z$Tpc <- z$Tfr * 100
    return(z)
  } else if(is(x, "source.spct")) {
    z <- copy(x)
    z$s.e.irrad <- log10(z$s.e.irrad)
    z$s.q.irrad <- log10(z$s.q.irrad)
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
    z$Tfr <- sqrt(z$Tfr)
    z$Tpc <- z$Tfr * 100
    return(z)
  } else if(is(x, "source.spct")) {
    z <- copy(x)
    z$s.e.irrad <- sqrt(z$s.e.irrad)
    z$s.q.irrad <- sqrt(z$s.q.irrad)
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
    z$Tfr <- exp(z$Tfr)
    z$Tpc <- z$Tfr * 100
    return(z)
  } else if(is(x, "source.spct")) {
    z <- copy(x)
    z$s.e.irrad <- exp(z$s.e.irrad)
    z$s.q.irrad <- exp(z$s.q.irrad)
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
#' A function to obtain the labels of a spectrum. Currently returns NA.
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

#' "gneric_spct" function
#'
#' Function that coverts absorbance into transmittance (fraction).
#'
#' @param x a "generic.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export A2T.generic.spct
#'
A2T.generic.spct <- function(x, action="add", byref=FALSE) {
  if (is(x, "filter.spct")) {
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
  } else {
    return(NA)
  }
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

#' "gneric_spct" function
#'
#' Function that coverts transmittance into absorbance (fraction).
#'
#' @param x a "generic.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export T2A.generic.spct
#'
T2A.generic.spct <- function(x, action="add", byref=FALSE) {
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

#' "generic.spct" function
#'
#' Function that coverts spectral energy irradiance into spectral photon irradiance (molar).
#'
#' @param x a "generic.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export e2q.generic.spct
#'
e2q.generic.spct <- function(x, action="add", byref=FALSE) {
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
#' @param x a "generic.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export q2e.generic.spct
#'
q2e.generic.spct <- function(x, action="add", byref=FALSE) {
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
