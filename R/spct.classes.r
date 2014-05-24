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
  if (!is(x, "sourcer.spct")) {
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

#' "range" function for spectra
#' 
#' Range function for spectra, returning wavelength range.
#' 
#' @param x an object of class "generic.spct"
#' @param ... not used in current version
#' @param na.rm a logical indicating whether missing values should be removed.
#' @export
#'
'range.generic.spct' <- function(x, ..., na.rm) {
  if(is(x, "filter.spct") || is(x, "source.spct")) {
    return(range(x$w.length, na.rm))
  } else {                                               
    return(NA)
  }
}  

#' "max" function for spectra
#' 
#' Maximun function for spectra, returning wavelength maximum.
#' 
#' @param x an object of class "generic.spct"
#' @param ... not used in current version
#' @param na.rm a logical indicating whether missing values should be removed.
#' @export
#'
'max.generic.spct' <- function(x, ..., na.rm) {
  if(is(x, "filter.spct") || is(x, "source.spct")) {
    return(max(x$w.length, na.rm))
  } else {                                               
    return(NA)
  }
}  

#' "min" function for spectra
#' 
#' Minimun function for spectra, returning wavelength minimum.
#' 
#' @param x an object of class "generic.spct"
#' @param ... not used in current version
#' @param na.rm a logical indicating whether missing values should be removed.
#' @export
#'
min.generic.spct <- function(x, ..., na.rm) {
  if(is(x, "filter.spct") || is(x, "source.spct")) {
    return(min(x$w.length, na.rm))
  } else {                                               
    return(NA)
  }
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
  return(NA)
}

