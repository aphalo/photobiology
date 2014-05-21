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
#' Sets the class attibute of a data.frame or data.table object to "generic.spct" an object to store spectra
#' if the object is a data.frame is is mane a data.table
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
    z$s.q.irrad <- NA
    return(z)
  }else {                                               
    return(NA)
  }
}  

#' "/" operator for spectra
#' 
#' Sets the class attibute of a data.frame or data.table object to "generic.spct" an object to store spectra
#' if the object is a data.frame is is mane a data.table
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
    z$s.q.irrad <- NA
    return(z)
  } else {                                               
    return(NA)
  }
}  

#' "+" operator for spectra
#' 
#' Sets the class attibute of a data.frame or data.table object to "generic.spct" an object to store spectra
#' if the object is a data.frame is is mane a data.table
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
#' Sets the class attibute of a data.frame or data.table object to "generic.spct" an object to store spectra
#' if the object is a data.frame is is mane a data.table
#' 
#' @param e1 an object of class "generic.spct"
#' @param e2 an object of class "generic.spct"
#' @name minus-.generic.spct
#' @export
#'
'-.generic.spct' <- function(e1, e2) {
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

