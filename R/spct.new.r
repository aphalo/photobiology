#' Create a new constant value source.spct
#'
#' @export
new.flat.source.spct <- function(w.length, irrad.unit = "energy", irrad.value = 1, time.unit="second") {
  if (irrad.unit == "energy") {
    out.spct <- data.table(w.length=w.length, s.e.irrad=irrad.value)
  } else if (irrad.unit == "quantum" || irrad.unit == "photon") {
    out.spct <- data.table(w.length=w.length, s.q.irrad=irrad.value)
  } else {
    out.spct <- data.table(w.length=w.length, s.e.irrad=NA)
  }
  setSourceSpct(out.spct)
  setattr(out.spct, "time.unit", time.unit)
  return(out.spct)
}

#' Create a new constant value filter.spct
#'
#' @export
new.flat.filter.spct <- function(w.length, Tfr.value = 1, T.type=c("total","internal")) {
  out.spct <- data.table(w.length=w.length, Tfr=Tfr.value)
  setFilterSpct(out.spct)
  setattr(out.spct, "T.type", T.type)
  return(out.spct)
}

#' Create a new constant value reflector.spct
#'
#' @export
new.flat.reflector.spct <- function(w.length, Rfr.value = 1, R.type=c("total","specular")) {
  out.spct <- data.table(w.length=w.length, Rfr=Rfr.value)
  setReflectorSpct(out.spct)
  setattr(out.spct, "R.type", R.type)
  return(out.spct)
}

