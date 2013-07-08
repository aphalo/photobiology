as_energy <- function(w.length, s.qmol.irrad){
  return(s.qmol.irrad * e2q_multipliers(w.length))
}
