PG <- function(norm=300) {
  new_waveband(w.low=250,w.high=390,weight="SWF",SWF.fun=PG.q.fun,norm=norm)
}