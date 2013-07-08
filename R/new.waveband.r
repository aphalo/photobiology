new_waveband <- function(w.low, w.high, weight=NULL, SWF.fun=NULL, norm=NULL, h, hinges=NULL){
  return(list(low=w.low, high=w.high, weight=weight, SWF.fun=SWF.fun, norm=norm, hinges=hinges))
} 