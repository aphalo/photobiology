insert_hinges <- function(w.length, s.irrad, hinges)
{
  # discard repeated values
  hinges <- unique(sort(hinges))
  
  # we will use a logical vector to select those hinges that we will actually use
  hinges.selector <- rep(FALSE, length(hinges))
  
  # select those hinges within the data range 
  hinges.selector[which(hinges > min(w.length) & (hinges < max(w.length)))] <- TRUE
  
  # deselect hinge wavelengths already included in the data
  for (i in 1:length(w.length))
  {
    hinges.selector[which(hinges==w.length[i])] <- FALSE  
  }
  
  # now we do the actual selection of hinges to insert
  hinges.used <- hinges[hinges.selector]
  
  # now we insert them one by one
  for (h in hinges.used)
  {
    i <- which(w.length > h)[1]
    # spectral irradiance 
    # (beware that the indexes used below are dependent on the order of statements)
    # (we need to interpolate the value of the inserted value between adjacent ones)
    delta.irrad <- (s.irrad[i] - s.irrad[i-1]) / 
      (w.length[i] - w.length[i-1]) * (h - w.length[i-1])
    irrad <- s.irrad[i-1] + delta.irrad
    s.irrad <- append(s.irrad, irrad, i-1)
    # w.length
    w.length <- append(w.length, h, i-1)
  }
  return(data.frame(w.length=w.length, s.irrad=s.irrad))
}
