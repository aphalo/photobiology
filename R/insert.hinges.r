#' Insert new wavelengths into the data, interpolating the corresponding spectral irradiance values.
#'
#' Inserting wavelengths values immediately bafore and after a discontinuity in the SWF,
#' gretaly reduces the errors caused by interpolating the weighted irradiance during
#' integration of the effective spectral irradiance. This is specially true when data
#' has a large wavelength step size.
#'  
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral irradiance values
#' @param hinges a numeric array giving the wavelengths at which the s.irrad should be inserted by
#' interpolation, no interpolation is indicated by an empty array (numeric(0))
#' 
#' @return a data.frame with variables \code{w.length} and \code{s.e.irrad}
#' 
#' @keywords manip misc
#' @export
#' @examples
#' data(sun.data)
#' with(sun.data, insert_hinges(w.length, s.e.irrad, c(399.99,400.00,699.99,700.00))

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
