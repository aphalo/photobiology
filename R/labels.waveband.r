#' Name of a "waveband" object.
#' 
#' A function to obtain the name of objects of class "waveband".
#' 
#' @param object an object of class "waveband"
#' @param ... not used in current version
#' 
#' @export
#' 
labels.waveband <- function(object, ...) {
  return(list(label = object$label, name = object$name))
}
