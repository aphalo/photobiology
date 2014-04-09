#' Print a "waveband" object that can be used as imput when calculating irradiances.
#' 
#' @usage print(w_band)
#' 
#' @param w_band an object of class "waveband"
#' 
#' @export
#' 
print.waveband <- function(w_band) {
  cat(w_band$name, "\n")
  cat("low (nm)", round(w_band$low, 0), "\n")
  cat("high (nm)", round(w_band$high, 0), "\n")
  if (!is.null(w_band$weight)) cat("weighted", w_band$weight, "\n")
  if (!is.null(w_band$norm)) cat("normalized at", w_band$norm, "nm \n")
}