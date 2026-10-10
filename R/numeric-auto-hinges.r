#' Guess whether insertion of hinges is needed or not
#'
#' Assess from wavelength vector whether insertion of hinges before further
#' calculations is needed.
#'
#' @param w.length numeric vector of wavelengths (nm).
#' @param step.limit numeric value (nm) for step, so that larger wavelength
#'    steps trigger insertion of hinges before computations on spectra.
#'
#' @return A logical value, TRUE if insertion of hinges is deemed necessary.
#'
#' @keywords internal
#'
auto_hinges <- function(w.length,
                        step.limit = 0.25) {
  # this uses average step.size because it is "cheap" to compute
  ((w.length[length(w.length)] - w.length[1]) / length(w.length)) >
    (step.limit * 0.75)
 }
