
# using options -----------------------------------------------------------

#' Use photobiology options
#'
#' Execute an R expression, possibly compound, using a certain setting for
#' spectral data related options.
#'
#' @param expr an R expression to execute.
#'
#' @return The value returned by the execution of \code{expression}.
#'
#' @references Based on \code{withOptions()} as offered by Thomas Lumley, and
#'   listed in \url{https://www.burns-stat.com/the-options-mechanism-in-r/},
#'   section Deep End, of "The Options mechanism in R" by Patrick Burns.
#'
#' @export
#'
using_Tfr <- function(expr) {
  old <- options(photobiology.filter.qty = "transmittance")
  on.exit(options(old))
  expr <- substitute(expr)
  eval.parent(expr)
}

#' @rdname using_Tfr
#'
#' @export
#'
using_Afr <- function(expr) {
  old <- options(photobiology.filter.qty = "absorptance")
  on.exit(options(old))
  expr <- substitute(expr)
  eval.parent(expr)
}

#' @rdname using_Tfr
#'
#' @export
#'
using_A <- function(expr) {
  old <- options(photobiology.filter.qty = "absorbance")
  on.exit(options(old))
  expr <- substitute(expr)
  eval.parent(expr)
}

#' @rdname using_Tfr
#'
#' @export
#'
using_energy <- function(expr) {
  old <- options(photobiology.radiation.unit = "energy")
  on.exit(options(old))
  expr <- substitute(expr)
  eval.parent(expr)
}

#' @rdname using_Tfr
#'
#' @export
#'
using_photon <- function(expr) {
  old <- options(photobiology.radiation.unit = "photon")
  on.exit(options(old))
  expr <- substitute(expr)
  eval.parent(expr)
}

#' @rdname using_Tfr
#'
#' @export
#'
using_quantum <- using_photon

# Set options -----------------------------------------------------------

#' Set spectral-data options
#'
#' Set spectral-data related options easily.
#'
#' @return Previous value of the modified option.
#'
#' @export
#'
energy_as_default <- function() {
  options(photobiology.radiation.unit = "energy")
}

#' @rdname energy_as_default
#'
#' @export
#'
photon_as_default <- function() {
  options(photobiology.radiation.unit = "photon")
}

#' @rdname energy_as_default
#'
#' @export
#'
quantum_as_default <- photon_as_default

#' @rdname energy_as_default
#'
#' @export
#'
Tfr_as_default <- function() {
  options(photobiology.filter.qty = "transmittance")
}

#' @rdname energy_as_default
#'
#' @export
#'
Afr_as_default <- function() {
  options(photobiology.filter.qty = "absorptance")
}

#' @rdname energy_as_default
#'
#' @export
#'
A_as_default <- function() {
  options(photobiology.filter.qty = "absorbance")
}

#' Set computation options
#'
#' Set computation related options easily.
#'
#' @param flag logical.
#'
#' @return Previous value of the modified option.
#'
#' @export
#'
wb_trim_as_default <- function(flag = TRUE) {
  options(photobiology.waveband.trim = flag)
}

#' @rdname wb_trim_as_default
#'
#' @export
#'
use_cached_mult_as_default <- function(flag = TRUE) {
  options(photobiology.use.cached.mult = flag)
}

#' @rdname energy_as_default
#'
#' @export
#'
unset_radiation_unit_default <- function() {
  options(photobiology.radiation.unit = NULL)
}

#' @rdname energy_as_default
#'
#' @export
#'
unset_filter_qty_default <- function() {
  options(photobiology.filter.qty = NULL)
}

#' @rdname energy_as_default
#'
#' @export
#'
unset_user_defaults <- function() {
  options(photobiology.filter.qty = NULL,
          photobiology.radiation.unit = NULL,
          photobiology.verbose = getOption("verbose"),
          photobiology.strict.range = NULL,
          photobiology.waveband.trim = NULL,
          photobiology.use.cached.mult = NULL)
}

#' Set error reporting options
#'
#' Set error reporting related options easily.
#'
#' @param flag logical.
#'
#' @return Previous value of the modified option.
#'
#' @export
#'
verbose_as_default <- function(flag = TRUE) {
  if (is.null(flag)) {
    flag <- getOption("verbose")
  }
  options(photobiology.verbose = flag)
}

#' @rdname verbose_as_default
#'
#' @export
#'
strict_range_as_default <- function(flag = TRUE) {
  options(photobiology.strict.range = flag)
}
