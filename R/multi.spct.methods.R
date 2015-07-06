#' Multi-sppct summary methods
#'
#' Functions
#'
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param f a function
#' @param idx logical whether to add a column with the names of the elements of mspct
#' @param ... other arguments passed to irrad.source.spct
#'
#' @export
#'
f_multi_spct <- function(mspct, f, ..., idx = !is.null(names(mspct))) {
  z <- NULL

  z <- lapply(mspct, f, ...)
  namesz <- names(z[[1]])
  z <- unlist(z, recursive = FALSE, use.names = FALSE)

  nspct <- length(mspct)
  nz <- length(z)
  nqty <- nz %/% nspct
  ncol <- attr(mspct, "ncol", exact = TRUE)
  nrow <- nspct / ncol

  stopifnot(nspct %% ncol == 0)

  rows <- rep(1:nrow, ncol)
  cols <- rep(1:ncol, rep(nrow, ncol))

  if (nz / nspct > 1) {
    z <- matrix(z, ncol = nqty, byrow = TRUE)
  }
  df <- data.frame(row = rows, col = cols, z = z)
  if (idx) {
    df[["idx"]] <- names(mspct)
  }

  qty.names <- paste(as.character(substitute(f)), gsub(" ", "", namesz), sep = "_")

  setnames(df, (1:nqty) + 2, qty.names)
  df
}

# generic_multi_spct methods -----------------------------------------------

#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#' @rdname  range.generic_spct
#'
range.generic_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, range, ..., idx = idx)
}

#' @param mspct an object of class generric_multi_spct or a derived class
#' @param idx logical whether to add a column with the names of the elements of mspct
#' @param ... other arguments passed to the underlaying \code{_spct} method
#'
#' @export
#' @rdname  min.generic_spct
#'
min.generic_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, min, ..., idx = idx)
}

#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#' @rdname  max.generic_spct
#'
max.generic_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, max, ..., idx = idx)
}

#' @describeIn stepsize  Method for "generic_multi_spct" objects for generic function.
#'
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
stepsize.generic_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, stepsize, ..., idx = idx)
}

#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#' @rdname  labels.generic_spct
#'
labels.generic_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, labels, ..., idx = idx)
}

# source_multi_spct methods -----------------------------------------------

#' @describeIn irrad  Calculates irradiance from a \code{source_multi_spct}
#'   object.
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
irrad.source_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, irrad, ..., idx = idx)
}

#' @describeIn q_irrad  Calculates photon (quantum) irradiance from a \code{source_multi_spct}
#'   object.
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
q_irrad.source_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, q_irrad, ..., idx = idx)
}

#' @describeIn e_irrad  Calculates energy irradiance from a \code{source_multi_spct}
#'   object.
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
e_irrad.source_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, e_irrad, ..., idx = idx)
}

#' @describeIn fluence Calculates fluence from a \code{source_multi_spct}
#'   object.
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
fluence.source_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, fluence, ..., idx = idx)
}

#' @describeIn e_fluence Calculates energy fluence from a \code{source_multi_spct}
#'   object.
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
e_fluence.source_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, e_fluence, ..., idx = idx)
}

#' @describeIn q_fluence Calculates photon (quantum) fluence from a \code{source_multi_spct}
#'   object.
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
q_fluence.source_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, q_fluence, ..., idx = idx)
}

#' @describeIn q_ratio Calculates photon:photon from a \code{source_multi_spct}
#'   object.
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
q_ratio.source_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, q_ratio, ..., idx = idx)
}

#' @describeIn e_ratio Calculates energy:energy ratio from a \code{source_multi_spct}
#'   object.
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
e_ratio.source_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, e_ratio, ..., idx = idx)
}

#' @describeIn eq_ratio Calculates energy:photon from a \code{source_multi_spct}
#'   object.
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
eq_ratio.source_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, eq_ratio, ..., idx = idx)
}

#' @describeIn qe_ratio Calculates photon:energy ratio from a \code{source_multi_spct}
#'   object.
#' @export
#'
qe_ratio.source_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, qe_ratio, ..., idx = idx)
}

#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#' @rdname color.source_spct
#'
color.source_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, color, ..., idx = idx)
}

# filter_multi_spct methods -----------------------------------------------

#' @describeIn transmittance Calculates transmittance from a \code{filter_multi_spct}
#'
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
transmittance.filter_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, transmittance, ..., idx = idx)
}

#' @describeIn absorptance Calculates absorptance from a \code{filter_multi_spct}
#'
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
absorptance.filter_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, absorptance, ..., idx = idx)
}

#' @describeIn absorbance Calculates absorbance from a \code{filter_multi_spct}
#'
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
absorbance.filter_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, absorbance, ..., idx = idx)
}

# reflector_multi_spct methods -----------------------------------------------

#' @describeIn reflectance Calculates reflectance from a \code{reflector_multi_spct}
#'
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
reflectance.reflector_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, reflectance, ..., idx = idx)
}

# object_multi_spct methods -----------------------------------------------

#' @describeIn transmittance Calculates transmittance from a \code{object_multi_spct}
#'
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
transmittance.object_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, transmittance, ..., idx = idx)
}

#' @describeIn absorptance Calculates absorptance from a \code{object_multi_spct}
#'
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
absorptance.object_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, absorptance, ..., idx = idx)
}

#' @describeIn reflectance Calculates reflectance from a \code{object_multi_spct}
#'
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
reflectance.object_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, reflectance, ..., idx = idx)
}

#' @describeIn absorbance Calculates absorbance from a \code{object_multi_spct}
#'
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
absorbance.object_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, absorbance, ..., idx = idx)
}

# response_multi_spct methods -----------------------------------------------

#' @describeIn response Calculates response from a \code{response_multi_spct}
#'
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
response.response_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, response, ..., idx = idx)
}

#' @describeIn q_response Calculates photon (quantum) response from a \code{response_multi_spct}
#'
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
q_response.response_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, q_response, ..., idx = idx)
}

#' @describeIn e_response Calculates energy response from a \code{response_multi_spct}
#'
#' @param mspct an object of class generric_multi_spct or a derived class
#' @param ... other arguments passed to the underlaying \code{_spct} method
#' @param idx logical whether to add a column with the names of the elements of mspct
#'
#' @export
#'
e_response.response_multi_spct <- function(mspct, ..., idx = !is.null(names(mspct))) {
  f_multi_spct(mspct, e_response, ..., idx = idx)
}

