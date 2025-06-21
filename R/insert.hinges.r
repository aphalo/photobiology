#' Insert wavelength values into spectral data.
#'
#' Inserting wavelengths values immediately before and after a discontinuity in
#' the SWF, greatly reduces the errors caused by interpolating the weighted
#' irradiance during integration of the effective spectral irradiance. This is
#' specially true when data have a large wavelength step size.
#'
#' @param x numeric vector (sorted in increasing order)
#' @param y numeric vector
#' @param h a numeric vector giving the wavelengths at which the y values
#'   should be inserted by interpolation, no interpolation is indicated by an
#'   empty vector (numeric(0))
#'
#' @return a data.frame with variables \code{x} and \code{y}. Unless the hinge
#'   values were already present in \code{y}, each inserted hinge, expands the
#'   vectors returned in the data frame by one value.
#'
#' @export
#' @examples
#' with(sun.data,
#'     insert_hinges(w.length, s.e.irrad,
#'        c(399.99, 400.00, 699.99, 700.00)))
#'
#' @note Insertion is a costly operation but I have tried to optimize this
#' function as much as possible by avoiding loops. Earlier this function was
#' implemented in C++, but a bug was discovered and I have now rewritten it
#' using R.
#'
#' @family low-level functions operating on numeric vectors.
#'
insert_hinges <- function(x, y, h) {
  # sanitize 'hinges'
  h <- setdiff(h, x)
  h <- h[h > x[1] & h < x[length(x)]]
  h <- sort(h)
  # save lengths
  j <- length(h)
  k <- length(x)
  # compute the new idx values after shift caused by insertions
  idxs.in <- findInterval(h, x)
  idxs.out <- idxs.in + 1:j
  idxs.diff <- diff(c(0, idxs.in, k))
  idxs.map <- 1:k + rep(0:j, idxs.diff)
  # we use recycling to interpolate all values and insert them into the gaps
  if (is.numeric(y)) {
    # allocate vectors for the output
    x.out <- numeric(j + k)
    y.out <- numeric(j + k)
    # we copy everything that does not require interpolation
    x.out[idxs.map] <- x
    y.out[idxs.map] <- y
    x.out[idxs.out] <- h
    # we fill the hinges by interpolation
    y.out[idxs.out] <- y[idxs.in + 1] -
      (x[idxs.in + 1] - x.out[idxs.out]) /
      (x[idxs.in + 1] - x[idxs.in]) *
      (y[idxs.in + 1] - y[idxs.in])
  } else if (length(unique(y)) == 1) {
    y.out <- rep(y[1], j + k)
  } else {
    y.out <- rep(y[NA], j + k)
    y.out[idxs.map] <- y
  }
  data.frame(x = x.out, y = y.out)
}

#' Insert wavelength values into spectral data.
#'
#' Inserting wavelengths values immediately before and after a discontinuity in
#' the SWF, greatly reduces the errors caused by interpolating the weighted
#' irradiance during integration of the effective spectral irradiance. This is
#' specially true when data have a large wavelength step size.
#'
#' @param x numeric vector (sorted in increasing order)
#' @param y numeric vector
#' @param h a numeric vector giving the wavelengths at which the y values
#'   should be inserted by interpolation, no interpolation is indicated by an
#'   empty vector (numeric(0))
#'
#' @return a list with variables \code{x} and \code{y}. Unless the hinge values
#'   were already present in \code{y}, each inserted hinge, expands the vectors
#'   by two values.
#'
#' @keywords internal
#'
#' @family low-level functions operating on numeric vectors.
#'
l_insert_hinges <- function(x, y, h) {
  # sanitize 'hinges'
  h <- setdiff(h, x)
  h <- h[h > x[1] & h < x[length(x)]]
  h <- sort(h)
  # save lengths
  j <- length(h)
  k <- length(x)
  # compute the new idx values after shift caused by insertions
  idxs.in <- findInterval(h, x)
  idxs.out <- idxs.in + 1:j
  idxs.diff <- diff(c(0, idxs.in, k))
  idxs.map <- 1:k + rep(0:j, idxs.diff)
  # we use recycling to interpolate all values and insert them into the gaps
  if (is.numeric(y)) {
    # allocate vectors for the output
    x.out <- numeric(j + k)
    y.out <- numeric(j + k)
    # we copy everything that does not require interpolation
    x.out[idxs.map] <- x
    y.out[idxs.map] <- y
    x.out[idxs.out] <- h
    # we fill the hinges by interpolation
    y.out[idxs.out] <- y[idxs.in + 1] -
      (x[idxs.in + 1] - x.out[idxs.out]) /
      (x[idxs.in + 1] - x[idxs.in]) *
      (y[idxs.in + 1] - y[idxs.in])
  } else if (length(unique(y)) == 1) {
    y.out <- rep(y[1], j + k)
  } else {
    y.out <- rep(y[NA], j + k)
    y.out[idxs.map] <- y
  }
  list(x = x.out, y = y.out)
}

#' Insert spectral data values at new wavelength values.
#'
#' Inserting wavelengths values immediately before and after a discontinuity in
#' the SWF, greatly reduces the errors caused by interpolating the weighted
#' irradiance during integration of the effective spectral irradiance. This is
#' specially true when data have a relatively large wavelength step size and/or
#' when the weighting function used has discontinuities in its value or slope.
#' This function differs from \code{insert_hinges()} in that it returns a vector
#' of \code{y} values instead of a \code{tibble}.
#'
#' @param x numeric vector (sorted in increasing order).
#' @param y numeric vector.
#' @param h a numeric vector giving the wavelengths at which the y values should
#'   be inserted by interpolation, no interpolation is indicated by an empty
#'   numeric vector (\code{numeric(0)}).
#'
#' @return A numeric vector with the numeric values of \code{y}, but longer.
#'   Unless the hinge values were already present in \code{y}, each inserted
#'   hinge, expands the vector by two values.
#'
#' @keywords internal.
#'
#' @family low-level functions operating on numeric vectors.
#'
v_insert_hinges <- function(x, y, h) {
  # sanitize 'hinges'
  h <- setdiff(h, x)
  h <- h[h > x[1] & h < x[length(x)]]
  h <- sort(h)
  # save lengths
  j <- length(h)
  k <- length(x)
  # compute the new idx values after shift caused by insertions
  idxs.in <- findInterval(h, x)
  idxs.out <- idxs.in + 1:j
  idxs.diff <- diff(c(0,idxs.in,k))
  idxs.map <- 1:k + rep(0:j, idxs.diff)
  # we use recycling to interpolate all values and insert them into the gaps
  if (is.numeric(y)) {
    # allocate vectors for the output
    x.out <- numeric(j + k)
    y.out <- numeric(j + k)
    # we copy everything that does not require interpolation
    x.out[idxs.map] <- x
    y.out[idxs.map] <- y
    x.out[idxs.out] <- h
    # we fill the hinges by interpolation
    y.out[idxs.out] <- y[idxs.in + 1] -
      (x[idxs.in + 1] - x.out[idxs.out]) /
      (x[idxs.in + 1] - x[idxs.in]) *
      (y[idxs.in + 1] - y[idxs.in])
  } else if (length(unique(y)) == 1) {
    y.out <- rep(y[1], j + k)
  } else {
    y.out <- rep(y[NA], j + k)
    y.out[idxs.map] <- y
  }
  y.out
}

#' Overwrite spectral data values at existing wavelength values.
#'
#' Overwriting spectral data with interpolated values at wavelengths values
#' containing bad data is needed when cleaning spectral data.
#' This function differs from \code{insert_hinges()} in that it returns a vector
#' of \code{y} values instead of a \code{tibble}.
#'
#' @param x numeric vector (sorted in increasing order).
#' @param y numeric vector.
#' @param h a numeric vector giving the wavelengths at which the y values should
#'   be replaced by interpolation, no interpolation is indicated by an empty
#'   numeric vector (\code{numeric(0)}).
#'
#' @return A numeric vector with the numeric values of \code{y} with values
#'   at the hinges replaced by interpolation of neighbours.
#'
#' @keywords internal.
#'
#' @family low-level functions operating on numeric vectors.
#'
v_replace_hinges <- function(x, y, h) {
  # sanitize 'hinges'
  h.idxs <- which(x %in% h)
  y.out <- y
  # we fill the hinges by interpolation
  if (is.numeric(y)) {
    y.out[h.idxs] <- y[h.idxs + 1] -
      (x[h.idxs + 1] - x[h.idxs - 1]) /
        (x[h.idxs + 1] - x[h.idxs]) *
        (y[h.idxs + 1] - y[h.idxs - 1])
  } else if (length(unique(y)) == 1) {
    y.out[h.idxs] <- y[1]
  } else {
    y.out[h.idxs] <- NA
  }
  y.out
}
