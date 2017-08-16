#' Join all spectra in a collection
#'
#' Join all the spectra contained in a homegenous collection, returning a data
#' frame with spectral-data columns named according to the names of the spectra
#' in the collection.
#'
#' @param x A generic_mspct object, or an object of a class derived from
#'   generic_mspct.
#' @param unit.out character Allowed values "energy", and "photon", or its alias
#'   "quantum".
#' @param qty.out character Allowed values "transmittance", and "absorbance".
#' @param ... currently ignored.
#'
#' @return A object of class dataframe, with the spectra joined by wavelength.
#'
#' @note Currently only source_mspct, response_mspct, filter_mspct, and
#'   reflector_mspct classes have this method implemented.
#'
#' @export
#'
join_mspct <- function(x, ...) UseMethod("join_mspct")

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.default <- function(x, ...) {
  stop("'join_mspct()' is only implemented for collections of spectra, ",
       "use 'plyr::join_all()' for lists of data frames.")
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.generic_spct <- function(x, ...) {
  warning("'join_mspct()' only implemented for homogeneous collections of spectra.")
  data.frame()
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.source_mspct <- function(x,
                                    unit.out = "energy",
                                    ...) {
  # if needed could be added as additional formal parameters
  by <- "w.length"
  type <- "full"
  match <- "first"

  if (length(x) == 0L) {
    return(data.frame())
  }
  names <- names(x)
  stopifnot(length(names) == length(x))
  if (unit.out == "energy") {
    x <- q2e(x, action = "replace")
    rmDerivedMspct(x)
    for (i in names) {
      x[[i]] <- as.data.frame(x[[i]])[c("w.length", "s.e.irrad")]
      x[[i]] <- plyr::rename(x[[i]], c(s.e.irrad = i))
    }
  } else if (unit.out %in% c("photon", "quantum")) {
    x <- e2q(x, action = "replace")
    rmDerivedMspct(x)
    for (i in names) {
      x[[i]] <- as.data.frame(x[[i]])[c("w.length", "s.q.irrad")]
      x[[i]] <- plyr::rename(x[[i]], c(s.q.irrad = i))
    }
  } else {
    stop("Unit out '", unit.out, "' unknown")
  }
  if (length(x) == 1L) {
    as.data.frame(x[[i]])
  } else {
    plyr::join_all(dfs = x, by = by, type = type, match = match)
  }
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.response_mspct <- function(x,
                                      unit.out = "energy",
                                      ...) {
  # if needed could be added as additional formal parameters
  by <- "w.length"
  type <- "full"
  match <- "first"

  if (length(x) == 0L) {
    return(data.frame())
  }
  names <- names(x)
  stopifnot(length(names) == length(x))
  if (unit.out == "energy") {
    x <- q2e(x, action = "replace")
    rmDerivedMspct(x)
    for (i in names) {
      x[[i]] <- as.data.frame(x[[i]])[c("w.length", "s.e.response")]
      x[[i]] <- plyr::rename(x[[i]], c(s.e.response = i))
    }
  } else if (unit.out %in% c("photon", "quantum")) {
    x <- e2q(x, action = "replace")
    rmDerivedMspct(x)
    for (i in names) {
      x[[i]] <- as.data.frame(x[[i]])[c("w.length", "s.q.response")]
      x[[i]] <- plyr::rename(x[[i]], c(s.q.response = i))
    }
  } else {
    stop("Unit out '", unit.out, "' unknown")
  }
  if (length(x) == 1L) {
    as.data.frame(x[[i]])
  } else {
    plyr::join_all(dfs = x, by = by, type = type, match = match)
  }
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.filter_mspct <- function(x,
                                    qty.out = "transmittance",
                                    ...) {
  # if needed could be added as additional formal parameters
  by <- "w.length"
  type <- "full"
  match <- "first"

  if (length(x) == 0L) {
    return(data.frame())
  }
  names <- names(x)
  stopifnot(length(names) == length(x))
  if (qty.out == "transmittance") {
    x <- A2T(x, action = "replace")
    rmDerivedMspct(x)
    for (i in names) {
      x[[i]] <- as.data.frame(x[[i]])[c("w.length", "Tfr")]
      x[[i]] <- plyr::rename(x[[i]], c(Tfr = i))
    }
  } else if (qty.out == "absorbance") {
    x <- T2A(x, action = "replace")
    rmDerivedMspct(x)
    for (i in names) {
      x[[i]] <- as.data.frame(x[[i]])[c("w.length", "A")]
      x[[i]] <- plyr::rename(x[[i]], c(A = i))
    }
  } else {
    stop("Unit out '", qty.out, "' unknown")
  }
  if (length(x) == 1L) {
    as.data.frame(x[[i]])
  } else {
    plyr::join_all(dfs = x, by = by, type = type, match = match)
  }
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.reflector_mspct <- function(x,
                                       ...) {
  # if needed could be added as additional formal parameters
  by <- "w.length"
  type <- "full"
  match <- "first"

  if (length(x) == 0L) {
    return(data.frame())
  }
  names <- names(x)
  stopifnot(length(names) == length(x))
  rmDerivedMspct(x)
  for (i in names) {
    x[[i]] <- as.data.frame(x[[i]])[c("w.length", "Rfr")]
    x[[i]] <- plyr::rename(x[[i]], c(Rfr = i))
  }
  if (length(x) == 1L) {
    as.data.frame(x[[i]])
  } else {
    plyr::join_all(dfs = x, by = by, type = type, match = match)
  }
}

