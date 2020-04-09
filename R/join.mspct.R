#' Join all spectra in a collection
#'
#' Join all the spectra contained in a homogenous collection, returning a data
#' frame with spectral-data columns named according to the names of the spectra
#' in the collection. By default a full join is done, filling the spectral
#' data for missing wave lengths in individual spectra with \code{NA}.
#'
#' @param x A generic_mspct object, or an object of a class derived from
#'   generic_mspct.
#' @param unit.out character Allowed values "energy", and "photon", or its alias
#'   "quantum".
#' @param qty.out character Allowed values "transmittance", and "absorbance".
#' @param type character Type of join: "left", "right", "inner" or "full"
#'   (default). See details for more information.
#' @param ... ignored (possibly used by derived methods).
#'
#' @return An object of class dataframe, with the spectra joined by wave length,
#'   with rows in addition sorted by wave length (variable \code{w.length}).
#'
#' @note Currently only generic_spct, source_mspct, response_mspct,
#'   filter_mspct, reflector_mspct and object_mspct classes have this method
#'   implemented.
#'
#' @export
#'
#' @family conversion of collections of spectra
#'
join_mspct <- function(x, type, ...) UseMethod("join_mspct")

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.default <- function(x, type = "full", ...) {
  stop("'join_mspct()' is only implemented for collections of spectra, ",
       "use 'plyr::join_all()' for lists of data frames.")
}

#' @describeIn join_mspct
#'
#' @param col.name character, name of the column in the spectra to be preserved,
#'   in addition to "w.length".
#'
#' @export
#'
join_mspct.generic_mspct <- function(x,
                                     type = "full",
                                     col.name, ...) {
  # if needed could be added as additional formal parameters
  by <- "w.length"
  match <- "first"

  if (length(x) == 0L) {
    return(data.frame())
  }
  names <- names(x)
  stopifnot(length(names) == length(x))
  rmDerivedMspct(x)
  col.selector <- c("w.length", col.name)
  for (i in names) {
    x[[i]] <- as.data.frame(x[[i]])[col.selector]
#    x[[i]] <- plyr::rename(x[[i]], replace = c(parse(col.name) = i))) silently failing!!
    col.names <- names(x[[i]])
    names(x[[i]])[col.names == col.name] <- i
  }
  if (length(x) == 1L) {
    z <- as.data.frame(x[[i]])
  } else {
    z <- plyr::join_all(dfs = x, by = by, type = type, match = match)
  }
  dplyr::arrange(z, .data$w.length)
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.source_mspct <- function(x,
                                    type = "full",
                                    unit.out = "energy",
                                    ...) {
  # if needed could be added as additional formal parameters
  by <- "w.length"
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
    z <- as.data.frame(x[[i]])
  } else {
    z <- plyr::join_all(dfs = x, by = by, type = type, match = match)
  }
  dplyr::arrange(z, .data$w.length)
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.response_mspct <- function(x,
                                      type = "full",
                                      unit.out = "energy",
                                      ...) {
  # if needed could be added as additional formal parameters
  by <- "w.length"
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
    z <- as.data.frame(x[[i]])
  } else {
    z <- plyr::join_all(dfs = x, by = by, type = type, match = match)
  }
  dplyr::arrange(z, .data$w.length)
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.filter_mspct <- function(x,
                                    type = "full",
                                    qty.out = "transmittance",
                                    ...) {
  # if needed could be added as additional formal parameters
  by <- "w.length"
  match <- "first"

  if (length(x) == 0L) {
    return(data.frame())
  }
  names <- names(x)
  stopifnot(length(names) == length(x))
  if (qty.out == "transmittance") {
    x <- any2T(x, action = "replace")
    rmDerivedMspct(x)
    for (i in names) {
      x[[i]] <- as.data.frame(x[[i]])[c("w.length", "Tfr")]
      x[[i]] <- plyr::rename(x[[i]], c(Tfr = i))
    }
  } else if (qty.out == "absorbance") {
    x <- any2A(x, action = "replace")
    rmDerivedMspct(x)
    for (i in names) {
      x[[i]] <- as.data.frame(x[[i]])[c("w.length", "A")]
      x[[i]] <- plyr::rename(x[[i]], c(A = i))
    }
  } else if (qty.out == "absorptance") {
    x <- any2Afr(x, action = "replace")
    rmDerivedMspct(x)
    for (i in names) {
      x[[i]] <- as.data.frame(x[[i]])[c("w.length", "Afr")]
      x[[i]] <- plyr::rename(x[[i]], c(Afr = i))
    }
  } else {
    stop("Unit out '", qty.out, "' unknown")
  }
  if (length(x) == 1L) {
    z <- as.data.frame(x[[i]])
  } else {
    z <- plyr::join_all(dfs = x, by = by, type = type, match = match)
  }
  dplyr::arrange(z, .data$w.length)
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.reflector_mspct <- function(x,
                                       type = "full",
                                       ...) {
  # if needed could be added as additional formal parameters
  by <- "w.length"
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
    z <- as.data.frame(x[[i]])
  } else {
    z <- plyr::join_all(dfs = x, by = by, type = type, match = match)
  }
  dplyr::arrange(z, .data$w.length)
}

#' @describeIn join_mspct
#'
#' @export
#'
join_mspct.object_mspct <- function(x,
                                    type = "full",
                                    qty.out,
                                    ...) {
  # # if needed could be added as additional formal parameters
  # by <- "w.length"
  # match <- "first"
  switch(qty.out,
         "transmittance" = join_mspct(as.filter_mspct(x), type = type, qty.out = qty.out, ...),
         "absorbance" = join_mspct(as.filter_mspct(x), type = type, qty.out = qty.out, ...),
         "absorbtance" = join_mspct(as.filter_mspct(x), type = type, qty.out = qty.out, ...),
         "reflectance" = join_mspct(as.reflector_mspct(x), type = type, ...),
         stop("'qty.out = ", qty.out, " not implemented.")
           )
}
