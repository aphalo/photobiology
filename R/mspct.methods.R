# Apply -------------------------------------------------------------------

#' Multi-spct transform methods
#'
#' Apply a function or operator to a collection of spectra.
#'
#' @param mspct an object of class generic_mspct or a derived class
#' @param .fun a function
#' @param ... other arguments passed to .fun
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @return a collection of spectra in the case of \code{msmsply}
#'
#' @export
#'
msmsply <- function(mspct, .fun, ...,
                    .parallel = FALSE, .paropts = NULL) {
  stopifnot(is.any_mspct(mspct))
  mspct.class <- class(mspct)
  byrow <- attr(mspct, "mspct.byrow", exact = TRUE)
  dim <- dim(mspct)
  ncol <- ncol(mspct)
  # llply returns a matrix for classes derived from list
  #
  rmDerivedMspct(mspct)
  y <- plyr::llply(.data = mspct,
                   .fun = .fun,
                   ...,
                   .parallel = .parallel,
                   .paropts = .paropts)

  stopifnot(length(y) == length(mspct))

  if (length(y) > 1) {
    result.class <- shared_member_class(y)[1]
  } else {
    result.class <- mspct.class[1]
  }
  stopifnot(length(result.class) == 1)

  generic_mspct(l = y,
                class = result.class,
                byrow = byrow,
                ncol = ncol,
                dim = dim)
}

#' @rdname  msmsply
#'
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @param col.names character Names to be used for data columns.
#'
#' @return a data frame in the case of \code{msdply}
#'
#' @export
#'
msdply <- function(mspct, .fun, ..., idx = NULL, col.names = NULL,
                   .parallel = FALSE, .paropts = NULL) {
  stopifnot(is.any_mspct(mspct))

  if ( (is.logical(idx) && idx) ||
       (is.null(idx) && !any(is.null(names(mspct)))) ) {
    .idx <- "spct.idx"
  } else if (is.logical(idx) && !idx) {
    .idx <- NULL
  } else {
    .idx <- idx
  }

  z <- plyr::ldply(.data = mspct,
                   .fun = .fun,
                   ...,
                   .id = .idx,
                   .parallel = .parallel,
                   .paropts = .paropts)

  f.name <- as.character(substitute(.fun))

  if (f.name %in% c("min", "max", "range",
                    "min_wl", "max_wl", "range_wl",
                    "spread", "midpoint", "stepsize",
                    "spread_wl", "midpoint_wl", "stepsize_wl",
                    "getWhenMeasured", "getWhereMeasured",
                    "irrad", "q_irrad", "e_irrad",
                    "fluence", "q_fluence", "e_fluence",
                    "q_ratio", "e_ratio", "eq_ratio", "qe_ratio",
                    "response", "q_response", "e_response",
                    "absorbance")) {
    qty.names <-
      switch(f.name,
             min = "min.wl", min_wl = "min.wl",
             max = "max.wl", max_wl = "max.wl",
             range = c("min.wl", "max.wl"), range_wl = c("min.wl", "max.wl"),
             spread = "spread.wl", spread_wl = "spread.wl",
             midpoint = "midpoint.wl", midpoint_wl = "midpoint.wl",
             stepsize = c("min.step.wl", "max.step.wl"),
             stepsize_wl = c("min.step.wl", "max.step.wl"),
             getWhenMeasured = "when.measured",
             NULL # default: use name of returned numeric values
    )
  } else if (!is.null(col.names) &&
             !any(col.names == "") &&
             !any(is.na(col.names)) &&
             length(col.names) == length(names(z)) - 1) {
    qty.names <- col.names
#  } else if (any(c("total", "mean", "contrib", "particip") %in% tolower(names(z)))) {
  } else {# make new names using function name
    if (is.null(.idx)) qty.names <- names(z) else qty.names <- names(z)[-1]
    qty.names <- paste(f.name,
                       gsub(" ", "", qty.names),
                       sep = "_")
#    qty.names <- NULL
  }

  if (!is.null(qty.names)) {
    if (is.null(.idx)) {
      names(z) <- qty.names
    } else {
      names(z)[-1] <- qty.names
    }
  }

  comment(z) <- paste("Applied function: '", f.name, "'.\n", sep = "", comment(mspct))

  mspct.nrow <- nrow(mspct)
  mspct.ncol <- ncol(mspct)
  mspct.byrow <- attr(mspct, "mspct.byrow", exact = TRUE)
  if (is.null(mspct.byrow)) {
    mspct.nrow <- FALSE
  }
  if (mspct.ncol > 1) {
    if (mspct.byrow) {
      z$col <- rep(1:mspct.ncol, mspct.nrow)
      z$row <- rep(1:mspct.nrow, rep(mspct.ncol, mspct.nrow))
    } else {
      z$col <- rep(1:mspct.ncol, rep(mspct.nrow, mspct.ncol))
      z$row <- rep(1:mspct.nrow, mspct.ncol)
    }
  }
  tibble::as_tibble(z)
}

#' @rdname  msmsply
#'
#' @return a list in the case of \code{mslply}
#'
#' @export
#'
mslply <- function(mspct, .fun, ...,
                   .parallel = FALSE, .paropts = NULL) {
  stopifnot(is.any_mspct(mspct))

  # llply returns a matrix for classes derived from list
  #
  rmDerivedMspct(mspct)
  z <- plyr::llply(.data = mspct,
                   .fun = .fun,
                   ...,
                   .parallel = .parallel,
                   .paropts = .paropts )

  names(z) <- names(mspct)

  f.name <- as.character( substitute(.fun))

  comment(z) <- paste("Applied function: '", f.name, "'.\n", sep = "", comment(mspct))

  z
}


#' @rdname  msmsply
#'
#' @param .drop should extra dimensions of length 1 in the output be dropped,
#'   simplifying the output. Defaults to TRUE
#' @return an vector in the case of \code{msaply}
#'
#' @export
#'
msaply <- function(mspct, .fun, ..., .drop = TRUE,
                   .parallel = FALSE, .paropts = NULL) {
  stopifnot(is.any_mspct(mspct))

  # As many of our summary functions return numeric values with names and other
  # attributes they need to be removed for dply::lapply to accept them.
  .ffun <- function(mspct, ...) {
    z <- .fun(mspct, ...)
    if (is.numeric(z)) {
      z <- as.numeric(z)
    }
    z
  }

  z <- plyr::laply(.data = mspct,
                   .fun = .ffun,
                   ...,
                   .drop = .drop,
                   .parallel = .parallel,
                   .paropts = .paropts)

  f.name <- as.character(substitute(.fun))

  comment(z) <- paste("Applied function: '", f.name, "'.\n", sep = "", comment(mspct))

  z
}



#' Get the "mspct.version" attribute
#'
#' Function to read the "mspct.version" attribute of an existing generic_mspct
#' object.
#'
#' @param x a generic_mspct object
#'
#' @return numeric value
#'
#' @note if x is not a \code{generic_mspct} object, \code{NA} is returned,
#'   and if it the attribute is missing, zero is returned with a warning.
#'
#' @export
#'
getMspctVersion <- function(x) {
  if (is.any_mspct(x)) {
    version <- attr(x, "mspct.version", exact = TRUE)
    if (is.null(version)) {
      # need to handle objects created with old versions
      version <- 0L
    }
  } else {
    version <- NA
  }
  version
}

#' Check that the "mspct.version" attribute is set
#'
#' Function to check the "mspct.version" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_mspct object
#'
#' @return numeric value
#'
#' @note if x is not a \code{generic_mspct} object, \code{NA} is returned,
#'   and if it the attribute is missing, zero is returned with a warning.
#'
#' @keywords internal
#'
checkMspctVersion <- function(x) {
  version <- getMspctVersion(x)
  stopifnot(!is.na(version))
  if (version < 1L) {
    warning("The object '", as.character(substitute(x)),
            "' is corrupted")
  }
}



# convolution ---------------------------------------------------------------

#' Convolve function for collections of spectra
#'
#' Convolve function for collections of spectra which applies an operation on
#' all the individual members of the collection(s) of spectra.
#'
#' @param e1 an object of class \code{generic_mspct} or \code{generic_scpt} or
#'   \code{numeric}
#' @param e2 an object of class \code{generic_mspct} or \code{generic_scpt} or
#'   \code{numeric}
#' @param oper function, usually but not necessarily an operator with two
#'   arguments.
#' @param sep character Used when pasting the names of members of \code{e1} and
#'   \code{e2} to form the names of members of the returned collection of
#'   spectra.
#' @param ... additional arguments passed to \code{oper} if present.
#'
#' @note At least one of e1 and e2 must be a \code{generic_mspct} object or
#'   derived.
#'
#' @export
#'
#' @family math operators and functions
#'
convolve_each <- function(e1, e2, oper = `*`, sep = "_", ...) {
  e3 <- list()
  if (is.any_mspct(e1) & !is.any_mspct(e2)) {
    for (spct.name in names(e1)) {
      e3[[spct.name]] <- oper(e1[[spct.name]], e2, ...)
    }
    z <- generic_mspct(e3, class = shared_member_class(e3),
                       ncol = ncol(e1),
                       byrow = attr(e1, "mspct.byrow", exact = TRUE))
  } else if (!is.any_mspct(e1) & is.any_mspct(e2)) {
    for (spct.name in names(e2)) {
      e3[[spct.name]] <- oper(e1, e2[[spct.name]], ...)
    }
    z <- generic_mspct(e3, class = shared_member_class(e3),
                       ncol = ncol(e2),
                       byrow = attr(e2, "mspct.byrow", exact = TRUE))
  } else if (is.any_mspct(e1) & is.any_mspct(e2)) {
    for (spct.name1 in names(e1)) {
      for (spct.name2 in names(e2)) {
        combined.name <- paste(spct.name1, spct.name2, sep = sep)
        e3[[combined.name]] <- oper(e1[[spct.name1]], e2[[spct.name2]], ...)
      }
     }
    z <- generic_mspct(e3, class = shared_member_class(e3),
                       ncol = nrow(e2),
                       byrow = FALSE)
    attr(z, "mspct.dimnames")  <- list(names(e1), names(e2))
  } else {
    stop("At least one of 'e1' and 'e2' should be a collection of spectra.")
  }
  z
}


# utility functions for attributes ----------------------------------------

#' Copy attributes from members of a generic_mspct
#'
#' Copy the when.measured, where.measured or what.measured attribute from
#' members of a generic_mspct object into a tibble or data.frame.
#'
#' @param mspct generic_mspct Any collection of spectra.
#' @param tb tibble or data.frame to which to add the data (optional).
#' @param col.names named character vector Name(s) of column(s) to create.
#'   Values are the names of the attributes to copy, while if named, the names
#'   provide the name for the column.
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#'
#' @return A tibble With the metadata attributes in separate new variables.
#'
#' @details The attributes are copied to a column in a tibble or data frame. If
#'   the \code{tb} formal parameter receives \code{NULL} as argument, a new
#'   \code{tibble} will be created. If an existing \code{data.frame} or
#'   \code{tibble} is passed as argument, new columns are added to it. However,
#'   the number of rows in the argument passed to \code{tb} must match the
#'   number of spectra in the argument passed to \code{mspct}. If the argument
#'   to \code{col.names} is aa named vector, with the names of members matching
#'   the names of attributes, then the values are used as names for the columns
#'   created. This permits setting any valid name for the new columns. If the
#'   vector passed to \code{col.names} has no names, then the values are
#'   interpreted as the names of the attributes to add, and also used as names
#'   for the new columns.
#'
#' @note Currently supported attributes are \code{"when.measured"},
#'   \code{"what.measured"} and \code{"where.measured"}. In the case of
#'   \code{"where.measured"}, which has different components the name
#'   \code{"where.measured"} is ignored, but instead the following names are
#'   recognized: \code{"lon"} and \code{"lat"} for creating numeric columns of
#'   longitudes and latitudes respectively, and \code{"geocode"} for creating a
#'   column of data frames, in which case, if \code{tb} is not already a
#'   \code{tibble} it is converted into one before adding the new column.  The
#'   order of the first two arguments is reversed in \code{add_attr2tb()}
#'   compared to the other functions. This is to allow its use in 'pipes', while
#'   the functions for single attributes are expected to be used mostly to
#'   create new tibbles.
#'
#' @examples
#'
#' library(dplyr)
#'
#' my.mspct <- source_mspct(list(sun1 = sun.spct, sun2 = sun.spct * 2))
#' q_irrad(my.mspct) %>%
#'   add_attr2tb(my.mspct, c(lat = "latitude",
#'                           lon = "longitude",
#'                           when.measured = "time"))
#'
#' when_measured2tb(my.mspct)
#'
#' @export
#'
add_attr2tb <- function(tb,
                        mspct,
                        col.names = NULL,
                        idx = "spct.idx") {
  if (length(col.names) < 1L) {
    return(tb)
  }
  if (length(names(col.names)) < 1L) {
    names(col.names) <- col.names
  }
  attributes <- intersect(names(col.names),
                          c("geocode",
                            "lon",
                            "lat",
                            "when.measured",
                            "what.measured"))
  if (length(attributes) < length(col.names)) {
    warning("Unrecognized attribute(s) '",
            paste(setdiff(names(col.names), attributes), collapse = ", "),
            "' where skipped.")
  }
  force(tb)
  for (a in attributes) {
    tb <-
      switch(a,
             geocode = geocode2tb(mspct = mspct,
                                  tb = tb,
                                  col.names = col.names["geocode"],
                                  idx = idx),
             lon = lon2tb(mspct = mspct,
                          tb = tb,
                          col.names = col.names["lon"],
                          idx = idx),
             lat = lat2tb(mspct = mspct,
                          tb = tb,
                          col.names = col.names["lat"],
                          idx = idx),
             when.measured = when_measured2tb(mspct = mspct,
                                              tb = tb,
                                              col.names = col.names["when.measured"],
                                              idx = idx),
             what.measured = what_measured2tb(mspct = mspct,
                                              tb = tb,
                                              col.names = col.names["what.measured"],
                                              idx = idx),
             tb)
  }
  tb
}

#' @rdname add_attr2tb
#'
#' @export
#'
when_measured2tb <- function(mspct,
                             tb = NULL,
                             col.names = "when.measured",
                             idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  when.tb <- getWhenMeasured(mspct, idx = idx)
  if (col.names == "when.measured") {
    if (length(names(col.names)) == 1L) {
      # syntax like in dplyr::rename()
      names(when.tb)[2L] <- names(col.names)
    }
  } else {
    # use value directly for column name
    names(when.tb)[2L] <- col.names
  }
  if (is.null(tb)) {
    when.tb
  } else {
    dplyr::full_join(tb, when.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
lonlat2tb <- function(mspct,
                      tb = NULL,
                      col.names = c("lon", "lat"),
                      idx = "spct.idx") {
  stopifnot(length(col.names) == 2L)
  lonlat.tb <- getWhereMeasured(mspct, idx = idx)[c(idx, "lon", "lat")]
  if (col.names[1L] == "lon" && col.names[2L] == "lat") {
    if (length(names(col.names)) == 2L) {
      # syntax like in dplyr::rename()
      names(lonlat.tb)[2L:3L] <- names(col.names)
    }
  } else {
    # use value directly for column name
    names(lonlat.tb)[2L:3L] <- col.names
  }
  if (is.null(tb)) {
    lonlat.tb
  } else {
    dplyr::full_join(tb, lonlat.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
lon2tb <- function(mspct,
                   tb = NULL,
                   col.names = "lon",
                   idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  lon.tb <- getWhereMeasured(mspct, idx = idx)[c(idx, "lon")]
  if (col.names == "lon") {
    if (length(names(col.names)) == 1L) {
      # syntax like in dplyr::rename()
      names(lon.tb)[2L] <- names(col.names)
    }
  } else {
    # use value directly for column name
    names(lon.tb)[2L] <- col.names
  }
  if (is.null(tb)) {
    lon.tb
  } else {
    dplyr::full_join(tb, lon.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
lat2tb <- function(mspct,
                   tb = NULL,
                   col.names = "lat",
                   idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  lat.tb <- getWhereMeasured(mspct, idx = idx)[c(idx, "lat")]
  if (col.names == "lat") {
    if (length(names(col.names)) == 1L) {
      # syntax like in dplyr::rename()
      names(lat.tb)[2L] <- names(col.names)
    }
  } else {
    # use value directly for column name
    names(lat.tb)[2L] <- col.names
  }
  if (is.null(tb)) {
    lat.tb
  } else {
    dplyr::full_join(tb, lat.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
geocode2tb <- function(mspct,
                       tb = NULL,
                       col.names = "geocode",
                       idx = "spct.idx") {
  if (is.null(tb)) {
    tb <- tibble::tibble(spct.idx = factor(names(mspct)))
  } else {
    stopifnot(nrow(tb) == length(mspct))
  }
  if (length(names(col.names)) < 1L) {
    names(col.names) <- col.names
  }
  geocodes <- list()
  row <- 1L
  for (x in mspct) {
    geocodes[[row]] <- getWhereMeasured(x)
    row <- row + 1L
  }
  if (!tibble::is_tibble(tb)) {
    tb <- tibble::as_tibble(tb)
  }
  tb[[col.names["geocode"]]] <- geocodes
  tb
}

#' @rdname add_attr2tb
#'
#' @export
#'
what_measured2tb <- function(mspct,
                             tb = NULL,
                             col.names = "what.measured",
                             idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  what.tb <- getWhatMeasured(mspct, idx = idx)
  if (col.names == "what.measured") {
    if (length(names(col.names)) == 1L) {
      # syntax like in dplyr::rename()
      names(what.tb)[2L] <- names(col.names)
    }
  } else {
    # use value directly for column name
    names(what.tb)[2L] <- col.names
  }
  if (is.null(tb)) {
    what.tb
  } else {
    dplyr::full_join(tb, what.tb, by = idx)
  }
}
