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
#' @return a collection of spectra in the case of \code{msmsply}, belonging to a
#'   different class than \code{mspct} if \code{.fun} modifies the class of the
#'   member spectra.
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
  rmDerivedMspct(mspct)
  y <- plyr::llply(.data = mspct,
                   .fun = .fun,
                   ...,
                   .parallel = .parallel,
                   .paropts = .paropts)

  stopifnot(length(y) == length(mspct))

  result.class <- shared_member_class(y)[1]
  stopifnot(length(result.class) == 1) # not a spectrum

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

  if ((is.logical(idx) && idx) ||
        (is.null(idx) && !any(is.null(names(mspct))))) {
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
  f.name <- gsub(paste0("\\.source_spct$|\\.response_spct$|\\.filter_spct$|",
                        "\\.reflector_spct$|\\.generic_spct$|\\.solute_spct$"),
                 "", x = f.name)

  if (f.name %in% c("min", "max", "range",
                    "min_wl", "max_wl", "range_wl",
                    "spread", "midpoint", "stepsize",
                    "spread_wl", "midpoint_wl", "stepsize_wl",
                    "getWhenMeasured", "getWhereMeasured",
                    "irrad", "q_irrad", "e_irrad",
                    "fluence", "q_fluence", "e_fluence",
                    "q_ratio", "e_ratio", "eq_ratio", "qe_ratio",
                    "q_fraction", "e_fraction",
                    "response", "q_response", "e_response",
                    "absorbance", "absorptance", "transmittance",
                    "reflectance")) {
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
  } else {# make new names using function name
    if (is.null(.idx)) qty.names <- names(z) else qty.names <- names(z)[-1]
    qty.names <- paste(f.name,
                       gsub(" ", "", qty.names),
                       sep = "_")
  }

  if (!is.null(qty.names)) {
    if (is.null(.idx)) {
      names(z) <- qty.names
    } else {
      names(z)[-1] <- qty.names
    }
  }

  comment(z) <-
    paste("Applied function: '", f.name, "'.\n", sep = "", comment(mspct))

  mspct.ncol <- ncol(mspct)
  if (is.na(mspct.ncol)) {
    mspct.ncol <- 1
  }
  mspct.nrow <- nrow(mspct)
  if (is.na(mspct.nrow)) {
    mspct.nrow <- length(mspct) %/% mspct.ncol
  }
  mspct.byrow <- attr(mspct, "mspct.byrow", exact = TRUE)
  if (is.null(mspct.byrow)) {
    mspct.nrow <- FALSE
  }
  if (mspct.ncol > 1) {
    if (mspct.byrow) {
      z[["col"]] <- rep(1:mspct.ncol, mspct.nrow)
      z[["row"]] <- rep(1:mspct.nrow, rep(mspct.ncol, mspct.nrow))
    } else {
      z[["col"]] <- rep(1:mspct.ncol, rep(mspct.nrow, mspct.ncol))
      z[["row"]] <- rep(1:mspct.nrow, mspct.ncol)
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
                   .paropts = .paropts)

  names(z) <- names(mspct)

  f.name <- as.character(substitute(.fun))
  f.name <- gsub(paste0("\\.source_spct$|\\.response_spct$|\\.filter_spct$|",
                        "\\.reflector_spct$|\\.generic_spct$|\\.solute_spct$"),
                 "", x = f.name)

  comment(z) <-
    paste("Applied function: '", f.name, "'.\n", sep = "", comment(mspct))

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
  f.name <- gsub(paste0("\\.source_spct$|\\.response_spct$|\\.filter_spct$|",
                        "\\.reflector_spct$|\\.generic_spct$|\\.solute_spct$"),
                 "", x = f.name)

  comment(z) <-
    paste("Applied function: '", f.name, "'.\n", sep = "", comment(mspct))

  z
}
