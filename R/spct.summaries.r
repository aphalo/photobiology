# print -------------------------------------------------------------------

#' Print spectral objects
#'
#' Print methods for objects of spectral classes, including collections of
#' spectra.
#'
#' @param x An object of one of the summary classes for spectra.
#' @param ... not used in current version.
#' @param attr.simplify logical If all members share the same attribute value
#'   return one copy instead of a data.frame, list or vector.
#' @param n	Number of rows to show. If NULL, the default, will print all rows if
#'   less than option \code{dplyr.print_max}. Otherwise, will print
#'   \code{dplyr.print_min} rows.
#' @param width	Width of text output to generate. This defaults to NULL, which
#'   means use getOption("width") and only display the columns that fit on one
#'   screen. You can also set option(dplyr.width = Inf) to override this default
#'   and always print all columns.
#'
#' @return Returns \code{x} invisibly.
#'
#' @export
#'
#' @details This is a wrapper on the print method for tibbles, with
#' additional information in the header. Currently, \code{width} applies only to
#' the table of data. To print only the header and a subset of data rows
#' starting from the shortest wavelengths pass a positive integer to \code{n}.
#'
#' Objects are printed as is, ignoring the current settings of R options
#' \code{photobiology.radiation.unit} and \code{photobiology.filter.qty}.
#'
#' @note Methods \code{\link[utils]{head}()}, \code{tail()} and
#'   \code{\link{head_tail}()} give additional flexibility on the selection
#'   of rows to print, while preserving the metadata. The information shown
#'   for wavelengths is in contrast to when using print that for the displayed
#'   rows
#'
#' @seealso \code{\link[tibble]{formatting}}.
#'
#' @name print.generic_spct
#'
#' @examples
#'
#' print(sun.spct)
#' print(sun.spct, n = 5)
#'
#' print(head(sun.spct, n = 5))
#' print(tail(sun.spct, n = 5))
#' print(head_tail(sun.spct, n = 5))
#'
#' print(q2e(sun.spct, action = "replace"))
#' print(e2q(sun.spct, action = "replace"))
#'
#' print(polyester.spct)
#' print(any2A(polyester.spct))
#' print(any2Afr(polyester.spct))
#'
#' print(two_filters.spct)
#'
print.generic_spct <- function(x,
                               ...,
                               attr.simplify = TRUE,
                               n = NULL,
                               width = NULL)
{
  # Skip checks of validity as we are only printing
  prev_state <- disable_check_spct()
  on.exit(set_check_spct(prev_state), add = TRUE)

  cat("Object: ", class_spct(x)[1], " ", dplyr::dim_desc(x), "\n", sep = "")
  if (nrow(x)) {
    m.wl <- getMultipleWl(x)
    if (m.wl > 1) {
      cat("containing ", m.wl, " spectra in long form\n", sep = "")
    }
    cat("Wavelength range ",
        paste(signif(wl_range(x), 8), sep = "", collapse = "-"), " nm, step ",
        paste(unique(signif(stepsize(x), 7)), sep = "", collapse = "-"),
        " nm \n", sep = "")
  }
  what.measured <- getWhatMeasured(x, simplify = attr.simplify)
  if (!any(is.na(what.measured))) {
    if (!is.list(what.measured)) {
      what.measured <- list(what.measured)
    }
    names <- names(what.measured)
    if (is.null(names)) {
      names <- "Label: "
    } else {
      names <- paste(names, "label: ")
    }
    cat(paste(names,
              what.measured,
              sep = "", collapse = "\n"), "\n")
  }
  when.measured <- getWhenMeasured(x, simplify = attr.simplify)
  if (!any(is.na(when.measured))) {
    if (!is.list(when.measured)) {
      when.measured <- list(when.measured)
    }
    names <- names(when.measured)
    if (is.null(names)) {
      names <- "Measured on "
    } else {
      names <- paste(names, "measured on ")
    }
    cat(paste(names,
              sapply(when.measured, as.character), " UTC",
              sep = "", collapse = "\n"), "\n")
  }
  where.measured <- getWhereMeasured(x, simplify = attr.simplify)
  if (!any(is.na(where.measured))) {
    if (is.data.frame(where.measured)) {
      where.measured <- list(where.measured)
    }
    names <- names(where.measured)
    if (is.null(names)) {
      names <- "Measured at "
    } else {
      names <- paste(names, "measured at ")
    }
    if (all(sapply(where.measured,
                   function(x) {"address" %in% names(x) && !any(is.na(x[["address"]]))}))) {
      cat(paste(names,
                sapply(where.measured, `[[`, i = "lat"), " N, ",
                sapply(where.measured, `[[`, i = "lon"), " E; ",
                stringr::str_trunc(sapply(where.measured, `[[`, i = "address"), width = 30, side = "right"),
                sep = "", collapse = "\n"),"\n")
    } else {
      cat(paste(names,
                sapply(where.measured, `[[`, i = "lat"), " N, ",
                sapply(where.measured, `[[`, i = "lon"), " E",
                sep = "", collapse = "\n"), "\n")
    }
  }
  if (class_spct(x)[1] %in% c("raw_spct", "cps_spct") &&
      getMultipleWl(x) == 1) {
    if (!all(is.na(getInstrDesc(x)))) {
      print(getInstrDesc(x))
      cat("\n")
    }
    if (!all(is.na(getInstrSettings(x)))) {
      print(getInstrSettings(x))
      cat("\n")
    }
  }
  if (class_spct(x)[1] == "filter_spct") {
    properties <- filter_properties(x, return.null = TRUE)
    if (!is.null(properties)) {
      print(properties)
      cat("\n")
    }
  }
  if (class_spct(x)[1] == "solute_spct") {
    properties <- solute_properties(x, return.null = TRUE)
    if (!is.null(properties)) {
      print(properties)
      cat("\n")
    }
  }
  if (is_scaled(x)) {
    scaling <- getScaling(x)
    if (all(is.na(scaling[["cols"]]))) {
      cat("Rescaled to '", scaling[["f"]], "' = ", scaling[["target"]],
          " for wavelengths in ", paste(scaling[["range"]], collapse = "-"), " nm\n", sep = "")
    } else {
      cat("Rescaled to '", scaling[["f"]], " of ",
          scaling[["cols"]], "' = ", scaling[["target"]],
          " for wavelengths in ", paste(scaling[["range"]], collapse = "-"), " nm\n", sep = "")
    }
  }
  if (is_normalized(x)) {
    norm <- getNormalized(x)
    normalization <- getNormalization(x)
    if (!anyNA(normalization[["norm.wl"]]) &&
        !anyNA(normalization[["norm.cols"]]) &&
        !anyNA(normalization[["norm.type"]])) {
      if (normalization[["norm.type"]][1] == "wavelength" ||
          all(is.na(normalization[["norm.range"]]))) {
        cat("Spectral data in ",
            paste(normalization[["norm.cols"]], collapse = ", "),
            " normalized to 1 at ",
            paste(round(normalization[["norm.wl"]], digits = 1),
                  collapse = " nm, "), " nm (",
            normalization[["norm.type"]][1], ")\n", sep = "")
      } else {
        cat("Spectral data in ",
            paste(normalization[["norm.cols"]], collapse = ", "),
            " normalized to 1 at ",
            paste(round(normalization[["norm.wl"]], digits = 1),
                  collapse = " nm, "), " nm (",
            normalization[["norm.type"]][1], " in ",
            paste(round(normalization[["norm.range"]], 2),
                  collapse = "-"),
            " nm)\n", sep = "")
      }
    } else {
      if (is.numeric(norm)) {
        cat("Spectral data normalized to 1 at ", norm, " nm \n", sep = "")
      } else {
        cat("Spectral data normalized to 1\n")
      }
    }
  }
  if (is_effective(x)) {
    BSWF <- getBSWFUsed(x)
    cat("Data weighted using '", BSWF, "' BSWF\n", sep = "")
  }
  var_labels <- make_var_labels(x)
  cat("Variables:\n", paste(names(var_labels), var_labels, sep = ": ", collapse = "\n "), "\n--\n")
  print(tibble::as_tibble(x), n = n, width = width)
  invisible(x)
}

# print method ------------------------------------------------------------

#' @describeIn print.generic_spct
#'
#' @param n.members	numeric Number of members of the collection to print.
#'
#' @export
#'
print.generic_mspct <- function(x,
                                ...,
                                attr.simplify = TRUE,
                                n = NULL,
                                width = NULL,
                                n.members = 10)  {
  cat("Object: ", class(x)[1], " ", dplyr::dim_desc(x), "\n", sep = "")
  member.names <- names(x)
  if (length(member.names) > n.members) {
    skipped.members <- length(member.names) - n.members
    member.names <- member.names[1:n.members]
  } else {
    skipped.members <- 0
  }
  for (name in member.names) {
    cat("--- Member:", name, "---\n")
    print(x[[name]], ..., attr.simplify = attr.simplify, n = n, width = width)
  }
  if (skipped.members > 0) {
    cat("..........................\n",
        skipped.members, " other member spectra not shown\n", sep = "")
  }
  cat("\n--- END ---")
  invisible(x)
}

# names of all spectral summary classes -----------------------------------

#' Function that returns a vector containing the names of spectral summary
#' classes.
#'
#' @export
#'
#' @return A \code{character} vector of class names.
#'
summary_spct_classes <- function() {
  c("summary_raw_spct", "summary_cps_spct",
    "summary_filter_spct", "summary_reflector_spct",
    "summary_source_spct", "summary_object_spct",
    "summary_response_spct", "summary_chroma_spct",
    "summary_solute_spct", "summary_generic_spct")
}

# is functions for spct summary classes --------------------------------------------

#' Query class of spectrum summary objects
#'
#' Functions to check if an object is of a given type of spectrum, or coerce it if
#' possible.
#'
#' @param x an R object.
#'
#' @return These functions return \code{TRUE} if its argument is a of the queried type
#'   of spectrum and \code{FALSE} otherwise.
#'
#' @note Derived types also return TRUE for a query for a base type such as
#' \code{generic_spct}.
#'
#' @export
#'
#' @examples
#' sm <- summary(sun.spct)
#' is.summary_source_spct(sm)
#'
is.summary_generic_spct <- function(x) inherits(x, "summary_generic_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_raw_spct <- function(x) inherits(x, "summary_raw_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_cps_spct <- function(x) inherits(x, "summary_cps_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_source_spct <- function(x) inherits(x, "summary_source_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_response_spct <- function(x) inherits(x, "summary_response_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_filter_spct <- function(x) inherits(x, "summary_filter_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_reflector_spct <- function(x) inherits(x, "summary_reflector_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_object_spct <- function(x) inherits(x, "summary_object_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_solute_spct <- function(x) inherits(x, "summary_solute_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_chroma_spct <- function(x) inherits(x, "summary_chroma_spct")

#' @rdname is.summary_generic_spct
#'
#' @export
#'
is.any_summary_spct <- function(x) {
  inherits(x, summary_spct_classes())
}

# summary -----------------------------------------------------------------

#' Summary of one or more spectra
#'
#' Methods of generic function summary for objects of spectral classes and
#' of classes for collections of spectra.
#'
#' @param object An object of one of the spectral classes for which a summary is
#'   desired.
#' @param maxsum integer Indicates how many levels should be shown for factors.
#' @param digits integer Used for number formatting with \code{\link{format}()}.
#' @param ... additional arguments affecting the summary produced, ignored in
#'   current version.
#' @param expand character One of \code{"none"}, \code{"collection"},
#'   \code{"each"} or \code{"auto"} indicating if multiple spectra in long form
#'   should be summarized as a collection or individually.
#'
#' @details
#' Objects are summarized as is, ignoring the current settings of R options
#' \code{photobiology.radiation.unit} and \code{photobiology.filter.qty}. Unlike
#' R's summary, these methods can optionally summarize each spectrum stored in
#' long form returning a list of summaries. Although this is frequently the most
#' informative approach, the default remains similar to \code{summary()} method
#' from R: to summarize \code{object} as a whole. Alternatively, multiple
#' spectra stored in long form, can optionally be summarized also as a
#' collection of spectra. Passing \code{"auto"} in the call, is equivalent to
#' passing \code{"each"} or \code{"collection"} depending on the number of
#' spectra contained in the object.
#'
#' @return A summary object matching the class of \code{object}, or a list of
#'   such objects or a summary object for a matching collection of spectra.
#'   Metadata stored in attributes are copied to identical attributes in the
#'   returned summary objects except when \code{object} is a collection
#'   of spectra or if \code{expand = "collection"} is passed in the call. In
#'   this two cases, a condensed summary is returned as a data frame and
#'   attributes from each member can be copied to variables in it.
#'
#' @seealso \code{\link[photobiology]{print.summary_generic_spct}}
#'
#' @export
#'
#' @method summary generic_spct
#'
#' @name summary.generic_spct
#'
#' @examples
#' summary(sun.spct)
#' class(summary(sun.spct))
#'
#' summary(two_filters.spct)
#' class(summary(two_filters.spct))
#'
#' summary(sun_evening.spct)
#' summary(two_filters.spct, expand = "none")
#' summary(two_filters.spct, expand = "each")
#' summary(two_filters.spct, expand = "collection")
#' summary(two_filters.spct, expand = "auto") # <= 4 spectra
#' summary(sun_evening.spct, expand = "auto") # > 4 spectra
#'
#' where_measured(sun.spct)
#' where_measured(summary(sun.spct))
#' what_measured(summary(two_filters.spct))
#' what_measured(summary(two_filters.spct, expand = "each")[[1]])
#'
#' summary(sun_evening.mspct)
#' summary(sun_evening.mspct, which.metadata = "when.measured")
#' summary(two_filters.mspct, which.metadata = "what.measured")
#' summary(two_filters.mspct, expand = "each")
#'
summary.generic_spct <- function(object,
                                 maxsum = 7,
                                 digits = max(3, getOption("digits") - 3),
                                 ...,
                                 expand = "none") {

  if (expand == "auto") {
    if (getMultipleWl(object) > 1 && getMultipleWl(object) <= 4L) {
      expand <- "each"
    } else {
      expand <- "collection"
    }
  }

  # optionally convert from long form into collection derived from generic_mspct
  if (expand %in% c("each", "collection")) {
    return(summary(subset2mspct(object),
                   maxsum = maxsum,
                   digits = digits,
                   expand = expand,
                   ...))
  }

  stopifnot(expand == "none")

  # summary of a generic_spct
  z <- list()
  class(z) <- c(paste("summary_", class_spct(object), sep = ""), class(z))

  object.name <- substitute(object)
  z[["orig.name"]] <- if (is.name(object.name)) as.character(object.name) else "anonymous"
  z[["orig.class"]] <- class_spct(object)[1]
  z[["orig.dim_desc"]] <- dplyr::dim_desc(object)
  z[["wl.range"]] <- wl_range(object)
  z[["wl.stepsize"]] <- wl_stepsize(object)
  z[["var.labels"]] <- make_var_labels(object)
  z[["summary"]] <- summary(as.data.frame(object), maxsum = maxsum, digits = digits, ...)

  z <- copy_attributes(object, z,
                       copy.class = FALSE,
                       which = c(all_spct_attr.ls[["generic_spct"]],
                                 all_spct_attr.ls[[class(object)[1]]]))
  z
}


# Print spectral summaries ------------------------------------------------

#' Print spectral summary
#'
#' A function to nicely print objects of classes "summary...spct".
#'
#' @param x An object of one of the summary classes for spectra
#' @param ... not used in current version
#'
#' @method print summary_generic_spct
#'
#' @name print.summary_generic_spct
#'
#' @export
#'
#' @examples
#' print(summary(sun.spct))
#'
print.summary_generic_spct <- function(x, ..., attr.simplify = TRUE) {
  # ensure backwards compatibility with summary objects created by versions < 0.9.30
  if (!exists("orig.name", x) || is.na(x[["orig.name"]])) {
    x[["orig.name"]] <- "'unknown name'"
  }
  cat("Summary of ", x[["orig.class"]], " ", x[["orig.dim_desc"]], " object: ", x[["orig.name"]] ,"\n", sep = "")
  m.wl <- getMultipleWl(x)
  if (m.wl > 1) {
    cat("containing ", m.wl, " spectra in long form\n")
  }
  cat("Wavelength range ",
      paste(signif(x[["wl.range"]], 8), sep = "", collapse = "-"), " nm, step ",
      paste(unique(signif(x[["wl.stepsize"]], 7)), sep = "", collapse = "-"),
      " nm\n", sep = "")
  what.measured <- getWhatMeasured(x, simplify = attr.simplify)
  if (!any(is.na(what.measured))) {
    if (length(what.measured) > 1) {
      if (!is.list(what.measured)) {
        what.measured <- list(what.measured)
      }
      names <- names(what.measured)
      if (is.null(names)) {
        names <- "Label: "
      } else {
        names <- paste(names, "label: ")
      }
      cat(paste(names,
                what.measured,
                sep = "", collapse = "\n"), "\n")
    } else {
      cat("Label: ", what.measured, "\n")
    }
  }
  when.measured <- getWhenMeasured(x, simplify = attr.simplify)
  if (!any(is.na(when.measured))) {
    if (length(when.measured) > 1) {
      if (!is.list(when.measured)) {
        when.measured <- list(when.measured)
      }
      names <- names(when.measured)
      if (is.null(names)) {
        names <- "Measured on "
      } else {
        names <- paste(names, "measured on ")
      }
      cat(paste(names,
                sapply(when.measured, as.character), " UTC",
                sep = "", collapse = "\n"), "\n")
    } else {
      cat("Measured on ", as.character(when.measured), " UTC\n")
    }
  }
  where.measured <- getWhereMeasured(x, simplify = attr.simplify)
  if (!any(is.na(where.measured))) {
    if (is.data.frame(where.measured)) {
      where.measured <- list(where.measured)
    }
    names <- names(where.measured)
    if (is.null(names)) {
      names <- "Measured at "
    } else {
      names <- paste(names, "measured at ")
    }
    if (all(sapply(where.measured,
                   function(x) {"address" %in% names(x) && !any(is.na(x[["address"]]))}))) {
      cat(paste(names,
                sapply(where.measured, `[[`, i = "lat"), " N, ",
                sapply(where.measured, `[[`, i = "lon"), " E; ",
                stringr::str_trunc(sapply(where.measured, `[[`, i = "address"), width = 30, side = "right"),
                sep = "", collapse = "\n"),"\n")
    } else {
      cat(paste(names,
                sapply(where.measured, `[[`, i = "lat"), " N, ",
                sapply(where.measured, `[[`, i = "lon"), " E",
                sep = "", collapse = "\n"), "\n")
    }
  }
  if (class(x)[1] %in% c("summary_raw_spct", "summary_cps_spct") &&
      getMultipleWl(x) == 1) {
    print(getInstrDesc(x))
    cat("\n")
    print(getInstrSettings(x))
  }
  if (class(x)[1] == "summary_solute_spct") {
    properties <- solute_properties(x, return.null = TRUE)
    if (!is.null(properties)) {
      print(properties)
      cat("\n")
    }
  }
  if (is_scaled(x)) {
    scaling <- getScaling(x)
    if (all(is.na(scaling[["cols"]]))) {
      cat("Rescaled to '", scaling[["f"]], "' = ", scaling[["target"]],
          " for wavelengths in ", paste(scaling[["range"]], collapse = "-"), " nm\n", sep = "")
    } else {
      cat("Rescaled to '", scaling[["f"]], " of ",
          scaling[["cols"]], "' = ", scaling[["target"]],
          " for wavelengths in ", paste(scaling[["range"]], collapse = "-"), " nm\n", sep = "")
    }
  }
  if (is_normalized(x)) {
    norm <- getNormalized(x)
    normalization <- getNormalization(x)
    if (!anyNA(normalization[["norm.wl"]]) &&
        !anyNA(normalization[["norm.cols"]]) &&
        !anyNA(normalization[["norm.type"]])) {
      if (normalization[["norm.type"]][1] == "wavelength" ||
          all(is.na(normalization[["norm.range"]]))) {
        cat("Spectral data in ",
            paste(normalization[["norm.cols"]], collapse = ", "),
            " normalized to 1 at ",
            paste(round(normalization[["norm.wl"]], digits = 1),
                  collapse = " nm, "), " nm (",
            normalization[["norm.type"]][1], ")\n", sep = "")
      } else {
        cat("Spectral data in ",
            paste(normalization[["norm.cols"]], collapse = ", "),
            " normalized to 1 at ",
            paste(round(normalization[["norm.wl"]], digits = 1),
                  collapse = " nm, "), " nm (",
            normalization[["norm.type"]][1], " in ",
            paste(round(normalization[["norm.range"]], 2), collapse = "-"),
            " nm)\n", sep = "")
      }
    } else {
      if (is.numeric(norm)) {
        cat("Spectral data normalized to 1 at ",
            round(norm, digits = 1), " nm \n", sep = "")
      } else {
        cat("Spectral data normalized to 1\n")
      }
    }
  }
  if (is_effective(x)) {
    BSWF <- getBSWFUsed(x)
    cat("Data weighted using '", BSWF, "' BSWF\n", sep = "")
  }
  if (exists("var.labels", x)) {
    cat("Variables:\n",
        paste(names(x[["var.labels"]]), x[["var.labels"]],
              sep = ": ", collapse = "\n "), "\n--\n")
  } else {
    cat("\n--\n")
  }
  print(x[["summary"]])
  invisible(x)
}

# summary -----------------------------------------------------------------

#' @describeIn summary.generic_spct
#'
#' @param idx character Name of the column with the names of the members of the
#' collection of spectra.
#' @param which.metadata character vector Names of attributes to retrieve, or
#'   "none" or "all". Obeyed if \code{expand = "collection"}, its default.
#'
#' @export
#'
#' @method summary generic_mspct
#'
summary.generic_mspct <- function(object,
                                  maxsum = 7,
                                  digits = max(3, getOption("digits") - 3),
                                  idx = "spct.idx",
                                  which.metadata = NULL,
                                  expand = "none",
                                  ...) {

  if (expand == "each") {
    return(lapply(subset2mspct(object), #
                  FUN = summary,
                  maxsum = maxsum,
                  digits = digits,
                  expand = "none", # avoid endless recursion
                  ...))
  }

  stopifnot(expand %in% c("none", "collection"))

  if (is.null(which.metadata)) {
    which.metadata <- switch(class(object)[1],
                             filter_mspct = c("multiple.wl", "Tfr.type"),
                             reflector_mspct = c("multiple.wl", "Rfr.type"),
                             object_mspct = c("multiple.wl", "Tfr.type", "Rfr.type"),
                             source_mspct = c("multiple.wl", "time.unit"),
                             generic_mspct = "multiple.wl",
                             "multiple.wl")
  }
  object.name <- substitute(object)
  z <- list()
  z[["orig.name"]] <- if (is.name(object.name)) as.character(object.name) else "anonymous"
  z[["orig.class"]] <- class(object)[1]
  z[["orig.dim_desc"]] <- dplyr::dim_desc(object)

  summary.tb <- tibble::tibble(.rows = length(object))
  summary.tb[[idx]] <- names(object)
  summary.tb[["class"]] <- sapply(object, function(x) {class(x)[1]})
  summary.tb[["dim"]] <- sapply(object, dplyr::dim_desc)
  summary.tb[["w.length.min"]] <- sapply(object, wl_min)
  summary.tb[["w.length.max"]] <- sapply(object, wl_max)
  summary.tb[["colnames"]] <- unname(lapply(object, colnames))

  if (length(which.metadata) != 0L && which.metadata[1] != "none") {
    if (length(which.metadata) == 1L && which.metadata[1] == "all") {
     which.metadata <-
       setdiff(spct_attributes(class(object)[1]),
               c("spct.tags", "spct.version", "comment",
                 "straylight.corrected",
                 "slit.corrected",
                 "QC_dark_pass",
                 "idfactor",
                 "spct.idx",
                 "normalization")
               )
    }
    summary.tb <-
      add_attr2tb(summary.tb, object, col.names = which.metadata)
  }

  z[["summary"]] <- summary.tb

  class(z) <- c(paste("summary_",
                      grep("_mspct$", class(object), value = TRUE),
                      sep = ""), class(z))

  comment(z) <- paste("Comment from '", z[["object.name"]], "'.\n",
                      comment(object), sep = "")
  z
}

#' @describeIn print.summary_generic_spct
#'
#' @param width integer Width of text output to generate. This defaults to NULL,
#'   which means use the width option.
#' @param ... named arguments passed to the \code{print()} method for class
#'   \code{"tbl_df"}.
#' @param attr.simplify logical If all members share the same attribute value
#'   return one copy instead of a data.frame, list or vector.
#' @param n integer Number of member spectra for which information is printed.
#'
#' @seealso \code{\link[tibble]{formatting}}
#'
#' @method print summary_generic_mspct
#'
#' @export
#'
#' @examples
#' print(summary(sun_evening.mspct))
#'
print.summary_generic_mspct <- function(x,
                                        width = NULL,
                                        ...,
                                        n = NULL) {
  cat("Summary of ", x[["orig.class"]], " ", x[["orig.dim_desc"]],
      " object: ", x[["orig.name"]] ,"\n", sep = "")
  print(x[["summary"]], width = width, ..., n = n)
}

# Print properties ---------------------------------------------------------
#
# These methods are called when printing spectra and their summaries.

#' Print methods for metadata records
#'
#' Print methods for objects of classes used to store different meta data
#' properties in the classes for different types of spectra.
#'
#' @param x An object of one of the summary classes for spectra.
#' @param ... not used in current version.
#'
#' @details These methods print an abbreviated representaion of objects used
#' to store metadata in attributes. They are similar to \emph{records} and
#' formatted printing is useful both on its own and in the print methods for
#' spectra and their summaries.
#'
#' @export
#'
#' @examples
#'
#' print(getInstrDesc(sun_evening.spct))
#' str(getInstrDesc(sun_evening.spct))
#'
#' print(getInstrSettings(sun_evening.spct))
#' str(getInstrSettings(sun_evening.spct))
#'
#' print(filter_properties(polyester.spct))
#' str(filter_properties(polyester.spct))
#'
#' print(solute_properties(phenylalanine.spct))
#' str(solute_properties(phenylalanine.spct))
#'
#' @name print.metadata
#'
print.instr_desc <- function(x, ...) {
  if (!length(x) || !is.list(x)) {
    warning("'x' is not an instrument descriptor")
  } else {
    if (!("entrance.optics" %in% names(x)) ||
        all(is.na(x[["entrance.optics"]]))) {
      diffuser <- "unknown"
    } else {
      diffuser <- x[["entrance.optics"]][["geometry"]]
    }
    cat("Data acquired with '",
        x[["spectrometer.name"]], "' s.n. ", x[["spectrometer.sn"]],
        "\ngrating '", x[["bench.grating"]],
        "', slit '", x[["bench.slit"]], "'",
        "\ndiffuser '", diffuser, "'",
        ...,
        sep = ""
    )
  }
  invisible(x)
}

#' @rdname print.metadata
#'
#' @export
#'
print.instr_settings <- function(x, ...) {
  if (!length(x) || !is.list(x)) {
    warning("'x' is not an instrument settings record")
  } else {
    cat("integ. time (s): ",
        paste(signif(as.numeric(x[["integ.time"]]) * 1e-6, digits = 3),
              collapse = ", "),
        "\ntotal time (s): ",
        paste(signif(as.numeric(x[["tot.time"]]) * 1e-6, digits = 3),
              collapse = ", "),
        "\ncounts @ peak (% of max): ",
        signif(as.numeric(x[["rel.signal"]]) * 100, digits = 3),
        sep = "",
        ...
    )
  }
  invisible(x)
}

#' @rdname print.metadata
#'
#' @export
#'
print.filter_properties <- function(x, ...) {
  if (!length(x) || !is.list(x)) {
    warning("'x' is not a filter properties record")
  } else {
    cat("Rfr (/1): ",
        paste(sprintf("%#.3f", x[["Rfr.constant"]]), collapse = " + "), ", ",
        "thickness (mm): ",
        paste(sprintf("%#.3g", x[["thickness"]] * 1e3), collapse = " + "), ", ",
        "attenuation mode: ", x[["attenuation.mode"]], ".",
        sep = "",
        ...
    )
  }
  invisible(x)
}

#' @rdname print.metadata
#'
#' @export
#'
print.solute_properties <- function(x, ...) {
  if (!length(x) || !is.list(x)) {
    warning("'x' is not a solute properties record")
  } else {
    cat("Name: ", x[["name"]][1], ", ",
      "Molar mass (Da): ", round(x[["mass"]], digits = 2), ", ",
      "Formula: ", x[["formula"]][1], ".",
      sep = "",
      ...
    )
  }
  invisible(x)
}
