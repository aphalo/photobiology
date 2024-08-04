# summary -----------------------------------------------------------------

#' Summary of a collection of spectra
#'
#' Method of generic function summary for objects of spectral collection
#' classes.
#'
#' @param object An object of one of the spectral collection classes for which a
#'   summary is desired
#' @param maxsum integer Indicates how many levels should be shown for factors.
#' @param digits integer Used for number formatting with \code{\link{format}()}.
#' @param idx character Name of the column with the names of the members of the
#' collection of spectra.
#' @param which.metadata character vector Names of attributes to retrieve, or
#'   "none" or "all".
#' @param ... additional arguments affecting the summary produced, ignored in
#'   current version
#'
#' @return A summary object matching the class of \code{object}.
#'
#' @export
#'
#' @method summary generic_mspct
#'
#' @name summary
#'
#' @examples
#' summary(sun_evening.mspct)
#' summary(sun_evening.mspct, which.metadata = "when.measured")
#' summary(two_filters.mspct, which.metadata = "what.measured")
#'
summary.generic_mspct <- function(object,
                                  maxsum = 7,
                                  digits = max(3, getOption("digits") - 3),
                                  idx = "spct.idx",
                                  which.metadata = NULL,
                                  ...) {
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

  if (length(which.metadata) != 1L || which.metadata != "none") {
    if (length(which.metadata) == 1L && which.metadata == "all") {
      which.metadata <- c("-", "names", "row.names", "spct.tags", "spct.version", "comment")
    }
    metadata.tb <- msdply(object, spct_attr2tb, which = which.metadata, idx = idx)
    names(metadata.tb) <- gsub("spct_attr2tb_", "", names(metadata.tb))
    summary.tb <- dplyr::left_join(summary.tb, metadata.tb, by = idx)
  }

  z[["summary"]] <- summary.tb

  class(z) <- c(paste("summary_",
                      grep("_mspct$", class(object), value = TRUE),
                      sep = ""), class(z))

  comment(z) <- paste("Comment from '", z[["object.name"]], "'.\n",
                      comment(object), sep = "")
  z
}

#' Print spectral collection summary
#'
#' A function to nicely print objects of classes "summary...mspct".
#'
#' @param x An object of one of the summary classes for collections of spectra.
#' @param width integer Width of text output to generate. This defaults to NULL,
#'   which means use the width option.
#' @param ... named arguments passed to the \code{print()} method for class
#'   \code{"tbl_df"}.
#' @param n integer Number of member spectra for which information is printed.
#'
#' @seealso \code{\link[tibble]{formatting}}
#'
#' @method print summary_generic_mspct
#'
#' @name print
#'
#' @export
#'
#' @examples
#' print(summary(sun_evening.mspct))
#'
print.summary_generic_mspct <- function(x, width = NULL, ..., n = NULL) {
  cat("Summary of ", x[["orig.class"]], " ", x[["orig.dim_desc"]], " object: ", x[["orig.name"]] ,"\n", sep = "")
  print(x[["summary"]], width = width, ..., n = n)
}
