# names of all multi spectral classes -------------------------------------------

#' Names of multi-spectra classes
#'
#' Function that returns a vector containing the names of multi-spectra classes
#' using for collections of spectra.
#'
#' @export
#'
#' @return A \code{character} vector of class names.
#'
#' @examples
#' mspct_classes()
#'
mspct_classes <- function() {
  c("calibration_mspct",
    "raw_mspct", "cps_mspct",
    "filter_mspct", "reflector_mspct",
    "source_mspct", "object_mspct",
    "response_mspct", "chroma_mspct", "generic_mspct")
}

# remove mspct class attributes --------------------------------------------

#' Remove "generic_mspct" and derived class attributes.
#'
#' Removes from an spectrum object the class attributes "generic_mspct" and any
#' derived class attribute such as "source_mspct". \strong{This operation is done
#' by reference!}
#'
#' @param x an R object.
#' @export
#'
#' @note If \code{x} is an object of any of the multi spectral classes defined
#'   in this package, this function changes by reference the multi spectrum
#'   object into the underlying list object. Otherwise, it just leaves \code{x}
#'   unchanged. The modified \code{x} is also returned invisibly.
#'
#' @return A character vector containing the removed class attribute values.
#'   This is different to the behaviour of function \code{unlist} in base R!
#'
#' @family set and unset 'multi spectral' class functions
#'
rmDerivedMspct <- function(x) {
  name <- substitute(x)
  mspctclasses <- mspct_classes()
  allclasses <- class(x)
  attr(x, "mspct.dim") <- NULL
  attr(x, "mspct.byrow") <- NULL
  attr(x, "mspct.version") <- NULL
  class(x) <- setdiff(allclasses, mspctclasses)
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(setdiff(allclasses, class(x)))
}


# query member classes ----------------------------------------------------

#' Classes common to all collection members.
#'
#' Finds the set intersection among the class attributes of all collection
#' member as a target set of class names.
#'
#' @param l a list or a generic_mscpt object or of a derived class.
#' @param target.set character The target set of classes within which to search
#'   for classes common to all members.
#' @export
#'
#' @return A character vector containing the class attribute values.
#'
#' @family set and unset 'multi spectral' class functions
#'
shared_member_class <- function(l, target.set = spct_classes()) {
  l.class <- target.set
  for (i in seq_along(l)) {
    member_class <- class(l[[i]])
    l.class <- intersect(l.class, member_class)
  }
  l.class
}

# Constructors ------------------------------------------------------------

#' @title Collection-of-spectra constructor
#'
#' @description Converts a list of spectral objects into a "multi spectrum"
#'   object by setting the class attribute of the list of spectra to the
#'   corresponding multi-spct class, check that components of the list belong to
#'   the expected class.
#'
#' @param l list of generic_spct or derived classes
#' @param class character The multi spectrum object class or the expected class
#'   for the elements of l
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#' @param dim integer vector of dimensions
#' @param ... ignored
#'
#' @export
#' @exportClass generic_mspct
#'
#' @note Setting class = source_spct or class = source_mspct makes no difference
#'
#' @family collections of spectra classes family
#' @examples
#' filter_mspct(list(polyester.spct, yellow_gel.spct))
#'
generic_mspct <- function(l = NULL, class = "generic_spct",
                          ncol = 1, byrow = FALSE,
                          dim = c(length(l) %/% ncol, ncol)) {
  if (is.generic_spct(l)) {
    l <- list(l)
  }
  if (is.null(l)) {
    l <- list()
  }
  stopifnot(is.list(l))

  class <- class[1]
  if (class %in% mspct_classes()) {
    multi_class <- class
    spct_class <- paste(sub("_mspct", "_spct", class))
  } else if (class %in% spct_classes()) {
    multi_class <- paste(sub("_spct", "_mspct", class))
    spct_class <- class
  } else {
    stop("'class' argument '", class, "' is not recognized as a spectral class")
  }

  if (class(l)[1] != multi_class) {
    if (is.any_mspct(l)) {
      rmDerivedMspct(l)
    }
    for (spct in l) {
      stopifnot(spct_class %in% class_spct(spct))
    }
    if (multi_class != "generic_mspct") {
      multi_class <- c(multi_class, "generic_mspct")
    }
    multi_class <- c(multi_class, class(l))
    class(l) <- multi_class
  }
  if (length(l) > 0 && is.null(names(l))) {
    attr(l, "names") <- paste("spct", seq_along(l), sep = "_")
  }
  attr(l, "mspct.version") <- 2

  dim(l) <- dim
  attr(l, "mspct.byrow") <- as.logical(byrow)
  l
}

#' @describeIn generic_mspct Specialization for collections of \code{calibration_spct} objects.
#'
#' @export
#' @exportClass calibration_mspct
#'
raw_mspct <- function(l = NULL, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "calibration_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{raw_spct} objects.
#'
#' @export
#' @exportClass raw_mspct
#'
raw_mspct <- function(l = NULL, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "raw_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{cps_spct} objects.
#'
#' @export
#' @exportClass cps_mspct
#'
cps_mspct <- function(l = NULL, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "cps_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{source_spct} objects.
#'
#' @export
#' @exportClass source_mspct
#'
source_mspct <- function(l = NULL, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "source_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{filter_spct} objects.
#'
#' @export
#' @exportClass filter_mspct
#'
filter_mspct <- function(l = NULL, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "filter_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{reflector_spct} objects.
#'
#' @export
#' @exportClass reflector_mspct
#'
reflector_mspct <- function(l = NULL, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "reflector_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{object_spct} objects.
#'
#' @export
#' @exportClass object_mspct
#'
object_mspct <- function(l = NULL, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "object_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{response_spct} objects.
#'
#' @export
#' @exportClass response_mspct
#'
response_mspct <- function(l = NULL, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "response_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{chroma_spct} objects.
#'
#' @export
#' @exportClass chroma_mspct
#'
chroma_mspct <- function(l = NULL, ncol = 1, byrow = FALSE, ...) {
  generic_mspct(l, class = "chroma_spct", ncol = ncol, byrow = byrow)
}

# is functions for mspct classes --------------------------------------------

#' Query class of spectrum objects
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
#' \code{generic_mspct}.
#'
#' @export
#' @rdname is.generic_mspct
#' @examples
#' my.mspct <- filter_mspct(list(polyester.spct, yellow_gel.spct))
#' is.any_mspct(my.mspct)
#' is.filter_mspct(my.mspct)
#' is.source_mspct(my.mspct)
#'
is.generic_mspct <- function(x) inherits(x, "generic_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.calibration_mspct <- function(x) inherits(x, "calibration_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.raw_mspct <- function(x) inherits(x, "raw_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.cps_mspct <- function(x) inherits(x, "cps_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.source_mspct <- function(x) inherits(x, "source_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.response_mspct <- function(x) inherits(x, "response_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.filter_mspct <- function(x) inherits(x, "filter_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.reflector_mspct <- function(x) inherits(x, "reflector_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.object_mspct <- function(x) inherits(x, "object_mspct")

#' @rdname is.generic_mspct
#' @export
#'
is.chroma_mspct <- function(x) inherits(x, "chroma_mspct")

#' @rdname is.generic_mspct
#'
#' @export
#'
is.any_mspct <- function(x) {
  inherits(x, "generic_mspct")
}

# as functions for mspct classes --------------------------------------------

#' @title Collection-of-spectra copy-constructor
#'
#' @description Return a copy of an R object with its class set to a given type
#'   of spectrum.
#'
#' @param x a list of spectral objects or a list of objects such as data frames
#'   that can be converted into spectral objects.
#' @param force.spct.class logical indicating whether to change the class
#'   of members to \code{generic_spct} or retain the existing class.
#'
#' @return These functions return a copy of \code{x} converted into a given
#'   class of spectral collection object, if \code{x} is a valid argument to the
#'   corresponding set function.
#'
#' @note Members of \code{generic_mspct} objects can be heterogeneous: they can
#'   belong any class derived from \code{generic_spct} and class is not
#'   enforced. In this case when \code{x} is a list of data frames,
#'   \code{force.spct.class = TRUE} needs to be supplied.
#'
#' @export
#'
#' @family conversion of collections of spectra
#'
#' @rdname as.generic_mspct
#'
as.generic_mspct <- function(x, force.spct.class = FALSE) {
  y <- x
  rmDerivedMspct(y)
  if (force.spct.class) {
    y <- plyr::llply(y, setGenericSpct)
  }
  generic_mspct(y)
}

#' @rdname as.generic_mspct
#'
#' @export
#'
as.calibration_mspct <- function(x) {
  y <- x
  rmDerivedMspct(y)
  z <- plyr::llply(y, setCalibrationSpct)
  raw_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @export
#'
as.raw_mspct <- function(x) {
  y <- x
  rmDerivedMspct(y)
  z <- plyr::llply(y, setRawSpct)
  raw_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @export
#'
as.cps_mspct <- function(x) {
  y <- x
  rmDerivedMspct(y)
  z <- plyr::llply(y, setCpsSpct)
  cps_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @param time.unit character A string, "second", "day" or "exposure"
#' @param bswf.used character
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning
#'
#' @export
#'
as.source_mspct <- function(x,
                            time.unit=c("second", "day", "exposure"),
                            bswf.used=c("none", "unknown"),
                            strict.range = getOption("photobiology.strict.range", default = FALSE)) {
  y <- x
  rmDerivedMspct(y)
  z <- plyr::llply(y, setSourceSpct, time.unit = time.unit,
                   strict.range = strict.range, bswf.used = bswf.used)
  source_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @export
#'
as.response_mspct <- function(x, time.unit = "second") {
  y <- x
  rmDerivedMspct(y)
  z <- plyr::llply(y, setResponseSpct, time.unit = time.unit)
  response_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @param Tfr.type a character string, either "total" or "internal"
#'
#' @export
#'
as.filter_mspct <- function(x,
                            Tfr.type=c("total", "internal"),
                            strict.range = TRUE) {
  y <- x
  rmDerivedMspct(y)
  z <- plyr::llply(y, setFilterSpct, Tfr.type = Tfr.type,
                   strict.range = strict.range)
  filter_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @param Rfr.type a character string, either "total" or "specular"
#'
#' @export
#'
as.reflector_mspct <- function(x,
                               Rfr.type = c("total", "specular"),
                               strict.range = TRUE) {
  y <- x
  rmDerivedMspct(y)
  z <- plyr::llply(y, setReflectorSpct, Rfr.type = Rfr.type,
                   strict.range = strict.range)
  reflector_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @export
#'
as.object_mspct <- function(x,
                            Tfr.type=c("total", "internal"),
                            Rfr.type=c("total", "specular"),
                            strict.range = TRUE) {
  y <- x
  rmDerivedMspct(y)
  z <- plyr::llply(y, setObjectSpct, Tfr.type = Tfr.type, Rfr.type = Rfr.type,
                   strict.range = strict.range)
  object_mspct(z)
}

#' @rdname as.generic_mspct
#'
#' @export
#'
as.chroma_mspct <- function(x) {
  y <- x
  rmDerivedMspct(y)
  z <- plyr::llply(y, setChromaSpct)
  chroma_mspct(z)
}

# constructor methods for data frames --------------------------------------

#' @title Convert a 'wide' or untidy data frame into a collection of spectra
#'
#' @description Convert a data frame object into a "multi spectrum" object by
#'   constructing a an object of a multi-spct class, converting numeric columns
#'   other than wavelength into individual spct objects.
#'
#' @param x data frame
#' @param member.class character Class of the collection members
#' @param spct.data.var character Name of the spectral data argument in the
#'   object constructor for \code{member.class}
#' @param w.length.var character Name of column containing wavelength data in
#'   nanometres
#' @param idx.var character Name of column containing data to be copied unchanged
#'   to each spct object
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#' @param ... additional named arguments passed to the member constructor
#'   function.
#'
#' @export
#'
#' @family conversion of collections of spectra
#'
split2mspct <- function(x,
                        member.class = NULL,
                        spct.data.var = NULL,
                        w.length.var = "w.length", idx.var = NULL,
                        ncol = 1, byrow = FALSE, ...) {
  stopifnot(!is.null(member.class) || !is.character(member.class))
  stopifnot(!is.null(spct.data.var) || !is.character(spct.data.var))
  if (is.generic_spct(x) && getMultipleWl(x) != 1) {
    stop("'split2mspct()' is for slicing vertically wide data in data frames ",
         "'subset2mspct()' is used in the case of tidy data in long form.")
  }
  if (!is.numeric(x[[w.length.var]])) {
    stop("Non-numeric variable '", w.length.var, "' is bad for wavelengths.")
  }
  collection.class <- sub("_spct", "_mspct", member.class, fixed = TRUE)
  member.constr <- member.class
  collection.constr <- collection.class
  col_names <- names(x)
  data.cols <- setdiff(col_names, c(w.length.var, idx.var))
  l <- list()
  for (col in data.cols) {
    if (!is.numeric(x[[col]])) {
      warning("Skipping non-numeric column in x: ", col)
      next
    }
    args <- list(w.length = x[[w.length.var]])
    args[[spct.data.var]] <- x[[col]]
    args.ellipsis <- list(...)
    l[[col]] <- do.call(member.constr, c(args, args.ellipsis))
    if (!is.null(idx.var)) {
      l[[col]][[idx.var]] <- x[[idx.var]]
    }
  }
  margs <- list(l = l, ncol = ncol, byrow = byrow)
  do.call(collection.constr, margs)
}

#' @rdname split2mspct
#' @export
#'
split2source_mspct <- function(x,
                               spct.data.var = "s.e.irrad",
                               w.length.var = "w.length", idx.var = NULL,
                               ncol = 1, byrow = FALSE, ...) {
  split2mspct(x = x,
              member.class = "source_spct",
              spct.data.var = spct.data.var,
              w.length.var = w.length.var,
              idx.var = idx.var,
              ncol = ncol, byrow = byrow,
              ...)
}

#' @rdname split2mspct
#' @export
#'
split2response_mspct <- function(x,
                                 spct.data.var = "s.e.response",
                                 w.length.var = "w.length", idx.var = NULL,
                                 ncol = 1, byrow = FALSE, ...) {
  split2mspct(x = x,
              member.class = "response_spct",
              spct.data.var = spct.data.var,
              w.length.var = w.length.var,
              idx.var = idx.var,
              ncol = ncol, byrow = byrow,
              ...)
}

#' @rdname split2mspct
#' @export
#'
split2filter_mspct <- function(x,
                               spct.data.var = "Tfr",
                               w.length.var = "w.length", idx.var = NULL,
                               ncol = 1, byrow = FALSE, ...) {
  split2mspct(x = x,
              member.class = "filter_spct",
              spct.data.var = spct.data.var,
              w.length.var = w.length.var,
              idx.var = idx.var,
              ncol = ncol, byrow = byrow,
              ...)
}

#' @rdname split2mspct
#' @export
#'
split2reflector_mspct <- function(x,
                                  spct.data.var = "Rfr",
                                  w.length.var = "w.length", idx.var = NULL,
                                  ncol = 1, byrow = FALSE, ...) {
  split2mspct(x = x,
              member.class = "reflector_spct",
              spct.data.var = spct.data.var,
              w.length.var = w.length.var,
              idx.var = idx.var,
              ncol = ncol, byrow = byrow,
              ...)
}

#' @rdname split2mspct
#' @export
#'
split2cps_mspct <- function(x,
                            spct.data.var = "cps",
                            w.length.var = "w.length", idx.var = NULL,
                            ncol = 1, byrow = FALSE, ...) {
  split2mspct(x = x,
              member.class = "cps_spct",
              spct.data.var = spct.data.var,
              w.length.var = w.length.var,
              idx.var = idx.var,
              ncol = ncol, byrow = byrow,
              ...)
}

#' @rdname split2mspct
#' @export
#'
split2raw_mspct <- function(x,
                            spct.data.var = "count",
                            w.length.var = "w.length", idx.var = NULL,
                            ncol = 1, byrow = FALSE, ...) {
  split2mspct(x = x,
              member.class = "raw_spct",
              spct.data.var = spct.data.var,
              w.length.var = w.length.var,
              idx.var = idx.var,
              ncol = ncol, byrow = byrow,
              ...)
}

#' @rdname split2mspct
#' @export
#'
split2calibration_mspct <- function(x,
                            spct.data.var = "irrad.mult",
                            w.length.var = "w.length", idx.var = NULL,
                            ncol = 1, byrow = FALSE, ...) {
  split2mspct(x = x,
              member.class = "calibration_spct",
              spct.data.var = spct.data.var,
              w.length.var = w.length.var,
              idx.var = idx.var,
              ncol = ncol, byrow = byrow,
              ...)
}

#' @title Convert 'long' or tidy spectral data into a collection of spectra
#'
#' @description Convert a data frame object or spectral object into a collection
#'   of spectra object of the corresponding class. For data frames converting
#'   numeric columns other than wavelength into individual spct objects.
#'
#' @param x a generic_spct object or a derived class, or a data frame
#' @param member.class character string
#' @param idx.var character Name of column containing data to be copied
#'   unchanged to each spct object
#' @param drop.idx logical Flag indicating whether to drop or keep idx.var in
#'   the collection members.
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#' @param ... additional named arguments passed to the member constructor
#'   function.
#'
#' @note A non-null value for \code{member.class} is mandatory only when
#'   \code{x} is a data frame.
#'
#' @return A collection of spectral objects, each with attributes set if x is a
#'   spectral object in long form with metadata attributes. If this object
#'   was created by row binding with 'photobiology' 0.9.14 or later then
#'   all metadata for each individual spectrum will be preserved, except for
#'   comments which are merged.
#'
#' @export
#'
#' @family conversion of collections of spectra
#'
subset2mspct <- function(x,
                         member.class = NULL,
                         idx.var = attr(x, "idfactor"),
                         drop.idx = TRUE,
                         ncol = 1, byrow = FALSE, ...) {
  stopifnot(is.data.frame(x))
  if (is.generic_spct(x) && is.null(member.class)) {
    member.class <- class(x)[1]
  }
  stopifnot(is.character(member.class))
  if (is.null(idx.var)) {
    idx.var <- "spct.idx"
  }
  stopifnot(idx.var %in% names(x))
  collection.class <- sub("_spct", "_mspct", member.class, fixed = TRUE)
  stopifnot(collection.class %in% mspct_classes())
  member.constr <- paste("as", member.class, sep = ".")
  collection.constr <- collection.class
  if (is.factor(x[[idx.var]])) {
    groups <- levels(x[[idx.var]])
    idx <- idx.var
  } else {
    # would hang or slowdown to a crawl if indexing by dates
    # could try benchmarking with as.numeric() to see how much faster it is
    if (lubridate::is.instant(x[[idx.var]])) {
    x[["tmp.idx"]] <- as.character(x[[idx.var]], tz = "UTC")
    idx <- "tmp.idx"
    } else {
      idx <- idx.var
    }
    groups <- unique(x[[idx]])
  }
  l <- list()
  for (grp in groups) {
    slice <- subset(x, x[[idx]] == grp)
    if (drop.idx) {
      slice[[idx.var]] <- NULL
    }
    if (idx != idx.var) {
      slice[[idx]] <- NULL
    }
    args <- list(x = slice)
    args.ellipsis <- list(...)
    l[[grp]] <- do.call(member.constr, c(args, args.ellipsis))
  }
  margs <- list(l = l, ncol = ncol, byrow = byrow)
  z <- do.call(collection.constr, margs)
  # copy metadata
  comment <- comment(x)
  if (!is.null(comment)) {
    z <- msmsply(z, `comment<-`, value = comment)
  }
  if (!is.generic_spct(x)) {
    return(z)
  }
  if (is_scaled(x)) {
    z <- msmsply(z, setScaled, scaled = TRUE)
  }
  if (is_normalized(x)) {
     z <- msmsply(z, setNormalized, norm = TRUE)
  }
  if (member.class == "source_spct" && is_effective(x)) {
    bswf.used <- getBSWFUsed(x)
    z <- msmsply(z, setBSWFUsed, bswf.used = bswf.used)
  }
  if (member.class %in% c("source_spct", "response_spct")) {
    time.unit <- getTimeUnit(x)
    z <- msmsply(z, setTimeUnit, time.unit = time.unit, override.ok = TRUE)
  }
  if (member.class %in% c("filter_spct", "object_spct")) {
    Tfr.type <- getTfrType(x)
    z <- msmsply(z, setTfrType, Tfr.type = Tfr.type)
  }
  if (member.class %in% c("reflector_spct", "object_spct")) {
    Rfr.type <- getRfrType(x)
    z <- msmsply(z, setRfrType, Rfr.type = Rfr.type)
  }
  # these methods return NA if attribute is not set
  when.measured <- getWhenMeasured(x)
  where.measured <- getWhereMeasured(x)
  what.measured <- getWhatMeasured(x)
  # these methods may return an empty list
  instr.desc <- getInstrDesc(x)
  instr.settings <- getInstrSettings(x)
  for (i in seq(along.with = z)) {
    if (!all(is.na(when.measured))) {
      if (is.list(when.measured) && length(when.measured) == length(groups)) {
        z[[i]] <- setWhenMeasured(z[[i]], when.measured[[i]])
      } else {
        z[[i]] <- setWhenMeasured(z[[i]], when.measured)
      }
    }
    if (!all(is.na(where.measured))) {
      if (is.list(where.measured) && !is.data.frame(where.measured) &&
            length(where.measured) == length(groups)) {
        z[[i]] <- setWhereMeasured(z[[i]], where.measured[[i]])
      } else {
        z[[i]] <- setWhereMeasured(z[[i]], where.measured)
      }
    }
    if (!all(is.na(what.measured))) {
      if (is.list(what.measured) && length(what.measured) == length(groups)) {
        z[[i]] <- setWhatMeasured(z[[i]], what.measured[[i]])
      } else {
        z[[i]] <- setWhatMeasured(z[[i]], what.measured)
      }
    }
    if (length(instr.desc) > 0) {
      if (is.list(instr.desc) &&
          !inherits(instr.desc, "instr_desc") &&
          length(instr.desc) == length(groups)) {
        z[[i]] <- setInstrDesc(z[[i]], instr.desc[[i]])
      } else {
        z[[i]] <- setInstrDesc(z[[i]], instr.desc)
      }
    }
    if (length(instr.settings) > 0) {
      if (is.list(instr.settings) &&
          !inherits(instr.settings, "instr_setting") &&
          length(instr.settings) == length(groups)) {
        z[[i]] <- setInstrSettings(z[[i]], instr.settings[[i]])
      } else {
        z[[i]] <- setInstrSettings(z[[i]], instr.settings)
      }
    }
  }
  z
}

#' Dimensions of an Object
#'
#' Retrieve or set the dimension of an object.
#'
#' @param x A \code{generic_mscpt} object or of a derived class.
#'
#' @return Either NULL or a numeric vector, which is coerced to integer (by
#'   truncation).
#'
#' @export
#'
dim.generic_mspct <- function(x) {
  z <- attr(x, "mspct.dim", exact = TRUE)
  if (!is.null(z)) {
    z <- as.integer(z)
  }
  z
}

#' @rdname dim.generic_mspct
#'
#' @param value Either NULL or a numeric vector, which is coerced to integer (by truncation).
#'
#' @export
#'
`dim<-.generic_mspct` <- function(x, value) {
  if (! is.null(value)) {
    value <- as.integer(value)
  }
  attr(x, "mspct.dim") <- value
  x
}
