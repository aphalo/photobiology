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
#' Removes from a spectrum object the class attributes "generic_mspct" and any
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
#' @param l a list or a generic_mspct object or of a derived class.
#' @param target.set character The target set of classes within which to search
#'   for classes common to all members.
#' @export
#'
#' @return A character vector containing the class attribute values.
#'
#' @family set and unset 'multi spectral' class functions
#'
shared_member_class <- function(l, target.set = spct_classes()) {
  if (length(l) == 0) {
    if (is.generic_mspct(l)) {
      # we return the least derived member class allowed
      gsub("_mspct$", "_spct", class(l)[1])
    } else {
      character()
    }
  } else {
    # we inspect the classes of members
    l.class <- target.set
    for (i in seq_along(l)) {
      member_class <- class(l[[i]])
      l.class <- intersect(l.class, member_class)
    }
    l.class
  }
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
#'
#'
#' @note Setting class = source_spct or class = source_mspct makes no difference
#'
#' @family collections of spectra classes family
#' @examples
#' filter_mspct(list(polyester.spct, yellow_gel.spct))
#'
generic_mspct <- function(l = NULL,
                          class = "generic_spct",
                          ncol = 1,
                          byrow = FALSE,
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
    spct_class <- sub("_mspct$", "_spct", class)
  } else if (class %in% spct_classes()) {
    multi_class <- sub("_spct$", "_mspct", class)
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
#'
#'
calibration_mspct <- function(l = NULL,
                              ncol = 1,
                              byrow = FALSE,
                              ...) {
  generic_mspct(l, class = "calibration_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{raw_spct} objects.
#'
#' @export
#'
#'
raw_mspct <- function(l = NULL,
                      ncol = 1,
                      byrow = FALSE,
                      ...) {
  generic_mspct(l, class = "raw_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{cps_spct} objects.
#'
#' @export
#'
#'
cps_mspct <- function(l = NULL,
                      ncol = 1,
                      byrow = FALSE, ...) {
  generic_mspct(l, class = "cps_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{source_spct} objects.
#'
#' @export
#'
#'
source_mspct <- function(l = NULL,
                         ncol = 1,
                         byrow = FALSE,
                         ...) {
  generic_mspct(l, class = "source_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{filter_spct} objects.
#'
#' @export
#'
#'
filter_mspct <- function(l = NULL,
                         ncol = 1,
                         byrow = FALSE,
                         ...) {
  generic_mspct(l, class = "filter_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{reflector_spct} objects.
#'
#' @export
#'
#'
reflector_mspct <- function(l = NULL,
                            ncol = 1,
                            byrow = FALSE,
                            ...) {
  generic_mspct(l, class = "reflector_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{object_spct} objects.
#'
#' @export
#'
#'
object_mspct <- function(l = NULL,
                         ncol = 1,
                         byrow = FALSE,
                         ...) {
  generic_mspct(l, class = "object_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{response_spct} objects.
#'
#' @export
#'
#'
response_mspct <- function(l = NULL,
                           ncol = 1,
                           byrow = FALSE,
                           ...) {
  generic_mspct(l, class = "response_spct", ncol = ncol, byrow = byrow)
}

#' @describeIn generic_mspct Specialization for collections of \code{chroma_spct} objects.
#'
#' @export
#'
#'
chroma_mspct <- function(l = NULL,
                         ncol = 1,
                         byrow = FALSE,
                         ...) {
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

# as Coercion methods for mspct classes -------------------------------------

#' @title Coerce to a collection-of-spectra
#'
#' @description Return a copy of an R object with its class set to a given type
#'   of spectrum.
#'
#' @param x a list of spectral objects or a list of objects such as data frames
#'   that can be converted into spectral objects.
#' @param force.spct.class logical indicating whether to change the class of
#'   members to \code{generic_spct} or retain the existing class.
#' @param ... passed to individual spectrum object constructor
#' @param w.length numeric A vector of wavelengthvalues sorted in strictly
#'   ascending order (nm).
#' @param spct.data.var character The name of the variable that will contain the
#'   spectral data. This indicates what physical quantity is stored in the
#'   matrix and the units of expression used.
#' @param member.class character The name of the class of the individual spectra
#'   to be constructed.
#' @param multiplier numeric A multiplier to be applied to the values in
#'   \code{x} to do unit or scale conversion.
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#' @param spct.names character Vector of names to be assigned to collection
#'   members, either of length 1, or with length equal to the number of spectra.
#'
#' @return A copy of \code{x} converted into a \code{generic_mspct} object.
#'
#' @note Members of \code{generic_mspct} objects can be heterogeneous: they can
#'   belong to any class derived from \code{generic_spct} and class is not
#'   enforced. When \code{x} is a list of data frames \code{force.spct.class =
#'   TRUE} needs to be supplied. When \code{x} is a square matrix an explicit
#'   argument is needed for \code{byrow} to indicate how data in \code{x} should
#'   be read. In every case the length of the \code{w.length} vector must match
#'   one of the dimensions of \code{x}.
#'
#' @export
#'
#' @family Coercion methods for collections of spectra
#'
as.generic_mspct <- function(x, ...) UseMethod("as.generic_mspct")

#' @describeIn as.generic_mspct
#'
#' @export
#'
as.generic_mspct.default <- function(x, ...) {
  message("'as.generic_mspct' not implemented for class: ", class(x)[1])
  generic_mspct()
}

#' @describeIn as.generic_mspct
#'
#' @export
#'
as.generic_mspct.data.frame <- function(x, force.spct.class = FALSE, ...) {
  as.generic_mspct(list(x),
                   force.spct.class = force.spct.class,
                   ...)
}

#' @describeIn as.generic_mspct
#'
#' @export
#'
as.generic_mspct.generic_spct <- function(x, force.spct.class = FALSE, ...) {
  generic_mspct(list(x),
                class = "generic_spct",
                ...)
}

#' @describeIn as.generic_mspct
#'
#' @export
#'
as.generic_mspct.list <- function(x,
                                  force.spct.class = FALSE,
                                  ...,
                                  ncol = 1,
                                  byrow = FALSE) {
  y <- x
  rmDerivedMspct(y)
  if (force.spct.class) {
    y <- plyr::llply(y, setGenericSpct, ...)
  }
  generic_mspct(y, ncol = ncol, byrow = byrow)
}

#' @describeIn as.generic_mspct
#'
#' @export
#'
as.generic_mspct.matrix <- function(x,
                                    w.length,
                                    member.class,
                                    spct.data.var,
                                    multiplier = 1,
                                    byrow = NULL,
                                    spct.names = "spct_",
                                    ...) {
  mat2mspct(x = x,
            w.length = w.length,
            member.class = member.class,
            spct.data.var = spct.data.var,
            multiplier = multiplier,
            byrow = byrow,
            spct.names = spct.names,
            ...)
}

#' @rdname as.generic_mspct
#'
#' @export
#'
mat2mspct <- function(x,
                      w.length,
                      member.class,
                      spct.data.var,
                      multiplier = 1,
                      byrow = NULL,
                      spct.names = "spct_",
                      ...) {
  stopifnot(is.matrix(x))
  if (length(spct.names) == 0) {
    spct.names = "spct"
  }
  if (is.null(byrow)) {
    if (nrow(x) == ncol(x)) {
      stop("For square matrices an argument for 'byrow' is mandatory")
    } else if (nrow(x) == length(w.length)) {
      byrow <- FALSE
    } else if (ncol(x) == length(w.length)) {
      byrow <- TRUE
    } else {
      stop("Length of 'w.length' vector is different to that of spectral data.")
    }
  }
  # spc data (spectra) can be stored as rows or as columns in a matrix,
  # consequently if stored by rows we transpose the matrix.
  if (byrow) {
    x <- t(x)
  }
  stopifnot(ncol(x) >= 1L)
  stopifnot(nrow(x) == length(w.length))
  # compatibility with as_tibble() >= 2.0.0
  if (is.null(colnames(x))) {
    colnames(x) <- as.character(1:ncol(x))
  }

  if (multiplier != 1) {
    x <- x * multiplier
  }

  y <- tibble::as_tibble(cbind(w.length, x))

  if (length(spct.names) == ncol(x)) {
    colnames(y) <- c("w.length", spct.names)
  } else {
    colnames(y) <- c("w.length", paste(spct.names[1], seq_len(ncol(x)), sep = ""))
  }

  # y contains the spectra as columns
  z <- split2mspct(x = y,
                   member.class = member.class,
                   spct.data.var = spct.data.var,
                   ncol = ncol(y),
                   ...)
  comment(z) <- paste('Converted from an R "matrix" object\n',
                      'with ', length(z), ' spectra stored ',
                      ifelse(byrow, "in rows.", "in columns."),
                      sep = "")
  z
}

#' @title Coerce to a collection-of-spectra
#'
#' @description Return a copy of an R object with its class set to a given type
#'   of spectrum.
#'
#' @param x a list of spectral objects or a list of objects such as data frames
#'   that can be converted into spectral objects.
#' @param ... passed to individual spectrum object constructor
#' @param w.length numeric A vector of wavelengthvalues sorted in strictly
#'   ascending order (nm).
#' @param spct.data.var character The name of the variable that will contain the
#'   spectral data. This indicates what physical quantity is stored in the
#'   matrix and the units of expression used.
#' @param multiplier numeric A multiplier to be applied to the values in
#'   \code{x} to do unit or scale conversion.
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#' @param spct.names character Vector of names to be assigned to collection
#'   members, either of length 1, or with length equal to the number of spectra.
#'
#' @note When \code{x} is a square matrix an explicit argument is needed for
#'   \code{byrow} to indicate how data in \code{x} should be read. In every case
#'   the length of the \code{w.length} vector must match one of the dimensions
#'   of \code{x}.
#'
#' @return A copy of \code{x} converted into a \code{calibration_mspctt} object.
#'
#' @export
#'
#' @family Coercion methods for collections of spectra
#'
as.calibration_mspct <- function(x, ...) UseMethod("as.calibration_mspct")

#' @describeIn as.calibration_mspct
#'
#' @export
#'
as.calibration_mspct.default <- function(x, ...) {
  message("'as.calibration_mspct' not implemented for class: ", class(x)[1])
  calibration_mspct()
}

#' @describeIn as.calibration_mspct
#'
#' @export
#'
as.calibration_mspct.data.frame <- function(x, ...) {
  as.calibration_mspct(x = list(x), ...)
}

#' @describeIn as.calibration_mspct
#'
#' @export
#'
as.calibration_mspct.calibration_spct <- function(x,
                                                  ...) {
  calibration_mspct(list(x))
}

#' @describeIn as.calibration_mspct
#'
#' @export
#'
as.calibration_mspct.list <- function(x,
                                      ...,
                                      ncol = 1,
                                      byrow = FALSE) {
  y <- x
  rmDerivedMspct(y)
  stopifnot(all(sapply(y, FUN = is.list)))
  z <- plyr::llply(y, setCalibrationSpct, ...)
  calibration_mspct(z, ncol = ncol, byrow = byrow)
}

#' @describeIn as.calibration_mspct
#'
#' @export
#'
as.calibration_mspct.matrix <- function(x,
                                        w.length,
                                        spct.data.var = "irrad.mult",
                                        multiplier = 1,
                                        byrow = NULL,
                                        spct.names = "spct_",
                                        ...) {
  mat2mspct(x = x,
            w.length = w.length,
            member.class = "calibration_spct",
            spct.data.var = spct.data.var,
            multiplier = multiplier,
            byrow = byrow,
            spct.names = spct.names,
            ...)
}

#' @title Coerce to a collection-of-spectra
#'
#' @description Return a copy of an R object with its class set to a given type
#'   of spectrum.
#'
#' @param x a list of spectral objects or a list of objects such as data frames
#'   that can be converted into spectral objects.
#' @param ... passed to individual spectrum object constructor
#' @param w.length numeric A vector of wavelengthvalues sorted in strictly
#'   ascending order (nm).
#' @param spct.data.var character The name of the variable that will contain the
#'   spectral data. This indicates what physical quantity is stored in the
#'   matrix and the units of expression used.
#' @param multiplier numeric A multiplier to be applied to the values in
#'   \code{x} to do unit or scale conversion.
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#' @param spct.names character Vector of names to be assigned to collection
#'   members, either of length 1, or with length equal to the number of spectra.
#'
#' @note When \code{x} is a square matrix an explicit argument is needed for
#'   \code{byrow} to indicate how data in \code{x} should be read. In every case
#'   the length of the \code{w.length} vector must match one of the dimensions
#'   of \code{x}.
#'
#'
#' @return A copy of \code{x} converted into a \code{raw_mspct} object.
#'
#' @export
#'
#' @family Coercion methods for collections of spectra
#'
as.raw_mspct <- function(x, ...) UseMethod("as.raw_mspct")

#' @describeIn as.raw_mspct
#'
#' @export
#'
as.raw_mspct.default <- function(x, ...) {
  message("'as.raw_mspct' not implemented for class: ", class(x)[1])
  raw_mspct()
}

#' @describeIn as.raw_mspct
#'
#' @export
#'
as.raw_mspct.data.frame <- function(x, ...) {
  as.raw_mspct(x = list(x), ...)
}

#' @describeIn as.raw_mspct
#'
#' @export
#'
as.raw_mspct.raw_spct <- function(x, ...) {
  raw_mspct(list(x), ...)
}

#' @describeIn as.raw_mspct
#'
#' @export
#'
as.raw_mspct.list <- function(x,
                              ...,
                              ncol = 1,
                              byrow = FALSE) {
  y <- x
  rmDerivedMspct(y)
  stopifnot(all(sapply(y, FUN = is.list)))
  z <- plyr::llply(y, setRawSpct, ...)
  raw_mspct(z, ncol = ncol, byrow = byrow)
}

#' @describeIn as.raw_mspct
#'
#' @export
#'
as.raw_mspct.matrix <- function(x,
                                w.length,
                                spct.data.var = "counts",
                                multiplier = 1,
                                byrow = NULL,
                                spct.names = "spct_",
                                ...) {
  mat2mspct(x = x,
            w.length = w.length,
            member.class = "raw_spct",
            spct.data.var = spct.data.var,
            multiplier = multiplier,
            byrow = byrow,
            spct.names = spct.names,
            ...)
}

#' @title Coerce to a collection-of-spectra
#'
#' @description Return a copy of an R object with its class set to a given type
#'   of spectrum.
#'
#' @param x a list of spectral objects or a list of objects such as data frames
#'   that can be converted into spectral objects.
#' @param ... passed to individual spectrum object constructor
#' @param w.length numeric A vector of wavelengthvalues sorted in strictly
#'   ascending order (nm).
#' @param spct.data.var character The name of the variable that will contain the
#'   spectral data. This indicates what physical quantity is stored in the
#'   matrix and the units of expression used.
#' @param multiplier numeric A multiplier to be applied to the values in
#'   \code{x} to do unit or scale conversion.
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#' @param spct.names character Vector of names to be assigned to collection
#'   members, either of length 1, or with length equal to the number of spectra.
#'
#' @note When \code{x} is a square matrix an explicit argument is needed for
#'   \code{byrow} to indicate how data in \code{x} should be read. In every case
#'   the length of the \code{w.length} vector must match one of the dimensions
#'   of \code{x}.
#'
#' @return A copy of \code{x} converted into a \code{cps_mspct} object.
#'
#' @export
#'
#' @family Coercion methods for collections of spectra
#'
as.cps_mspct <- function(x, ...) UseMethod("as.cps_mspct")

#' @describeIn as.cps_mspct
#'
#' @export
#'
as.cps_mspct.default <- function(x, ...) {
  message("'as.cps_mspct' not implemented for class: ", class(x)[1])
  cps_mspct()
}

#' @describeIn as.cps_mspct
#'
#' @export
#'
as.cps_mspct.data.frame <- function(x, ...) {
  as.cps_mspct(x = list(x), ...)
}

#' @describeIn as.cps_mspct
#'
#' @export
#'
as.cps_mspct.cps_spct <- function(x, ...) {
  cps_mspct(list(x), ...)
}

#' @describeIn as.cps_mspct
#'
#' @export
#'
as.cps_mspct.list <- function(x,
                              ...,
                              ncol = 1,
                              byrow = FALSE) {
  y <- x
  rmDerivedMspct(y)
  stopifnot(all(sapply(y, FUN = is.list)))
  z <- plyr::llply(y, setCpsSpct, ...)
  cps_mspct(z, ncol = ncol, byrow = byrow)
}

#' @describeIn as.cps_mspct
#'
#' @export
#'
as.cps_mspct.matrix <- function(x,
                                w.length,
                                spct.data.var = "cps",
                                multiplier = 1,
                                byrow = NULL,
                                spct.names = "spct_",
                                ...) {
  mat2mspct(x = x,
            w.length = w.length,
            member.class = "cps_spct",
            spct.data.var = spct.data.var,
            multiplier = multiplier,
            byrow = byrow,
            spct.names = spct.names,
            ...)
}

#' @title Coerce to a collection-of-spectra
#'
#' @description Return a copy of an R object with its class set to a given type
#'   of spectrum.
#'
#' @param x a list of spectral objects or a list of objects such as data frames
#'   that can be converted into spectral objects.
#' @param time.unit character A string, "second", "day" or "exposure"
#' @param bswf.used character
#' @param strict.range logical Flag indicating how off-range values are handled
#' @param ... passed to individual spectrum object constructor
#' @param w.length numeric A vector of wavelengthvalues sorted in strictly
#'   ascending order (nm).
#' @param spct.data.var character The name of the variable that will contain the
#'   spectral data. This indicates what physical quantity is stored in the
#'   matrix and the units of expression used.
#' @param multiplier numeric A multiplier to be applied to the values in
#'   \code{x} to do unit or scale conversion.
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#' @param spct.names character Vector of names to be assigned to collection
#'   members, either of length 1, or with length equal to the number of spectra.
#'
#' @note When \code{x} is a square matrix an explicit argument is needed for
#'   \code{byrow} to indicate how data in \code{x} should be read. In every case
#'   the length of the \code{w.length} vector must match one of the dimensions
#'   of \code{x}.
#'
#'
#' @return A copy of \code{x} converted into a \code{source_mspct} object.
#'
#' @export
#'
#' @family Coercion methods for collections of spectra
#'
as.source_mspct <- function(x, ...) UseMethod("as.source_mspct")

#' @describeIn as.source_mspct
#'
#' @export
#'
as.source_mspct.default <- function(x, ...) {
  message("'as.source_mspct' not implemented for class: ", class(x)[1])
  source_mspct()
}

#' @describeIn as.source_mspct
#'
#' @export
#'
as.source_mspct.data.frame <-
  function(x,
           time.unit=c("second", "day", "exposure"),
           bswf.used=c("none", "unknown"),
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           ...) {
    as.source_mspct(x = list(x),
                    time.unit = time.unit,
                    bswf.used = bswf.used,
                    strict.range = strict.range,
                    ...)
  }

#' @describeIn as.source_mspct
#'
#' @export
#'
as.source_mspct.source_spct <- function(x, ...) {
  source_mspct(list(x), ...)
}

#' @describeIn as.source_mspct
#'
#' @export
#'
as.source_mspct.list <-
  function(x,
           time.unit=c("second", "day", "exposure"),
           bswf.used=c("none", "unknown"),
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           ...,
           ncol = 1,
           byrow = FALSE) {
    y <- x
    rmDerivedMspct(y)
    stopifnot(all(sapply(y, FUN = is.list)))
    z <- plyr::llply(y, setSourceSpct, time.unit = time.unit,
                     strict.range = strict.range, bswf.used = bswf.used, ...)
    source_mspct(z, ncol = ncol, byrow = byrow)
  }

#' @describeIn as.source_mspct
#'
#' @export
#'
as.source_mspct.matrix <- function(x,
                                   w.length,
                                   spct.data.var = "s.e.irrad",
                                   multiplier = 1,
                                   byrow = NULL,
                                   spct.names = "spct_",
                                   ...) {
  mat2mspct(x = x,
            w.length = w.length,
            member.class = "source_spct",
            spct.data.var = spct.data.var,
            multiplier = multiplier,
            byrow = byrow,
            spct.names = spct.names,
            ...)
}

#' @title Coerce to a collection-of-spectra
#'
#' @description Return a copy of an R object with its class set to a given type
#'   of spectrum.
#'
#' @param x a list of spectral objects or a list of objects such as data frames
#'   that can be converted into spectral objects.
#' @param time.unit character A string, "second", "day" or "exposure"
#' @param ... passed to individual spectrum object constructor
#' @param w.length numeric A vector of wavelengthvalues sorted in strictly
#'   ascending order (nm).
#' @param spct.data.var character The name of the variable that will contain the
#'   spectral data. This indicates what physical quantity is stored in the
#'   matrix and the units of expression used.
#' @param multiplier numeric A multiplier to be applied to the values in
#'   \code{x} to do unit or scale conversion.
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#' @param spct.names character Vector of names to be assigned to collection
#'   members, either of length 1, or with length equal to the number of spectra.
#'
#' @note When \code{x} is a square matrix an explicit argument is needed for
#'   \code{byrow} to indicate how data in \code{x} should be read. In every case
#'   the length of the \code{w.length} vector must match one of the dimensions
#'   of \code{x}.
#'
#' @return A copy of \code{x} converted into a \code{response_mspct} object.
#'
#' @export
#'
#' @family Coercion methods for collections of spectra
#'
as.response_mspct <- function(x, ...) UseMethod("as.response_mspct")

#' @describeIn as.response_mspct
#'
#' @export
#'
as.response_mspct.default <- function(x, ...) {
  message("'as.response_mspct' not implemented for class: ", class(x)[1])
  response_mspct()
}

#' @describeIn as.response_mspct
#'
#' @export
#'
as.response_mspct.data.frame <-
  function(x,
           time.unit= "second",
           ...) {
    as.source_mspct(x = list(x),
                    time.unit = time.unit,
                    ...)
  }

#' @describeIn as.response_mspct
#'
#' @export
#'
as.response_mspct.response_spct <- function(x, ...) {
  response_mspct(list(x), ...)
}

#' @describeIn as.response_mspct
#'
#' @export
#'
as.response_mspct.list <- function(x,
                                   time.unit = "second",
                                   ...,
                                   ncol = 1,
                                   byrow = FALSE) {
  y <- x
  rmDerivedMspct(y)
  stopifnot(all(sapply(y, FUN = is.list)))
  z <- plyr::llply(y, setResponseSpct, time.unit = time.unit, ...)
  response_mspct(z, ncol = ncol, byrow = byrow)
}

#' @describeIn as.response_mspct
#'
#' @export
#'
as.response_mspct.matrix <- function(x,
                                     w.length,
                                     spct.data.var = "s.e.response",
                                     multiplier = 1,
                                     byrow = NULL,
                                     spct.names = "spct_",
                                     ...) {
  mat2mspct(x = x,
            w.length = w.length,
            member.class = "response_spct",
            spct.data.var = spct.data.var,
            multiplier = multiplier,
            byrow = byrow,
            spct.names = spct.names,
            ...)
}


#' @title Coerce to a collection-of-spectra
#'
#' @description Return a copy of an R object with its class set to a given type
#'   of spectrum.
#'
#' @param x a list of spectral objects or a list of objects such as data frames
#'   that can be converted into spectral objects.
#' @param Tfr.type a character string, either "total" or "internal"
#' @param strict.range logical Flag indicating how off-range values are handled
#' @param ... passed to individual spectrum object constructor
#' @param w.length numeric A vector of wavelengthvalues sorted in strictly
#'   ascending order (nm).
#' @param spct.data.var character The name of the variable that will contain the
#'   spectral data. This indicates what physical quantity is stored in the
#'   matrix and the units of expression used.
#' @param multiplier numeric A multiplier to be applied to the values in
#'   \code{x} to do unit or scale conversion.
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#' @param spct.names character Vector of names to be assigned to collection
#'   members, either of length 1, or with length equal to the number of spectra.
#'
#' @note When \code{x} is a square matrix an explicit argument is needed for
#'   \code{byrow} to indicate how data in \code{x} should be read. In every case
#'   the length of the \code{w.length} vector must match one of the dimensions
#'   of \code{x}.
#'
#' @return A copy of \code{x} converted into a \code{filter_mspct} object.
#'
#' @export
#'
#' @family Coercion methods for collections of spectra
#'
as.filter_mspct <- function(x, ...) UseMethod("as.filter_mspct")

#' @describeIn as.filter_mspct
#'
#' @export
#'
as.filter_mspct.default <- function(x, ...) {
  message("'as.filter_mspct' not implemented for class: ", class(x)[1])
  filter_mspct()
}

#' @describeIn as.filter_mspct
#'
#' @export
#'
as.filter_mspct.data.frame <-
  function(x,
           Tfr.type = c("total", "internal"),
           strict.range = TRUE,
           ...) {
    as.filter_mspct(x = list(x),
                    Tfr.type = Tfr.type,
                    strict.range = strict.range,
                    ...)
  }

#' @describeIn as.filter_mspct
#'
#' @export
#'
as.filter_mspct.filter_spct <- function(x, ...) {
  filter_mspct(list(x), ...)
}

#' @describeIn as.filter_mspct
#'
#' @export
#'
as.filter_mspct.list <- function(x,
                            Tfr.type = c("total", "internal"),
                            strict.range = TRUE,
                            ...,
                            ncol = 1,
                            byrow = FALSE) {
  y <- x
  rmDerivedMspct(y)
  stopifnot(all(sapply(y, FUN = is.list)))
  z <- plyr::llply(y, setFilterSpct,
                   Tfr.type = Tfr.type,
                   strict.range = strict.range,
                   ...)
  filter_mspct(z, ncol = ncol, byrow = byrow)
}

#' @describeIn as.filter_mspct
#'
#' @export
#'
as.filter_mspct.matrix <- function(x,
                                   w.length,
                                   spct.data.var = "Tfr",
                                   multiplier = 1,
                                   byrow = NULL,
                                   spct.names = "spct_",
                                   ...) {
  mat2mspct(x = x,
            w.length = w.length,
            member.class = "filter_spct",
            spct.data.var = spct.data.var,
            multiplier = multiplier,
            byrow = byrow,
            spct.names = spct.names,
            ...)
}

#' @title Coerce to a collection-of-spectra
#'
#' @description Return a copy of an R object with its class set to a given type
#'   of spectrum.
#'
#' @param x a list of spectral objects or a list of objects such as data frames
#'   that can be converted into spectral objects.
#' @param Rfr.type a character string, either "total" or "specular"
#' @param strict.range logical Flag indicating how off-range values are handled
#' @param ... passed to individual spectrum object constructor
#' @param w.length numeric A vector of wavelengthvalues sorted in strictly
#'   ascending order (nm).
#' @param spct.data.var character The name of the variable that will contain the
#'   spectral data. This indicates what physical quantity is stored in the
#'   matrix and the units of expression used.
#' @param multiplier numeric A multiplier to be applied to the values in
#'   \code{x} to do unit or scale conversion.
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#' @param spct.names character Vector of names to be assigned to collection
#'   members, either of length 1, or with length equal to the number of spectra.
#'
#' @note When \code{x} is a square matrix an explicit argument is needed for
#'   \code{byrow} to indicate how data in \code{x} should be read. In every case
#'   the length of the \code{w.length} vector must match one of the dimensions
#'   of \code{x}.
#'
#' @return A copy of \code{x} converted into a \code{reflector_mspct} object.
#'
#' @export
#'
#' @family Coercion methods for collections of spectra
#'
as.reflector_mspct <- function(x, ...) UseMethod("as.reflector_mspct")

#' @describeIn as.reflector_mspct
#'
#' @export
#'
as.reflector_mspct.default <- function(x, ...) {
  message("'as.reflector_mspct' not implemented for class: ", class(x)[1])
  reflector_mspct()
}

#' @describeIn as.reflector_mspct
#'
#' @export
#'
as.reflector_mspct.data.frame <-
  function(x,
           Rfr.type = c("total", "specular"),
           strict.range = TRUE,
           ...) {
    as.filter_mspct(x = list(x),
                    Rfr.type = Rfr.type,
                    strict.range = strict.range,
                    ...)
  }

#' @describeIn as.reflector_mspct
#'
#' @export
#'
as.reflector_mspct.reflector_spct <- function(x, ...) {
  reflector_mspct(list(x), ...)
}

#' @describeIn as.reflector_mspct
#'
#' @export
#'
as.reflector_mspct.list <- function(x,
                                    Rfr.type = c("total", "specular"),
                                    strict.range = TRUE,
                                    ...,
                                    ncol = 1,
                                    byrow = FALSE) {
  y <- x
  rmDerivedMspct(y)
  z <- plyr::llply(y,
                   setReflectorSpct,
                   Rfr.type = Rfr.type,
                   strict.range = strict.range,
                   ...)
  reflector_mspct(z, ncol = ncol, byrow = byrow)
}

#' @describeIn as.reflector_mspct
#'
#' @export
#'
as.reflector_mspct.matrix <- function(x,
                                      w.length,
                                      spct.data.var = "Rfr",
                                      multiplier = 1,
                                      byrow = NULL,
                                      spct.names = "spct_",
                                      ...) {
  mat2mspct(x = x,
            w.length = w.length,
            member.class = "reflector_spct",
            spct.data.var = spct.data.var,
            multiplier = multiplier,
            byrow = byrow,
            spct.names = spct.names,
            ...)
}

#' @title Coerce to a collection-of-spectra
#'
#' @description Return a copy of an R object with its class set to a given type
#'   of spectrum.
#'
#' @param x a list of spectral objects or a list of objects such as data frames
#'   that can be converted into spectral objects.
#' @param Tfr.type a character string, either "total" or "internal"
#' @param Rfr.type a character string, either "total" or "specular"
#' @param strict.range logical Flag indicating how off-range values are handled
#' @param ... passed to individual spectrum object constructor
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#'
#'
#' @return A copy of \code{x} converted into a \code{object_mspct} object.
#'
#' @export
#'
#' @family Coercion methods for collections of spectra
#'
as.object_mspct <- function(x, ...) UseMethod("as.object_mspct")

#' @describeIn as.object_mspct
#'
#' @export
#'
as.object_mspct.default <- function(x, ...) {
  message("'as.object_mspct' not implemented for class: ", class(x)[1])
  object_mspct()
}

#' @describeIn as.object_mspct
#'
#' @export
#'
as.object_mspct.data.frame <-
  function(x,
           Tfr.type=c("total", "internal"),
           Rfr.type=c("total", "specular"),
           strict.range = TRUE,
           ...) {
    as.object_mspct(x = list(x),
                    Tfr.type = Tfr.type,
                    Rfr.type = Rfr.type,
                    strict.range = strict.range,
                    ...)
  }

#' @describeIn as.object_mspct
#'
#' @export
#'
as.object_mspct.object_spct <- function(x, ...) {
  object_mspct(list(x), ...)
}

#' @describeIn as.object_mspct
#'
#' @export
#'
as.object_mspct.list <- function(x,
                            Tfr.type=c("total", "internal"),
                            Rfr.type=c("total", "specular"),
                            strict.range = TRUE,
                            ...,
                            ncol = 1,
                            byrow = FALSE) {
  y <- x
  rmDerivedMspct(y)
  z <- plyr::llply(y,
                   setObjectSpct,
                   Tfr.type = Tfr.type,
                   Rfr.type = Rfr.type,
                   strict.range = strict.range,
                   ...)
  object_mspct(z, ncol = ncol, byrow = byrow)
}

#' @title Coerce to a collection-of-spectra
#'
#' @description Return a copy of an R object with its class set to a given type
#'   of spectrum.
#'
#' @param x a list of spectral objects or a list of objects such as data frames
#'   that can be converted into spectral objects.
#' @param ... passed to individual spectrum object constructor
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#'
#'
#' @return A copy of \code{x} converted into a \code{chroma_mspct} object.
#'
#' @export
#'
#' @family Coercion methods for collections of spectra
#'
as.chroma_mspct <- function(x, ...) UseMethod("as.chroma_mspct")

#' @describeIn as.chroma_mspct
#'
#' @export
#'
as.chroma_mspct.default <- function(x, ...) {
  message("'as.chroma_mspct' not implemented for class: ", class(x)[1])
  chroma_mspct()
}

#' @describeIn as.chroma_mspct
#'
#' @export
#'
as.chroma_mspct.data.frame <- function(x, ...) {
  as.chroma_mspct(x = list(x), ...)
}

#' @describeIn as.chroma_mspct
#'
#' @export
#'
as.chroma_mspct.chroma_spct <- function(x, ...) {
  chroma_mspct(list(x), ...)
}

#' @describeIn as.chroma_mspct
#'
#' @export
#'
as.chroma_mspct.list <- function(x,
                                 ...,
                                 ncol = 1,
                                 byrow = FALSE) {
  y <- x
  rmDerivedMspct(y)
  z <- plyr::llply(y, setChromaSpct, ...)
  chroma_mspct(z, ncol = ncol, byrow = byrow)
}


# coerce to matrix --------------------------------------------------------

#' Coerce a collection of spectra into a matrix
#'
#' Convert an object of class \code{generic_mspct} or a derived class into an R
#' matrix with wavelengths saved as an attribute and spectral data in rows or
#' columns.
#'
#' @note Only collections of spectra containing spectra with exactly the same
#'   \code{w.length} values can by converted. If needed, the spectra can be
#'   re-expressed before attempting the conversion to a matrix.
#'
#' @param x generic_mspct object.
#' @param spct.data.var character The name of the variable containing the
#'   spectral data.
#' @param byrow logical. If FALSE (the default) the matrix is filled with the
#'   spectra stored by columns, otherwise the matrix is filled by rows.
#' @param ... currently ignored.
#'
#' @section Warning!: This conversion preserves the spectral data but discards
#'   almost all the metadata contained in the spectral objects. In other words a
#'   matrix created with this function cannot be used to recreate the original
#'   object unless the same metadata is explicitly supplied when converting the
#'   matrix into new collection of spectra.
#'
#' @export
#'
#' @name as.matrix-mspct
#'
as.matrix.generic_mspct <- function(x,
                                    spct.data.var,
                                    byrow = attr(x, "mspct.byrow"),
                                    ...) {
  mspct2mat(x = x,
            spct.data.var = spct.data.var,
            byrow = byrow,
            ...)
}

#' @rdname as.matrix-mspct
#'
#' @export
#'
mspct2mat <- function(x,
                      spct.data.var,
                      byrow = attr(x, "mspct.byrow"),
                      ...) {
  stopifnot(is.any_mspct(x))
  if (length(x) == 0L) {
    return(matrix(numeric()))
  }
  spct.names <- names(x)
  spct.selector <- rep(TRUE, length(x))
  mat <- numeric()
  for (i in seq_along(x)) {
    temp <- x[[i]]
    s.column <- temp[[spct.data.var]]
    wl.current <- temp[["w.length"]]
    if (i == 1L) {
      wl.prev <- wl.current
    }
    if (!all(wl.current == wl.prev) || length(s.column) == 0L) {
      spct.selector[i] <- FALSE
      next()
    }
    mat <- c(mat, s.column) # one long numeric vector
  }
  if (any(!spct.selector)) {
    warning("Spectra dropped: ", sum(!spct.selector), " out of ",
            length(spct.selector), ".")
  }
  if (byrow) {
    z <- matrix(mat, nrow = sum(spct.selector), byrow = byrow,
                dimnames = list(spct = c(spct.names[spct.selector]),
                                w.length = wl.prev))
  } else {
    z <- matrix(mat, ncol = sum(spct.selector), byrow = byrow,
                dimnames = list(w.length = wl.prev,
                                spct = c(spct.names[spct.selector])))
  }
  attr(z, "w.length") <- wl.prev
  comment(z) <- comment(x)
  z
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
#' @param idx.var character Name of column containing data to be copied
#'   unchanged to each spct object
#' @param ncol integer Number of 'virtual' columns in data
#' @param byrow logical If \code{ncol > 1} how to read in the data
#' @param ... additional named arguments passed to the member constructor
#'   function.
#'
#' @export
#'
#' @family Coercion methods for collections of spectra
#'
split2mspct <- function(x,
                        member.class = NULL,
                        spct.data.var = NULL,
                        w.length.var = "w.length",
                        idx.var = NULL,
                        ncol = 1,
                        byrow = FALSE,
                        ...) {
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
#' @family Coercion methods for collections of spectra
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
  what.measured <- getWhatMeasured(x)
  # these methods return a data.frame
  where.measured <- getWhereMeasured(x)
  # these methods may return an empty list
  instr.desc <- getInstrDesc(x)
  instr.settings <- getInstrSettings(x)
  filter.properties <- getFilterProperties(x, return.null = TRUE)
  if (is.null(filter.properties)) {
    filter.properties <- list()
  }
  for (i in seq(along.with = z)) {
    if (!all(is.na(when.measured))) {
      if (is.list(when.measured) && length(when.measured) == length(groups)) {
        z[[i]] <- setWhenMeasured(z[[i]], when.measured[[i]])
      } else {
        z[[i]] <- setWhenMeasured(z[[i]], when.measured)
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
    if (length(filter.properties) > 0) {
      if (is.list(filter.properties) &&
          !inherits(filter.properties, "filter_properties") &&
          length(filter.properties) == length(groups)) {
        z[[i]] <- setFilterProperties(z[[i]], filter.properties[[i]])
      } else {
        z[[i]] <- setFilterProperties(z[[i]], filter.properties)
      }
    }
  }
  z <- setWhereMeasured(z, where.measured)
  z
}

#' Dimensions of an Object
#'
#' Retrieve or set the dimension of an object.
#'
#' @param x A \code{generic_mspct} object or of a derived class.
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
#' @param value Either NULL or a numeric vector, which is coerced to integer (by
#'   truncation).
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
