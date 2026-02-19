# transmittance and absorbance --------------------------------------------

# A2T ---------------------------------------------------------------------

#' Convert absorbance into transmittance
#'
#' Function that converts absorbance (a.u.) into transmittance (fraction).
#'
#' @param x an R object.
#' @param action a character string "add" or "replace".
#' @param byref logical indicating if new object will be created by reference or
#'   by copy of \code{x}.
#' @param ... not used in current version.
#'
#' @return A copy of \code{x} with a column \code{Tfr} added and \code{A} and
#'   \code{Afr} possibly deleted except for \code{w.length}. If \code{action =
#'   "replace"}, in all cases, the additional columns are removed, even if no
#'   column needs to be added.
#'
#' @export
#' @family quantity conversion functions
#'
A2T <- function(x, action, byref, ...) UseMethod("A2T")

#' @describeIn A2T Default method for generic function
#'
#' @export
#'
A2T.default <- function(x, action = NULL, byref = FALSE, ...) {
  warning("'A2T()' not implemented for class \"", class(x)[1], "\".")
  return(x)
}

#' @describeIn A2T method for numeric vectors
#'
#' @export
#'
A2T.numeric <- function(x, action = NULL, byref = FALSE, ...) {
  return(10^-x)
}

#' @describeIn A2T Method for filter spectra
#'
#' @export
#'
A2T.filter_spct <- function(x, action = "add", byref = FALSE, ...) {

  if (byref) {
    name <- substitute(x)
  }

  action <- check_action_arg(x, action)
  # nothing to do
  if (action == "replace" &&
      "Tfr" %in% colnames(x) && !any(c("A", "Afr") %in% colnames(x))) {
    return(x)
  }
  if (action == "add" && "Tfr" %in% colnames(x)) {
    return(x)
  }
  # we remove normalization
  ## should skip this when only removing s.e.irrad by removing col from normalization
  if (is_normalised(x)) {
    old.normalization.ls <- getNormalization(x)
    x <- denormalize_spct(x) # remove normalization
    norm.action <- "update"
  } else {
    norm.action <- "skip"
  }

  if (!exists("Tfr", x, inherits = FALSE)) {
    if (exists("A", x, inherits = FALSE)) {
      x[["Tfr"]] <- 10^-x[["A"]]
    } else {
      x[["Tfr"]] <- NA_real_
      action <- "add"
      warning("'A' missing in 'A2T()', filling 'Tfr' with 'NA'.")
    }
  }

  if (action == "replace" && exists("A", x, inherits = FALSE)) {
    x[["A"]] <- NULL
  }
  if (action == "replace" && exists("Afr", x, inherits = FALSE)) {
    x[["Afr"]] <- NULL
  }

  if (norm.action == "update") {
    x <- restore_normalization(x, old.normalization.ls)
  }

  if (byref && is.name(name)) { # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

#' @describeIn A2T Method for collections of filter spectra
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
A2T.filter_mspct <- function(x,
                             action = "add",
                             byref = FALSE,
                             ...,
                             .parallel = FALSE,
                             .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          .fun = A2T.filter_spct,
          action = action,
          byref = byref,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}


# T2A ---------------------------------------------------------------------


#' Convert transmittance into absorbance.
#'
#' Function that converts transmittance (fraction) into \eqn{\log_{10}}-based
#' absorbance (a.u.).
#'
#' @details
#' Absorbance, \eqn{A}, is frequently used in chemistry as it is linearly
#' related to the concentration of a solute dissolved in a solvent.
#'
#' \deqn{A = -\log_{10} \tau}
#'
#' where, \eqn{A} absorbance and \eqn{\tau} is internal transmittance. By
#' default, if total transmittance, \eqn{T}, is stored in \code{x}, the
#' returned value computed as
#'
#' \deqn{A = -\log_{10} T}
#'
#' is not strictly absorbance. In this
#' case and in cases when the measured light attenuation is the result of
#' scattering, or when part of measured light is re-emitted after absorption
#' the use of \emph{attenuance} is the IUPAC-recommended name for this quantity.
#'
#' If \code{strict.A = TRUE} is passed in the call and total
#' transmittance, \eqn{T}, and total
#' reflectance, \eqn{\rho}, are both available, absorbance is computed as:
#'
#' \deqn{A = -\log_{10} (T - \rho) / (1 - \rho)}
#'
#' where \eqn{\rho} can be either spectral total reflectance stored in \code{x}
#' as data or a single approximate \code{Rfr.constant} value stored as part
#' of the metadata.
#'
#' @note The default \code{A.strict = FALSE} ensures indentical behaviour
#'   as in 'photobiology' (<= 0.11.0).
#'
#' @param x an R object.
#' @param action character Allowed values \code{"replace"} and \code{"add"}.
#' @param byref logical indicating if new object will be created by reference
#'   or by copy of \code{x}.
#' @param clean logical replace off-boundary values before conversion
#' @param strict.A logical Attempt to compute a true internal absorbance even
#'   if \code{"total"} transmittance is stored in \code{x}.
#' @param ... not used in current version
#'
#' @return A copy of \code{x} with a column \code{A} added and other columns
#'   possibly deleted except for \code{w.length}. If \code{action = "replace"},
#'   in all cases, the additional columns are removed, even if no column needs
#'   to be added.
#'
#' @export
#' @family quantity conversion functions
#'
T2A <- function(x, action, byref, clean, ...) UseMethod("T2A")

#' @describeIn T2A Default method for generic function
#'
#' @export
#'
T2A.default <- function(x, action = NULL, byref = FALSE, ...) {
  warning("'T2A()' not implemented for class \"", class(x)[1], "\".")
  return(x)
}

#' @describeIn T2A Method for numeric vectors
#'
#' @export
#'
T2A.numeric <- function(x,
                        action = NULL,
                        byref = FALSE,
                        clean = TRUE,
                        ...) {
  if (clean) {
    Tfr.zero <- getOption("photobiology.Tfr.zero", default = 0)
    if (any(x < Tfr.zero)) {
      warning("Replacing ", sum(x < Tfr.zero), " values < ",
              Tfr.zero, " with ", Tfr.zero)
      x <- ifelse(x <= 0, Tfr.zero, x)
    }
  }
  return(-log10(x))
}

#' @describeIn T2A Method for filter spectra
#'
#' @export
#'
T2A.filter_spct <- function(x,
                            action = "add",
                            byref = FALSE,
                            clean = TRUE,
                            strict.A = FALSE,
                            ...) {
  if (byref) {
    name <- substitute(x)
  }

  action <- check_action_arg(x, action)
  # nothing to do
  if (action == "replace" &&
      "A" %in% colnames(x) && !any(c("Tfr", "Afr") %in% colnames(x))) {
    return(x)
  }
  if (action == "add" && "A" %in% colnames(x)) {
    return(x)
  }
  # we remove normalization
  ## should skip this when only removing s.e.irrad by removing col from normalization
  if (is_normalised(x)) {
    old.normalization.ls <- getNormalization(x)
    x <- denormalize_spct(x) # remove normalization
    norm.action <- "update"
  } else {
    norm.action <- "skip"
  }

  if (!exists("A", x, inherits = FALSE)) {
    if (exists("Tfr", x, inherits = FALSE)) {
      x[["A"]] <- -log10(x[["Tfr"]])
    } else {
      x[["A"]] <- NA_real_
      action <- "add"
      warning("'Tfr' missing in 'T2A()', filling 'A' with 'NA'.")
    }
  }

  if (action == "replace" && exists("Tfr", x, inherits = FALSE)) {
    x[["Tfr"]] <- NULL
  }
  if (action == "replace" && exists("Afr", x, inherits = FALSE)) {
    x[["Afr"]] <- NULL
  }

  if (norm.action == "update") {
    x <- restore_normalization(x, old.normalization.ls)
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  if (any(is.infinite(x[["A"]]))) {
    warning("'Inf' absorbance values generated as some Tfr values ",
            "were equal to zero!")
  }
  return(x)
}

#' @describeIn T2A Method for collections of filter spectra
#'
#' @param .parallel	if \code{TRUE}, apply function in parallel, using parallel
#'   backend provided by foreach.
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
T2A.filter_mspct <- function(x,
                             action = "add",
                             byref = FALSE,
                             clean = TRUE,
                             strict.A = TRUE,
                             ...,
                             .parallel = FALSE,
                             .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          .fun = T2A.filter_spct,
          action = action,
          byref = byref,
          clean = clean,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

# T2Afr ---------------------------------------------------------------------

#' Convert transmittance into absorptance.
#'
#' Function that converts transmittance (fraction) into absorptance (fraction).
#' If reflectance (fraction) is available, it also allows conversions between
#' internal and total absorptance.
#'
#' @details
#' Absorptance, internal transmittance and total reflectance when expressed as
#' fractions, add up to one:
#'
#' \deqn{1 = \alpha + \tau + \rho}
#'
#' where, \eqn{\alpha} is absorptance, \eqn{\tau} is internal transmittance and
#' \eqn{\rho} is total reflectance. If any two of these quantities are known,
#' the third one can be computed from them.
#'
#' On the other hand:
#'
#' \deqn{1 = \alpha\prime + T}
#'
#' where, \eqn{\alpha\prime = \alpha + \rho}, measured together. In this case,
#' there is not enough information available to compute \eqn{\alpha}.
#'
#' Thus, method \code{T2Afr()} computes
#' either \eqn{\alpha} or \eqn{\alpha\prime},
#' depending on whether \eqn{\tau} or \eqn{T} are contained in the argument
#' passed to \code{x}, but neither of them when only \eqn{\tau} is known. To
#' know which quantity has been computed, use \code{getTfrType()} to query
#' whether the computations were based on \eqn{\tau} or \eqn{T}.
#'
#' The R names used are: \code{Tfr} for \eqn{\tau} and \eqn{T} are \code{Tfr},
#' \code{Afr} for \eqn{\alpha} and \eqn{\alpha\prime}, and \code{Rfr} for
#' \eqn{rho}. The distinction between \eqn{\tau} and \eqn{T} and
#' between \eqn{\alpha} and \eqn{\alpha\prime} is made based on metadata
#' attributes.
#'
#' @param x an R object.
#' @param action character Allowed values \code{"replace"} and \code{"add"}.
#' @param byref logical indicating if new object will be created by reference
#' or by copy of \code{x}.
#' @param clean logical replace off-boundary values before conversion.
#' @param ... not used in current version.
#'
#' @return A copy of \code{x} with a column \code{Afr} added and other columns
#'   possibly deleted except for \code{w.length}. If \code{action = "replace"},
#'   in all cases, the redundant columns are removed, even when
#'   column \code{Afr} was present in the argument passed to \code{x}.
#'
#' @export
#' @family quantity conversion functions
#'
#' @examples
#' T2Afr(Ler_leaf.spct)
#'
T2Afr <- function(x, action, byref, clean, ...) UseMethod("T2Afr")

#' @describeIn T2Afr Default method for generic function
#'
#' @export
#'
T2Afr.default <- function(x,
                          action = NULL,
                          byref = FALSE,
                          clean = FALSE,
                          ...) {
  warning("'T2Afr()' not implemented for class \"", class(x)[1], "\".")
  x
}

#' @describeIn T2Afr Default method for generic function
#'
#' @param Rfr numeric vector. Spectral reflectance o reflectance factor.
#'   Set to zero if \code{x} is internal reflectance,
#' @export
#'
T2Afr.numeric <- function(x,
                          action = NULL,
                          byref = FALSE,
                          clean = FALSE,
                          Rfr = NA_real_,
                          ...) {
  if (byref) {
    stop("Conversion by reference not supported for \"numeric\" objects.")
  }
  if (is.na(Rfr)) {
    warning("Convertion requires 'Rfr' to be known.")
  }
  if (any(Rfr > 1) || any(Rfr < 0)) {
    warning("Bad 'Rfr' input valies.")
  }
  if (any(x > 1) || any(x < 0)) {
    warning("Bad 'Tfr' input valies.")
  }
  Tfr.internal <- x * (1 - Rfr)
  1 - Tfr.internal
}

#' @describeIn T2Afr Method for filter spectra
#'
#' @export
#'
T2Afr.filter_spct <- function(x,
                              action = "add",
                              byref = FALSE,
                              clean = FALSE,
                              ...) {
  if (byref) {
    name <- substitute(x)
  }

  action <- check_action_arg(x, action)
  # nothing to do
  if (action == "replace" &&
      "Afr" %in% colnames(x) && !any(c("A", "Tfr") %in% colnames(x))) {
    return(x)
  }
  if (action == "add" && "Afr" %in% colnames(x)) {
    return(x)
  }
  # we remove normalization
  ## should skip this when only removing s.e.irrad by removing col from normalization
  if (is_normalised(x)) {
    old.normalization.ls <- getNormalization(x)
    x <- denormalize_spct(x) # remove normalization
    norm.action <- "update"
  } else {
    norm.action <- "skip"
  }

  if (exists("Tfr", x, inherits = FALSE)) {
    if (clean) {
      x <- using_Tfr(clean(x))
    }
    current.Tfr.type <- getTfrType(x)
    properties <- getFilterProperties(x)

    if (current.Tfr.type == "total") {
      if (exists("Rfr", x, inherits = FALSE)) {
        x[["Afr"]] <- 1 - x[["Tfr"]] - x[["Rfr"]]
      } else if (!is.na(properties[["Rfr.constant"]])) {
        x[["Afr"]] <- 1 - x[["Tfr"]] - properties[["Rfr.constant"]]
      } else {
        x[["Afr"]] <- NA_real_
        action <- "add" # avoid loss of information
        warning("Conversion from total Tfr to Afr possible only ",
                "if Rfr or Rfr.constant are known.")
      }
    } else if (current.Tfr.type == "internal") {
      x[["Afr"]] <- 1 - x[["Tfr"]]
    } else {
      stop("Missing or bad 'Tfr.type' attribute setting.")
    }
  } else {
    x[["Afr"]] <- NA_real_
    action <- "add" # avoid loss of information
    warning("'Tfr' missing in 'T2Afr()', filling 'Afr' with 'NA'.")
  }

  if (action == "replace" && exists("A", x, inherits = FALSE)) {
    x[["A"]] <- NULL
  }
  if (action == "replace" && exists("Tfr", x, inherits = FALSE)) {
    x[["Tfr"]] <- NULL
  }

  if (current.Tfr.type == "total") {
    if (action == "add") {
      x <- convertTfrType(x, Tfr.type = "total")
    } else {
      # no Tfr stored in object, but keep for future conversion operations
      x <- setTfrType(x, "total")
    }
  }

  if (norm.action == "update") {
    x <- restore_normalization(x, old.normalization.ls)
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  x
}

#' @describeIn T2Afr Method for object spectra
#'
#' @export
#'
T2Afr.object_spct <- T2Afr.filter_spct

#' @describeIn T2Afr Method for collections of filter spectra
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
T2Afr.filter_mspct <- function(x,
                               action = "add",
                               byref = FALSE,
                               clean = FALSE,
                               ...,
                               .parallel = FALSE,
                               .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          .fun = T2Afr.filter_spct,
          action = action,
          byref = byref,
          clean = FALSE,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn T2Afr Method for collections of object spectra
#'
#' @export
#'
T2Afr.object_mspct <- T2Afr.filter_mspct

# Afr2T ---------------------------------------------------------------------

#' Convert transmittance into absorptance.
#'
#' Function that converts transmittance (fraction) into absorptance (fraction).
#' If reflectance (fraction) is available, it allows conversions between
#' internal and total absorptance.
#'
#' @param x an R object
#' @param action character Allowed values "replace" and "add"
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @param clean logical replace off-boundary values before conversion
#' @param ... not used in current version
#'
#' @return A copy of \code{x} with a column \code{Tfr} added and other columns
#'   possibly deleted except for \code{w.length}. If \code{action = "replace"},
#'   in all cases, the additional columns are removed, even if no column needs
#'   to be added.
#'
#' @export
#'
#' @family quantity conversion functions
#'
#' @examples
#' T2Afr(Ler_leaf.spct)
#'
Afr2T <- function(x, action, byref, clean, ...) UseMethod("Afr2T")

#' @describeIn Afr2T Default method for generic function
#'
#' @export
#'
Afr2T.default <- function(x,
                          action = NULL,
                          byref = FALSE,
                          clean = FALSE,
                          ...) {
  warning("'Afr2T()' not implemented for class \"", class(x)[1], "\".")
  x
}

#' @describeIn Afr2T Default method for generic function
#'
#' @param Rfr numeric vector. Spectral reflectance o reflectance factor.
#'   Set to zero if \code{x} is internal reflectance,
#' @export
#'
Afr2T.numeric <- function(x,
                          action = NULL,
                          byref = FALSE,
                          clean = FALSE,
                          Rfr = NA_real_,
                          ...) {
  if (byref) {
    stop("Conversion by reference not supported for \"numeric\" objects.")
  }
  if (is.na(Rfr)) {
    warning("Convertion requires 'Rfr'.")
  }
  if (any(Rfr > 1) || any(Rfr < 0)) {
    warning("Bad 'Rfr' input valies.")
  }
  if (any(x > 1) || any(x < 0)) {
    warning("Bad 'Afr' input valies.")
  }
  Afr.internal <- x / (1 - Rfr)
  1 - Afr.internal
}

#' @describeIn Afr2T Method for filter spectra
#'
#' @export
#'
Afr2T.filter_spct <- function(x,
                              action = "add",
                              byref = FALSE,
                              clean = FALSE,
                              ...) {
  if (byref) {
    name <- substitute(x)
  }

  action <- check_action_arg(x, action)
  # nothing to do
  if (action == "replace" &&
      "Tfr" %in% colnames(x) && !any(c("A", "Afr") %in% colnames(x))) {
    return(x)
  }
  if (action == "add" && "Tfr" %in% colnames(x)) {
    return(x)
  }
  # we remove normalization
  ## should skip this when only removing s.e.irrad by removing col from normalization
  if (is_normalised(x)) {
    old.normalization.ls <- getNormalization(x)
    x <- denormalize_spct(x) # remove normalization
    norm.action <- "update"
  } else {
    norm.action <- "skip"
  }

  if (exists("Afr", x, inherits = FALSE)) {
    if (clean) {
      x <- using_Afr(clean(x))
    }
    current.Tfr.type <- getTfrType(x)
    properties <- getFilterProperties(x)

    if (current.Tfr.type == "total") {
      if (exists("Rfr", x, inherits = FALSE)) {
        x[["Tfr"]] <- 1 - x[["Afr"]] - x[["Rfr"]]
      } else if (!is.na(properties[["Rfr.constant"]])) {
        properties <- getFilterProperties(x, return.null = FALSE)
        x[["Tfr"]] <- 1 - x[["Afr"]] - properties[["Rfr.constant"]]
      } else {
        x[["Tfr"]] <- NA_real_
        action <- "add" # avoid loss of information
        warning("Conversion from total Afr to Tfr possible only ",
                "if Rfr or Rfr.constant are known.")
      }
    } else if (current.Tfr.type == "internal") {
      x[["Tfr"]] <- 1 - x[["Afr"]]
    } else {
      stop("Missing or bad 'Tfr.type' attribute setting.")
    }
  } else {
    x[["Tfr"]] <- NA_real_
    action <- "add" # avoid loss of information
    warning("'Afr' missing in 'Afr2T()', filling 'Tfr' with 'NA'.")
  }

  if (action == "replace" && exists("Afr", x, inherits = FALSE)) {
    x[["Afr"]] <- NULL
  }
  if (action == "replace" && exists("A", x, inherits = FALSE)) {
    x[["A"]] <- NULL
  }

  if (norm.action == "update") {
    x <- restore_normalization(x, old.normalization.ls)
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  check_spct(x)
}

#' @describeIn Afr2T Method for object spectra
#'
#' @export
#'
Afr2T.object_spct <- Afr2T.filter_spct

#' @describeIn Afr2T Method for collections of filter spectra
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
Afr2T.filter_mspct <- function(x,
                               action = "add",
                               byref = FALSE,
                               clean = FALSE,
                               ...,
                               .parallel = FALSE,
                               .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          .fun = T2Afr.filter_spct,
          action = action,
          byref = byref,
          clean = FALSE,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn Afr2T Method for collections of object spectra
#'
#' @export
#'
Afr2T.object_mspct <- Afr2T.filter_mspct

# "any" filter conversions ------------------------------------------------

#' Convert filter quantities.
#'
#' Functions that convert or add related physical quantities to
#' \code{filter_spct} or  \code{object_spct} objects. transmittance (fraction)
#' into absorptance (fraction).
#'
#' @param x an filter_spct or a filter_mspct object.
#' @param action character Allowed values "replace" and "add".
#' @param clean logical replace off-boundary values before conversion
#'
#' @details These functions are dispatchers for \code{\link{A2T}},
#'   \code{\link{Afr2T}}, \code{\link{T2A}}, and \code{\link{T2Afr}}. The
#'   dispatch is based on the names of the variables stored in \code{x}. They
#'   do not support in-place modification of \code{x}.
#'
#' @return A copy of \code{x} with the columns for the different quantities
#'   added or replaced. If \code{action = "replace"}, in all cases, the
#'   additional columns are removed, even if no column needs to be added.
#'
#' @family quantity conversion functions
#'
#' @export
#'
#' @examples
#' any2Afr(Ler_leaf.spct)
#' any2T(Ler_leaf.spct)
#' any2T(polyester.spct)
#'
any2T <- function(x, action = "add", clean = FALSE) {
  if (is.filter_mspct(x) || is.object_mspct(x)) {
    if (!length(x)) return(x)
    return(msmsply(mspct = x,
                   .fun =  any2T,
                   action = action,
                   clean = clean))
  }
  stopifnot(is.filter_spct(x) || is.object_spct(x))
  if (any(c("A", "Tfr") %in% colnames(x))) {
    A2T(x, action = action, clean = clean, byref = FALSE)
  } else {
    Afr2T(x, action = action, clean = clean, byref = FALSE)
  }
}

#' @rdname any2T
#'
#' @export
#'
any2A <- function(x, action = "add", clean = FALSE) {
  if (is.filter_mspct(x) || is.object_mspct(x)) {
    if (!length(x)) return(x)
    return(msmsply(mspct = x,
                   .fun =  any2A,
                   action = action,
                   clean = clean))
  }
  stopifnot(is.filter_spct(x) || is.object_spct(x))
  if (any(c("A", "Tfr") %in% colnames(x))) {
    T2A(x, action = action, clean = clean, byref = FALSE)
  } else {
    Afr2T(x, action = action)
    T2A(x, action = action, clean = clean, byref = FALSE)
  }
}

#' @rdname any2T
#'
#' @export
#'
any2Afr <- function(x, action = "add", clean = FALSE) {
  if (is.filter_mspct(x) || is.object_mspct(x)) {
    if (!length(x)) return(x)
    return(msmsply(mspct = x,
                   .fun =  any2Afr,
                   action = action,
                   clean = clean))
  }
  stopifnot(is.filter_spct(x) || is.object_spct(x))
  if (any(c("Afr", "Tfr") %in% colnames(x))) {
    T2Afr(x, action = action, clean = clean, byref = FALSE)
  } else {
    A2T(x, action = action)
    T2Afr(x, action = action, clean = clean, byref = FALSE)
  }
}


# energy - photon and photon - energy conversions -------------------------

# energy to photon ---------------------------------------------------------------------


#' Convert energy-based quantities into photon-based quantities.
#'
#' Conversion methods for spectral energy irradiance into spectral photon
#' irradiance and for spectral energy response into spectral photon
#' response.
#'
#' @details The converted spectral values are added to or replace the existing
#'   spectral values depending on the argument passed to parameter
#'   \code{action}. Addition is currently not supported for normalized spectra.
#'   If the spectrum has been normalized with a recent version of package
#'   'photobiology' the spectrum will be renormalized after conversion using the
#'   same arguments as previously.
#'
#' @param x an R object.
#' @param action a character string, one of "add", or "replace".
#' @param byref logical indicating if a new object will be created by reference
#'   or a new object returned.
#' @param ... not used in current version.
#'
#' @export
#' @family quantity conversion functions
#'
e2q <- function(x, action, byref, ...) UseMethod("e2q")

#' @describeIn e2q Default method
#'
#' @export
#'
e2q.default <- function(x, action = "add", byref = FALSE, ...) {
  return(NA)
}

#' @describeIn e2q Method for spectral irradiance
#'
#' @export
#'
e2q.source_spct <- function(x,
                            action = NULL,
                            byref = FALSE,
                            ...) {

  if (byref) {
    name <- substitute(x)
  }

  action <- check_action_arg(x, action)
  # nothing to do
  if (action == "add" &&
      all(c("s.e.irrad", "s.q.irrad") %in% colnames(x))) {
    return(x)
  }
  if (action == "replace" &&
      "s.q.irrad" %in% colnames(x) && !"s.e.irrad" %in% colnames(x)) {
    return(x)
  }
  # we remove normalization
  ## should skip this when only removing s.e.irrad by removing col from normalization
  if (is_normalised(x)) {
    old.normalization.ls <- getNormalization(x)
    x <- denormalize_spct(x) # remove normalization
    norm.action <- "update"
  } else {
    norm.action <- "skip"
  }

  if (!exists("s.q.irrad", x, inherits = FALSE) &&
      exists("s.e.irrad", x, inherits = FALSE)) {
    x[["s.q.irrad"]] <-
      x[["s.e.irrad"]] * e2qmol_multipliers(x[["w.length"]])
  }

  if (action == "replace" &&
        exists("s.e.irrad", x, inherits = FALSE)) {
    x[["s.e.irrad"]] <- NULL
  }

  if (norm.action == "update") {
    x <- restore_normalization(x, old.normalization.ls)
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

#' @describeIn e2q Method for spectral responsiveness
#'
#' @export
#'
e2q.response_spct <- function(x,
                              action = "add",
                              byref = FALSE,
                              ...) {
  if (byref) {
    name <- substitute(x)
  }

  action <- check_action_arg(x, action)
  # nothing to do
  if (action == "add" &&
      all(c("s.e.response", "s.q.response") %in% colnames(x))) {
    return(x)
  }
  if (action == "replace" &&
      "s.q.response" %in% colnames(x) && !"s.e.response" %in% colnames(x)) {
    return(x)
  }
  # we remove normalization
  ## should skip this when only removing s.e.irrad by removing col from normalization
  if (is_normalised(x)) {
    old.normalization.ls <- getNormalization(x)
    x <- denormalize_spct(x) # remove normalization
    norm.action <- "update"
  } else {
    norm.action <- "skip"
  }

  if (!exists("s.q.response", x, inherits = FALSE) &&
      exists("s.e.response", x, inherits = FALSE)) {
    x[["s.q.response"]] <-
      x[["s.e.response"]] / e2qmol_multipliers(x[["w.length"]])
  }

  if (action == "replace" &&
      exists("s.e.response", x, inherits = FALSE)) {
    x[["s.e.response"]] <- NULL
  }

  if (norm.action == "update") {
    x <- restore_normalization(x, old.normalization.ls)
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

#' @describeIn e2q Method for collections of (light) source spectra
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
e2q.source_mspct <- function(x,
                             action = "add",
                             byref = FALSE,
                             ...,
                             .parallel = FALSE,
                             .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          .fun = e2q.source_spct,
          action = action,
          byref = byref,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn e2q Method for collections of response spectra
#'
#' @export
#'
e2q.response_mspct <- function(x,
                               action = "add",
                               byref = FALSE,
                               ...,
                               .parallel = FALSE,
                               .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          .fun = e2q.response_spct,
          action = action,
          byref = byref,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

# photon to energy ------------------------------------------------------------

#' Convert photon-based quantities into energy-based quantities
#'
#' Conversion methods for spectral photon irradiance into spectral energy
#' irradiance and for spectral photon response into spectral energy
#' response.
#'
#' @details The converted spectral values are added to or replace the existing
#'   spectral values depending on the argument passed to parameter
#'   \code{action}. Addition is currently not supported for normalized spectra.
#'   If the spectrum has been normalized with a recent version of package
#'   'photobiology' the spectrum will be renormalized after conversion using the
#'   same arguments as previously.
#'
#' @param x an R object.
#' @param action a character string, one of "add", or "replace".
#' @param byref logical indicating if a new object will be created by reference
#'   or a new object returned.
#' @param ... not used in current version.
#'
#' @export
#' @family quantity conversion functions
#'
q2e <- function(x, action, byref, ...) UseMethod("q2e")

#' @describeIn q2e Default method
#'
#' @export
#'
q2e.default <- function(x, action = "add", byref = FALSE, ...) {
  return(NA)
}

#' @describeIn q2e Method for spectral irradiance
#'
#' @export
#'
q2e.source_spct <- function(x,
                            action = "add",
                            byref = FALSE,
                            ...) {
  if (byref) {
    name <- substitute(x)
  }

  action <- check_action_arg(x, action)
  # nothing to do
  if (action == "add" &&
      all(c("s.e.irrad", "s.q.irrad") %in% colnames(x))) {
    return(x)
  }
  if (action == "replace" &&
      "s.e.irrad" %in% colnames(x) && !"s.q.irrad" %in% colnames(x)) {
    return(x)
  }
  # we remove normalization
  ## should skip this when only removing s.e.irrad by removing col from normalization
  if (is_normalised(x)) {
    old.normalization.ls <- getNormalization(x)
    x <- denormalize_spct(x) # remove normalization
    norm.action <- "update"
  } else {
    norm.action <- "skip"
  }

  if (!exists("s.e.irrad", x, inherits = FALSE) &&
      exists("s.q.irrad", x, inherits = FALSE)) {
    x[["s.e.irrad"]] <-
      x[["s.q.irrad"]] / e2qmol_multipliers(x[["w.length"]])
  }

  if (action == "replace" &&
      exists("s.q.irrad", x, inherits = FALSE)) {
    x[["s.q.irrad"]] <- NULL
  }

  if (norm.action == "update") {
    x <- restore_normalization(x, old.normalization.ls)
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

#' @describeIn q2e Method for spectral responsiveness
#'
#' @export
#'
q2e.response_spct <- function(x,
                              action = "add",
                              byref = FALSE,
                              ...) {
  if (byref) {
    name <- substitute(x)
  }

  action <- check_action_arg(x, action)
  # nothing to do
  if (action == "add" &&
      all(c("s.e.response", "s.q.response") %in% colnames(x))) {
    return(x)
  }
  if (action == "replace" &&
      "s.e.response" %in% colnames(x) && !"s.q.response" %in% colnames(x)) {
    return(x)
  }
  # we remove normalization
  ## should skip this when only removing s.e.irrad by removing col from normalization
  if (is_normalised(x)) {
    old.normalization.ls <- getNormalization(x)
    x <- denormalize_spct(x) # remove normalization
    norm.action <- "update"
  } else {
    norm.action <- "skip"
  }

  if (!exists("s.e.response", x, inherits = FALSE) &&
      exists("s.q.response", x, inherits = FALSE)) {
    x[["s.e.response"]] <-
      x[["s.q.response"]] * e2qmol_multipliers(x[["w.length"]])
  }

  if (action == "replace" &&
      exists("s.q.response", x, inherits = FALSE)) {
    x[["s.q.response"]] <- NULL
  }

  if (norm.action == "update") {
    x <- restore_normalization(x, old.normalization.ls)
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

#' @describeIn q2e Method for collections of (light) source spectra
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
q2e.source_mspct <- function(x,
                             action = "add",
                             byref = FALSE,
                             ...,
                             .parallel = FALSE,
                             .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          q2e.source_spct,
          action = action,
          byref = byref,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}


#' @describeIn q2e Method for collections of response spectra
#'
#' @export
#'
q2e.response_mspct <- function(x,
                               action = "add",
                               byref = FALSE,
                               ...,
                               .parallel = FALSE,
                               .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          q2e.response_spct,
          action = action,
          byref = byref,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @keywords internal
#'
check_action_arg <-
  function(x, action, accepted = c("add", "replace")) {
    if (is.null(action) || is.na(action[1])) {
      if (is_normalized(x)) {
        action <- "replace"
      } else {
        action <- "add"
      }
    } else if (!is.character(action)) {
      stop("Argument to 'action' is \"", class(action)[1], "\" but one of \"",
           paste(accepted, collapse = "\", \""), "\" expected.")
    } else if (!action[1] %in% accepted) {
      stop("'action = \"", action, "\"' but one of \"",
           paste(accepted, collapse = "\", \""), "\" expected.")
    } else {
      action <- action[1]
    }
    action
  }
