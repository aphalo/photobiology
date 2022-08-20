
# names of all spectral classes -------------------------------------------

#' Function returning a vector containing the names of spectra classes.
#'
#' @export
#'
#' @return A \code{character} vector of class names.
#' @examples
#' spct_classes()
#'
spct_classes <- function() {
  c("calibration_spct",
    "raw_spct", "cps_spct",
    "filter_spct", "reflector_spct",
    "source_spct", "object_spct",
    "response_spct", "chroma_spct",
    "solute_spct", "generic_spct")
}

# check -------------------------------------------------------------------

#' Enable or disable checks
#'
#' Choose between protection against errors or faster performance by enabling
#' (the default) or disabling data-consistency and sanity checks.
#'
#' @family data validity check functions
#'
#' @return The previous value of the option, which can be passed as argument
#'   to function \code{set_check_spct()} to restore the previous state of the
#'   option.
#'
#' @export
#'
enable_check_spct <- function() {
  options(photobiology.check.spct = TRUE)[[1]]
}

#' @rdname enable_check_spct
#'
#' @export
#'
disable_check_spct <- function() {
  options(photobiology.check.spct = FALSE)[[1]]
}

#' @rdname enable_check_spct
#'
#' @param x logical Flag to enable (TRUE), disable (FALSE) or unset (NULL)
#'   option.
#'
#' @export
#'
set_check_spct <- function(x) {
  options(photobiology.check.spct = x)[[1]]
}

#' Check validity of spectral objects
#'
#' Check that an R object contains the expected data members.
#'
#' @param x An R object
#' @param byref logical indicating if new object will be created by reference or
#'   by copy of \code{x}
#' @param strict.range logical indicating whether off-range values result in an
#'   error instead of a warning, \code{NA} disables the test.
#' @param force logical If \code{TRUE} check is done even if checks are
#'   disabled.
#' @param ... additional param possible in derived methods
#'
#' @export
#' @examples
#' check_spct(sun.spct)
#'
#' @family data validity check functions
#'
#' @examples
#' check_spct(sun.spct)
#' # try(check_spct(-sun.spct))
#' # try(check_spct((sun.spct[1, "w.length"] <- 1000)))
#'
check_spct <- function(x, byref, strict.range, force = FALSE, ...) {
  if (force || getOption("photobiology.check.spct", TRUE)) {
    UseMethod("check_spct")
  } else {
    x
  }
}

#' @describeIn check_spct Default for generic function.
#' @export
check_spct.default <-
  function(x, byref = FALSE, strict.range = NA, force = FALSE, ...) {
  x
}

#' @describeIn check_spct Specialization for generic_spct.
#'
#' @param multiple.wl numeric Maximum number of repeated w.length entries with
#'   same value.
#'
#' @export
check_spct.generic_spct <-
  function(x,
           byref = TRUE,
           strict.range = NA,
           force = FALSE,
           multiple.wl = getMultipleWl(x),
           ...)
  {
    # assert that option is set so that we can keep remaining code simpler.
    # defensive code in case the option has been unset by the user
    if (is.null(getOption("photobiology.verbose"))) {
      options(photobiology.verbose = getOption("verbose"))
    }

    # fix old class attributes
    class.x <- class_spct(x)
    if (!("tbl_df") %in% class(x)) {
      x <- tibble::as_tibble(x)
    }
    class(x) <- union(class.x, class(x)) # can change order!! BUG PRONE
    # check variables
    if (exists("wl", x, mode = "numeric", inherits = FALSE)) {
      dots <- list(~wl)
      x <- dplyr::rename_(x, .dots = stats::setNames(dots, "w.length"))
    } else if (exists("wavelength", x, mode = "numeric", inherits = FALSE)) {
      dots <- list(~wavelength)
      x <- dplyr::rename_(x, .dots = stats::setNames(dots, "w.length"))
    } else if (exists("Wavelength", x, mode = "numeric", inherits = FALSE)) {
      dots <- list(~Wavelength)
      x <- dplyr::rename_(x, .dots = stats::setNames(dots, "w.length"))
    }
    if (!exists("w.length", x, mode = "numeric", inherits = FALSE)) {
      stop("No wavelength data found in generic_spct")
    }

    if (nrow(x)) {
      wl.min <- min(x[["w.length"]], na.rm = TRUE)
      #  wl.max <- max(x[["w.length"]], na.rm = TRUE)
      if (wl.min == Inf) {
        warning("No valid 'w.length' values found") # could be stop()
      } else if (wl.min < 1) {
        stop("Negative or zero 'w.length' values found: aborting!")
      } else if ((wl.min < 99.999 || wl.min > 2.8e3) &&
                 getOption("photobiology.verbose")) { # catch use of Angstrom
        warning("Possibly off-range w.length values, minimum = ", signif(wl.min, 4), " nm. (Nanometers expected.)")
      }
      # we use run length encoding to find the maximum number of copies of any w.length value
      # this be needed. This redundancy needs to be fixed.
      if (multiple.wl == 1) {
        if (is.unsorted(x[["w.length"]], na.rm = TRUE, strictly = TRUE)) {
          if (is.unsorted(-x[["w.length"]], na.rm = TRUE, strictly = TRUE)) {
            stop("'w.length' must be sorted and have unique values")
          } else {
            # if unsorted is TRUE, then nrow >= 1 is TRUE
            # w.length in decreasing order, which we reverse
            x <- x[nrow(x):1, ]
          }
        }
      } else if (multiple.wl > 1) {
        runs <- rle(sort(x[["w.length"]]))
        num.copies <- max(runs[["lengths"]])
        if (num.copies > multiple.wl) {
          stop("Too many copies of w.length values: ", num.copies)
        }
      } else {
        stop("ASSERTION FAILED: invalid 'multiple.wl' value: ", multiple.wl)
      }
    }
    x
  }

#' @describeIn check_spct Specialization for calibration_spct.
#' @export
check_spct.calibration_spct <-
  function(x,
           byref = TRUE,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           force = FALSE,
           multiple.wl = getMultipleWl(x),
           ...) {

  x <- check_spct.generic_spct(x, multiple.wl = multiple.wl)

  mult.cols <- grep("^irrad.mult$", names(x))

  if (length(mult.cols) == 1 &&
      is.numeric(x[["irrad.mult"]]) &&
      all(na.omit(x[["irrad.mult"]]) >= 0)) {
    if (getOption("photobiology.verbose") && anyNA(x[["irrad.mult"]])) {
      warning("At least one NA in 'irrad.mult'")
    }
    return(x)
  } else {
    warning("No valid 'irrad.mult' data found in calibration_spct")
    x[["irrad.mult"]] = NA_real_
    return(x)
  }
}

#' @describeIn check_spct Specialization for raw_spct.
#' @export
check_spct.raw_spct <-
  function(x,
           byref = TRUE,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           force = FALSE,
           multiple.wl = getMultipleWl(x),
           ...) {

  x <- check_spct.generic_spct(x, multiple.wl = multiple.wl)

  counts.cols <- grep("^counts", names(x))

  if (length(counts.cols) >= 1) {
    return(x)
  } else {
    warning("No raw 'counts' data found in raw_spct")
    x[["counts"]] = NA_real_
    return(x)
  }
}

#' @describeIn check_spct Specialization for cps_spct.
#' @export
check_spct.cps_spct <-
  function(x,
           byref = TRUE,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           force = FALSE,
           multiple.wl = getMultipleWl(x),
           ...) {

  range_check <- function(x, cps.cols) {
    for (col in cps.cols) {
      stopifnot(is.numeric(x[[col]]))
      if (all(is.na(x[[col]]))) {
        next()
      }
      # we need to include zero and a reasonably high number as otherwise dark scans may not pass the test
      cps.range <- range(0, x[[col]], 1e3, na.rm = TRUE)
      # we need to be very lax here as during processing of scans we can get negative
      # values due to subtraction of dark scans
      if (cps.range[2] < 0 || abs(cps.range[1]) > (cps.range[2] / 10)) {
        message.text <- paste("Possible off-range cps values [",
                              formatted_range(cps.range),
                              "]", sep = "")
        if (is.null(strict.range) || is.na(strict.range)) {
          message(message.text)
        } else if (strict.range) {
          stop(message.text)
        } else if (!strict.range) {
          warning(message.text)
        } else {
          stop ("Bad argument for 'strict.range': ", strict.range)
        }
      }
    }
  }

  x <- check_spct.generic_spct(x, multiple.wl = multiple.wl)

  cps.cols <- grep("^cps", names(x), value = TRUE)

  if (length(cps.cols) >= 1) {
    range_check(x, cps.cols)
    return(x)
  } else {
    warning("No counts per second data found in cps_spct")
    x[["cps_1"]] <- NA_real_
    return(x)
  }
}

#' @describeIn check_spct Specialization for filter_spct.
#' @export
check_spct.filter_spct <-
  function(x,
           byref = TRUE,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           force = FALSE,
           multiple.wl = getMultipleWl(x),
           ...)
  {

    range_check_Tfr <- function(x, strict.range) {
      if (!all(is.na(x[["Tfr"]]))) {
        Tfr.min <- min(x[["Tfr"]], na.rm = TRUE)
        Tfr.max <- max(x[["Tfr"]], na.rm = TRUE)
        if (Tfr.min < -1e-4 || Tfr.max > 1 + 1e-6) {
          message.text <- paste("Off-range transmittance values [",
                                formatted_range(c(Tfr.min, Tfr.max)),
                                "] instead of  [0..1]", sep = "")
          if (is.null(strict.range) || is.na(strict.range)) {
            message(message.text)
          } else if (strict.range) {
            stop(message.text)
          } else if (!strict.range) {
            warning(message.text)
          } else {
            stop ("Bad argument for 'strict.range': ", strict.range)
          }
        }
      }
    }

    range_check_Afr <- function(x, strict.range) {
      if (!all(is.na(x[["Afr"]]))) {
        Afr.min <- min(x[["Afr"]], na.rm = TRUE)
        Afr.max <- max(x[["Afr"]], na.rm = TRUE)
        if (Afr.min < -1e-4 || Afr.max > 1 + 1e-6) {
          message.text <-
            paste0(
              "Off-range absorptance values [",
              formatted_range(c(Afr.min, Afr.max)),
              "] instead of  [0..1]",
              sep = ""
            )
          if (is.null(strict.range) || is.na(strict.range)) {
            message(message.text)
          } else if (strict.range) {
            stop(message.text)
          } else if (!strict.range) {
            warning(message.text)
          } else {
            stop ("Bad argument for 'strict.range': ", strict.range)
          }
        }
      }
    }

    range_check_A <- function(x, strict.range) {
      if (!all(is.na(x[["A"]]))) {
        A.min <- min(x[["A"]], na.rm = TRUE)
        A.max <- max(x[["A"]], na.rm = TRUE)
        if (A.min < -1e-7 || A.max > 20) {
          message.text <- paste("Off-range absorbance values [",
                                formatted_range(c(A.min, A.max)),
                                "] instead of  [0..20]", sep = "")
          if (is.null(strict.range) || is.na(strict.range)) {
            message(message.text)
          } else if (strict.range) {
            stop(message.text)
          } else if (!strict.range) {
            warning(message.text)
          } else {
            stop ("Bad argument for 'strict.range': ", strict.range)
          }
        }
      }
    }

    x <- check_spct.generic_spct(x, multiple.wl = multiple.wl)

    if (is.null(getTfrType(x))) {
      setTfrType(x, "total")
      warning("Missing Tfr.type attribute replaced by 'total'")
    }
    if (is.null(getTfrType(x))) {
      setTfrType(x, "total")
      warning("Missing Tfr.type attribute replaced by 'total'")
    }
    # check and replace 'other' quantity names
    if (exists("transmittance", x, mode = "numeric", inherits = FALSE)) {
      dots <- list(~transmittance)
      x <- dplyr::rename_(x, .dots = stats::setNames(dots, "Tpc"))
      warning("Found variable 'transmittance', I am assuming it is expressed as percent")
    }
    if (exists("absorbance", x, mode = "numeric", inherits = FALSE)) {
      dots <- list(~absorbance)
      x <- dplyr::rename_(x, .dots = stats::setNames(dots, "A"))
      warning("Found variable 'absorbance', I am assuming it is in log10-based absorbance units")
    } else if (exists("Absorbance", x, mode = "numeric", inherits = FALSE)) {
      dots <- list(~Absorbance)
      x <- dplyr::rename_(x, .dots = stats::setNames(dots, "A"))
      warning("Found variable 'Absorbance', I am assuming it is in log10-based absorbance units")
    }
    # look for percentages and change them into fractions of one
    if (exists("Tfr", x, mode = "numeric", inherits = FALSE)) {
      range_check_Tfr(x, strict.range = strict.range)
    } else if (exists("Tpc", x, mode = "numeric", inherits = FALSE)) {
      x[["Tfr"]] <- x[["Tpc"]] / 100
      x[["Tpc"]] <-  NULL
      range_check_Tfr(x, strict.range = strict.range)
    } else if (exists("A", x, mode = "numeric", inherits = FALSE)) {
      range_check_A(x, strict.range = strict.range)
    } else if (exists("Afr", x, mode = "numeric", inherits = FALSE)) {
      range_check_Afr(x, strict.range = strict.range)
    } else {
      warning("No transmittance, absortance or absorbance data found in filter_spct")
      x[["Tfr"]] <- NA_real_
    }
    if (getOption("photobiology.verbose")) {
      if (exists("Tfr", x, mode = "numeric", inherits = FALSE) && anyNA(x[["Tfr"]])) {
        warning("At least one NA in 'Tfr'")
      }
      if (exists("A", x, mode = "numeric", inherits = FALSE) && anyNA(x[["A"]])) {
        warning("At least one NA in 'A'")
      }
    }
    x
  }

#' @describeIn check_spct Specialization for solute_spct.
#' @export
check_spct.solute_spct <-
  function(x,
           byref = TRUE,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           force = FALSE,
           multiple.wl = getMultipleWl(x),
           ...) {

    range_check <- function(x, strict.range, k.base) {
      K.col <- paste("K", k.base, sep = ".")
      stopifnot(K.col %in% names(x))

      if (!all(is.na(x[[K.col]]))) {
        K.col.min <- min(x[[K.col]], na.rm = TRUE)
        K.col.max <- max(x[[K.col]], na.rm = TRUE)
        if (K.col.min < -1e-4) {
          message.text <-
            paste0(
              "Negative attenuation values [",
              formatted_range(c(K.col.min, K.col.max)),
              "]",
              sep = ""
            )
          if (is.null(strict.range) || is.na(strict.range)) {
            message(message.text)
          } else if (strict.range) {
            stop(message.text)
          } else if (!strict.range) {
            warning(message.text)
          } else {
            stop ("Bad argument for 'strict.range': ", strict.range)
          }
        }
      }
    }

    x <- check_spct.generic_spct(x, multiple.wl = multiple.wl)
    if (is.null(getKType(x))) {
      setKType(x, "attenaution")
      warning("Missing K.type attribute replaced by 'attenuation'")
    }

    if (exists("K.mole", x, mode = "numeric", inherits = FALSE)) {
      range_check(x, strict.range = strict.range, k.base = "mole")
      if (getOption("photobiology.verbose") && anyNA(x[["K.mole"]])) {
        warning("At least one NA in 'K.mole'")
      }
    } else if (exists("K.mass", x, mode = "numeric", inherits = FALSE)) {
      range_check(x, strict.range = strict.range, k.base = "mass")
      if (getOption("photobiology.verbose") && anyNA(x[["K.mass"]])) {
        warning("At least one NA in 'K.mass'")
      }
    } else {
      warning("No coefficient of attenuation data found in solute_spct")
      x[["K.mole"]] <- NA_real_
    }
    x
  }

#' @describeIn check_spct Specialization for reflector_spct.
#' @export
check_spct.reflector_spct <-
  function(x,
           byref = TRUE,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           force = FALSE,
           multiple.wl = getMultipleWl(x),
           ...) {

    range_check <- function(x, strict.range) {
      if (!all(is.na(x[["Rfr"]]))) {
        Rfr.min <- min(x[["Rfr"]], na.rm = TRUE)
        Rfr.max <- max(x[["Rfr"]], na.rm = TRUE)
        if (Rfr.min < -1e-4 ||  Rfr.max > 1 + 1e-6) {
          message.text <-
            paste0(
              "Off-range reflectance values [",
              formatted_range(c(Rfr.min, Rfr.max)),
              "] instead of  [0..1]",
              sep = ""
            )
          if (is.null(strict.range) || is.na(strict.range)) {
            message(message.text)
          } else if (strict.range) {
            stop(message.text)
          } else if (!strict.range) {
            warning(message.text)
          } else {
            stop ("Bad argument for 'strict.range': ", strict.range)
          }
        }
      }
    }

    x <- check_spct.generic_spct(x, multiple.wl = multiple.wl)

    if (is.null(getRfrType(x))) {
      setRfrType(x, "total")
      warning("Missing Rfr.type attribute replaced by 'total'")
    }
    if (exists("reflectance", x, mode = "numeric", inherits=FALSE)) {
      dots <- list(~reflectance)
      x <- dplyr::rename_(x, .dots = stats::setNames(dots, "Rpc"))
      warning("Found variable 'reflectance', I am assuming it is expressed as percent")
    }
    if (exists("Rfr", x, mode = "numeric", inherits=FALSE)) {
      range_check(x, strict.range=strict.range)
    } else if (exists("Rpc", x, mode = "numeric", inherits=FALSE)) {
      x[["Rfr"]] <- x[["Rpc"]] / 100
      x[["Rpc"]] <- NULL
      range_check(x, strict.range=strict.range)
    } else {
      warning("No reflectance data found in reflector_spct")
      x[["Rfr"]] <- NA_real_
    }
    if (getOption("photobiology.verbose") &&
        exists("Rfr", x, mode = "numeric", inherits = FALSE) && anyNA(x[["Rfr"]])) {
      warning("At least one NA in 'Tfr'")
    }
    x
  }

#' @describeIn check_spct Specialization for object_spct.
#' @export

check_spct.object_spct <-
  function(x,
           byref = TRUE,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           force = FALSE,
           multiple.wl = getMultipleWl(x),
           ...) {

    range_check_Tfr <- function(x, strict.range) {
      if (!all(is.na(x[["Tfr"]]))) {
        Tfr.min <- min(x[["Tfr"]], na.rm = TRUE)
        Tfr.max <- max(x[["Tfr"]], na.rm = TRUE)
        if (Tfr.min < -1e-4 || Tfr.max > 1 + 1e-6) {
          message.text <-
            paste0(
              "Off-range transmittance values [",
              formatted_range(c(Tfr.min, Tfr.max)),
              "] instead of  [0..1]",
              sep = ""
            )
          if (is.null(strict.range) || is.na(strict.range)) {
            message(message.text)
          } else if (strict.range) {
            stop(message.text)
          } else if (!strict.range) {
            warning(message.text)
          } else {
            stop ("Bad argument for 'strict.range': ", strict.range)
          }
        }
      }
    }

    range_check_Afr <- function(x, strict.range) {
      if (!all(is.na(x[["Afr"]]))) {
        Afr.min <- min(x[["Afr"]], na.rm = TRUE)
        Afr.max <- max(x[["Afr"]], na.rm = TRUE)
        if (Afr.min < -1e-4 || Afr.max > 1 + 1e-6) {
          message.text <-
            paste0(
              "Off-range absorptance values [",
              formatted_range(c(Afr.min, Afr.max)),
              "] instead of  [0..1]",
              sep = ""
            )
          if (is.null(strict.range) || is.na(strict.range)) {
            message(message.text)
          } else if (strict.range) {
            stop(message.text)
          } else if (!strict.range) {
            warning(message.text)
          } else {
            stop ("Bad argument for 'strict.range': ", strict.range)
          }
        }
      }
    }

    range_check_Rfr <- function(x, strict.range) {
      if (!all(is.na(x[["Rfr"]]))) {
        Rfr.min <- min(x[["Rfr"]], na.rm = TRUE)
        Rfr.max <- max(x[["Rfr"]], na.rm = TRUE)
        if (!is.na(Rfr.min) && !is.na(Rfr.max)) {
          if (Rfr.min < -1e-4 ||  Rfr.max > 1 + 1e-6) {
            message.text <-
              paste0(
                "Off-range reflectance values [",
                formatted_range(c(Rfr.min, Rfr.max)),
                "] instead of  [0..1]",
                sep = ""
              )
            if (is.null(strict.range) || is.na(strict.range)) {
              message(message.text)
            } else if (strict.range) {
              stop(message.text)
            } else if (!strict.range) {
              warning(message.text)
            } else {
              stop ("Bad argument for 'strict.range': ", strict.range)
            }
          }
        }
      }
    }

    x <- check_spct.generic_spct(x, multiple.wl = multiple.wl)

    if (is.null(getTfrType(x))) {
      setTfrType(x, "total")
      warning("Missing Tfr.type attribute replaced by 'total'")
    }
    if (is.null(getRfrType(x))) {
      setRfrType(x, "total")
      warning("Missing Rfr.type attribute replaced by 'total'")
    }
    if (exists("reflectance", x, mode = "numeric", inherits=FALSE)) {
      dots <- list(~reflectance)
      x <- dplyr::rename_(x, .dots = stats::setNames(dots, "Rpc"))
      warning("Found variable 'reflectance', I am assuming it is expressed as percent")
    }
    if (exists("Rfr", x, mode = "numeric", inherits=FALSE)) {
      range_check_Rfr(x, strict.range=strict.range)
    } else if (exists("Rpc", x, mode = "numeric", inherits=FALSE)) {
      x[["Rfr"]] <- x[["Rpc"]] / 100
      x[["Rpc"]] <- NULL
      range_check_Rfr(x, strict.range=strict.range)
    } else {
      warning("No reflectance data found in object_spct")
      x[["Rfr"]] <- NA_real_
    }

    if (exists("transmittance", x, mode = "numeric", inherits=FALSE)) {
      dots <- list(~transmittance)
      x <- dplyr::rename_(x, .dots = stats::setNames(dots, "Tpc"))
      warning("Found variable 'transmittance', I am assuming it expressed as percent")
    }
    if (exists("Afr", x, mode = "numeric", inherits=FALSE)) {
      range_check_Afr(x, strict.range=strict.range)
    }

    if (exists("Tfr", x, mode = "numeric", inherits=FALSE)) {
      range_check_Tfr(x, strict.range=strict.range)
    } else if (exists("Tpc", x, mode = "numeric", inherits=FALSE)) {
      x[["Tfr"]] <- x[["Tpc"]] / 100
      x[["Tpc"]] <- NULL
      range_check_Tfr(x, strict.range=strict.range)
    }

    quantities <- colnames(x)
    if (!"Rfr" %in% quantities) {
      x[["Rfr"]] <- NA_real_
    }
    if (! any(c("Tfr", "Afr") %in% quantities)) {
      x[["Tfr"]] <- NA_real_
    }
    ### Creates an endless recursive call!!
    # if ("Afr" %in% quantities) {
    #   x <- Afr2T(x, action = "add")
    # }

    if (getOption("photobiology.verbose")) {
      if (exists("Tfr", x, mode = "numeric", inherits = FALSE) && anyNA(x[["Tfr"]])) {
        warning("At least one NA in 'Tfr'")
      }
      if (exists("Rfr", x, mode = "numeric", inherits = FALSE) && anyNA(x[["Rfr"]])) {
        warning("At least one NA in 'Rfr'")
      }
    }
    x
  }

#' @describeIn check_spct Specialization for response_spct.
#' @export
check_spct.response_spct <-
  function(x,
           byref = TRUE,
           strict.range = NA,
           force = FALSE,
           multiple.wl = getMultipleWl(x),
           ...) {

    x <- check_spct.generic_spct(x, multiple.wl = multiple.wl)

    x <- checkTimeUnit(x)

    if (exists("s.e.response", x, mode = "numeric", inherits=FALSE) ||
        exists("s.q.response", x, mode = "numeric", inherits=FALSE)) {
      NULL # nothing to do
    } else if (exists("response", x, mode = "numeric", inherits=FALSE)) {
      x[["s.e.response"]] <- x[["response"]]
      x[["response"]] <- NULL
      warning("Found variable 'response', I am assuming it is expressed on an energy basis")
    } else if (exists("signal", x, mode = "numeric", inherits=FALSE)) {
      x[["s.e.response"]] <- x[["signal"]]
      x[["signal"]] <- NULL
      warning("Found variable 'signal', I am assuming it is expressed on an energy basis")
    } else {
      warning("No response data found in response_spct")
      x[["s.e.response"]] <- NA_real_
      return(x)
    }
    if (getOption("photobiology.verbose")) {
      if (exists("s.e.response", x, mode = "numeric", inherits = FALSE) && anyNA(x[["s.e.response"]])) {
        warning("At least one NA in 's.e.response'")
      }
      if (exists("s.q.response", x, mode = "numeric", inherits = FALSE) && anyNA(x[["s.q.response"]])) {
        warning("At least one NA in 's.q.response'")
      }
    }
    x
  }

#' @describeIn check_spct Specialization for source_spct.
#' @export
check_spct.source_spct <-
  function(x,
           byref = TRUE,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           force = FALSE,
           multiple.wl = getMultipleWl(x),
           ...) {

    range_check <- function(x, strict.range) {
      min.limit <- -0.10 # we accept small negative values
      if (exists("s.e.irrad", x, inherits = FALSE) &&
          !all(is.na(x[["s.e.irrad"]]))) {
        s.e.range <- range(x[["s.e.irrad"]], na.rm = TRUE)
        s.e.spread <- s.e.range[2] # for irradiance zero is meaningful
        # we need to be fairly lax as dark reference spectra may have
        # proportionally lots of noise.
        if (s.e.range[1] < (min.limit * max(s.e.spread, 0.04) )) {
          message.text <-
            paste(
              "Negative spectral energy irradiance values; minimum s.e.irrad =",
              format(s.e.range[1], digits = 3, nsmall = 2)
            )
          if (is.null(strict.range) || is.na(strict.range)) {
            message(message.text)
          } else if (strict.range) {
            stop(message.text)
          } else if (!strict.range) {
            warning(message.text)
          } else {
            stop ("Bad argument for 'strict.range': ", strict.range)
          }
        }
      }
      if (exists("s.q.irrad", x, inherits = FALSE) &&
          !all(is.na(x[["s.q.irrad"]]))) {
        s.q.range <- range(x[["s.q.irrad"]], na.rm = TRUE)
        s.q.spread <- s.q.range[2] # zero is meaningful
        # we need to be fairly lax as dark reference spectra may have
        # proportionally lots of noise.
        if (s.q.range[1] < (min.limit * (max(s.q.spread, 1e-5)) )) {
          message.text <-
            paste(
              "Negative spectral photon irradiance values; minimum s.q.irrad =",
              signif(s.q.range[1], 2)
            )
          if (is.null(strict.range) || is.na(strict.range)) {
            message(message.text)
          } else if (strict.range) {
            stop(message.text)
          } else if (!strict.range) {
            warning(message.text)
          } else {
            stop ("Bad argument for 'strict.range': ", strict.range)
          }
        }
      }
    }

    x <- check_spct.generic_spct(x, multiple.wl = multiple.wl)
    x <- checkTimeUnit(x)

    if (is.null(is_effective(x))) {
      setBSWFUsed(x, "none")
      warning("Missing attribute 'bswf.used' set to 'none'")
    }
    if (exists("s.e.irrad", x, mode = "numeric", inherits=FALSE) ||
        exists("s.q.irrad", x, mode = "numeric", inherits=FALSE)) {
      NULL
    } else if (exists("irradiance", x, mode = "numeric", inherits=FALSE)) {
      x[["s.e.irradiance"]] <- x[["irradiance"]]
      x[["irradiance"]] <- NULL
      warning("Found variable 'irradiance', I am assuming it is expressed on an energy basis")
    } else {
      warning("No spectral irradiance data found in source_spct")
      x[["s.e.irrad"]] <- NA_real_
      return(x)
    }
    if (!is.null(strict.range) && !is.na(strict.range)) {
      range_check(x, strict.range = strict.range)
    }
    if (getOption("photobiology.verbose")) {
      if (exists("s.e.irrad", x, mode = "numeric", inherits = FALSE) && anyNA(x[["s.e.irrad"]])) {
        warning("At least one NA in 's.e.irrad'")
      }
      if (exists("s.q.irrad", x, mode = "numeric", inherits = FALSE) && anyNA(x[["s.q.irrad"]])) {
        warning("At least one NA in 's.q.irrad'")
      }
    }
    x
  }

#' @describeIn check_spct Specialization for chroma_spct.
#' @export

check_spct.chroma_spct <-
  function(x,
           byref = TRUE,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           force = FALSE,
           multiple.wl = getMultipleWl(x),
           ...) {

    names_x <- names(x)

    x <- check_spct.generic_spct(x, multiple.wl = multiple.wl)

    idxs <- grep("[XYZ]", names_x)
    names(x)[idxs] <- tolower(names_x[idxs])
    if (!exists("x", x, mode="numeric", inherits=FALSE)) {
      warning("Chromaticity coordinate 'x' data missing")
      x[["x"]] <- NA_real_
    }
    if (!exists("y", x, mode="numeric", inherits=FALSE)) {
      warning("Chromaticity coordinate 'y' data missing")
      x[["y"]] <- NA_real_
    }
    if (!exists("z", x, mode="numeric", inherits=FALSE)) {
      warning("Chromaticity coordinate 'z' data missing")
      x[["z"]] <- NA_real_
    }
    if (getOption("photobiology.verbose") && (anyNA(x[["x"]]) || anyNA(x[["y"]]) || anyNA(x[["z"]]))) {
      warning("One or more NAs in chromaticity coordinates")
    }
    return(x)
  }


# set class ---------------------------------------------------------------

#' Remove "generic_spct" and derived class attributes.
#'
#' Removes from a spectrum object the class attributes \code{"generic_spct"} and
#' any derived class attribute such as \code{"source_spct"}. \strong{This
#' operation is done by reference!}
#'
#' @param x an R object.
#' @param keep.classes character vector Names of classes to keep. Can be used
#'   to retain base class \code{"generic_spct"}.
#'
#' @export
#'
#' @note If \code{x} is an object of any of the spectral classes defined in this
#'   package, this function changes by reference the spectrum object into the
#'   underlying data.frame object. Otherwise, it just leaves \code{x} unchanged.
#'
#' @details This function alters \code{x} itself by reference. If \code{x} is
#'   not a \code{generic_spct} object, \code{x} is not modified. This function
#'   behaves similarly to \code{setdiff()} but preserving the original order of
#'   the character vector of the S3 class names.
#'
#' @return A character vector containing the removed class attribute values.
#'   This is different to the behaviour of function \code{unlist} in base R!
#'
#' @family set and unset spectral class functions
#'
#' @examples
#' my.spct <- sun.spct
#' removed <- rmDerivedSpct(my.spct)
#' removed
#' class(sun.spct)
#' class(my.spct)
#'
rmDerivedSpct <- function(x, keep.classes = NULL) {
  name <- substitute(x)
  allclasses <- class(x)
  classes2remove <- setdiff(spct_classes(), keep.classes)
  class(x) <- setdiff(allclasses, classes2remove)
  if (!is.generic_spct(x)) {
    attr(x, "spct.version") <- NULL
  }
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(setdiff(allclasses, class(x)))
}

#' Convert an R object into a spectrum object.
#'
#' Sets the class attribute of a data.frame or an object of a derived
#' class to "generic_spct".
#'
#' @param x data.frame, list or generic_spct and derived classes
#' @param multiple.wl numeric Maximum number of repeated \code{w.length} entries
#'   with same value.
#' @param idfactor character Name of factor distinguishing multiple spectra when
#'   stored longitudinally (required if \code{mulitple.wl} > 1).
#'
#' @export
#'
#' @return x
#'
#' @details This method alters \code{x} itself by reference and in addition
#'   returns the modified \code{x} invisibly. The wavelength values and data are
#'   checked for validity and out-of-range values trigger warnings. These checks
#'   are done during construction by means of the matching
#'   \code{\link{check_spct}} methods, unless checks have been disabled by
#'   setting the corresponding option (see \code{\link{enable_check_spct}}).
#'
#' @family set and unset spectral class functions
#'
#' @examples
#' my.df <- data.frame(w.length = 300:309, s.e.irrad = rep(100, 10))
#' is.source_spct(my.df)
#' setSourceSpct(my.df)
#' is.source_spct(my.df)
#'
setGenericSpct <-
  function(x,
           multiple.wl = 1L,
           idfactor = NULL) {
    name <- substitute(x)
    if (is.null(multiple.wl) && is.any_spct(x)) {
      multiple.wl <- attr(x, "multiple.wl", exact = TRUE)
    }
    if (is.null(idfactor) && is.any_spct(x)) {
      idfactor <- attr(x, "idfactor", exact = TRUE)
    }
    rmDerivedSpct(x)
    if (!is.data.frame(x) || inherits(x, "data.table")) {
      x <- tibble::as_tibble(x)
    }
    if (!is.generic_spct(x)) {
      class(x) <- c("generic_spct", class(x))
      attr(x, "spct.tags") <- NA
      x <- setMultipleWl(x, multiple.wl = multiple.wl)
    }
    x <- check_spct(x)
    attr(x, "idfactor") <- idfactor
    attr(x, "spct.version") <- 2
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

#' @describeIn setGenericSpct Set class of a an object to "calibration_spct".
#'
#' @export
#'
#'
setCalibrationSpct <-
  function(x,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           multiple.wl = 1L,
           idfactor = NULL) {
    name <- substitute(x)
    setGenericSpct(x, multiple.wl = multiple.wl, idfactor = idfactor)
    class(x) <- c("calibration_spct", class(x))
    x <- check_spct(x, strict.range = strict.range)
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

#' @describeIn setGenericSpct Set class of a an object to "raw_spct".
#'
#' @export
#'
#'
setRawSpct <-
  function(x,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           multiple.wl = 1L,
           idfactor = NULL) {
    name <- substitute(x)
    setGenericSpct(x, multiple.wl = multiple.wl, idfactor = idfactor)
    class(x) <- c("raw_spct", class(x))
    x <- check_spct(x, strict.range = strict.range)
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

#' @describeIn setGenericSpct Set class of a an object to "cps_spct".
#'
#' @export
#'
#'
setCpsSpct <-
  function(x,
           time.unit="second",
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           multiple.wl = 1L,
           idfactor = NULL) {
    name <- substitute(x)
    setGenericSpct(x, multiple.wl = multiple.wl, idfactor = idfactor)
    class(x) <- c("cps_spct", class(x))
    setTimeUnit(x, time.unit)
    x <- check_spct(x, strict.range = strict.range)
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

#' @describeIn setGenericSpct Set class of an object to "filter_spct".
#'
#' @param Tfr.type character Either "total" or "internal".
#' @param Rfr.constant numeric The value of the reflection factor [\eqn{/1}].
#' @param thickness numeric The thickness of the material.
#' @param attenuation.mode character One of \code{"reflection"},
#'   \code{"absorption"} or \code{"mixed"}.
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning.
#' @export
#'
#' @section Warning!: Not entering metadata when creating an object will limit
#'   the available operations!
#'
#' @note \code{"internal"} \strong{transmittance} is defined as the
#'   transmittance of the material body itself, while \code{"total"}
#'   transmittance includes the effects of surface reflectance on the amount of
#'   light transmitted. For non-diffusing materials like glass an approximate
#'   \code{Rfr.constant} value can be used to inter-convert total and internal
#'   transmittance values. Use \code{NA} if the the mode is not known, or not
#'   applicable, e.g., for materials subject to internal scattering. The
#'   validity of computations related to thickness of the material or length of
#'   the light path depends on the availability and accuracy of the metadata.
#'
setFilterSpct <-
  function(x,
           Tfr.type = c("total", "internal"),
           Rfr.constant = NA_real_,
           thickness = NA_real_,
           attenuation.mode = NA_character_,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           multiple.wl = 1L,
           idfactor = NULL) {
    name <- substitute(x)
    if ((is.object_spct(x) || is.filter_spct(x)) &&
               getTfrType(x) != "unknown") {
      if (length(Tfr.type) > 1) {
        Tfr.type <- getTfrType(x)
      } else if (Tfr.type != getTfrType(x)) {
        warning("Overwriting attribute 'Tfr.type' from ", getTfrType(x),
                " into ", Tfr.type)
      }
    }
    setGenericSpct(x, multiple.wl = multiple.wl, idfactor = idfactor)
    class(x) <- c("filter_spct", class(x))
    setTfrType(x, Tfr.type[1])
    setFilterProperties(x,
                        Rfr.constant = Rfr.constant,
                        thickness = thickness,
                        attenuation.mode = attenuation.mode)
    x <- check_spct(x, strict.range = strict.range)
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

#' @describeIn setGenericSpct Set class of an object to "solute_spct".
#'
#' @param K.type character A string, either "attenuation", "absorption" or
#'   "scattering".
#' @param name,solvent.name character The names of the substance and of the
#'   solvent. A named character vector, with member names such as "IUPAC" for
#'   the authority.
#' @param mass numeric The mass in Dalton (Da = g/mol).
#' @param formula character The molecular formula.
#' @param structure raster A bitmap of the structure.
#' @param ID,solvent.ID character The IDs of the substance and of the solvent. A
#'   named character vector, with member names such as "ChemSpider" or "PubChen"
#'   for the authority.
#'
#' @export
#'
#' @note Particles in suspension unlike disolved \strong{solutes} scatter light.
#'   Thus two different processes can attenuate light in liquid media:
#'   absorption and scattering. Coefficients of attenuation are always based on
#'   measurements of internal absorbance or internal transmittance. In practice
#'   this is achieved by using as reference pure solvent in a vessel, such as a
#'   spectrometer cuvette, called \emph{blank}. The measurement of the blank is
#'   done sequentially, before or after the \emph{sample} of interest in single
#'   beam spectrophotometers and concurrently in double beam spectrophotometers.
#'   \code{K.type} describes the process of attenuation: \code{"attenuation"},
#'   \code{"absorption"} or \code{"scattering"}, with \code{"attenuation"} used
#'   for cases of mixed modes of attenuation. Set \code{K.type = NA} if not
#'   available or unknown, or not applicable.
#'
setSoluteSpct <-
  function(x,
           K.type = c("attenuation", "absorption", "scattering"),
           name = NA_character_,
           mass = NA_character_,
           formula = NA_character_,
           structure = grDevices::as.raster(matrix()),
           ID = NA_character_,
           solvent.name = NA_character_,
           solvent.ID = NA_character_,
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           multiple.wl = 1L,
           idfactor = NULL) {
    name <- substitute(x)
    if (is.solute_spct(x) &&
        getKType(x) != "unknown") {
      if (length(K.type) > 1) {
        K.type <- getKType(x)
      } else if (K.type != getKType(x)) {
        warning("Overwriting attribute 'K.type' from ", getKType(x),
                " into ", K.type)
      }
    }
    setGenericSpct(x, multiple.wl = multiple.wl, idfactor = idfactor)
    class(x) <- c("solute_spct", class(x))
    setKType(x, K.type[1])
    setSoluteProperties(x,
                        pass.null = FALSE,
                        mass = mass,
                        formula = formula,
                        structure = structure,
                        name = name,
                        ID = ID,
                        solvent.name = solvent.name,
                        solvent.ID = solvent.ID)
    x <- check_spct(x, strict.range = strict.range)
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

#' @describeIn setGenericSpct Set class of a an object to "reflector_spct".
#'
#' @param Rfr.type character A string, either "total" or "specular".
#'
#' @export
#'
#' @note \code{"specular"} \strong{reflectance} is defined as that measured by
#'   collecting the light reflected by the surface at the \dQuote{mirror} of the
#'   angle of incidence; i.e., using a probe with a narrow angle of aperture.
#'   Usually measured close to normal angle of incidence. \code{"total"}
#'   \strong{reflectance} is defined as that measured by collecting all the
#'   light reflected by the surface; i.e., using an integrating sphere. In a
#'   mirror, reflectance is mostly specular, while on the white surface of a
#'   sheet of paper scattering predominates. In the first case the value for
#'   total reflectance is not much more than for specular reflectance, while in
#'   the second case the difference is much larger as the "specular" component
#'   is much smaller.
#'
setReflectorSpct <-
  function(x,
           Rfr.type = c("total", "specular"),
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           multiple.wl = 1L,
           idfactor = NULL) {
    name <- substitute(x)
    if ((is.object_spct(x) || is.reflector_spct(c)) && getRfrType(x) != "unknown") {
      if (length(Rfr.type) > 1) {
        Rfr.type <- getRfrType(x)
      } else if (Rfr.type != getRfrType(x)) {
        warning("Overwriting attribute 'Rfr.type' from ", getRfrType(x),
                " into ", Rfr.type)
      }
    }
    setGenericSpct(x, multiple.wl = multiple.wl, idfactor = idfactor)
    class(x) <- c("reflector_spct", class(x))
    setRfrType(x, Rfr.type[1])
    x <- check_spct(x, strict.range = strict.range)
    #  setkey_spct(x, w.length)
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

#' @describeIn setGenericSpct Set class of an object to "object_spct".
#'
#' @export
#'
setObjectSpct <-
  function(x,
           Tfr.type = c("total", "internal"),
           Rfr.type = c("total", "specular"),
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           multiple.wl = 1L,
           idfactor = NULL) {
    name <- substitute(x)
    if ((is.filter_spct(x) || is.object_spct(x)) && getTfrType(x) != "unknown") {
      if (length(Tfr.type) > 1) {
        Tfr.type <- getTfrType(x)
      } else if (Tfr.type != getTfrType(x)) {
        warning("Overwriting attribute 'Tfr.type' from ", getTfrType(x),
                " into ", Tfr.type)
      }
    } else {
      Tfr.type <- Tfr.type[1]
    }
    if ((is.reflector_spct(x) || is.object_spct(x)) && getRfrType(x) != "unknown") {
      if (length(Rfr.type) > 1) {
        Rfr.type <- getRfrType(x)
      } else if (Rfr.type != getRfrType(x)) {
        warning("Overwriting attribute 'Rfr.type' from ", getRfrType(x),
                " into ", Rfr.type)
      }
    } else {
      Rfr.type <- Rfr.type[1]
    }
    if (Tfr.type == "total" && Rfr.type != "total") {
      message("Rfr is not \"total\", making conversions between Afr and Tfr impossible.")
    }
    setGenericSpct(x, multiple.wl = multiple.wl, idfactor = idfactor)
    class(x) <- c("object_spct", class(x))
    setTfrType(x, Tfr.type)
    setRfrType(x, Rfr.type)
    x <- check_spct(x, strict.range = strict.range)
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

#' @describeIn setGenericSpct Set class of an object to "response_spct".
#'
#' @param time.unit character string indicating the time unit used for spectral
#'   irradiance or exposure (\code{"second"}, \code{"day"} or \code{"exposure"})
#'   or an object of class duration as defined in package lubridate.
#' @param response.type a character string, either \code{"response"} or \code{"action"}.
#' @export
#'
setResponseSpct <-
  function(x,
           time.unit = "second",
           response.type = "response",
           multiple.wl = 1L,
           idfactor = NULL) {
    name <- substitute(x)
    setGenericSpct(x, multiple.wl = multiple.wl, idfactor = idfactor)
    class(x) <- c("response_spct", class(x))
    setResponseType(x, response.type)
    setTimeUnit(x, time.unit)
    x <- check_spct(x)
    #  setkey_spct(x, w.length)
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

#' @describeIn setGenericSpct Set class of an object to "source_spct".
#'
#' @param bswf.used character A string, either \code{"none"} or the name of a
#'   BSWF. (Users seldom need to change the default, as this metadata value
#'   is in normal use set by operators or functions that apply a BSWF.)
#'
#' @export
#'
#'
setSourceSpct <-
  function(x,
           time.unit = "second",
           bswf.used = c("none", "unknown"),
           strict.range = getOption("photobiology.strict.range", default = FALSE),
           multiple.wl = 1L,
           idfactor = NULL) {
    name <- substitute(x)
    setGenericSpct(x, multiple.wl = multiple.wl, idfactor = idfactor)
    class(x) <- c("source_spct", class(x))
    setTimeUnit(x, time.unit)
    setBSWFUsed(x, bswf.used = bswf.used)
    x <- check_spct(x, strict.range = strict.range)
    #  setkey_spct(x, w.length)
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

#' @describeIn setGenericSpct Set class of an object to "chroma_spct".
#'
#' @export
#'
#'
setChromaSpct <-
  function(x,
           multiple.wl = 1L,
           idfactor = NULL) {
    name <- substitute(x)
    setGenericSpct(x, multiple.wl = multiple.wl, idfactor = idfactor)
    class(x) <- c("chroma_spct", class(x))
    x <- check_spct(x)
    #  setkey_spct(x, w.length)
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

# is functions for spct classes --------------------------------------------

#' Query class of spectrum objects
#'
#' Functions to query whether an object is of a given type of spectrum.
#'
#' @param x an R object.
#'
#' @return A logical value, \code{TRUE} if the argument passed to \code{x} is an
#'   object of the queried type of spectrum and \code{FALSE} otherwise.
#'
#' @note Derived types also return \code{TRUE} for a query for a base type such
#'   as \code{generic_spct}, following R's practice.
#'
#' @examples
#' is.source_spct(sun.spct)
#' is.filter_spct(sun.spct)
#' is.generic_spct(sun.spct)
#' is.generic_spct(sun.spct)
#'
#' @export
#'
#' @examples
#' is.source_spct(sun.spct)
#' is.filter_spct(sun.spct)
#' is.generic_spct(sun.spct)
#' is.generic_spct(sun.spct)
#'
is.generic_spct <- function(x) inherits(x, "generic_spct")

#' @rdname is.generic_spct
#' @export
#'
is.raw_spct <- function(x) inherits(x, "raw_spct")

#' @rdname is.generic_spct
#' @export
#'
is.calibration_spct <- function(x) inherits(x, "calibration_spct")

#' @rdname is.generic_spct
#' @export
#'
is.cps_spct <- function(x) inherits(x, "cps_spct")

#' @rdname is.generic_spct
#' @export
#'
is.source_spct <- function(x) inherits(x, "source_spct")

#' @rdname is.generic_spct
#' @export
#'
is.response_spct <- function(x) inherits(x, "response_spct")

#' @rdname is.generic_spct
#' @export
#'
is.filter_spct <- function(x) inherits(x, "filter_spct")

#' @rdname is.generic_spct
#' @export
#'
is.reflector_spct <- function(x) inherits(x, "reflector_spct")

#' @rdname is.generic_spct
#' @export
#'
is.object_spct <- function(x) inherits(x, "object_spct")

#' @rdname is.generic_spct
#' @export
#'
is.solute_spct <- function(x) inherits(x, "solute_spct")

#' @rdname is.generic_spct
#' @export
#'
is.chroma_spct <- function(x) inherits(x, "chroma_spct")

#' @rdname is.generic_spct
#'
#' @export
#'
is.any_spct <- function(x) {
  inherits(x, "generic_spct")
}

#' Query which is the class of a spectrum
#'
#' Extract class information from a generic spectrum.
#'
#' @details The value returned is equivalent to the set intersection of the
#'   value returned by \code{class(x)} and the value returned by
#'   \code{\link{spct_classes}}, but preserving the order of the members of
#'   the character vector.
#'
#' @param x any R object
#'
#' @return A character vector containing all matching xxxx.spct S3 classes.
#'
#' @export
#'
#' @examples
#' class_spct(sun.spct)
#' class(sun.spct)
#'
class_spct <- function(x) {
  #  intersect(spct_classes(), class(x)) # alters order!
  class(x)[class(x) %in% spct_classes()] # maintains order
}

#' Query if a spectrum is tagged
#'
#' Functions to check if an spct object contains tags.
#'
#' @param x any R object
#'
#' @return \code{is_tagged} returns a logical value, \code{TRUE} if its argument is a a spectrum that contains
#'   tags and \code{FALSE} if it is an untagged spectrum, but returns \code{NA}
#'   for any other R object.
#'
#' @export
#'
#' @family tagging and related functions
#' @examples
#' is_tagged(sun.spct)
#'
is_tagged <- function(x) {
  if (!is.generic_spct(x)) {
    return(NA)
  } else {
    tags <- attr(x, "spct.tags", exact=TRUE)
    return(!is.null(tags) && length(tags) > 0 && !is.na(tags[[1]]))
  }
}

# is_photon_based ---------------------------------------------------------

#' Query if a spectrum contains photon- or energy-based data.
#'
#' Functions to query if \code{source_spct} and \code{response_spct} objects
#' contain photon-based or energy-based data.
#'
#' @param x any R object
#'
#' @return \code{is_photon_based} returns a logical value, \code{TRUE} if its argument is a \code{source_spct}
#'   or a \code{response_spct} object that contains photon base data and
#'   \code{FALSE} otherwise, but returns \code{NA} for any other R object,
#'   including those belonging other \code{generic_spct}-derived classes.
#'
#' @export
#' @family query units functions
#'
#' @rdname is_photon_based
#' @examples
#' colnames(sun.spct)
#' is_photon_based(sun.spct)
#' my.spct <- sun.spct[ , c("w.length", "s.e.irrad")]
#' is.source_spct(my.spct)
#' is_photon_based(my.spct)
#'
is_photon_based <- function(x) {
  if (is.source_spct(x) || is.summary_source_spct(x)) {
    return("s.q.irrad" %in% names(x))
  } else if (is.response_spct(x)) {
    return("s.q.response" %in% names(x))
  } else {
    return(NA_integer_)
  }
}

# is_energy_based ---------------------------------------------------------

#' @rdname is_photon_based
#'
#' @return \code{is_energy_based} returns a logical value, \code{TRUE} if its argument is a \code{source_spct}
#'   or a \code{response_spct} object that contains energy base data and
#'   \code{FALSE} otherwise, but returns \code{NA} for any other R object,
#'   including those belonging other \code{generic_spct}-derived classes
#'
#' @export
#' @examples
#' colnames(sun.spct)
#' is_energy_based(sun.spct)
#' my.spct <- sun.spct[ , c("w.length", "s.q.irrad")]
#' is.source_spct(my.spct)
#' is_energy_based(my.spct)
#'
is_energy_based <- function(x) {
  if (is.source_spct(x) || is.summary_source_spct(x)) {
    return("s.e.irrad" %in% names(x))
  } else if (is.response_spct(x)) {
    return("s.e.response" %in% names(x))
  } else {
    return(NA_integer_)
  }
}

# is_absorbance_based ---------------------------------------------------------

#' Query if a spectrum contains absorbance or transmittance data
#'
#' Functions to query if an filter spectrum contains spectral absorbance data or
#' spectral transmittance data.
#'
#' @param x an R object
#'
#' @return \code{is_absorbance_based} returns a logical value, \code{TRUE} if its argument is a \code{filter_spct}
#'   object that contains spectral absorbance data and \code{FALSE} otherwise,
#'   but returns \code{NA} for any other R object, including those belonging
#'   other \code{generic_spct}-derived classes.
#'
#' @export
#' @family query units functions
#'
#' @rdname is_absorbance_based
#' @examples
#' is_absorbance_based(polyester.spct)
#' my.spct <- T2A(polyester.spct)
#' is.filter_spct(my.spct)
#' is_absorbance_based(my.spct)
#'
is_absorbance_based <- function(x) {
  if (is.filter_spct(x) || is.summary_filter_spct(x)) {
    return("A" %in% names(x))
  } else {
    return(NA_integer_)
  }
}

# is_absorptance_based ---------------------------------------------------------

#' @rdname is_absorbance_based
#'
#' @return \code{is_absorptance_based} returns a logical value, if its argument
#'   is a \code{filter_spct} object, \code{TRUE} if it contains data as spectral
#'   absorptance and \code{FALSE} otherwise, but returns \code{NA} for any other
#'   R object, including those belonging other \code{generic_spct}-derived
#'   classes.
#'
#' @export
#' @examples
#' is_absorptance_based(polyester.spct)
#'
is_absorptance_based <- function(x) {
  if (is.filter_spct(x) || is.summary_source_spct(x)) {
    return("Afr" %in% names(x))
  } else {
    return(NA_integer_)
  }
}

# is_transmittance_based ---------------------------------------------------------

#' @rdname is_absorbance_based
#'
#' @return \code{is_transmittance_based} returns TRUE if its argument is a
#'   \code{filter_spct} object that contains spectral transmittance data and
#'   FALSE if it does not contain such data, but returns NA for any other R
#'   object, including those belonging other \code{generic_spct}-derived
#'   classes.
#'
#' @export
#' @examples
#' is_transmittance_based(polyester.spct)
#'
is_transmittance_based <- function(x) {
  if (is.filter_spct(x) || is.summary_source_spct(x)) {
    return("Tfr" %in% names(x))
  } else {
    return(NA_integer_)
  }
}

# is_mole_based ---------------------------------------------------------

#' Query if a spectrum contains mole or mass based data
#'
#' Functions to check if an solute attenuation spectrum contains coefficients
#' on expressed on mole of mass base.
#'
#' @param x an R object
#'
#' @return \code{is_mole_based} returns TRUE if its argument is a
#'   \code{solute_spct} object that contains spectral \code{K.mole} data and
#'   \code{FALSE} if it contains \code{K.mass} data, but returns NA for any
#'   other R object, including those belonging other \code{generic_spct}-derived
#'   classes. \code{is_mass_based} returns the complement of
#'   \code{is_mole_based}.
#'
#' @export
#' @family query units functions
#'
#' @examples
#' print("missing example")
#'
is_mole_based <- function(x) {
  if (is.solute_spct(x) || is.summary_solute_spct(x)) {
    return("K.mole" %in% names(x))
  } else {
    return(NA_integer_)
  }
}

#' @rdname is_mole_based
#'
is_mass_based <- function(x) {
  !is_mole_based(x)
}

# time.unit attribute -----------------------------------------------------

#' Set the "time.unit" attribute of an existing source_spct object
#'
#' Function to set by reference the "time.unit" attribute
#'
#' @param x a source_spct object
#' @param time.unit character string indicating the time unit used for spectral
#'   irradiance or exposure ("second" , "day" or "exposure") or an object of
#'   class duration as defined in package lubridate.
#' @param override.ok logical Flag that can be used to silence warning when
#'   overwriting an existing attribute value (used internally)
#'
#' @return x
#' @note This function alters x itself by reference and in addition
#'   returns x invisibly. If x is not a source_spct or response_spct object, x is not modified.
#'   The behaviour of this function is 'unusual' in that the default for
#'   parameter \code{time.unit} is used only if \code{x} does not already have
#'   this attribute set. \code{time.unit = "hour"} is currently not fully
#'   supported.
#'
#' @export
#' @family time attribute functions
#' @examples
#' my.spct <- sun.spct
#' setTimeUnit(my.spct, time.unit = "second")
#' setTimeUnit(my.spct, time.unit = lubridate::duration(1, "seconds"))
#'
setTimeUnit <- function(x,
                        time.unit = c("second", "hour", "day", "exposure", "none"),
                        override.ok = FALSE) {
  if (!(class(x)[1] %in%
        c("source_spct", "summary_source_spct",
          "response_spct", "response_spct",
          "raw_spct", "cps_spct"))) {
     return(invisible(x))
  }

  if (is.character(time.unit)) {
    time.unit <- time.unit[1]
  }
  name <- substitute(x)
  old.time.unit <- getTimeUnit(x)
  override.ok <- ifelse(is.na(old.time.unit) ||
                          (is.character(old.time.unit) && old.time.unit == "unknown"),
                        TRUE, override.ok)
  if (override.ok) {
    if (is.character(time.unit)) {
      if (!(time.unit %in% c("second", "hour", "day", "none", "exposure", "unknown"))) {
        warning("Unrecognized 'time.unit' argument ", time.unit, " set to 'unknown'.")
        time.unit <- "unknown"
      }
    } else if (lubridate::is.duration(time.unit)) {
      if (time.unit <= lubridate::duration(0, "seconds")) {
        stop("When 'time.unit' is a duration, it must be > 0")
      }
    }
    attr(x, "time.unit") <- time.unit
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' Get the "time.unit" attribute of an existing source_spct object
#'
#' Function to read the "time.unit" attribute
#'
#' @param x a source_spct object
#' @param force.duration logical If TRUE a lubridate::duration is returned even
#'   if the object attribute is a character string, if no conversion is possible
#'   NA is returned.
#'
#' @return character string or a lubridate::duration
#'
#' @note if x is not a \code{source_spct} or a \code{response_spct} object, NA
#' is returned
#'
#' @export
#' @family time attribute functions
#' @examples
#' getTimeUnit(sun.spct)
#'
getTimeUnit <- function(x,
                        force.duration = FALSE) {
  if (class(x)[1] %in%
      c("source_spct", "summary_source_spct",
        "response_spct", "response_spct",
        "raw_spct", "cps_spct")) {
    time.unit <- attr(x, "time.unit", exact = TRUE)
    # need to handle objects created with old versions
    if (!length(time.unit)) {
      time.unit <- "unknown"
    }
    if (is.character(time.unit)) {
      time.unit <- time.unit[[1]]
    }
    # this is safe in case class attribute is lost, as duration is stored as seconds
    if (!lubridate::is.duration(time.unit) && is.numeric(time.unit)) {
      time.unit <- lubridate::duration(seconds = time.unit)
    }
    # convert strings to durations
    if (force.duration && is.character(time.unit)) {
      time.unit <- char2duration(time.unit)
    }
    return(time.unit)
  } else {
    if (force.duration) {
      lubridate::duration(NA_character_)
    } else {
      NA_character_
    }
  }
}

#' Convert the "time.unit" attribute of an existing source_spct object
#'
#' Function to set the "time.unit" attribute and simultaneously rescaling the
#' spectral data to be expressed using the new time unit as basis of expression.
#' The change is done by reference ('in place').
#'
#' @param x source_spct or response_spct object
#' @param time.unit a character string, either "second", "hour", "day",
#'   "exposure" or "none", or a lubridate::duration
#' @param ... (currently ignored)
#'
#' @return x possibly with the \code{time.unit} attribute modified
#'
#' @note if x is not a \code{source_spct} or a \code{response_spct} object, or
#'   time.unit is NULL x is returned unchanged, if the existing or new time.unit
#'   cannot be converted to a duration, then the returned spectrum will contain
#'   NAs.
#'
#' @export
#' @family time attribute functions
#' @examples
#'
#' my.spct <- sun.spct
#' my.spct
#' convertTimeUnit(my.spct, "day")
#' my.spct
#'
convertTimeUnit <- function(x, time.unit = NULL, ...) {
  if (!is.generic_spct(x)) {
    warning("'convertTimeUnit()' mot applicable to class '", class(x)[1], "'. Skipping!")
    return(invisible(x))
  }
  if (is.null(time.unit)) {
    # nothing to do
    return(invisible(x))
  }
  columns <- intersect(names(x), c("s.e.irrad", "s.q.irrad", "s.e.response", "s.q.response") )
  if (length(columns) == 0) {
    warning("No column to convert to new time unit.")
    return(invisible(x))
  }

  x.out <- checkTimeUnit(x)

  new.time.unit <- char2duration(time.unit)
  old.time.unit <- getTimeUnit(x.out, force.duration = TRUE)

  multiplier <- as.numeric(new.time.unit) / as.numeric(old.time.unit)

  for (col in columns) {
    x.out[[col]] <- x.out[[col]] * multiplier
  }

  setTimeUnit(x.out, time.unit, override.ok = TRUE)
}


#' Check the "time.unit" attribute of an existing source_spct object
#'
#' Function to read the "time.unit" attribute
#'
#' @param x a source_spct object
#'
#' @return x possibly with the \code{time.unit} attribute modified
#'
#' @note if x is not a \code{source_spct} or a \code{response_spct} object, NA
#' is returned
#'
#' @export
#' @family time attribute functions
#'
checkTimeUnit <- function(x) {
  if (is.source_spct(x) || is.response_spct(x) || is.cps_spct(x)) {
    time.unit <- getTimeUnit(x)
    ## Handled already in getTimeUnit()
    # if (!length(time.unit)) {
    #   setTimeUnit(x, "second")
    #   warning("Missing attribute 'time.unit' set to 'second'")
    # }

    if (is.character(time.unit)) {
      if (!(time.unit %in% c("second", "minute", "hour", "day", "exposure", "none", "unknown"))) {
        stop("'time.unit' ",  time.unit, " is unknown")
      }
    } else if (lubridate::is.duration(time.unit)) {
      if (time.unit <= lubridate::duration(0, "seconds")) {
        stop("When 'time.unit' is a duration, it must be > 0")
      }
    } else {
      stop("'time.unit' must be of class character or lubridate::duration, but found class '",
           class(time.unit), "' instead.")
    }
  }
  invisible(x)
}

# private
char2duration <- function(time.unit) {
  if (is.character(time.unit)) {
    time.duration <- switch(time.unit,
                            second  = lubridate::duration(1, "seconds"),
                            minute  = lubridate::duration(1, "minutes"),
                            hour    = lubridate::duration(1, "hours"),
                            day     = lubridate::duration(1, "days"),
                            exposure = lubridate::duration(NA_character_),
                            none    = lubridate::duration(NA_character_),
                            unknown = lubridate::duration(NA_character_)
    )
  } else if (lubridate::is.duration(time.unit)) {
    time.duration <- time.unit
  }
  return(time.duration)
}

# response.type attribute ------------------------------------------------------

#' Set the "response.type" attribute
#'
#' Function to set by reference the "response.type" attribute of an existing
#' response_spct object.
#'
#' Objects of class \code{response_spct()} can contain data for a response
#' spectrum or an action spectrum. Response spectra are measured using the
#' same photon (or energy) irradiance at each wavelength. Action spectra are
#' derived from dose response curves at each wavelength, and responsivity
#' at each wavelength is expressed as the reciprocal of the photon fluence
#' required to obtain a fixed level of response.
#'
#' @param x a response_spct object
#' @param response.type a character string, either "response" or "action"
#'
#' @return x
#' @note This function alters x itself by reference and in addition returns x
#'   invisibly. If x is not a response_spct object, x is not modified The
#'   behaviour of this function is 'unusual' in that the default for parameter
#'   \code{response.type} is used only if \code{x} does not already have this
#'   attribute set.
#'
#' @export
#' @family response type attribute functions
#' @examples
#' my.spct <- ccd.spct
#' setResponseType(my.spct, "action")
#'
setResponseType <- function(x,
                            response.type = c("response", "action")) {
  name <- substitute(x)
  if (length(response.type) > 1) {
    if (getResponseType(x) != "unknown") {
      response.type <- getResponseType(x)
    } else {
      response.type <- response.type[[1]]
    }
  }
  if (is.response_spct(x) || is.summary_response_spct(x)) {
    if  (!(response.type %in% c("response", "action", "unknown"))) {
      warning("Invalid 'response.type' argument, only 'response' and 'action' supported.")
      return(x)
    }
    attr(x, "response.type") <- response.type
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' Get the "response.type" attribute
#'
#' Function to read the "response.type" attribute of an existing response_spct
#' object.
#'
#' Objects of class \code{response_spct()} can contain data for a response
#' spectrum or an action spectrum. Response spectra are measured using the
#' same photon (or energy) irradiance at each wavelength. Action spectra are
#' derived from dose response curves at each wavelength, and responsivity
#' at each wavelength is expressed as the reciprocal of the photon fluence
#' required to obtain a fixed level of response.
#'
#' @param x a response_spct object
#'
#' @return character string
#'
#' @note If x is not a \code{response_spct} object, \code{NA} is returned.
#'
#' @export
#' @family response.type attribute functions
#' @examples
#' getResponseType(ccd.spct)
#' getResponseType(sun.spct)
#'
getResponseType <- function(x) {
  if (is.response_spct(x) || is.summary_response_spct(x)) {
    response.type <- attr(x, "response.type", exact = TRUE)
    if (is.null(response.type) || is.na(response.type)) {
      # need to handle objects created with old versions
      response.type <- "unknown"
    }
    return(response.type[[1]])
  } else {
    return(NA_character_)
  }
}

# bswf attribute -----------------------------------------------------

#' Set the "bswf.used" attribute
#'
#' Function to set by reference the "time.unit" attribute of an existing
#' source_spct object
#'
#' @param x a source_spct object
#' @param bswf.used a character string, either "none" or the name of a BSWF
#'
#' @return x
#' @note This function alters x itself by reference and in addition
#'   returns x invisibly. If x is not a source_spct, x is not modified. The behaviour of this
#'   function is 'unusual' in that the default for parameter \code{bswf.used} is
#'   used only if \code{x} does not already have this attribute set.
#'   \code{time.unit = "hour"} is currently not fully supported.
#'
#' @export
#' @family BSWF attribute functions
#'
setBSWFUsed <- function(x, bswf.used=c("none", "unknown")) {
  if (is.null(bswf.used) || length(bswf.used) < 1) {
    bswf.used <- "none"
  }
  if (length(bswf.used) > 1) {
    if (is_effective(x)) {
      bswf.used <- getBSWFUsed(x)
    } else {
      bswf.used <- bswf.used[[1]]
    }
  }
  if (is.source_spct(x) || is.summary_source_spct(x)) {
    name <- substitute(x)
    if  (!(is.character(bswf.used))) {
      warning("Only character strings are valid values for 'bswf.used' argument")
      bswf.used <- "unknown"
    }
    attr(x, "bswf.used") <- bswf.used
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' Get the "bswf.used" attribute
#'
#' Function to read the "time.unit" attribute of an existing source_spct object
#'
#' @param x a source_spct object
#'
#' @return character string
#'
#' @note if x is not a \code{source_spct} object, NA is returned
#'
#' @export
#' @family BSWF attribute functions
#' @examples
#' getBSWFUsed(sun.spct)
#'
getBSWFUsed <- function(x) {
  if (is.source_spct(x) || is.summary_source_spct(x)) {
    bswf.used <- attr(x, "bswf.used", exact = TRUE)
    if (is.null(bswf.used) || length(bswf.used) < 1) {
      # need to handle objects created with old versions
      bswf.used <- "none"
    }
    return(bswf.used[[1]])
  } else {
    return(NA_character_)
  }
}

# is_effective.source_spct defined in file "waveband.class.r" to avoid the need
# of using collate to get the documentation in the correct order.

# Tfr.type attribute ------------------------------------------------------

#' Set the "Tfr.type" attribute
#'
#' Function to set by reference the "Tfr.type" attribute of an existing
#' filter_spct or object_spct object
#'
#' @param x a filter_spct or an object_spct object
#' @param Tfr.type a character string, either "total" or "internal"
#'
#' @return x
#' @note This function alters x itself by reference and in addition
#'   returns x invisibly. If x is not a filter_spct or an object_spct object, x is not modified
#'   The behaviour of this function is 'unusual' in that the default for
#'   parameter \code{Tfr.type} is used only if \code{x} does not already have
#'   this attribute set.
#'
#' @export
#' @family Tfr attribute functions
#' @examples
#' my.spct <- polyester.spct
#' getTfrType(my.spct)
#' setTfrType(my.spct, "internal")
#' getTfrType(my.spct)
#'
setTfrType <- function(x, Tfr.type=c("total", "internal")) {
  name <- substitute(x)
  if (length(Tfr.type) > 1) {
    if (getTfrType(x) != "unknown") {
      Tfr.type <- getTfrType(x)
    } else {
      Tfr.type <- Tfr.type[[1]]
    }
  }
  if (is.filter_spct(x) || is.object_spct(x) ||
      is.summary_filter_spct(x) || is.summary_object_spct(x)) {
    if  (!(Tfr.type %in% c("total", "internal", "unknown"))) {
      warning("Invalid 'Tfr.type' argument, only 'total' and 'internal' supported.")
      return(x)
    }
    attr(x, "Tfr.type") <- Tfr.type
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' Get the "Tfr.type" attribute
#'
#' Function to read the "Tfr.type" attribute of an existing filter_spct or
#' object_spct object.
#'
#' @param x a filter_spct or object_spct object
#'
#' @return character string
#'
#' @note If x is not a \code{filter_spct} or an \code{object_spct} object,
#'   \code{NA} is returned.
#'
#' @export
#' @family Tfr attribute functions
#' @examples
#' getTfrType(polyester.spct)
#'
getTfrType <- function(x) {
  if (is.filter_spct(x) || is.object_spct(x) ||
      is.summary_filter_spct(x) || is.summary_object_spct(x)) {
    Tfr.type <- attr(x, "Tfr.type", exact = TRUE)
    if (is.null(Tfr.type) || is.na(Tfr.type)) {
      # need to handle objects created with old versions
      Tfr.type <- "unknown"
    }
    return(Tfr.type[[1]])
  } else {
    return(NA_character_)
  }
}

# Rfr.type attribute ------------------------------------------------------

#' Set the "Rfr.type" attribute
#'
#' Function to set by reference the "Rfr.type" attribute  of an existing
#' reflector_spct or object_spct object.
#'
#' @param x a reflector_spct or an object_spct object
#' @param Rfr.type a character string, either "total" or "specular"
#'
#' @return x
#' @note This function alters x itself by reference and in addition
#'   returns x invisibly. If x is not a reflector_spct or object_spct object, x is not modified.
#'   The behaviour of this function is 'unusual' in that the default for
#'   parameter Rfr.type is used only if \code{x} does not already have this
#'   attribute set.
#'
#' @export
#' @family Rfr attribute functions
#' @examples
#' my.spct <- reflector_spct(w.length = 400:409, Rfr = 0.1)
#' getRfrType(my.spct)
#' setRfrType(my.spct, "specular")
#' getRfrType(my.spct)
#'
setRfrType <- function(x, Rfr.type=c("total", "specular")) {
  name <- substitute(x)
  if (length(Rfr.type) > 1) {
    if (getRfrType(x) != "unknown") {
      Rfr.type <- getRfrType(x)
    } else {
      Rfr.type <- Rfr.type[[1]]
    }
  }
  if (is.reflector_spct(x) || is.object_spct(x) ||
      is.summary_reflector_spct(x) || is.summary_object_spct(x)) {
    if  (!(Rfr.type %in% c("total", "specular", "unknown"))) {
      warning("Invalid 'Rfr.type' argument, only 'total' and 'internal' supported.")
      return(x)
    }
    attr(x, "Rfr.type") <- Rfr.type
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' Get the "Rfr.type" attribute
#'
#' Function to read the "Rfr.type" attribute of an existing reflector_spct
#' object or object_spct object.
#'
#' @param x a source_spct object
#'
#' @return character string
#'
#' @note if x is not a \code{filter_spct} object, \code{NA} is returned
#'
#' @export
#' @family Rfr attribute functions
#'
getRfrType <- function(x) {
  if (is.reflector_spct(x) || is.object_spct(x) ||
      is.summary_reflector_spct(x) || is.summary_object_spct(x)) {
    Rfr.type <- attr(x, "Rfr.type", exact = TRUE)
    if (is.null(Rfr.type) || is.na(Rfr.type)) {
      # need to handle objects created with old versions
      Rfr.type <- "unknown"
    }
    return(Rfr.type[[1]])
  } else {
    return(NA_character_)
  }
}

# Tfr.type attribute ------------------------------------------------------

#' Set the "K.type" attribute
#'
#' Function to set by reference the "K.type" attribute of an existing
#' solute_spct object
#'
#' @param x a solute_spct or a summary_solute_spct object.
#' @param K.type character A string, either "attenuation", "absorption" or
#'   "scattering".
#'
#' @return x
#' @note This function alters x itself by reference and in addition
#'   returns x invisibly. If x is not a solute_spct object, x is not modified
#'   The behaviour of this function is 'unusual' in that the default for
#'   parameter \code{K.type} is used only if \code{x} does not already have
#'   this attribute set.
#'
#' @export
#' @family K attribute functions
#' @examples
#' print("missing example")
#'
setKType <- function(x,
                     K.type = c("attenuation", "absorption", "scattering")) {
  name <- substitute(x)
  if (length(K.type) > 1) {
    if (getKType(x) != "unknown") {
      K.type <- getKType(x)
    } else {
      K.type <- K.type[[1]]
    }
  }
  if (is.solute_spct(x) || is.summary_solute_spct(x)) {
    if  (!(K.type %in% c("attenuation", "absorption", "scattering"))) {
      warning("Invalid 'K.type' argument, only 'attenuation', 'absorption' and 'scattering' supported.")
      return(x)
    }
    attr(x, "K.type") <- K.type
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' Get the "K.type" attribute
#'
#' Function to read the "K.type" attribute of an existing solute_spct object.
#'
#' @param x a solute_spct object
#'
#' @return character string
#'
#' @note If x is not a \code{solute_spct} or a \code{summary_solute_spct} object,
#'   \code{NA} is returned.
#'
#' @export
#' @family K attribute functions
#' @examples
#' print("missing example")
#'
getKType <- function(x) {
  if (is.solute_spct(x) || is.summary_solute_spct(x)) {
    K.type <- attr(x, "K.type", exact = TRUE)
    if (is.null(K.type) || is.na(K.type)) {
      # need to handle corrupted objects
      K.type <- "unknown"
    }
    return(K.type[[1]])
  } else {
    return(NA_character_)
  }
}

# spct.version ------------------------------------------------------------

#' Get the "spct.version" attribute
#'
#' Function to read the "spct.version" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#'
#' @return integer value
#'
#' @note if x is not a \code{generic_spct} object, \code{NA} is returned,
#'   and if it the attribute is missing, zero is returned with a warning.
#'
#' @export
#'
getSpctVersion <- function(x) {
  if (is.generic_spct(x) || is.old_spct(x)) {
    version <- attr(x, "spct.version", exact = TRUE)
    if (is.null(version)) {
      # need to handle objects created with old versions
      version <- 0L
    }
  } else {
    version <- NA_integer_
  }
  version
}

#' Check that the "spct.version" attribute is set
#'
#' Function to check the "spct.version" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#'
#' @return numeric value
#'
#' @note if x is not a \code{generic_spct} object, \code{NA} is returned,
#'   and if it the attribute is missing, zero is returned with a warning.
#'
#' @keywords internal
#'
checkSpctVersion <- function(x) {
  version <- getSpctVersion(x)
  stopifnot(!is.na(version))
  if (version < 1L) {
    warning("The object '", as.character(substitute(x)),
            "' was created in a version (< 0.7.0) or has become corrupted")
  }
}


# multiple wl -------------------------------------------------------------


#' Find repeated w.length values
#'
#' @param x a generic_spct object
#' @param same.wls logical If TRUE all spectra spected to share same w.length
#'   values.
#'
#' @return integer Number of spectra, guessed from the number of copies of each
#'   individual w.length value.
#'
findMultipleWl <- function(x, same.wls = TRUE) {
  stopifnot(is.generic_spct(x))
  if (nrow(x) == 0L) {
    return(0L)
  }
  runs <- rle(sort(x[["w.length"]]))
  if (same.wls) {
    num.copies <- unique(runs[["lengths"]])
    stopifnot(length(num.copies) %in% c(0L, 1L))
  } else {
    num.copies <- max(runs[["lengths"]])
  }
  num.copies
}

#' Set the "multiple.wl" attribute
#'
#' Function to set by reference the "multiple.wl" attribute  of an existing
#' generic_spct or an object of a class derived from generic_spct.
#'
#' @param x a generic_spct object
#' @param multiple.wl numeric >= 1 If \code{multiple.wl} is \code{NULL}, the
#'   default, the attribute is not modified if it is already present and valid,
#'   and set to 1 otherwise.
#'
#' @return x
#' @note This function alters x itself by reference and in addition
#'   returns x invisibly. If x is not a generic_spct or an object of a class derived from
#'   generic_spct, x is not modified. If \code{multiple.wl}
#'
#' @export
#' @family multiple.wl attribute functions
#'
setMultipleWl <- function(x, multiple.wl = NULL) {
  stopifnot(is.generic_spct(x) || is.summary_generic_spct(x))
  name <- substitute(x)
  if (is.null(multiple.wl)) {
    multiple.wl <- findMultipleWl(x)
  } else {
    multiple.wl <- trunc(multiple.wl)
    stopifnot(multiple.wl >= 0) # 0L only for empty spectral objects
  }
  attr(x, "multiple.wl") <- multiple.wl
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' Get the "multiple.wl" attribute
#'
#' Function to read the "multiple.wl" attribute of an existing generic_spct.
#'
#' @param x a generic_spct object
#'
#' @return integer
#'
#' @note If x is not a \code{generic_spct} or an object of a derived class
#'   \code{NA} is returned.
#'
#' @export
#' @family multiple.wl attribute functions
#' @examples
#' getMultipleWl(sun.spct)
#'
getMultipleWl <- function(x) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    multiple.wl <- attr(x, "multiple.wl", exact = TRUE)
    if (is.null(multiple.wl) || is.na(multiple.wl) || !is.numeric(multiple.wl)) {
      # need to handle objects created with old versions
      multiple.wl <- 1
    }
    return(multiple.wl)
  } else {
    return(NA_integer_)
  }
}


# idfactor -------------------------------------------------------------

#' Set the "idfactor" attribute
#'
#' Function to set by reference the "idfactor" attribute  of an existing
#' generic_spct or an object of a class derived from generic_spct.
#'
#' @param x a generic_spct object
#' @param idfactor character The name of a factor identifying multiple
#'    spectra stored longitudinally.
#'
#' @return x
#'
#' @note This function alters x itself by reference and in addition
#'   returns x invisibly. If x is not a generic_spct or an object of a class derived from
#'   generic_spct, x is not modified.
#'
#' @export
#' @family idfactor attribute functions
#'
setIdFactor <- function(x, idfactor) {
  stopifnot(is.generic_spct(x) || is.summary_generic_spct(x))
  stopifnot(is.null(idfactor) || is.character(idfactor))
  name <- substitute(x)
  if (is.null(idfactor) || exists(idfactor, x, inherits = FALSE)) {
    attr(x, "idfactor") <- idfactor
  } else {
    stop("'idfactor' points to a non-existant variable")
  }
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' Get the "idfactor" attribute
#'
#' Function to read the "idfactor" attribute of an existing generic_spct.
#'
#' @param x a generic_spct object
#'
#' @return character
#'
#' @note If x is not a \code{generic_spct} or an object of a derived class
#'   \code{NA} is returned.
#'
#' @export
#' @family idfactor attribute functions
#' @examples
#' getMultipleWl(sun.spct)
#'
getIdFactor <- function(x) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    idfactor <- attr(x, "idfactor", exact = TRUE)
    if (is.null(idfactor) || is.na(idfactor) || !is.character(idfactor)) {
      # need to handle objects created with old versions
      idfactor <- NA_character_
    }
  } else {
    idfactor <- NA_character_
  }
  idfactor
}

# "filter.properties" attribute ----------------------------------------------

#' Set the "filter.properties" attribute
#'
#' Function to set by reference the "filter.properties" attribute  of an existing
#' filter_spct object.
#'
#' @param x a filter_spct object
#' @param filter.properties,value a list with fields named \code{"Rfr.constant"},
#'   \code{"thickness"} and \code{"attenuation.mode"}.
#' @param pass.null logical If TRUE, the parameters to the next three
#'    parameters will be always ignored, otherwise they will be used to
#'    build an object of class \code{"filter.properties"} when the argument
#'    passed to parameter \code{filter.properties} is \code{NULL}.
#' @param Rfr.constant numeric The value of the reflection factor [/1].
#' @param thickness numeric The thickness of the material [\eqn{m}].
#' @param attenuation.mode character One of \code{"reflection"}, \code{"absorption"},
#'    \code{"absorption.layer"}, \code{"mixed"} or \code{"stack"}.
#'
#' @details Storing filter properties allows inter-conversion between internal
#'   and total transmittance, as well as computation of transmittance for
#'   arbitrary thickness of the material. Whether computations are valid depend
#'   on the homogeneity of the material. The parameter \code{pass.null} makes
#'   it possible to remove the attribute.
#'
#' @return \code{x}
#' @note This function alters \code{x} itself by reference and in addition
#'   returns \code{x} invisibly. If \code{x} is not a \code{filter_spct} object,
#'   \code{x} is not modified.
#'
#'   The values of \code{attenuation.mode} \code{"reflection"} and
#'   \code{"absorption"} should be used when one of these processes is clearly
#'   the main one; \code{"mixed"} is for cases when they both play a role, i.e.,
#'   when a simple correction using a single value of \code{Rfr} across
#'   wavelengths is not possible; \code{"absorption.layer"} is for cases when a
#'   thin absorbing layer is deposited on the surface of a transparent support
#'   or enclosed between two sheets of glass or other transparent material. If
#'   in doubt, set this to \code{NA} to ensure that computation of spectra for a
#'   different thickness remains disabled.
#'
#' @export
#'
#' @family measurement metadata functions
#'
#' @examples
#'
#' my.spct <- polyester.spct
#' filter_properties(my.spct)
#' filter_properties(my.spct) <- NULL
#' filter_properties(my.spct)
#' filter_properties(my.spct, return.null = TRUE)
#' filter_properties(my.spct) <- list(Rfr.constant = 0.01,
#'                                    thickness = 125e-6,
#'                                    attenuation.mode = "absorption")
#' filter_properties(my.spct)
#'
setFilterProperties <- function(x,
                                filter.properties = NULL,
                                pass.null = FALSE,
                                Rfr.constant = NA_real_,
                                thickness = NA_real_,
                                attenuation.mode = NA_character_) {
  name <- substitute(x)
  if (is.filter_spct(x) || is.object_spct(x)) {
    if (!(pass.null && is.null(filter.properties))) {
      if (is.null(filter.properties)) {
        filter.properties <- list(Rfr.constant = Rfr.constant,
                                  thickness = thickness,
                                  attenuation.mode = attenuation.mode)
        class(filter.properties) <-
          c("filter_properties", class(filter.properties))
      } else {
        stopifnot(setequal(names(filter.properties),
                           c("Rfr.constant", "thickness", "attenuation.mode")))
        if (class(filter.properties)[1] != "filter_properties") {
          class(filter.properties) <-
            c("filter_properties", class(filter.properties))
        }
      }
      if (!is.numeric(filter.properties[["Rfr.constant"]])) {
        filter.properties[["Rfr.constant"]] <-
          as.numeric(filter.properties[["Rfr.constant"]])
      }
      if (!is.na(filter.properties[["Rfr.constant"]]) &&
          (filter.properties[["Rfr.constant"]] < 0 ||
          filter.properties[["Rfr.constant"]] > 0.2)) {
        warning("Off-range value '",
                filter.properties[["Rfr.constant"]],
                "' for \"Rfr.constant\" set to NA")
        filter.properties[["Rfr.constant"]] <- NA_real_
      }
      if (!is.numeric(filter.properties[["thickness"]])) {
        filter.properties[["thickness"]] <-
          as.numeric(filter.properties[["thickness"]])
      }
      if (any(!is.na(filter.properties[["thickness"]]) &
          filter.properties[["thickness"]] <= 0)) {
        warning("'thickness' (m) <= 0 set to NA")
        filter.properties[["thickness"]][!is.na(filter.properties[["thickness"]]) &
                                           filter.properties[["thickness"]] <= 0] <- NA_real_
      }
      # one could have a list with the properties of the stacked filters as an additional field
      # but would require surgery of the code in several other places and careful thought
      if (length(filter.properties[["thickness"]]) > 1L &&
          filter.properties[["attenuation.mode"]] != "stack") {
        filter.properties[["attenuation.mode"]] <- "stack"
      }
      if (!is.character(filter.properties[["attenuation.mode"]])) {
        # handle NA which is logical or numeric
        filter.properties[["attenuation.mode"]] <-
          as.character(filter.properties[["attenuation.mode"]])
       }
      if (!is.na(filter.properties[["attenuation.mode"]]) &&
                  !filter.properties[["attenuation.mode"]] %in%
             c("reflection", "absorption", "absorption.layer", "mixed", "stack")) {
        warning("Bad value(s) '",
                filter.properties[["attenuation.mode"]],
                "' for \"attenuation.mode\" set to NA")
        filter.properties[["attenuation.mode"]][!filter.properties[["attenuation.mode"]] %in%
                                                  c("reflection", "absorption", "absorption.layer", "mixed", "stack")] <- NA_character_
      }
      if (any(!is.na(filter.properties[["Rfr.constant"]])) &&
          filter.properties[["attenuation.mode"]] == "stack") {
        warning("Setting 'Rfr.constant' to 'NA' for filter stack")
        filter.properties[["Rfr.constant"]] <- NA_real_
      }
    }
    attr(x, "filter.properties") <- filter.properties
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  } else {
    warning("'setFilterProperties()' not applicable to objects of class ",
            class(x)[1], ", skipping.")
  }
  invisible(x)
}

#' @rdname setFilterProperties
#'
#' @export
#'
`filter_properties<-` <- function(x,
                                  value = NULL) {
  setFilterProperties(x = x,
                      filter.properties = value,
                      pass.null = TRUE)
}

#' Get the "filter.properties" attribute
#'
#' Function to read the "filter.properties" attribute of an existing filter_spct
#' or a filter_mspct.
#'
#' @param x a filter_spct object
#' @param return.null logical If true, \code{NULL} is returned if the attribute
#'   is not set, otherwise the expected list is returned with all fields set to
#'   \code{NA}.
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return a list with fields named \code{"Rfr.constant"} [\eqn{/1}],
#'   \code{"thickness"} [\eqn{m}] and \code{"attenuation.mode"}. If the
#'   attribute is not set, and \code{return.null} is FALSE, a list with fields
#'   set to \code{NA} is returned, otherwise, \code{NULL}.
#'
#' @export
#' @family measurement metadata functions
#'
#' @examples
#' filter_properties(polyester.spct)
#'
getFilterProperties <- function(x, return.null, ...) UseMethod("getFilterProperties")

#' @rdname getFilterProperties
#'
#' @export
#'
filter_properties <- getFilterProperties

#' @describeIn getFilterProperties default
#' @export
getFilterProperties.default <- function(x,
                                        return.null = FALSE,
                                        ...) {
  if (!is.any_spct(x) && !is.any_summary_spct(x)) {
    warning("Methods 'getFilterProperties()' not implemented for class: ",
            class(x)[1])
  }
  if (return.null) {
    NULL
  } else {
    # we return an NA
    filter.properties <- list(Rfr.constant = NA_real_,
                              thickness = NA_real_,
                              attenuation.mode = NA)
    class(filter.properties) <-
      c("filter_properties", class(filter.properties))
    filter.properties
  }
}

#' @describeIn getFilterProperties generic_spct
#' @export
getFilterProperties.filter_spct <- function(x,
                                            return.null = FALSE,
                                            ...) {
  filter.properties <- attr(x, "filter.properties", exact = TRUE)
  if (is.null(filter.properties)) {
    if (!return.null) {
      # need to handle objects created with old versions
      filter.properties <- list(Rfr.constant = NA_real_,
                                thickness = NA_real_,
                                attenuation.mode = NA)
      class(filter.properties) <-
        c("filter_properties", class(filter.properties))
    }
  } else {
    stopifnot(setequal(names(filter.properties),
                       c("Rfr.constant", "thickness", "attenuation.mode")))
  }
  filter.properties
}

#' @describeIn getFilterProperties summary_generic_spct
#'
#' @export
#'
getFilterProperties.summary_filter_spct <- getFilterProperties.filter_spct

#' @describeIn getFilterProperties filter_mspct
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @note The method for collections of spectra returns the a tibble with a
#'   column of lists.
#' @export
#'
getFilterProperties.generic_mspct <- function(x,
                                              return.null = FALSE,
                                              ...,
                                              idx = "spct.idx") {
  l <- mslply(mspct = x, .fun = getFilterProperties, ...)
  comment(l) <- NULL
  z <- list(filter.properties = l)
  z[[idx]] <- factor(names(l), levels = names(l))
  tibble::as_tibble(z[c(2, 1)])
}

# Modify filter properties -----------------------------------------------

#' Convert the "thickness" attribute of an existing filter_spct object.
#'
#' Function to set the "thickness" attribute and simultaneously converting the
#' spectral data to correspond to the new thickness.
#'
#' @details For spectral transmittance at a different thickness to be exactly
#'   computed, it needs to be based on internal transmittance. This function
#'   will apply \code{converTfrType()} to \code{x} if needed, but to succeed
#'   metadata should be available. Please, see \code{\link{convertTfrType}}.
#'
#' @param x a filter_spct, object_spct, filter_mspct or object_mspct object.
#' @param thickness numeric [\eqn{m}].
#'
#' @return \code{x} possibly with the \code{"thickness"} field of the
#'   \code{"filter.properties"} attribute modified and \code{Tfr} or \code{A}
#'   computed for the requested thickness.
#'
#' @note if \code{x} is not a \code{filter_spct}, \code{object_spct},
#'   \code{filter_mspct} or \code{object_mspct} object or a collection of such
#'   objects, \code{x} is returned unchanged. If \code{x} does not have the
#'   \code{"filter.properties"} attribute set or has it with missing member
#'   data, \code{x} is returned with \code{Tfr} set to \code{NA} values.
#'
#' @export
#' @family time attribute functions
#' @examples
#'
#' my.spct <- polyester.spct
#' filter_properties(my.spct)
#' convertThickness(my.spct, thickness = 250e-6)
#'
convertThickness <- function(x, thickness = NULL) {
  if (is.filter_mspct(x) || is.object_mspct(x)) {
    return(msmsply(mspct = x,
                   .fun =  convertThickness,
                   thickness = thickness))
  }
  if (!(is.filter_spct(x) || is.object_spct(x))) {
    warning("'convertThickness()' not applicable to class '", class(x)[1], "'. Skipping!")
    return(invisible(x))
  }
  if (is.null(thickness)) {
    # nothing to do
    return(invisible(x))
  }

  properties <- filter_properties(x)
  if (properties[["attenuation.mode"]] == "mixed") {
    warning("Conversion not possible for non-absorbent materials.")
    return(x * NA_real_)
  } else if (properties[["attenuation.mode"]] == "absorption.layer") {
    warning("Conversion is undefined when absorbing material is layered.")
    return(x * NA_real_)
  } else if (properties[["attenuation.mode"]] == "stack") {
    warning("Conversion is undefined for stacks of filters.")
    return(x * NA_real_)
  } else if(properties[["attenuation.mode"]] == "reflection") {
    warning("Transmittance remains unchanged for purely reflective materials.")
    properties[["thickness"]] <- thickness
    setFilterProperties(x, properties)
    return(x)
  } else if (properties[["attenuation.mode"]] == "absorption") {
    columns <- intersect(colnames(x), c("Tfr", "Afr", "A") )
    if (length(columns) == 0) {
      warning("No column to convert to new thickness.")
      return(invisible(x))
    }
    if ("Tfr" %in% columns || "A" %in% columns) {
      if (!"Tfr" %in% columns) {
        .fun <- T2A
      } else {
        .fun <- NULL
      }
      # "A" column converted or deleted as needed
      z <- A2T(x, action = "replace")
    } else if ("Afr" %in% columns) {
      .fun <- T2Afr
      z <- Afr2T(x, action = "replace")
    } else {
      stop("conversion failed")
    }

    current.Tfr.type <- getTfrType(x)
    if (current.Tfr.type == "total") {
      z <- convertTfrType(z, "internal")
    }
    # convert Tfr, formula is valid only for internal transmittance
    z <- using_Tfr(z^(thickness / properties[["thickness"]]))
    properties[["thickness"]] <- thickness
    setFilterProperties(z, properties)
    if (current.Tfr.type == "total") {
      z <- convertTfrType(z, "total")
    }
    if (!is.null(.fun)) {
      z <- .fun(z, action = "replace")
    }
    z
  }
}

#' Convert the "Tfr.type" attribute
#'
#' Function to set the \code{"Tfr.type"} attribute and simultaneously converting
#' the spectral data to correspond to the new type.
#'
#' @details Internal transmittance uses as reference the light entering the
#'   object while total transmittance takes the incident light as reference. The
#'   conversion is possible only if reflectance is known. Either as spectral
#'   data in an \code{object_spct} object, a \code{filter_spct} object that is
#'   "under-the-hood" an \code{object_spct}, or if a fixed reflectance factor
#'   applicable to all wavelengths is stored in the \code{filter.properties}
#'   attribute of the \code{filter_spct} object.
#'
#' @param x a filter_spct, object_spct, filter_mspct or object_mspct object.
#' @param Tfr.type character One of \code{"internal"} or \code{"total"}.
#'
#' @return \code{x} possibly with the \code{"thickness"} field of the
#'   \code{"filter.properties"} attribute modified
#'
#' @note if \code{x} is not a \code{filter_spct} object, \code{x} is returned
#'   unchanged. If \code{x} does not have the \code{"filter.properties"}
#'   attribute set if it is missing data, \code{x} is returned with
#'   \code{Tfr} set to \code{NA} values.
#'
#' @export
#' @family time attribute functions
#' @examples
#'
#' my.spct <- polyester.spct
#' filter_properties(my.spct) <- list(Rfr.constant = 0.07,
#'                                    thickness = 125e-6,
#'                                    attenuation.mode = "absorption")
#' convertTfrType(my.spct, Tfr.type = "internal")
#'
convertTfrType <- function(x, Tfr.type = NULL) {
  if (is.filter_mspct(x) || is.object_mspct(x)) {
    return(msmsply(mspct = x,
                   .fun =  convertTfrType,
                   Tfr.type = Tfr.type))
  }
  if (!(is.filter_spct(x) || is.object_spct(x))) {
    warning("'convertTfrType()' mot applicable to class '", class(x)[1L], "'. Skipping!")
    return(invisible(x))
  }

  if (is.null(Tfr.type) || Tfr.type[1] == getTfrType(x)[1]) {
    # nothing to do
    return(invisible(x))
  }

  columns <- intersect(colnames(x), c("Tfr", "Afr", "A", "Rfr") )
  if (length(setdiff(columns, "Rfr")) == 0L) {
    warning("No column to convert to new Tfr.type")
    return(invisible(x))
  }

  # we keep columns "A" or "Afr" as their values do not depend on "Tfr.type"
  if (is.filter_spct(x)) {
    if ("Rfr" %in% columns) {
      z <- as.object_spct(x)
    } else if ("Tfr" %in% columns) {
      z <- x
    } else if ("A" %in% columns) {
      # "A" column converted as needed
      z <- A2T(x, action = "add")
    } else if ("Afr" %in% columns) {
      z <- Afr2T(x, action = "add")
    } else {
      stop("conversion of input failed")
    }
  } else { # user passed an object_spct
    z <- x
  }

  current.Tfr.type <- getTfrType(x)
  if (is.filter_spct(z)) {
    # no spectral Rfr available, we use a factor
    properties <- filter_properties(x)
    if (is.na(properties[["attenuation.mode"]]) ||
        properties[["attenuation.mode"]] != "absorption") {
      properties[["Rfr.constant"]] <- NA_real_
      if (!is.na(properties[["attenuation.mode"]])) {
        warning("Setting internal Tfr to NA as filter is not absorptive")
      }
    }
    if (is.na(current.Tfr.type)) {
      warning("Current Tfr type is not set, returning NAs.")
    }
    if (current.Tfr.type == "internal" && Tfr.type == "total") {
      z[["Tfr"]] <- z[["Tfr"]] * (1 - properties[["Rfr.constant"]])
    } else if (current.Tfr.type == "total" && Tfr.type == "internal") {
      z[["Tfr"]] <- z[["Tfr"]] / (1 - properties[["Rfr.constant"]])
    }
    setTfrType(z, Tfr.type)
  } else if (is.object_spct(z)) {
    if (current.Tfr.type == "internal" && Tfr.type == "total") {
      z[["Tfr"]] <- z[["Tfr"]] * (1 - z[["Rfr"]])
    } else if (current.Tfr.type == "total" && Tfr.type == "internal") {
      z[["Tfr"]] <- z[["Tfr"]] / (1 - z[["Rfr"]])
    }
    setTfrType(z, Tfr.type)
    if (is.filter_spct(x)) {
      z <- as.filter_spct(z)
    }
  }
  z
}

# "solute.properties" attribute ----------------------------------------------

#' Set the "solute.properties" attribute
#'
#' Function to set by reference the \code{"solute.properties"} attribute of an
#' existing \code{solute_spct} object.
#'
#' @param x solute_spct A spectrum of coefficients of attenuation.
#' @param solute.properties,value a list with fields named \code{"mass"}, \code{"formula"},
#'   \code{"structure"}, \code{"name"} and \code{"ID"}.
#' @param pass.null logical If \code{TRUE}, the parameters to the next three
#'    parameters will be always ignored, otherwise they will be used to
#'    build an object of class \code{"solute.properties"} when the argument to
#'    \code{solute.properties} is \code{NULL}.
#' @param mass numeric The mass in Dalton [\eqn{Da = g\,mol^{-1}}{Da = g/mol}].
#' @param formula character The molecular formula.
#' @param structure raster A bitmap of the structure.
#' @param name,solvent.name character The name of the substance and the name of the solvent. A named character
#'     vector, with member names such as "IUPAC" for the authority.
#' @param ID,solvent.ID character The names of the substance and of the solvent. A named character
#'     vector, with member names such as "ChemSpider" or "PubChen" for the
#'     authority.
#'
#' @details Storing solute properties allows inter-conversion between bases of
#'   expression, and ensures the unambiguous identification of the substances to
#'   which the spectral data refer. These properties make it possible to compute
#'   \code{filter_spct} objects for solutions of the solute, i.e., absorption
#'   spectra of liquid filters. The parameter \code{pass.null} makes it possible
#'   to remove the attribute. The solvent used for the determination of the
#'   attenuation coefficient is important metadata as the solvent can alter
#'   the spectral ansorption properties of the solute.
#'
#' @return \code{x}
#'
#' @note This function alters \code{x} itself by reference and in addition
#'   returns \code{x} invisibly. If \code{x} is not a filter_spct object,
#'   \code{x} is not modified.
#'
#' @export
#' @family measurement metadata functions
#'
#' @examples
#'
#' solute.properties <-
#'   list(formula = c(text = "H2O", html = "H<sub>2</sub>", TeX = "$H_2O$"),
#'        name = c("water", IUPAC = "oxidane"),
#'        structure = grDevices::as.raster(matrix()),
#'        mass = 18.015, # Da
#'        ID = c(ChemSpider = "917", CID = "962"),
#'        solvent.name = NA_character_,
#'        solvent.ID = NA_character_)
#' my.spct <- solute_spct()
#' solute_properties(my.spct) <- solute.properties
#' solute_properties(my.spct)
#' solute_properties(my.spct) <- NULL
#' solute_properties(my.spct)
#' solute_properties(my.spct, return.null = TRUE)
#' solute_properties(my.spct)
#'
setSoluteProperties <- function(x,
                                solute.properties = NULL,
                                pass.null = FALSE,
                                mass = NA_real_,
                                formula = NULL,
                                structure = grDevices::as.raster(matrix()),
                                name = NA_character_,
                                ID = NA_character_,
                                solvent.name = NA_character_,
                                solvent.ID = NA_character_) {
  obj.name <- substitute(x)
  if (is.solute_spct(x)) {
    if (!(pass.null && is.null(solute.properties))) {
      if (is.null(solute.properties)) {
        solute.properties <- list(mass = mass,
                                  formula = formula,
                                  structure = structure,
                                  name = name,
                                  ID = ID,
                                  solvent.name = solvent.name,
                                  solvent.ID = solvent.ID)
        class(solute.properties) <-
          c("solute_properties", class(solute.properties))
      } else {
        stopifnot(setequal(names(solute.properties),
                           c("mass", "formula", "structure", "name", "ID", "solvent.name", "solvent.ID")))
        if (class(solute.properties)[1] != "solute_properties") {
          class(solute.properties) <-
            c("solute_properties", class(solute.properties))
        }
      }
      if (!is.numeric(solute.properties[["mass"]])) {
        solute.properties[["mass"]] <-
          as.numeric(solute.properties[["mass"]])
      }
      if (!is.na(solute.properties[["mass"]]) &&
          solute.properties[["mass"]] <= 1) {
        warning("Found 'mass' (Da = g/mol) <= 1 and set it to NA")
        solute.properties[["mass"]] <- NA_real_
      }
      if (!is.character(solute.properties[["formula"]])) {
        solute.properties[["formula"]] <-
          as.character(solute.properties[["formula"]])
      }
      if (!grDevices::is.raster(solute.properties[["structure"]])) {
        solute.properties[["structure"]] <-
          grDevices::as.raster(solute.properties[["structure"]])
      }
      if (!is.character(solute.properties[["name"]])) {
        solute.properties[["name"]] <-
          as.character(solute.properties[["name"]])
      }
      if (!is.character(solute.properties[["ID"]])) {
        solute.properties[["ID"]] <-
          as.character(solute.properties[["ID"]])
      }
    }
    attr(x, "solute.properties") <- solute.properties
    if (is.name(obj.name)) {
      obj.name <- as.character(obj.name)
      assign(obj.name, x, parent.frame(), inherits = TRUE)
    }
  } else {
    warning("'setSoluteProperties()' not applicable to objects of class ",
            class(x)[1], ", skipping.")
  }
  invisible(x)
}

#' @rdname setSoluteProperties
#'
#' @export
#'
`solute_properties<-` <- function(x,
                                  value = NULL) {
  setSoluteProperties(x = x,
                      solute.properties = value,
                      pass.null = TRUE)
}

#' Get the "solute.properties" attribute
#'
#' Function to read the \code{"solute.properties"} attribute of an existing
#' \code{solute_spct} or a \code{solute_mspct} objects.
#'
#' @param x solute_spct A spectrum of coefficients of attenuation.
#' @param return.null logical If true, \code{NULL} is returned if the attribute
#'   is not set, otherwise the expected list is returned with all fields set to
#'   \code{NA}.
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return a \code{list} with fields named \code{"mass"}, \code{"formula"},
#'   \code{"structure"}, \code{"name"} and \code{"ID"}. If the attribute is not
#'   set, and \code{return.null} is \code{FALSE}, a list with fields set to
#'   \code{NA} is returned, otherwise, \code{NULL}.
#'
#' @export
#' @family measurement metadata functions
#'
#' @examples
#' solute_properties(water.spct)
#'
getSoluteProperties <- function(x, return.null, ...) UseMethod("getSoluteProperties")

#' @rdname getSoluteProperties
#'
#' @export
#'
solute_properties <- getSoluteProperties

#' @describeIn getSoluteProperties default
#' @export
getSoluteProperties.default <- function(x,
                                        return.null = FALSE,
                                        ...) {
  if (!is.any_spct(x) && !is.any_summary_spct(x)) {
    warning("Methods 'getSoluteProperties()' not implemented for class: ",
            class(x)[1])
  }
  if (return.null) {
    NULL
  } else {
    # we return an NA
    solute.properties <- list(mass = NA_real_,
                              formula = NA_character_,
                              structure = grDevices::as.raster(matrix()),
                              name = NA_character_,
                              ID = NA_character_,
                              solvent.name =  NA_character_,
                              solvent.ID = NA_character_)
    class(solute.properties) <-
      c("solute_properties", class(solute.properties))
    solute.properties
  }
}

#' @describeIn getSoluteProperties solute_spct
#' @export
getSoluteProperties.solute_spct <- function(x,
                                            return.null = FALSE,
                                            ...) {
  solute.properties <- attr(x, "solute.properties", exact = TRUE)
  if (is.null(solute.properties)) {
    if (!return.null) {
      # need to handle objects created with old versions
      solute.properties <- list(mass = NA_real_,
                                formula = NA_character_,
                                structure = grDevices::as.raster(matrix()),
                                name = NA_character_,
                                ID = NA_character_,
                                solvent.name =  NA_character_,
                                solvent.ID = NA_character_)
      class(solute.properties) <-
        c("solute_properties", class(solute.properties))
    }
  } else {
    stopifnot(setequal(names(solute.properties),
                       c("mass", "formula", "structure", "name", "ID", "solvent.name", "solvent.ID")))
  }
  solute.properties
}

#' @describeIn getSoluteProperties summary_solute_spct
#'
#' @export
#'
getSoluteProperties.summary_solute_spct <- getSoluteProperties.solute_spct

#' @describeIn getSoluteProperties solute_mspct
#'
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#'
#' @note The method for collections of spectra returns the a tibble with a
#'   column of lists.
#'
#' @export
#'
getSoluteProperties.solute_mspct <- function(x,
                                             return.null = FALSE,
                                             ...,
                                             idx = "spct.idx") {
  l <- mslply(mspct = x, .fun = getSoluteProperties, ...)
  comment(l) <- NULL
  z <- list(solute.properties = l)
  z[[idx]] <- factor(names(l), levels = names(l))
  tibble::as_tibble(z[c(2, 1)])
}
