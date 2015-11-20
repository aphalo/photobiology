#' Clean (=replace) off-range values in a spectrum
#'
#' These functions implement the equivalent of replace() but for spectral
#' objects instead of vectors.
#'
#' @param x an R object
#' @param range numeric vector of wavelengths
#' @param range.s.data numeric vector of length two giving the allowable
#'     range for the spectral data.
#' @param fill numeric vector of length 1 or 2, giving the replacement
#'     values to use at each extreme of the range.
#' @param ... currently ignored
#'
#' @keywords manip misc
#' @export
#'
clean <- function(x, range, range.s.data, fill, ...) UseMethod("clean")

#' @describeIn clean Default for generic function
#'
#' @export
#'
clean.default <- function(x, range, range.s.data, fill, ...) {
  warning("'clean' is not defined for objects of class ", class(x)[1])
  return(x)
}

#' @describeIn clean Replace off-range values in a source spectrum
#' @param unit.out character string with allowed values "energy", and "photon",
#'   or its alias "quantum"
#'
#' @export
#'
clean.source_spct <-
  function(x,
           range = x,
           range.s.data = c(0,NA),
           fill = range.s.data,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...) {
    stopifnot(length(range) >= 2L &&
                length(range.s.data) == 2L &&
                length(fill) <= 2L)
    # wavelength range
    if (is.any_spct(range) || is.numeric(range) && length(range) > 2L) {
      range <- range(range, na.rm = TRUE)
    } else {
      if (is.na(range[1])) {
        range[1] <- min(x)
      }
      if (is.na(range[2])) {
        range[2] <- max(x)
      }
    }
    selector <- x[["w.length"]] >= range[1] & x[["w.length"]] <= range[2]

    if (unit.out == "quantum") {
      unit.out <- "photon"
    }
    if (unit.out == "photon") {
      x <- e2q(x, action = "replace")
      if (is.na(range.s.data[1])) {
        range.s.data[1] <- min(x[["s.q.irrad"]], na.rm = TRUE)
      }
      if (is.na(range.s.data[2])) {
        range.s.data[2] <- max(x[["s.q.irrad"]], na.rm = TRUE)
      }
      x[selector, "s.q.irrad"] <- with(x[selector, ],
                                       ifelse(s.q.irrad < range.s.data[1],
                                              fill[1],
                                              ifelse(s.q.irrad  > range.s.data[2],
                                                    fill[2],
                                                    s.q.irrad)))
    } else if (unit.out == "energy") {
      x <- q2e(x, action = "replace")
      if (is.na(range.s.data[1])) {
        range.s.data[1] <- min(x[["s.e.irrad"]], na.rm = TRUE)
      }
      if (is.na(range.s.data[2])) {
        range.s.data[2] <- max(x[["s.e.irrad"]], na.rm = TRUE)
      }
      x[selector, "s.e.irrad"] <- with(x[selector, ],
                                       ifelse(s.e.irrad < range.s.data[1],
                                              fill[1],
                                              ifelse(s.e.irrad  > range.s.data[2],
                                              fill[2],
                                              s.e.irrad)))
    } else {
      stop("unit.out: '", unit.out, "' not supported")
    }
   x
  }

#' @describeIn clean Replace off-range values in a filter spectrum
#' @param qty.out character string with allowed values "energy", and "photon",
#'   or its alias "quantum"
#'
#' @export
#'
clean.filter_spct <-
  function(x,
           range = x,
           range.s.data = NULL,
           fill = range.s.data,
           qty.out = getOption("photobiology.filter.qty", default = "transmittance"),
           ...) {
    if (is.null(range.s.data)) {
      if (qty.out == "transmittance") {
        range.s.data <- c(0,1)
      } else {
        range.s.data <- c(0,NA)
      }
      fill <- range.s.data
    }
    stopifnot(length(range) >= 2L &&
                length(range.s.data) == 2L &&
                length(fill) <= 2L)
    # wavelength range
    if (is.any_spct(range) || is.numeric(range) && length(range) > 2L) {
      range <- range(range, na.rm = TRUE)
    } else {
      if (is.na(range[1])) {
        range[1] <- min(x)
      }
      if (is.na(range[2])) {
        range[2] <- max(x)
      }
    }
    selector <- x[["w.length"]] >= range[1] & x[["w.length"]] <= range[2]

    if (qty.out == "transmittance") {
      x <- A2T(x, action = "replace")
      if (is.na(range.s.data[1])) {
        range.s.data[1] <- min(x[["Tfr"]], na.rm = TRUE)
      }
      if (is.na(range.s.data[2])) {
        range.s.data[2] <- max(x[["Tfr"]], na.rm = TRUE)
      }
      x[selector, "Tfr"] <- with(x[selector, ],
                                 ifelse(Tfr < range.s.data[1],
                                        fill[1],
                                        ifelse(Tfr  > range.s.data[2],
                                               fill[2],
                                               Tfr)))
    } else if (qty.out == "absorbance") {
      x <- T2A(x, action = "replace")
      if (is.na(range.s.data[1])) {
        range.s.data[1] <- min(x[["A"]], na.rm = TRUE)
      }
      if (is.na(range.s.data[2])) {
        range.s.data[2] <- max(x[["A"]], na.rm = TRUE)
      }
      x[selector, "A"] <- with(x[selector, ],
                                 ifelse(A < range.s.data[1],
                                        fill[1],
                                        ifelse(A  > range.s.data[2],
                                               fill[2],
                                               A)))
    } else {
      stop("qty.out: '", qty.out, "' not supported")
    }
    x
  }

#' @describeIn clean Replace off-range values in a reflector spectrum
#'
#' @export
#'
clean.reflector_spct <-
  function(x,
           range = x,
           range.s.data = c(0, 1),
           fill = range.s.data,
           ...) {
    stopifnot(length(range) >= 2L &&
                length(range.s.data) == 2L &&
                length(fill) <= 2L)
    # wavelength range
    if (is.any_spct(range) || is.numeric(range) && length(range) > 2L) {
      range <- range(range, na.rm = TRUE)
    } else {
      if (is.na(range[1])) {
        range[1] <- min(x)
      }
      if (is.na(range[2])) {
        range[2] <- max(x)
      }
    }
    selector <- x[["w.length"]] >= range[1] & x[["w.length"]] <= range[2]

    if (is.na(range.s.data[1])) {
      range.s.data[1] <- min(x[["Rfr"]], na.rm = TRUE)
    }
    if (is.na(range.s.data[2])) {
      range.s.data[2] <- max(x[["Rfr"]], na.rm = TRUE)
    }
    x[selector, "Rfr"] <- with(x[selector, ],
                                        ifelse(Rfr < range.s.data[1],
                                               fill[1],
                                               ifelse(Rfr  > range.s.data[2],
                                                      fill[2],
                                                      Rfr)))
    x
  }


#' @describeIn clean Replace off-range values in a response spectrum
#'
#' @export
#'
clean.response_spct <-
  function(x,
           range = x,
           range.s.data = c(0,NA),
           fill = range.s.data,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {
    stopifnot(length(range) >= 2L &&
                length(range.s.data) == 2L &&
                length(fill) <= 2L)
    # wavelength range
    if (is.any_spct(range) || is.numeric(range) && length(range) > 2L) {
      range <- range(range, na.rm = TRUE)
    } else {
      if (is.na(range[1])) {
        range[1] <- min(x)
      }
      if (is.na(range[2])) {
        range[2] <- max(x)
      }
    }
    selector <- x[["w.length"]] >= range[1] & x[["w.length"]] <= range[2]

    if (unit.out == "quantum") {
      unit.out <- "photon"
    }
    if (unit.out == "photon") {
      x <- e2q(x, action = "replace")
      if (is.na(range.s.data[1])) {
        range.s.data[1] <- min(x[["s.q.response"]], na.rm = TRUE)
     }
      if (is.na(range.s.data[2])) {
        range.s.data[2] <- max(x[["s.q.response"]], na.rm = TRUE)
      }
      x[selector, "s.q.response"] <- with(x[selector, ],
                                       ifelse(s.q.response < range.s.data[1],
                                              fill[1],
                                              ifelse(s.q.response  > range.s.data[2],
                                                     fill[2],
                                                     s.q.response)))
    } else if (unit.out == "energy") {
      x <- q2e(x, action = "replace")
      if (is.na(range.s.data[1])) {
        range.s.data[1] <- min(x[["s.e.response"]], na.rm = TRUE)
      }
      if (is.na(range.s.data[2])) {
        range.s.data[2] <- max(x[["s.e.response"]], na.rm = TRUE)
      }
      x[selector, "s.e.response"] <- with(x[selector, ],
                                       ifelse(s.e.response < range.s.data[1],
                                              fill[1],
                                              ifelse(s.e.response  > range.s.data[2],
                                                     fill[2],
                                                     s.e.response)))
    } else {
      stop("unit.out: '", unit.out, "' not supported")
    }
    x
  }

#' @describeIn clean
#'
#' @export
#'
clean.source_mspct <-
  function(x,
           range = x,
           range.s.data = c(0,NA),
           fill = range.s.data,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...) {
    msmsply(x,
            clean,
            range = range,
            f = f,
            range.s.data = range.s.data,
            fill = fill,
            unit.out = unit.out,
            ...)
  }

#' @describeIn clean
#'
#' @export
#'
clean.filter_mspct <-
  function(x,
           range = x,
           range.s.data = NULL,
           fill = range.s.data,
           qty.out = getOption("photobiology.filter.qty",
                               default = "transmittance"),
           ...) {
    msmsply(x,
            clean,
            range = range,
            f = f,
            range.s.data = range.s.data,
            fill = fill,
            qty.out = qty.out,
            ...)
  }

#' @describeIn clean
#'
#' @export
#'
clean.reflector_mspct <-
  function(x,
           range = x,
           range.s.data = c(0, 1),
           fill = range.s.data,
           ...) {
    msmsply(x,
            clean,
            range = range,
            f = f,
            range.s.data = range.s.data,
            fill = fill,
            ...)
  }

#' @describeIn clean
#'
#' @export
#'
clean.response_mspct <-
  function(x,
           range = x,
           range.s.data = c(0,NA),
           fill = range.s.data,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...) {
    msmsply(x,
            clean,
            range = range,
            f = f,
            range.s.data = range.s.data,
            fill = fill,
            unit.out = unit.out,
            ...)
  }
