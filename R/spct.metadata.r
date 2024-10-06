# This file contains definitions for all methods related to setting and
# accessing metadata that are not tightly tied to how computations are
# performed or data are plotted. In other words, ancillary matadata.

# when.measured ---------------------------------------------------------------

#' Set the "when.measured" attribute
#'
#' Function to set by reference the "when" attribute  of an existing
#' generic_spct or an object of a class derived from generic_spct.
#'
#' @param x a generic_spct object
#' @param when.measured,value POSIXct to add as attribute, or a list of POSIXct.
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return x
#' @note This method alters x itself by reference and in addition
#'   returns x invisibly. If x is not a generic_spct or an object of a class derived from
#'   generic_spct, x is not modified. If \code{when} is not a POSIXct object
#'   or \code{NULL} an error is triggered. A \code{POSIXct} describes an
#'   instant in time (date plus time-of-day plus time zone).
#'
#' @export
#' @family measurement metadata functions
#' @examples
#' my.spct <- sun.spct
#' when_measured(my.spct)
#' when_measured(my.spct) <- lubridate::ymd_hms("2020-01-01 08:00:00")
#' when_measured(my.spct)
#'
setWhenMeasured <- function(x, when.measured, ...) UseMethod("setWhenMeasured")

#' @rdname setWhenMeasured
#'
#' @export
#'
`when_measured<-` <- function(x, value) {
  setWhenMeasured(x, when.measured = value)
}

#' @describeIn setWhenMeasured default
#' @export
setWhenMeasured.default <- function(x, when.measured, ...) {
  warning("Default dummy method called.")
  invisible(x)
}

#' @describeIn setWhenMeasured generic_spct
#' @export
setWhenMeasured.generic_spct <-
  function(x,
           when.measured = lubridate::now(tzone = "UTC"),
           ...) {
    name <- substitute(x)
    if (!is.null(when.measured)) {
      if (!is.list(when.measured)) {
        when.measured <- list(when.measured)
      } else if (!length(when.measured) %in% c(1L, getMultipleWl(x))) {
        warning("Length of 'when.measured' does not match spectrum object")
      }
      if (all(sapply(when.measured, lubridate::is.instant))) {
        when.measured <-
          lapply(when.measured, lubridate::with_tz, tzone = "UTC")
      }
      if (is.list(when.measured) && length(when.measured) == 1) {
        when.measured <- when.measured[[1]]
      }
    }
    attr(x, "when.measured") <- when.measured
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

#' @describeIn setWhenMeasured summary_generic_spct
#' @export
#'
setWhenMeasured.summary_generic_spct <- setWhenMeasured.generic_spct

#' @describeIn setWhenMeasured data.frame
#' @export
#'
setWhenMeasured.data.frame <- setWhenMeasured.generic_spct

#' @describeIn setWhenMeasured generic_mspct
#' @export
setWhenMeasured.generic_mspct <-
  function(x,
           when.measured = lubridate::now(tzone = "UTC"),
           ...) {
    name <- substitute(x)
    stopifnot((lubridate::is.POSIXct(when.measured) && length(when.measured) == 1) ||
                is.list(when.measured))
    if (lubridate::is.POSIXct(when.measured) || length(when.measured) == 1) {
      if (is.list(when.measured)) {
        when.measured <- when.measured[[1]]
        stopifnot(lubridate::is.POSIXct(when.measured))
      }
      when <- lubridate::with_tz(when.measured, "UTC")
      x <- msmsply(mspct = x, .fun = setWhenMeasured, when.measured = when)
    } else if (length(when.measured) == length(x)) {
      for (i in seq_along(x)) {
        when <- when.measured[[i]]
        stopifnot(lubridate::is.POSIXct(when))
        when <- lubridate::with_tz(when, "UTC")
        x[[i]] <- setWhenMeasured(x[[i]], when.measured = when)
      }
    }
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)
  }

#' Get the "when.measured" attribute
#'
#' Function to read the "when.measured" attribute of an existing generic_spct
#' or a generic_mspct.
#'
#' @param x a generic_spct object
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return POSIXct An object with date and time.
#'
#' @note If x is not a \code{generic_spct} or an object of a derived class
#'   \code{NA} is returned.
#'
#' @export
#' @family measurement metadata functions
#' @examples
#'
#' when_measured(sun.spct)
#'
getWhenMeasured <- function(x, ...) UseMethod("getWhenMeasured")

#' @rdname getWhenMeasured
#'
#' @export
#'
when_measured <- getWhenMeasured

#' @describeIn getWhenMeasured default
#' @export
getWhenMeasured.default <- function(x, ...) {
  # we return an NA of class POSIXct
  suppressWarnings(lubridate::ymd_hms(NA_character_, tz = "UTC"))
}

#' @describeIn getWhenMeasured generic_spct
#'
#' @param as.df logical If TRUE return a data frame instead of a list, when
#'   the value stored in the attribute is a list.
#'
#' @export
#'
getWhenMeasured.generic_spct <- function(x, as.df = FALSE, ...) {
  when.measured <- attr(x, "when.measured", exact = TRUE)
  if (is.null(when.measured)) {
    when.measured <- lubridate::NA_POSIXct_
  } else if (lubridate::is.POSIXct(when.measured)) {
    when.measured <-
      as.POSIXct(when.measured, tz = "UTC", origin = lubridate::origin)
  } else if (as.df && is.list(when.measured)) {
    if (all(sapply(when.measured, lubridate::is.instant))) {
     when.measured <-
        tibble::tibble(spct.idx = names(when.measured),
                       when.measured = as.POSIXct(unlist(when.measured, use.names = FALSE),
                                                  tz = lubridate::tz(when.measured[[1]])))
    } else {
      when.measured <-
        tibble::tibble(spct.idx = names(when.measured),
                       when.measured = rep(lubridate::NA_POSIXct_,
                                           length(when.measured)))
    }
  }
  when.measured
}

#' @describeIn getWhenMeasured summary_generic_spct
#' @export
getWhenMeasured.summary_generic_spct <- getWhenMeasured.generic_spct

#' @describeIn getWhenMeasured data.frame
#' @export
getWhenMeasured.data.frame <- getWhenMeasured.generic_spct

#' @describeIn getWhenMeasured generic_mspct
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @note The method for collections of spectra returns the
#'   a tibble with the correct times in TZ = "UTC".
#' @export
getWhenMeasured.generic_mspct <- function(x,
                                          ...,
                                          idx = "spct.idx") {
  z <- msdply(mspct = x, .fun = getWhenMeasured, ..., idx = idx, col.names = "when.measured")
  z[["when.measured"]] <- lubridate::with_tz(z[["when.measured"]], "UTC")
  z
}

# where.measured ---------------------------------------------------------------

#' Set the "where.measured" attribute
#'
#' Function to set by reference the "where.measured" attribute  of an existing
#' generic_spct or an object of a class derived from generic_spct.
#'
#' @param x a generic_spct object
#' @param where.measured,value A one row data.frame such as returned by
#'   function \code{geocode} from package 'ggmap' for a location search.
#' @param lat numeric Latitude in decimal degrees North
#' @param lon numeric Longitude in decimal degrees West
#' @param address character Human readable address
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return x
#' @note This method alters x itself by reference and in addition
#'   returns x invisibly. If x is not a generic_spct or an object of a class derived from
#'   generic_spct, x is not modified. If \code{where} is not a POSIXct object
#'   or \code{NULL} an error is triggered. A \code{POSIXct} describes an
#'   instant in time (date plus time-of-day plus time zone). As expected
#'   passing \code{NULL} as argument for \code{where.measured} unsets the
#'   attribute.
#'
#' @export
#'
#' @family measurement metadata functions
#'
#' @examples
#'
#' my.spct <- sun.spct
#' where_measured(my.spct)
#' where_measured(my.spct) <- data.frame(lon = 0, lat = -60)
#' where_measured(my.spct)
#'
setWhereMeasured <-
  function(x, where.measured, lat, lon, address, ...) UseMethod("setWhereMeasured")

#' @rdname setWhereMeasured
#'
#' @export
#'
`where_measured<-` <- function(x, value) {
  setWhereMeasured(x, where.measured = value)
}

#' @describeIn setWhereMeasured default
#' @export
setWhereMeasured.default <- function(x,
                                     where.measured,
                                     lat,
                                     lon,
                                     address,
                                     ...) {
  x
}

#' @describeIn setWhereMeasured generic_spct
#' @export
setWhereMeasured.generic_spct <- function(x,
                                          where.measured = NA,
                                          lat = NA,
                                          lon = NA,
                                          address = NA,
                                          ...) {
  name <- substitute(x)
  if (!is.null(where.measured)) {
    if (is.atomic(where.measured) && all(is.na(where.measured))) {
      # replace missing geocode with a valid one
      # type conversion needed for NA
      where.measured <-
        validate_geocode(data.frame(lon = as.numeric(lon),
                                    lat = as.numeric(lat),
                                    address = as.character(address),
                                    stringsAsFactors = FALSE))
      stopifnot(is_valid_geocode(where.measured))
    } else if (is.list(where.measured) && !is.data.frame(where.measured)) {
      where.measured <- sapply(where.measured, validate_geocode)
      stopifnot(all(sapply(where.measured, is_valid_geocode)))
    } else {
      where.measured <- validate_geocode(where.measured)
      stopifnot(is_valid_geocode(where.measured))
    }
  }
  attr(x, "where.measured") <- where.measured
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' @describeIn setWhereMeasured summary_generic_spct
#'
#' @export
setWhereMeasured.summary_generic_spct <- setWhereMeasured.generic_spct

#' @describeIn setWhereMeasured data.frame
#'
#' @export
setWhereMeasured.data.frame <- setWhereMeasured.generic_spct

#' @describeIn setWhereMeasured generic_mspct
#' @note Method for collections of spectra recycles the location information
#'   only if it is of length one.
#' @export
setWhereMeasured.generic_mspct <- function(x,
                                           where.measured = NA,
                                           lat = NA,
                                           lon = NA,
                                           address = NA,
                                           ...) {
  name <- substitute(x)
  if (!is.null(where.measured)) {
    if (is.atomic(where.measured) && all(is.na(where.measured))) {
      # replace missing geocode with a valid one
      # type conversion needed for NA
      where.measured <- data.frame(lon = as.numeric(lon),
                                   lat = as.numeric(lat),
                                   address = as.character(address),
                                   stringsAsFactors = FALSE)
    } else if (!is_valid_geocode(where.measured)) {
      stop("Bad 'where.measured' argument of class: ", class(where.measured))
    }
  }
  if (is.null(where.measured) ||
      (is.data.frame(where.measured) && nrow(where.measured) == 1)) {
    x <- msmsply(mspct = x,
                 .fun = setWhereMeasured,
                 where.measured = where.measured)
  } else if (is.data.frame(where.measured) &&
             nrow(where.measured) == length(x)) {
    if (exists("spct.idx", where.measured)) {
      if (setequal(where.measured[["spct.idx"]], names(x))) {
        # we use name matching
        j <- which(colnames(where.measured) != "spct.idx")
        for (i in names(x)) {
          wm <- where.measured[where.measured[["spct.idx"]] == i, j]
          x[[i]] <- setWhereMeasured(x[[i]],
                                     where.measured = wm)
        }
      } else {
        stop("'spct-idx' values '", where.measured[["spct.idx"]],
             "' do not match names of spectra in collection.")
      }
    } else {
      # we match by position
      for (i in seq_along(x)) {
        x[[i]] <- setWhereMeasured(x[[i]], where.measured = where.measured[i, ])
      }
    }
  } else if (is.list(where.measured) && length(where.measured) == length(x)) {
    for (i in seq_along(x)) {
      x[[i]] <- setWhereMeasured(x[[i]], where.measured = where.measured[[i]])
    }
  } else {
    stop("Length of geocode must be either 1, or equal to the number of spectra.")
  }
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' Get the "where.measured" attribute
#'
#' Function to read the "where.measured" attribute of an existing generic_spct.
#'
#' @param x a generic_spct object
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return a data.frame with a single row and at least columns "lon" and "lat",
#'    unless expand is set to \code{FALSE}.
#'
#' @note If x is not a \code{generic_spct} or an object of a derived class
#'   \code{NA} is returned.
#'
#' @export
#'
#' @family measurement metadata functions
#'
#' @examples
#' where_measured(sun.spct)
#'
getWhereMeasured <- function(x, ...) UseMethod("getWhereMeasured")

#' @rdname getWhereMeasured
#'
#' @export
#'
where_measured <- getWhereMeasured

#' @describeIn getWhereMeasured default
#' @export
#'
getWhereMeasured.default <- function(x, ...) {
  na_geocode()
}

#' @describeIn getWhereMeasured generic_spct
#' @export
#'
getWhereMeasured.generic_spct <- function(x, ...) {
  where.measured <- attr(x, "where.measured", exact = TRUE)
  if (is.null(where.measured)) return(na_geocode())

  if (is.list(where.measured) && !is.data.frame(where.measured)) {
    x <- dplyr::bind_rows(where.measured)
  }
  if (!is.data.frame(where.measured)) {
    # need to handle invalid or missing attribute values
    where.measured <- na_geocode()
  }
  # needed to clean inconsistent values from previous versions
  validate_geocode(where.measured)
}

#' @describeIn getWhereMeasured summary_generic_spct
#' @export
getWhereMeasured.summary_generic_spct <- getWhereMeasured.generic_spct

#' @describeIn getWhereMeasured generic_mspct
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @param .bind.geocodes logical In the case of collections of spectra if
#'    \code{.bind.geocodes = TRUE}, the default, the returned value is a single
#'    geocode with one row for each member spectrum. Otherwise the individual
#'    geocode data frames are returned in a list column within a tibble.
#'
#' @export
#'
getWhereMeasured.generic_mspct <- function(x,
                                           ...,
                                           idx = "spct.idx",
                                           .bind.geocodes = TRUE) {
  if (.bind.geocodes) {
    msdply(mspct = x, .fun = getWhereMeasured, idx = idx, ...)
  } else {
    l <- mslply(mspct = x, .fun = getWhereMeasured, ...)
    comment(l) <- NULL
    z <- list(where.measured = l)
    z[[idx]] <- factor(names(l), levels = names(l))
    tibble::as_tibble(z[c(2, 1)])
  }
}

#' @describeIn getWhereMeasured data.frame
#' @export
#'
getWhereMeasured.data.frame <- function(x, ...) {
  where.measured <- attr(x, "where.measured", exact = TRUE)
  if (is.null(where.measured)) return(na_geocode())
  stopifnot("The value of 'geocode' attribute is invalid" =
              is_valid_geocode(where.measured))

  if (is.list(where.measured) && !is.data.frame(where.measured)) {
    x <- dplyr::bind_rows(where.measured)
  }
  # needed to clean inconsistent values from previous versions
  validate_geocode(where.measured)
}

# how.measured attributes -------------------------------------------------

#' Set the "how.measured" attribute
#'
#' Function to set by reference the "how.measured" attribute  of an existing
#' generic_spct or derived-class object.
#'
#' @param x a generic_spct object
#' @param how.measured,value a list
#'
#' @return x
#' @note This function alters x itself by reference and in addition
#'   returns x invisibly. If x is not a generic_spct object, x is not
#'   modified.
#'
#' @export
#' @family measurement metadata functions
#'
#' @examples
#'
#' my.spct <- sun.spct
#' how_measured(my.spct)
#' how_measured(my.spct) <- "simulated with a radiation transfer model"
#' how_measured(my.spct)
#'
setHowMeasured <- function(x, how.measured) {
  name <- substitute(x)
  if (is.generic_spct(x) || is.summary_generic_spct(x) || is.data.frame(x)) {
    attr(x, "how.measured") <- how.measured
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' @rdname setHowMeasured
#'
#' @export
#'
`how_measured<-` <- function(x, value) {
  setHowMeasured(x, how.measured = value)
}

#' Get the "how.measured" attribute
#'
#' Function to read the "how.measured" attribute of an existing generic_spct
#' or a generic_mspct.
#'
#' @param x a generic_spct object
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return character vector An object containing a description of the data.
#'
#' @export
#' @family measurement metadata functions
#'
#' @examples
#' how_measured(sun.spct)
#'
getHowMeasured <- function(x, ...) UseMethod("getHowMeasured")

#' @rdname getHowMeasured
#'
#' @export
#'
how_measured <- getHowMeasured

#' @describeIn getHowMeasured default
#' @export
getHowMeasured.default <- function(x, ...) {
  # we return an NA of class character
  NA_character_
}

#' @describeIn getHowMeasured generic_spct
#' @export
getHowMeasured.generic_spct <- function(x, ...) {
  how.measured <- attr(x, "how.measured", exact = TRUE)
  if (is.null(how.measured) || (is.atomic(how.measured) && all(is.na(how.measured)))) {
    # need to handle objects created with old versions
    NA_character_
  } else {
    how.measured
  }
}

#' @describeIn getHowMeasured summary_generic_spct
#' @export
getHowMeasured.summary_generic_spct <- getHowMeasured.generic_spct

#' @describeIn getHowMeasured data.frame
#' @export
getHowMeasured.data.frame <- getHowMeasured.generic_spct

#' @describeIn getHowMeasured generic_mspct
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @note The method for collections of spectra returns the
#'   a tibble with a column of character strings.
#' @export
#'
getHowMeasured.generic_mspct <- function(x,
                                          ...,
                                          idx = "spct.idx") {
  msdply(mspct = x, .fun = getHowMeasured, ..., idx = idx, col.names = "how.measured")
}

##

#' Set the "instr.desc" attribute
#'
#' Function to set by reference the "instr.desc" attribute  of an existing
#' generic_spct or derived-class object.
#'
#' @param x a generic_spct object
#' @param instr.desc,value a list
#'
#' @return x
#' @note This function alters x itself by reference and in addition
#'   returns x invisibly. If x is not a generic_spct object, x is not
#'   modified.
#'
#' @note
#' The fields to be passed in the list \code{instr.desc} in part vary
#' depending on the instrument brand and model.
#'
#' @export
#'
#' @family measurement metadata functions
#'
setInstrDesc <- function(x, instr.desc) {
  name <- substitute(x)
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    attr(x, "instr.desc") <- instr.desc
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' @rdname setInstrDesc
#'
#' @export
#'
`instr_descriptor<-` <- function(x, value) {
  setInstrDesc(x, instr.desc = value)
}

#' Get the "instr.desc" attribute
#'
#' Function to read the "instr.desc" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#'
#' @return list (depends on instrument type)
#'
#'
#' @export
#' @family measurement metadata functions
#'
getInstrDesc <- function(x) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    if (isValidInstrDesc(x)) {
      instr.desc <- attr(x, "instr.desc", exact = TRUE)
    } else {
      instr.desc <- list(spectrometer.name = NA_character_,
                         spectrometer.sn = NA_character_,
                         bench.grating = NA_character_,
                         bench.slit = NA_character_)
    }
    if (!inherits(instr.desc, "instr_desc") &&
        !inherits(instr.desc[[1]], "instr_desc")) {
      class(instr.desc) <- c("instr_desc", class(instr.desc))
    }
    instr.desc
  } else {
    list()
  }
}

#' @rdname getInstrDesc
#'
#' @export
#'
instr_descriptor <- getInstrDesc

#' Trim the "instr.desc" attribute
#'
#' Function to trim the "instr.desc" attribute of an existing generic_spct
#' object, discarding all fields except for `spectrometer.name`,
#' `spectrometer.sn`, `bench.grating`, `bench.slit`, and calibration name.
#'
#' @param x a generic_spct object
#' @param fields a character vector with the names of the fields to keep,
#'   or if first member is `"-"`, the names of fields to delete; "*" as
#'   first member of the vector makes the function a no-op, leaving the spectrum
#'   object unaltered.
#'
#' @return x
#' @note This function alters x itself by reference and in addition
#'   returns x invisibly. If x is not a generic_spct object, x is not
#'   modified.
#'
#' @export
#' @family measurement metadata functions
#'
trimInstrDesc <- function(x,
                          fields = c("time",
                                     "spectrometer.name",
                                     "spectrometer.sn",
                                     "bench.grating",
                                     "bench.slit",
                                     "entrance.optics")
) {
  name <- substitute(x)
  if ((is.generic_spct(x) || is.summary_generic_spct(x)) &&
      fields[1] != "*") {
    instr.desc <- attr(x, "instr.desc", exact = TRUE)
    if (inherits(instr.desc, "instr_desc") ||
        "spectrometer.name" %in% names(instr.desc)) {
      instr.desc <- list(instr.desc)
    }
    for (i in seq(along.with = instr.desc)) {
      if (!(is.null(instr.desc[[i]]) || all(is.na(instr.desc[[i]])))) {
        if (fields[1] == "-") {
          fields.tmp <- setdiff(names(instr.desc[[i]]), fields[-1])
        } else if (fields[1] == "=") {
          fields.tmp <- fields[-1]
        } else {
          fields.tmp <- fields
        }
        instr.desc[[i]] <- instr.desc[[i]][fields.tmp]
        if (!inherits(instr.desc[[i]], "instr_desc")) {
          class(instr.desc[[i]]) <- c("instr_desc", class(instr.desc[[i]]))
        }
      }
    }
    if (length(instr.desc) == 1) {
      instr.desc <- instr.desc[[1]]
    }
    attr(x, "instr.desc") <- instr.desc
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' Check the "instr.desc" attribute
#'
#' Function to validate the "instr.settings" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#'
#' @return logical TRUE if at least instrument name and serial number is found.
#'
#' @export
#'
#' @family measurement metadata functions
#'
isValidInstrDesc <- function(x) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    instr.desc <- attr(x, "instr.desc", exact = TRUE)
    if (is.null(instr.desc)) {
      return(FALSE)
    }
    if (inherits(instr.desc, "instr_desc") ||
        "spectrometer.name" %in% names(instr.desc)) {
      # need to guard in case of objects created with earlier
      # versions
      instr.desc <- list(instr.desc)
    }
    valid <- TRUE
    for (desc in instr.desc) {
      if (length(desc) == 0 || (length(desc) == 1 && is.na(desc))) {
        # need to handle objects created with old versions
        valid <- FALSE
      } else if (is.list(desc)) {
        valid <- valid &&
          length(intersect(names(desc),
                           c("spectrometer.name", "spectrometer.sn"))) != 0
      } else {
        valid <- FALSE
      }
    }
  } else {
    valid <- NA_integer_
  }
  valid
}

#' Set the "instr.settings" attribute
#'
#' Function to set by reference the "what.measured" attribute  of an existing
#' generic_spct or derived-class object.
#'
#' @param x a generic_spct object
#' @param instr.settings,value a list
#'
#' @return x
#' @note This function alters x itself by reference and in addition
#'   returns x invisibly. If x is not a generic_spct object, x is not
#'   modified.
#'
#' @export
#' @family measurement metadata functions
#'
setInstrSettings <- function(x, instr.settings) {
  name <- substitute(x)
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    attr(x, "instr.settings") <- instr.settings
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' @rdname setInstrSettings
#'
#' @export
#'
`instr_settings<-` <- function(x, value) {
  setInstrSettings(x, instr.settings = value)
}

#' Get the "instr.settings" attribute
#'
#' Function to read the "instr.settings" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#'
#' @return list
#'
#'
#' @export
#'
#' @family measurement metadata functions
#'
getInstrSettings <- function(x) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    if (isValidInstrSettings(x)) {
      instr.settings <- attr(x, "instr.settings", exact = TRUE)
    } else {
      instr.settings <- list(integ.time = NA_real_,
                             tot.time = NA_real_,
                             num.scans = NA_integer_,
                             rel.signal = NA_real_)
    }
    if (!inherits(instr.settings, "instr_settings") &&
        !inherits(instr.settings[[1]], "instr_settings")) {
      class(instr.settings) <- c("instr_settings", class(instr.settings))
    }
    instr.settings
  } else {
    list()
  }
}

#' @rdname getInstrSettings
#'
#' @export
#'
instr_settings <- getInstrSettings

#' Trim the "instr.settings" attribute
#'
#' Function to trim the "instr.settings" attribute of an existing generic_spct
#' object, by discarding some fields.
#'
#' @param x a generic_spct object
#' @param fields a character vector with the names of the fields to keep,
#'   or if first member is `"-"`, the names of fields to delete; "*" as
#'   first member of the vector makes the function a no-op, leaving the spectrum
#'   object unaltered.
#'
#' @return x
#' @note This function alters x itself by reference and in addition
#'   returns x invisibly. If x is not a generic_spct object, x is not
#'   modified.
#'
#' @export
#' @family measurement metadata functions
#'
trimInstrSettings <- function(x,
                              fields = "*" ) {
  name <- substitute(x)
  if ((is.generic_spct(x) || is.summary_generic_spct(x)) &&
      fields[1] != "*") {
    instr.settings <- attr(x, "instr.settings", exact = TRUE)
    if (inherits(instr.settings, "instr_settings") ||
        "integ.time" %in% names(instr.settings)) {
      instr.settings <- list(instr.settings)
    }
    for (i in seq(along.with = instr.settings)) {
      if (!(length(instr.settings[[i]]) == 0 || all(is.na(instr.settings[[i]])))) {
        if (fields[1] == "-") {
          fields.tmp <- setdiff(names(instr.settings[[i]]), fields[-1])
        } else if (fields[1] == "=") {
          fields.tmp <- fields[-1]
        } else {
          fields.tmp <- fields
        }
        instr.settings[[i]] <- instr.settings[[i]][fields.tmp]
        if (!inherits(instr.settings[[i]], "instr_settings")) {
          class(instr.settings[[i]]) <- c("instr_settings", class(instr.settings[[i]]))
        }
      }
    }
    if (length(instr.settings) == 1) {
      instr.settings <- instr.settings[[1]]
    }
    attr(x, "instr.settings") <- instr.settings
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' Check the "instr.settings" attribute
#'
#' Function to validate the "instr.settings" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#'
#' @return logical TRUE if at least integration time data is found.
#'
#' @export
#'
#' @family measurement metadata functions
#'
isValidInstrSettings <- function(x) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    instr.settings <- attr(x, "instr.settings", exact = TRUE)
    if (is.null(instr.settings)) {
      return(FALSE)
    }
    if (inherits(instr.settings, "instr_settings") ||
        "integ.time" %in% names(instr.settings)) {
      # need to guard in case of objects created with earlier
      # versions
      instr.settings <- list(instr.settings)
    }
    valid <- TRUE
    for (setting in instr.settings) {
      if (length(setting) == 0 || (length(setting) == 1 && is.na(setting))) {
        # need to handle objects created with old versions
        valid <- FALSE
      } else if (is.list(setting)) {
        integ.time <- setting[["integ.time"]]
        if (is.null(integ.time) || any(is.na(integ.time)) || !is.numeric(integ.time)) {
          valid <- FALSE
        } # else we keep valid unchanged
      } else {
        valid <- FALSE
      }
    }
  } else {
    valid <- NA_integer_
  }
  valid
}

# what measured attributes -------------------------------------------------

#' Set the "what.measured" attribute
#'
#' Function to set by reference the "what.measured" attribute  of an existing
#' generic_spct or derived-class object.
#'
#' @param x a generic_spct object
#' @param what.measured,value a list
#'
#' @return x
#' @note This function alters x itself by reference and in addition
#'   returns x invisibly. If x is not a generic_spct object, x is not
#'   modified.
#'
#' @export
#'
#' @examples
#' my.spct <- sun.spct
#' what_measured(my.spct)
#' what_measured(my.spct) <- "Sun"
#' what_measured(my.spct)
#'
#' @family measurement metadata functions
#'
setWhatMeasured <- function(x, what.measured) {
  name <- substitute(x)
  if (is.generic_spct(x) || is.summary_generic_spct(x) || is.data.frame(x)) {
    attr(x, "what.measured") <- what.measured
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

#' @rdname setWhatMeasured
#'
#' @export
#'
`what_measured<-` <- function(x, value) {
  setWhatMeasured(x, what.measured = value)
}

#' Get the "what.measured" attribute
#'
#' Function to read the "what.measured" attribute of an existing generic_spct
#' or a generic_mspct.
#'
#' @param x a generic_spct object
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return character vector An object containing a description of the data.
#'
#' @export
#'
#' @family measurement metadata functions
#'
#' @examples
#'
#' what_measured(sun.spct)
#'
getWhatMeasured <- function(x, ...) UseMethod("getWhatMeasured")

#' @rdname getWhatMeasured
#'
#' @export
#'
what_measured <- getWhatMeasured

#' @describeIn getWhatMeasured default
#' @export
getWhatMeasured.default <- function(x, ...) {
  # we return an NA of class character
  NA_character_
}

#' @describeIn getWhatMeasured generic_spct
#' @export
getWhatMeasured.generic_spct <- function(x, ...) {
  what.measured <- attr(x, "what.measured", exact = TRUE)
  if (is.null(what.measured) || (is.atomic(what.measured) && all(is.na(what.measured)))) {
    # need to handle objects created with old versions
    NA_character_
  } else {
    what.measured
  }
}

#' @describeIn getWhatMeasured summary_generic_spct
#' @export
getWhatMeasured.summary_generic_spct <- getWhatMeasured.generic_spct

#' @describeIn getWhatMeasured data.frame
#' @export
getWhatMeasured.data.frame <- getWhatMeasured.generic_spct

#' @describeIn getWhatMeasured generic_mspct
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @note The method for collections of spectra returns the
#'   a tibble with a column of character strings.
#' @export
#'
getWhatMeasured.generic_mspct <- function(x,
                                          ...,
                                          idx = "spct.idx") {
  msdply(mspct = x, .fun = getWhatMeasured, ..., idx = idx, col.names = "what.measured")
}

# utility functions for attributes ----------------------------------------

#' Copy attributes from members of a generic_mspct
#'
#' Copy metadata attributes from members of a generic_mspct object into a tibble
#' or data.frame.
#'
#' @param mspct generic_mspct Any collection of spectra.
#' @param tb tibble or data.frame to which to add the data (optional).
#' @param col.names named character vector Name(s) of metadata attributes
#'   to copy, while if named, the names provide the name for the column.
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @param unnest logical Flag controlling if metadata attributes that are lists
#'   of values should be returned in a list column or in separate columns.
#'
#' @return A tibble With the metadata attributes in separate new variables.
#'
#' @details The attributes are copied to a column in a tibble or data frame. If
#'   the \code{tb} formal parameter receives \code{NULL} as argument, a new
#'   \code{tibble} will be created. If an existing \code{data.frame} or
#'   \code{tibble} is passed as argument, new columns are added to it. However,
#'   the number of rows in the argument passed to \code{tb} must match the
#'   number of spectra in the argument passed to \code{mspct}. Only in the case
#'   of method \code{add_attr2tb()} if the argument
#'   to \code{col.names} is a named vector, the names of members are used as names for the columns
#'   created. This permits setting any valid name for the new columns. If the
#'   vector passed to \code{col.names} has no names the names of the attributes
#'   are used for the new columns. If the fields of the attributes are unnested
#'   their names are used as names for the columns.
#'
#'   Valid accepted as argument to \code{col.names} are \code{NULL},
#'   \code{"lon"}, \code{"lat"}, \code{"address"}, \code{"geocode"},
#'   \code{"where.measured"}, \code{"when.measured"}, \code{"what.measured"},
#'   \code{"how.measured"}, \code{"comment"}, \code{"normalised"},
#'   \code{"normalized"}, \code{"scaled"}, \code{"bswf.used"},
#'   \code{"instr.desc"}, \code{"instr.settings"},
#'   \code{solute.properties}, \code{"filter.properties"},
#'   \code{"Tfr.type"}, \code{"Rfr.type"}, \code{"time.unit"}.
#'
#' @note The order of the first two arguments
#'   is reversed in \code{add_attr2tb()} compared to the other functions. This
#'   is to allow its use in 'pipes', while the functions for single attributes
#'   are expected to be used mostly to create new tibbles.
#'
#' @family measurement metadata functions
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
add_attr2tb <- function(tb = NULL,
                        mspct,
                        col.names = NULL,
                        idx = "spct.idx",
                        unnest = FALSE) {
  stopifnot(is.generic_mspct(mspct))

  force(col.names)
  if (length(col.names) < 1L) {
    return(tb)
  }
  if (all(is.na(col.names))) {
    return(tb)
  } else {
    col.names <- na.omit(col.names)
  }
  names.out <- names(col.names)
  if (is.null(names.out)) {
    # set names
    names(col.names) <- col.names
  } else {
    # fill-in only missing names
    selector <- names.out == ""
    names(col.names)[selector] <- col.names[selector]
  }
  if (unnest && any(c("geocode", "where.measured") %in% col.names)) {
    # setdiff removes names from the vector!
    col.names <- col.names[!col.names %in% c("lat", "lon")]
  }
  # We walk the list of attributes adding columns
  tb.cols <- names(tb)
  for (a in names(col.names)) {
    tb <-
      switch(a,
             lon = lon2tb(mspct = mspct,
                          tb = tb,
                          col.names = col.names["lon"],
                          idx = idx),
             lat = lat2tb(mspct = mspct,
                          tb = tb,
                          col.names = col.names["lat"],
                          idx = idx),
             address = address2tb(mspct = mspct,
                                  tb = tb,
                                  col.names = col.names["address"],
                                  idx = idx),
             geocode = geocode2tb(mspct = mspct,
                                  tb = tb,
                                  col.names = col.names["geocode"],
                                  idx = idx),
             where.measured = geocode2tb(mspct = mspct,
                                         tb = tb,
                                         col.names = col.names["where.measured"],
                                         idx = idx),
             when.measured = when_measured2tb(mspct = mspct,
                                              tb = tb,
                                              col.names = col.names["when.measured"],
                                              idx = idx),
             what.measured = what_measured2tb(mspct = mspct,
                                              tb = tb,
                                              col.names = col.names["what.measured"],
                                              idx = idx),
             how.measured = how_measured2tb(mspct = mspct,
                                            tb = tb,
                                            col.names = col.names["how.measured"],
                                            idx = idx),
             comment = comment2tb(mspct = mspct,
                                  tb = tb,
                                  col.names = col.names["comment"],
                                  idx = idx),
             normalized = normalized2tb(mspct = mspct,
                                        tb = tb,
                                        col.names = col.names["normalized"],
                                        idx = idx),
             normalised = normalized2tb(mspct = mspct,
                                        tb = tb,
                                        col.names = col.names["normalised"],
                                        idx = idx),
             scaled = scaled2tb(mspct = mspct,
                                tb = tb,
                                col.names = col.names["scaled"],
                                idx = idx),
             instr.desc = instr_desc2tb(mspct = mspct,
                                        tb = tb,
                                        col.names = col.names["instr.desc"],
                                        idx = idx),
             instr.settings = instr_settings2tb(mspct = mspct,
                                                tb = tb,
                                                col.names = col.names["instr.settings"],
                                                idx = idx),
             filter.properties = filter_properties2tb(mspct = mspct,
                                                      tb = tb,
                                                      col.names = col.names["filter.properties"],
                                                      idx = idx),
             solute.properties = solute_properties2tb(mspct = mspct,
                                                      tb = tb,
                                                      col.names = col.names["solute.properties"],
                                                      idx = idx),
             Tfr.type = Tfr_type2tb(mspct = mspct,
                                    tb = tb,
                                    col.names = col.names["Tfr.type"],
                                    idx = idx),
             Rfr.type = Rfr_type2tb(mspct = mspct,
                                    tb = tb,
                                    col.names = col.names["Rfr.type"],
                                    idx = idx),
             time.unit = time_unit2tb(mspct = mspct,
                                      tb = tb,
                                      col.names = col.names["time.unit"],
                                      idx = idx),
             bswf.used = BSWF_used2tb(mspct = mspct,
                                      tb = tb,
                                      col.names = col.names["bswf.used"],
                                      idx = idx),
             {warning("Skipping unknown metada name: ", a);
               tb})
  }
  if (unnest) {
    list.cols <- colnames(tb)[sapply(tb, is.list)]
    # do not expand preexisting list columns
    list.cols <- setdiff(list.cols, tb.cols)
    # expand metadata fields into columns
    for (col in list.cols) {
      # handles lists of lists or lists of dataframes
      tb <- tidyr::unnest_wider(tb, tidyr::all_of(col))
    }
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
                             idx = NULL) {
  stopifnot((is.generic_spct(mspct) || is.generic_mspct(mspct)) &&
              length(col.names) == 1L)
  if (is.generic_spct(mspct)) {
    if (is.null(idx)) {
      idx <- getIdFactor(mspct)
      if (is.na(idx)) {
        idx <- "spct.idx"
      }
    }
    when.ls  <- getWhenMeasured(mspct, idx = idx)
    if (getMultipleWl(mspct) == 1) {
      when.ls <- list(spct.nn = when.ls)
    }
    when.tb <-
      tibble::tibble(names(when.ls),
                     as.POSIXct(unlist(when.ls, use.names = FALSE),
                                tz = "UTC", origin = lubridate::origin))
    names(when.tb) <- c(idx, col.names)
  } else if (is.generic_mspct(mspct)) {
    if (is.null(idx)) {
      idx <- "spct.idx"
    }
    when.tb <- getWhenMeasured(mspct, idx = idx)
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
geocode2tb <- function(mspct,
                       tb = NULL,
                       col.names = "geocode",
                       idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  where.tb <- getWhereMeasured(mspct, idx = idx, .bind.geocodes = FALSE)
  names(where.tb)[2L] <- col.names
  if (is.null(tb)) {
    where.tb
  } else {
    dplyr::full_join(tb, where.tb, by = idx)
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
  names(lonlat.tb)[2L:3L] <- col.names
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
  names(lon.tb)[2L] <- col.names
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
  names(lat.tb)[2L] <- col.names
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
address2tb <- function(mspct,
                       tb = NULL,
                       col.names = "address",
                       idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  address.tb <- getWhereMeasured(mspct, idx = idx)[c(idx, "address")]
  names(address.tb)[2L] <- col.names
  if (is.null(tb)) {
    address.tb
  } else {
    dplyr::full_join(tb, address.tb, by = idx)
  }
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
  names(what.tb)[2L] <- col.names
  if (is.null(tb)) {
    what.tb
  } else {
    dplyr::full_join(tb, what.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
how_measured2tb <- function(mspct,
                            tb = NULL,
                            col.names = "how.measured",
                            idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  how.tb <- getHowMeasured(mspct, idx = idx)
  names(how.tb)[2L] <- col.names
  if (is.null(tb)) {
    how.tb
  } else {
    dplyr::full_join(tb, how.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
normalized2tb <- function(mspct,
                          tb = NULL,
                          col.names = "normalized",
                          idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  # method not implemented yet for collections
  normalized.tb <- msdply(mspct = mspct,
                          .fun = getNormalized,
                          idx = idx,
                          .force.numeric = TRUE)
  names(normalized.tb)[2L] <- col.names
  if (is.null(tb)) {
    normalized.tb
  } else {
    dplyr::full_join(tb, normalized.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
scaled2tb <- function(mspct,
                      tb = NULL,
                      col.names = "scaled",
                      idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  # method not implemented yet for collections
  l <- mslply(mspct = mspct, .fun = getScaled, .force.list = TRUE)
  comment(l) <- NULL
  z <- list(instr.desc = l)
  z[[idx]] <- factor(names(l), levels = names(l))
  fscaled.tb <- tibble::as_tibble(z[c(2, 1)])
  names(fscaled.tb)[2L] <- col.names
  if (is.null(tb)) {
    fscaled.tb
  } else {
    dplyr::full_join(tb, fscaled.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
instr_desc2tb <- function(mspct,
                          tb = NULL,
                          col.names = "instr.desc",
                          idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  # method not implemented yet for collections
  l <- mslply(mspct = mspct, .fun = getInstrDesc)
  comment(l) <- NULL
  z <- list(instr.desc = l)
  z[[idx]] <- factor(names(l), levels = names(l))
  desc.tb <- tibble::as_tibble(z[c(2, 1)])
  names(desc.tb)[2L] <- col.names
  if (is.null(tb)) {
    desc.tb
  } else {
    dplyr::full_join(tb, desc.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
instr_settings2tb <- function(mspct,
                              tb = NULL,
                              col.names = "instr.settings",
                              idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  # method not implemented yet for collections
  l <- mslply(mspct = mspct, .fun = getInstrSettings)
  comment(l) <- NULL
  z <- list(instr.settings = l)
  z[[idx]] <- factor(names(l), levels = names(l))
  settings.tb <- tibble::as_tibble(z[c(2, 1)])
  names(settings.tb)[2L] <- col.names
  if (is.null(tb)) {
    settings.tb
  } else {
    dplyr::full_join(tb, settings.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
BSWF_used2tb <- function(mspct,
                         tb = NULL,
                         col.names = "BSWF.used",
                         idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  # method not implemented yet for collections
  bswf.tb <- msdply(mspct = mspct, .fun = getBSWFUsed)
  names(bswf.tb)[2L] <- col.names
  if (is.null(tb)) {
    bswf.tb
  } else {
    dplyr::full_join(tb, bswf.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
filter_properties2tb <- function(mspct,
                                 tb = NULL,
                                 col.names = "filter.properties",
                                 idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  properties.tb <- getFilterProperties(mspct, idx = idx)
  names(properties.tb)[2L] <- col.names
  if (is.null(tb)) {
    properties.tb
  } else {
    dplyr::full_join(tb, properties.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
solute_properties2tb <- function(mspct,
                                 tb = NULL,
                                 col.names = "solute.properties",
                                 idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  properties.tb <- getSoluteProperties(mspct, idx = idx)
  names(properties.tb)[2L] <- col.names
  if (is.null(tb)) {
    properties.tb
  } else {
    dplyr::full_join(tb, properties.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
Tfr_type2tb <- function(mspct,
                        tb = NULL,
                        col.names = "Tfr.type",
                        idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  # method not implemented yet for collections
  tfr_type.tb <-  msdply(mspct = mspct, .fun = getTfrType)
  names(tfr_type.tb)[2L] <- col.names
  if (is.null(tb)) {
    tfr_type.tb
  } else {
    dplyr::full_join(tb, tfr_type.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
Rfr_type2tb <- function(mspct,
                        tb = NULL,
                        col.names = "Rfr.type",
                        idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  # method not implemented yet for collections
  rfr_type.tb <- msdply(mspct = mspct, .fun = getRfrType)
  names(rfr_type.tb)[2L] <- col.names
  if (is.null(tb)) {
    rfr_type.tb
  } else {
    dplyr::full_join(tb, rfr_type.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
time_unit2tb <- function(mspct,
                         tb = NULL,
                         col.names = "time.unit",
                         idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  # method not implemented yet for collections
  time_unit.tb <- msdply(mspct = mspct, .fun = getTimeUnit, force.duration = TRUE)
  names(time_unit.tb)[2L] <- col.names
  if (is.null(tb)) {
    time_unit.tb
  } else {
    dplyr::full_join(tb, time_unit.tb, by = idx)
  }
}

#' @rdname add_attr2tb
#'
#' @export
#'
comment2tb <- function(mspct,
                       tb = NULL,
                       col.names = "comment",
                       idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  # method not implemented yet for collections
  l <- mslply(mspct = mspct, .fun = comment)
  comment(l) <- NULL
  z <- list(instr.settings = l)
  z[[idx]] <- factor(names(l), levels = names(l))
  comments.tb <- tibble::as_tibble(z[c(2, 1)])
  #  settings.tb <- getInstrSettings(mspct, idx = idx)
  names(comments.tb)[2L] <- col.names
  if (is.null(tb)) {
    comments.tb
  } else {
    dplyr::full_join(tb, comments.tb, by = idx)
  }
}

# get all metadata --------------------------------------------------------

#' Access metadata
#'
#' Return metadata attributes from a single spectrum or a collection of spectra
#' as a tibble.
#'
#' @param x generic_mspct or generic_spct Any collection of spectra or spectrum.
#' @param col.names named character vector Name(s) of column(s) to create.
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @param na.rm logical Flag controlling deletion of columns containing only NA
#'   values.
#' @param unnest logical Flag controlling if metadata attributes that are lists
#'   of values should be returned in a list column or in separate columns.
#'
#' @return A tibble With the metadata attributes and an index column.
#'
#' @details Attributes are returned as columns in a tibble. If the argument to
#'   \code{col.names} is a named vector, with the names of members matching the
#'   names of attributes, then the values are used as names for the columns
#'   created. This permits setting any valid name for the new columns. If the
#'   vector passed to \code{col.names} has no names, then the values are
#'   interpreted as the names of the attributes to add, and also used as names
#'   for the new columns.
#'
#'   Some metadata values are stored in lists or data frames, these can be
#'   returned as a list columns or the individual fields unnested into separate
#'   columns.
#'
#' @seealso \code{\link{add_attr2tb}} for more details.
#'
#' @family measurement metadata functions
#'
#' @examples
#'
#' my.mspct <- source_mspct(list(sun1 = sun.spct, sun2 = sun.spct * 2))
#'
#' spct_metadata(my.mspct)
#'
#' spct_metadata(sun.spct)
#'
#' spct_metadata(my.mspct, na.rm = TRUE)
#'
#' spct_metadata(sun.spct, na.rm = TRUE)
#'
#' spct_metadata(my.mspct, col.names = c(geocode = "geo", "instr.desc"))
#'
#' spct_metadata(sun.spct, col.names = c(geocode = "geo", "instr.desc"))
#'
#' spct_metadata(sun.spct, col.names = "where.measured")$where.measured
#'
#' @export
#'
spct_metadata <- function(x,
                          col.names = NULL,
                          idx = "spct.idx",
                          na.rm = is.null(col.names),
                          unnest = TRUE) {
  force(na.rm) # compute default before assignment to col.names
  if (length(col.names) < 1L) {
    col.names <- c("where.measured",
                   "when.measured",
                   "what.measured",
                   "how.measured",
                   "normalized",
                   "scaled",
                   "time.unit",
                   "bswf.used",
                   "Tfr.type",
                   "Rfr.type")
  }
  if (is.any_spct(x)) {
    # ensure we operate on a collection of spectra
    name <- substitute(x)
    l <- list()
    l[[name]] <- x
    x <- generic_mspct(l, class = class(x)[1])
  }
  z <- add_attr2tb(tb = NULL,
                   mspct = x,
                   col.names = col.names,
                   idx = idx,
                   unnest = unnest)
  if (na.rm) {
    # omit columns with no data
    col.has.data <- sapply(X = as.list(z),
                           FUN = function(x) {
                             (is.atomic(x) & !all(is.na(x))) |
                               (is.list(x) & !all(is.na(unlist(x))))
                           })

    z <- z[ , col.has.data]
  }
  z
}

