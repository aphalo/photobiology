# This file contains definitions for all methods related to setting and
# accessing metadata that are not tightly tied to how computations are
# performed or data are plotted. In other words, ancillary metadata.

# when.measured ---------------------------------------------------------------

#' Set the "when.measured" attribute
#'
#' Method to set by reference the \code{"when.measured"} attribute  of an R
#' object.
#'
#' @param x an R object
#' @param when.measured,value POSIXct to add as attribute, or a list of POSIXct.
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return \code{x}, with its \code{"when.measured"} set.
#' @details This method alters \code{x} itself by reference and in addition
#'   returns \code{x} invisibly. If \code{x} is not an object of a supported
#'   class, \code{x} is not modified. If the arguments to \code{"when.measured"}
#'   or \code{value} are not a \code{POSIXct} object or \code{NULL} an error is
#'   triggered. A \code{POSIXct} describes an instant in time (date plus
#'   time-of-day plus time zone).
#'
#'   Be aware that \code{lubridate::ymd()} returns an incompatible \code{Date}
#'   object while \code{lubridate::ymd_h()}, \code{lubridate::ymd_hm()} and
#'   \code{lubridate::ymd_hms()} and similar functions return objects of class
#'   \code{POSIXct} acceptable as arguments for parameter \code{when.measured}.
#'
#' @export
#'
#' @family measurement metadata functions
#'
#' @examples
#' my.spct <- sun.spct
#' when_measured(my.spct)
#' when_measured(my.spct) <- lubridate::ymd_hms("2020-01-01 08:00:00")
#' when_measured(my.spct)
#' when_measured(my.spct) <- NULL
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
#' Method to read the \code{"when.measured"} attribute of an R object.
#'
#' @param x an R object
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return a \code{POSIXct} object with date and time, or named list of such
#'   objects, or, on user request, a data frame.
#'
#' @note If \code{x} is not an object of one of the supported classes,
#'   \code{NA} is returned.
#'
#' @export
#'
#' @inherit setWhenMeasured examples
#'
#' @family measurement metadata functions
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
#' @param as.df logical If \code{TRUE} return a data frame instead of a list,
#'   when the value stored in the attribute is a list.
#'
#' @export
#'
getWhenMeasured.generic_spct <- function(x, as.df = FALSE, ..., simplify = FALSE) {
  when.measured <- attr(x, "when.measured", exact = TRUE)
  if (is.null(when.measured)) {
    when.measured <- lubridate::NA_POSIXct_
  } else if (lubridate::is.POSIXct(when.measured)) {
    if (simplify && length(unique(when.measured)) == 1) {
      when.measured <- when.measured[1]
    }
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
    if (simplify &&
        sum(!duplicated(when.measured[ , -which(names(when.measured) == getIdFactor(x))])) == 1) {
      when.measured <- when.measured[1 , -which(names(when.measured) == getIdFactor(x))]
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
#' @param simplify logical If all members share the same attribute value return
#'   one copy instead of a data.frame.
#' @note The method for collections of spectra returns a tibble with the
#'   times expressed in TZ = "UTC".
#' @export
#' @examples
#' my.spct <- sun.spct
#' when_measured(my.spct)
#' when_measured(my.spct) <- lubridate::ymd_hms("2020-01-01 08:00:00")
#' when_measured(my.spct)
#' when_measured(my.spct) <- NULL
#' when_measured(my.spct)
#'
getWhenMeasured.generic_mspct <- function(x,
                                          ...,
                                          idx = "spct.idx",
                                          simplify = FALSE) {
  z <- msdply(mspct = x,
              .fun = getWhenMeasured,
              ...,
              idx = idx,
              col.names = "when.measured")
  z[["when.measured"]] <- lubridate::with_tz(z[["when.measured"]], "UTC")

  if (simplify) {
    zz <- unique(z[["when.measured"]])
    if (length(zz) <= 1) {
      z <- zz
    } else {
      z <- z[["when.measured"]]
    }
  }
  z
}

# where.measured ---------------------------------------------------------------

#' Set the "where.measured" attribute
#'
#' Method to set by reference the \code{"where.measured"} attribute  of an R
#' object.
#'
#' @param x an R object
#' @param where.measured,value A one row \code{data.frame} with the same format
#'   as returned by function \code{geocode} from package 'ggmap' for a location
#'   search.
#' @param lat numeric Latitude in decimal degrees North.
#' @param lon numeric Longitude in decimal degrees West.
#' @param address character Human readable address.
#' @param idFactor character Name of the column with IDs of the spectra stored
#'   in long form or ID column name in bound geocodes to use for IDs of
#'   collection of spectra members.
#' @param simplify logical If all members share the same geocode value set as
#'   attribute value a one row geocode instead of a named list of data frames.
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @details Code \code{setWhereMeasured()} methods alter \code{x} itself by
#'   reference and in addition return \code{x} invisibly. If \code{x} is not an
#'   object of a supported class, \code{x} is not modified. If the argument to
#'   \code{where.measured} is not a \code{data.frame} or \code{tibble} object or
#'   \code{NULL} an error is triggered as a validation test is applied. A
#'   geocode describes a geographic location based on longitude (\code{lon}) and
#'   latitude (\code{lat}) as \code{numeric} values and can optionally contain
#'   an address (\code{address}) as a single \code{character} string. Passing
#'   \code{NULL} as argument for parameter \code{where.measured} unsets the
#'   attribute. Parameters \code{lon}, \code{lat} and \code{address} provide an
#'   alternative to passing a ready constructed geocode data frame as input.
#'
#'   By default, when setting the geocode attribute for multiple spectra stored
#'   in long form, geocodes are stored as named lists of data frames unless they
#'   are identical and can be simplified. It is possible to
#'   disable simplification and force the use of a named list.
#'
#'   If the argument passed to parameter \code{geocode} is a data frame with one
#'   row per spectrum and the \code{idFactor} name matches, it will be split
#'   into a named list and, if possibly, simplified. It is also possible but
#'   deprecated to set the attribute to an indexed geocode with multiple rows.
#'
#' @return x, with the \code{"where.measured"} attribute set or unset.
#'
#' @export
#'
#' @family measurement metadata functions
#'
#' @examples
#' my.spct <- sun.spct
#' where_measured(my.spct)
#' where_measured(my.spct) <- data.frame(lon = 0, lat = -60)
#' where_measured(my.spct)
#' where_measured(my.spct) <- NULL
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
                                          idFactor = getIdFactor(x),
                                          simplify = TRUE,
                                          ...) {
  name <- substitute(x)
  if (!is.null(where.measured)) {
    if (is.atomic(where.measured) && all(is.na(where.measured))) {
      # replace missing geocode with a valid one
      # type conversion needed for NA
      where.measured <-
        SunCalcMeeus::validate_geocode(
          data.frame(lon = as.numeric(lon),
                     lat = as.numeric(lat),
                     address = as.character(address),
                     stringsAsFactors = FALSE))
      stopifnot(SunCalcMeeus::is_valid_geocode(where.measured))
    } else if (is.list(where.measured) && !is.data.frame(where.measured)) {
      where.measured <- lapply(where.measured, SunCalcMeeus::validate_geocode)
      stopifnot(all(sapply(where.measured, SunCalcMeeus::is_valid_geocode)))
      if (simplify && length(where.measured) == 1L) {
        where.measured <- where.measured[[1]]
      }
    } else {
      where.measured <-
        SunCalcMeeus::validate_geocode(where.measured)
#      stopifnot(SunCalcMeeus::is_valid_geocode(where.measured))
      if (getMultipleWl(x) > 1L &&
#          nrow(where.measured) > 1L &&
          !is.na(idFactor)) {
        where.measured <-
          SunCalcMeeus::split_geocodes(geocode = where.measured,
                                       idx = idFactor,
                                       simplify = simplify)
      }
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
#'   only if it is a one row \code{data.frame}.
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
      # replace missing geocode with a valid one built from other arguments
      # type conversion needed for NA
      where.measured <- data.frame(lon = as.numeric(lon),
                                   lat = as.numeric(lat),
                                   address = as.character(address),
                                   stringsAsFactors = FALSE)
    }
    if (!SunCalcMeeus::is_valid_geocode(where.measured)) {
      stop("Bad 'where.measured' argument. ",
           "Class: ", class(where.measured),
           "; named: ", names(where.measured),
           "; length: ", length(where.measured))
    }
  }
  if (is.null(where.measured) ||
      (is.data.frame(where.measured) && nrow(where.measured) == 1)) {
    # recycle and apply to each member
    x <- msmsply(mspct = x,
                 .fun = setWhereMeasured,
                 where.measured = where.measured)
  } else if (is.data.frame(where.measured) &&
             nrow(where.measured) == length(x)) {
    # match and split metadata to members
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
#' Method to read the "where.measured" attribute of
#' generic_spct, generic_mspct, summary_generic_spct, data.frame or a
#' derived-class object.
#'
#' @param x a generic_spct object
#' @param ... Allows use of additional arguments in methods for other classes.
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @param simplify logical If all members share the same attribute value return
#'   one copy instead of a data.frame.
#' @param .bind.geocodes logical In the case of collections of spectra if
#'    \code{.bind.geocodes = TRUE}, the default, the returned value is a single
#'    geocode with one row for each member spectrum. Otherwise the individual
#'    geocode data frames are returned in a list column within a tibble.
#'
#' @return a data.frame with a single row and at least columns "lon" and "lat",
#'    unless expand is set to \code{FALSE}.
#'
#' @note If x is not a \code{generic_spct} or an object of a derived class
#'   \code{NA} is returned.
#'
#' @export
#'
#' @examples
#' getWhereMeasured(sun.spct)
#' getWhereMeasured(sun_evening.spct)
#' getWhereMeasured(sun_evening.spct, simplify = TRUE)
#' getWhereMeasured(sun_evening.mspct)
#' getWhereMeasured(sun_evening.mspct, .bind.geocodes = FALSE)
#' getWhereMeasured(sun_evening.mspct, simplify = TRUE)
#' getWhereMeasured(sun_evening.mspct, simplify = FALSE)
#'
#' getWhereMeasured(summary(sun.spct))
#' getWhereMeasured(summary(sun_evening.spct))
#' getWhereMeasured(summary(sun_evening.mspct, expand = "each")[[1]])
#'
#' getWhereMeasured(polyester.spct)
#'
#' @family measurement metadata functions
#'
#' @inherit setWhereMeasured examples
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
  SunCalcMeeus::na_geocode()
}

#' @describeIn getWhereMeasured generic_spct
#' @export
#'
getWhereMeasured.generic_spct <- function(x,
                                          ...,
                                          idx = getIdFactor(x),
                                          simplify = TRUE,
                                          .bind.geocodes = TRUE) {
  where.measured <- attr(x, "where.measured", exact = TRUE)
  # attribute not set
  if (is.null(where.measured)) return(SunCalcMeeus::na_geocode())
  # single spectrum and not returning a list
  if ((getMultipleWl(x) == 1L) &&
      (simplify || .bind.geocodes || is.na(idx))) return(where.measured)
  # bind list of geocodes into a single data frame
  if (.bind.geocodes && !is.data.frame(where.measured)) {
    if (is.na(idx)) {
      idx <- "spct.idx"
    }
    SunCalcMeeus::bind_geocodes(where.measured,
                                idx = idx)
  } else if (!.bind.geocodes && is.data.frame(where.measured)) {
    # split a multi-row data frame into list of data frames
    if (idx %in% colnames(where.measured)) {
      SunCalcMeeus::split_geocodes(where.measured,
                                   idx = idx,
                                   simplify = FALSE)
    # single row geocode is already simplified
    } else if (simplify && nrow(where.measured) == 1L) {
      where.measured
    } else {
      warning("No 'idx' column \"", idx, "\" found! Returning 'where.measured' unchanged")
      where.measured
    }
  } else if (.bind.geocodes && is.data.frame(where.measured)) {
    # if a simplified geocode has an idx column, remove it
    if (nrow(where.measured) == 1L && idx %in% colnames(where.measured)) {
      where.measured[ , -which(colnames(where.measured) == idx)]
    } else if (nrow(where.measured) == 1L &&
               "spct.idx" %in% colnames(where.measured)) {
      where.measured[ , -which(colnames(where.measured) == "spct.idx")]
    } else {
      where.measured
    }
  } else if (!.bind.geocodes && is.list(where.measured)) {
    if (simplify && length(where.measured) == 1L) {
      where.measured[[1]]
    } else {
      if (!all(names(where.measured) %in% unique(x[[idx]]))) {
        warning("Unexpected 'where.measured' value! Returning it unchanged")
      }
      where.measured
    }
  }
}

#' @describeIn getWhereMeasured summary_generic_spct
#' @export
getWhereMeasured.summary_generic_spct <- getWhereMeasured.generic_spct

#' @describeIn getWhereMeasured generic_mspct
#'
#' @export
#'
getWhereMeasured.generic_mspct <- function(x,
                                           ...,
                                           idx = "spct.idx",
                                           .bind.geocodes = TRUE,
                                           simplify = FALSE) {
  if (.bind.geocodes) {
    z <- msdply(mspct = x, .fun = getWhereMeasured, idx = idx, ...)
    if (simplify) {
      unique.rows <- !duplicated(z[ , -which(names(z) == idx)])
      if (sum(unique.rows) <= 1) {
        z <- z[1, -which(names(z) == idx)]
      }
    }
  } else {
    z <- mslply(mspct = x, .fun = getWhereMeasured, ...)
    comment(z) <- NULL
    if (simplify) {
      if (length[z] == 1L) {
        z <- z[[1]]
      } else {
        num.unique <- 1L
        for (i in 2:length(z)) {
          num.unique <- num.unique + isFALSE(all.equal(z[[1]], z[[i]]))
        }
        if (num.unique == 1L) {
          z <- z[[1]]
        }
      }
    }
  }
  z
}

#' @describeIn getWhereMeasured data.frame
#' @export
#'
getWhereMeasured.data.frame <- function(x, ...) {
  where.measured <- attr(x, "where.measured", exact = TRUE)
  if (is.null(where.measured)) return(SunCalcMeeus::na_geocode())
  stopifnot("The value of 'geocode' attribute is invalid" =
              SunCalcMeeus::is_valid_geocode(where.measured))

  if (is.list(where.measured) && !is.data.frame(where.measured)) {
    x <- dplyr::bind_rows(where.measured)
  }
  # needed to clean inconsistent values from previous versions
  SunCalcMeeus::validate_geocode(where.measured)
}

# how.measured attributes -------------------------------------------------

#' Set the "how.measured" attribute
#'
#' Method to set the \code{"how.measured"} attribute of an R object.
#'
#' @param x a R object.
#' @param how.measured,value a list or a character string.
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return x modified by reference.
#'
#' @note This function alters \code{x} itself by reference and in addition
#'   returns \code{x} invisibly. If \code{x} is not an object of a supported
#'   class, \code{x} is silently returned unchanged.
#'
#' @export
#' @family measurement metadata functions
#'
#' @examples
#' my.spct <- sun.spct
#' how_measured(my.spct)
#' how_measured(my.spct) <- "Simulated with a radiation transfer model"
#' how_measured(my.spct)
#' how_measured(my.spct) <- NULL
#' how_measured(my.spct)
#'
setHowMeasured <- function(x, ...) {UseMethod("setHowMeasured")}


#' @rdname setHowMeasured
#'
#' @export
#'
`how_measured<-` <- function(x, value) {
  setHowMeasured(x, how.measured = value)
}

#' @describeIn setHowMeasured default
#' @export
setHowMeasured.default <- function(x,
                                   how.measured,
                                   ...) {
  x
}

#' @describeIn setHowMeasured generic_spct
#' @export
setHowMeasured.generic_spct <- function(x, how.measured, ...) {
  name <- substitute(x)
  attr(x, "how.measured") <- how.measured
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' @describeIn setHowMeasured summary_generic_spct
#' @export
setHowMeasured.summary_generic_spct <- setHowMeasured.generic_spct

#' @describeIn setHowMeasured data.frame
#' @export
setHowMeasured.data.frame <- setHowMeasured.generic_spct

#' @describeIn setHowMeasured generic_mspct
#'
#' @export
setHowMeasured.generic_mspct <- function(x,
                                         how.measured,
                                         ...) {
  msmsply(mspct = x, .fun = setHowMeasured, ..., how.measured = how.measured)
}

#' Get the "how.measured" attribute
#'
#' Method to read the \code{"how.measured"} attribute of an R object.
#'
#' @param x an R object.
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return character vector An object containing a verbal description of the
#'   data.
#'
#' @export
#'
#' @inherit setHowMeasured examples
#'
#' @family measurement metadata functions
#'
getHowMeasured <- function(x, ...) UseMethod("getHowMeasured")

#' @rdname getHowMeasured
#'
#' @export
#'
how_measured <- getHowMeasured

#' @describeIn getHowMeasured default
#'
#' @export
getHowMeasured.default <- function(x, ...) {
  # we return an NA of class character
  NA_character_
}

#' @describeIn getHowMeasured generic_spct
#'
#' @export
getHowMeasured.generic_spct <- function(x, ..., simplify = FALSE) {
  z <- attr(x, "how.measured", exact = TRUE)
  if (is.null(z) || (is.atomic(z) && all(is.na(z)))) {
    # need to handle objects created with old versions
    z <- NA_character_
  } else if (simplify && is.list(z)) {
    if (length(z) == 1L || sum(!duplicated(z) == 1L)) {
      z <- z[[1]]
    }
  }
  z
}

#' @describeIn getHowMeasured summary_generic_spct
#'
#' @export
getHowMeasured.summary_generic_spct <- getHowMeasured.generic_spct

#' @describeIn getHowMeasured data.frame
#'
#' @export
getHowMeasured.data.frame <- getHowMeasured.generic_spct

#' @describeIn getHowMeasured generic_mspct
#'
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @param simplify logical If all members share the same attribute value return
#'   one copy instead of a data.frame.
#' @note The method for collections of spectra returns the
#'   a data frame with a column of character strings.
#' @export
#'
getHowMeasured.generic_mspct <- function(x,
                                         ...,
                                         idx = "spct.idx",
                                         simplify = FALSE) {
  z <- msdply(mspct = x,
              .fun = getHowMeasured,
              ...,
              idx = idx,
              col.names = "how.measured")
  if (simplify) {
    zz <- unique(z[["how.measured"]])
    if (length(zz) <= 1) {
      z <- zz
    } else {
      z <- paste(z[["how.measured"]], "\n", collapse = "", sep = "")
    }
  }
  z
}


# Instrument descriptor attribute -----------------------------------------

#' Set the "instr.desc" attribute
#'
#' Function to set by reference the \code{"instr.desc"} attribute of an existing
#' \code{generic_spct} or derived-class object, or of a
#' \code{summary_generic_spct} or derived-class object.
#'
#' @param x a \code{generic_spct} object or a \code{summary_generic_spct}
#'  object.
#' @param instr.desc,value a \code{list}, \code{instr_desc} object, or
#'   \code{NULL}.
#'
#' @return \code{x}, with the value of its \code{"instr.desc"} attribute set to
#'   the value of the argument passed to \code{instr.desc} or to \code{value}.
#'
#' @details This function alters \code{x} itself by reference and in addition
#'   returns \code{x} invisibly. If \code{x} is not a \code{generic_spct} object,
#'   \code{x} is not modified, silently. If \code{inst.desc = NULL} is passed
#'   in the call, the attribute \code{"instr.desc"} is removed.
#'   \emph{This function is very rarely called from user code.}
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
    if (!is.null(instr.desc)) {
      if (length(instr.desc) == 1L) {
        if (is.na(instr.desc)) {
          instr.desc <- list()
        }
        stopifnot("The argument passed to 'instr.desc' must be a list" =
                    is.list(instr.desc))
        minimal.desc <- list(spectrometer.name = NA_character_,
                             spectrometer.sn = NA_character_,
                             bench.grating = NA_character_,
                             bench.slit = NA_character_)
        missing <- !c("spectrometer.name",
                      "spectrometer.sn",
                      "bench.grating",
                      "bench.slit") %in% names(instr.desc)
        if (any(missing)) {
          instr.desc <- c(instr.desc, minimal.desc[missing])
        }
        if (!inherits(instr.desc, "instr_desc") &&
            !inherits(instr.desc[[1]], "instr_desc")) {
          class(instr.desc) <- c("instr_desc", class(instr.desc))
        }
      # when all descriptors are the same the attr is not expanded into a list
      } else if (!is.list(instr.desc) ||
                 (length(instr.desc) != getMultipleWl(x) &&
                  !any(c("spectrometer.name", "spectrometer.sn") %in%
                       names(instr.desc)))) {
        warning("Length of 'instr.desc' list different to number of spectra: ",
                length(instr.desc), " != ", getMultipleWl(x))
      }
    }
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
#' Function to query the \code{"instr.desc"} attribute of an existing
#' \code{generic_spct} or derived-class object, or of a
#' \code{summary_generic_spct} or derived-class object.
#'
#' @param x a \code{generic_spct} object or a \code{summary_generic_spct}
#'   object.
#'
#' @return an object of class \code{"instr_desc"} derived from \code{"list"}.
#'   The fields \code{spectrometer.name}, \code{spectrometer.sn},
#'   \code{bench.grating} and \code{bench.slit} are always present, although may
#'   be set to \code{NA}. Additional fields can be present depending on the
#'   origin of the data.
#'
#' @export
#'
#' @examples
#' valid.descriptor <- getInstrDesc(white_led.cps_spct)
#' class(valid.descriptor)
#' print(valid.descriptor)
#' print(str(valid.descriptor))
#'
#' missing.descriptor <- getInstrDesc(white_body.spct)
#' class(missing.descriptor)
#' print(missing.descriptor)
#' print(str(missing.descriptor))
#'
#' @family measurement metadata functions
#'
getInstrDesc <- function(x) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    if (isValidInstrDesc(x)) {
      instr.desc <- attr(x, "instr.desc", exact = TRUE)
      # earlier bug created objects with extra empty fields with bad names
      bad.fields <- which(is.na(names(instr.desc)) |
                            names(instr.desc) %in% c("NA", ""))
      if (length(bad.fields)) {
        instr.desc <- instr.desc[-bad.fields]
      }
    } else {
      instr.desc <- list()
    }
    if (getMultipleWl(x) == 1) {
      minimal.desc <- list(spectrometer.name = NA_character_,
                           spectrometer.sn = NA_character_,
                           bench.grating = NA_character_,
                           bench.slit = NA_character_)
      missing <- !c("spectrometer.name",
                    "spectrometer.sn",
                    "bench.grating",
                    "bench.slit") %in% names(instr.desc)
      if (any(missing)) {
        instr.desc <- c(instr.desc, minimal.desc[missing])
      }
      # very old records lack class attribute
      if (!inherits(instr.desc, "instr_desc") &&
          !inherits(instr.desc[[1]], "instr_desc")) {
        class(instr.desc) <- c("instr_desc", class(instr.desc))
      }
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
#' Function to trim the \code{"instr.desc"} attribute of a \code{generic_spct}
#' or a \code{summary_generic_spct} object, by default discarding all fields
#' except for \code{spectrometer.name}, \code{spectrometer.sn},
#' \code{bench.grating}, \code{bench.slit}, and \code{entrance.optics}.
#'
#' @param x a \code{generic_spct} object or a \code{summary_generic_spct}
#'   object.
#' @param fields a character vector with the names of the fields to keep,
#'   or if first member is \code{"-"}, the names of fields to delete; \code{"*"}
#'   as the first member of the vector makes the function a no-op, leaving the
#'   spectrum object unaltered.
#'
#' @details This function alters \code{x} itself by reference and in addition
#'   returns \code{x} invisibly. If \code{x} is not a \code{generic_spct} object
#'   or a \code{summary_generic_spct} object, or if the \code{"instr.desc"}
#'   attribute is not present in a \code{generic_spct} object, \code{x} is not
#'   modified.
#'
#'   Attempts to remove or keep fields that are not present in the attribute are
#'   ignored silently. The value of fields in the attribute is never modified,
#'   fields are either kept unchanged or removed.
#'
#' @note Some of the spectrometer-specific metadata can be large, as they can
#'   include calibration coefficients. In the case of R package 'ooacquire' also
#'   pointers to Java objects may need to be deleted.
#'
#' @return \code{x}, possibly with the \code{"instr.desc"} attribute
#'   modified.
#'
#' @export
#'
#' @examples
#' my.spct <- white_led.cps_spct
#' names(instr_descriptor(my.spct))
#' trimInstrDesc(my.spct) # modified by reference!
#' names(instr_descriptor(my.spct))
#'
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
      # in earlier versions of 'photobiology' the record had fewer fields
      # thus, we silently ignore missing fields, avoiding NAs
      if (!(is.null(instr.desc[[i]]) || all(is.na(instr.desc[[i]])))) {
        if (fields[1] == "-") {
          fields.tmp <- setdiff(names(instr.desc[[i]]), fields[-1])
        } else if (fields[1] == "=") {
          fields.tmp <- intersect(fields[-1], names(instr.desc[[i]]))
        } else {
          fields.tmp <- intersect(fields, names(instr.desc[[i]]))
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
#' Function to validate the \code{"instr.settings"} attribute of an existing
#' \code{generic_spct} object or \code{summary_generic_spct} object.
#'
#' @param x a \code{generic_spct} object or a \code{summary_generic_spct}
#'   object.
#'
#' @details Test if at least one of instrument name (field
#'  \code{spectrometer.name}) or serial number (field \code{spectrometer.sn})
#'  is found in the value of the R attribute \code{"instr.desc"} of \code{x}.
#'  \code{FALSE} is silently returned if \code{x} does not belong to a class
#'  derived from class \code{generic_spct} or from class
#'  \code{summary_generic_spct}, or if it is derived from these classes but the
#'  attribute is not set.
#'
#' @return A \code{logical} vector of length one.
#'
#' @export
#'
#' @examples
#' isValidInstrDesc(white_led.cps_spct)
#' isValidInstrDesc(white_body.spct)
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

# Instrument settings attribute -------------------------------------------

#' Set the "instr.settings" attribute
#'
#' Function to set by reference the \code{"what.measured"} attribute  of a
#' \code{generic_spct}, or of a \code{summary_generic_spct} object.
#'
#' @param x a \code{generic_spct} object or a \code{summary_generic_spct}
#'   object.
#' @param instr.settings,value a \code{list} or a \code{instr_settings} object.
#'
#' @return x
#'
#' @details This function alters \code{x} itself by reference and in addition
#'   returns \code{x} invisibly. If \code{x} is not a \code{generic_spct} object
#'   or a \code{summary_generic_spct} object, \code{x} is not modified,
#'   silently. If \code{inst.desc = NULL} is passed in the call, the attribute
#'   \code{instr.settings} is removed. \emph{This function is very rarely called
#'   from user code.}
#'
#' @export
#' @family measurement metadata functions
#'
setInstrSettings <- function(x, instr.settings) {
  name <- substitute(x)
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    if (!is.null(instr.settings)) {
      if (getMultipleWl(x) == 1L) {
        if (length(instr.settings) == 1L && is.na(instr.settings)) {
          instr.settings <- list(integ.time = NA_real_,
                                 tot.time = NA_real_,
                                 num.scans = NA_integer_,
                                 rel.signal = NA_real_)
        }
        stopifnot("The argument passed to 'instr.settings' must be a list" =
                    is.list(instr.settings))
        if (!inherits(instr.settings, "instr_settings") &&
            !inherits(instr.settings[[1]], "instr_settings")) {
          class(instr.settings) <- c("instr_settings", class(instr.settings))
        }
      # when all settings are the same the attr is not expanded into a list
      } else if (!is.list(instr.settings) ||
                 (length(instr.settings) != getMultipleWl(x) &&
                  !any(c("integ.time", "tot.time", "num.scans") %in%
                       names(instr.settings)))) {
        warning("Length of 'instr.settings' list different to number of spectra: ",
                length(instr.settings), " != ", getMultipleWl(x))
      }
    }
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
#' Function to extract the \code{"instr.settings"} attribute from
#' \code{generic_spct} object or from a \code{summary_generic_spct}.
#'
#' @param x a \code{generic_spct} object or a \code{summary_generic_spct}
#'   object.
#'
#' @return an object of class \code{"instr_settings"} derived from \code{"list"}.
#'
#' @details If \code{x} is derived from \code{generic_spct} or from
#' \code{summary_generic_spct}, the value of attribute \code{"instr.settings"}
#' is returned (\code{NULL}, if missing). Otherwise \code{list()} is returned.
#'
#' @export
#'
#' @examples
#' settings <- getInstrSettings(white_led.cps_spct)
#' class(settings)
#' print(settings)
#' print(str(settings))
#'
#' @family measurement metadata functions
#'
getInstrSettings <- function(x) {
  if (is.generic_spct(x) || is.summary_generic_spct(x)) {
    if (getMultipleWl(x) == 1) {
      if (isValidInstrSettings(x)) {
        instr.settings <- attr(x, "instr.settings", exact = TRUE)
        # earlier bug created objects with extra empty fields with bad names
        bad.fields <- which(is.na(names(instr.settings)) |
                              names(instr.settings) %in% c("NA", ""))
        if (length(bad.fields)) {
          instr.settings <- instr.settings[-bad.fields]
        }
      } else {
        instr.settings <- list(integ.time = NA_real_,
                               tot.time = NA_real_,
                               num.scans = NA_integer_,
                               rel.signal = NA_real_)
      }
      # very old records lack class attribute
      if (!inherits(instr.settings, "instr_settings") &&
          !inherits(instr.settings[[1]], "instr_settings")) {
        class(instr.settings) <- c("instr_settings", class(instr.settings))
      }
      instr.settings
    } else {
      attr(x, "instr.settings", exact = TRUE)
    }
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
#' Trim the \code{"instr.settings"} attribute of an existing \code{generic_spct}
#' object or of a \code{summary_generic_spct} object, by discarding some fields.
#'
#' @param x a \code{generic_spct} object or a \code{summary_generic_spct}
#'   object.
#' @param fields a character vector with the names of the fields to keep, or if
#'   first member is \code{"-"}, the names of fields to delete; \code{"*"} as
#'   first member of the vector makes the function a no-op, leaving the spectrum
#'   object unaltered.
#'
#' @details This function alters \code{x} itself by reference and in addition
#'   returns \code{x} invisibly. If \code{x} is not a \code{generic_spct} object
#'   or a \code{summary_generic_spct} object, or if the \code{"instr.settings"}
#'   attribute is not present in \code{x}, \code{x} is not modified.
#'
#'   Attempts to remove or keep fields that are not present in the attribute are
#'   ignored silently. The value of fields in the attribute is never modified,
#'   fields are either kept unchanged or removed.
#'
#' @return \code{x}, possibly with the \code{"instr.settings"} attribute
#'   modified.
#'
#' @export
#'
#' @examples
#'
#' my.spct <- white_led.cps_spct
#' names(instr_settings(my.spct))
#' trimInstrSettings(my.spct, fields = c("-", "pix.selector")) # by reference!
#' names(instr_settings(my.spct))
#'
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
      # silently ignore missing fields, avoiding NAs
      if (!(length(instr.settings[[i]]) == 0 || all(is.na(instr.settings[[i]])))) {
        if (fields[1] == "-") {
          fields.tmp <- setdiff(names(instr.settings[[i]]), fields[-1])
        } else if (fields[1] == "=") {
          fields.tmp <- intersect(fields[-1], names(instr.settings[[i]]))
        } else {
          fields.tmp <- intersect(fields, names(instr.settings[[i]]))
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
#' Function to validate the \code{"instr.settings"} attribute of an existing
#' \code{generic_spct} or \code{summary_generic_spct} object.
#'
#' @param x a \code{generic_spct} object or a \code{summary_generic_spct}
#'   object.
#'
#' @return logical TRUE if at least the integration time is found in the
#'   metadata attribute. If \code{x} is not a \code{generic_spct} or
#'   a \code{summary_generic_spct} object, \code{NA} is returned.
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
#' Method to set by reference the \code{"what.measured"} attribute of an R
#' object.
#'
#' @param x an R object.
#' @param what.measured,value a list
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return x
#' @details This function alters \code{x} itself by reference and in addition
#'   returns \code{x} invisibly. If \code{x} does not belong to one of the
#'   supported classes, \code{x} is not modified.
#'
#' @export
#'
#' @examples
#' my.spct <- sun.spct
#' what_measured(my.spct)
#' what_measured(my.spct) <- "Sun"
#' what_measured(my.spct)
#' what_measured(my.spct) <- NULL
#' what_measured(my.spct)
#'
#' @family measurement metadata functions
#'
setWhatMeasured <- function(x, ...) {UseMethod("setWhatMeasured")}

#' @rdname setWhatMeasured
#'
#' @export
#'
`what_measured<-` <- function(x, value) {
  setWhatMeasured(x, what.measured = value)
}

#' @describeIn setWhatMeasured default
#' @export
setWhatMeasured.default <- function(x, what.measured, ...) {
  x
}

#' @describeIn setWhatMeasured generic_spct
#' @export
setWhatMeasured.generic_spct <- function(x, what.measured, ...) {
  name <- substitute(x)
  attr(x, "what.measured") <- what.measured
  if (is.name(name)) {
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  invisible(x)
}

#' @describeIn setWhatMeasured summary_generic_spct
#' @export
setWhatMeasured.summary_generic_spct <- setWhatMeasured.generic_spct

#' @describeIn setWhatMeasured data.frame
#' @export
setWhatMeasured.data.frame <- setWhatMeasured.generic_spct

#' @describeIn setWhatMeasured generic_mspct
#'
#' @export
setWhatMeasured.generic_mspct <- function(x,
                                          what.measured,
                                         ...) {
  msmsply(mspct = x, .fun = setWhatMeasured, ..., what.measured = what.measured)
}

#' Get the \code{"what.measured"} attribute
#'
#' Method to read the \code{"what.measured"} attribute of an R object.
#'
#' @param x an R object.
#' @param ... Allows use of additional arguments in methods for other classes.
#'
#' @return \code{character} vector An object containing a description of the
#'   data. If \code{x} does not belong to a supported class \code{NA} is
#'   returned.
#'
#' @export
#'
#' @family measurement metadata functions
#'
#' @inherit setWhatMeasured examples
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
getWhatMeasured.generic_spct <- function(x,
                                         ...,
                                         simplify = FALSE) {
  z <- attr(x, "what.measured", exact = TRUE)
  if (is.null(z) ||
      (is.atomic(z) && all(is.na(z)))) {
    # need to handle objects created with old versions
    z <- NA_character_
  } else {
    if (simplify) {
      if (is.list(z)) {
        if (length(z) == 1L || sum(!duplicated(z) == 1L)) {
          z <- z[[1]]
        } else {
          z <- unlist(z, use.names = TRUE)
        }
      }
    }
  }
  z
}

#' @describeIn getWhatMeasured summary_generic_spct
#' @export
getWhatMeasured.summary_generic_spct <- getWhatMeasured.generic_spct

#' @describeIn getWhatMeasured data.frame
#' @export
#'
getWhatMeasured.data.frame <- getWhatMeasured.generic_spct

#' @describeIn getWhatMeasured generic_mspct
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @param simplify logical If all members share the same attribute value return
#'   one copy instead of a \code{data.frame}.
#' @note The method for collections of spectra returns the
#'   a \code{data.frame} with a column of character strings.
#' @export
#'
getWhatMeasured.generic_mspct <- function(x,
                                          ...,
                                          idx = "spct.idx",
                                          simplify = FALSE) {
  z <- msdply(mspct = x,
              .fun = getWhatMeasured,
              ...,
              idx = idx,
              col.names = "what.measured")
  if (simplify) {
    zz <- unique(z[["what.measured"]])
    if (length(zz) <= 1) {
      z <- zz
    } else {
      z <- paste(z[["what.measured"]], "\n", collapse = "", sep = "")
    }
  }
  z
}

# copy attributes to data.frame -------------------------------------------

#' Copy attributes from members of a \code{generic_mspct}
#'
#' Copy metadata attributes from members of a \code{generic_mspct} object into
#' a \code{data.frame} or a \code{tibble}.
#'
#' @param mspct generic_mspct or generic_spct Any collection of spectra or one
#'   or more spectra in long form.
#' @param tb tibble or \code{data.frame} to which to add the data (optional).
#' @param col.names named \code{character} vector Name(s) of metadata attributes
#'   to copy. If named, the names provide the name for the columns.
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @param unnest logical Flag controlling if metadata attributes that are lists
#'   of values should be returned in a list column or in separate columns.
#'
#' @return A \code{data.frame} or a \code{tibble} With the metadata attributes
#'   in separate new variables.
#'
#' @details Each attribute is by default copied to a column in a \code{tibble}
#'   or a \code{data.frame}. If the argument for \code{tb} is \code{NULL}, as by
#'   default, a new \code{tibble} will be created. If an existing
#'   \code{data.frame} or \code{tibble} is passed as argument, new columns are
#'   added to it. However, the number of rows in the argument passed to
#'   \code{tb} must match the number of spectra in the argument passed to
#'   \code{mspct}. Only in the case of methods \code{add_attr2tb()} and
#'   \code{spct_metadata()} if the argument to \code{col.names} is a named
#'   vector, the names of members are used as names for the columns created.
#'   This permits setting any valid name for the new columns. If the members of
#'   the vector passed to \code{col.names} have no names, then the value is
#'   interpreted as the name of the attributes to add, and also used as name for
#'   the new column.
#'
#'   Valid values accepted as argument to \code{col.names} are \code{NULL}, or a
#'   vector containing one or more of the following \code{character} strings:
#'   \code{"lon"}, \code{"lat"}, \code{"address"}, \code{"geocode"},
#'   \code{"where.measured"}, \code{"when.measured"}, \code{"what.measured"},
#'   \code{"how.measured"}, \code{"comment"}, \code{"normalised"},
#'   \code{"normalized"}, \code{"scaled"}, \code{"bswf.used"},
#'   \code{"instr.desc"}, \code{"instr.sn"}, \code{solute.properties},
#'   \code{"filter.properties"}, \code{"Tfr.type"}, \code{"Rfr.type"},
#'   \code{"time.unit"}, \code{bswf.used}, \code{multiple.wl}. Invalid character
#'   values are ignored with a warning.
#'
#' @note The order of the first two arguments is reversed in
#'   \code{add_attr2tb()}, \code{when_measured2tb()}, \code{what_measured2tb()},
#'   etc., compared to attribute query functions, such as \code{spct_metadata},
#'   \code{when_measured()}, \code{what_measured()}, \code{how_measured()}, etc.
#'   This is to allow the use of \code{add_attr2tb()} and related functions in
#'   'pipes' to add metadata to summaries computed at earlier steps in the pipe.
#'
#' @family measurement metadata functions
#'
#' @examples
#' # Add attributes to irradiance
#' ## from collection of spectra
#' e_irrad(sun_evening.mspct) |>
#'   add_attr2tb(sun_evening.mspct,
#'               c(when.measured = "time"))
#'
#' ## from spectra in long form
#' e_irrad(sun_evening.spct) |>
#'   add_attr2tb(sun_evening.spct,
#'               c(when.measured = "time"))
#'
#' # Add attributes to transmittance
#' ## from collection of spectra
#' transmittance(two_filters.mspct) |>
#'   add_attr2tb(two_filters.mspct, col.names = "what.measured")
#'
#' transmittance(two_filters.mspct) |>
#'   add_attr2tb(two_filters.mspct,
#'               col.names = c("filter.properties", "what.measured"),
#'               unnest = TRUE)
#'
#' # Create a new data frame
#' add_attr2tb(mspct = two_filters.mspct,
#'             idx = "filter",
#'             col.names = c("filter.properties", "what.measured"),
#'             unnest = TRUE)
#'
#' @export
#'
add_attr2tb <- function(tb = NULL,
                        mspct,
                        col.names = NULL,
                        idx = "spct.idx",
                        unnest = FALSE) {
  # ensure we operate on a collection of spectra
  if (is.any_spct(mspct)) {
    mspct <- subset2mspct(mspct)
  }
  # we accept NULL, list(), generic_mspct() as input
  if (length(mspct) == 0L) {
    return(data.frame())
  }
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
    # setdiff removes names from the vector to avoid duplicated columns!
    col.names <- col.names[!col.names %in% c("lat", "lon")]
  }
  # We walk the list of attributes adding columns
  tb.cols <- names(tb)
  for (a in names(col.names)) {
    tb <-
      switch(a,
             multiple.wl =
               multiple_wl2tb(mspct = mspct,
                              tb = tb,
                              col.names = col.names["multiple.wl"],
                              idx = idx),
             lon =
               lon2tb(mspct = mspct,
                      tb = tb,
                      col.names = col.names["lon"],
                      idx = idx),
             lat =
               lat2tb(mspct = mspct,
                      tb = tb,
                      col.names = col.names["lat"],
                      idx = idx),
             address =
               address2tb(mspct = mspct,
                          tb = tb,
                          col.names = col.names["address"],
                          idx = idx),
             geocode =
               geocode2tb(mspct = mspct,
                          tb = tb,
                          col.names = col.names["geocode"],
                          idx = idx),
             where.measured =
               geocode2tb(mspct = mspct,
                          tb = tb,
                          col.names = col.names["where.measured"],
                          idx = idx),
             when.measured =
               when_measured2tb(mspct = mspct,
                                tb = tb,
                                col.names = col.names["when.measured"],
                                idx = idx),
             what.measured =
               what_measured2tb(mspct = mspct,
                                tb = tb,
                                col.names = col.names["what.measured"],
                                idx = idx),
             how.measured =
               how_measured2tb(mspct = mspct,
                               tb = tb,
                               col.names = col.names["how.measured"],
                               idx = idx),
             comment =
               comment2tb(mspct = mspct,
                          tb = tb,
                          col.names = col.names["comment"],
                          idx = idx),
             normalized =
               normalized2tb(mspct = mspct,
                             tb = tb,
                             col.names = col.names["normalized"],
                             idx = idx),
             normalised =
               normalized2tb(mspct = mspct,
                             tb = tb,
                             col.names = col.names["normalised"],
                             idx = idx),
             scaled =
               scaled2tb(mspct = mspct,
                         tb = tb,
                         col.names = col.names["scaled"],
                         idx = idx),
             instr.desc =
               instr_desc2tb(mspct = mspct,
                             tb = tb,
                             col.names = col.names["instr.desc"],
                             idx = idx),
             instr.sn =
               instr_desc2tb(mspct = mspct,
                             tb = tb,
                             fields = "spectrometer.sn",
                             col.names = col.names["instr.sn"],
                             idx = idx),
             instr.settings =
               instr_settings2tb(mspct = mspct,
                                 tb = tb,
                                 col.names = col.names["instr.settings"],
                                 idx = idx),
             filter.properties =
               filter_properties2tb(mspct = mspct,
                                    tb = tb,
                                    col.names = col.names["filter.properties"],
                                    idx = idx),
             solute.properties =
               solute_properties2tb(mspct = mspct,
                                    tb = tb,
                                    col.names = col.names["solute.properties"],
                                    idx = idx),
             Tfr.type =
               Tfr_type2tb(mspct = mspct,
                           tb = tb,
                           col.names = col.names["Tfr.type"],
                           idx = idx),
             Rfr.type =
               Rfr_type2tb(mspct = mspct,
                           tb = tb,
                           col.names = col.names["Rfr.type"],
                           idx = idx),
             time.unit =
               time_unit2tb(mspct = mspct,
                            tb = tb,
                            col.names = col.names["time.unit"],
                            idx = idx),
             bswf.used =
               BSWF_used2tb(mspct = mspct,
                            tb = tb,
                            col.names = col.names["bswf.used"],
                            idx = idx),
             {warning("Skipping unknown metada name: ", a);
               tb})
  }
  if (unnest) {
    list.cols <- colnames(tb)[sapply(tb, is.list)]
    # do not expand pre-existing list columns
    list.cols <- setdiff(list.cols, tb.cols)
    # expand metadata fields into columns
    for (col in list.cols) {
      # avoid duplicate names
      tb[[col]] <- lapply(tb[[col]], function(x) x[setdiff(names(x), colnames(tb))])
      # handles lists of lists or lists of data frames
      tb <- tidyr::unnest_wider(tb,
                                tidyr::all_of(col),
                                names_repair = "check_unique")
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
  where.ls <- getWhereMeasured(mspct, idx = idx, .bind.geocodes = FALSE)
  where.tb <- tibble::tibble(factor(names(where.ls)), where.ls)
  names(where.tb) <- c(idx, col.names)
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
  lonlat.tb <- getWhereMeasured(mspct,
                                idx = idx,
                                .bind.geocodes = TRUE)[c(idx, "lon", "lat")]
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
  lon.tb <- getWhereMeasured(mspct,
                             idx = idx,
                             .bind.geocodes = TRUE)[c(idx, "lon")]
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
  lat.tb <- getWhereMeasured(mspct,
                             idx = idx,
                             .bind.geocodes = TRUE)[c(idx, "lat")]
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
  address.tb <- getWhereMeasured(mspct,
                                 idx = idx,
                                 .bind.geocodes = TRUE)[c(idx, "address")]
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
#' @param fields character vector or logical Names of fields to extract from
#'   each descriptor record.
#'
#' @export
#'
instr_desc2tb <- function(mspct,
                          tb = NULL,
                          col.names = "instr.desc",
                          fields = TRUE,
                          idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  # allow partial extraction
  getInstrDescFields <- function(x, fields) {
    z <- getInstrDesc(x)
    selector <- ifelse(is.character(fields),
                       intersect(names(z), fields),
                       fields)
    z[selector]
  }
  l <- mslply(mspct = mspct, .fun = getInstrDescFields, fields = fields)
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
                              fields = TRUE,
                              idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  # allow partial extraction
  getInstrSettingsFields <- function(x, fields) {
    z <- getInstrSettings(x)
    selector <- ifelse(is.character(fields),
                       intersect(names(z), fields),
                       fields)
    z[selector]
  }
  l <- mslply(mspct = mspct, .fun = getInstrSettingsFields, fields = fields)
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
                         col.names = "bswf.used",
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

#' @rdname add_attr2tb
#'
#' @export
#'
multiple_wl2tb <- function(mspct,
                           tb = NULL,
                           col.names = "multiple.wl",
                           idx = "spct.idx") {
  stopifnot(length(col.names) == 1L)
  # method not implemented yet for collections
  multiple_wl.tb <- msdply(mspct = mspct, .fun = getMultipleWl)
  names(multiple_wl.tb)[2L] <- col.names
  if (is.null(tb)) {
    multiple_wl.tb
  } else {
    dplyr::full_join(tb, multiple_wl.tb, by = idx)
  }
}

# extract all metadata ----------------------------------------------------

#' Access metadata
#'
#' Return metadata attributes from a single spectrum or a collection of spectra
#' as a \code{data.frame}. A wrapper on \code{add_attr2tb} providing an
#' alternative order of formal parameters and constrained functionality.
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
#' @inherit add_attr2tb details return note
#'
#' @seealso \code{\link{add_attr2tb}} for more details.
#'
#' @family measurement metadata functions
#'
#' @examples
#' # collection of spectra
#' spct_metadata(sun_evening.mspct)
#'
#' spct_metadata(sun_evening.mspct, na.rm = FALSE)
#'
#' spct_metadata(sun_evening.mspct,
#'               col.names = "geocode",
#'               unnest = FALSE)
#'
#' spct_metadata(sun_evening.mspct,
#'               col.names = c(when.measured = "time", "what.measured"))
#'
#' # multiple spectra in long form
#' spct_metadata(sun_evening.spct,
#'               col.names = c("geocode", "when.measured"))
#'
#' # single spectrum
#' spct_metadata(sun.spct,
#'               col.names = c("geocode", "when.measured"))
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
  # ensure we operate on a collection of spectra
  if (is.any_spct(x)) {
    x <- subset2mspct(x)
  }
  if (!is.generic_mspct(x) || length(x) == 0L) {
    return(data.frame())
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
