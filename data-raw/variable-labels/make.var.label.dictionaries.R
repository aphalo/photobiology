#' Column or variable labels
#'
#' Create a named list of character strings suitable to be used to set variable
#' labels with methods from package 'labelled' or similar.
#'
#' @param x An object of a class derived from \code{generic_spct}.
#' @param ... Currently ignored.
#'
#' @details Objects of classes derived from \code{generic_spct} are used to
#'   store different types of spectral data. The data stored in some of the
#'   classes needs to be interpreted differently depending on how they were
#'   measured or are expressed and this information is stored in attributes of
#'   the objects. In other cases, even if consistent across different objects,
#'   the units of expression may not be obvious to users. The names of the
#'   variables are concise, thus using variable labels makes it possible to
#'   make these features visible when exploring the data. The methods provided
#'   do not add the labels, only supply the character strings.
#'
#' @return A named list of character strings with one member for each recognized
#'   column in \code{x}.
#'
#' @export
#'
#' @examples
#' col_labels(sun.spct)
#' col_labels(sun_daily.spct)
#'
col_labels <- function(x, ...) {
    UseMethod("col_labels")
}

#' @describeIn col_labels
#'
#' @export
#'
col_labels.default <-
  function(x, ...) {
    list()
  }

#' @describeIn col_labels
#'
#' @export
#'
col_labels.source_spct <- function(x, ...) {
  time.unit <- getTimeUnit(x, force.duration = FALSE)
  if (lubridate::is.duration(time.unit)) {
    if (time.unit == lubridate::seconds(1)) {
      time.unit <- "second"
    } else if (time.unit == lubridate::hours(1)) {
      time.unit <- "hour"
    } else if (time.unit == lubridate::days(1)) {
      time.unit <- "day"
    } else {
      return(
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy exposure [J m-2 nm-1]",
             s.q.irrad = "Spectral photon exposure [mol m-2 nm-1]")[colnames(x)]
      )
    }
  }
  if (is.character(time.unit)) {
    if (time.unit == "second") {
      return(
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy irradiance [W m-2 nm-1]",
             s.q.irrad = "Spectral photon irradiance [mol s-1 m-2 nm-1]")[colnames(x)]
      )
    } else if (time.unit == "hour") {
      return(
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy irradiance [J h-1 m-2 nm-1]",
             s.q.irrad = "Spectral photon irradiance [mol h-1 m-2 nm-1]")[colnames(x)]
      )
    } else if (time.unit == "day") {
      return(list(w.length = "Wavelength [nm]",
                  s.e.irrad = "Spectral energy exposure [J d-1 m-2 nm-1]",
                  s.q.irrad = "Spectral photon exposure [mol d-1 m-2 nm-1]")[colnames(x)]
      )
    } else if (time.unit == "exposure") {
      return(
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy exposure [J m-2 nm-1]",
             s.q.irrad = "Spectral photon exposure [mol m-2 nm-1]")[colnames(x)]
      )
    } else if (time.unit == "none") {
      return(
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy exposure",
             s.q.irrad = "Spectral photon exposure")[colnames(x)]
      )
    }
  }
}

#' @describeIn col_labels
#'
#' @export
#'
col_labels.response_spct <- function(x, ...) {
  time.unit <- getTimeUnit(x, force.duration = FALSE)
  if (lubridate::is.duration(time.unit)) {
    if (time.unit == lubridate::seconds(1)) {
      time.unit <- "second"
    } else if (time.unit == lubridate::hours(1)) {
      time.unit <- "hour"
    } else if (time.unit == lubridate::days(1)) {
      time.unit <- "day"
    } else {
      return(
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy response [J-1 m2 nm]",
             s.q.irrad = "Spectral photon response [mol-1 m nm1]")[colnames(x)]
      )
    }
  }
  if (is.character(time.unit)) {
    if (time.unit == "second") {
      return(
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy response [W-1 m2 nm]",
             s.q.irrad = "Spectral photon response [mol-1 s m2 nm]")[colnames(x)]
      )
    } else if (time.unit == "hour") {
      return(
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy response [J-1 h m2 nm]",
             s.q.irrad = "Spectral photon response [mol-1 h m2 nm]")[colnames(x)]
      )
    } else if (time.unit == "day") {
      return(list(w.length = "Wavelength [nm]",
                  s.e.irrad = "Spectral energy response [J-1 d m2 nm]",
                  s.q.irrad = "Spectral photon response [mol-1 d m2 nm]")[colnames(x)]
      )
    } else if (time.unit == "exposure") {
      return(
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy response [J-1 m2 nm]",
             s.q.irrad = "Spectral photon response [mol-1 m2 nm]")[colnames(x)]
      )
    } else if (time.unit == "none") {
      return(
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy response",
             s.q.irrad = "Spectral photon response")[colnames(x)]
      )
    }
  }
}

#' @describeIn col_labels
#'
#' @export
#'
col_labels.filter_spct <- function(x, ...) {
  Tfr.type <- getTfrType(x)
  if (Tfr.type == "total") {
    return(
      list(w.length = "Wavelength [nm]",
           Tfr = "Total spectral transmittance [/1]",
           Afr = "Spectral absorptance [/1]",
           A = "Spectral absorbance log10 based [a.u.]")[colnames(x)]
    )
  } else if (Tfr.type == "internal") {
    return(
      list(w.length = "Wavelength [nm]",
           Tfr = "Internal spectral transmittance [/1]",
           Afr = "Spectral absorptance [/1]",
           A = "Spectral absorbance log10 based [a.u.]")[colnames(x)]
    )
  } else if (Tfr.type == "unknown") {
    return(list(w.length = "Wavelength [nm]",
                Tfr = "Spectral transmittance [/1]",
                Afr = "Spectral absorptance [/1]",
                A = "Spectral absorbance log10 based [a.u.]")[colnames(x)]
    )
  } else {
    return(
      list(w.length = "Wavelength [nm]")[colnames(x)]
    )
  }
}

#' @describeIn col_labels
#'
#' @export
#'
col_labels.reflector_spct <- function(x, ...) {
  Rfr.type <- getRfrType(x)
  if (Rfr.type == "total") {
    return(
      list(w.length = "Wavelength [nm]",
           Rfr = "Total spectral reflectance [/1]")[colnames(x)]
    )
  } else if (Rfr.type == "specular") {
    return(
      list(w.length = "Wavelength [nm]",
           Rfr = "Specular spectral reflectance [/1]")[colnames(x)]
    )
  } else if (Rfr.type == "unknown") {
    return(list(w.length = "Wavelength [nm]",
                Rfr = "Spectral reflectance [/1]")[colnames(x)]
    )
  } else {
    return(
      list(w.length = "Wavelength [nm]")[colnames(x)]
    )
  }
}

#' @describeIn col_labels
#'
#' @export
#'
col_labels.object_spct <- function(x, ...) {
  Tfr.type <- getTfrType(x)
  Rfr.type <- getRfrType(x)
  Rfr.label <- c(total = "Total spectral reflectance [/1]",
                 specular = "Specular spectral reflectance [/1]")
  if (Tfr.type == "total") {
    return(
      list(w.length = "Wavelength [nm]",
           Tfr = "Total spectral transmittance [/1]",
           Rfr = Rfr.label[Rfr.type],
           Afr = "Spectral absorptance [/1]",
           A = "Spectral absorbance log10 based [a.u.]")[colnames(x)]
    )
  } else if (Tfr.type == "internal") {
    return(
      list(w.length = "Wavelength [nm]",
           Tfr = "Internal spectral transmittance [/1]",
           Rfr = Rfr.label[Rfr.type],
           Afr = "Spectral absorptance [/1]",
           A = "Spectral absorbance log10 based [a.u.]")[colnames(x)]
    )
  } else if (Tfr.type == "unknown") {
    return(list(w.length = "Wavelength [nm]",
                Tfr = "Spectral transmittance [/1]",
                Rfr = Rfr.label[Rfr.type],
                Afr = "Spectral absorptance [/1]",
                A = "Spectral absorbance log10 based [a.u.]")[colnames(x)]
    )
  } else {
    return(
      list(w.length = "Wavelength [nm]")[colnames(x)]
    )
  }
}

# summary.sec.variable.labels <-
#   list(e.irrad = "Energy irradiance [W m-2]",
#        q.irrad = "Photon irradiance [umol m-2 s-1]",
#        Tfr = "Transmittance [/1]",
#        Rfr = "Reflectance [/1]",
#        Afr = "Absorptance [/1]",
#        A = "Absorbance log10 based [a.u.]",
#        e.response = "Energy responsivity",
#        q.response = "Photon responsivity")
#
# summary.day.variable.labels <-
#   list(e.irrad = "Energy exposure [W m-2]",
#        q.irrad = "Photon exposure [umol m-2 s-1]",
#        Tfr = "Transmittance [/1]",
#        Rfr = "Reflectance [/1]",
#        Afr = "Absorptance [/1]",
#        A = "Absorbance log10 based [a.u.]",
#        e.response = "Energy responsivity",
#        q.response = "Photon responsivity")
