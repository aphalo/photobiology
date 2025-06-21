#' Column or variable labels
#'
#' Create a named list of character strings describing the variables contained
#' in a spectrum object.
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
#'   do not add the labels, only supply the character strings. Variable labels
#'   are implemented in packages 'labelled' by setting the \code{label}
#'   attribute in each variable (= column) of a data frame or tibble. This is
#'   compatible with the approach used by package 'haven'.
#'
#' @note These methods are still under development and the text of the labels
#'   may change. Not all classes derived from \code{generic_spct} are yet
#'   supported.
#'
#' @return A named list of character strings with one member for each recognized
#'   column in \code{x}. This list can be used to set variable labels with
#'   methods from package 'labelled'. However, package 'photobiology' does not
#'   natively support variable labels stored in attribute \code{label}.
#'
#' @export
#'
#' @examples
#' make_var_labels(sun.spct)
#' # str() prints more compactly than print()
#' str(make_var_labels(sun.spct))
#' str(make_var_labels(normalize(sun.spct)))
#' str(make_var_labels(fscale(sun.spct)))
#'
#' str(make_var_labels(sun_daily.spct))
#'
#' str(make_var_labels(polyester.spct))
#' str(make_var_labels(normalize(polyester.spct)))
#' str(make_var_labels(fscale(polyester.spct)))
#'
#' str(make_var_labels(white_led.cps_spct))
#' str(make_var_labels(white_led.raw_spct))
#'
make_var_labels <- function(x, ...) {
  UseMethod("make_var_labels")
}

#' @describeIn make_var_labels
#'
#' @export
#'
make_var_labels.default <-
  function(x, ...) {
    list()
  }

#' @describeIn make_var_labels
#'
#' @export
#'
make_var_labels.source_spct <- function(x, ...) {
  time.unit <- getTimeUnit(x, force.duration = FALSE)
  if (lubridate::is.duration(time.unit)) {
    if (time.unit == lubridate::seconds(1)) {
      time.unit <- "second"
    } else if (time.unit == lubridate::hours(1)) {
      time.unit <- "hour"
    } else if (time.unit == lubridate::days(1)) {
      time.unit <- "day"
    } else {
      labels <-
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy exposure [J m-2 nm-1]",
             s.q.irrad = "Spectral photon exposure [mol m-2 nm-1]")[colnames(x)]
    }
  }
  if (is.character(time.unit)) {
    if (time.unit == "second") {
      labels <-
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy irradiance [W m-2 nm-1]",
             s.q.irrad = "Spectral photon irradiance [mol s-1 m-2 nm-1]")
    } else if (time.unit == "hour") {
      labels <-
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy irradiance [J h-1 m-2 nm-1]",
             s.q.irrad = "Spectral photon irradiance [mol h-1 m-2 nm-1]")
    } else if (time.unit == "day") {
      labels <-
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy exposure [J d-1 m-2 nm-1]",
             s.q.irrad = "Spectral photon exposure [mol d-1 m-2 nm-1]")
    } else if (time.unit == "exposure") {
      labels <-
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy exposure [J m-2 nm-1]",
             s.q.irrad = "Spectral photon exposure [mol m-2 nm-1]")
    } else if (time.unit == "none") {
      labels <-
        list(w.length = "Wavelength [nm]",
             s.e.irrad = "Spectral energy exposure",
             s.q.irrad = "Spectral photon exposure")[colnames(x)]
    }
  }
  # lapply preserves the list and member names
  sub.pattern <-
    "J m-2 nm-1|mol m-2 nm-1|W m-2 nm-1|J [shd]-1 m-2 nm-1|mol [shd]-1 m-2 nm-1"
  if (is_normalized(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "normalized")
  } else if (is_scaled(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "scaled")
  }

  labels[intersect(colnames(x), names(labels))]

}

#' @describeIn make_var_labels
#'
#' @export
#'
make_var_labels.response_spct <- function(x, ...) {
  time.unit <- getTimeUnit(x, force.duration = FALSE)
  if (lubridate::is.duration(time.unit)) {
    if (time.unit == lubridate::seconds(1)) {
      time.unit <- "second"
    } else if (time.unit == lubridate::hours(1)) {
      time.unit <- "hour"
    } else if (time.unit == lubridate::days(1)) {
      time.unit <- "day"
    } else {
      labels <-
        list(w.length = "Wavelength [nm]",
             s.e.response = "Spectral energy response [J-1 m2 nm]",
             s.q.response = "Spectral photon response [mol-1 m nm]")
    }
  }
  if (is.character(time.unit)) {
    if (time.unit == "second") {
      labels <-
        list(w.length = "Wavelength [nm]",
             s.e.response = "Spectral energy response [W-1 m2 nm]",
             s.q.response = "Spectral photon response [mol-1 s m2 nm]")
    } else if (time.unit == "hour") {
      labels <-
        list(w.length = "Wavelength [nm]",
             s.e.response = "Spectral energy response [J-1 h m2 nm]",
             s.q.response = "Spectral photon response [mol-1 h m2 nm]")
    } else if (time.unit == "day") {
      labels <-
        list(w.length = "Wavelength [nm]",
             s.e.response = "Spectral energy response [J-1 d m2 nm]",
             s.q.response = "Spectral photon response [mol-1 d m2 nm]")
    } else if (time.unit == "exposure") {
      labels <-
        list(w.length = "Wavelength [nm]",
             s.e.response = "Spectral energy response [J-1 m2 nm]",
             s.q.response = "Spectral photon response [mol-1 m2 nm]")
    } else if (time.unit == "none") {
      labels <-
        list(w.length = "Wavelength [nm]",
             s.e.response = "Spectral energy response",
             s.q.response = "Spectral photon response")
    }
  }
  sub.pattern <- "J-1 m2 nm|mol-1 m nm|J-1 [shd] m2 nm|mol-1 [shd] m2 nm"
  if (is_normalized(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "normalized")
  } else if (is_scaled(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "scaled")
  }

  labels[intersect(colnames(x), names(labels))]

}

#' @describeIn make_var_labels
#'
#' @export
#'
make_var_labels.filter_spct <- function(x, ...) {

  Tfr.type <- getTfrType(x)
  Tfr.label <- c(total = "Total spectral transmittance [/1]",
                 internal = "Internal spectral transmittance [/1]",
                 unknown = "Spectral transmittance [/1]")
  labels <-
    list(w.length = "Wavelength [nm]",
         Tfr = Tfr.label[Tfr.type],
         Afr = "Spectral absorptance [/1]",
         A = "Spectral absorbance log10 based [a.u.]")
  sub.pattern <- "/1|a\\.u\\."
  if (is_normalized(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "normalized")
  } else if (is_scaled(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "scaled")
  }

  labels[intersect(colnames(x), names(labels))]

}

#' @describeIn make_var_labels
#'
#' @export
#'
make_var_labels.reflector_spct <- function(x, ...) {

  Rfr.type <- getRfrType(x)
  Rfr.label <- c(total = "Total spectral reflectance [/1]",
                 specular = "Specular spectral reflectance [/1]",
                 unknown = "Spectral reflectance [/1]")
  labels <-
    list(w.length = "Wavelength [nm]",
         Rfr = Rfr.label[Rfr.type])
  sub.pattern <- "/1"
  if (is_normalized(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "normalized")
  } else if (is_scaled(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "scaled")
  }

  labels[intersect(colnames(x), names(labels))]

}

#' @describeIn make_var_labels
#'
#' @export
#'
make_var_labels.object_spct <- function(x, ...) {
  Tfr.type <- getTfrType(x)
  Tfr.label <- c(total = "Total spectral transmittance [/1]",
                 internal = "Internal spectral transmittance [/1]",
                 unknown = "Spectral transmittance [/1]")
  Rfr.type <- getRfrType(x)
  Rfr.label <- c(total = "Total spectral reflectance [/1]",
                 specular = "Specular spectral reflectance [/1]")
  labels <-
    list(w.length = "Wavelength [nm]",
         Tfr = Tfr.label[Tfr.type],
         Rfr = Rfr.label[Rfr.type],
         Afr = "Spectral absorptance [/1]",
         A = "Spectral absorbance log10 based [a.u.]")
  # scaling and normalization not supported by class object_spct

  labels[intersect(colnames(x), names(labels))]

}

#' @describeIn make_var_labels
#'
#' @export
#'
make_var_labels.solute_spct <- function(x, ...) {
  # To follow IUPAC "attenuation" could be used in all cases
  K.type <- getKType(x)
  labels <-
    list(w.length = "Wavelength [nm]",
         K.mole = paste("Molar",  K.type, "coefficient [m2 mol-1]"),
         K.mass = paste("Mass", K.type, "coefficient [m2 g-1]"))
  sub.pattern <- "m2 mol-1|m2 g-1"
  if (is_normalized(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "normalized")
  } else if (is_scaled(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "scaled")
  }

  labels[intersect(colnames(x), names(labels))]

}

#' @describeIn make_var_labels
#'
#' @export
#'
make_var_labels.chroma_spct <- function(x, ...) {
  labels <-
    list(w.length = "Wavelength [nm]",
         x = "Numeric colour coordinate X",
         y = "Numeric colour coordinates Y",
         z = "Numeric colour coordinates Z")
  # scaling and normalization not supported by class chroma_spct

  labels[intersect(colnames(x), names(labels))]

}

#' @describeIn make_var_labels
#'
#' @export
#'
make_var_labels.calibration_spct <- function(x, ...) {
  labels <-
    list(w.length = "Wavelength [nm]",
         irrad.mult = "Irradiance multipliers [J m-2 nm-1 n-1]")
  # scaling and normalization not supported by class calibration_spct

  labels[intersect(colnames(x), names(labels))]

}

#' @describeIn make_var_labels
#'
#' @export
#'
make_var_labels.raw_spct <- function(x, ...) {
  column.names <- colnames(x)
  count.cols <- grepl("^counts", column.names)
  labels <- c("Wavelength [nm]",
              rep("Raw detector counts [number]", sum(count.cols)))
  names(labels) <- c("w.length", column.names[count.cols])
  sub.pattern <- "number"
  if (is_normalized(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "normalized")
  } else if (is_scaled(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "scaled")
  }

  as.list(labels)

}

#' @describeIn make_var_labels
#'
#' @export
#'
make_var_labels.cps_spct <- function(x, ...) {
  column.names <- colnames(x)
  count.cols <- grepl("^cps", column.names)
  labels <- c("Wavelength [nm]",
              rep("Detector counts [number s-1]", sum(count.cols)))
  names(labels) <- c("w.length", column.names[count.cols])
  sub.pattern <- "number  s-1"
  if (is_normalized(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "normalized")
  } else if (is_scaled(x)) {
    labels <-
      lapply(labels, gsub, pattern = sub.pattern, replacement = "scaled")
  }

  as.list(labels)

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
