# color_of ------------------------------------------------------------------

#' Color of an object
#'
#' Equivalent RGB color of an object such as a spectrum, wavelength or waveband.
#'
#' @param x an R object.
#' @param ... ignored (possibly used by derived methods).
#'
#' @export color_of
#'
#' @return A color definition in hexadecimal format as a \code{character} string
#'   of 7 characters, "#" followed by the red, blue, and green values in
#'   hexadecimal (scaled to 0 ... 255). In the case of the specialization for
#'   \code{list}, a list of such definitions is returned. In the case of a
#'   collection of spectra, a \code{data.frame} with one column with such
#'   definitions and by default an additional column with names of the spectra
#'   as index. In case of missing input the returned value is \code{NA}.
#'
#' @examples
#' wavelengths <- c(300, 420, 500, 600, NA) # nanometres
#' color_of(wavelengths)
#' color_of(waveband(c(300,400)))
#' color_of(list(blue = waveband(c(400,480)), red = waveband(c(600,700))))
#' color_of(numeric())
#' color_of(NA_real_)
#'
#' color_of(sun.spct)
#'
color_of <- function(x, ...) UseMethod("color_of")

#' @describeIn color_of Default method (returns always "black").
#'
#' @export
#'
color_of.default <- function(x, ...) {
  if (length(x) == 0) {
    character()
  } else if (is.vector(x)) {
    warning("'color_of()' not implemented for class \"", class(x)[1], "\".")
    rep(NA_character_, length(x))
  } else if (is.any_spct(x)) {
    warning("To obtain the colour of a filter or reflector, first multiply",
            " by spectral irradiance under which it is observed.")
    NA_character_
  }
}

#' @describeIn color_of Method that returns Color definitions corresponding to
#'   numeric values representing a wavelengths in nm.
#'
#' @param type,chroma.type character telling whether "CMF", "CC", or "both" should be returned
#'   for human vision, or an object of class \code{chroma_spct} for any other
#'   trichromic visual system.
#'
#' @export
#'
color_of.numeric <- function(x,
                             type = "CMF",
                             chroma.type = type,
                             ...) {
  if (length(x) == 0) {
    return(character())
  }
  stopifnot(all(ifelse(is.na(x), Inf, x) > 0)) # w.length > 0 or NA
  if (is.character(chroma.type)) {
    if (chroma.type == "CMF") {
      color.out <-
        w_length2rgb(x, sens = photobiology::ciexyzCMF2.spct,
                     color.name = NULL)
    } else if (chroma.type == "CC") {
      color.out <-
        w_length2rgb(x, sens = photobiology::ciexyzCC2.spct, color.name = NULL)
    } else {
      warning("Chroma type = ", chroma.type, " not implemented for 'numeric'.")
      color.out <- rep(NA_character_, length(x))
    }
    if (!is.null(names(x))) {
      names(color.out) <- paste(names(x), chroma.type, sep = ".")
    } else {
      names(color.out) <- paste(names(color.out), chroma.type, sep = ".")
    }
  } else if (is.chroma_spct(chroma.type)) {
    #    warning("Using a CC or CMF definition, RGB may not denote red, green, blue!")
    color.out <-
      w_length2rgb(x, sens = chroma.type,
                   color.name = NULL)
    if (!is.null(names(x))) {
      names(color.out) <- paste(names(x), "chroma", sep = ".")
    } else {
      names(color.out) <- paste(names(color.out), "chroma", sep = ".")
    }
  } else {
    warning("Chroma type of class \"", class(chroma.type)[1], "\" not supported.")
    color.out <- rep(NA_character_, length(x))
  }
  color.out
}

#' @describeIn color_of Method that returns Color of elements in a list.
#'
#' @param short.names logical indicating whether to use short or long names for
#'   wavebands
#'
#' @note When \code{x} is a list but not a waveband, if a method  \code{color_of}
#'   is not available for the class of each element of the list, then
#'   \code{color_of.default} will be called.
#' @export
#'
color_of.list <- function(x,
                          short.names = TRUE,
                          type = "CMF",
                          chroma.type = type,
                          ...) {
  color.out <- character(0)
  for (xi in x) {
    color.out <- c(color.out,
                   color_of(xi,
                            short.names = short.names,
                            chroma.type = chroma.type,
                            ...))
  }
  if (!is.null(names(x))) {
    names(color.out) <- paste(names(x), chroma.type, sep = ".")
  }
  return(color.out)
}

#' @describeIn color_of Color at midpoint of a \code{\link{waveband}} object.
#'
#' @export
#'
color_of.waveband <-  function(x,
                               short.names = TRUE,
                               type = "CMF",
                               chroma.type = type,
                               ...) {
  idx <- ifelse(!short.names, "name", "label")
  name <- labels(x)[[idx]]
  if (is.character(chroma.type)) {
    if (chroma.type == "both") {
      color <- list(CMF = w_length_range2rgb(range(x),
                                             sens = photobiology::ciexyzCMF2.spct,
                                             color.name = paste(name, "CMF", sep = ".")),
                    CC  = w_length_range2rgb(range(x),
                                             sens = photobiology::ciexyzCC2.spct,
                                             color.name = paste(name, "CC", sep = ".")))
    } else if (chroma.type == "CMF") {
      color <- w_length_range2rgb(range(x),
                                  sens = photobiology::ciexyzCMF2.spct,
                                  color.name = paste(name, "CMF", sep = "."))
    } else if (chroma.type == "CC") {
      color <- w_length_range2rgb(range(x),
                                  sens = photobiology::ciexyzCC2.spct,
                                  color.name = paste(name, "CC", sep = "."))
    } else {
      color <- NA_character_
    }
  }  else if (is.chroma_spct(chroma.type)) {
    #    warning("Using a CC or CMF definition, RGB may not denote red, green, blue!")
    color <- w_length_range2rgb(range(x),
                                sens = chroma.type,
                                color.name = paste(name, "chroma", sep = "."))
  } else {
    warning("Chroma type of class \"", class(chroma.type)[1], "\" not supported.")
    color <- NA_character_
  }
  return(color)
}

#' @describeIn color_of
#'
#' @export
#'
color_of.source_spct <- function(x,
                                 type = "CMF",
                                 chroma.type = type,
                                 ...) {
  if (length(x) == 0) {
    return(character())
  }
  x.name <- "source"
  q2e(x, byref = TRUE)
  if (is.character(chroma.type)) {
    if (chroma.type == "CMF") {
      s_e_irrad2rgb(x[["w.length"]], x[["s.e.irrad"]],
                    sens = photobiology::ciexyzCMF2.spct,
                    color.name = paste(x.name, "CMF", sep = "."))
    } else if (chroma.type == "CC") {
      s_e_irrad2rgb(x[["w.length"]], x[["s.e.irrad"]],
                    sens = photobiology::ciexyzCC2.spct,
                    color.name = paste(x.name, "CC", sep = "."))
    } else if (chroma.type == "both") {
      c(s_e_irrad2rgb(x[["w.length"]], x[["s.e.irrad"]],
                      sens = photobiology::ciexyzCMF2.spct,
                      color.name = paste(x.name, "CMF", sep = ".")),
        s_e_irrad2rgb(x[["w.length"]], x[["s.e.irrad"]],
                      sens = photobiology::ciexyzCC2.spct,
                      color.name = paste(x.name, "CC", sep = ".")))
    } else {
      warning("Chroma type = ", chroma.type, " not implemented for 'source_spct'.")
      color.out <- NA_character_
      names(color.out) <- paste(x.name, chroma.type, sep = ".")
      color.out
    }
  } else if (is.chroma_spct(chroma.type)) {
#    warning("Using a CC or CMF definition, RGB may not denote red, green, blue!")
    s_e_irrad2rgb(x[["w.length"]],
                  x[["s.e.irrad"]],
                  sens = chroma.type,
                  color.name = paste(x.name, "chroma", sep = "."))
  } else {
    warning("Chroma type of class \"", class(chroma.type)[1], "\" not supported.")
    NA_character_
  }
}

#' @describeIn color_of
#'
#' @export
#'
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#'
#' @export
#'
color_of.source_mspct <- function(x,
                                  ...,
                                  idx = "spct.idx") {
  msdply(mspct = x, color_of, ..., idx = idx, col.names = "color")
}

# compatibility -----------------------------------------------------------

#' @rdname color_of
#'
#' @export
#'
colour_of <- function(x, ...) {
  color_of(x, ...)
}

#' @rdname color_of
#'
#' @export
#'
#' @section Deprecated: Use of color() is deprecated as this wrapper function
#'   may be removed in future versions of the package because of name clashes.
#'   Use color_of() instead.
#'
color <- function(x, ...) {
  message("Method photobiology::color() has been renamed color_of(), color() is deprecated and will be removed.")
  color_of(x, ...)
}
