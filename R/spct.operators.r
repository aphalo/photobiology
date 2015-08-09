# Private function which takes the operator as its third argument
# The operator must enclosed in backticks to be recognized as such
#
# This avoids the repetition of very similar code for each operator
# as in the older versions.

oper.e.generic_spct <- function(e1, e2, oper) {
  if (is.object_spct(e1) || is.object_spct(e2)) {
    stop("Operators are not defined for object_spct objects")
  }
  if (is.logical(e1)) {
    e1 <- as.integer(e1)
  }
  if (is.logical(e2)) {
    e2 <- as.integer(e2)
  }
  if (is.waveband(e1)) {
    e_temp <- e2
    e2 <- e1
    e1 <- e_temp
  }
  if (is.numeric(e1)) {
    class1 <- "numeric"
  } else {
    class1 <- class_spct(e1)[1]
  }
  if (is.numeric(e2)) {
    class2 <- "numeric"
  } else if (is.waveband(e2)) {
    class2 <- "waveband"
  } else {
    class2 <- class_spct(e2)[1]
  }
  if (class1 == "cps_spct") {
    if (is.numeric(e2)) {
    z <- e1
    z[["cps"]] <- oper(z[["cps"]], e2)
    return(z)
    } else if (class2 == "cps_spct") {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$cps, e2$cps, bin.oper=oper, trim="intersection")
    names(z)[2] <- "cps"
    setCpsSpct(z, strict.range = FALSE)
    return(z)
    } else {
      warning("operation between 'cps_spct' and ", class(e2)[1], " objects not implemented")
    }
  } else if (class1 == "source_spct") {
    q2e(e1, action = "add", byref = TRUE)
    if (is.waveband(e2)) {
      if (!identical(oper, `*`)) {
        warning("Only '*' is allowed between source_spct and waveband objects")
        return(NA)
      }
      if (is_effective(e1) && !is_effective(e2)) {
        bwswf.used <- getBSWFUsed(e1)
      } else if (!is_effective(e1) && is_effective(e2)) {
        bswf.used <- labels(e2)[["name"]]
      } else if (is_effective(e1) && is_effective(e2)) {
        bswf.used <- paste(getBSWFUsed(e1), "*", labels(e2)[["name"]])
      } else if (!is_effective(e1) && !is_effective(e2)) {
        bswf.used <- "none"
      } else {
        stop("Failed assertion! BUG IN PACKAGE CODE")
      }
      e1 <- trim_spct(e1, low.limit=min(e2), high.limit=max(e2) - 1e-12, verbose=FALSE,
                      use.hinges = getOption("photobiology.use.hinges", default=NULL))
      mult <-
        calc_multipliers(w.length=e1$w.length, w.band=e2, unit.out="energy",
                         unit.in="energy",
                         use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE))
      return(source_spct(w.length=e1$w.length, s.e.irrad = e1$s.e.irrad * mult,
                           time.unit=getTimeUnit(e1), bswf.used=bswf.used, strict.range = FALSE))
    }
    if (is.numeric(e2)) {
      z <- e1
      if (exists("s.q.irrad", z, inherits = FALSE)) {
        z[["s.q.irrad"]] <- NULL
      }
      z[["s.e.irrad"]] <- oper(z[["s.e.irrad"]], e2)
      return(z)
    } else if (class2 == "source_spct") {
      q2e(e2, action = "replace", byref = TRUE)
      if (is_effective(e1) || is_effective(e2)) {
        warning("One or both operands are effective spectral irradiances")
        bswf.used <- paste(getBSWFUsed(e1), getBSWFUsed(e2))
      } else {
        bswf.used <- "none"
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      names(z)[2] <- "s.e.irrad"
      setSourceSpct(z, time.unit=getTimeUnit(e1), bswf.used = bswf.used, strict.range = FALSE)
      return(z)
    } else if (class2 == "filter_spct") {
      filter.quantity <- getOption("photobiology.filter.qty", default="transmittance")
      if (filter.quantity=="transmittance") {
        if (!identical(oper, `*`) && !identical(oper, `/`)) {
          warning("Only '*' and '/' are allowed between source_spct and filter_spct objects")
          return(NA)
        }
        A2T(e2, action = "replace", byref = TRUE)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$Tfr, bin.oper=oper, trim="intersection")
        names(z)[2] <- "s.e.irrad"
        setSourceSpct(z, time.unit=getTimeUnit(e1), bswf.used = getBSWFUsed(e1), strict.range = FALSE)
      } else if (filter.quantity=="absorbance") {
        if (!identical(oper, `*`) && !identical(oper, `/`)) {
          warning("Only '*' and '/' are allowed between source_spct and filter_spct objects")
          return(NA)
        }
        T2A(e2, action = "add", byref = TRUE)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$A, bin.oper=oper, trim="intersection")
        names(z)[2] <- "s.e.response"
        setResponseSpct(z, time.unit=getTimeUnit(e1))
      }
      return(z)
    } else if (class2 == "reflector_spct") {
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        warning("Only '*' and '/' are allowed between source_spct and reflector_spct objects")
        return(NA)
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      names(z)[2] <- "s.e.irrad"
      setSourceSpct(z, time.unit=getTimeUnit(e1), bswf.used = getBSWFUsed(e1), strict.range = FALSE)
      return(z)
    } else if (class2 == "response_spct") {
      q2e(e2, action = "replace", byref = TRUE)
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        warning("Only '*' and '/' are allowed between source_spct and response_spct objects")
        return(NA)
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.response, bin.oper=oper, trim="intersection")
      names(z)[2] <- "s.e.response"
      setResponseSpct(z, time.unit=getTimeUnit(e1))
      return(z)
    } else if (class2 == "chroma_spct") {
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        warning("Only '*' and '/' are allowed between source_spct and chroma_spct objects")
        return(NA)
      }
      x <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$x, bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$y, bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$z, bin.oper=oper, trim="intersection")
      z <- dplyr::data_frame(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]])
      setChromaSpct(z)
      return(z)
    } else { # this traps also e2 == "generic_spct"
      warning("Operations involving generic_spct are undefined and always return NA")
      return(NA)
    }
  } else if (class1 == "filter_spct") {
    filter.quantity <- getOption("photobiology.filter.qty", default="transmittance")
    if (filter.quantity=="transmittance") {
      e1 <- A2T(e1)
      if (is.numeric(e2)) {
        z <- e1
        if (exists("A", z, inherits = FALSE)) {
          z[["A"]] <- NULL
        }
        z[["Tfr"]] <- oper(z[["Tfr"]], e2)
        return(z)
      }
      else if (class2 == "source_spct") {
        q2e(e2, action = "add", byref = TRUE)
        if (!identical(oper, `*`) && !identical(oper, `/`)) {
          warning("Only '*' and '/' are allowed between filter_spct and source_spct objects")
          return(NA)
        }
        z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.e.irrad, bin.oper=oper, trim="intersection")
        names(z)[2] <- "s.e.irrad"
        setSourceSpct(z, time.unit=getTimeUnit(e2), bswf.used = getBSWFUsed(e2), strict.range = FALSE)
      } else if (class2 == "filter_spct") {
        e2 <- A2T(e2)
        if (!identical(oper, `*`) && !identical(oper, `/`)) {
          warning("Only '*' and '/' are allowed for two filter_spct objects (transmittance)")
          return(NA)
        }
        z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=oper, trim="intersection")
        names(z)[2] <- "Tfr"
        setFilterSpct(z, strict.range = FALSE)
        return(z)
      } else { # this traps optically illegal operations
        warning("The operation attempted in undefined according to Optics laws or the input is malformed")
        return(NA)
      }
    } else if (filter.quantity=="absorbance") {
      T2A(e1, action = "add", byref = TRUE)
      if (is.numeric(e2)) {
        z <- e1
        if (exists("Tfr", z, inherits = FALSE)) {
          z[["Tfr"]] <- NULL
        }
        z[["A"]] <- oper(z[["A"]], e2)
        return(z)
      } else if (class2 == "source_spct") {
        q2e(e2, action = "add", byref = TRUE)
        if (!identical(oper, `*`) && !identical(oper, `/`)) {
          warning("Only '*' and '/' are allowed between filter_spct and source_spct objects")
        }
        z <- oper_spectra(e1$w.length, e2$w.length, e1$A, e2$s.e.irrad, bin.oper=oper, trim="intersection")
        names(z)[2] <- "s.e.response"
        setResponseSpct(z, time.unit=getTimeUnit(e2))
        return(z)
      } else if (class2 == "filter_spct") {
        T2A(e2, action = "add", byref = TRUE)
        if (!identical(oper, `+`) && !identical(oper, `-`)) {
          warning("Only '+' and '-' are allowed for two filter_spct objects (absorbance)")
          return(NA)
        }
        z <- oper_spectra(e1$w.length, e2$w.length, e1$A, e2$A, bin.oper=oper, trim="intersection")
        names(z)[2] <- "A"
        setFilterSpct(z, strict.range = FALSE)
        return(z)
      } else { # this traps optically illegal operations
        warning("The operation attempted in undefined according to Optics laws or the input is malformed")
        return(NA)
      }
    }
  } else if (class1 == "reflector_spct") {
    if (is.numeric(e2)) {
      z <- e1
      z[["Rfr"]] <- oper(z[["Rfr"]], e2)
      return(z)
    } else if (class2 == "reflector_spct") {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$Rfr, bin.oper=oper, trim="intersection")
      names(z)[2] <- "Rfr"
      setReflectorSpct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "source_spct") {
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        warning("Only '*' and '/' are allowed between reflector_spct and source_spct objects")
        return(NA)
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      names(z)[2] <- "s.e.irrad"
      setSourceSpct(z, time.unit=getTimeUnit(e2), bswf.used = getBSWFUsed(e2), strict.range = FALSE)
      return(z)
    } else { # this traps optically illegal operations
      warning("The operation attempted in undefined according to Optics laws or the input is malformed")
      return(NA)
    }
  } else if (class1 == "response_spct") {
    q2e(e1, action = "replace", byref = TRUE)
    if (is.numeric(e2)) {
      z <- e1
      z[["s.e.response"]] <- oper(z[["s.e.response"]], e2)
      return(z)
    } else if (class2 == "response_spct") {
      q2e(e2, action = "replace", byref = TRUE)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.response, e2$s.e.response, bin.oper=oper, trim="intersection")
      names(z)[2] <- "s.e.response"
      setResponseSpct(z, time.unit=getTimeUnit(e1))
      return(z)
    } else if (class2 == "source_spct") {
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        warning("Only '*' and '/' are allowed between response_spct and source_spct objects")
        return(NA)
      }
      q2e(e2, action = "replace", byref = TRUE)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.response, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      names(z)[2] <- "s.e.response"
      setResponseSpct(z, time.unit=getTimeUnit(e2))
      return(z)
    } else { # this traps optically illegal operations
      warning("The operation attempted in undefined according to Optics laws or the input is malformed")
      return(NA)
    }
  } else if (class1 == "chroma_spct") {
    if (is.numeric(e2)) {
      if (length(e2) == 3 && names(e2) == c("x", "y", "z")) {
        return(chroma_spct(w.length = e1$w.length,
                           x = oper(e1$x, e2["x"]),
                           y = oper(e1$y, e2["y"]),
                           z = oper(e1$z, e2["z"])))
        } else {
          return(chroma_spct(w.length = e1$w.length,
                             x = oper(e1$x, e2),
                             y = oper(e1$y, e2),
                             z = oper(e1$z, e2)))
         }
    } else if (class2 == "chroma_spct") {
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$x, bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$y, bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$z, bin.oper=oper, trim="intersection")
      return(chroma_spct(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]]))
    } else if (class2 == "source_spct") {
      q2e(e2, action = "replace", byref = TRUE)
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      return(chroma_spct(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]]))
    }
  } else if (is.numeric(e1)) {
    if (class2 == "cps_spct") {
      z <- e2
      z[["cps"]] <- oper(e1, z[["cps"]])
      return(z)
    } else if (class2 == "source_spct") {
      q2e(e2, action = "replace", byref = TRUE)
      z <- e2
      z[["s.e.irrad"]] <- oper(e1, z[["s.e.irrad"]])
      return(z)
    } else if (class2 == "filter_spct") {
      filter.quantity <- getOption("photobiology.filter.qty", default="transmittance")
      if (filter.quantity=="transmittance") {
        A2T(e2, action = "add", byref = TRUE)
        return(filter_spct(w.length=e2$w.length, Tfr=oper(e1, e2$Tfr), strict.range = FALSE))
      } else if (filter.quantity=="absorbance") {
        T2A(e2, action = "add", byref = TRUE)
        return(filter_spct(w.length=e2$w.length, A=oper(e1, e2$A), strict.range = FALSE))
      } else {
        stop("Assertion failed: bug in code!")
      }
    } else if (class2 == "reflector_spct") {
      return(reflector_spct(w.length=e2$w.length, Rfr=oper(e1, e2$Rfr), strict.range = FALSE))
    } else if (class2 == "response_spct") {
      q2e(e2, action = "replace", byref = TRUE)
      return(response_spct(w.length=e2$w.length, s.e.response=oper(e1, e2$s.e.response)))
    } else if (class2 == "chroma_spct") {
      if (length(e2) == 3 && names(e2) == c("x", "y", "z")) {
        return(chroma_spct(w.length = e2$w.length,
                           x = oper(e1["x"], e2$x),
                           y = oper(e1["y"], e2$y),
                           z = oper(e1["z"], e2$z)))
      } else {
        return(chroma_spct(w.length = e2$w.length,
                           x = oper(e1, e2$x),
                           y = oper(e1, e2$y),
                           z = oper(e1, e2$z)))
      }
    }
  } else {
    warning("The operation attempted in undefined according to Optics laws or the input is malformed")
    return(NA)
  }
}

# Private function which takes the operator as its third argument
# The operator must enclosed in backticks to be recognized as such
#
# This avoids the repetition of very similar code for each operator
# as in the older versions.

oper.q.generic_spct <- function(e1, e2, oper) {
  if (is.object_spct(e1) || is.object_spct(e2)) {
    stop("Operators are not defined for object_spct objects")
  }
  if (is.logical(e1)) {
    e1 <- as.integer(e1)
  }
  if (is.logical(e2)) {
    e2 <- as.integer(e2)
  }
  if (is.waveband(e1)) {
    e_temp <- e2
    e2 <- e1
    e1 <- e_temp
  }
  if (is.numeric(e1)) {
    class1 <- "numeric"
  } else {
    class1 <- class_spct(e1)[1]
  }
  if (is.numeric(e2)) {
    class2 <- "numeric"
  } else if (is.waveband(e2)) {
    class2 <- "waveband"
  } else {
    class2 <- class_spct(e2)[1]
  }
  if (class1 == "cps_spct") {
    if (is.numeric(e2)) {
      z <- e1
      z[["cps"]] <- oper(z[["cps"]], e2)
      return(z)
    } else if (class2 == "cps_spct") {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$cps, e2$cps, bin.oper=oper, trim="intersection")
      names(z)[2] <- "cps"
      setCpsSpct(z, strict.range = FALSE)
      return(z)
    } else {
      warning("operation between 'cps_spct' and ", class(e2)[1], " objects not implemented")
    }
  } else if (class1 == "source_spct") {
    e2q(e1, action = "replace", byref = TRUE)
    if (is.waveband(e2)) {
      if (!identical(oper, `*`)) {
        warning("The operation attempted in undefined according to Optics laws or the input is malformed")
        return(NA)
      }
      if (is_effective(e1) && !is_effective(e2)) {
        bwswf.used <- getBSWFUsed(e1)
      } else if (!is_effective(e1) && is_effective(e2)) {
        bswf.used <- labels(e2)[["name"]]
      } else if (is_effective(e1) && is_effective(e2)) {
        bswf.used <- paste(getBSWFUsed(e1), "*", labels(e2)[["name"]])
      } else if (!is_effective(e1) && !is_effective(e2)) {
        bswf.used <- "none"
      } else {
        stop("Failed assertion! BUG IN PACKAGE CODE")
      }
      e1 <- trim_spct(e1, low.limit = min(e2), high.limit = max(e2) - 1e-12, verbose=FALSE,
                      use.hinges = getOption("photobiology.use.hinges", default=NULL))
      mult <- calc_multipliers(w.length=e1$w.length, w.band=e2, unit.out="photon",
                               unit.in="photon",
                               use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE))
      return(source_spct(w.length=e1$w.length, s.q.irrad = e1$s.q.irrad * mult,
                         bswf.used = bswf.used, time.unit=getTimeUnit(e1),
                         strict.range = FALSE))
     }
    if (is.numeric(e2)) {
      z <- e1
      z[["s.q.irrad"]] <- oper(z[["s.q.irrad"]], e2)
      return(z)
    } else if (class2 == "source_spct") {
      e2q(e2, action = "replace", byref = TRUE)
      if (getTimeUnit(e1) != getTimeUnit(e2)) {
        warning("operands have different value for 'time.unit' attribute")
        return(NA)
      }
      if (is_effective(e1) || is_effective(e2)) {
        warning("One or both operands are effective spectral irradiances")
        bswf.used <- paste(getBSWFUsed(e1), getBSWFUsed(e2))
      } else {
        bswf.used <- "none"
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      names(z)[2] <- "s.q.irrad"
      setSourceSpct(z, time.unit = getTimeUnit(e1), bswf.used = bswf.used, strict.range = FALSE)
      return(z)
    } else if (class2 == "filter_spct") {
      filter.quantity <- getOption("photobiology.filter.qty", default="transmittance")
      if (filter.quantity=="transmittance") {
        if (!identical(oper, `*`) && !identical(oper, `/`)) {
          warning("Only '*' and '/' are allowed between source_spct and filter_spct objects")
          return(NA)
        }
        A2T(e2, action = "add", byref = TRUE)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$Tfr, bin.oper=oper, trim="intersection")
        names(z)[2] <- "s.q.irrad"
        setSourceSpct(z, time.unit = getTimeUnit(e1), bswf.used = getBSWFUsed(e1), strict.range = FALSE)
      } else if (filter.quantity=="absorbance") {
        if (!identical(oper, `*`)) return(NA)
        T2A(e2, action = "add", byref = TRUE)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$A, bin.oper=oper, trim="intersection")
        names(z)[2] <- "s.q.response"
        setResponseSpct(z)
      }
      return(z)
    } else if (class2 == "reflector_spct") {
      if (!identical(oper, `*`)) return(NA)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      names(z)[2] <- "s.q.irrad"
      setSourceSpct(z, time.unit = getTimeUnit(e1), bswf.used = getBSWFUsed(e1), strict.range = FALSE)
      return(z)
    } else if (class2 == "response_spct") {
      if (!identical(oper, `*`)) return(NA)
      e2q(e2, action = "replace", byref = TRUE)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$s.q.response, bin.oper=oper, trim="intersection")
      names(z)[2] <- "s.q.response"
      setResponseSpct(z)
      return(z)
    } else if (class2 == "chroma_spct") {
      if (!identical(oper, `*`)) return(NA)
      x <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$x, bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$y, bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$z, bin.oper=oper, trim="intersection")
      z <- dplyr::data_frame(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]])
      setChromaSpct(z)
      return(z)
    } else { # this traps also e2 == "generic_spct"
      return(NA)
    }
  } else if (class1 == "filter_spct") {
    filter.quantity <- getOption("photobiology.filter.qty", default="transmittance")
    if (filter.quantity=="transmittance") {
      e1 <- A2T(e1)
      if (is.numeric(e2)) {
        z <- e1
        if (exists("A", z, inherits = FALSE)) {
          z[["A"]] <- NULL
        }
        z[["Tfr"]] <- oper(z[["Tfr"]], e2)
        return(z)
      } else if (class2 == "source_spct") {
        e2q(e2, action = "replace", byref = TRUE)
        if (!identical(oper, `*`)) return(NA)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.q.irrad, bin.oper=oper, trim="intersection")
        names(z)[2] <- "s.q.irrad"
        setSourceSpct(z, time.unit = getTimeUnit(e2), bswf.used = getBSWFUsed(e2), strict.range = FALSE)
      } else if (class2 == "filter_spct") {
        e2 <- A2T(e2)
        if (!identical(oper, `*`)) return(NA)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=oper, trim="intersection")
        names(z)[2] <- "Tfr"
        setFilterSpct(z, strict.range = FALSE)
        return(z)
      } else { # this traps optically illegal operations
        warning("The operation attempted in undefined according to Optics laws or the input is malformed")
        return(NA)
      }
    } else if (filter.quantity=="absorbance") {
      T2A(e1, action = "add", byref = TRUE)
      if (is.numeric(e2)) {
        if (exists("Tfr", z, inherits = FALSE)) {
          z[["Tfr"]] <- NULL
        }
        z[["A"]] <- oper(z[["A"]], e2)
        return(z)
      } else if (class2 == "source_spct") {
        e2q(e2, action = "add", byref = TRUE)
        if (!identical(oper, `*`)) return(NA)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$A, e2$s.q.irrad, bin.oper=oper, trim="intersection")
        names(z)[2] <- "s.q.response"
        setResponseSpct(z)
        return(z)
      } else if (class2 == "filter_spct") {
        T2A(e2, action = "add", byref = TRUE)
        if (!identical(oper, `+`) || !identical(oper, `-`)) return(NA)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$A, e2$A, bin.oper=oper, trim="intersection")
        names(z)[2] <- "A"
        setFilterSpct(z, strict.range = FALSE)
        return(z)
      } else { # this traps optically illegal operations
        warning("The operation attempted in undefined according to Optics laws or the input is malformed")
        return(NA)
      }
    }
  } else if (class1 == "reflector_spct") {
    if (is.numeric(e2)) {
      z <- e1
      z[["Rfr"]] <- oper(z[["Rfr"]], e2)
      return(z)
    } else if (class2 == "reflector_spct") {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$Rfr, bin.oper=oper, trim="intersection")
      names(z)[2] <- "Rfr"
      setReflectorSpct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "source_spct") {
      if (!identical(oper, `*`)) return(NA)
      e2q(e2, action = "replace", byref = TRUE)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      names(z)[2] <- "s.q.irrad"
      setSourceSpct(z, time.unit = getTimeUnit(e2), bswf.used = getBSWFUsed(e2), strict.range = FALSE)
      return(z)
    } else { # this traps optically illegal operations
      warning("The operation attempted is undefined according to Optics laws or the input is malformed")
      return(NA)
    }
  } else if (class1 == "response_spct") {
    e2q(e1, action = "replace", byref = TRUE)
    if (is.numeric(e2)) {
      z <- e1
      if (exists("s.e.response", z, inherits = FALSE)) {
        z[["s.e.response"]] <- NULL
      }
      z[["s.q.response"]] <- oper(z[["s.q.response"]], e2)
      return(z)
    }  else if (class2 == "response_spct") {
      e2q(e2, action = "replace", byref = TRUE)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.response, e2$s.q.response, bin.oper=oper, trim="intersection")
      names(z)[2] <- "s.q.response"
      setResponseSpct(z, time.unit=getTimeUnit(e1))
      return(z)
    } else if (class2 == "source_spct") {
      if (!identical(oper, `*`)) return(NA)
      e2q(e2, action = "replace", byref = TRUE)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.response, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      names(z)[2] <- "s.q.response"
      setResponseSpct(z)
      return(z)
    } else { # this traps optically illegal operations
      warning("The operation attempted in undefined according to Optics laws or the input is malformed")
      return(NA)
    }
  } else if (class1 == "chroma_spct") {
    if (is.numeric(e2)) {
      if (length(e2) == 3 && names(e2) == c("x", "y", "z")) {
        return(chroma_spct(w.length = e1$w.length,
                           x = oper(e1$x, e2["x"]),
                           y = oper(e1$y, e2["y"]),
                           z = oper(e1$z, e2["z"])))
      } else {
        return(chroma_spct(w.length = e1$w.length,
                           x = oper(e1$x, e2),
                           y = oper(e1$y, e2),
                           z = oper(e1$z, e2)))
      }
    } else if (class2 == "chroma_spct") {
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$x, bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$y, bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$z, bin.oper=oper, trim="intersection")
      return(chroma_spct(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]]))
    } else if (class2 == "source_spct") {
      e2q(e2, action = "replace", byref = TRUE)
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      return(chroma_spct(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]]))
    }
  } else if (is.numeric(e1)) {
    if (class2 == "cps_spct") {
      z <- e2
      z[["cps"]] <-  oper(e1, z[["cps"]])
      return(z)
    } else if (class2 == "source_spct") {
      e2q(e2, action = "replace", byref = TRUE)
      z <- e2
      z[["s.q.irrad"]] <- oper(e1, z[["s.q.irrad"]])
      return(z)
    } else if (class2 == "filter_spct") {
      filter.quantity <- getOption("photobiology.filter.qty", default="transmittance")
      if (filter.quantity=="transmittance") {
        A2T(e2, action = "add", byref = TRUE)
        return(filter_spct(w.length=e2$w.length, Tfr=oper(e1, e2$Tfr), strict.range = FALSE))
      } else if (filter.quantity=="absorbance") {
        T2A(e2, action = "add", byref = TRUE)
        return(filter_spct(w.length=e2$w.length, A=oper(e1, e2$A), strict.range = FALSE))
      } else {
        stop("Assertion failed: bug in code!")
      }
    } else if (class2 == "reflector_spct") {
      return(reflector_spct(w.length=e2$w.length, Rfr=oper(e1, e2$Rfr), strict.range = FALSE))
    } else if (class2 == "response_spct") {
      e2q(e2, action = "replace", byref = TRUE)
      return(response_spct(w.length=e2$w.length, s.q.response=oper(e1, e2$s.q.response)))
    } else if (class2 == "chroma_spct") {
      if (length(e2) == 3 && names(e2) == c("x", "y", "z")) {
        return(chroma_spct(w.length = e2$w.length,
                           x = oper(e1["x"], e2$x),
                           y = oper(e1["y"], e2$y),
                           z = oper(e1["z"], e2$z)))
      } else {
        return(chroma_spct(w.length = e2$w.length,
                           x = oper(e1, e2$x),
                           y = oper(e1, e2$y),
                           z = oper(e1, e2$z)))
      }
    }
  } else  {
    warning("The operation attempted in undefined according to Optics laws or the input is malformed")
    return(NA)
  }
}

# multiplication ----------------------------------------------------------

#' "*" operator for spectra
#'
#' Multiplication operator for spectra.
#'
#' @param e1 an object of class "generic_spct"
#' @param e2 an object of class "generic_spct"
#' @name times-.generic_spct
#' @export
#' @family math operators and functions
#'
'*.generic_spct' <- function(e1, e2) {
  unit <- getOption("photobiology.radiation.unit", default="energy")
  if (unit == "energy") {
    return(oper.e.generic_spct(e1, e2, `*`))
  } else if (unit == "photon" || unit == "quantum") {
    return(oper.q.generic_spct(e1, e2, `*`))
  } else {
    return(NA)
  }
}

# division ----------------------------------------------------------------


#' "/" operator for spectra
#'
#' Division operator for generic spectra.
#'
#' @param e1 an object of class "generic_spct"
#' @param e2 an object of class "generic_spct"
#' @name slash-.generic_spct
#' @export
#' @family math operators and functions
#'
'/.generic_spct' <- function(e1, e2) {
  unit <- getOption("photobiology.radiation.unit", default="energy")
  if (unit == "energy") {
    return(oper.e.generic_spct(e1, e2, `/`))
  } else if (unit == "photon" || unit == "quantum") {
    return(oper.q.generic_spct(e1, e2, `/`))
  } else {
    return(NA)
  }
}

# Sum ---------------------------------------------------------------

#' "+" operator for spectra
#'
#' Division operator for generic spectra.
#'
#' @param e1 an object of class "generic_spct"
#' @param e2 an object of class "generic_spct"
#' @name plus-.generic_spct
#' @export
#' @family math operators and functions
#'
'+.generic_spct' <- function(e1, e2 = NULL) {
  unit <- getOption("photobiology.radiation.unit", default="energy")
  if (is.null(e2)) {
    e2 <- 0.0
  }
  if (unit == "energy") {
    return(oper.e.generic_spct(e1, e2, `+`))
  } else if (unit == "photon" || unit == "quantum") {
    return(oper.q.generic_spct(e1, e2, `+`))
  } else {
    return(NA)
  }
}

# Minus -------------------------------------------------------------------

#' "-" operator for spectra
#'
#' Substraction operator for generic spectra.
#'
#' @param e1 an object of class "generic_spct"
#' @param e2 an object of class "generic_spct"
#' @name minus-.generic_spct
#' @export
#' @family math operators and functions
#'
'-.generic_spct' <- function(e1, e2 = NULL) {
  unit <- getOption("photobiology.radiation.unit", default="energy")
  if (is.null(e2)) {
    e2 <- e1
    e1 <- 0.0
  }
  if (unit == "energy") {
    return(oper.e.generic_spct(e1, e2, `-`))
  } else if (unit == "photon" || unit == "quantum") {
    return(oper.q.generic_spct(e1, e2, `-`))
  } else {
    return(NA)
  }
}

# other operators  ---------------------------------------------------------------------

#' "^" operator for spectra
#'
#' Power operator for spectra.
#'
#' @param e1 an object of class "generic_spct"
#' @param e2 a numeric vector. possibly of length one.
#' @export
#' @family math operators and functions
#'
'^.generic_spct' <- function(e1, e2) {
  unit <- getOption("photobiology.radiation.unit", default="energy")
  if (unit == "energy") {
    return(oper.e.generic_spct(e1, e2, `^`))
  } else if (unit %in% c("photon", "quantum")) {
    return(oper.q.generic_spct(e1, e2, `^`))
  } else {
    return(NA)
  }
}


# math functions ----------------------------------------------------------

#' math function dispatcher for spectra
#'
#' Function that dispatches the function supplied as argument using different variables depending
#' on the class of the spectrum argument.
#'
#' @param x an object of class "generic_spct"
#' @param f an R function with signature function(x, ...)
#' @param ... additional arguments passed to f
#'
#' @keywords internal
#'
f_dispatcher_spct <- function(x, f, ...) {
  if (is.cps_spct(x)) {
    z <- x
    z[["cps"]] <- f(z[["cps"]], ...)
    return(z)
  } else if (is.filter_spct(x)) {
    filter.qty <- getOption("photobiology.filter.qty", default="transmittance")
    if (filter.qty == "transmittance") {
      z <- A2T(x, action = "replace", byref = FALSE)
      z[["Tfr"]] <- f(z[["Tfr"]], ...)
    } else if (filter.qty == "absorbance") {
      z <- T2A(x, action = "replace", byref = FALSE)
      z[["A"]] <- f(z[["A"]], ...)
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    return(z)
  } else if(is.reflector_spct(x)) {
    z <- x
    z[["Rfr"]] <- f(z[["Rfr"]], ...)
    return(z)
  } else if(is.source_spct(x)) {
    unit <- getOption("photobiology.radiation.unit", default="energy")
    if (unit == "energy") {
      z <- q2e(x, action = "replace", byref = FALSE)
      z[["s.e.irrad"]] <- f(z[["s.e.irrad"]], ...)
      return(z)
    } else if (unit == "photon" || unit == "quantum") {
      z <- e2q(x, action = "replace", byref = FALSE)
      z[["s.q.irrad"]] <- f(z[["s.q.irrad"]], ...)
      return(z)
    } else {
      return(NA)
    }
  } else if (is.response_spct(x)) {
    unit <- getOption("photobiology.radiation.unit", default="energy")
    if (unit == "energy") {
      z <- q2e(x, action = "replace", byref = FALSE)
      z[["s.e.response"]] <- f(z[["s.e.response"]], ...)
      return(z)
    } else if (unit == "photon" || unit == "quantum") {
      z <- e2q(x, action = "replace", byref = FALSE)
      z[["s.q.response"]] <- f(z[["s.q.response"]], ...)
      return(z)
    }
  } else if (is.chroma_spct(x)) {
    z <- x
    z[["x"]] <- f(z[["x"]], ...)
    z[["y"]] <- f(z[["y"]], ...)
    z[["z"]] <- f(z[["z"]], ...)
    return(z)
  } else {
      warning("Function not implemented for ", class(x)[1], " objects.")
    return(NA)
  }
}

#' "log" function for spectra
#'
#' Logarirthm function for spectra.
#'
#' @param x an object of class "generic_spct"
#' @param base a positive number: the base with respect to which logarithms are computed. Defaults to e=exp(1).
#' @export
#' @family math operators and functions
#'
log.generic_spct <- function(x, base = exp(1)) {
  f_dispatcher_spct(x, log, base)
}

#' "log10" function for spectra
#'
#' Base 10 logarithm function for spectra.
#'
#' @param x an object of class "generic_spct"
#' @export
#' @family math operators and functions
#'
log10.generic_spct <- function(x) {
  f_dispatcher_spct(x, log, base = 10)
}

#' "sqrt" function for spectra
#'
#' Square root function for spectra.
#'
#' @param x an object of class "generic_spct"
#' @export
#' @family math operators and functions
#'
sqrt.generic_spct <- function(x) {
  f_dispatcher_spct(x, sqrt)
}

#' "exp" function for spectra
#'
#' Exponential function for spectra.
#'
#' @param x an object of class "generic_spct"
#' @export
#' @family math operators and functions
#'
exp.generic_spct <- function(x) {
  f_dispatcher_spct(x, exp)
}


# transmittance and absorbance --------------------------------------------


# A2T ---------------------------------------------------------------------


#' Generic function
#'
#' Function that coverts absorbance into transmittance (fraction).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @param ... not used in current version
#'
#' @export A2T
#' @family quantity conversion functions
#'
A2T <- function(x, action, byref, ...) UseMethod("A2T")

#' @describeIn A2T Default method for generic function
#'
#' @export
#'
A2T.default <- function(x, action=NULL, byref=FALSE, ...) {
  return(10^-x)
}

#' @describeIn A2T Method for filter spectra
#'
#' @export
#'
A2T.filter_spct <- function(x, action="add", byref=FALSE, ...) {
  if (byref) {
    name <- substitute(x)
  }
  if (exists("Tfr", x, inherits=FALSE)) {
    NULL
  } else if (exists("A", x, inherits=FALSE)) {
    x[["Tfr"]] <- 10^-x[["A"]]
  } else {
    x[["Tfr"]] <- NA
  }
  if (action=="replace" && exists("A", x, inherits=FALSE)) {
    x[["A"]] <- NULL
  }
  if (byref && is.name(name)) { # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

#' @describeIn A2T Method for collections of filter spectra
#'
#' @export
#'
A2T.filter_mspct <- function(x, action = "add", byref = FALSE, ...) {
  mutate_mspct(x, A2T, action = action, byref = byref, ...)
}


# T2A ---------------------------------------------------------------------


#' Convert transmittance into absorbance.
#'
#' Function that coverts transmittance into absorbance (fraction).
#'
#' @param x an R object
#' @param action character Allowed values "replace" and "add"
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @param ... not used in current version
#'
#' @export T2A
#' @family quantity conversion functions
#'
T2A <- function(x, action, byref, ...) UseMethod("T2A")

#' @describeIn T2A Default method for generic function
#'
#' @export
#'
T2A.default <- function(x, action=NULL, byref=FALSE, ...) {
  if (any(x < 0)) {
    Tfr.zero <- getOption("photobiology.Tfr.zero", default = 1e-10)
    warning("Replacing zeros by", Tfr.zero)
    x <- ifelse(x <= 0, Tfr.zero, x)
  }
  return(-log10(x))
}

#' @describeIn T2A Method for filter spectra
#'
#' @export
#'
T2A.filter_spct <- function(x, action="add", byref=FALSE, ...) {
  if (byref) {
    name <- substitute(x)
  }
  if (exists("A", x, inherits=FALSE)) {
    NULL
  } else if (exists("Tfr", x, inherits=FALSE)) {
    x[["A"]] <- -log10(x[["Tfr"]])
  } else {
    x[["A"]] <- NA
  }
  if (action=="replace" && exists("Tfr", x, inherits=FALSE)) {
    x[["Tfr"]] <- NULL
  }
  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  if (any(is.infinite(x$A))) {
    warning("'Inf' absorbance values generated as some Tfr values were equal to zero!")
  }
  return(x)
}

#' @describeIn T2A Method for collections of filter spectra
#'
#' @export
#'
T2A.filter_mspct <- function(x, action = "add", byref = FALSE, ...) {
  mutate_mspct(x, T2A, action = action, byref = byref, ...)
}


# energy - photon and photon - energy conversions -------------------------

# energy to photon ---------------------------------------------------------------------


#' Convert energy-based spectra into photon-based spectra.
#'
#' Function that coverts spectral energy irradiance into spectral photon irradiance (molar).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @param ... not used in current version
#'
#' @export e2q
#' @family quantity conversion functions
#'
e2q <- function(x, action, byref, ...) UseMethod("e2q")

#' @describeIn e2q Default method
#'
#' @export
#'
e2q.default <- function(x, action="add", byref=FALSE, ...) {
  return(NA)
}

#' @describeIn e2q Method for spectral irradiance
#'
#' @export
#'
e2q.source_spct <- function(x, action="add", byref=FALSE, ...) {
  if (byref) {
    name <- substitute(x)
  }
  if (exists("s.q.irrad", x, inherits=FALSE)) {
    NULL
  } else if (exists("s.e.irrad", x, inherits=FALSE)) {
    x[["s.q.irrad"]] <- x[["s.e.irrad"]] * e2qmol_multipliers(x[["w.length"]])
  } else {
    x[["s.q.irrad"]] <- NA
  }
  if (action=="replace" && exists("s.e.irrad", x, inherits=FALSE)) {
    x[["s.e.irrad"]] <- NULL
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
e2q.response_spct <- function(x, action="add", byref=FALSE, ...) {
  if (byref) {
    name <- substitute(x)
  }
  if (exists("s.q.response", x, inherits=FALSE)) {
    NULL
  } else if (exists("s.e.response", x, inherits=FALSE)) {
    x[["s.q.response"]] <- x[["s.e.response"]] / e2qmol_multipliers(x[["w.length"]])
  } else {
    x[["s.q.response"]] <- NA
  }
  if (action=="replace" && exists("s.e.response", x, inherits=FALSE)) {
    x[["s.e.response"]] <- NULL
  }
  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

#' @describeIn e2q Method for collections of (light) source spectra
#'
#' @export
#'
e2q.source_mspct <- function(x, action = "add", byref = FALSE, ...) {
  mutate_mspct(x, e2q, action = action, byref = byref, ...)
}

#' @describeIn e2q Method for for collections of response spectra
#'
#' @export
#'
e2q.response_mspct <- function(x, action = "add", byref = FALSE, ...) {
  mutate_mspct(x, e2q, action = action, byref = byref, ...)
}

# photon to energy ---------------------------------------------------------------------

#' Convert photon-based spectra into energy-based spectra.
#'
#' Function that coverts spectral photon irradiance (molar) into spectral energy irradiance.
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @param ... not used in current version
#'
#' @export q2e
#' @family quantity conversion functions
#'
q2e <- function(x, action, byref, ...) UseMethod("q2e")

#' @describeIn q2e Default method
#'
#' @export
#'
q2e.default <- function(x, action="add", byref=FALSE, ...) {
  return(NA)
}

#' @describeIn q2e Method for spectral irradiance
#'
#' @export
#'
q2e.source_spct <- function(x, action="add", byref=FALSE, ...) {
  if (byref) {
    name <- substitute(x)
  }
  if (exists("s.e.irrad", x, inherits=FALSE)) {
    NULL
  } else if (exists("s.q.irrad", x, inherits=FALSE)) {
    x[["s.e.irrad"]] <- x[["s.q.irrad"]] / e2qmol_multipliers(x[["w.length"]])
  } else {
    x[["s.e.irrad"]] <- NA
  }
  if (action=="replace" && exists("s.q.irrad", x, inherits=FALSE)) {
    x[["s.q.irrad"]] <- NULL
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
q2e.response_spct <- function(x, action="add", byref=FALSE, ...) {
  if (byref) {
    name <- substitute(x)
  }
  if (exists("s.e.response", x, inherits=FALSE)) {
    NULL
  } else if (exists("s.q.response", x, inherits=FALSE)) {
    x[["s.e.response"]] <- x[["s.q.response"]] * e2qmol_multipliers(x[["w.length"]])
  } else {
    x[["s.e.response"]] <- NA
  }
  if (action=="replace" && exists("s.q.response", x, inherits=FALSE)) {
    x[["s.q.response"]] <- NULL
  }
  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

#' @describeIn q2e Method for collections of (light) source spectra
#'
#' @export
#'
q2e.source_mspct <- function(x, action = "add", byref = FALSE, ...) {
  mutate_mspct(x, q2e, action = action, byref = byref, ...)
}


#' @describeIn q2e Method for collections of response spectra
#'
#' @export
#'
q2e.response_mspct <- function(x, action = "add", byref = FALSE, ...) {
  mutate_mspct(x, q2e, action = action, byref = byref, ...)
}

