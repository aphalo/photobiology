# Private function which takes the operator as its third argument
# The operator must enclosed in backticks to be recognized as such
#
# This avoids the repetition of very similar code for each operator
# as in the older versions.

oper.e.generic.spct <- function(e1, e2, oper) {
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
    class1 <- class.spct(e1)[1]
  }
  if (is.numeric(e2)) {
    class2 <- "numeric"
  } else if (is.waveband(e2)) {
    class2 <- "waveband"
  } else {
    class2 <- class.spct(e2)[1]
  }
  if (class1 == "source.spct") {
    q2e(e1, action = "add", byref = TRUE)
    if (is.waveband(e2)) {
      if (!identical(oper, `*`)) {
        warning("Only '*' is allowed between source.spct and waveband objects")
        return(NA)
      }
      e1 <- trim_spct(e1, low.limit=min(e2), high.limit=max(e2) - 1e-3, verbose=FALSE, use.hinges = TRUE)
      mult <- calc_multipliers(w.length=e1$w.length, w.band=e2, unit.out="energy",
                               unit.in="energy", use.cached.mult=FALSE)
      if (is_effective(e2)) {
        return(response.spct(w.length=e1$w.length, s.e.response = e1$s.e.irrad * mult, time.unit=attr(e1, "time.unit", exact="TRUE")))
      } else {
        return(source.spct(w.length=e1$w.length, s.e.irrad = e1$s.e.irrad * mult,
                           time.unit=attr(e1, "time.unit", exact="TRUE"), strict.range = FALSE))
      }
    }
    if (is.numeric(e2)) {
      out.spct <- copy(e1)
      if (exists("s.q.irrad", out.spct, inherits = FALSE)) {
        out.spct[ , s.q.irrad := NULL]
      }
      out.spct[ , s.e.irrad := oper(s.e.irrad, e2)]
      return(out.spct)
    } else if (class2 == "source.spct") {
      q2e(e2, action = "add", byref = TRUE)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z, time.unit=attr(e1, "time.unit", exact="TRUE"), strict.range = FALSE)
      return(z)
    } else if (class2 == "filter.spct") {
      filter.quantity <- getOption("photobiology.filter.qty", default="transmittance")
      if (filter.quantity=="transmittance") {
        if (!identical(oper, `*`) && !identical(oper, `/`)) {
          warning("Only '*' and '/' are allowed between source.spct and filter.spct objects")
          return(NA)
        }
        A2T(e2, action = "add", byref = TRUE)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$Tfr, bin.oper=oper, trim="intersection")
        setnames(z, 2, "s.e.irrad")
        setSourceSpct(z, time.unit=attr(e1, "time.unit", exact="TRUE"), strict.range = FALSE)
      } else if (filter.quantity=="absorbance") {
        if (!identical(oper, `*`) && !identical(oper, `/`)) {
          warning("Only '*' and '/' are allowed between source.spct and filter.spct objects")
          return(NA)
        }
        T2A(e2, action = "add", byref = TRUE)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$A, bin.oper=oper, trim="intersection")
        setnames(z, 2, "s.e.response")
        setResponseSpct(z, time.unit=attr(e1, "time.unit", exact="TRUE"))
      }
      return(z)
    } else if (class2 == "reflector.spct") {
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        warning("Only '*' and '/' are allowed between source.spct and reflector.spct objects")
        return(NA)
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z, time.unit=attr(e1, "time.unit", exact="TRUE"), strict.range = FALSE)
      return(z)
    } else if (class2 == "response.spct") {
      q2e(e2, action = "add", byref = TRUE)
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        warning("Only '*' and '/' are allowed between source.spct and response.spct objects")
        return(NA)
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$s.e.response, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.e.response")
      setResponseSpct(z, time.unit=attr(e1, "time.unit", exact="TRUE"))
      return(z)
    } else if (class2 == "chroma.spct") {
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        warning("Only '*' and '/' are allowed between source.spct and chroma.spct objects")
        return(NA)
      }
      x <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$x, bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$y, bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.irrad, e2$z, bin.oper=oper, trim="intersection")
      out.spct <- data.table(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]])
      setChromaSpct(out.spct)
      return(out.spct)
    } else { # this traps also e2 == "generic.spct"
      warning("Operations involving generic.spct are undefined and always return NA")
      return(NA)
    }
  } else if (class1 == "filter.spct") {
    filter.quantity <- getOption("photobiology.filter.qty", default="transmittance")
    if (filter.quantity=="transmittance") {
      e1 <- A2T(e1)
      if (is.numeric(e2)) {
        out.spct <- copy(e1)
        if (exists("A", out.spct, inherits = FALSE)) {
          out.spct[ , A := NULL]
        }
        out.spct[ , Tfr := oper(Tfr, e2)]
        return(out.spct)
      }
      else if (class2 == "source.spct") {
        q2e(e2, action = "add", byref = TRUE)
        if (!identical(oper, `*`) && !identical(oper, `/`)) {
          warning("Only '*' and '/' are allowed between filter.spct and source.spct objects")
          return(NA)
        }
        z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.e.irrad, bin.oper=oper, trim="intersection")
        setnames(z, 2, "s.e.irrad")
        setSourceSpct(z, time.unit=attr(e2, "time.unit", exact="TRUE"), strict.range = FALSE)
      } else if (class2 == "filter.spct") {
        e2 <- A2T(e2)
        if (!identical(oper, `*`) && !identical(oper, `/`)) {
          warning("Only '*' and '/' are allowed for two filter.spct objects (transmittance)")
          return(NA)
        }
        z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=oper, trim="intersection")
        setnames(z, 2, "Tfr")
        setFilterSpct(z, strict.range = FALSE)
        return(z)
      } else { # this traps optically illegal operations
        warning("The operation attempted in undefined according to Optics laws or the input is malformed")
        return(NA)
      }
    } else if (filter.quantity=="absorbance") {
      T2A(e1, action = "add", byref = TRUE)
      if (is.numeric(e2)) {
        out.spct <- copy(e1)
        if (exists("Tfr", out.spct, inherits = FALSE)) {
          out.spct[ , Tfr := NULL]
        }
        out.spct[ , A := oper(A, e2)]
        return(out.spct)
      } else if (class2 == "source.spct") {
        q2e(e2, action = "add", byref = TRUE)
        if (!identical(oper, `*`) && !identical(oper, `/`)) {
          warning("Only '*' and '/' are allowed between filter.spct and source.spct objects")
        }
        z <- oper_spectra(e1$w.length, e2$w.length, e1$A, e2$s.e.irrad, bin.oper=oper, trim="intersection")
        setnames(z, 2, "s.e.response")
        setResponseSpct(z, time.unit=attr(e2, "time.unit", exact="TRUE"))
        return(z)
      } else if (class2 == "filter.spct") {
        T2A(e2, action = "add", byref = TRUE)
        if (!identical(oper, `+`) && !identical(oper, `-`)) {
          warning("Only '+' and '-' are allowed for two filter.spct objects (absorbance)")
          return(NA)
        }
        z <- oper_spectra(e1$w.length, e2$w.length, e1$A, e2$A, bin.oper=oper, trim="intersection")
        setnames(z, 2, "A")
        setFilterSpct(z, strict.range = FALSE)
        return(z)
      } else { # this traps optically illegal operations
        warning("The operation attempted in undefined according to Optics laws or the input is malformed")
        return(NA)
      }
    }
  } else if (class1 == "reflector.spct") {
    if (is.numeric(e2)) {
      out.spct <- copy(e1)
      out.spct[ , Rfr := oper(Rfr, e2)]
      return(out.spct)
    } else if (class2 == "reflector.spct") {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$Rfr, bin.oper=oper, trim="intersection")
      setnames(z, 2, "Rfr")
      setReflectorSpct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "source.spct") {
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        warning("Only '*' and '/' are allowed between reflector.spct and source.spct objects")
        return(NA)
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.e.irrad")
      setSourceSpct(z, time.unit=attr(e2, "time.unit", exact="TRUE"), strict.range = FALSE)
      return(z)
    } else { # this traps optically illegal operations
      warning("The operation attempted in undefined according to Optics laws or the input is malformed")
      return(NA)
    }
  } else if (class1 == "response.spct") {
    q2e(e1, action = "add", byref = TRUE)
    if (is.numeric(e2)) {
      out.spct <- copy(e1)
      if (exists("s.q.response", out.spct, inherits = FALSE)) {
        out.spct[ , s.q.response := NULL]
      }
      out.spct[ , s.e.response := oper(s.e.response, e2)]
      return(out.spct)
    } else if (class2 == "response.spct") {
      e2 <- q2e(e2, action = "add", byref = TRUE)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.response, e2$s.e.response, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.e.response")
      setResponseSpct(z, time.unit=attr(e1, "time.unit", exact="TRUE"))
      return(z)
    } else if (class2 == "source.spct") {
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        warning("Only '*' and '/' are allowed between response.spct and source.spct objects")
        return(NA)
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.e.response, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.e.response")
      setResponseSpct(z, time.unit=attr(e2, "time.unit", exact="TRUE"))
      return(z)
    } else { # this traps optically illegal operations
      warning("The operation attempted in undefined according to Optics laws or the input is malformed")
      return(NA)
    }
  } else if (class1 == "chroma.spct") {
    if (is.numeric(e2)) {
      if (length(e2) == 3 && names(e2) == c("x", "y", "z")) {
        return(chroma.spct(w.length = e1$w.length,
                           x = oper(e1$x, e2["x"]),
                           y = oper(e1$y, e2["y"]),
                           z = oper(e1$z, e2["z"])))
        } else {
          return(chroma.spct(w.length = e1$w.length,
                             x = oper(e1$x, e2),
                             y = oper(e1$y, e2),
                             z = oper(e1$z, e2)))
         }
    } else if (class2 == "chroma.spct") {
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$x, bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$y, bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$z, bin.oper=oper, trim="intersection")
      return(chroma.spct(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]]))
    } else if (class2 == "source.spct") {
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$s.e.irrad, bin.oper=oper, trim="intersection")
      return(chroma.spct(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]]))
    }
  } else if (is.numeric(e1)) {
    if (class2 == "source.spct") {
      out.spct <- copy(e2)
      if (exists("s.q.irrad", out.spct, inherits = FALSE)) {
        out.spct[ , s.q.irrad := NULL]
      }
      out.spct[ , s.e.irrad := oper(e1, s.e.irrad)]
      return(out.spct)
    } else if (class2 == "filter.spct") {
      filter.quantity <- getOption("photobiology.filter.qty", default="transmittance")
      if (filter.quantity=="transmittance") {
        A2T(e2, action = "add", byref = TRUE)
        return(filter.spct(w.length=e2$w.length, Tfr=oper(e1, e2$Tfr), strict.range = FALSE))
      } else if (filter.quantity=="absorbance") {
        T2A(e2, action = "add", byref = TRUE)
        return(filter.spct(w.length=e2$w.length, A=oper(e1, e2$A), strict.range = FALSE))
      } else {
        stop("Assertion failed: bug in code!")
      }
    } else if (class2 == "reflector.spct") {
      return(reflector.spct(w.length=e2$w.length, Rfr=oper(e1, e2$Rfr), strict.range = FALSE))
    } else if (class2 == "response.spct") {
      return(response.spct(w.length=e2$w.length, s.e.response=oper(e1, e2$s.e.response)))
    } else if (class2 == "chroma.spct") {
      if (length(e2) == 3 && names(e2) == c("x", "y", "z")) {
        return(chroma.spct(w.length = e2$w.length,
                           x = oper(e1["x"], e2$x),
                           y = oper(e1["y"], e2$y),
                           z = oper(e1["z"], e2$z)))
      } else {
        return(chroma.spct(w.length = e2$w.length,
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

oper.q.generic.spct <- function(e1, e2, oper) {
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
    class1 <- class.spct(e1)[1]
  }
  if (is.numeric(e2)) {
    class2 <- "numeric"
  } else if (is.waveband(e2)) {
    class2 <- "waveband"
  } else {
    class2 <- class.spct(e2)[1]
  }
  if (class1 == "source.spct") {
    e2q(e1, action = "add", byref = TRUE)
    if (is.waveband(e2)) {
      if (!identical(oper, `*`)) {
        warning("The operation attempted in undefined according to Optics laws or the input is malformed")
        return(NA)
      }
      e1 <- trim_spct(e1, low.limit = min(e2), high.limit = max(e2) - 1e-3, verbose=FALSE, use.hinges = TRUE)
      mult <- calc_multipliers(w.length=e1$w.length, w.band=e2, unit.out="photon",
                               unit.in="photon", use.cached.mult=FALSE)
      if (is_effective(e2)) {
        return(response.spct(w.length=e1$w.length, s.q.response = e1$s.q.irrad * mult, time.unit=attr(e1, "time.unit")))
      } else {
        return(source.spct(w.length=e1$w.length, s.q.irrad = e1$s.q.irrad * mult,
                           time.unit=attr(e1, "time.unit"), strict.range = FALSE))
      }
    }
    if (is.numeric(e2)) {
      out.spct <- copy(e1)
      if (exists("s.e.irrad", out.spct, inherits = FALSE)) {
        out.spct[ , s.e.irrad := NULL]
      }
      out.spct[ , s.q.irrad := oper(s.q.irrad, e2)]
      return(out.spct)
    } else if (class2 == "source.spct") {
      e2q(e2, action = "add", byref = TRUE)
      if (attr(e1, "time.unit") != attr(e2, "time.unit")) {
        warning("operands have different value for 'time.unit' attribute")
        return(NA)
      }
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.q.irrad")
      setSourceSpct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "filter.spct") {
      filter.quantity <- getOption("photobiology.filter.qty", default="transmittance")
      if (filter.quantity=="transmittance") {
        if (!identical(oper, `*`) && !identical(oper, `/`)) {
          warning("Only '*' and '/' are allowed between source.spct and filter.spct objects")
          return(NA)
        }
        A2T(e2, action = "add", byref = TRUE)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$Tfr, bin.oper=oper, trim="intersection")
        setnames(z, 2, "s.q.irrad")
        setSourceSpct(z, strict.range = FALSE)
      } else if (filter.quantity=="absorbance") {
        if (!identical(oper, `*`)) return(NA)
        T2A(e2, action = "add", byref = TRUE)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$A, bin.oper=oper, trim="intersection")
        setnames(z, 2, "s.q.response")
        setResponseSpct(z)
      }
      return(z)
    } else if (class2 == "reflector.spct") {
      if (!identical(oper, `*`)) return(NA)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.q.irrad")
      setSourceSpct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "response.spct") {
      e2q(e2, action = "add", byref = TRUE)
      if (!identical(oper, `*`)) return(NA)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$s.q.response, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.q.response")
      setResponseSpct(z)
      return(z)
    } else if (class2 == "chroma.spct") {
      if (!identical(oper, `*`)) return(NA)
      x <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$x, bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$y, bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.irrad, e2$z, bin.oper=oper, trim="intersection")
      out.spct <- data.table(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]])
      setChromaSpct(out.spct)
      return(out.spct)
    } else { # this traps also e2 == "generic.spct"
      return(NA)
    }
  } else if (class1 == "filter.spct") {
    filter.quantity <- getOption("photobiology.filter.qty", default="transmittance")
    if (filter.quantity=="transmittance") {
      e1 <- A2T(e1)
      if (is.numeric(e2)) {
        out.spct <- copy(e1)
        if (exists("A", out.spct, inherits = FALSE)) {
          out.spct[ , A := NULL]
        }
        out.spct[ , Tfr := oper(Tfr, e2)]
        return(out.spct)
      } else if (class2 == "source.spct") {
        e2q(e2, action = "add", byref = TRUE)
        if (!identical(oper, `*`)) return(NA)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$s.q.irrad, bin.oper=oper, trim="intersection")
        setnames(z, 2, "s.q.irrad")
        setSourceSpct(z, strict.range = FALSE)
      } else if (class2 == "filter.spct") {
        e2 <- A2T(e2)
        if (!identical(oper, `*`)) return(NA)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$Tfr, e2$Tfr, bin.oper=oper, trim="intersection")
        setnames(z, 2, "Tfr")
        setFilterSpct(z, strict.range = FALSE)
        return(z)
      } else { # this traps optically illegal operations
        warning("The operation attempted in undefined according to Optics laws or the input is malformed")
        return(NA)
      }
    } else if (filter.quantity=="absorbance") {
      T2A(e1, action = "add", byref = TRUE)
      if (is.numeric(e2)) {
        if (exists("Tfr", out.spct, inherits = FALSE)) {
          out.spct[ , Tfr := NULL]
        }
        out.spct[ , A := oper(A, e2)]
        return(out.spct)
      } else if (class2 == "source.spct") {
        e2q(e2, action = "add", byref = TRUE)
        if (!identical(oper, `*`)) return(NA)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$A, e2$s.q.irrad, bin.oper=oper, trim="intersection")
        setnames(z, 2, "s.q.response")
        setResponseSpct(z)
        return(z)
      } else if (class2 == "filter.spct") {
        T2A(e2, action = "add", byref = TRUE)
        if (!identical(oper, `+`) || !identical(oper, `-`)) return(NA)
        z <- oper_spectra(e1$w.length, e2$w.length, e1$A, e2$A, bin.oper=oper, trim="intersection")
        setnames(z, 2, "A")
        setFilterSpct(z, strict.range = FALSE)
        return(z)
      } else { # this traps optically illegal operations
        warning("The operation attempted in undefined according to Optics laws or the input is malformed")
        return(NA)
      }
    }
  } else if (class1 == "reflector.spct") {
    if (is.numeric(e2)) {
      out.spct <- copy(e1)
      out.spct[ , Rfr := oper(Rfr, e2)]
      return(out.spct)
    } else if (class2 == "reflector.spct") {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$Rfr, bin.oper=oper, trim="intersection")
      setnames(z, 2, "Rfr")
      setReflectorSpct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "source.spct") {
      if (!identical(oper, `*`)) return(NA)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$Rfr, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.q.irrad")
      setSourceSpct(z, strict.range = FALSE)
      return(z)
    } else { # this traps optically illegal operations
      warning("The operation attempted in undefined according to Optics laws or the input is malformed")
      return(NA)
    }
  } else if (class1 == "response.spct") {
    e2q(e1, action = "add", byref = TRUE)
    if (is.numeric(e2)) {
      out.spct <- copy(e1)
      if (exists("s.e.response", out.spct, inherits = FALSE)) {
        out.spct[ , s.e.response := NULL]
      }
      out.spct[ , s.q.response := oper(s.q.response, e2)]
      return(out.spct)
    }  else if (class2 == "response.spct") {
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.response, e2$s.q.response, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.q.response")
      setResponseSpct(z, time.unit=attr(e1, "time.unit", exact="TRUE"))
      return(z)
    } else if (class2 == "source.spct") {
      if (!identical(oper, `*`)) return(NA)
      z <- oper_spectra(e1$w.length, e2$w.length, e1$s.q.response, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      setnames(z, 2, "s.q.response")
      setResponseSpct(z)
      return(z)
    } else { # this traps optically illegal operations
      warning("The operation attempted in undefined according to Optics laws or the input is malformed")
      return(NA)
    }
  } else if (class1 == "chroma.spct") {
    if (is.numeric(e2)) {
      if (length(e2) == 3 && names(e2) == c("x", "y", "z")) {
        return(chroma.spct(w.length = e1$w.length,
                           x = oper(e1$x, e2["x"]),
                           y = oper(e1$y, e2["y"]),
                           z = oper(e1$z, e2["z"])))
      } else {
        return(chroma.spct(w.length = e1$w.length,
                           x = oper(e1$x, e2),
                           y = oper(e1$y, e2),
                           z = oper(e1$z, e2)))
      }
    } else if (class2 == "chroma.spct") {
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$x, bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$y, bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$z, bin.oper=oper, trim="intersection")
      return(chroma.spct(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]]))
    } else if (class2 == "source.spct") {
      x <- oper_spectra(e1$w.length, e2$w.length, e1$x, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1$w.length, e2$w.length, e1$y, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1$w.length, e2$w.length, e1$z, e2$s.q.irrad, bin.oper=oper, trim="intersection")
      return(chroma.spct(w.length=x$w.length, x=x[["s.irrad"]], y=y[["s.irrad"]], z=z[["s.irrad"]]))
    }
  } else if (is.numeric(e1)) {
    if (class2 == "source.spct") {
      out.spct <- copy(e2)
      if (exists("s.e.irrad", out.spct, inherits = FALSE)) {
        out.spct[ , s.e.irrad := NULL]
      }
      out.spct[ , s.q.irrad := oper(e1, s.q.irrad)]
      return(out.spct)
    } else if (class2 == "filter.spct") {
      filter.quantity <- getOption("photobiology.filter.qty", default="transmittance")
      if (filter.quantity=="transmittance") {
        A2T(e2, action = "add", byref = TRUE)
        return(filter.spct(w.length=e2$w.length, Tfr=oper(e1, e2$Tfr), strict.range = FALSE))
      } else if (filter.quantity=="absorbance") {
        T2A(e2, action = "add", byref = TRUE)
        return(filter.spct(w.length=e2$w.length, A=oper(e1, e2$A), strict.range = FALSE))
      } else {
        stop("Assertion failed: bug in code!")
      }
    } else if (class2 == "reflector.spct") {
      return(reflector.spct(w.length=e2$w.length, Rfr=oper(e1, e2$Rfr), strict.range = FALSE))
    } else if (class2 == "response.spct") {
      return(response.spct(w.length=e2$w.length, s.q.response=oper(e1, e2$s.q.response)))
    } else if (class2 == "chroma.spct") {
      if (length(e2) == 3 && names(e2) == c("x", "y", "z")) {
        return(chroma.spct(w.length = e2$w.length,
                           x = oper(e1["x"], e2$x),
                           y = oper(e1["y"], e2$y),
                           z = oper(e1["z"], e2$z)))
      } else {
        return(chroma.spct(w.length = e2$w.length,
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

#' "*" operator for generic spectra
#'
#' Multiplication operator for generic spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 an object of class "generic.spct"
#' @name times-.generic.spct
#' @export
#'
'*.generic.spct' <- function(e1, e2) {
  unit <- getOption("photobiology.radiation.unit", default="energy")
  if (unit == "energy") {
    return(oper.e.generic.spct(e1, e2, `*`))
  } else if (unit == "photon" || unit == "quantum") {
    return(oper.q.generic.spct(e1, e2, `*`))
  } else {
    return(NA)
  }
}

# division ----------------------------------------------------------------


#' "/" operator for generic spectra
#'
#' Division operator for generic spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 an object of class "generic.spct"
#' @name slash-.generic.spct
#' @export
#'
'/.generic.spct' <- function(e1, e2) {
  unit <- getOption("photobiology.radiation.unit", default="energy")
  if (unit == "energy") {
    return(oper.e.generic.spct(e1, e2, `/`))
  } else if (unit == "photon" || unit == "quantum") {
    return(oper.q.generic.spct(e1, e2, `/`))
  } else {
    return(NA)
  }
}

# Sum ---------------------------------------------------------------

#' "+" operator for generic spectra
#'
#' Division operator for generic spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 an object of class "generic.spct"
#' @name plus-.generic.spct
#' @export
#'
'+.generic.spct' <- function(e1, e2 = NULL) {
  unit <- getOption("photobiology.radiation.unit", default="energy")
  if (is.null(e2)) {
    e2 <- 0.0
  }
  if (unit == "energy") {
    return(oper.e.generic.spct(e1, e2, `+`))
  } else if (unit == "photon" || unit == "quantum") {
    return(oper.q.generic.spct(e1, e2, `+`))
  } else {
    return(NA)
  }
}

# Minus -------------------------------------------------------------------

#' "-" operator for generic spectra
#'
#' Substraction operator for generic spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 an object of class "generic.spct"
#' @name minus-.generic.spct
#' @export
#'
'-.generic.spct' <- function(e1, e2 = NULL) {
  unit <- getOption("photobiology.radiation.unit", default="energy")
  if (is.null(e2)) {
    e2 <- e1
    e1 <- 0.0
  }
  if (unit == "energy") {
    return(oper.e.generic.spct(e1, e2, `-`))
  } else if (unit == "photon" || unit == "quantum") {
    return(oper.q.generic.spct(e1, e2, `-`))
  } else {
    return(NA)
  }
}

# other operators  ---------------------------------------------------------------------

#' "^" operator for spectra
#'
#' Power operator for spectra.
#'
#' @param e1 an object of class "generic.spct"
#' @param e2 a numeric vector. possibly of length one.
#' @export
#'
'^.generic.spct' <- function(e1, e2) {
  unit <- getOption("photobiology.radiation.unit", default="energy")
  if (unit == "energy") {
    return(oper.e.generic.spct(e1, e2, `^`))
  } else if (unit %in% c("photon", "quantum")) {
    return(oper.q.generic.spct(e1, e2, `^`))
  } else {
    return(NA)
  }
}


# math functions ----------------------------------------------------------

#' math function dispatcher for spectra
#'
#' Function that dispatches the fucntion supplied as argument using different variables depending
#' on the class of the spectrum argument.
#'
#' @param x an object of class "generic.spct"
#' @param f an R function with signature function(x, ...)
#' @param ... additional arguments passed to f
#'
f_dispatcher_spct <- function(x, f, ...) {
  if (is(x, "filter.spct")) {
    z <- A2T(x, action = "replace", byref = FALSE)
    z[ , Tfr := f(Tfr, ...)]
    return(z)
  } else if(is(x, "reflector.spct")) {
    z <- copy(x)
    z[ , Rfr := f(Rfr, ...)]
    return(z)
  } else if(is(x, "source.spct")) {
    unit <- getOption("photobiology.radiation.unit", default="energy")
    if (unit == "energy") {
      z <- q2e(x, action = "replace", byref = FALSE)
      z[ , s.e.irrad := f(s.e.irrad, ...)]
      return(z)
    } else if (unit == "photon" || unit == "quantum") {
      z <- e2q(x, action = "replace", byref = FALSE)
      z[ , s.q.irrad := f(s.q.irrad, ...)]
      return(z)
    } else {
      return(NA)
    }
  } else if(is(x, "response.spct")) {
    unit <- getOption("photobiology.radiation.unit", default="energy")
    if (unit == "energy") {
      z <- q2e(x, action = "replace", byref = FALSE)
      z[ , s.e.response := f(s.e.response, ...)]
      return(z)
    } else if (unit == "photon" || unit == "quantum") {
      z <- e2q(x, action = "replace", byref = FALSE)
      z[ , s.q.response := f(s.q.response, ...)]
      return(z)
    }
  } else if (is.chroma.spct(x)) {
    z <- copy(x)
    z[ , `:=` (x = f(x, ...), y = f(y, ...), z = f(z, ...)) ]
    return(z)
  } else {
      stop("Assertion failed! Bug in package photobology's code")
  }
}

#' "log" function for spectra
#'
#' Logarirthm function for spectra.
#'
#' @param x an object of class "generic.spct"
#' @param base a positive number: the base with respect to which logarithms are computed. Defaults to e=exp(1).
#' @export
#'
log.generic.spct <- function(x, base = exp(1)) {
  f_dispatcher_spct(x, log, base)
}

#' "log10" function for spectra
#'
#' Base 10 logarithm function for spectra.
#'
#' @param x an object of class "generic.spct"
#' @export
#'
log10.generic.spct <- function(x) {
  f_dispatcher_spct(x, log, base = 10)
}

#' "sqrt" function for spectra
#'
#' Square root function for spectra.
#'
#' @param x an object of class "generic.spct"
#' @export
#'
sqrt.generic.spct <- function(x) {
  f_dispatcher_spct(x, sqrt)
}

#' "exp" function for spectra
#'
#' Exponential function for spectra.
#'
#' @param x an object of class "generic.spct"
#' @export
#'
exp.generic.spct <- function(x) {
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
#' @export A2T
A2T <- function(x, action, byref) UseMethod("A2T")

#' Default for generic function
#'
#' Function that coverts absorbance into transmittance (fraction).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export A2T.default
A2T.default <- function(x, action=NULL, byref=FALSE) {
  return(10^-x)
}

#' "generic.spct" function
#'
#' Function that coverts absorbance into transmittance (fraction).
#'
#' @param x a "filter.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export A2T.filter.spct
#'
A2T.filter.spct <- function(x, action="add", byref=FALSE) {
  if (byref) {
    name <- substitute(x)
  } else {
    x <- copy(x)
  }
  if (exists("Tfr", x, inherits=FALSE)) {
    NULL
  } else if (exists("A", x, inherits=FALSE)) {
    x[ , Tfr := 10^-A]
  } else {
    x[ , Tfr := NA]
  }
  if (action=="replace" && exists("A", x, inherits=FALSE)) {
    x[ , A := NULL]
  }
  if (byref && is.name(name)) { # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}


# T2A ---------------------------------------------------------------------


#' Generic function
#'
#' Function that coverts transmittance into absorbance (fraction).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export T2A
T2A <- function(x, action, byref) UseMethod("T2A")

#' Default for generic function
#'
#' Function that coverts transmittance into absorbance (fraction).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export T2A.default
T2A.default <- function(x, action=NULL, byref=FALSE) {
  return(-log10(x))
}

#' "filter.spct" function
#'
#' Function that coverts transmittance into absorbance (fraction).
#'
#' @param x a "filter.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export T2A.filter.spct
#'
T2A.filter.spct <- function(x, action="add", byref=FALSE) {
  if (byref) {
    name <- substitute(x)
  } else {
    x <- copy(x)
  }
  if (exists("A", x, inherits=FALSE)) {
    NULL
  } else if (exists("Tfr", x, inherits=FALSE)) {
    x[ , A := -log10(Tfr)]
  } else {
    x[ , A := NA]
  }
  if (action=="replace" && exists("Tfr", x, inherits=FALSE)) {
    x[ , Tfr := NULL]
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


# energy - photon and photon - energy conversions -------------------------

# energy to photon ---------------------------------------------------------------------


#' Generic function
#'
#' Function that coverts spectral energy irradiance into spectral photon irradiance (molar).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export e2q
e2q <- function(x, action, byref) UseMethod("e2q")

#' Default for generic function
#'
#' Function that coverts spectral energy irradiance into spectral photon irradiance (molar).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export e2q.default
e2q.default <- function(x, action="add", byref=FALSE) {
  return(NA)
}

#' "source.spct" function
#'
#' Function that coverts spectral energy irradiance into spectral photon irradiance (molar).
#'
#' @param x a "source.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export e2q.source.spct
#'
e2q.source.spct <- function(x, action="add", byref=FALSE) {
  if (byref) {
    name <- substitute(x)
  } else {
    x <- copy(x)
  }
  if (exists("s.q.irrad", x, inherits=FALSE)) {
    NULL
  } else if (exists("s.e.irrad", x, inherits=FALSE)) {
    x[ , s.q.irrad := s.e.irrad * e2qmol_multipliers(w.length)]
  } else {
    x[ , s.q.irrad := NA]
  }
  if (action=="replace" && exists("s.e.irrad", x, inherits=FALSE)) {
    x[ , s.e.irrad := NULL]
  }
  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

#' "response.spct" function
#'
#' Function that coverts response to spectral energy irradiance into response to spectral photon irradiance (molar).
#'
#' @param x a "response.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export e2q.response.spct
#'
e2q.response.spct <- function(x, action="add", byref=FALSE) {
  if (byref) {
    name <- substitute(x)
  } else {
    x <- copy(x)
  }
  if (exists("s.q.response", x, inherits=FALSE)) {
    NULL
  } else if (exists("s.e.response", x, inherits=FALSE)) {
    x[ , s.q.response := s.e.response / e2qmol_multipliers(w.length)]
  } else {
    x[ , s.q.response := NA]
  }
  if (action=="replace" && exists("s.e.response", x, inherits=FALSE)) {
    x[ , s.e.response := NULL]
  }
  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

# photon to energy ---------------------------------------------------------------------


#' Generic function
#'
#' Function that coverts spectral photon irradiance (molar) into spectral energy irradiance.
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export q2e
q2e <- function(x, action, byref) UseMethod("q2e")

#' Default for generic function
#'
#' Function that coverts spectral photon irradiance (molar) into spectral energy irradiance.
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export q2e.default
q2e.default <- function(x, action="add", byref=FALSE) {
  return(NA)
}

#' "source.spct" function
#'
#' Function that coverts spectral photon irradiance (molar) into spectral energy irradiance.
#'
#' @param x a "source.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export q2e.source.spct
#'
q2e.source.spct <- function(x, action="add", byref=FALSE) {
  if (byref) {
    name <- substitute(x)
  } else {
    x <- copy(x)
  }
  if (exists("s.e.irrad", x, inherits=FALSE)) {
    NULL
  } else if (exists("s.q.irrad", x, inherits=FALSE)) {
    x[ , s.e.irrad := s.q.irrad / e2qmol_multipliers(w.length)]
  } else {
    x[ , s.e.irrad := NA]
  }
  if (action=="replace" && exists("s.q.irrad", x, inherits=FALSE)) {
    x[ , s.q.irrad := NULL]
  }
  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

#' "response.spct" function
#'
#' Function that coverts response to spectral photon irradiance (molar) into response to spectral energy irradiance.
#'
#' @param x a "response.spct"  object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @export q2e.response.spct
#'
q2e.response.spct <- function(x, action="add", byref=FALSE) {
  if (byref) {
    name <- substitute(x)
  } else {
    x <- copy(x)
  }
  if (exists("s.e.response", x, inherits=FALSE)) {
    NULL
  } else if (exists("s.q.response", x, inherits=FALSE)) {
    x[ , s.e.response := s.q.response * e2qmol_multipliers(w.length)]
  } else {
    x[ , s.e.response := NA]
  }
  if (action=="replace" && exists("s.q.response", x, inherits=FALSE)) {
    x[ , s.q.response := NULL]
  }
  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

