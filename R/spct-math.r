# Private function which takes the operator as its third argument
# The operator must enclosed in back ticks to be recognized as such
#
# This avoids the repetition of very similar code for each operator
# as in the older versions.

apply_oper <- function(e1, e2, oper) {

  # adjust to current option settings using function variables
  # photon vs. energy
  .unit.irrad <- getOption("photobiology.radiation.unit",
                           default = "energy")
  if (.unit.irrad == "energy") {
    .fun.irrad <- q2e
    .name.irrad <- "s.e.irrad"
    .name.response <- "s.e.response"
  } else if (.unit.irrad == "photon") {
    .fun.irrad <- e2q
    .name.irrad <- "s.q.irrad"
    .name.response <- "s.q.response"
  }
  if (is.source_spct(e1) || is.response_spct(e1)) {
    e1 <- .fun.irrad(e1, action = "replace")
  }
  if (is.source_spct(e2) || is.response_spct(e2)) {
    e2 <- .fun.irrad(e2, action = "replace")
  }
  # filters
  .qty.filter <- getOption("photobiology.filter.qty", default = "transmittance")
  if (.qty.filter == "transmittance") {
    .fun.filter.qty <- any2T
    .name.filter.qty <- "Tfr"
  } else if (.qty.filter == "absorbance") {
    .fun.filter.qty <- any2A
    .name.filter.qty <- "A"
  } else if (.qty.filter == "absorptance") {
    .fun.filter.qty <- any2Afr
    .name.filter.qty <- "Afr"
  }
  # solutes
  if (is.solute_spct(e1)) {
    .name.solute.qty <- intersect(c("K.mole", "K.mass"), names(e1))
    stopifnot(length(.name.filter.qty) == 1)
  } else {
    .name.solute.qty <- NULL
  }
  if (is.solute_spct(e2)) {
    col.name <- intersect(c("K.mole", "K.mass"), names(e2))
    stopifnot(length(.name.filter.qty) == 1)
    if (is.null(.name.solute.qty)) {
      .name.solute.qty <- col.name
    } else if (.name.solute.qty != col.name) {
      stop("Both arguments must be on the same base of expression, mole or mass.")
    }
  }

  # conversion and/or deletion of other quantities
  if (is.filter_spct(e1)) {
    e1 <- .fun.filter.qty(e1, action = "replace")
  }
  if (is.filter_spct(e2)) {
    e2 <- .fun.filter.qty(e2, action = "replace")
  }

  if (is.object_spct(e1) || is.object_spct(e2)) {
    stop("Operators are not defined for \"object_spct\" objects")
  }

  # convert logical to integer
  if (is.logical(e1)) {
    e1 <- as.integer(e1)
  }
  if (is.logical(e2)) {
    e2 <- as.integer(e2)
  }

  # only product defined so transitive
  if (is.waveband(e1)) {
    e_temp <- e2
    e2 <- e1
    e1 <- e_temp
  }

  # save class to avoid multiple calls
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

  # get out of the way cases requiring special handling
  if (is.numeric(e1) && is.any_spct(e2)) {
    if (!length(e2)) return(e2)
    else if (!length(e1)) return(e2[FALSE, ])
  } else if (is.numeric(e2) && is.any_spct(e1)) {
    if (!length(e1)) return(e1)
    else if (!length(e2)) return(e1[FALSE, ])
  }

  # now we walk through all valid combinations of classes
  if (class1 == "calibration_spct") {
    if (is.numeric(e2)) {
      z <- e1
      z[["irrad.mult"]] <- oper(z[["irrad.mult"]], e2)
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "cps_spct" && identical(oper, `*`)) {
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[["irrad.mult"]], e2[["cps"]],
                        bin.oper=oper, trim="intersection")
      names(z)[2] <- .name.irrad
      setSourceSpct(z, strict.range = getOption("photobiology.strict.range",
                                                default = FALSE))
      return(merge_attributes(e1, e2, z))
    } else {
      if (identical(oper, `*`)) {
        stop("operation between 'calibration_spct' and ", class(e2)[1],
             " objects not implemented")
      } else {
        stop("only multiplication between 'calibration_spct' and ",
             "'cps_spct' objects is implemented")
      }
    }

  } else if (class1 == "raw_spct") {
    if (is.numeric(e2)) {
      z <- e1
      z[["counts"]] <- oper(z[["counts"]], e2)
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "raw_spct") {
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[["counts"]], e2[["counts"]],
                        bin.oper=oper, trim="intersection")
      names(z)[2] <- "counts"
      setRawSpct(z, strict.range = getOption("photobiology.strict.range",
                                             default = FALSE))
      return(merge_attributes(e1, e2, z))
    } else {
      stop("operation between 'raw_spct' and ", class(e2)[1],
           " objects not implemented")
    }

  } else if (class1 == "cps_spct") {
    if (is.numeric(e2)) {
      z <- e1
      z[["cps"]] <- oper(z[["cps"]], e2)
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "calibration_spct" && identical(oper, `*`)) {
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[["cps"]], e2[["irrad.mult"]],
                        bin.oper=oper, trim="intersection")
      names(z)[2] <- .name.irrad
      setSourceSpct(z, strict.range = getOption("photobiology.strict.range",
                                                default = FALSE))
      return(merge_attributes(e1, e2, z))
    } else if (class2 == "cps_spct") {
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[["cps"]], e2[["cps"]],
                        bin.oper=oper, trim="intersection")
      names(z)[2] <- "cps"
      setCpsSpct(z, strict.range = getOption("photobiology.strict.range",
                                             default = FALSE))
      return(z)
    } else {
      stop("operation between 'cps_spct' and ", class(e2)[1],
           " objects not implemented")
    }

  } else if (class1 == "source_spct") {
    if (is.waveband(e2)) {
      if (!identical(oper, `*`)) {
        stop("Only '*' is allowed between source_spct and waveband objects")
      }
      if (is_effective(e1) && !is_effective(e2)) {
        bswf.used <- getBSWFUsed(e1)
      } else if (!is_effective(e1) && is_effective(e2)) {
        bswf.used <- labels(e2)[["name"]]
      } else if (is_effective(e1) && is_effective(e2)) {
        bswf.used <- paste(getBSWFUsed(e1), "*", labels(e2)[["name"]])
      } else if (!is_effective(e1) && !is_effective(e2)) {
        bswf.used <- "none"
      } else {
        stop("Failed assertion! BUG IN PACKAGE CODE")
      }
      e1 <- trim_spct(e1, low.limit = min(e2), high.limit = max(e2) - 1e-12,
                      verbose = FALSE,
                      use.hinges = getOption("photobiology.use.hinges",
                                             default = NULL))
      mult <-
        calc_multipliers(w.length = e1[["w.length"]], w.band = e2, unit.out = .unit.irrad,
                         unit.in = .unit.irrad,
                         use.cached.mult = getOption("photobiology.use.cached.mult",
                                                     default = FALSE))
      z <- e1
      z[[.name.irrad]] <- z[[.name.irrad]] * mult
      setBSWFUsed(z, bswf.used = bswf.used)
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (is.numeric(e2)) {
      z <- e1
      z[[.name.irrad]] <- oper(z[[.name.irrad]], e2)
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "source_spct") {
      if (is_effective(e1) || is_effective(e2)) {
        if (getBSWFUsed(e1) != getBSWFUsed(e2)) {
          warning("Operands are weighted quantities, but use a different BSWF.")
          bswf.used <- paste(getBSWFUsed(e1), getBSWFUsed(e2))
        } else {
          bswf.used <- getBSWFUsed(e1)
        }
      } else {
        bswf.used <- "none"
      }
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[[.name.irrad]], e2[[.name.irrad]],
                        bin.oper = oper, trim="intersection")
      names(z)[2] <- .name.irrad
      setSourceSpct(z, time.unit=getTimeUnit(e1), bswf.used = bswf.used,
                    strict.range = getOption("photobiology.strict.range",
                                             default = FALSE))
      return(merge_attributes(e1, e2, z))
    } else if (class2 == "filter_spct") {
      if (!(identical(oper, `*`) || identical(oper, `/`))) {
        stop("Only '*' and '/' are allowed between irradiance and \"Tfr\", \"A\" or \"Afr\" values.")
      }
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[[.name.irrad]], e2[[.name.filter.qty]],
                        bin.oper = oper, trim="intersection")
      names(z)[2] <- .name.irrad
      setSourceSpct(z, time.unit = getTimeUnit(e1), bswf.used = getBSWFUsed(e1),
                    strict.range = getOption("photobiology.strict.range",
                                             default = FALSE))
      return(merge_attributes(e1, e2, z))
    } else if (class2 == "reflector_spct") {
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        stop("Only '*' and '/' are allowed between source_spct and reflector_spct objects")
      }
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[[.name.irrad]], e2[["Rfr"]],
                        bin.oper=oper, trim="intersection")
      names(z)[2] <- .name.irrad
      setSourceSpct(z, time.unit=getTimeUnit(e1), bswf.used = getBSWFUsed(e1),
                    strict.range = getOption("photobiology.strict.range",
                                             default = FALSE))
      return(merge_attributes(e1, e2, z))
    } else if (class2 == "response_spct") {
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        stop("Only '*' and '/' are allowed between source_spct and response_spct objects")
      }
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[[.name.irrad]], e2[[.name.response]],
                        bin.oper=oper, trim="intersection")
      names(z)[2] <- .name.response
      setResponseSpct(z, time.unit=getTimeUnit(e1))
      return(merge_attributes(e1, e2, z))
    } else if (class2 == "chroma_spct") {
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        stop("Only '*' and '/' are allowed between source_spct and chroma_spct objects")
      }
      x <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[[.name.irrad]], e2[["x"]],
                        bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[[.name.irrad]], e2[["y"]],
                        bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[[.name.irrad]], e2[["z"]],
                        bin.oper=oper, trim="intersection")
      z <- tibble::tibble(w.length = x[[1L]],
                          x = x[[2L]], y = y[[2L]], z = z[[2L]])
      setChromaSpct(z)
      return(merge_attributes(e1, e2, z))
    } else { # this traps also e2 == "generic_spct"
      stop("Operations involving generic_spct are undefined and always return NA")
    }

  } else if (class1 == "filter_spct") {
    if (is.numeric(e2)) {
      z <- e1
      z[[.name.filter.qty]] <- oper(z[[.name.filter.qty]], e2)
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "source_spct") {
      if (!(identical(oper, `*`) || identical(oper, `/`))) {
        stop("Only '*' and '/' are allowed between filter_spct and source_spct objects")
      }
      if (identical(oper, `/`)) { # only for verbose!!
        warning("Dividing a 'filter.spct' by a 'source.spct' is a very unusual operation!")
      }
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[[.name.filter.qty]], e2[[.name.irrad]],
                        bin.oper=oper, trim="intersection")
      names(z)[2] <- .name.irrad
      setSourceSpct(z, time.unit = getTimeUnit(e2), bswf.used = getBSWFUsed(e2),
                    strict.range = getOption("photobiology.strict.range",
                                             default = FALSE))
      return(merge_attributes(e1, e2, z))
    } else if (class2 == "filter_spct") {
      if (.qty.filter %in% c("transmittance", "absorptance") &&
          !(identical(oper, `*`) || identical(oper, `/`))) {
        stop("Only '*' and '/' are allowed between \"Tfr\", or \"Afr\" values.")
      } else if (.qty.filter == "absorbance" &&
                 !(identical(oper, `+`) || identical(oper, `-`))) {
        stop("Only '+' and '-' are allowed between \"A\" values.")
      }
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[[.name.filter.qty]], e2[[.name.filter.qty]],
                        bin.oper=oper, trim="intersection")
      names(z)[2] <- .name.filter.qty
      setFilterSpct(z, strict.range = getOption("photobiology.strict.range",
                                                default = FALSE))
      return(merge_attributes(e1, e2, z))
    } else if (class2 == "source_spct") {
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        stop("Only '*' and '/' are allowed between filter_spct and source_spct objects")
      }
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[[.name.filter.qty]], e2[[.name.irrad]],
                        bin.oper=oper, trim="intersection")
      if (.qty.filter == "Tfr") {
        names(z)[2] <- .name.irrad
        setSourceSpct(z, time.unit = getTimeUnit(e2))
      } else {
        names(z)[2] <- .name.response
        setResponseSpct(z, time.unit = getTimeUnit(e2))
      }
      return(merge_attributes(e1, e2, z))
    }

  } else if (class1 == "reflector_spct") {
    if (is.numeric(e2)) {
      z <- e1
      z[["Rfr"]] <- oper(z[["Rfr"]], e2)
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "reflector_spct") {
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[["Rfr"]], e2[["Rfr"]],
                        bin.oper=oper, trim="intersection")
      names(z)[2] <- "Rfr"
      setReflectorSpct(z, strict.range = getOption("photobiology.strict.range",
                                                   default = FALSE))
      return(merge_attributes(e1, e2, z))
    } else if (class2 == "source_spct") {
      if (!identical(oper, `*`) && !identical(oper, `/`)) {
        stop("Only '*' and '/' are allowed between reflector_spct and source_spct objects")
      }
      if (identical(oper, `/`)) { # if verbose!!
        warning("Dividing a 'reflector.spct' by a 'source.spct' is a very unusual operation!")
      }
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[["Rfr"]], e2[[.name.irrad]],
                        bin.oper = oper, trim = "intersection")
      names(z)[2] <- .name.irrad
      setSourceSpct(z, time.unit=getTimeUnit(e2), bswf.used = getBSWFUsed(e2),
                    strict.range = getOption("photobiology.strict.range",
                                             default = FALSE))
      return(merge_attributes(e1, e2, z))
    } else { # this traps optically illegal operations
      stop("The operation attempted is undefined according to Optics laws or the input is malformed")
    }

  } else if (class1 == "solute_spct") {
    if (is.numeric(e2)) {
      z <- e1
      z[[.name.solute.qty]] <- oper(z[[.name.solute.qty]], e2)
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "solute_spct") {
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[[.name.solute.qty]], e2[[.name.solute.qty]],
                        bin.oper = oper, trim = "intersection")
      names(z)[2] <- .name.solute.qty
      setSoluteSpct(z, strict.range = getOption("photobiology.strict.range",
                                                default = FALSE))
      return(merge_attributes(e1, e2, z))
    } else { # this traps optically illegal operations
      stop("The operation attempted is undefined according to Optics laws or the input is malformed")
    }

  } else if (class1 == "response_spct") {
    if (is.numeric(e2)) {
      z <- e1
      z[[.name.response]] <- oper(z[[.name.response]], e2)
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "response_spct") {
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[[.name.response]], e2[[.name.response]],
                        bin.oper=oper, trim="intersection")
      names(z)[2] <- .name.response
      setResponseSpct(z, time.unit=getTimeUnit(e1))
      return(merge_attributes(e1, e2, z))
    } else if (class2 == "source_spct") {
      if (!(identical(oper, `*`) || identical(oper, `/`))) {
        stop("Only '*' and '/' are allowed between response_spct and source_spct objects")
      }
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[[.name.response]], e2[[.name.irrad]],
                        bin.oper=oper, trim="intersection")
      names(z)[2] <- .name.response
      setResponseSpct(z, time.unit=getTimeUnit(e2))
      return(merge_attributes(e1, e2, z))
    } else { # this traps optically illegal operations
      stop("The operation attempted is undefined according to Optics laws or the input is malformed")
    }

  } else if (class1 == "chroma_spct") {
    if (is.numeric(e2)) {
      if (length(e2) == 3 && names(e2) == c("x", "y", "z")) {
        return(chroma_spct(w.length = e1[["w.length"]],
                           x = oper(e1[["x"]], e2["x"]),
                           y = oper(e1[["y"]], e2["y"]),
                           z = oper(e1[["z"]], e2["z"])))
      } else {
        return(chroma_spct(w.length = e1[["w.length"]],
                           x = oper(e1[["x"]], e2),
                           y = oper(e1[["y"]], e2),
                           z = oper(e1[["z"]], e2)))
      }
    } else if (class2 == "chroma_spct") {
      x <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1["x"], e2["x"],
                        bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1["y"], e2["y"],
                        bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1["z"], e2["z"],
                        bin.oper=oper, trim="intersection")
      zz <- chroma_spct(w.length = x[["w.length"]],
                        x = x[["s.irrad"]], y = y[["s.irrad"]], z = z[["s.irrad"]])
      return(merge_attributes(e1, e2, zz))
    } else if (class2 == "source_spct") {
      x <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[["x"]], e2[[.name.irrad]],
                        bin.oper=oper, trim="intersection")
      y <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[["y"]], e2[[.name.irrad]],
                        bin.oper=oper, trim="intersection")
      z <- oper_spectra(e1[["w.length"]], e2[["w.length"]],
                        e1[["z"]], e2[[.name.irrad]],
                        bin.oper=oper, trim="intersection")
      zz <- chroma_spct(w.length = x[["w.length"]],
                        x = x[[2L]], y = y[[2L]], z = z[[2L]])
      return(merge_attributes(e1, e2, zz))
    }

  } else if (is.numeric(e1)) {
    if (class2 == "calibration_spct") {
      z <- e2
      z[["irrad.mult"]] <- oper(e1, z[["irrad.mult"]])
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "raw_spct") {
      z <- e2
      z[["counts"]] <- oper(e1, z[["counts"]])
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "cps_spct") {
      z <- e2
      z[["cps"]] <- oper(e1, z[["cps"]])
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "source_spct") {
      z <- e2
      z[[.name.irrad]] <- oper(e1, z[[.name.irrad]])
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "filter_spct") {
      z <- e2
      z[[.name.filter.qty]] <- oper(e1, z[[.name.filter.qty]])
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "reflector_spct") {
      z <- e2
      z[["Rfr"]] <- oper(e1, e2[["Rfr"]])
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "solute_spct") {
      z <- e2
      z[[.name.solute.qty]] <- oper(e1, e2[[.name.solute.qty]])
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "response_spct") {
      z <- e2
      z[[.name.response]] <- oper(e1, e2[[.name.response]])
      check_spct(z, strict.range = FALSE)
      return(z)
    } else if (class2 == "chroma_spct") {
      z <- e2
      if (length(e1) == 3 && names(e1) == c("x", "y", "z")) {
        z[["x"]] <- oper(e1["x"], z["x"])
        z[["y"]] <- oper(e1["y"], z["y"])
        z[["z"]] <- oper(e1["z"], z["z"])
        check_spct(z, strict.range = FALSE)
        return(z)
      } else {
        z[["x"]] <- oper(e1, z["x"])
        z[["y"]] <- oper(e1, z["y"])
        z[["z"]] <- oper(e1, z["z"])
        check_spct(z, strict.range = FALSE)
        return(z)
      }
    }
  } else {
    stop("The operation attempted is undefined according to Optics laws or the input is malformed")
  }
}

# multiplication ----------------------------------------------------------

#' Arithmetic Operators
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
  apply_oper(e1, e2, `*`)
}

# division ----------------------------------------------------------------


#' Arithmetic Operators
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
  apply_oper(e1, e2, `/`)
}

#' Arithmetic Operators
#'
#' Integer-division operator for generic spectra.
#'
#' @param e1 an object of class "generic_spct"
#' @param e2 an object of class "generic_spct"
#' @name div-.generic_spct
#' @export
#' @family math operators and functions
#'
'%/%.generic_spct' <- function(e1, e2) {
  apply_oper(e1, e2, `%/%`)
}

#' Arithmetic Operators
#'
#' Reminder operator for generic spectra.
#'
#' @param e1 an object of class "generic_spct"
#' @param e2 an object of class "generic_spct"
#' @name mod-.generic_spct
#' @export
#' @family math operators and functions
#'
'%%.generic_spct' <- function(e1, e2) {
  apply_oper(e1, e2, `%%`)
}

# Sum ---------------------------------------------------------------

#' Arithmetic Operators
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
  if (is.null(e2)) {
    # needed to ensure same conversion as if other operator used
    apply_oper(e1, 1, `*`)
  } else {
    apply_oper(e1, e2, `+`)
  }
}

# Minus -------------------------------------------------------------------

#' Arithmetic Operators
#'
#' Subtraction operator for generic spectra.
#'
#' @param e1 an object of class "generic_spct"
#' @param e2 an object of class "generic_spct"
#' @name minus-.generic_spct
#' @export
#' @family math operators and functions
#'
'-.generic_spct' <- function(e1, e2 = NULL) {
  if (is.null(e2)) {
    apply_oper(e1, -1, `*`)
  } else {
    apply_oper(e1, e2, `-`)
  }
}

# other operators  ---------------------------------------------------------------------

#' Arithmetic Operators
#'
#' Power operator for spectra.
#'
#' @param e1 an object of class "generic_spct"
#' @param e2 a numeric vector. possibly of length one.
#' @export
#' @family math operators and functions
#'
'^.generic_spct' <- function(e1, e2) {
  apply_oper(e1, e2, `^`)
}

# math functions ----------------------------------------------------------

#' Math function dispatcher for spectra
#'
#' Function that dispatches the function supplied as argument using different variables depending
#' on the class of the spectrum argument.
#'
#' @param x an object of class "generic_spct"
#' @param .fun an R function with signature function(x, ...)
#' @param ... additional arguments passed to f
#'
#' @keywords internal
#'
f_dispatcher_spct <- function(x, .fun, ...) {
  # Skip checks for intermediate results
  prev_state <- disable_check_spct()
  on.exit(set_check_spct(prev_state), add = TRUE)
  class.x <- class_spct(x)[1]

  # radiation unit
  if (class.x %in% c("source_spct", "response_spct")) {
    .unit.irrad <- getOption("photobiology.radiation.unit",
                             default = "energy")
    if (.unit.irrad == "energy") {
      .fun.irrad <- q2e
      .name.irrad <- "s.e.irrad"
      .name.response <- "s.e.response"
    } else if (.unit.irrad == "photon") {
      .fun.irrad <- e2q
      .name.irrad <- "s.q.irrad"
      .name.response <- "s.q.response"
    }
  }

  # filter qty
  if (class.x == "filter_spct") {
    .qty.filter <- getOption("photobiology.filter.qty", default = "transmittance")
    if (.qty.filter == "transmittance") {
      .fun.filter.qty <- any2T
      .name.filter.qty <- "Tfr"
    } else if (.qty.filter == "absorbance") {
      .fun.filter.qty <- any2A
      .name.filter.qty <- "A"
    } else if (.qty.filter == "absorptance") {
      .fun.filter.qty <- any2Afr
      .name.filter.qty <- "Afr"
    }
  }

  # coeff base
  if (class.x == "solute_spct") {
    .name.solute.qty <- intersect(c("K.mole", "K.mass"), names(x))
    stopifnot(length(.name.solute.qty) == 1)
  }

  # conversion and/or deletion of other quantities
  if (is.source_spct(x) || is.response_spct(x)) {
    z <- .fun.irrad(x, action = "replace")
  } else if (is.filter_spct(x)) {
    z <- .fun.filter.qty(x, action = "replace")
  } else {
    z <- x
  }

  # dispatch
  if (is.calibration_spct(x)) {
    z[["irrad.mult"]] <- .fun(z[["irrad.mult"]], ...)
  } else if (is.raw_spct(x)) {
    z[["counts"]] <- .fun(z[["counts"]], ...)
  } else if (is.cps_spct(x)) {
    z[["cps"]] <- .fun(z[["cps"]], ...)
  } else if (is.filter_spct(x)) {
    z[[.name.filter.qty]] <- .fun(z[[.name.filter.qty]], ...)
  } else if (is.reflector_spct(x)) {
    z[["Rfr"]] <- .fun(z[["Rfr"]], ...)
  } else if (is.solute_spct(x)) {
    z[[.name.solute.qty]] <- .fun(z[[.name.solute.qty]], ...)
  } else if (is.source_spct(x)) {
    z[[.name.irrad]] <- .fun(z[[.name.irrad]], ...)
  } else if (is.response_spct(x)) {
    z[[.name.response]] <- .fun(z[[.name.response]], ...)
  } else if (is.chroma_spct(x)) {
    z[["x"]] <- .fun(z[["x"]], ...)
    z[["y"]] <- .fun(z[["y"]], ...)
    z[["z"]] <- .fun(z[["z"]], ...)
  } else if (is.generic_spct(x)) {
    numeric.cols <- sapply(x, is.numeric)
    col.names <- setdiff(colnames(x)[numeric.cols], "w.length")
    if (length(col.names > 1L)) {
      .name.response <- col.names[1]
      warning("Function applied only to numeric column \"",
              .name.response, "\".")
    } else if (!length(col.names == 0L)) {
      z[[.name.response]] <- .fun(z[[.name.response]], ...)
    }
  } else {
      stop("Function not implemented for ", class(x)[1], " objects.")
  }
  check_spct(z, force = TRUE)
}

#' Logarithms and Exponentials
#'
#' Logarithms and Exponentials for Spectra. The functions are applied to the
#' spectral data, not the wavelengths. The quantity in the spectrum to which the
#' function is applied depends on the class of \code{x} and the current value of
#' output options
#'
#' @name log
#'
#' @param x an object of class "generic_spct"
#' @param base a positive number: the base with respect to which logarithms are
#'   computed. Defaults to e=exp(1).
#'
#' @return An object of the same class as \code{x}.
#'
#' @note In most cases a logarithm of an spectral quantity will yield off-range
#'   values. For this reason unless \code{x} is an object of base class
#'   \code{generic_spct}, checks will not be passed, resulting in warnings or
#'   errors.
#'
#' @export
#' @family math operators and functions
#'
log.generic_spct <- function(x, base = exp(1)) {
  f_dispatcher_spct(x, log, base)
}

#' @rdname log
#'
#' @method log2 generic_spct
#' @export
#'
log2.generic_spct <- function(x) {
  f_dispatcher_spct(x, log, base = 2)
}

#' @rdname log
#'
#' @method log10 generic_spct
#' @export
#'
log10.generic_spct <- function(x) {
  f_dispatcher_spct(x, log, base = 10)
}

#' @rdname log
#'
#' @export
#'
exp.generic_spct <- function(x) {
  f_dispatcher_spct(x, exp)
}

#' Miscellaneous Mathematical Functions
#'
#' \code{abs(x)} computes the absolute value of \code{x}, \code{sqrt(x)}
#' computes the (principal) square root of \code{x}. The functions are applied
#' to the spectral data, not the wavelengths. The quantity in the spectrum to
#' which the function is applied depends on the class of \code{x} and the current
#' value of output options.
#'
#' @name MathFun
#'
#' @param x an object of class "generic_spct"
#' @export
#' @family math operators and functions
#'
sqrt.generic_spct <- function(x) {
  f_dispatcher_spct(x, sqrt)
}

#' @rdname MathFun
#'
#' @export
#'
abs.generic_spct <- function(x) {
  f_dispatcher_spct(x, abs)
}

#' Sign
#'
#' \code{sign} returns a vector with the signs of the corresponding elements of
#' x (the sign of a real number is 1, 0, or -1 if the number is positive, zero,
#' or negative, respectively).
#'
#' @name sign
#'
#' @param x an object of class "generic_spct"
#' @export
#' @family math operators and functions
#'
sign.generic_spct <- function(x) {
  f_dispatcher_spct(x, sign)
}

#' @title
#' Rounding of Numbers
#'
#' @description
#' \code{ceiling} takes a single numeric argument x and returns a numeric vector
#' containing the smallest integers not less than the corresponding elements of
#' x. \\
#' \code{floor} takes a single numeric argument x and returns a numeric vector
#' containing the largest integers not greater than the corresponding elements
#' of x. \\
#' \code{trunc} takes a single numeric argument x and returns a numeric vector
#' containing the integers formed by truncating the values in x toward 0. \\
#' \code{round} rounds the values in its first argument to the specified number of
#' decimal places (default 0). \\
#' \code{signif} rounds the values in its first argument to the specified number of
#' significant digits. \\
#' The functions are applied to the spectral data, not the wavelengths. The
#' quantity in the spectrum to which the function is applied depends on the class
#' of \code{x} and the current value of output options.
#'
#' @param x an object of class "generic_spct" or a derived class.
#' @param digits integer indicating the number of decimal places (round) or
#'   significant digits (signif) to be used. Negative values are allowed (see
#'   'Details').
#'
#' @name round
#' @export
#' @family math operators and functions
#'
round.generic_spct <- function(x, digits = 0) {
  f_dispatcher_spct(x, sign, digits = digits)
}

#' @rdname round
#'
#' @export
#'
signif.generic_spct <- function(x, digits = 6) {
  f_dispatcher_spct(x, signif, digits = digits)
}

#' @rdname round
#'
#' @export
#'
ceiling.generic_spct <- function(x) {
  f_dispatcher_spct(x, ceiling)
}

#' @rdname round
#'
#' @export
#'
floor.generic_spct <- function(x) {
  f_dispatcher_spct(x, floor)
}

#' @rdname round
#'
#' @param ...	arguments to be passed to methods.
#'
#' @export
#'
trunc.generic_spct <- function(x, ...) {
  f_dispatcher_spct(x, trunc, ...)
}

#' Trigonometric Functions
#'
#' Trigonometric functions for object of \code{generic_spct} and derived
#' classes.  \\
#' The functions are applied to the spectral data, not the wavelengths. The
#' quantity in the spectrum to which the function is applied depends on the class
#' of \code{x} and the current value of output options.
#'
#' @name Trig
#'
#' @param x an object of class "generic_spct" or a derived class.
#'
#' @export
#'
cos.generic_spct <- function(x) {
  f_dispatcher_spct(x, cos)
}

#' @rdname Trig
#'
#' @export
#'
sin.generic_spct <- function(x) {
  f_dispatcher_spct(x, sin)
}

#' @rdname Trig
#'
#' @export
#'
tan.generic_spct <- function(x) {
  f_dispatcher_spct(x, tan)
}

#' @rdname Trig
#'
#' @export
#'
acos.generic_spct <- function(x) {
  f_dispatcher_spct(x, acos)
}

#' @rdname Trig
#'
#' @export
#'
asin.generic_spct <- function(x) {
  f_dispatcher_spct(x, asin)
}

#' @rdname Trig
#'
#' @export
#'
atan.generic_spct <- function(x) {
  f_dispatcher_spct(x, atan)
}

