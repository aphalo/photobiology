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

# transmittance and absorbance --------------------------------------------


# A2T ---------------------------------------------------------------------


#' Convert absorbance into transmittance
#'
#' Function that converts absorbance (a.u.) into transmittance (fraction).
#'
#' @param x an R object
#' @param action a character string
#' @param byref logical indicating if new object will be created by reference or
#'   by copy of \code{x}
#' @param ... not used in current version
#'
#' @return A copy of \code{x} with a column \code{Tfr} added and \code{A} and
#'   \code{Afr} possibly deleted except for \code{w.length}. If \code{action =
#'   "replace"}, in all cases, the additional columns are removed, even if no
#'   column needs to be added.
#'
#' @export
#' @family quantity conversion functions
#'
A2T <- function(x, action, byref, ...) UseMethod("A2T")

#' @describeIn A2T Default method for generic function
#'
#' @export
#'
A2T.default <- function(x, action=NULL, byref = FALSE, ...) {
  warning("'A2T()' not implemented for class \"", class(x)[1], "\".")
  return(x)
}

#' @describeIn A2T method for numeric vectors
#'
#' @export
#'
A2T.numeric <- function(x, action=NULL, byref = FALSE, ...) {
  return(10^-x)
}

#' @describeIn A2T Method for filter spectra
#'
#' @export
#'
A2T.filter_spct <- function(x, action="add", byref = FALSE, ...) {
  if (byref) {
    name <- substitute(x)
  }
  if (is_normalised(x) && !action %in% c("add.raw", "replace.raw")) {
    x <- normalise(x, norm = "update", qty.out = "transmittance")
  } else {
    if (exists("Tfr", x, inherits = FALSE)) {
      NULL
    } else if (exists("A", x, inherits = FALSE)) {
      x[["Tfr"]] <- 10^-x[["A"]]
    } else {
      x[["Tfr"]] <- NA_real_
      action <- "add"
      warning("'A' not available in 'A2T()', ignoring \"replace\" action.")
    }
    if (action=="replace" && exists("A", x, inherits = FALSE)) {
      x[["A"]] <- NULL
    }
    if (action=="replace" && exists("Afr", x, inherits = FALSE)) {
      x[["Afr"]] <- NULL
    }
  }
  if (byref && is.name(name)) { # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

#' @describeIn A2T Method for collections of filter spectra
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
A2T.filter_mspct <- function(x,
                             action = "add",
                             byref = FALSE,
                             ...,
                             .parallel = FALSE,
                             .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          .fun = A2T,
          action = action,
          byref = byref,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}


# T2A ---------------------------------------------------------------------


#' Convert transmittance into absorbance.
#'
#' Function that converts transmittance (fraction) into absorbance (a.u.).
#'
#' @param x an R object
#' @param action character Allowed values "replace" and "add"
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @param clean logical replace off-boundary values before conversion
#' @param ... not used in current version
#'
#' @return A copy of \code{x} with a column \code{A} added and other columns
#'   possibly deleted except for \code{w.length}. If \code{action = "replace"},
#'   in all cases, the additional columns are removed, even if no column needs
#'   to be added.
#'
#' @export
#' @family quantity conversion functions
#'
T2A <- function(x, action, byref, clean, ...) UseMethod("T2A")

#' @describeIn T2A Default method for generic function
#'
#' @export
#'
T2A.default <- function(x, action=NULL, byref = FALSE, ...) {
  warning("'T2A()' not implemented for class \"", class(x)[1], "\".")
  return(x)
}

#' @describeIn T2A Method for numeric vectors
#'
#' @export
#'
T2A.numeric <- function(x, action=NULL, byref = FALSE, clean = TRUE, ...) {
  if (clean) {
    Tfr.zero <- getOption("photobiology.Tfr.zero", default = 0)
    if (any(x < Tfr.zero)) {
      warning("Replacing ", sum(x < Tfr.zero), " values < ",
              Tfr.zero, " with ", Tfr.zero)
      x <- ifelse(x <= 0, Tfr.zero, x)
    }
  }
  return(-log10(x))
}

#' @describeIn T2A Method for filter spectra
#'
#' @export
#'
T2A.filter_spct <- function(x, action="add", byref = FALSE, clean = TRUE, ...) {
  if (byref) {
    name <- substitute(x)
  }
  if (is_normalised(x) && !action %in% c("add.raw", "replace.raw")) {
    x <- normalise(x, norm = "update", qty.out = "absorbance")
  } else {
    if (exists("A", x, inherits = FALSE)) {
      NULL
    } else if (exists("Tfr", x, inherits = FALSE)) {
      if (clean) {
        # we need to avoid infinite recursion
        using_Tfr(x <- clean(x))
      }
      x[["A"]] <- -log10(x[["Tfr"]])
    } else {
      x[["A"]] <- NA_real_
      warning("'Tfr' not available in 'T2A()', filling 'A' with 'NA'.")
      action <- "add"
    }
    if (action=="replace" && exists("Tfr", x, inherits = FALSE)) {
      x[["Tfr"]] <- NULL
    }
    if (action=="replace" && exists("Afr", x, inherits = FALSE)) {
      x[["Afr"]] <- NULL
    }
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  if (any(is.infinite(x[["A"]]))) {
    warning("'Inf' absorbance values generated as some Tfr values were equal to zero!")
  }
  return(x)
}

#' @describeIn T2A Method for collections of filter spectra
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
T2A.filter_mspct <- function(x,
                             action = "add",
                             byref = FALSE,
                             clean = TRUE,
                             ...,
                             .parallel = FALSE,
                             .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          .fun = T2A,
          action = action,
          byref = byref,
          clean = clean,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

# T2Afr ---------------------------------------------------------------------

#' Convert transmittance into absorptance.
#'
#' Function that converts transmittance (fraction) into absorptance (fraction).
#' If reflectance (fraction) is available, it allows conversions between
#' internal and total absorptance.
#'
#' @param x an R object
#' @param action character Allowed values "replace" and "add"
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @param clean logical replace off-boundary values before conversion
#' @param ... not used in current version
#'
#' @return A copy of \code{x} with a column \code{Afr} added and other columns
#'   possibly deleted except for \code{w.length}. If \code{action = "replace"},
#'   in all cases, the additional columns are removed, even if no column needs
#'   to be added.
#'
#' @export
#' @family quantity conversion functions
#'
#' @examples
#' T2Afr(Ler_leaf.spct)
#'
T2Afr <- function(x, action, byref, clean, ...) UseMethod("T2Afr")

#' @describeIn T2Afr Default method for generic function
#'
#' @export
#'
T2Afr.default <- function(x,
                          action = NULL,
                          byref = FALSE,
                          clean = FALSE,
                          ...) {
  warning("'T2Afr()' not implemented for class \"", class(x)[1], "\".")
  x
  }

#' @describeIn T2Afr Default method for generic function
#'
#' @param Rfr numeric vector. Spectral reflectance o reflectance factor.
#'   Set to zero if \code{x} is internal reflectance,
#' @export
#'
T2Afr.numeric <- function(x,
                          action = NULL,
                          byref = FALSE,
                          clean = FALSE,
                          Rfr = NA_real_,
                          ...) {
  if (byref) {
    stop("Conversion by reference not supported for \"numeric\" objects.")
  }
  if (is.na(Rfr)) {
    warning("Convertion requires 'Rfr' to be known.")
  }
  if (any(Rfr > 1) || any(Rfr < 0)) {
    warning("Bad 'Tfr' input valies.")
  }
  if (any(x > 1) || any(x < 0)) {
    warning("Bad 'Tfr' input valies.")
  }
  Tfr.internal <- x * (1 - Rfr)
  1 - Tfr.internal
}

#' @describeIn T2Afr Method for filter spectra
#'
#' @export
#'
T2Afr.filter_spct <- function(x,
                              action = "add",
                              byref = FALSE,
                              clean = FALSE,
                              ...) {
  if (byref) {
    name <- substitute(x)
  }
  current.Tfr.type <- getTfrType(x)
  if (is_normalised(x) && !action %in% c("add.raw", "replace.raw")) {
    x <- normalise(x, norm = "update", qty.out = "absorptance")
  } else {
    if (exists("Afr", x, inherits = FALSE)) {
      NULL
    } else if (exists("Tfr", x, inherits = FALSE)) {
      if (clean) {
        x <- using_Tfr(clean(x))
      }
      if (current.Tfr.type == "total") {
        if (exists("Rfr", x, inherits = FALSE)) {
          x[["Afr"]] <- 1 - x[["Tfr"]] - x[["Rfr"]]
        } else {
          x <- convertTfrType(x, "internal")
          x[["Afr"]] <- 1 - x[["Tfr"]]
          if (all(is.na(x[["Afr"]]))) {
            action <- "add"
            warning("'Tfr.type' or 'Rfr.constant' not available in ')'.")
          }
        }
      } else if (current.Tfr.type == "internal") {
        x[["Afr"]] <- 1 - x[["Tfr"]]
      }
    } else {
      x[["Afr"]] <- NA_real_
      action <- "add"
      warning("'Tfr' not available in 'T2Afr()'.")
    }
    if (action == "replace" && exists("A", x, inherits = FALSE)) {
      x[["A"]] <- NULL
    }
    if (action == "replace" && exists("Tfr", x, inherits = FALSE)) {
      x[["Tfr"]] <- NULL
    }
  }
  if (current.Tfr.type == "total") {
    if (action == "add") {
      x <- convertTfrType(x, Tfr.type = "total")
    } else {
      # no Tfr stored in object, but keep for future conversion operations
      x <- setTfrType(x, "total")
    }
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  x
}

#' @describeIn T2Afr Method for object spectra
#'
#' @export
#'
T2Afr.object_spct <- function(x,
                              action = "add",
                              byref = FALSE,
                              clean = FALSE, ...) {
  if (byref) {
    name <- substitute(x)
  }
  if (exists("Afr", x, inherits = FALSE)) {
    NULL
  } else if (exists("Tfr", x, inherits = FALSE)) {
    if (clean) {
      x <- using_Tfr(clean(x))
    }
    current.Tfr.type <- getTfrType(x)
    if (current.Tfr.type == "internal") {
      x <- convertTfrType(x, "total")
    }
    x[["Afr"]] <- 1 - x[["Tfr"]] - x[["Rfr"]]
    if (current.Tfr.type == "internal") {
      x <- convertTfrType(x, Tfr.type = "internal")
    }
  } else {
    x[["Afr"]] <- NA_real_
    action <- "add"
    warning("'Tfr' not available in 'T2Afr()', ignoring \"replace\" action.")
  }
  if (action == "replace" && exists("Tfr", x, inherits = FALSE)) {
    x[["Tfr"]] <- NULL
  }
  if (action == "replace" && exists("A", x, inherits = FALSE)) {
    x[["A"]] <- NULL
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  x
}

#' @describeIn T2Afr Method for collections of filter spectra
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
T2Afr.filter_mspct <- function(x,
                               action = "add",
                               byref = FALSE,
                               clean = FALSE,
                               ...,
                               .parallel = FALSE,
                               .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          .fun = T2Afr,
          action = action,
          byref = byref,
          clean = FALSE,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn T2Afr Method for collections of object spectra
#'
#' @export
#'
T2Afr.object_mspct <- T2Afr.filter_mspct

# Afr2T ---------------------------------------------------------------------

#' Convert transmittance into absorptance.
#'
#' Function that converts transmittance (fraction) into absorptance (fraction).
#' If reflectance (fraction) is available, it allows conversions between
#' internal and total absorptance.
#'
#' @param x an R object
#' @param action character Allowed values "replace" and "add"
#' @param byref logical indicating if new object will be created by reference or by copy of x
#' @param clean logical replace off-boundary values before conversion
#' @param ... not used in current version
#'
#' @return A copy of \code{x} with a column \code{Tfr} added and other columns
#'   possibly deleted except for \code{w.length}. If \code{action = "replace"},
#'   in all cases, the additional columns are removed, even if no column needs
#'   to be added.
#'
#' @export
#'
#' @family quantity conversion functions
#'
#' @examples
#' T2Afr(Ler_leaf.spct)
#'
Afr2T <- function(x, action, byref, clean, ...) UseMethod("Afr2T")

#' @describeIn Afr2T Default method for generic function
#'
#' @export
#'
Afr2T.default <- function(x,
                          action = NULL,
                          byref = FALSE,
                          clean = FALSE,
                          ...) {
  warning("'Afr2T()' not implemented for class \"", class(x)[1], "\".")
  x
}

#' @describeIn Afr2T Default method for generic function
#'
#' @param Rfr numeric vector. Spectral reflectance o reflectance factor.
#'   Set to zero if \code{x} is internal reflectance,
#' @export
#'
Afr2T.numeric <- function(x,
                          action = NULL,
                          byref = FALSE,
                          clean = FALSE,
                          Rfr = NA_real_,
                          ...) {
  if (byref) {
    stop("Conversion by reference not supported for \"numeric\" objects.")
  }
  if (is.na(Rfr)) {
    warning("Convertion requires 'Rfr'.")
  }
  if (any(Rfr > 1) || any(Rfr < 0)) {
    warning("Bad 'Tfr' input valies.")
  }
  if (any(x > 1) || any(x < 0)) {
    warning("Bad 'Afr' input valies.")
  }
  Afr.internal <- x / (1 - Rfr)
  1 - Afr.internal
}

#' @describeIn Afr2T Method for filter spectra
#'
#' @export
#'
Afr2T.filter_spct <- function(x,
                              action = "add",
                              byref = FALSE,
                              clean = FALSE,
                              ...) {
  if (byref) {
    name <- substitute(x)
  }
  current.Tfr.type <- getTfrType(x)
  if (is_normalised(x) && !action %in% c("add.raw", "replace.raw")) {
    x <- normalise(x, norm = "update", qty.out = "transmittance")
  } else {
    if (exists("Tfr", x, inherits = FALSE)) {
      NULL
    } else if (current.Tfr.type == "internal") {
      # we assume this is what is desired
      x[["Tfr"]] <- 1 - x[["Afr"]]
    } else if (current.Tfr.type == "total") {
      if (exists("Rfr", x, inherits = FALSE)) {
        x[["Tfr"]] <- 1 - x[["Afr"]] - x[["Rfr"]]
      } else {
        properties <- getFilterProperties(x, return.null = FALSE)
        x[["Tfr"]] <- 1 - x[["Afr"]] - properties[["Rfr.constant"]]
      }
    } else {
      stop("Invalid 'Tfr.type' attribute: ", current.Tfr.type)
    }
    if (action == "replace" && exists("Tfr", x, inherits = FALSE)) {
      x[["Afr"]] <- NULL
    }
    if (action == "replace" && exists("A", x, inherits = FALSE)) {
      x[["A"]] <- NULL
    }
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  check_spct(x)
}

#' @describeIn Afr2T Method for object spectra
#'
#' @export
#'
Afr2T.object_spct <- function(x,
                              action = "add",
                              byref = FALSE,
                              clean = FALSE, ...) {
  if (byref) {
    name <- substitute(x)
  }
  current.Tfr.type <- getTfrType(x)
  if (exists("Tfr", x, inherits = FALSE)) {
    NULL
  } else if (current.Tfr.type == "internal") {
    x[["Tfr"]] <- 1 - x[["Afr"]]
  } else if (current.Tfr.type == "total") {
      x[["Tfr"]] <- 1 - x[["Afr"]] - x[["Rfr"]]
  } else {
    stop("Invalid 'Tfr.type' attribute: ", current.Tfr.type)
  }
  if (action == "replace" && exists("Tfr", x, inherits = FALSE)) {
    x[["Afr"]] <- NULL
  }
  if (action == "replace" && exists("A", x, inherits = FALSE)) {
    x[["A"]] <- NULL
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  check_spct(x)
}

#' @describeIn Afr2T Method for collections of filter spectra
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
Afr2T.filter_mspct <- function(x,
                               action = "add",
                               byref = FALSE,
                               clean = FALSE,
                               ...,
                               .parallel = FALSE,
                               .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          .fun = T2Afr,
          action = action,
          byref = byref,
          clean = FALSE,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn Afr2T Method for collections of object spectra
#'
#' @export
#'
Afr2T.object_mspct <- Afr2T.filter_mspct

# "any" filter conversions ------------------------------------------------

#' Convert filter quantities.
#'
#' Functions that convert or add related physical quantities to
#' \code{filter_spct} or  \code{object_spct} objects. transmittance (fraction)
#' into absorptance (fraction).
#'
#' @param x an filter_spct or a filter_mspct object.
#' @param action character Allowed values "replace" and "add".
#' @param clean logical replace off-boundary values before conversion
#'
#' @details These functions are dispatchers for \code{\link{A2T}},
#'   \code{\link{Afr2T}}, \code{\link{T2A}}, and \code{\link{T2Afr}}. The
#'   dispatch is based on the names of the variables stored in \code{x}. They
#'   do not support in-place modification of \code{x}.
#'
#' @return A copy of \code{x} with the columns for the different quantities
#'   added or replaced. If \code{action = "replace"}, in all cases, the
#'   additional columns are removed, even if no column needs to be added.
#'
#' @family quantity conversion functions
#'
#' @export
#'
#' @examples
#' any2Afr(Ler_leaf.spct)
#' any2T(Ler_leaf.spct)
#' any2T(polyester.spct)
#'
any2T <- function(x, action = "add", clean = FALSE) {
  if (is.filter_mspct(x) || is.object_mspct(x)) {
    if (!length(x)) return(x)
    return(msmsply(mspct = x,
                   .fun =  any2T,
                   action = action,
                   clean = clean))
  }
  stopifnot(is.filter_spct(x) || is.object_spct(x))
  if (any(c("A", "Tfr") %in% colnames(x))) {
    A2T(x, action = action, clean = clean, byref = FALSE)
  } else {
    Afr2T(x, action = action, clean = clean, byref = FALSE)
  }
}

#' @rdname any2T
#'
#' @export
#'
any2A <- function(x, action = "add", clean = FALSE) {
  if (is.filter_mspct(x) || is.object_mspct(x)) {
    if (!length(x)) return(x)
    return(msmsply(mspct = x,
                   .fun =  any2A,
                   action = action,
                   clean = clean))
  }
  stopifnot(is.filter_spct(x) || is.object_spct(x))
  if (any(c("A", "Tfr") %in% colnames(x))) {
    T2A(x, action = action, clean = clean, byref = FALSE)
  } else {
    Afr2T(x, action = action)
    T2A(x, action = action, clean = clean, byref = FALSE)
  }
}

#' @rdname any2T
#'
#' @export
#'
any2Afr <- function(x, action = "add", clean = FALSE) {
  if (is.filter_mspct(x) || is.object_mspct(x)) {
    if (!length(x)) return(x)
    return(msmsply(mspct = x,
                   .fun =  any2Afr,
                   action = action,
                   clean = clean))
  }
  stopifnot(is.filter_spct(x) || is.object_spct(x))
  if (any(c("Afr", "Tfr") %in% colnames(x))) {
    T2Afr(x, action = action, clean = clean, byref = FALSE)
  } else {
    A2T(x, action = action)
    T2Afr(x, action = action, clean = clean, byref = FALSE)
  }
}


# energy - photon and photon - energy conversions -------------------------

# energy to photon ---------------------------------------------------------------------


#' Convert energy-based quantities into photon-based quantities.
#'
#' Conversion methods for spectral energy irradiance into spectral photon
#' irradiance and for spectral energy response into spectral photon
#' response.
#'
#' @details The converted spectral values are added to or replace the existing
#'   spectral values depending on the argument passed to parameter
#'   \code{action}. Addition is currently not supported for normalized spectra.
#'   If the spectrum has been normalized with a recent version of package
#'   'photobiology' the spectrum will be renormalized after conversion using the
#'   same arguments as previously. \code{"add.raw"} and \code{"replace.raw"}
#'   prevent the re-normalization, are included for completeness and as a way
#'   of restoring previous behaviour.
#'
#' @param x an R object.
#' @param action a character string, one of "add", "replace", "add.raw" or
#'   "replace.raw".
#' @param byref logical indicating if a new object will be created by reference
#'   or a new object returned.
#' @param ... not used in current version.
#'
#' @export
#' @family quantity conversion functions
#'
e2q <- function(x, action, byref, ...) UseMethod("e2q")

#' @describeIn e2q Default method
#'
#' @export
#'
e2q.default <- function(x, action = "add", byref = FALSE, ...) {
  return(NA)
}

#' @describeIn e2q Method for spectral irradiance
#'
#' @export
#'
e2q.source_spct <- function(x,
                            action = "add",
                            byref = FALSE,
                            ...) {
  if (byref) {
    name <- substitute(x)
  }

  if (is_normalised(x) && !action %in% c("add.raw", "replace.raw")) {
    x <- normalise(x, norm = "update", unit.out = "photon")
  } else {
    if (exists("s.q.irrad", x, inherits = FALSE)) {
      NULL
    } else if (exists("s.e.irrad", x, inherits = FALSE)) {
      x[["s.q.irrad"]] <- x[["s.e.irrad"]] * e2qmol_multipliers(x[["w.length"]])
    } else {
      x[["s.q.irrad"]] <- NA
    }
    if (action %in% c("replace", "replace.raw") &&
        exists("s.e.irrad", x, inherits = FALSE)) {
      x[["s.e.irrad"]] <- NULL
    }
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
e2q.response_spct <- function(x,
                              action = "add",
                              byref = FALSE,
                              ...) {
  if (byref) {
    name <- substitute(x)
  }

  if (is_normalised(x) && !action %in% c("add.raw", "replace.raw")) {
    x <- normalise(x, norm = "update", unit.out = "photon")
  } else {
    if (exists("s.q.response", x, inherits = FALSE)) {
      NULL
    } else if (exists("s.e.response", x, inherits = FALSE)) {
      x[["s.q.response"]] <- x[["s.e.response"]] / e2qmol_multipliers(x[["w.length"]])
    } else {
      x[["s.q.response"]] <- NA
    }
    if (action %in% c("replace", "replace.raw") &&
        exists("s.e.response", x, inherits = FALSE)) {
      x[["s.e.response"]] <- NULL
    }
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

#' @describeIn e2q Method for collections of (light) source spectra
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
e2q.source_mspct <- function(x,
                             action = "add",
                             byref = FALSE,
                             ...,
                             .parallel = FALSE,
                             .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          .fun = e2q,
          action = action,
          byref = byref,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

#' @describeIn e2q Method for collections of response spectra
#'
#' @export
#'
e2q.response_mspct <- function(x,
                               action = "add",
                               byref = FALSE,
                               ...,
                               .parallel = FALSE,
                               .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          .fun = e2q,
          action = action,
          byref = byref,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}

# photon to energy ---------------------------------------------------------------------

#' Convert photon-based quantities into energy-based quantities
#'
#' Conversion methods for spectral photon irradiance into spectral energy
#' irradiance and for spectral photon response into spectral energy
#' response.
#'
#' @details The converted spectral values are added to or replace the existing
#'   spectral values depending on the argument passed to parameter
#'   \code{action}. Addition is currently not supported for normalized spectra.
#'   If the spectrum has been normalized with a recent version of package
#'   'photobiology' the spectrum will be renormalized after conversion using the
#'   same arguments as previously. \code{"add.raw"} and \code{"replace.raw"}
#'   prevent the re-normalization, are included for completeness and as a way of
#'   restoring previous behaviour.
#'
#' @param x an R object.
#' @param action a character string, one of "add", "replace", "add.raw" or
#'   "replace.raw".
#' @param byref logical indicating if a new object will be created by reference
#'   or a new object returned.
#' @param ... not used in current version.
#'
#' @export
#' @family quantity conversion functions
#'
q2e <- function(x, action, byref, ...) UseMethod("q2e")

#' @describeIn q2e Default method
#'
#' @export
#'
q2e.default <- function(x, action = "add", byref = FALSE, ...) {
  return(NA)
}

#' @describeIn q2e Method for spectral irradiance
#'
#' @export
#'
q2e.source_spct <- function(x,
                            action = "add",
                            byref = FALSE,
                            ...) {
  if (byref) {
    name <- substitute(x)
  }

  if (is_normalised(x) && !action %in% c("add.raw", "replace.raw")) {
    x <- normalise(x, norm = "update", unit.out = "energy")
  } else {
    if (exists("s.e.irrad", x, inherits = FALSE)) {
      NULL
    } else if (exists("s.q.irrad", x, inherits = FALSE)) {
      x[["s.e.irrad"]] <- x[["s.q.irrad"]] / e2qmol_multipliers(x[["w.length"]])
    } else {
      x[["s.e.irrad"]] <- NA
    }
    if (action %in% c("replace", "replace.raw") &&
        exists("s.q.irrad", x, inherits = FALSE)) {
      x[["s.q.irrad"]] <- NULL
    }
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
q2e.response_spct <- function(x,
                              action = "add",
                              byref = FALSE,
                              ...) {
  if (byref) {
    name <- substitute(x)
  }

  if (is_normalised(x) && !action %in% c("add.raw", "replace.raw")) {
    x <- normalise(x, norm = "update", unit.out = "energy")
  } else {
    if (exists("s.e.response", x, inherits = FALSE)) {
      NULL
    } else if (exists("s.q.response", x, inherits = FALSE)) {
      x[["s.e.response"]] <- x[["s.q.response"]] * e2qmol_multipliers(x[["w.length"]])
    } else {
      x[["s.e.response"]] <- NA
    }
    if (action %in% c("replace", "replace.raw") &&
        exists("s.q.response", x, inherits = FALSE)) {
      x[["s.q.response"]] <- NULL
    }
  }

  if (byref && is.name(name)) {  # this is a temporary safe net
    name <- as.character(name)
    assign(name, x, parent.frame(), inherits = TRUE)
  }
  return(x)
}

#' @describeIn q2e Method for collections of (light) source spectra
#'
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
q2e.source_mspct <- function(x,
                             action = "add",
                             byref = FALSE,
                             ...,
                             .parallel = FALSE,
                             .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          q2e,
          action = action,
          byref = byref,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}


#' @describeIn q2e Method for collections of response spectra
#'
#' @export
#'
q2e.response_mspct <- function(x,
                               action = "add",
                               byref = FALSE,
                               ...,
                               .parallel = FALSE,
                               .paropts = NULL) {
  if (!length(x)) return(x)
  msmsply(x,
          q2e,
          action = action,
          byref = byref,
          ...,
          .parallel = .parallel,
          .paropts = .paropts)
}


# using options -----------------------------------------------------------

#' Use photobiology options
#'
#' Execute an R expression, possibly compound, using a certain setting for
#' spectral data related options.
#'
#' @param expr an R expression to execute.
#'
#' @return The value returned by the execution of \code{expression}.
#'
#' @references Based on \code{withOptions()} as offered by Thomas Lumley, and
#'   listed in \url{https://www.burns-stat.com/the-options-mechanism-in-r/},
#'   section Deep End, of "The Options mechanism in R" by Patrick Burns.
#'
#' @export
#'
using_Tfr <- function(expr) {
  old <- options(photobiology.filter.qty = "transmittance")
  on.exit(options(old))
  expr <- substitute(expr)
  eval.parent(expr)
}

#' @rdname using_Tfr
#'
#' @export
#'
using_Afr <- function(expr) {
  old <- options(photobiology.filter.qty = "absorptance")
  on.exit(options(old))
  expr <- substitute(expr)
  eval.parent(expr)
}

#' @rdname using_Tfr
#'
#' @export
#'
using_A <- function(expr) {
  old <- options(photobiology.filter.qty = "absorbance")
  on.exit(options(old))
  expr <- substitute(expr)
  eval.parent(expr)
}

#' @rdname using_Tfr
#'
#' @export
#'
using_energy <- function(expr) {
  old <- options(photobiology.radiation.unit = "energy")
  on.exit(options(old))
  expr <- substitute(expr)
  eval.parent(expr)
}

#' @rdname using_Tfr
#'
#' @export
#'
using_photon <- function(expr) {
  old <- options(photobiology.radiation.unit = "photon")
  on.exit(options(old))
  expr <- substitute(expr)
  eval.parent(expr)
}

#' @rdname using_Tfr
#'
#' @export
#'
using_quantum <- using_photon

# Set options -----------------------------------------------------------

#' Set spectral-data options
#'
#' Set spectral-data related options easily.
#'
#' @return Previous value of the modified option.
#'
#' @export
#'
energy_as_default <- function() {
  options(photobiology.radiation.unit = "energy")
}

#' @rdname energy_as_default
#'
#' @export
#'
photon_as_default <- function() {
  options(photobiology.radiation.unit = "photon")
}

#' @rdname energy_as_default
#'
#' @export
#'
quantum_as_default <- photon_as_default

#' @rdname energy_as_default
#'
#' @export
#'
Tfr_as_default <- function() {
  options(photobiology.filter.qty = "transmittance")
}

#' @rdname energy_as_default
#'
#' @export
#'
Afr_as_default <- function() {
  options(photobiology.filter.qty = "absorptance")
}

#' @rdname energy_as_default
#'
#' @export
#'
A_as_default <- function() {
  options(photobiology.filter.qty = "absorbance")
}

#' Set computation options
#'
#' Set computation related options easily.
#'
#' @param flag logical.
#'
#' @return Previous value of the modified option.
#'
#' @export
#'
wb_trim_as_default <- function(flag = TRUE) {
  options(photobiology.waveband.trim = flag)
}

#' @rdname wb_trim_as_default
#'
#' @export
#'
use_cached_mult_as_default <- function(flag = TRUE) {
  options(photobiology.use.cached.mult = flag)
}

#' @rdname energy_as_default
#'
#' @export
#'
unset_radiation_unit_default <- function() {
  options(photobiology.radiation.unit = NULL)
}

#' @rdname energy_as_default
#'
#' @export
#'
unset_filter_qty_default <- function() {
  options(photobiology.filter.qty = NULL)
}

#' @rdname energy_as_default
#'
#' @export
#'
unset_user_defaults <- function() {
  options(photobiology.filter.qty = NULL,
          photobiology.radiation.unit = NULL,
          photobiology.verbose = getOption("verbose"),
          photobiology.strict.range = NULL,
          photobiology.waveband.trim = NULL,
          photobiology.use.cached.mult = NULL)
}

#' Set error reporting options
#'
#' Set error reporting related options easily.
#'
#' @param flag logical.
#'
#' @return Previous value of the modified option.
#'
#' @export
#'
verbose_as_default <- function(flag = TRUE) {
  if (is.null(flag)) {
    flag <- getOption("verbose")
  }
  options(photobiology.verbose = flag)
}

#' @rdname verbose_as_default
#'
#' @export
#'
strict_range_as_default <- function(flag = TRUE) {
  options(photobiology.strict.range = flag)
}
