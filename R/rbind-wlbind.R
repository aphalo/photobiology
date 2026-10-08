# rbind -------------------------------------------------------------------

#' Row-bind spectra
#'
#' A wrapper on \code{dplyr::rbind_fill} that preserves class and other
#' attributes of spectral objects.
#'
#' @param l A \code{source_mspct}, \code{filter_mspct}, \code{reflector_mspct},
#'   \code{response_mspct}, \code{chroma_mspct}, \code{cps_mspct},
#'   \code{generic_mspct} object or a list containing \code{source_spct},
#'   \code{filter_spct}, \code{reflector_spct}, \code{response_spct},
#'   \code{chroma_spct}, \code{cps_spct}, or \code{generic_spct} objects.
#'
#' @param use.names logical If \code{TRUE} items will be bound by matching
#'   column names. By default \code{TRUE} for \code{rbindspct}. Columns with
#'   duplicate names are bound in the order of occurrence, similar to base. When
#'   TRUE, at least one item of the input list has to have non-null column
#'   names.
#'
#' @param fill logical If \code{TRUE} fills missing columns with NAs. By default
#'   \code{TRUE}. When \code{TRUE}, \code{use.names} has also to be \code{TRUE},
#'   and all items of the input list have to have non-null column names.
#'
#' @param idfactor logical or character Generates an index column of
#'   \code{factor} type. Default is (\code{idfactor=TRUE}) for both lists and
#'   \code{_mspct} objects. If \code{idfactor=TRUE} then the column is auto
#'   named \code{spct.idx}. Alternatively the column name can be directly
#'   provided to \code{idfactor} as a character string.
#'
#' @param attrs.source integer Index into the members of the list from which
#'   attributes should be copied. If \code{NULL}, all attributes are collected
#'   into named lists, except that unique comments are pasted.
#'
#' @param attrs.simplify logical Flag indicating that when all values of an
#'   attribute are equal for all members, the named list will be replaced by
#'   a single copy of the value.
#'
#' @details Each item of \code{l} should be a spectrum, including \code{NULL}
#'   (skipped) or an empty object (0 rows). \code{rbindspc} is most useful when
#'   there are a variable number of (potentially many) objects to stack.
#'   \code{rbindspct} always returns at least a \code{generic_spct} as long as
#'   all elements in l are spectra.
#'
#' @note Note that any additional 'user added' attributes that might exist on
#'   individual items of the input list will not be preserved in the result.
#'   The attributes used by the \code{photobiology} package are preserved, and
#'   if they are not consistent across the bound spectral objects, a warning is
#'   issued.
#'
#' @return An spectral object of a type common to all bound items containing a
#'   concatenation of all the items passed in. If the argument 'idfactor' is
#'   TRUE, then a factor 'spct.idx' will be added to the returned spectral
#'   object.
#'
#' @export
#'
#' @details
#'  \code{dplyr::rbind_fill} is called internally and the result returned is
#'   the highest class in the inheritance hierarchy which is common to all
#'   elements in the list. If not all members of the list belong to one of the
#'   \code{_spct} classes, an error is triggered. The function sets all data in
#'   \code{source_spct} and \code{response_spct} objects supplied as arguments
#'   into energy-based quantities, and all data in \code{filter_spct} objects
#'   into transmittance before the row binding is done. If any member spectrum
#'   is tagged, it is untagged before row binding.
#'
#'   \code{spctbind()} is another name for \code{rbindspct()}.
#'
#' @family Methods for row-binding spectra
#'
#' @examples
#' # default, adds factor 'spct.idx' with letters as levels
#' spct <- rbindspct(list(sun.spct, sun.spct))
#' spct
#' class(spct)
#'
#' # adds factor 'spct.idx' with letters as levels
#' spct <- rbindspct(list(sun.spct, sun.spct), idfactor = TRUE)
#' head(spct)
#' class(spct)
#'
#' # adds factor 'spct.idx' with the names given to the spectra in the list
#' # supplied as formal argument 'l' as levels
#' spct <- rbindspct(list(one = sun.spct, two = sun.spct), idfactor = TRUE)
#' head(spct)
#' class(spct)
#'
#' # adds factor 'ID' with the names given to the spectra in the list
#' # supplied as formal argument 'l' as levels
#' spct <- rbindspct(list(one = sun.spct, two = sun.spct),
#'                   idfactor = "ID")
#' head(spct)
#' class(spct)
#'
rbindspct <- function(l,
                      use.names = TRUE,
                      fill = TRUE,
                      idfactor = TRUE,
                      attrs.source = NULL,
                      attrs.simplify = FALSE) {
  if (is.null(l) || !is.list(l) || length(l) < 1) {
    # _mspct classes are derived from "list"
    warning("Argument 'l' should be a non-empty list or ",
            "a collection of spectra.")
    return(generic_spct())
  }

  if ((is.null(idfactor) &&
       (!is.null(names(l)))) || is.logical(idfactor) && idfactor) {
    idfactor <- "spct.idx"
  } else if (is.logical(idfactor) && !idfactor) {
    idfactor <- NULL
  }

  # inefficient but simpler to implement, and ensures proper naming
  # make sure each member spct object contains a single spectrum
  if (any(sapply(l, getMultipleWl) > 1L)) {
    l <- subset2mspct(l)
  }

  # we skip spectra with no rows
  selector <- unname(sapply(l, nrow)) > 0

  if (use.names && !rlang::is_named(l)) {
    names(l) <- paste("spct", seq_along(l), sep = "_")
  }
  add.idfactor <- is.character(idfactor)

  # We find the most derived common class for spectra
  l.class <- shared_member_class(l)
  if (length(l.class) < 1L) {
    stop("Argument 'l' should contain spectra.")
  } else {
    l.class <- l.class[1L]
  }
  if (!any(selector)) {
    return(do.call(what = l.class, args = list()))
  }
  if (length(l[selector]) == 1L) {
    z <- l[selector][[1L]]
    if (add.idfactor) {
      z[[idfactor]] <- factor(rep(names(l[selector]), times = nrow(z)))
      setIdFactor(z, idfactor)
    }
    return(z)
  }

  # list may have members which already have multiple spectra in long form
  mltpl.wl <- sum(sapply(l, FUN = getMultipleWl))

  # we check that all spectral data contain consistent quantities
  if (l.class %in% c("source_spct", "response_spct")) {
    photon.based <- sapply(l, FUN = is_photon_based)
    energy.based <- sapply(l, FUN = is_energy_based)
    qe.consistent.based <-
      all(photon.based) && !any(energy.based) ||
      all(energy.based) && !any(photon.based) ||
      all(energy.based) && all(photon.based)
  } else {
    qe.consistent.based <- NA
  }

  if (l.class == "filter_spct") {
    absorbance.based <- sapply(l, FUN = is_absorbance_based)
    transmittance.based <- sapply(l, FUN = is_transmittance_based)
    absorptance.based <- sapply(l, FUN = is_absorptance_based)
    TA.consistent.based <- all(absorbance.based) ||
      all(absorptance.based) ||
      all(transmittance.based)
  } else {
    TA.consistent.based <- NA
  }

  # check for transformed data
  scaled.input <- sapply(l, FUN = is_scaled)
  normalized.input <- sapply(l, FUN = is_normalized)
  effective.input <- sapply(l, FUN = is_effective)

  if (any(scaled.input) && !all(scaled.input)) {
    warning("Spectra being row-bound have been differently re-scaled")
  }
  if (any(normalized.input) && length(unique(normalized.input)) > 1L) {
    warning("Spectra being row-bound have been differently normalized")
  }

  for (i in seq_along(l)) {
    class_spct <- class(l[[i]])[1]
    l.class <- intersect(l.class, class_spct)
    if (is_tagged(l[[i]])) {
      l[[i]] <- untag(l[[i]])
    }
    if (!is.na(qe.consistent.based) && !qe.consistent.based) {
      l[[i]] <- q2e(l[[i]], action = "replace", byref = FALSE)
    }
    if (!is.na(TA.consistent.based) && !TA.consistent.based) {
      l[[i]] <- A2T(l[[i]], action = "replace", byref = FALSE)
    }
  }

  # check class is same for all spectra
  #  print(l.class)
  if (length(l.class) != 1L) {
    stop("All spectra in 'l' should belong to the same spectral class.")
  }

  # Here we do the actual binding
  if (length(l) == 1) {
    ans <- l[[1]]
  } else {
    ans <- plyr::rbind.fill(l)
    ans <- tibble::as_tibble(ans)
  }
  if (is.null(ans)) {
    return(generic_spct())
  }

  names.spct <- names(l)
  if (is.null(names.spct) || anyNA(names.spct) || length(names.spct) < length(l)) {
    names.spct <- paste("spct", seq_along(l), sep = "_")
  } else {
    if (anyDuplicated(names.spct)) {
      warning("Duplicated member names have been de-ambiguated before binding spectra.")
      names.spct <- make.unique(names.spct, sep = "_")
      names(l) <- names.spct
    }
  }
  if (add.idfactor) {
    ans[[idfactor]] <- factor(rep(names.spct, times = sapply(l, FUN = nrow)),
                              levels = names.spct)
  }

  comment.ans <- "rbindspct: concatenated comments"
  comments.found <- FALSE

  if (length(attrs.source)) {
    idxs <- intersect(seq_along(l), attrs.source)
  } else {
    idxs <- seq_along(l)
  }

  # get methods and functions return NA if attr is not set
  if (length(idxs) == 1L) {
    comment.ans <- comment(l[[idxs]])
    instr.desc <- getInstrDesc(l[[idxs]])
    instr.settings <- getInstrSettings(l[[idxs]])
    when.measured <- getWhenMeasured(l[[idxs]])
    where.measured <- getWhereMeasured(l[[idxs]])
    what.measured <- getWhatMeasured(l[[idxs]])
    how.measured <- getHowMeasured(l[[idxs]])
    normalized <- getNormalized(l[[idxs]])
    normalization <- getNormalization(l[[idxs]])
  } else {
    # we avoid duplicating the attributes when possible
    comments <- lapply(l[idxs], comment)
    comment.ans <- paste(unique(comments))

    instr.desc <- lapply(l[idxs], getInstrDesc)
    if (attrs.simplify && length(unique(instr.desc)) == 1) {
      instr.desc <- instr.desc[[1]]
    } else {
      names(instr.desc) <- names.spct[idxs]
    }

    instr.settings <- lapply(l[idxs], getInstrSettings)
    if (attrs.simplify && length(unique(instr.settings)) == 1) {
      instr.settings <- instr.settings[[1]]
    } else {
      names(instr.settings) <- names.spct[idxs]
    }

    when.measured <- lapply(l[idxs], getWhenMeasured)
    names(when.measured) <- names.spct[idxs]

    where.measured <- lapply(l[idxs], getWhereMeasured)
    if (attrs.simplify &&
        (all(is.na(where.measured$lon)) ||
         length(unique(where.measured$lon)) == 1) &&
        (all(is.na(where.measured$lat)) ||
         length(unique(where.measured$lat)) == 1) &&
        (all(is.na(where.measured$address)) ||
         length(unique(where.measured$address)) == 1)) {
      where.measured <- where.measured[[1]]
    } else {
      names(where.measured) <- names.spct[idxs]
    }

    what.measured <- lapply(l[idxs], getWhatMeasured)
    if (attrs.simplify && length(unique(what.measured)) == 1) {
      what.measured <- what.measured[[1]]
    } else {
      names(what.measured) <- names.spct[idxs]
    }

    how.measured <- lapply(l[idxs], getHowMeasured)
    if (attrs.simplify && length(unique(how.measured)) == 1) {
      how.measured <- how.measured[[1]]
    } else {
      names(how.measured) <- names.spct[idxs]
    }

    normalized <- lapply(l[idxs], getNormalized)
    names(normalized) <- names.spct[idxs]

    normalization <- lapply(l[idxs], getNormalization)
    names(normalization) <- names.spct[idxs]

  }

  if (l.class == "source_spct") {
    time.unit <- sapply(l, FUN = getTimeUnit)
    names(time.unit) <- NULL
    time.unit <- unique(time.unit)
    if (length(time.unit) > 1L) {
      warning("Inconsistent time units among source spectra ",
              "passed to rbindspct")
      return(source_spct())
    }
    if (any(effective.input)) {
      bswfs.input <- sapply(l, FUN = getBSWFUsed)
      if (length(unique(bswfs.input)) > 1L) {
        bswf.used <- "multiple"
        ans[["BSWF"]] <-
          factor(rep(bswfs.input, times = sapply(l, FUN = nrow)),
                 levels = bswfs.input)
      } else {
        bswf.used <- bswfs.input[1]
      }
    } else {
      bswf.used <- "none"
    }
    setSourceSpct(ans,
                  time.unit = time.unit[1],
                  bswf.used = bswf.used,
                  multiple.wl = mltpl.wl)
    if (!qe.consistent.based) {
      e2q(ans, action = "add", byref = TRUE)
    }
  } else if (l.class == "filter_spct") {
    Tfr.type <- sapply(l, FUN = getTfrType)
    names(Tfr.type) <- NULL
    Tfr.type <- unique(Tfr.type)
    if (length(Tfr.type) > 1L) {
      warning("Inconsistent 'Tfr.type' among filter spectra ",
              "passed to rbindspct")
      return(filter_spct())
    }
    setFilterSpct(ans, Tfr.type = Tfr.type[1], multiple.wl = mltpl.wl)
    filter.descriptor <-
      lapply(l, FUN = getFilterProperties, return.null = TRUE)
    filter.descriptor <- unique(filter.descriptor)
    if (length(filter.descriptor) == 1L) {
      setFilterProperties(ans, filter.descriptor[[1]])
    } else if (length(filter.descriptor) != 0L) {
      message("Discarding heterogeous 'filter.descriptor' attributes!!")
    }
    if (!TA.consistent.based) {
      T2A(ans, action = "add", byref = TRUE)
    }
  } else if (l.class == "reflector_spct") {
    Rfr.type <- sapply(l, FUN = getRfrType)
    names(Rfr.type) <- NULL
    Rfr.type <- unique(Rfr.type)
    if (length(Rfr.type) > 1L) {
      warning("Inconsistent 'Rfr.type' among reflector spectra in rbindspct")
      return(reflector_spct())
    }
    setReflectorSpct(ans, Rfr.type = Rfr.type[1], multiple.wl = mltpl.wl)
    # filter.descriptor <-
    #   lapply(l, FUN = getFilterProperties, return.null = TRUE)
    # filter.descriptor <- unique(filter.descriptor)
    # if (length(filter.descriptor) == 1L) {
    #   setFilterProperties(ans, filter.descriptor[[1]])
    # } else if (length(filter.descriptor) != 0L) {
    #   message("Discarding heterogeous 'filter.descriptor' attributes!!")
    # }
  } else if (l.class == "object_spct") {
    Tfr.type <- sapply(l, FUN = getTfrType)
    names(Tfr.type) <- NULL
    Tfr.type <- unique(Tfr.type)
    Rfr.type <- sapply(l, FUN = getRfrType)
    names(Rfr.type) <- NULL
    Rfr.type <- unique(Rfr.type)
    if (length(Tfr.type) > 1L) {
      warning("Inconsistent 'Tfr.type' among object spectra ",
              "passed to rbindspct")
      return(filter_spct())
    }
    if (length(Rfr.type) > 1L) {
      warning("Inconsistent 'Rfr.type' among object spectra ",
              "passed to rbindspct")
      return(reflector_spct())
    }
    setObjectSpct(ans, Tfr.type = Tfr.type[1], Rfr.type = Rfr.type[1],
                  multiple.wl = mltpl.wl)
    filter.descriptor <-
      lapply(l, FUN = getFilterProperties, return.null = TRUE)
    filter.descriptor <- unique(filter.descriptor)
    if (length(filter.descriptor) == 1L) {
      setFilterProperties(ans, filter.descriptor[[1]])
    } else if (length(filter.descriptor) != 0L) {
      message("Discarding heterogeous 'filter.descriptor' attributes!!")
    }
  } else if (l.class == "solute_spct") {
    K.type <- sapply(l, FUN = getKType)
    names(K.type) <- NULL
    K.type <- unique(K.type)
    if (length(K.type) > 1L) {
      warning("Inconsistent 'K.type' among solute spectra in rbindspct")
      return(solute_spct())
    }
    setSoluteSpct(ans, K.type = K.type, multiple.wl = mltpl.wl)
    solute.descriptor <-
      lapply(l, FUN = getSoluteProperties, return.null = TRUE)
    solute.descriptor <- unique(solute.descriptor)
    if (length(solute.descriptor) == 1L) {
      setSoluteProperties(ans, solute.descriptor[[1]])
    } else if (length(solute.descriptor) != 0L) {
      message("Discarding heterogeous 'solute.descriptor' attributes!!")
    }
  } else if (l.class == "response_spct") {
    time.unit <- sapply(l, FUN = getTimeUnit)
    names(time.unit) <- NULL
    time.unit <- unique(time.unit)
    if (length(time.unit) > 1L) {
      warning("Inconsistent time units among response spectra in rbindspct")
      return(response_spct())
    }
    setResponseSpct(ans, time.unit = time.unit[1], multiple.wl = mltpl.wl)
    if (!qe.consistent.based) {
      e2q(ans, action = "add", byref = TRUE)
    }
    sensor.descriptor <-
      lapply(l, FUN = getSensorProperties, return.null = TRUE)
    sensor.descriptor <- unique(sensor.descriptor)
    if (length(sensor.descriptor) == 1L) {
      setSensorProperties(ans, sensor.descriptor[[1]])
    } else if (length(sensor.descriptor) != 0L) {
      message("Discarding heterogeous 'sensor.descriptor' attributes!!")
    }
  } else if (l.class == "chroma_spct") {
    setChromaSpct(ans, multiple.wl = mltpl.wl)
  } else if (l.class == "cps_spct") {
    setCpsSpct(ans, multiple.wl = mltpl.wl)
  } else if (l.class == "raw_spct") {
    setRawSpct(ans, multiple.wl = mltpl.wl)
  } else if (l.class == "generic_spct") {
    setGenericSpct(ans, multiple.wl = mltpl.wl)
  }
  if (any(scaled.input)) {
    attr(ans, "scaled") <- TRUE
  }
  if (!is.null(comment.ans)) {
    comment(ans) <- comment.ans
  }
  if (is.character(idfactor)) {
    setIdFactor(ans, idfactor)
  }
  setWhenMeasured(ans, when.measured)
  attr(ans, "where.measured") <- where.measured
  # setWhereMeasured(ans, where.measured, simplify = TRUE)
  setWhatMeasured(ans, what.measured)
  setHowMeasured(ans, how.measured)
  attr(ans, "normalized") <- normalized
  if (any(normalized.input)) {
    attr(ans, "normalization") <- normalization
  }
  if (!all(is.na(instr.desc))) {
    setInstrDesc(ans, instr.desc)
  }
  if (!all(is.na(instr.settings))) {
    setInstrSettings(ans, instr.settings)
  }
  ans
}

#' @rdname rbindspct
#'
spctbind <- rbindspct


# wlbind ------------------------------------------------------------------

#' Row bind spectra by wavelength
#'
#' Row bind spectra maintaining wavelength values sorted by interspersing the
#' rows as needed.
#'
#' @details Two objects belonging to the same class, and containing each data
#'   for a single spectrum are row bound keeping wavelengths sorted. Only
#'   columns present in both \code{x} and \code{y} are preserved if \code{fill =
#'   FALSE} and otherwise missing values are filled with \code{NA}.
#'
#' @param x,y generic_spct or of the same derived class, containing each
#'   data for a single spectrum.
#' @param ids character Named vector with the names to use to identify
#'   the origin of the rows.
#' @param strict.wls logical If \code{TRUE} the presence of the same
#'   wavelengths in \code{x} and \code{y} triggers an error. I \code{FALSE}
#'   a message is issued and when a wavelength is both \code{x} and in \code{y},
#'   the row from \code{y} prevails.
#' @inheritParams rbindspct fill
#'
#' @return An object of the same class as \code{x} and \code{y} with data for
#' columns shared by \code{x} and \code{y} based on names.
#'
#' @export
#'
#' @family Methods for row-binding spectra
#'
#' @examples
#' wlbind(sun.spct[-(1:3), ], sun.spct[1:3, ])
#' wlbind(peaks(white_led.source_spct, span = NULL),
#'        wls_at_target(white_led.source_spct),
#'        ids = c(x = "peak", y = "fwhm"))
#'
wlbind <- function(x,
                   y,
                   ids = c(x = "x", y = "y"),
                   strict.wls = FALSE,
                   fill = TRUE) {
  stopifnot(is.any_spct(x) && is.any_spct(y))
  stopifnot(getMultipleWl(x) == 1 && getMultipleWl(y) == 1)
  stopifnot(class_spct(x) == class_spct(y))
  stopifnot(length(ids) == 2 && all(sort(names(ids)) == c("x", "y")))

  shared.wls <- intersect(x[["w.length"]], y[["w.length"]])
  if (length(shared.wls)) {
    if (strict.wls) {
      stop("Same 'w.length' value(s) in 'x' and 'y': ",
           paste(round(shared.wls, 2), collapse = ", "))
    } else {
      x <- x[!x[["w.length"]] %in% shared.wls, ]
      message("Replacing ", length(shared.wls),
              " rows from 'x' with rows from 'y' at ",
              "wavelength(s): ",
              paste(round(shared.wls, 2), collapse = ", "))
    }
  }
  shared.cols <- intersect(colnames(x), colnames(y))
  if (!fill && length(colnames(x)) > length(shared.cols)) {
    x.droped.cols <- setdiff(colnames(x), shared.cols)
    x <- x[ , shared.cols]
    if (!is.na(id_factor(x)) && id_factor(x) %in% x.droped.cols) {
      id_factor(x) <- NULL
    }
    message("Columns dropped from 'x': ", paste(x.droped.cols, collapse = ", "))
  }
  if (!fill && length(colnames(y)) > length(shared.cols)) {
    y.droped.cols <- setdiff(colnames(y), shared.cols)
    y <- y[ , shared.cols]
    if (!is.na(id_factor(y)) && id_factor(y) %in% y.droped.cols) {
      id_factor(y) <- NULL
    }
    message("Columns dropped from 'y': ", paste(y.droped.cols, collapse = ", "))
  }
  x[["id"]] <- ids[["x"]]
  y[["id"]] <- ids[["y"]]
  old.setting <- disable_check_spct()
  z <- rbindspct(list(x, y),
                 fill = fill,
                 idfactor = FALSE,
                 attrs.simplify = TRUE)
  z <- z[order(z[["w.length"]]), ]
  when_measured(z) <- unique(when_measured(z))
  multiple_wl(z) <- 1
  set_check_spct(old.setting)
  check_spct(z)
}
