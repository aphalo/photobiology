
# illuminace --------------------------------------------------------------

#' Irradiance
#'
#' Computes illuminance (lux), or the luminous flux incident on a surface, from
#' spectral irradiance stored in a \code{source_spct} object.
#'
#' @param spct an R object.
#' @param std character The luminous efficiency function to use,
#'   \code{"CIE2deg"} or \code{"CIE10deg"}.
#' @param scale.factor numeric vector of length 1, or the character string
#'   \code{exposure}.
#' @param allow.scaled logical indicating whether scaled or normalized spectra
#'   as argument to spct are flagged as an error.
#' @param naming character one of \code{"long"}, \code{"default"},
#'   \code{"short"} or \code{"none"}. Used to select the type of names to assign
#'   to returned value.
#' @param ... other arguments (possibly ignored)
#'
#' @note Formal parameter \code{allow.scaled} is used internally for calculation
#'   of ratios, as rescaling and normalization do not invalidate the calculation
#'   of ratios within one spectrum.
#'
#' @return A named \code{numeric} vector of length one in the case of methods
#'   for individual spectra. A \code{data.frame} in the case of collections of
#'   spectra, containing one column with illuminance, an index column with the
#'   names of the spectra, and optionally additional columns with metadata
#'   values retrieved from the attributes of the member spectra.
#'
#'   The \code{time.unit} attribute is always second. Units are as follows: if
#'   time.unit of the argument passed to \code{spct} is second, [W m-2 nm-1] ->
#'   [lx], otherwise average value [lx] for the period unless
#'   \code{exposure = TRUE}.
#'
#' @export
#' @examples
#' illuminance(sun.spct)
#' illuminance(sun.daily.spct)
#' illuminance(sun.daily.spct, scale.factor = "exposure")
#' illuminance(sun.daily.spct, scale.factor = 1e-3)
#'
#' @references
#' Stockman, A. (2019) Cone fundamentals and CIE standards.
#' \emph{Current Opinion in Behavioral Sciences}, 30, 87-93.
#' \doi{10.1016/j.cobeha.2019.06.005}
#'
#' @family illumination functions
#'
illuminance <- function(spct, std, scale.factor, allow.scaled, ...) UseMethod("illuminance")

#' @describeIn illuminance Default for generic function
#'
#' @export
#'
illuminance.default <- function(spct, std, scale.factor, allow.scaled, ...) {
  warning("'illuminance' is not defined for objects of class ", class(spct)[1])
  return(NA_real_)
}

#' @describeIn illuminance  Calculates illuminance from a \code{source_spct}
#'   object.
#'
#' @method illuminance source_spct
#' @export
#'
illuminance.source_spct <-
  function(spct,
           std = "CIE2deg",
           scale.factor = 1,
           allow.scaled = FALSE,
           naming = "default",
           ...) {

    # we look for multiple spectra in long form
    if (getMultipleWl(spct) > 1) {
      # convert to a collection of spectra
      mspct <- subset2mspct(x = spct,
                            idx.var = getIdFactor(spct),
                            drop.idx = FALSE)
      # call method on the collection
      return(illuminance(spct = mspct,
                         std = std,
                         scale.factor = scale.factor,
                         allow.scaled = allow.scaled,
                         naming = naming,
                         ...))
    }

    if (!allow.scaled && is_normalized(spct)) {
      warning("The spectral data have been normalized, ",
              "preventing calculation of illuminance. ",
              "See 'setNormalised()' and 'normalise()'.")
      return(NA_real_)
    }
    if (!allow.scaled && is_scaled(spct)) {
      warning("The spectral data have been scaled, ",
              "preventing calculation of illuminance. ",
              "See 'setScaled()' and 'fscale()'.")
      return(NA_real_)
    }

    data.time.unit <- getTimeUnit(spct, force.duration = TRUE)

    value.name <- "Ev[lx]"
    if (data.time.unit != "second" & is.numeric(scale.factor)) {
      spct <- convertTimeUnit(spct, time.unit = "second")
      if (scale.factor == 1) {
        unit.attr <- "mean illuminance [lx];"
      } else {
        unit.attr <- "scaled mean illuminance [lx];"
      }
    } else if (is.numeric(scale.factor)) {
      unit.attr <- "scaled illuminance [lx];"
    } else if (is.character(scale.factor) && scale.factor == "exposure") {
      value.name <- "Hv[lx]"
      unit.attr <- "luminous exposure [lx];"
      scale.factor <- 1
    }

    std.label <- std
    if (wl_min(spct) > 390) { # limit of CIE luminous efficiency
      std.label <- paste("]", std.label, sep = "")
    }
    if (wl_max(spct) < 760) { # limit of VIS
      std.label <- paste(std.label, "[", sep = "")
    }
    if (std == "CIE2deg" || std == "CIEV2") {
      z <- response(using_energy(spct * photobiology::ciev2.spct), unit.out = "energy") * 683.002
    } else if (std == "CIE10deg" || std == "CIEV10") {
      z <- response(using_energy(spct * photobiology::ciev10.spct), unit.out = "energy") * 683.002
    } else {
      stop("Unkown standard:", std)
    }
    if (std.label != std) {
      warning("Bad value: spectral data does not cover CIE spectrum")
    }
    attr(z, "radiation.unit") <- paste(unit.attr, std.label)
    names(z) <- value.name
    z * scale.factor
  }

# source_mspct methods -----------------------------------------------

#' @describeIn illuminance Calculates illuminance from a \code{source_mspct} object.
#'
#' @param attr2tb character vector, see \code{\link{add_attr2tb}} for the syntax
#'   for \code{attr2tb} passed as is to formal parameter \code{col.names}.
#' @param idx character Name of the column with the names of the members of the
#'   collection of spectra.
#' @param .parallel	if TRUE, apply function in parallel, using parallel backend
#'   provided by foreach.
#' @param .paropts a list of additional options passed into the foreach function
#'   when parallel computation is enabled. This is important if (for example)
#'   your code relies on external data or packages: use the .export and
#'   .packages arguments to supply them so that all cluster nodes have the
#'   correct environment set up for computing.
#'
#' @export
#'
illuminance.source_mspct <-
  function(spct,
           std = "CIE2deg",
           scale.factor = 1,
           allow.scaled = FALSE,
           naming = "default",
           ...,
           attr2tb = NULL,
           idx = "spct.idx",
           .parallel = FALSE,
           .paropts = NULL) {

    spct <- subset2mspct(spct) # expand long form spectra within collection

    z <-
      msdply(
        mspct = spct,
        .fun = illuminance,
        std = std,
        scale.factor = scale.factor,
        allow.scaled = allow.scaled,
        naming = naming,
        idx = idx,
        col.names = "illuminace",
        .parallel = .parallel,
        .paropts = .paropts
      )
    add_attr2tb(tb = z,
                mspct = spct,
                col.names = attr2tb,
                idx = idx)
  }
