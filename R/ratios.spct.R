#' Calculate photon ratio from spectral irradiance.
#'
#' This function returns the photon ratio for a given pair of wavebands of a
#' light source spectrum.
#'
#' @usage q_ratio(spct, w.band.num=NULL, w.band.denom=NULL, use.cached.mult =
#'   getOption("photobiology.use.cached.mult", default = FALSE),
#'   use.hinges=getOption("photobiology.use.hinges", default=NULL))
#'
#' @param spct an object of class "source.spct"
#' @param w.band.num waveband definition created with new_waveband()
#' @param w.band.denom waveband definition created with new_waveband()
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#'
#' @return a single numeric nondimensional value giving a photon ratio between
#'   pairs of wavebands, with name attribute set to the name of the wavebands
#'   unless a named list of wavebands is supplied in which case the names of the
#'   list elements are used, with "(q:q)" appended.
#' @keywords manip misc
#' @export
#' @examples
#' q_ratio(sun.spct, new_waveband(400,500), new_waveband(400,700))
#'
#' @note Recycling for wavebans takes place when the number of denominator and
#'   denominator wavebands differ. The last two parameters control speed
#'   optimizations. The defaults should be suitable in mosts cases. If you will
#'   use repeatedly the same SWFs on many spectra measured at exactly the same
#'   wavelengths you may obtain some speed up by setting
#'   \code{use.cached.mult=TRUE}. However, be aware that you are responsible for
#'   ensuring that the wavelengths are the same in each call, as the only test
#'   done is for the length of the \code{w.length} vector.
#'
#' @aliases q_ratio
#'
#' @family photon and energy ratio functions
#'
q_ratio <-
  function(spct,
           w.band.num=NULL, w.band.denom=NULL,
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL) ) {
    q.irrad.num <- irrad_spct(spct, w.band=w.band.num,
                              unit.out="photon", quantity="total",
                              wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
                              use.cached.mult=use.cached.mult, use.hinges=use.hinges, allow.scaled=TRUE)
    q.irrad.denom <- irrad_spct(spct, w.band=w.band.denom,
                                unit.out="photon", quantity="total",
                                wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
                                use.cached.mult=use.cached.mult, use.hinges=use.hinges, allow.scaled=TRUE)
    ratio <- q.irrad.num / q.irrad.denom
    names(ratio) <- paste(names(q.irrad.num), ":", names(q.irrad.denom), "(q:q)", sep="")
    setattr(ratio, "time.unit", NULL)
    setattr(ratio, "radiation.unit", "q:q ratio")
    return(ratio)
  }

#' Calculate energy ratio from spectral irradiance.
#'
#' This function returns the photon ratio for a given pair of wavebands of a
#' light source spectrum.
#'
#' @usage e_ratio(spct, w.band.num=NULL, w.band.denom=NULL, use.cached.mult =
#'   getOption("photobiology.use.cached.mult", default = FALSE),
#'   use.hinges=getOption("photobiology.use.hinges", default=NULL))
#'
#' @param spct asource.spct
#' @param w.band.num waveband or list of waveband objects
#' @param w.band.denom waveband or list of waveband objects
#' @param use.cached.mult logical Flag telling whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical Flag telling whether to use hinges to reduce
#'   interpolation errors
#'
#' @return A single numeric nondimensional value giving an energy ratio between
#'   pairs of wavebands, with name attribute set to the name of each waveband
#'   unless a named list of wavebands is supplied in which case the names of the
#'   list elements are used, with "(e:e)" appended.
#'
#' @keywords manip misc
#' @export e_ratio
#' @examples
#' e_ratio(sun.spct, new_waveband(400,500), new_waveband(400,700))
#'
#' @note Recycling for wavebans takes place when the number of denominator and
#'   denominator wavebands differ. The last two parameters control speed
#'   optimizations. The defaults should be suitable in mosts cases. If you will
#'   use repeatedly the same SWFs on many spectra measured at exactly the same
#'   wavelengths you may obtain some speed up by setting
#'   \code{use.cached.mult=TRUE}. However, be aware that you are responsible for
#'   ensuring that the wavelengths are the same in each call, as the only test
#'   done is for the length of the \code{w.length} vector.
#'
#' @aliases e_ratio
#'
#' @family photon and energy ratio functions
#'
e_ratio <-
  function(spct,
           w.band.num=NULL, w.band.denom=NULL,
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL) ) {
    e.irrad.num <- irrad_spct(spct, w.band=w.band.num, unit.out="energy", quantity="total", wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
                              use.cached.mult=use.cached.mult, use.hinges=use.hinges, allow.scaled=TRUE)
    e.irrad.denom <- irrad_spct(spct, w.band=w.band.denom, unit.out="energy", quantity="total", wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
                                use.cached.mult=use.cached.mult, use.hinges=use.hinges, allow.scaled=TRUE)
    ratio <- e.irrad.num / e.irrad.denom
    names(ratio) <- paste(names(e.irrad.num), ":", names(e.irrad.denom), "(e:e)", sep="")
    setattr(ratio, "time.unit", NULL)
    setattr(ratio, "radiation.unit", "e:e ratio")
    return(ratio)
  }

#' Calculate the photon to energy ratio from spectral irradiance.
#'
#' This function returns the photon to energy ratio for each waveband of a light
#' source spectrum.
#'
#' @usage qe_ratio(spct, w.band=NULL, use.cached.mult =
#'   getOption("photobiology.use.cached.mult", default = FALSE),
#'   use.hinges=getOption("photobiology.use.hinges", default=NULL) )
#'
#' @param spct source.spct
#' @param w.band waveband or list of waveband objects
#' @param use.cached.mult logical Flag telling whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical Flag telling whether to use hinges to reduce
#'   interpolation errors
#'
#' @return A vector of \code{numeric} values giving number of moles of photons
#'   per Joule for each waveband, with name attribute set to the name of each
#'   waveband unless a named list of wavebands is supplied in which case the
#'   names of the list elements are used, with "q:e" prepended..
#'
#' @keywords manip misc
#' @export
#' @examples
#' qe_ratio(sun.spct, new_waveband(400,700))
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in mosts cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @family photon and energy ratio functions
#'
qe_ratio <-
  function(spct, w.band=NULL,
           use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL) ) {
    q.irrad <- irrad_spct(spct, w.band=w.band, unit.out="photon",
                          quantity="total",
                          wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
                          use.cached.mult=use.cached.mult, use.hinges=use.hinges, allow.scaled=TRUE)
    e.irrad <- irrad_spct(spct, w.band=w.band, unit.out="energy",
                          quantity="total",
                          wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
                          use.cached.mult=use.cached.mult, use.hinges=use.hinges, allow.scaled=TRUE)
    ratio <- q.irrad / e.irrad
    names(ratio) <- paste("q:e(", names(q.irrad), ")", sep="")
    setattr(ratio, "time.unit", NULL)
    setattr(ratio, "radiation.unit", "q:e ratio")
    return(ratio)
  }

#' Calculate energy to mol photon ratio from spectral irradiance.
#'
#' This function returns the energy to molle photn ratio for each waveband and a
#' light source spectrum.
#'
#' @usage eq_ratio(spct, w.band=NULL, use.cached.mult =
#'   getOption("photobiology.use.cached.mult", default = FALSE),
#'   use.hinges=getOption("photobiology.use.hinges", default=NULL) )
#'
#' @param spct source.spct
#' @param w.band waveband or list of waveband objects
#' @param use.cached.mult logical Flag telling whether multiplier values should
#'   be cached between calls
#' @param use.hinges logical Flag telling whether to use hinges to reduce
#'   interpolation errors
#'
#' @return a numeric value giving number of Joule per mol of photons for each
#'   waveband, with name attribute set to the name of each waveband unless a
#'   named list of wavebands is supplied in which case the names of the list
#'   elements are used, with "e:q" prepended..
#'
#' @keywords manip misc
#' @export
#' @examples
#' eq_ratio(sun.spct, new_waveband(400,700))
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in mosts cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @family photon and energy ratio functions
#'
eq_ratio <-
  function(spct, w.band=NULL, use.cached.mult = getOption("photobiology.use.cached.mult", default = FALSE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL) ) {
    ratio <- 1 / qe_ratio(spct, w.band, use.cached.mult, use.hinges)
    names(ratio) <- gsub("q:e", "e:q", names(ratio), fixed=TRUE )
    setattr(ratio, "time.unit", NULL)
    setattr(ratio, "radiation.unit", "e:q ratio")
    return(ratio)
  }
