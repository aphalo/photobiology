#' Spectrum to rgb color conversion
#'
#' Calculates rgb values from spectra based on human color matching functions
#'
#' @usage s_e_irrad2rgb(w.length, s.e.irrad, sens=ciexyz.data, colour.name=NULL)
#'
#' @param w.length numeric array of wavelengths (nm)
#' @param s.e.irrad numeric array of spectral irradiance values
#' @param sens a dataframe with variables w.length, x, y, and z, giving the chromaticity definition 
#' @param colour.name character string for naming the rgb colour definition
#' 
#' @return A colour defined using \code{rgb()}. The numeric values of the RGB components can be obtained 
#' using function \code{col2rgb()}.
#' 
#' @export
#' @examples
#' data(sun.data)
#' my.color <- with(sun.data, s_e_irrad2rgb(w.length, s.e.irrad, colour.name="sunWhite"))
#' col2rgb(my.color)
#' 
#' @author Pedro J. Aphalo (modified from Chad Eliason's \email{cme16@@zips.uakron.edu} spec2rgb function in package Pavo). 
#' @references CIE(1932). Commission Internationale de l'Eclairage Proceedings, 1931. Cambridge: Cambridge University Press.
#' @references Color matching functions obtained from Colour and Vision Research Laboratory 
#' online data respository at \url{http://www.cvrl.org/}.
#' @references \url{http://www.cs.rit.edu/~ncs/color/t_spectr.html}.

s_e_irrad2rgb <- function(w.length, s.e.irrad, sens=ciexyz.data, colour.name=NULL) {
  low.limit <- min(sens$w.length)
  high.limit <- max(sens$w.length)
  if (length(w.length) == 1) {
    if (w.length < low.limit || w.length > high.limit) {
      return(rgb(0, 0, 0, name=colour.name))
    } else {
      s.e.irrad = 1.0
    }
  } else {
    if (!check_spectrum(w.length, s.e.irrad)) {
      return(NA)
    } else {
      if (min(w.length) > low.limit | max(w.length) < high.limit) {
        warning('Wavelength range does not capture the full cromaticity range\nfilling missing values with zeros.')
      }
    }
  }
  
# will expand and fill with zeros when needed

sens$s.e.irrad <- interpolate_spectrum(w.length, s.e.irrad, sens$w.length, fill=0.0)
sens$s.e.irrad.norm <- with(sens, s.e.irrad / integrate_irradiance(w.length, s.e.irrad))

X <- with(sens, integrate_irradiance(w.length, s.e.irrad.norm * x))
Y <- with(sens, integrate_irradiance(w.length, s.e.irrad.norm * y))
Z <- with(sens, integrate_irradiance(w.length, s.e.irrad.norm * z))

XYZ <- rbind(X, Y, Z) / (X + Y + Z)

xyzmat <- rbind(c(3.240479, -1.537150, -0.498535),
                c(-0.969256, 1.875992, 0.041556),
                c(0.055648, -0.204043, 1.057311))

rgb1 <- xyzmat %*% as.matrix(XYZ)

# print(rgb1)

# normalization
rgb1[rgb1 < 0] <- 0
rgb1[rgb1 > 1] <- 1

rgb.color <- rgb(red=rgb1[1,1], green=rgb1[2,1], blue=rgb1[3,1], name=colour.name)

#class(colrs) <- c('spec2rgb', 'character')

rgb.color

}
