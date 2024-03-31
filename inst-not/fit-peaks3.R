# from https://pisquare.osisoft.com/s/Blog-Detail/a8r1I000000Gvd3QAC/curve-fitting-using-r-peak-asymmetry

# The most commonly used function for fitting chromatography peak is the
# Exponential Modified Gaussian (short EMG). I tried this function but
# unfortunately it didn't work well. The function is numerically unstable for
# very low asymmetries and I got a lot of problems during the optimization. A
# much better choice is the Exponential-Gaussian Hybrid (short EGH):
# http://acadine.physics.jmu.edu/group/technical_notes/GC_peak_fitting/X_lan_Jorgenson.pdf

# The function can be written in R as follows:
# See Lan2001 Lan and Jorgenson https://doi.org/10.1016/S0021-9673(01)00594-5

egh <- function(t, H, tr, sigmag, tau) {
  result <- rep(0, length(t))
  index <- which(2*sigmag^2 + tau*(t-tr)>0)
  result[index] <- H*exp(-(t[index]-tr)^2/(2*sigmag^2+tau*(t[index]-tr)))
  return(result)
}

# where
#
# t : time vector
#
# H : amplitude
#
# tr : peak retention
#
# sigmag : Gaussian standard deviation
#
# tau : exponential time constan
#
# The fitting requires an optimizer and I found the out-of-the-box to be
# sufficient. The code looks as follows:

egh_optim <- function(x, y, par, lower, upper) {
  dat=data.frame(x=x,y=y)
  # fitting routine
  egh_penalty <-function(data, par)(sum((egh(data$x,par[1],par[2],par[3],par[4])-data$y)^2))
  fitting <- optim(par = par, egh_penalty, data = dat,lower=lower, upper=upper,
                  control=list(maxit=1000,trace=0),method="L-BFGS-B")
  # calculate predicted values
  ypred <- egh(x,fitting$par[1],fitting$par[2],fitting$par[3],fitting$par[4])
   # find peak maximum
       H <- fitting$par[1]
       tr <- fitting$par[2]
       sigma <- fitting$par[3]
       tau <- fitting$par[4]
  inverse <-function(x,y)(sum((y-egh(x,H,tr,sigma,tau))^2))

  l10 <-optimize(inverse,c(tr-5*(sigma+tau),tr),tol = 0.0001, y = 0.1*H,maximum = FALSE)
  r10 <-optimize(inverse,c(tr,tr+5*(sigma+tau)),tol = 0.0001, y = 0.1*H,maximum = FALSE)
  asymmetry <- (tr-l10$minimum)/(r10$minimum-tr)

  # collect results and return
  result <- list(par=fitting$par,
                 value=fitting$value,
                 convergence=fitting$convergence,
                 ypred=ypred,
                 peak_max=fitting$par[1],
                 asymmetry=asymmetry)
}


###
# from: http://pd.chem.ucl.ac.uk/pdnn/peaks/loren.htm
#
# Lorentzian (or Cauchy)
#
# The Lorentzian is also a well-used peak function with the form:
#
#   I(2θ) = 	w2 w2 + (2θ − 2θ0)2
#
# where w is equal to half of the peak width (w = 0.5 H). The main features of the Lorentzian function are:
#
#   that it is also easy to calculate
# that, relative to the Gaussian function, it emphasises the tails of the peak
# its integral breadth β = π H / 2
# it has a convenient convolution property (see later)
# it is symmetrical
# instrumental peak shapes are not normally Lorentzian except at high angles where wavelength dispersion is dominant
#
# We note again that since peak intensity is identified with peak area, it is often convenient to also have a form of Lorentz function normalised so that the area is unity; i.e.
#
# L = 	2 / H π 1 + 4 (2θ − 2θ0)2 / H2

# from: https://stat.ethz.ch/pipermail/r-help/2001-October/015687.html

lorentz <- nls( y ~ 1 / (Q*sqrt((1-(x/f0)^2)^2 + ((x/f0) * 1/Q)^2)),
                data = testframe, alg = "plinear",
                start = list (f0=84321.6, Q=600000), trace = TRUE)
