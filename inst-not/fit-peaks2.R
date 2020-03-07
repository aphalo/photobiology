library(photobiology)
library(ggspectra)
library(RcppFaddeeva)
library(photobiologyLamps)

spct <- white_led.source_spct
autoplot(spct)

fit <- nls(s.e.irrad ~ alpha*Voigt(w.length, x0, sigma, gamma),
           data = spct,
           start = list(alpha = 1, x0 = 604, sigma = 100, gamma = 1))

summary(fit)

pre.spct <- source_spct(w.length = spct$w.length, s.e.irrad = predict(fit))

autoplot(pre.spct) +
  geom_line(data = spct, linetype = "dotted")

spct1 <- clip_wl(spct, c(500,700))
autoplot(spct1)

fit <- nls(s.e.irrad ~ alpha*Voigt(w.length, x0, sigma, gamma),
           data = spct1,
           start = list(alpha = 100, x0 = 603.6, sigma = 30, gamma = 30))

pre.spct <- source_spct(w.length = spct1$w.length, s.e.irrad = predict(fit))

autoplot(pre.spct) +
  geom_line(data = spct1, linetype = "dotted")

f <- splinefun(spct$w.length, spct$s.e.irrad)

f(spct$w.length)

optimize(f, interval=c(590, 615), maximum=TRUE)

wls <- seq(590, 615, by = 0.005)
wls[which.max(f(wls))]

wls <- seq(450, 465, by = 0.005)
wls[which.max(f(wls))]

spct <- lamps.mspct$philips.tl12
autoplot(spct)

f <- splinefun(spct$w.length, spct$s.e.irrad, method = "natural")
pre.spct <- source_spct(w.length = seq(300, 325, by = 0.005),
                        s.e.irrad = f(seq(300, 325, by = 0.005)))

ggplot(spct) +
  geom_point() +
  geom_line(data = pre.spct, linetype = "dotted") +
  xlim(300, 325)

wls <- seq(300, 325, by = 0.005)
wls[which.max(f(wls))]

wls <- seq(360, 380, by = 0.005)
wls[which.max(f(wls))]

autoplot(spct)

