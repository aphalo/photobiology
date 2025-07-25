## ----include=FALSE, echo=FALSE------------------------------------------------
knitr::opts_chunk$set(fig.width=7.2, fig.height=4.3)

## ----printing-spectra, eval=TRUE, include=FALSE-------------------------------
# library(tibble)
options(tibble.print_max = 6, tibble.print_min = 4)

## ----pkg-load, eval=TRUE, message = FALSE-------------------------------------
library(photobiology)
library(lubridate)
library(dplyr)

## ----example-1, eval=FALSE----------------------------------------------------
# # not run
# my.spct <- source_spct(w.length = wavelength/10, s.e.irrad = irrad/1000)

## ----query-class-1------------------------------------------------------------
is.any_spct(sun.spct)
is.generic_spct(sun.spct)
is.source_spct(sun.spct)
is.data.frame(sun.spct)

## ----query-class-2------------------------------------------------------------
class_spct(sun.spct)
class(sun.spct)

## ----construction-1-----------------------------------------------------------
my.df <- data.frame(w.length = 400:410, s.e.irrad = 1)
my.spct <- as.source_spct(my.df)
class(my.spct)
class(my.df)
my.spct

## ----construction-2-----------------------------------------------------------
my.g.spct <- as.generic_spct(my.spct)
class(my.g.spct)

## ----construction-3-----------------------------------------------------------
source_spct(w.length = 300:305, s.e.irrad = 1)

## ----construction-4-----------------------------------------------------------
z <- 300:305
y <- 2
source_spct(w.length = z, s.e.irrad = y)

## ----construction-5-----------------------------------------------------------
w.length <- 300:305
s.e.irrad <- 1
source_spct(w.length, s.e.irrad)

## ----construction-6-----------------------------------------------------------
my.d.spct <- as.source_spct(my.df, time.unit = "day")

## ----construction-7-----------------------------------------------------------
source_spct(w.length = 300:305, s.e.irrad = -1)
source_spct(w.length = 300:305, s.e.irrad = -1, strict.range = NULL)

## ----construction-8-----------------------------------------------------------
my.cm.spct <- source_spct(w.length = 300:305, s.e.irrad = 1,
                          comment = "This is a comment")
comment(my.cm.spct)

## ----attr-1-------------------------------------------------------------------
my.spct <- sun.spct
when_measured(my.spct) <-  NULL
when_measured(my.spct)
when_measured(my.spct) <- lubridate::ymd_hms("2015-10-31 22:55:00", tz = "Europe/Helsinki")
when_measured(my.spct)

## ----attr-2-------------------------------------------------------------------
where_measured(my.spct) <- NULL
where_measured(my.spct)
where_measured(my.spct) <- data.frame(lat = 60, lon = -10)
where_measured(my.spct)
where_measured(my.spct) <- data.frame(lat = 60, lon = -10, address = "Somewhere")
where_measured(my.spct)
my.spct

## ----attr-2a------------------------------------------------------------------
what_measured(my.spct) <- "something"
what_measured(my.spct)
my.spct

## -----------------------------------------------------------------------------
sun.spct

## ----attr-3-------------------------------------------------------------------
is_effective(sun.spct)
is_effective(sun.spct * waveband(c(400, 700)))

## ----attr-4-------------------------------------------------------------------
ten.minutes.spct <-
  convertTimeUnit(sun.spct, time.unit = duration(10, "minutes"))
ten.minutes.spct
getTimeUnit(ten.minutes.spct)

## -----------------------------------------------------------------------------
polyester.spct

## -----------------------------------------------------------------------------
convertThickness(polyester.spct, thickness = 2e-3)

## ----col-construction-1-------------------------------------------------------
two_suns.mspct <- source_mspct(list(sun1 = sun.spct, sun2 = sun.spct))
two_suns.mspct

## ----col-construction-2-------------------------------------------------------
mixed.mspct <- generic_mspct(list(filter = polyester.spct, source = sun.spct))
class(mixed.mspct)
lapply(mixed.mspct, class_spct)

## ----col-construction-3-------------------------------------------------------
two_gen.mspct <- as.generic_mspct(two_suns.mspct)
class(two_gen.mspct)
str(two_gen.mspct, max.level = 1, give.attr = FALSE)

## ----col-construction-4-------------------------------------------------------
one_sun.mspct <- as.source_mspct(sun.spct)
class(one_sun.mspct)
str(one_sun.mspct, max.level = 1, give.attr = FALSE)

## ----col-construction-5-------------------------------------------------------
x <- matrix(1:100, ncol = 2)
wl <- 501:550 # wavelengths in nanometres
as.filter_mspct(x, wl, "Tpc")

## ----col-construction-6-------------------------------------------------------
as.filter_mspct(x, wl, "Tpc", spct.names = c("A", "B"))

## ----col-construction-7-------------------------------------------------------
xrow <- matrix(1:100, nrow = 2, byrow = TRUE)
as.filter_mspct(xrow, wl, "Tpc")

## ----col-construction-8-------------------------------------------------------
two_suns.mat <- as.matrix(two_suns.mspct, "s.e.irrad")
class(two_suns.mat)
dim(two_suns.mat)
head(dimnames(two_suns.mat)$spct)
head(dimnames(two_suns.mat)$w.length)
head(attr(two_suns.mat, "w.length"))

## ----col-construction-9-------------------------------------------------------
two_suns.row_mat <- as.matrix(two_suns.mat, "s.e.irrad", byrow = TRUE)
class(two_suns.row_mat)
dim(two_suns.row_mat)
head(dimnames(two_suns.row_mat)$spct)
head(attr(two_suns.row_mat, "w.length"))

## ----bind-1-------------------------------------------------------------------
sun_evening_12.spct <- rbindspct(sun_evening.mspct[1:2])
sun_evening_12.spct

## ----bind-1a------------------------------------------------------------------
subset2mspct(sun_evening_12.spct)

## ----bind-2-------------------------------------------------------------------
test1.df <- data.frame(w.length = rep(200:210, 2),
                       s.e.irrad = rep(c(1, 2), c(11, 11)),
                       spectrum = factor(rep(c("A", "B"), c(11,11))))
subset2mspct(test1.df, member.class = "source_spct", idx.var = "spectrum")

## ----bind-3-------------------------------------------------------------------
subset2mspct(test1.df, member.class = "source_spct", idx.var = "spectrum",
             time.unit = "day")

## ----set-class-1--------------------------------------------------------------
test2.df <- test1.df
setSourceSpct(test2.df, multiple.wl = 2L, idfactor = "spectrum")
getMultipleWl(test2.df)
getIdFactor(test2.df)

## ----set-class-2--------------------------------------------------------------
test3.df <- test1.df
setSourceSpct(test3.df, multiple.wl = NULL, idfactor = "spectrum")
getMultipleWl(test3.df)
getIdFactor(test3.df)

## ----split-1a-----------------------------------------------------------------
test2.df <- data.frame(w.length = 200:210, A = 1, B = 2)
split2source_mspct(x = test2.df)

## ----split-1b-----------------------------------------------------------------
test3.df <- data.frame(w = 200:210, A = 1, B = 2, z = "Z")
split2source_mspct(x = test3.df, w.length.var = "w", idx.var = "z")

## -----------------------------------------------------------------------------
split2source_mspct(test2.df, spct.data.var = "s.q.irrad")

## ----split-2------------------------------------------------------------------
split2source_mspct(test2.df, spct.data.var = "s.q.irrad", time.unit = "day")

## ----join-mspct-01------------------------------------------------------------
my.df <- join_mspct(sun_evening.mspct)
head(my.df)

## ----col-query-class-1--------------------------------------------------------
is.source_mspct(sun_evening.mspct)
class(sun_evening.mspct)

## ----col-query-class-2--------------------------------------------------------
is.filter_mspct(sun_evening.mspct)
is.any_mspct(sun_evening.mspct)
class(sun_evening.mspct)
lapply(sun_evening.mspct, class_spct)
lapply(sun_evening.mspct, class)

## ----extract-1----------------------------------------------------------------
sun_evening.mspct[1]

## ----extract-1a---------------------------------------------------------------
sun_evening.mspct[1:3]

## ----extract-2----------------------------------------------------------------
# warning: this does not swap the names, even if it swaps the spectra
my.mspct <- sun_evening.mspct
summary(my.mspct, which.metadata = "when.measured")$summary[ , -(2:6)]

# member spectr swapped positions, but not the slot names
my.mspct[1:2] <- my.mspct[2:1]
summary(my.mspct, which.metadata = "when.measured")$summary[ , -(2:6)]

# of course, we can also swap the names if needed
names(my.mspct)[1:2] <- names(my.mspct)[2:1]
summary(my.mspct, which.metadata = "when.measured")$summary[ , -(2:6)]

## ----extract-3----------------------------------------------------------------
sun_evening.mspct[[1]]
sun_evening.mspct$time.01
sun_evening.mspct[["time.01"]]

## ----extract-4----------------------------------------------------------------
# local copy
my.mspct <- sun_evening.mspct
names(my.mspct)
# add computed member
my.mspct[["time.01x2"]] <- my.mspct[["time.01"]] * 2
names(my.mspct)
# delete a member
my.mspct[["time.01x2"]] <- NULL
names(my.mspct)

## ----extract-5----------------------------------------------------------------
new.spct <- c(my.mspct[5:4], my.mspct[3])
summary(new.spct, which.metadata = "when.measured")$summary[ , -(2:6)]

## -----------------------------------------------------------------------------
set.seed(1234564)
sampled.mspct <- pull_sample(sun_evening.mspct, size = 2)
summary(sampled.mspct, which.metadata = "when.measured")$summary[ , -(2:6)]

## -----------------------------------------------------------------------------
set.seed(1234564)
sampled.spct <- pull_sample(sun_evening.spct, size = 2)
summary(sampled.spct)

## ----apply-1------------------------------------------------------------------
two.mspct <- sun_evening.mspct[1:2]
msmsply(two.mspct, `+`, 0.1)

## ----apply-2------------------------------------------------------------------
msmsply(two.mspct, trim_wl, range = c(285, 500), fill = NA)

## ----apply-3------------------------------------------------------------------
msdply(two.mspct, wl_max)

## ----apply-4------------------------------------------------------------------
wl_ranges.df <- msdply(two.mspct, wl_range)
wl_ranges.df
cat(comment(wl_ranges.df))

## ----apply-5------------------------------------------------------------------
msdply(two.mspct, wl_range, na.rm = TRUE)

## ----apply-6------------------------------------------------------------------
str(mslply(two.mspct, colnames))

## ----apply-7------------------------------------------------------------------
str(msaply(two.mspct, wl_max))

## ----apply-8------------------------------------------------------------------
msaply(two.mspct, wl_range)

## -----------------------------------------------------------------------------
wl_range(two.mspct)

## -----------------------------------------------------------------------------
s_mean(sun_evening.mspct)

## -----------------------------------------------------------------------------
s_se(sun_evening.mspct)

## -----------------------------------------------------------------------------
s_mean_se(sun_evening.mspct)

## ----convolve-1---------------------------------------------------------------
convolve_each(two.mspct, yellow_gel.spct)

## ----convolve-2---------------------------------------------------------------
convolve_each(yellow_gel.spct, two.mspct)

## ----convolve-3---------------------------------------------------------------
another_two.mspct <- two.mspct
names(another_two.mspct) <- c("a", "b")
convolve_each(another_two.mspct, two.mspct)

## -----------------------------------------------------------------------------
convolve_each(two.mspct, sun.spct, oper = `+`)

## ----col-attr-1---------------------------------------------------------------
when_measured(two.mspct)
when_measured(two.mspct, simplify = TRUE)
when_measured(two.mspct) <- ymd("2015-10-31", tz = "Europe/Helsinki")
when_measured(two.mspct)
when_measured(two.mspct, simplify = TRUE)
when_measured(two.mspct) <- list(ymd_hm("2015-10-31 10:00", tz = "Europe/Helsinki"),
                                 ymd_hm("2015-10-31 11:00", tz = "Europe/Helsinki"))
when_measured(two.mspct) # UTC shown!
when_measured(two.mspct, simplify = TRUE)
two.mspct

## -----------------------------------------------------------------------------
when_measured2tb(sun_evening.mspct)

## -----------------------------------------------------------------------------
when_measured2tb(sun_evening.mspct, col.names = c(when.measured = "acquisition.time"))

## -----------------------------------------------------------------------------
spct_metadata(two.mspct)

## -----------------------------------------------------------------------------
spct_metadata(two.mspct, 
              col.names = c("when.measured" = "time", 
                            "where.measured" = "geocode"),
              unnest = FALSE)

## -----------------------------------------------------------------------------
q_irrad(two.mspct) %>%
  add_attr2tb(two.mspct, 
              col.names = c("lon", "lat", "when.measured"))

## -----------------------------------------------------------------------------
q_irrad(two.mspct) %>%
  add_attr2tb(two.mspct, 
              col.names = c(lon = "longitude", 
                            lat = "latitude", 
                            when.measured = "time"))

## ----wb-1---------------------------------------------------------------------
PAR.wb <- waveband(c(400, 700), wb.name = "PAR")
UVA.wb <- waveband(c(315, 400), wb.name = "UVA")
UVB.wb <- waveband(c(280, 315), wb.name = "UVB")
UVC.wb <- waveband(c(100, 280), wb.name = "UVC")
UV.wb  <- waveband(c(100, 400), wb.name =  "UV")
UV_bands.lst <- list(UVC.wb, UVB.wb, UVA.wb)

## ----wb-2---------------------------------------------------------------------
CIE_e_fun <-
function(w.length){
    CIE.energy <- numeric(length(w.length))
    CIE.energy[w.length <= 298] <- 1
    CIE.energy[(w.length > 298) & (w.length <= 328)] <-
      10^(0.094*(298-w.length[(w.length > 298) & (w.length <= 328)]))
    CIE.energy[(w.length > 328) & (w.length <= 400)] <-
      10^(0.015*(139-w.length[(w.length > 328) & (w.length <= 400)]))
    CIE.energy[w.length > 400] <- 0
    return(CIE.energy)
}

## ----wb-3---------------------------------------------------------------------
CIE.wb <- waveband(c(250, 400), weight = "SWF",
                   SWF.e.fun = CIE_e_fun, SWF.norm = 298)

## ----wb-4---------------------------------------------------------------------
waveband(sun.spct)

## -----------------------------------------------------------------------------
waveband()

## ----wb-5---------------------------------------------------------------------
is.waveband(PAR.wb)

## ----wb-6---------------------------------------------------------------------
is_effective(waveband(c(400,500)))

## ----wb-list-1----------------------------------------------------------------
wavebands <- list(waveband(c(300,400)), waveband(c(400,500)))
wavebands

## ----wb-split-1---------------------------------------------------------------
split_bands(c(200, 225, 300))
split_bands(c(200, 225, 300), length.out = 2)

## ----wb-split-2---------------------------------------------------------------
split_bands(sun.spct, length.out = 2)
split_bands(PAR.wb, length.out = 2)
split_bands(c(200, 800), length.out = 3)

## ----wb-split-3---------------------------------------------------------------
split_bands(list(A = c(200, 300), B = c(400, 500), C = c(250, 350)))
split_bands(list(c(100, 150, 200), c(800, 825)))

## ----wb-split-4---------------------------------------------------------------
split_bands(UV_bands.lst, length.out  =  2)
split_bands(list(c(100, 150, 200), c(800, 825)), length.out = 1)

## ----set-up-printing, eval=FALSE----------------------------------------------
# options(tibble.print_max = 4)
# options(tibble.print_min = 4)

## ----print-1------------------------------------------------------------------
print(sun.spct, n = 3)

## ----print-2------------------------------------------------------------------
str(summary(sun.spct))

## ----print-2a-----------------------------------------------------------------
summary(sun.spct)

## ----print-2b-----------------------------------------------------------------
summary(sun_evening.mspct)

## ----print-2c, eval = FALSE---------------------------------------------------
# summary(sun_evening.mspct, expand = "each")

## ----print-2d-----------------------------------------------------------------
summary(sun_evening.spct)

## ----print-2e-----------------------------------------------------------------
summary(sun_evening.spct, expand = "collection")

## -----------------------------------------------------------------------------
na.omit(sun.spct)
na.exclude(sun.spct)

## ----bin-oper-1---------------------------------------------------------------
sun.spct * sun.spct

## ----bin-oper-2---------------------------------------------------------------
sun.spct * polyester.spct

## ----bin-oper-3---------------------------------------------------------------
sun.spct * polyester.spct * polyester.spct
sun.spct * polyester.spct^2

## ----bin-oper-4---------------------------------------------------------------
sun.spct * 2
2 * sun.spct
sun.spct * c(0,1)

## ----bin-oper-5---------------------------------------------------------------
sun.spct * UVB.wb

## ----bin-oper-6, eval=FALSE---------------------------------------------------
# sun.spct * CIE.wb

## ----unary-oper-1-------------------------------------------------------------
-sun.spct
sqrt(sun.spct)

## ----options-1----------------------------------------------------------------
options(photobiology.radiation.unit = "photon")
sun.spct * UVB.wb
options(photobiology.radiation.unit = "energy")
sun.spct * UVB.wb

## ----options-2----------------------------------------------------------------
photon_as_default()
sun.spct * UVB.wb
energy_as_default()
sun.spct * UVB.wb
unset_radiation_unit_default()

## ----options-3----------------------------------------------------------------
using_photon(sun.spct * UVB.wb)
using_energy(sun.spct * UVB.wb)

## ----manip-1------------------------------------------------------------------
# STOPGAP
shade.spct <- sun.spct

## ----manip-2------------------------------------------------------------------
rbindspct(list(sun.spct, shade.spct))
rbindspct(list(A = sun.spct, B = shade.spct), idfactor = "site")

## ----manip-3------------------------------------------------------------------
sun.spct[1:10, ]
sun.spct[1:10, 1]
sun.spct[1:10, 1, drop = TRUE]
sun.spct[1:10, "w.length", drop = TRUE]

## ----manip-4------------------------------------------------------------------
subset(sun.spct, s.e.irrad > 0.2)
subset(sun.spct, w.length > 600)
subset(sun.spct, c(TRUE, rep(FALSE, 99)))

## ----manip-5------------------------------------------------------------------
e2q(sun.spct, "add")
e2q(sun.spct, "replace")

## -----------------------------------------------------------------------------
polyester.spct

## -----------------------------------------------------------------------------
any2Afr(polyester.spct, "add")

## -----------------------------------------------------------------------------
any2Afr(polyester.spct, "replace")

## ----manip-6------------------------------------------------------------------
normalize(sun.spct)

## ----manip-7------------------------------------------------------------------
normalize(sun.spct, range = PAR.wb, norm = "max")

## ----manip-8------------------------------------------------------------------
normalize(sun.spct, norm = 600.3)

## -----------------------------------------------------------------------------
my.spct <- normalize(sun.spct)
is_normalized(my.spct)
getNormalized(my.spct)

## ----manip-9------------------------------------------------------------------
fscale(sun.spct)
fscale(sun.spct, set.scaled = FALSE)
fscale(sun.spct, target = 100)
fscale(sun.spct, target = 100, set.scaled = TRUE)

## ----manip-9a-----------------------------------------------------------------
fscale(sun.spct, f = "integral")
fscale(sun.spct, range = PAR.wb, f = e_irrad)
fscale(sun.spct, range = PAR.wb, f = q_irrad, target = 800e-6)

## ----manip-9b-----------------------------------------------------------------
my.spct <- fscale(sun.spct)
is_scaled(my.spct)
getScaled(my.spct)

## ----manip-10-----------------------------------------------------------------
fshift(white_led.source_spct, range = UVB.wb, f = "mean")
fshift(sun.spct, range = c(280,290), f = "min")

## ----manip-11-----------------------------------------------------------------
clean(polyester.spct - 0.053)

## ----manip-11a----------------------------------------------------------------
clean(sun.spct - 0.01, range = c(280.5, 282), fill = NA)

## ----manip-12-----------------------------------------------------------------
spikes(sun.spct)

## ----manip-12a----------------------------------------------------------------
my_sun.spct <- despike(sun.spct)
spikes(my_sun.spct)

## -----------------------------------------------------------------------------
smooth_spct(sun.spct)

## -----------------------------------------------------------------------------
smooth_spct(polyester.spct, method = "supsmu", strength = 2)

## ----manip-13-----------------------------------------------------------------
interpolate_wl(sun.spct, seq(400, 500, by = 0.1))

## ----trim-1-------------------------------------------------------------------
clip_wl(sun.spct, range = c(400, 402))
clip_wl(sun.spct, range = c(400, NA))

## ----trim-2-------------------------------------------------------------------
clip_wl(sun.spct, range = UVA.wb)

## ----trim-3-------------------------------------------------------------------
clip_wl(sun.spct, range = c(100, 200))

## ----trim-4-------------------------------------------------------------------
trim_wl(sun.spct, c(282.5, NA))
clip_wl(sun.spct, c(282.5, NA))

## ----trim-5-------------------------------------------------------------------
trim_wl(sun.spct, PAR.wb)

## ----trim-6-------------------------------------------------------------------
trim_wl(sun.spct, c(281.5, NA), fill = NA)

## ----trim-7-------------------------------------------------------------------
trim_wl(sun.spct, c(275, NA), fill = 0)

## ----trim-8-------------------------------------------------------------------
trim_wl(sun.spct, c(281.5, NA), fill = NA)
trim_wl(sun.spct, c(281.5, NA), fill = NA, use.hinges = FALSE)

## -----------------------------------------------------------------------------
trim2overlap(two.mspct)

## -----------------------------------------------------------------------------
extend2extremes(two.mspct, fill = 0)

## -----------------------------------------------------------------------------
nrow(yellow_gel.spct)
wl_stepsize(yellow_gel.spct)
thinned.spct <- thin_wl(yellow_gel.spct)
nrow(thinned.spct)
wl_stepsize(thinned.spct)

## ----weights-1----------------------------------------------------------------
sun.spct * CIE.wb

## ----tag-1--------------------------------------------------------------------
tag(sun.spct, PAR.wb)
tag(sun.spct, UV_bands.lst)

## ----tag-2--------------------------------------------------------------------
tg.sun.spct <- tag(sun.spct, PAR.wb)
attr(tg.sun.spct, "spct.tags")

## ----tag-3--------------------------------------------------------------------
wb2tagged_spct(UV_bands.lst)
wb2rect_spct(UV_bands.lst)

## ----tag-4--------------------------------------------------------------------
tg.sun.spct
is_tagged(tg.sun.spct)
untg.sun.spct <- untag(tg.sun.spct)
is_tagged(untg.sun.spct)

## -----------------------------------------------------------------------------
is_tagged(untg.sun.spct)
untag(tg.sun.spct, byref = TRUE)
is_tagged(untg.sun.spct)


## ----summary-1----------------------------------------------------------------
summary(sun.spct)

## ----summary-2----------------------------------------------------------------
summary(two_filters.mspct)

## ----summary-2-each-----------------------------------------------------------
summary(two_filters.mspct, expand = "each")

## ----summary-2a---------------------------------------------------------------
summary(two_filters.spct)

## ----summary-2a-collection, eval=FALSE----------------------------------------
# summary(two_filters.spct, expand = "collection")

## ----summary-2a-each, eval=FALSE----------------------------------------------
# summary(two_filters.spct, expand = "each")

## ----summary-3----------------------------------------------------------------
wl_range(sun.spct)
wl_min(sun.spct)
wl_max(sun.spct)
wl_midpoint(sun.spct)
wl_expanse(sun.spct)
wl_stepsize(sun.spct)

## ----summary-4----------------------------------------------------------------
filters.mspct <- filter_mspct(list(none = clear.spct,
                                   pet = polyester.spct,
                                   yellow = yellow_gel.spct))
wl_range(filters.mspct)

## ----summary-5----------------------------------------------------------------
nrow(sun.spct)
nrow(peaks(sun.spct))
nrow(valleys(sun.spct))

## ----summary-5a---------------------------------------------------------------
peaks(sun.spct, span = 51)
valleys(sun.spct, span = 51)

## ----summary-8----------------------------------------------------------------
peaks(sun.spct, span = NULL)

## ----summary-6a---------------------------------------------------------------
peaks(sun.spct, 
      span = NULL,
      refine.wl = TRUE)

## ----summary-6----------------------------------------------------------------
peaks(sun.spct, 
      span = NULL, 
      unit.out = "energy")

peaks(sun.spct, 
      span = NULL, 
      unit.out = "photon")

## ----summary-10---------------------------------------------------------------
spikes(sun.spct)

## ----summary-10a--------------------------------------------------------------
spikes(sun.spct, z.threshold = 6)

## ----col-summary-1------------------------------------------------------------
peaks(sun_evening.mspct, span = NULL)

## ----col-summary-1a-----------------------------------------------------------
peaks(sun_evening.spct, span = NULL)

## ----find-wls-1---------------------------------------------------------------
wls_at_target(Ler_leaf_trns.spct, 
              target = "half.maximum")
wls_at_target(Ler_leaf_trns.spct, 
              target = "half.maximum", 
              interpolate = TRUE)

## ----find-wls-2---------------------------------------------------------------
wls_at_target(filters.mspct, target = "half.maximum")

## ----irrad-1------------------------------------------------------------------
e_irrad(sun.spct)

## ----irrad-2------------------------------------------------------------------
e_irrad(sun.spct, PAR.wb)

## ----irrad-3------------------------------------------------------------------
e_irrad(sun.spct, c(400, 700))

## ----irrad-4a-----------------------------------------------------------------
q_irrad(sun.spct, PAR.wb, scale.factor = 1e6) # umol s-1 m-2

## ----irrad-5------------------------------------------------------------------
q_irrad(white_led.source_spct, PAR.wb, time.unit = "hour")

## ----irrad-6------------------------------------------------------------------
e_irrad(sun_daily.spct, PAR.wb, time.unit = "second")

## ----irrad-7------------------------------------------------------------------
e_irrad(sun.spct, UV_bands.lst) # W m-2

## ----irrad-7a-----------------------------------------------------------------
q_irrad(sun.spct, UV_bands.lst) # mol s-1 m-2

## ----irrad-7b-----------------------------------------------------------------
q_irrad(sun.spct, UV_bands.lst, scale.factor = 1e6) # umol s-1 m-2

## ----irrad-8a-----------------------------------------------------------------
e_irrad(sun.spct, UV_bands.lst, quantity = "total") # watt m-2

## ----irrad-8b-----------------------------------------------------------------
e_irrad(sun.spct, UV_bands.lst, quantity = "average") # watt m-2 nm-1

## ----irrad-8c-----------------------------------------------------------------
e_irrad(sun.spct, UV_bands.lst, quantity = "contribution")

## ----irrad-8cc----------------------------------------------------------------
e_irrad(sun.spct, UV_bands.lst, quantity = "relative")

## ----irrad-8d-----------------------------------------------------------------
e_irrad(sun.spct, UV_bands.lst, quantity = "contribution.pc")

## ----irrad-8dd----------------------------------------------------------------
e_irrad(sun.spct, UV_bands.lst, quantity = "relative.pc")

## -----------------------------------------------------------------------------
e_irrad(sun.spct, PAR.wb, time.unit = duration(8, "hours"))

## -----------------------------------------------------------------------------
e_fluence(sun.spct, PAR.wb, exposure.time = duration(8, "hours"))

## -----------------------------------------------------------------------------
q_irrad(sun.spct, UV_bands.lst, naming = "short")

## -----------------------------------------------------------------------------
q_irrad(sun.spct, UV_bands.lst, naming = "none")

## -----------------------------------------------------------------------------
names(UV_bands.lst) <- c("UV-C", "UV-B", "UV-A")
q_irrad(sun.spct, UV_bands.lst, naming = "short", scale.factor = 1e6)

## -----------------------------------------------------------------------------
e_irrad(sun_evening.mspct, w.band = PAR.wb)

## -----------------------------------------------------------------------------
q_irrad(sun_evening.mspct, 
        w.band = PAR.wb,
        scale.factor = 1e6, # umol m-2 s-1
        attr2tb = c(when.measured = "time", lon = "lon", lat = "lat"))

## ----col-convolve-1-----------------------------------------------------------
filtered_sun <- convolve_each(filters.mspct, sun.spct)
q_irrad(filtered_sun,
        list(UVA.wb, PAR.wb),
        scale.factor = 1e6,
        idx = "Filter")

## ----col-convolve-2-----------------------------------------------------------
q_irrad(convolve_each(filters.mspct, sun.spct), 
        list("UV-A" = UVA.wb, PAR.wb),
        scale.factor = 1e6,  # umol m-2 s-1
        naming = "short",
        attr2tb = c(what.measured = "Filter type"))[ , c(4, 2, 3)]

## ----fluence-1----------------------------------------------------------------
fluence(sun.spct, exposure.time = duration(1, "hours"))

## ----fluence-1a---------------------------------------------------------------
fluence(sun.spct, exposure.time = 3600) # seconds

## ----fluence-2----------------------------------------------------------------
q_fluence(sun.spct, PAR.wb, exposure.time = duration(25, "minutes"))

## -----------------------------------------------------------------------------
q_ratio(sun.spct, UVB.wb, PAR.wb)

## -----------------------------------------------------------------------------
q_ratio(sun.spct, list(UVC.wb, UVB.wb, UVA.wb))

## -----------------------------------------------------------------------------
q_ratio(sun.spct, list(UVC.wb, UVB.wb, UVA.wb), PAR.wb)

## ----ratios-2-----------------------------------------------------------------
qe_ratio(sun.spct,
         list("UV-B" = UVB.wb, PAR.wb), 
         scale.factor = 1e6,
         name.tag = " (umol/J)")

## ----ratios-3-----------------------------------------------------------------
q_ratio(filtered_sun, 
        list(UVB.wb, UVA.wb),
        PAR.wb)

## ----ratios-3a----------------------------------------------------------------
q_ratio(filtered_sun, 
        list(UVB.wb, UVA.wb),
        PAR.wb, 
        scale.factor = 100, 
        name.tag = " (% photons)", 
        idx = "Filter")

## -----------------------------------------------------------------------------
transmittance(polyester.spct, list(UVB.wb, UVA.wb, PAR.wb))

## -----------------------------------------------------------------------------
transmittance(polyester.spct, 
              list(UVB.wb, UVA.wb, PAR.wb),
              naming = "none")

## -----------------------------------------------------------------------------
transmittance(polyester.spct, 
              list(UVB.wb, UVA.wb, PAR.wb),
              naming = "short")

## -----------------------------------------------------------------------------
reflectance(green_leaf.spct, waveband(c(600, 700)))

## -----------------------------------------------------------------------------
transmittance(filters.mspct, 
              w.band = list(UVA.wb, PAR.wb))

## -----------------------------------------------------------------------------
transmittance(filters.mspct, 
              w.band = list(UVA.wb, PAR.wb),
              naming = "short")

## -----------------------------------------------------------------------------
transmittance(filters.mspct, 
              w.band = list("UV-A" = UVA.wb, PAR = PAR.wb))

## -----------------------------------------------------------------------------
transmittance(filters.mspct, 
              w.band = list(UVA = UVA.wb, PAR = PAR.wb),
              naming = "short")

## -----------------------------------------------------------------------------
transmittance(filters.mspct, 
              w.band = UVA.wb,
              naming = "short",
              attr2tb = c("what.measured" = "Filter type"))

## -----------------------------------------------------------------------------
transmittance(filters.mspct,
              w.band = UVA.wb,
              attr2tb = "what.measured")[ , c(3, 2)]

## -----------------------------------------------------------------------------
normalized_diff_ind(Ler_leaf_rflt.spct,
                    waveband(c(740, 840)), waveband(c(590, 690)),
                    reflectance)

## -----------------------------------------------------------------------------
normalized_diff_ind(sun.spct,
                    waveband(c(600, 700)), waveband(c(400, 500)),
                    q_irrad)

## -----------------------------------------------------------------------------
response(photodiode.spct)

## -----------------------------------------------------------------------------
e_response(photodiode.spct, list(UVB.wb, UVA.wb))

## -----------------------------------------------------------------------------
sensors <- response_mspct(list(GaAsP = photodiode.spct,
                               CCD = ccd.spct))
response(sensors, list(UVB.wb, UVA.wb, PAR.wb), quantity = "contribution")

## -----------------------------------------------------------------------------
integrate_spct(sun.spct)

## -----------------------------------------------------------------------------
average_spct(sun.spct)

## -----------------------------------------------------------------------------
compare_spct(source_mspct(list(sun1 = sun.spct, sun2 = sun.spct * 2)))

## -----------------------------------------------------------------------------
compare_spct(filter_mspct(list(pet = polyester.spct,
                              yllw = yellow_gel.spct)),
             w.band = 50,
            .comparison.fun = `<`)

## -----------------------------------------------------------------------------
illuminance(sun.spct)

## -----------------------------------------------------------------------------
illuminance(sun.spct, std = "CIE10deg")

## -----------------------------------------------------------------------------
illuminance(sun_daily.spct)

## -----------------------------------------------------------------------------
color_of(550) # green
color_of(630) # red
color_of(c(550, 630, 380, 750)) # vectorized

## -----------------------------------------------------------------------------
color_of(sun.spct)
color_of(sun.spct * yellow_gel.spct)

## -----------------------------------------------------------------------------
color_of(waveband(c(400, 500), wb.name = "my_BL"))

## -----------------------------------------------------------------------------
color_of(sun.spct, chroma.type = "CC")
color_of(sun.spct, chroma.type = "CMF")
color_of(sun.spct, chroma.type = beesxyzCMF.spct)

