## ---- include=FALSE, echo=FALSE------------------------------------------
library(knitr)
opts_chunk$set(fig.path = 'figure/pos-', fig.align = 'center', fig.show = 'hold', fig.width = 7, fig.height = 4)

## ---- printing-spectra, eval=TRUE, include=FALSE-------------------------
library(tibble)
options(tibble.print_max = 6, tibble.print_min = 4)

## ---- pkg-load, eval=TRUE------------------------------------------------
library(photobiology)
library(lubridate)

## ---- example-1, eval=FALSE----------------------------------------------
#  # not run
#  my.spct <- source_spct(w.length = wavelength/10, s.e.irrad = irrad/1000)

## ---- query-class-1------------------------------------------------------
is.any_spct(sun.spct)
is.generic_spct(sun.spct)
is.source_spct(sun.spct)

## ---- query-class-2------------------------------------------------------
class_spct(sun.spct)
class(sun.spct)

## ---- construction-1-----------------------------------------------------
my.df <- data.frame(w.length = 400:410, s.e.irrad = 1)
my.spct <- as.source_spct(my.df)
class(my.spct)
class(my.df)
my.spct

## ---- construction-2-----------------------------------------------------
my.g.spct <- as.generic_spct(my.spct)
class(my.g.spct)

## ---- construction-3-----------------------------------------------------
source_spct(w.length = 300:305, s.e.irrad = 1)

## ---- construction-4-----------------------------------------------------
z <- 300:305
y <- 2
source_spct(w.length = z, s.e.irrad = y)

## ---- construction-5-----------------------------------------------------
w.length <- 300:305
s.e.irrad <- 1
source_spct(w.length, s.e.irrad)

## ---- construction-6-----------------------------------------------------
my.d.spct <- as.source_spct(my.df, time.unit = "day")

## ---- construction-7-----------------------------------------------------
source_spct(w.length = 300:305, s.e.irrad = -1)
source_spct(w.length = 300:305, s.e.irrad = -1, strict.range = NULL)

## ---- construction-8-----------------------------------------------------
my.cm.spct <- source_spct(w.length = 300:305, s.e.irrad = 1,
                          comment = "This is a comment")
comment(my.cm.spct)

## ---- attr-1-------------------------------------------------------------
my.spct <- sun.spct
setWhenMeasured(my.spct, NULL)
getWhenMeasured(my.spct)
setWhenMeasured(my.spct, ymd_hms("2015-10-31 22:55:00", tz = "EET"))
getWhenMeasured(my.spct)

## ---- attr-2-------------------------------------------------------------
setWhereMeasured(my.spct, NULL)
getWhereMeasured(my.spct)
setWhereMeasured(my.spct, lat = 60, lon = -10)
getWhereMeasured(my.spct)
getWhereMeasured(my.spct)$lon
my.spct

## ---- attr-3-------------------------------------------------------------
is_effective(sun.spct)
is_effective(sun.spct * waveband(c(400, 700)))

## ---- attr-4-------------------------------------------------------------
ten.minutes.spct <-
  convertTimeUnit(sun.spct, time.unit = duration(10, "minutes"))
ten.minutes.spct
getTimeUnit(ten.minutes.spct)

## ---- col-construction-1-------------------------------------------------
two_suns.mspct <- source_mspct(list(sun1 = sun.spct, sun2 = sun.spct))
two_suns.mspct

## ---- col-construction-2-------------------------------------------------
mixed.mspct <- generic_mspct(list(filter = clear.spct, source = sun.spct))
class(mixed.mspct)
lapply(mixed.mspct, class_spct)

## ---- col-construction-3-------------------------------------------------
two_gen.mscpt <- as.generic_mspct(two_suns.mspct)
class(two_gen.mscpt)
lapply(two_gen.mscpt, class_spct)

## ---- bind-1-------------------------------------------------------------
two_suns.spct <- rbindspct(list(a = sun.spct, b = sun.spct / 2))
subset2mspct(two_suns.spct)

## ---- bind-2-------------------------------------------------------------
test1.df <- data.frame(w.length = rep(200:210, 2),
                       s.e.irrad = rep(c(1, 2), c(11, 11)),
                       spectrum = factor(rep(c("A", "B"), c(11,11))))
subset2mspct(test1.df, member.class = "source_spct", idx.var = "spectrum")

## ---- set-class-1--------------------------------------------------------
setSourceSpct(test1.df, multiple.wl = 2L)
test1.df

## ---- split-1------------------------------------------------------------
test2.df <- data.frame(w.length = 200:210, A = 1, B = 2, z = "A")
split2source_mspct(test2.df)
split2source_mspct(test2.df, spct.data.var = "s.q.irrad")

## ---- col-query-class-1--------------------------------------------------
is.source_mspct(two_suns.mspct)
class(two_suns.mspct)

## ---- col-query-class-2--------------------------------------------------
is.filter_mspct(mixed.mspct)
is.any_mspct(mixed.mspct)
class(mixed.mspct)
lapply(mixed.mspct, class_spct)
lapply(mixed.mspct, class)

