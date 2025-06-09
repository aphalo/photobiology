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

