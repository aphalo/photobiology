library(ggplot2)
library(photobiology)
library(ggspectra)

my.df <- data.frame(w.length = 300:700, s.e.irrad = 1)
my.spct <- as.source_spct(my.df)
class(my.spct)

ggplot(my.spct)+geom_line() +
  stat_wl_strip(ymin = 0.9, ymax = 1.1) +
  scale_fill_identity() +
  expand_limits(y = 1.2)
