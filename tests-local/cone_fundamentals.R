library(photobiology)
library(ggspectra)
library(photobiologyWavebands)

ggplot(cone_fundamentals10.spct) +
  geom_line(aes(w.length, z, colour = "blue")) +
  geom_line(aes(w.length, y, colour = "green")) +
  geom_line(aes(w.length, x, colour = "red")) +
  scale_color_manual(values = c(green = "green", blue = "blue", red = "red")) +
  theme_bw()

autoplot(, 
         w.band = VIS_bands(), unit.out = "photon",
         annotations = c("-", "summaries", "peaks")) +
  theme_bw()
