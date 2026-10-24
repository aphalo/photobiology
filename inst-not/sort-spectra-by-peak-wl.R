library(photobiologyLEDs)
library(ggspectra)
library(dplyr)

peaks(leds.mspct, span = NULL) |> rbind_spct() |> arrange(w.length) -> z

