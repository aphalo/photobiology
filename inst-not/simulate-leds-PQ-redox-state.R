library(photobiologyPlants)
library(photobiologyLEDs)
library(dplyr)
library(ggpmisc)

peaks(leds.mspct, span = NULL) |> rbind_spct() |> arrange(w.length) -> led_peaks.tb
peaks(led_arrays.mspct, span = NULL) |> rbind_spct() |> arrange(w.length) -> array_peaks.tb

PQ_redox_state(leds.mspct$LedEngin_LZ4_40R208_660nm)
peaks(leds.mspct$LedEngin_LZ4_40R208_660nm, span = NULL)

PQ_redox_state(led_arrays.mspct$Weili_120W.array.12ch.custom.A) -> z

full_join(z, peaks(subset2mspct(led_arrays.mspct$Weili_120W.array.12ch.custom.A),
                   span = NULL) |>
            rbind_spct())

ref.spct <- leds.mspct$LedEngin_LZ4_40R208_660nm
leds_simul.mspct <- source_mspct()
for (i in 1:110) {
  temp.spct <- ref.spct
  temp.spct$w.length <- temp.spct$w.length - 70 + i
  spct.name <- paste("spct_", peaks(temp.spct, span = NULL)$w.length, "nm", sep = "")
  leds_simul.mspct[[spct.name]] <- temp.spct
}

PQ_rdx.tb <- PQ_redox_state(leds_simul.mspct)
PQ_rdx.tb$peak.wl <- as.numeric(gsub("^spct_|nm$", "", PQ_rdx.tb$spct.idx))

ggplot(PQ_rdx.tb, aes(peak.wl, PQ_redox.state)) +
  geom_line() +
  stat_peaks() +
  stat_peaks(geom = "text", vjust = -0.5) +
  stat_valleys(global.threshold = -0.9) +
  stat_valleys(geom = "text", vjust = 1.5, global.threshold = -0.9) +
  scale_x_continuous(name = "Wavelength at peak (nm)", n.breaks = 10) +
  scale_y_continuous(name = "Platoquinone redox state (%)") +
  labs(title = "Platoquinone redox state in light from red LEDs",
       subtitle = "LEDs simulated by shifting the measured spectrum of a 660 nm LED",
       caption = "Led Engin type LZ4-40R208 uses as starting point with peak at 662.71 nm") +
  theme_classic()

  
