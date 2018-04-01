library(ggplot2)
library(ggspectra)
library(photobiology)
library(hexSticker)

p <- ggplot(normalize(clip_wl(sun.spct, range = c(380,700)))) +
  labs(x = NULL, y = NULL) + theme_void() +
  theme_transparent() + stat_wl_strip(ymin = -Inf, ymax = Inf) +
#  geom_line(colour = "white") + geom_line(colour = "black", size = 0.2) +
  scale_fill_identity()
#p <- p + theme_sticker()

sticker(p, package="R for Photobiology", p_size=10, p_y = 0.5, s_x=1, s_y=0.65, s_width=1.85, s_height=0.25,
        h_fill="white", h_color="dodgerblue2", filename="test.png")
