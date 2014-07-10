library(photobiology)
library(data.table)
oldwd <- setwd("raw.data/bees")
Maxwell.data <- read.table(file="Maxwell.data", header=TRUE)
bee.coordinates <- Maxwell.data[, 1:4]
D65.data <- Maxwell.data[, c(1,5)]
flower.data <- Maxwell.data[, c("w.length", "Flower")]
leaf.data <- Maxwell.data[, c("w.length", "Leaf")]
1 / sum(Maxwell.data$UV * D65.data$D65 * leaf.data$Leaf)
sum(Maxwell.data$UV * D65.data$D65 * leaf.data$Leaf)

leaf.conv.spct <- with(Maxwell.data, data.table(w.length = w.length,
                                                U = UV * D65 * Leaf,
                                                B = Blue * D65 * Leaf,
                                                G = Green * D65 * Leaf))

flower.conv.spct <- with(Maxwell.data, data.table(w.length = w.length,
                                                U = UV * D65 * Flower,
                                                B = Blue * D65 * Flower,
                                                G = Green * D65 * Flower))

R.flower <- 1 / integrate_spct(flower.conv.spct)
R.leaf <- 1 / integrate_spct(leaf.conv.spct)

P.flower <- R.leaf * 1 / R.flower

E.flower <- P.flower / (P.flower + 1)

rel.q.abs.flower <- P.flower / sum(P.flower)

bee.xyz.coord <- c(x = 0.8667 * (P.flower["U"] - P.flower["G"]), y = P.flower["B"] - 0.5 * (P.flower["G"] + P.flower["U"]))

setwd(oldwd)
