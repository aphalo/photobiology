library(photobiology)

oldwd <- setwd("raw.data/humans")

# Human
ciexyzCMF2.spct <- read.csv(file="lin2012xyz2e_1_7sf.csv", comment.char = "#")
ciexyzCC2.spct <- read.csv(file="cc2012xyz2_1_5dp.csv", comment.char = "#")
ciexyzCMF10.spct <- read.csv(file="lin2012xyz10e_1_7sf.csv", comment.char = "#")
ciexyzCC10.spct <- read.csv(file="cc2012xyz10_1_5dp.csv", comment.char = "#")
setChromaSpct(ciexyzCMF2.spct)
setChromaSpct(ciexyzCC2.spct)
setChromaSpct(ciexyzCMF10.spct)
setChromaSpct(ciexyzCC10.spct)

setwd(oldwd)

olwd <- setwd("data")

# save(ciexyzCMF2.data, ciexyzCMF10.data, ciexyzCC2.data, ciexyzCC10.data, file="ciexyz2006.data.rda")
save(ciexyzCMF2.spct, ciexyzCMF10.spct, ciexyzCC2.spct, ciexyzCC10.spct, file="ciexyz2006.spct.rda")

setwd(olwd)

# Honeybee

oldwd <- setwd("raw.data/bees")

Maxwell.data <- read.table(file="Maxwell.data", header=TRUE)

beesxyzCMF.spct <- Maxwell.data[ , 1:4]
setnames(beesxyzCMF.spct, 1:4, c("w.length", "x", "y", "z"))
setChromaSpct(beesxyzCMF.spct)

setwd(oldwd)

olwd <- setwd("data")
save(beesxyzCMF.spct, file="beesxyz.spct.rda")
setwd(olwd)
