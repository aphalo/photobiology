library(photobiology)
library(dplyr)

energy_as_default()

oldwd <- setwd("data-raw/humans")

# Human
ciexyzCMF2.spct <- read.csv(file="lin2012xyz2e_1_7sf.csv", comment.char = "#")
ciexyzCC2.spct <- read.csv(file="cc2012xyz2_1_5dp.csv", comment.char = "#")
ciev2.spct <- read.csv(file="linCIE2008v2e_1.csv", comment.char = "#")
ciexyzCMF10.spct <- read.csv(file="lin2012xyz10e_1_7sf.csv", comment.char = "#")
ciexyzCC10.spct <- read.csv(file="cc2012xyz10_1_5dp.csv", comment.char = "#")
ciev10.spct <- read.csv(file="linCIE2008v10e_1.csv", comment.char = "#")
cone_fundamentals10.spct <- read.csv(file="linss10e_1.csv", col.names = c("w.length", "x", "y", "z"))
setChromaSpct(ciexyzCMF2.spct)
what_measured(ciexyzCMF2.spct) <- "CIE 2012 2 degrees CMF"
comment(ciexyzCMF2.spct) <- "CIE 2012 2 degrees CMF (color matching function) from lin2012xyz2e_1_7sf.csv"
setChromaSpct(ciexyzCC2.spct)
what_measured(ciexyzCC2.spct) <- "CIE 2012 2 degrees CC"
comment(ciexyzCC2.spct) <- "CIE 2012 2 degrees CC (color coordinates) from cc2012xyz2_1_5dp.csv"
setResponseSpct(ciev2.spct)
what_measured(ciev2.spct) <- "CIE 2008 2 degrees V"
comment(ciev2.spct) <- "CIE 2008 2 degrees V from linCIE2008v2e_1.csv"
setChromaSpct(ciexyzCMF10.spct)
what_measured(ciexyzCMF10.spct) <- "CIE 2008 10 degrees CMF"
comment(ciexyzCMF10.spct) <- "CIE 2012 10 degrees CMF (color matching function) from lin2012xyz10e_1_7sf.csv"
setChromaSpct(ciexyzCC10.spct)
what_measured(ciexyzCC10.spct) <- "CIE 2008 10 degrees CC"
comment(ciexyzCC10.spct) <- "CIE 2012 10 degrees CC (color coordinates) from cc2012xyz10_1_5dp.csv"
setResponseSpct(ciev10.spct)
what_measured(ciev10.spct) <- "CIE 2008 10 degrees V"
comment(ciev10.spct) <- "CIE 2008 10 degrees V from linCIE2008v10e_1.csv"
cone_fundamentals10.spct %>%
  split2response_mspct() %>%
  msmsply(`what_measured<-`, value = "Cone fundamentals 10 degrees.") %>%
  msmsply(normalize, norm = "max") %>%
  msmsply(na.omit) ->
  cone_fundamentals10.mspct
setChromaSpct(cone_fundamentals10.spct) %>% mutate(z = ifelse(is.na(z), 0, z)) -> cone_fundamentals10.spct
what_measured(cone_fundamentals10.spct) <- "Cone fundamentals 10 degrees"
comment(cone_fundamentals10.spct) <- "10-deg cone fundamentals based on the Stiles & Burch 10-deg CMFs; Stockman & Sharpe (2000)"
setwd(oldwd)

olwd <- setwd("data")

# save(ciexyzCMF2.data, ciexyzCMF10.data, ciexyzCC2.data, ciexyzCC10.data, file="ciexyz2006.data.rda")
save(ciexyzCMF2.spct, file="ciexyzCMF2.spct.rda")
save(ciexyzCMF10.spct, file="ciexyzCMF10.spct.rda")
save(ciexyzCC2.spct, file="ciexyzCC2.spct.rda")
save(ciexyzCC10.spct, file="ciexyzCC10.spct.rda")
save(ciev2.spct, file="ciev2.spct.rda")
save(ciev10.spct, file="ciev10.spct.rda")
save(cone_fundamentals10.spct, cone_fundamentals10.mspct, file="cone_fundamentals10.spct.rda")

setwd(olwd)

# Honeybee

oldwd <- setwd("data-raw/bees")

Maxwell.data <- read.table(file="Maxwell.data", header=TRUE)

beesxyzCMF.spct <- Maxwell.data[ , 1:4]
names(beesxyzCMF.spct)[1:4] <- c("w.length", "z", "y", "x")
setChromaSpct(beesxyzCMF.spct)
what_measured(beesxyzCMF.spct) <- "Honey bee CMF"
comment(beesxyzCMF.spct) <- "Maxwell's color matching function for honey bee, interpolated to 1 nm"
beesxyzCMF.spct <- interpolate_spct(beesxyzCMF.spct, 300:700)
setwd(oldwd)

olwd <- setwd("data")
save(beesxyzCMF.spct, file="beesxyz.spct.rda")
setwd(olwd)

