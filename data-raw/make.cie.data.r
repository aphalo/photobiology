library(readxl)
library(photobiology)

D65.illuminant.spct <- read_excel(path = "data-raw/CIE/204.xls", sheet = "D65", skip = 5, col_names = FALSE)
D65.illuminant.spct <- na.omit(D65.illuminant.spct)
names(D65.illuminant.spct)[1:2] <- c("w.length", "s.e.irrad")
D65.illuminant.spct <- as.source_spct(D65.illuminant.spct)
D65.illuminant.spct <- normalize(D65.illuminant.spct, norm = 560)
comment(D65.illuminant.spct) <- "CIE D65 standard illuminant, normalized to one at 560 nm"

save(D65.illuminant.spct, file = "./data/D65.illuminant.spct.rda")
rm(D65.illuminant.spct)

A.illuminant.spct <- read_excel(path = "data-raw/CIE/204.xls", sheet = "Ill. A", skip = 5, col_names = FALSE)
A.illuminant.spct <- na.omit(A.illuminant.spct)
names(A.illuminant.spct)[1:2] <- c("w.length", "s.e.irrad")
A.illuminant.spct <- as.source_spct(A.illuminant.spct)
A.illuminant.spct <- normalize(A.illuminant.spct, norm = 560)
comment(A.illuminant.spct) <- "CIE A standard illuminant, normalized to one at 560 nm"

save(A.illuminant.spct, file = "./data/A.illuminant.spct.rda")
rm(A.illuminant.spct)

