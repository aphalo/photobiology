library(readxl)
library(photobiology)

D65.illuminant.spct <- read.csv(file = "data-raw/CIE/CIE_std_illum_D65.csv", header = FALSE, col.names = c("w.length", "s.e.irrad"))
D65.illuminant.spct <- as.source_spct(D65.illuminant.spct)
D65.illuminant.spct <- normalize(D65.illuminant.spct, norm = 560)
comment(D65.illuminant.spct) <- "CIE D65 standard illuminant, normalized to one at 560 nm"
setWhatMeasured(D65.illuminant.spct, "CIE D65 standard illuminant (DOI: 10.25039/CIE.DS.hjfjmt59), normalized to one at 560 nm.")
setHowMeasured(D65.illuminant.spct, "Spectral data imported from a published .CSV file.")
attr(D65.illuminant.spct, "header") <- fromJSON("data-raw/CIE/CIE_std_illum_D65.csv_metadata.json")
save(D65.illuminant.spct, file = "./data/D65.illuminant.spct.rda")

D50.illuminant.spct <- read.csv(file = "data-raw/CIE/CIE_std_illum_D50.csv", header = FALSE, col.names = c("w.length", "s.e.irrad"))
D50.illuminant.spct <- as.source_spct(D50.illuminant.spct)
D50.illuminant.spct <- normalize(D50.illuminant.spct, norm = 560)
comment(D50.illuminant.spct) <- "CIE D50 standard illuminant, normalized to one at 560 nm"
setWhatMeasured(D50.illuminant.spct, "CIE D50 standard illuminant (DOI: 10.25039/CIE.DS.etgmuqt5), normalized to one at 560 nm.")
setHowMeasured(D50.illuminant.spct, "Spectral data imported from a published .CSV file.")
attr(D50.illuminant.spct, "header") <- fromJSON("data-raw/CIE/CIE_std_illum_D50.csv_metadata.json")
save(D50.illuminant.spct, file = "./data/D50.illuminant.spct.rda")

A.illuminant.spct <- read.csv(file = "data-raw/CIE/CIE_std_illum_A_1nm.csv", header = FALSE, col.names = c("w.length", "s.e.irrad"))
A.illuminant.spct <- as.source_spct(A.illuminant.spct)
A.illuminant.spct <- normalize(A.illuminant.spct, norm = 560)
comment(A.illuminant.spct) <- "CIE A standard illuminant, normalized to one at 560 nm"
setWhatMeasured(A.illuminant.spct, "CIE A standard illuminant (DOI: 10.25039/CIE.DS.8jsxjrsn), normalized to one at 560 nm.")
setHowMeasured(A.illuminant.spct, "Spectral data imported from a published .CSV file.")
attr(A.illuminant.spct, "header") <- fromJSON("data-raw/CIE/CIE_std_illum_A_1nm.csv_metadata.json")
save(A.illuminant.spct, file = "./data/A.illuminant.spct.rda")

# rm(D65.illuminant.spct)
# rm(D50.illuminant.spct)
# rm(A.illuminant.spct)
