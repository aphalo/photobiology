library(photobiology)

water.tb <- read.table("data-raw/solutes/water-pure-buiteveld94.txt",
                       skip = 6,
                       col.names = c("w.length", "K.mole"))
comment(water.tb) <-
  paste(readLines("data-raw/solutes/water-pure-buiteveld94.txt", n = 6),
        collapse = "\n")

water.spct <- as.solute_spct(water.tb)

water.properties <-
  list(formula = c(text = "H2O", html = "H<sub>2</sub>", TeX = "$H_2O$"),
       name = c("water", IUPAC = "oxidane"),
       structure = grDevices::as.raster(matrix()),
       mass = 18.015, # Da
       ID = c(ChemSpider = "917", CID = "962", CAS = "7732-18-5"),
       solvent.name = NA_character_,
       solvent.ID = NA_character_)

setSoluteProperties(water.spct, water.properties)
water.spct
solute_properties(water.spct)

phenylalanine.tb <- read.table("data-raw/solutes/phenylalanine-abs.txt",
                               skip = 23,
                               col.names = c("w.length", "K.mole"))
comment(phenylalanine.tb) <-
  paste(readLines("data-raw/solutes/phenylalanine-abs.txt", n = 23),
        collapse = "\n")

phenylalanine.spct <- as.solute_spct(phenylalanine.tb)

phenylalanine.properties <-
  list(formula = c(text = "C9H11NO2",
                   html = "C<sub>9</sub>H<sub>11</sub>NO<sub>2</sub>",
                   TeX = "$C_9H_{11}NO_2$"),
       name = c("phenylalanine", IUPAC = "(2S)-2-amino-3-phenylpropanoic acid"),
       structure = grDevices::as.raster(matrix()),
       mass = 165.19, # Da
       ID = c(ChemSpider = "969", CID = "6140", CAS = "63-91-2"),
       solvent.name = c("water", IUPAC = "oxidane"),
       solvent.ID = c(ChemSpider = "917", CID = "962", CAS = "7732-18-5"))

setSoluteProperties(phenylalanine.spct, phenylalanine.properties)
phenylalanine.spct
solute_properties(phenylalanine.spct)
str(solute_properties(phenylalanine.spct))

save(water.spct, phenylalanine.spct, file = "./data/solutes.rda")
