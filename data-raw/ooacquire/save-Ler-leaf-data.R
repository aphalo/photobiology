library(photobiology)
library(dplyr)

load("data-raw/ooacquire/Ler_06_black.Rda")

Ler_06_black.rfr_spct %>%
  trimInstrDesc() %>%
  trimInstrSettings() %>%
  trim_wl(c(250, 850)) %>%
  smooth_spct(method = "supsmu", strength = 3) %>%
  clean() -> Ler_leaf_rflt.spct

Ler_06_black.tfr_spct %>%
  trimInstrDesc() %>%
  trimInstrSettings() %>%
  trim_wl(c(250, 850)) %>%
  smooth_spct(method = "supsmu", strength = 3) %>%
  clean() -> Ler_leaf_trns.spct

merge2object_spct(Ler_leaf_rflt.spct,
                  Ler_leaf_trns.spct,
                  w.length.out = seq(from = 250, to = 850, by = 0.25),
                  Tfr.type.out = "total") %>%
  na.omit() %>%
  clean() -> Ler_leaf.spct

getTfrType(Ler_leaf.spct)
getRfrType(Ler_leaf.spct)

Ler_leaf.spct %>%
  convertTfrType(Tfr.type = "internal") %>%
  as.filter_spct()  -> Ler_leaf_trns_i.spct

save(Ler_leaf_rflt.spct,
     Ler_leaf_trns.spct,
     Ler_leaf_trns_i.spct,
     Ler_leaf.spct,
     file = "data/Ler-leaf-spct.rda")

