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

merge(Ler_leaf_rflt.spct,
      Ler_leaf_trns.spct,
      w.length.out = seq(from = 250, to = 850, by = 0.25),
      Tfr.type.out = "total") %>%
  na.omit() -> Ler_leaf.spct

merge(Ler_leaf_rflt.spct,
      Ler_leaf_trns.spct,
      w.length.out = seq(from = 250, to = 850, by = 0.25),
      Tfr.type.out = "internal") %>%
  na.omit() -> Ler_leaf_i.spct

Ler_leaf_trns_i.spct
  Ler_leaf_i.spct %>%
  as.filter_spct() %>%
  select(-Rfr) %>%
  clean() -> Ler_leaf_trns_i.spct

save(Ler_leaf_rflt.spct,
     Ler_leaf_trns.spct,
     Ler_leaf_trns_i.spct,
     Ler_leaf.spct,
     file = "data/Ler-leaf-spct.rda")

