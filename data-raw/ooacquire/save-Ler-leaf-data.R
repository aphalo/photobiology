library(photobiology)
library(dplyr)

load("data-raw/ooacquire/Ler_06_black.Rda")

Ler_06_black.rfr_spct %>%
  trimInstrDesc() %>%
  trimInstrSettings() -> Ler_leaf_rflt.spct

Ler_06_black.tfr_spct %>%
  trimInstrDesc() %>%
  trimInstrSettings() -> Ler_leaf_trns.spct

Ler_leaf.spct <-
  merge(smooth_spct(Ler_leaf_rflt.spct, method = "supsmu"),
        smooth_spct(Ler_leaf_trns.spct, method = "supsmu"),
        w.length.out = seq(from = 250, to = 850, by = 0.25),
        Tfr.type.out = "total")

Ler_leaf_i.spct <-
  merge(smooth_spct(Ler_leaf_rflt.spct, method = "supsmu"),
        smooth_spct(Ler_leaf_trns.spct, method = "supsmu"),
        w.length.out = seq(from = 250, to = 850, by = 0.25),
        Tfr.type.out = "internal")

Ler_leaf_trns_i.spct <- as.filter_spct(Ler_leaf_i.spct) %>% select(-Rfr)

save(Ler_leaf_rflt.spct,
     Ler_leaf_trns.spct,
     Ler_leaf_trns_i.spct,
     Ler_leaf.spct,
     file = "data/Ler-leaf-spct.rda")

