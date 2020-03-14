library(hyperSpec)
library(photobiologyInOut)
library(ggspectra)

load("./inst-not/Raman spikes/TestData.RData")
wl(test_data) <- list(wl = 1 / wl(test_data) * 1e-2 * 1e9,
                      label = expression(lambda / nm))
chk.hy(test_data)
Raman.spct <- hyperSpec2spct(test_data,
                             member.class = "response_spct",
                             spct.data.var = "s.e.response")
Raman.mspct <- hyperSpec2mspct(test_data,
                               member.class = "response_spct",
                               spct.data.var = "s.e.response")

rbindspct(Raman.mspct)
comment(Raman.mspct[[1]])

autoplot(Raman.mspct[1:5])

