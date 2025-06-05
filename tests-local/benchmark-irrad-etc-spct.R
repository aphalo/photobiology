# This are just benchmarks
# We need to check also if values are correct
library(microbenchmark)
library(photobiology)
library(photobiologyWavebands)
library(ggplot2)

caption <-
  paste(trimws(unname(Sys.info()["nodename"])),
        "; photobiology ", packageVersion("photobiology"),
        "; tibble ", packageVersion("tibble"),
        "; dplyr ", packageVersion("dplyr"),
        with(R.version, paste("; R ", major, ".", minor, sep = "")), sep = "")

out.folder <- paste("./tests-local/spct-irrad-test", format(lubridate::today()), sep = "-")

if (!dir.exists(out.folder)) {
  dir.create(out.folder)
}

mbm.irrad.spct <-
  microbenchmark(irrad(sun.spct),
                 irrad(sun.spct, range = c(400, 700)),
                 irrad(sun.spct, w.band = c(400, 700)),
                 irrad(sun.spct, w.band = VIS_bands()),
                 irrad(sun.spct, w.band = CIE()),
                 irrad(sun.spct, w.band = PAR()),
                 q_irrad(sun.spct),
                 q_irrad(sun.spct, range = c(400, 700)),
                 q_irrad(sun.spct, w.band = c(400, 700)),
                 q_irrad(sun.spct, w.band = VIS_bands()),
                 q_irrad(sun.spct, w.band = CIE()),
                 q_irrad(sun.spct, w.band = PAR()),
                 e_irrad(sun.spct),
                 e_irrad(sun.spct, range = c(400, 700)),
                 e_irrad(sun.spct, w.band = c(400, 700)),
                 times = 300, unit = "second", control = list(warmup = 5))

fig1 <-
autoplot(mbm.irrad.spct, y_max = 0.0015) +
  labs(caption = caption)
fig1

svglite(paste(out.folder, "irrad.spct.svg", sep = "/"), width = 6, height = 4)
print(fig1)
dev.off()

mbm.irrad.lspct <-
  microbenchmark(irrad(sun_evening.spct),
                 irrad(sun_evening.spct, range = c(400, 700)),
                 irrad(sun_evening.spct, w.band = c(400, 700)),
                 irrad(sun_evening.spct, w.band = VIS_bands()),
                 irrad(sun_evening.spct, w.band = CIE()),
                 irrad(sun_evening.spct, w.band = PAR()),
                 q_irrad(sun_evening.spct),
                 q_irrad(sun_evening.spct, range = c(400, 700)),
                 q_irrad(sun_evening.spct, w.band = c(400, 700)),
                 q_irrad(sun_evening.spct, w.band = VIS_bands()),
                 q_irrad(sun_evening.spct, w.band = CIE()),
                 q_irrad(sun_evening.spct, w.band = PAR()),
                 e_irrad(sun_evening.spct),
                 e_irrad(sun_evening.spct, range = c(400, 700)),
                 e_irrad(sun_evening.spct, w.band = c(400, 700)),
                 times = 300, unit = "second", control = list(warmup = 5))

fig2 <-
autoplot(mbm.irrad.lspct, y_max = 0.015) +
  labs(caption = caption)
fig2

svglite(paste(out.folder, "irrad.lspct.svg", sep = "/"), width = 6, height = 4)
print(fig2)
dev.off()

colnames(sun_evening.mspct[[1]])

mbm.irrad.mspct <-
  microbenchmark(irrad(sun_evening.mspct),
                 irrad(sun_evening.mspct, range = c(400, 700)),
                 irrad(sun_evening.mspct, w.band = c(400, 700)),
                 irrad(sun_evening.mspct, w.band = VIS_bands()),
                 irrad(sun_evening.mspct, w.band = CIE()),
                 irrad(sun_evening.mspct, w.band = PAR()),
                 q_irrad(sun_evening.mspct),
                 q_irrad(sun_evening.mspct, range = c(400, 700)),
                 q_irrad(sun_evening.mspct, w.band = c(400, 700)),
                 q_irrad(sun_evening.mspct, w.band = VIS_bands()),
                 q_irrad(sun_evening.mspct, w.band = CIE()),
                 q_irrad(sun_evening.mspct, w.band = PAR()),
                 e_irrad(sun_evening.mspct),
                 e_irrad(sun_evening.mspct, range = c(400, 700)),
                 e_irrad(sun_evening.mspct, w.band = c(400, 700)),
                 times = 300, unit = "second", control = list(warmup = 5))

fig3 <-
autoplot(mbm.irrad.mspct, y_max = 0.015) +
  labs(caption = caption)
fig3

svglite(paste(out.folder, "irrad.mspct.svg", sep = "/"), width = 6, height = 4)
print(fig3)
dev.off()


## Ratios

mbm.ratio.spct <-
  microbenchmark(q_ratio(sun.spct, UVA(), PAR()),
                 q_ratio(sun.spct, VIS_bands(), PAR()),
                 q_ratio(sun.spct, PAR(), VIS_bands()),
                 e_ratio(sun.spct, UVA(), PAR()),
                 e_ratio(sun.spct, VIS_bands(), PAR()),
                 e_ratio(sun.spct, PAR(), VIS_bands()),
                 qe_ratio(sun.spct, PAR()),
                 qe_ratio(sun.spct, VIS_bands()),
                 eq_ratio(sun.spct, PAR()),
                 eq_ratio(sun.spct, VIS_bands()),
                 times = 300, unit = "second", control = list(warmup = 5))

fig4 <-
autoplot(mbm.ratio.spct, y_max = 0.015) +
  labs(caption = caption)
fig4

svglite(paste(out.folder, "ratio.spct.svg", sep = "/"), width = 6, height = 4)
print(fig4)
dev.off()

mbm.ratio.lspct <-
  microbenchmark(q_ratio(sun_evening.spct, UVA(), PAR()),
                 q_ratio(sun_evening.spct, VIS_bands(), PAR()),
                 q_ratio(sun_evening.spct, PAR(), VIS_bands()),
                 e_ratio(sun_evening.spct, UVA(), PAR()),
                 e_ratio(sun_evening.spct, VIS_bands(), PAR()),
                 e_ratio(sun_evening.spct, PAR(), VIS_bands()),
                 qe_ratio(sun_evening.spct, PAR()),
                 qe_ratio(sun_evening.spct, VIS_bands()),
                 eq_ratio(sun_evening.spct, PAR()),
                 eq_ratio(sun_evening.spct, VIS_bands()),
                 times = 300, unit = "second", control = list(warmup = 5))

fig5 <-
autoplot(mbm.ratio.lspct, y_max = 0.04) +
  expand_limits(x = 0.002) +
  labs(caption = caption)
fig5

svglite(paste(out.folder, "ratio.lspct.svg", sep = "/"), width = 6, height = 4)
print(fig5)
dev.off()

colnames(sun_evening.mspct[[1]])

mbm.ratio.mspct <-
  microbenchmark(q_ratio(sun_evening.mspct, UVA(), PAR()),
                 q_ratio(sun_evening.mspct, VIS_bands(), PAR()),
                 q_ratio(sun_evening.mspct, PAR(), VIS_bands()),
                 e_ratio(sun_evening.mspct, UVA(), PAR()),
                 e_ratio(sun_evening.mspct, VIS_bands(), PAR()),
                 e_ratio(sun_evening.mspct, PAR(), VIS_bands()),
                 qe_ratio(sun_evening.mspct, PAR()),
                 qe_ratio(sun_evening.mspct, VIS_bands()),
                 eq_ratio(sun_evening.mspct, PAR()),
                 eq_ratio(sun_evening.mspct, VIS_bands()),
                 times = 300, unit = "second", control = list(warmup = 5))

fig6 <-
autoplot(mbm.ratio.mspct, y_max = 0.04) +
  expand_limits(x = 0.002) +
  labs(caption = caption)
fig6

svglite(paste(out.folder, "ratio.mspct.svg", sep = "/"), width = 6, height = 4)
print(fig2)
dev.off()

save(caption,
     mbm.irrad.spct, mbm.irrad.lspct, mbm.irrad.mspct,
     mbm.ratio.spct, mbm.ratio.lspct, mbm.ratio.mspct,
     file = paste(out.folder, "microbenchmark-results.rda", sep = "/"))
