# This are just benchmarks
# We need to check also if values are correct
library(microbenchmark)
library(photobiology)
library(tibble)
library(ggplot2)
library(svglite)

caption <-
  paste(trimws(unname(Sys.info()["nodename"])),
        "; photobiology ", packageVersion("photobiology"),
        "; tibble ", packageVersion("tibble"),
        "; dplyr ", packageVersion("dplyr"),
        with(R.version, paste("; R ", major, ".", minor, sep = "")), sep = "")

out.folder <- paste("./tests-local/spct-extract-test", format(lubridate::today()), sep = "-")

if (!dir.exists(out.folder)) {
  dir.create(out.folder)
}

x.spct <- sun.spct
x.df <- sun.data
x.tb <- as_tibble(sun.data)

mbm.spct <-
microbenchmark(z.clp.spct <- clip_wl(x.spct, range = c(400, 700)),
               z.trmwl.spct <- trim_wl(x.spct, range = c(400, 700)),
               z.trmspct.spct <- trim_spct(x.spct, range = c(400, 700)),
               z.xtr.spct <- x.spct[x.spct$w.length >= 400 & x.spct$w.length < 700, ],
               z.xtr.tb <- x.tb[x.tb$w.length >= 400 & x.tb$w.length < 700, ],
               z.xtr.df <- x.df[x.df$w.length >= 400 & x.df$w.length < 700, ])

fig1 <-
  autoplot(mbm.spct) +
  labs(caption = caption)

fig1

svglite(paste(out.folder, "spct.svg", sep = "/"), width = 6, height = 4)
print(fig1)
dev.off()

x.spct <- sun_evening.spct
x.df <- as.data.frame(sun_evening.spct)
class(x.df)
x.tb <- as_tibble(sun.data)
class(x.tb)

mbm.lspct <-
microbenchmark(z.clp.spct <- clip_wl(x.spct, range = c(400, 700)),
               z.trm.spct <- trim_wl(x.spct, range = c(400, 700)),
               z.xtr.spct <- x.spct[x.spct$w.length >= 400 & x.spct$w.length <= 700, ],
               z.xtr.tb <- x.tb[x.tb$w.length >= 400 & x.tb$w.length <= 700, ],
               z.xtr.df <- x.df[x.df$w.length >= 400 & x.df$w.length <= 700, ],
               times = 300, unit = "second", control = list(warmup = 5))

fig2 <-
autoplot(mbm.lspct) +
  labs(caption = caption)

svglite(paste(out.folder, "lspct.svg", sep = "/"), width = 6, height = 4)
print(fig2)
dev.off()

mbm.lmspct <-
  microbenchmark(z.clp.spct <- clip_wl(sun_evening.spct, range = c(400, 700)),
                 z.trm.spct <- trim_wl(sun_evening.spct, range = c(400, 700)),
                 z.clp.spct <- clip_wl(sun_evening.mspct, range = c(400, 700)),
                 z.trm.spct <- trim_wl(sun_evening.mspct, range = c(400, 700)),
                 times = 300, unit = "second", control = list(warmup = 5))

fig3 <-
  autoplot(mbm.lmspct) +
  labs(caption = caption)

fig3

svglite(paste(out.folder, "lmspct.svg", sep = "/"), width = 6, height = 4)
print(fig3)
dev.off()

## subset

x.spct <- sun_evening.spct

mbm.subset <-
  microbenchmark(z.sb.spct <- subset(x.spct, w.length >= 400 & x.spct$w.length <= 700),
                 z.xtr.spct <- x.spct[x.spct$w.length >= 400 & x.spct$w.length <= 700, ],
                 z.sb.tb <- subset(x.tb, w.length >= 400 & x.tb$w.length <= 700),
                 z.xtr.df <- subset(x.df, w.length >= 400 & x.df$w.length <= 700),
                 times = 1000, unit = "second", control = list(warmup = 5))

autoplot(mbm.subset) +
  labs(caption = caption)

mbm.pull <-
  microbenchmark(z.smp4.spct <- pull_sample(x.spct, 1),
                 z.smp2.spct <- pull_sample(x.spct, 2),
                 z.smp10.spct <- pull_sample(x.spct, 3),
                 z.smp1.spct <- pull_sample(x.spct, 4),
                 times = 1000, unit = "second", control = list(warmup = 5))

autoplot(mbm.pull) +
  labs(caption = caption)

save(caption,
     mbm.spct, mbm.lspct, mbm.lmspct, mbm.subset, mbm.pull,
     file = paste(out.folder, "microbenchmark-results.rda", sep = "/"))
