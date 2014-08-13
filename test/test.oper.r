library(ggplot2)
library(photobiologyAll)
library(profr)

test <- copy(sun.spct)

Rprof(interval = 0.0001)
z <- test * 2
Rprof(NULL)
summaryRprof()

prof <- profr({z <- test * 2}, interval = 0.001)
ggplot(prof)

prof <- profr({z <- test + test}, interval = 0.001)
ggplot(prof)

test1 <- interpolate_spct(test, length.out=1000)

prof <- profr({z <- test1 + test}, interval = 0.001)
ggplot(prof)
