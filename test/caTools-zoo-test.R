library(caTools)
library(zoo)

vec <- rnorm(2000)

## min

min.catools <- runmin(vec, k = 5)

length(min.catools)

min.zoo <- rollapply(vec, width = 5, FUN = min, partial = TRUE)

length(min.zoo)

all(min.catools == min.zoo)

## mad

mad.catools <- runmad(vec, k = 15)

length(mad.catools)

mad.zoo <- rollapply(vec, width = 15, FUN = mad, partial = TRUE)

length(mad.zoo)

mad.different <- which(mad.catools != mad.zoo)

mad.different

mad.catools[mad.different]

mad.zoo[mad.different]


