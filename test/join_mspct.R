library(dplyr)
library(photobiology)
library(photobiologyLEDs)
library(plyr)

test.mspct <- leds.mspct[1:9]
# join_all(test.mspct, by = "w.length")

join_mspct <- function(mspct, by = "w.length",
                       type = "full", match = "first",
                       unit_out = "energy") {
  stopifnot(is.any_mspct(mspct))
  names <- names(mspct)
  stopifnot(length(names) == length(mspct))
  if (unit_out == "energy") {
    mspct <- q2e(mspct, action = "replace")
    rmDerivedMspct(mspct)
    for (i in names) {
      mspct[[i]] <- as.data.frame(mspct[[i]])[c("w.length", "s.e.irrad")]
      mspct[[i]] %>%
        dplyr::rename(s.e.irrad = i) -> mspct[[i]]
    }
  } else if (unit_out %in% c("photon", "quantum")) {
    mspct <- e2q(mspct, action = "replace")
    rmDerivedMspct(mspct)
    for (i in names) {
      mspct[[i]] <- as.data.frame(mspct[[i]])[c("w.length", "s.q.irrad")]
      mspct[[i]] %>%
        dplyr::rename(s.q.irrad = i) -> mspct[[i]]
    }
  } else {
    stop("Unit out '", unit_out, "' unknown")
  }
  plyr::join_all(dfs = mspct, by = by, type = type, match = match)
}

energy.df <- print(join_mspct(test.mspct))
photon.df <- print(join_mspct(test.mspct, unit_out = "photon"))

