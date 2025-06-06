library("photobiology")

context("peaks")

test_that("source_spct single", {

  my.spct <- sun.spct

  peaks.spct <- peaks(sun.spct, span = NULL, strict = TRUE)
  expect_equal(nrow(peaks.spct), 1)
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, span = NULL, strict = FALSE)
  expect_equal(nrow(peaks.spct), 1)
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, span = NULL, strict = TRUE, global.threshold = 0.9)
  expect_equal(nrow(peaks.spct), 1)
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, span = NULL, strict = FALSE, global.threshold = -0.1)
  expect_equal(nrow(peaks.spct), 0)
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, span = NULL, strict = TRUE,
                      global.threshold = 0.9, threshold.range = c(0, 0.82))
  expect_equal(nrow(peaks.spct), 1)
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, span = 5)
  expect_equal(nrow(peaks.spct), 76)
  expect_equal(names(peaks.spct), c("w.length", "s.e.irrad"))
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, span = 51)
  expect_equal(nrow(peaks.spct), 3)
  expect_equal(names(peaks.spct), c("w.length", "s.e.irrad"))
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, span = 5) # default
  expect_equal(nrow(peaks.spct), 76)
  expect_equal(names(peaks.spct), c("w.length", "s.e.irrad"))
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct)
  expect_equal(nrow(peaks.spct), 76)
  expect_equal(names(peaks.spct), c("w.length", "s.e.irrad"))
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, global.threshold = 0.9)
  expect_equal(nrow(peaks.spct), 15)
  expect_equal(names(peaks.spct), c("w.length", "s.e.irrad"))
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, global.threshold = -0.1)
  expect_equal(nrow(peaks.spct), 0)
  expect_equal(names(peaks.spct), c("w.length", "s.e.irrad"))
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, local.threshold = 0.1)
  expect_equal(nrow(peaks.spct), 16)
  expect_equal(names(peaks.spct), c("w.length", "s.e.irrad"))
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, local.threshold = 0.1, local.reference = "minimum")
  expect_equal(nrow(peaks.spct), 16)
  expect_equal(names(peaks.spct), c("w.length", "s.e.irrad"))
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, local.threshold = 0.05, local.reference = "median")
  expect_equal(nrow(peaks.spct), 9)
  expect_equal(names(peaks.spct), c("w.length", "s.e.irrad"))
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, span = NULL, unit.out = "photon")
  expect_equal(nrow(peaks.spct), 1)

  peaks.spct <- peaks(my.spct, unit.out = "photon")

  expect_equal(nrow(peaks.spct), 77)
  expect_equal(names(peaks.spct), c("w.length", "s.q.irrad"))
  expect_is(peaks.spct, "source_spct")

  my_thn.spct <- thin_wl(sun.spct)
  expect_silent(peaks(my_thn.spct))
  expect_silent(peaks(my_thn.spct, span = 5L))
  expect_silent(peaks(my_thn.spct, span = NULL))
  expect_warning(peaks(my_thn.spct, span = 7L))
  expect_warning(peaks(my_thn.spct, span = 101L))
})

test_that("source_mspct", {

#  spct.l <- list(A = sun.spct, B = sun.spct)
#  my.mspct <- source_mspct(spct.l)
  my.mspct <- sun_evening.mspct

  peaks.mspct <- peaks(my.mspct)

  expect_equal(length(peaks.mspct), length(my.mspct))
  expect_equal(names(peaks.mspct[[1]]), c("w.length", "s.e.irrad"))
  expect_is(peaks.mspct[[1]], "source_spct")
  expect_is(peaks.mspct, "source_mspct")

  peaks.mspct <- peaks(my.mspct, unit.out = "photon")

  expect_equal(length(peaks.mspct), length(my.mspct))
  expect_equal(names(peaks.mspct[[1]]), c("w.length", "s.q.irrad"))
  expect_is(peaks.mspct[[1]], "source_spct")
  expect_is(peaks.mspct, "source_mspct")

})


context("valleys")

test_that("source_spct", {

  my.spct <- sun.spct

  valleys.spct <- valleys(sun.spct, span = NULL, strict = TRUE)
  expect_equal(nrow(valleys.spct), 0)
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct, span = NULL, strict = FALSE)
  expect_equal(nrow(valleys.spct), 14)
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct, span = NULL, strict = FALSE, global.threshold = -0.1)
  expect_equal(nrow(valleys.spct), 14)
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct, span = NULL, strict = FALSE, global.threshold = 0.1)
  expect_equal(nrow(valleys.spct), 0)
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct, span = NULL, strict = FALSE, global.threshold = -0.9)
  expect_equal(nrow(valleys.spct), 14)
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct, strict = FALSE, global.threshold = -0.1)
  expect_equal(nrow(valleys.spct), 12)
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct, strict = FALSE, global.threshold = 0.1)
  expect_equal(nrow(valleys.spct), 0)
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct, strict = FALSE, global.threshold = -0.9)
  expect_equal(nrow(valleys.spct), 87)
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct, strict = FALSE,
                          global.threshold = -0.9,
                          threshold.range = c(0, 0.82))
  expect_equal(nrow(valleys.spct), 87)
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct, strict = FALSE,
                          local.threshold = 0.1)
  expect_equal(nrow(valleys.spct), 16)
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct, span = 51)
  expect_equal(nrow(valleys.spct), 9)
  expect_equal(names(valleys.spct), c("w.length", "s.e.irrad"))
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct, span = 5)
  expect_equal(nrow(valleys.spct), 75)
  expect_equal(names(valleys.spct), c("w.length", "s.e.irrad"))
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct)
  expect_equal(nrow(valleys.spct), 75)
  expect_equal(names(valleys.spct), c("w.length", "s.e.irrad"))
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct)
  expect_equal(nrow(valleys.spct), 75)
  expect_equal(names(valleys.spct), c("w.length", "s.e.irrad"))
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct)
  expect_equal(nrow(valleys.spct), 75)
  expect_equal(names(valleys.spct), c("w.length", "s.e.irrad"))
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct, span = NULL,
                          unit.out = "photon")
  expect_equal(nrow(valleys.spct), 0)

  valleys.spct <- valleys(my.spct, unit.out = "photon")

  expect_equal(nrow(valleys.spct), 77)
  expect_equal(names(valleys.spct), c("w.length", "s.q.irrad"))
  expect_is(valleys.spct, "source_spct")

  my_thn.spct <- thin_wl(sun.spct)
  expect_silent(valleys(my_thn.spct))
  expect_silent(valleys(my_thn.spct, span = 5L))
  expect_silent(valleys(my_thn.spct, span = NULL))
  expect_warning(valleys(my_thn.spct, span = 7L))
  expect_warning(valleys(my_thn.spct, span = 101L))
})

test_that("source_mspct", {

  spct.l <- list(A = sun.spct, B = sun.spct)
  my.mspct <- source_mspct(spct.l)

  valleys.mspct <- valleys(my.mspct)

  expect_equal(length(valleys.mspct), length(my.mspct))
  expect_equal(names(valleys.mspct[[1]]), c("w.length", "s.e.irrad"))
  expect_is(valleys.mspct[[1]], "source_spct")
  expect_is(valleys.mspct, "source_mspct")

  valleys.mspct <- valleys(my.mspct, unit.out = "photon")

  expect_equal(length(valleys.mspct), length(my.mspct))
  expect_equal(names(valleys.mspct[[1]]), c("w.length", "s.q.irrad"))
  expect_is(valleys.mspct[[1]], "source_spct")
  expect_is(valleys.mspct, "source_mspct")

})

context("wls_at_target")

test_that("source_spct", {

  my.spct <- white_led.source_spct

  wls.spct <- wls_at_target(my.spct, interpolate = TRUE)
  expect_equal(wls.spct[["w.length"]], c(541.0686, 661.0016), tolerance = 1e-6)
  expect_equal(nrow(wls.spct), 2)
  expect_equal(names(wls.spct), c("w.length", "s.e.irrad"))
  expect_is(wls.spct, "source_spct")

  wls.spct <- wls_at_target(my.spct, idfactor = TRUE, interpolate = TRUE)
  expect_equal(wls.spct[["w.length"]], c(541.0686, 661.0016), tolerance = 1e-6)
  expect_equal(nrow(wls.spct), 2)
  expect_equal(names(wls.spct), c("w.length", "s.e.irrad", "target.idx"))
  expect_is(wls.spct, "source_spct")

  wls.spct <- wls_at_target(my.spct, idfactor = "TARGET", interpolate = TRUE)
  expect_equal(wls.spct[["w.length"]], c(541.0686, 661.0016), tolerance = 1e-5)
  expect_equal(nrow(wls.spct), 2)
  expect_equal(names(wls.spct), c("w.length", "s.e.irrad", "TARGET"))
  expect_is(wls.spct, "source_spct")

  wls.spct <- wls_at_target(my.spct, target = c("half.maximum", "half.range"))
  expect_equal(nrow(wls.spct), 4)
  expect_equal(names(wls.spct), c("w.length", "s.e.irrad", "target.idx"))
  expect_is(wls.spct[["target.idx"]], "factor")
  expect_equal(levels(wls.spct[["target.idx"]]), c("0.5max", "0.5range"))
  expect_is(wls.spct, "source_spct")

  wls.spct <- wls_at_target(my.spct, target = "half.maximum")
  expect_equal(nrow(wls.spct), 2)
  expect_equal(names(wls.spct), c("w.length", "s.e.irrad"))
  expect_is(wls.spct, "source_spct")

  expect_equal(wls_at_target(my.spct, target = "HM"), wls.spct)
  expect_equal(wls_at_target(my.spct, target = "0.5max"), wls.spct)
  expect_equal(wls_at_target(my.spct, target = "0.5 max"), wls.spct)

  wls.spct <- wls_at_target(my.spct, target = "half.range", interpolate = TRUE)
  expect_equal(wls.spct[["w.length"]], c(541.0686, 661.0016), tolerance = 1e-5)
  expect_equal(nrow(wls.spct), 2)
  expect_equal(names(wls.spct), c("w.length", "s.e.irrad"))
  expect_is(wls.spct, "source_spct")

  expect_equal(wls_at_target(my.spct, target = "HR", interpolate = TRUE), wls.spct)
  expect_equal(wls_at_target(my.spct, target = "0.5range", interpolate = TRUE), wls.spct)
  expect_equal(wls_at_target(my.spct, target = "0.5 range", interpolate = TRUE), wls.spct)

  wls.spct <- wls_at_target(my.spct, target = "half.range")
  expect_equal(wls_at_target(my.spct, target = "HR"), wls.spct)
  expect_equal(wls_at_target(my.spct, target = "0.5range"), wls.spct)
  expect_equal(wls_at_target(my.spct, target = "0.5 range"), wls.spct)

  expect_lt(max(abs(wls_at_target(my.spct, target = "0.5 range")$w.length -
                             wls_at_target(my.spct, target = "0.5 range",
                                           interpolate = TRUE)$w.length)), 0.225)

  wls.spct <- wls_at_target(my.spct, unit.out = "photon")
  expect_equal(nrow(wls.spct), 2)
  expect_equal(names(wls.spct), c("w.length", "s.q.irrad"))
  expect_is(wls.spct, "source_spct")

  wls.spct <- wls_at_target(my.spct, target = "half.maximum",
                            unit.out = "photon")
  expect_equal(nrow(wls.spct), 2)
  expect_equal(names(wls.spct), c("w.length", "s.q.irrad"))
  expect_is(wls.spct, "source_spct")

  wls.spct <- find_wls(my.spct)
  expect_equal(nrow(wls.spct), 0)
  expect_equal(names(wls.spct), c("w.length", "s.e.irrad"))
  expect_is(wls.spct, "source_spct")
})

test_that("source_mspct", {

  spct.l <- list(A = sun.spct, B = sun.spct)
  my.mspct <- source_mspct(spct.l)

  wls.mspct <- wls_at_target(my.mspct, 0.5)

  expect_equal(length(wls.mspct), length(my.mspct))
  expect_true(all(c("w.length", "s.e.irrad") %in% names(wls.mspct[[1]])))
  expect_is(wls.mspct[[1]], "source_spct")
  expect_is(wls.mspct, "source_mspct")

  wls.mspct <- wls_at_target(my.mspct, target = 1e-6, unit.out = "photon")

  expect_equal(length(wls.mspct), length(my.mspct))
  expect_true(all(c("w.length", "s.q.irrad") %in% names(wls.mspct[[1]])))
  expect_is(wls.mspct[[1]], "source_spct")
  expect_is(wls.mspct, "source_mspct")

  wls.mspct <- wls_at_target(my.mspct, target = 1e-6, interpolate = TRUE)

  expect_equal(length(wls.mspct), length(my.mspct))
  expect_equal(names(wls.mspct[[1]]), c("w.length", "s.e.irrad"))
  expect_is(wls.mspct[[1]], "source_spct")
  expect_is(wls.mspct, "source_mspct")

  wls.mspct <- wls_at_target(my.mspct, target = 1e-6, interpolate = TRUE, unit.out = "photon")

  expect_equal(length(wls.mspct), length(my.mspct))
  expect_equal(names(wls.mspct[[1]]), c("w.length", "s.q.irrad"))
  expect_is(wls.mspct[[1]], "source_spct")
  expect_is(wls.mspct, "source_mspct")

})

test_that("data.frame", {

  my.df <- data.frame(x = as.double(4:10), y = as.double(c(1:4, 3:1)))

  wls.df <- find_wls(my.df, 2.5, col.name.x = "x")

  expect_equal(nrow(wls.df), 2)
  expect_true(all(c("x", "y") %in% colnames(wls.df)))
  expect_equal(wls.df[["x"]], c(6, 9))
  expect_is(wls.df, "data.frame")
  expect_type(wls.df[["x"]], "double")
  expect_type(wls.df[["y"]], "double")

  wls.df <- find_wls(my.df, 2.5, col.name.x = "x", interpolate = TRUE)

  expect_equal(nrow(wls.df), 2)
  expect_true(all(c("x", "y") %in% colnames(wls.df)))
  expect_equal(wls.df[["x"]], c(5.5, 8.5))
  expect_is(wls.df, "data.frame")
  expect_type(wls.df[["x"]], "double")
  expect_type(wls.df[["y"]], "double")

  wls.df <- find_wls(my.df, "HM", col.name.x = "x")

  expect_equal(nrow(wls.df), 2)
  expect_true(all(c("x", "y") %in% colnames(wls.df)))
  expect_equal(wls.df[["x"]], c(5, 9))
  expect_is(wls.df, "data.frame")
  expect_type(wls.df[["x"]], "double")
  expect_type(wls.df[["y"]], "double")

  wls.df <- find_wls(my.df, "HM", col.name.x = "x", interpolate = TRUE)

  expect_equal(nrow(wls.df), 2)
  expect_true(all(c("x", "y") %in% colnames(wls.df)))
  expect_equal(wls.df[["x"]], c(5, 9))
  expect_is(wls.df, "data.frame")
  expect_type(wls.df[["x"]], "double")
  expect_type(wls.df[["y"]], "double")

  wls.df <- find_wls(my.df, "HR", col.name.x = "x")

  expect_equal(nrow(wls.df), 2)
  expect_true(all(c("x", "y") %in% colnames(wls.df)))
  expect_equal(wls.df[["x"]], c(6, 9))
  expect_is(wls.df, "data.frame")
  expect_type(wls.df[["x"]], "double")
  expect_type(wls.df[["y"]], "double")

  wls.df <- find_wls(my.df, "HR", col.name.x = "x", interpolate = TRUE)

  expect_equal(nrow(wls.df), 2)
  expect_true(all(c("x", "y") %in% colnames(wls.df)))
  expect_equal(wls.df[["x"]], c(5.5, 8.5))
  expect_is(wls.df, "data.frame")
  expect_type(wls.df[["x"]], "double")
  expect_type(wls.df[["y"]], "double")

})

context("spikes")

test_that("source_spct", {

  my.spct <- sun.spct

  spikes.spct <- spikes(sun.spct, max.spike.width = 2)
  expect_equal(nrow(spikes.spct), 2)

  spikes.spct <- spikes(sun.spct, max.spike.width = 5, z.threshold = 3.5)
  expect_equal(nrow(spikes.spct), 13)

  spikes.spct <- spikes(sun.spct, max.spike.width = 2)
  expect_equal(nrow(spikes.spct), 2)

  spikes.spct <- spikes(sun.spct, max.spike.width = 2)
  expect_equal(nrow(spikes.spct), 2)


  spikes.spct <- spikes(sun.spct)

  expect_equal(nrow(spikes.spct), 2)
  expect_equal(names(spikes.spct), c("w.length", "s.e.irrad"))
  expect_is(spikes.spct, "source_spct")

  spikes.spct <- spikes(sun.spct, max.spike.width = 2, unit.out = "photon")
  expect_equal(nrow(spikes.spct), 1)

  spikes.spct <- spikes(my.spct, unit.out = "photon")

  expect_equal(nrow(spikes.spct), 1)
  expect_equal(names(spikes.spct), c("w.length", "s.q.irrad"))
  expect_is(spikes.spct, "source_spct")

})

test_that("source_mspct", {

  # spct.l <- list(A = sun.spct, B = sun.spct)
  # my.mspct <- source_mspct(spct.l)
  my.mspct <- sun_evening.mspct

  spikes.mspct <- spikes(my.mspct)

  expect_equal(length(spikes.mspct), length(my.mspct))
  expect_equal(names(spikes.mspct[[1]]), c("w.length", "s.e.irrad"))
  expect_is(spikes.mspct[[1]], "source_spct")
  expect_is(spikes.mspct, "source_mspct")

  spikes.mspct <- spikes(my.mspct, unit.out = "photon")

  expect_equal(length(spikes.mspct), length(my.mspct))
  expect_equal(names(spikes.mspct[[1]]), c("w.length", "s.q.irrad"))
  expect_is(spikes.mspct[[1]], "source_spct")
  expect_is(spikes.mspct, "source_mspct")

})



