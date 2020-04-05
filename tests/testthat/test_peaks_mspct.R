library("photobiology")

context("peaks")

test_that("source_spct", {

  my.spct <- sun.spct

  peaks.spct <- peaks(sun.spct, span = NULL, strict = TRUE)
  expect_equal(nrow(peaks.spct), 1)

  peaks.spct <- peaks(sun.spct, span = NULL, strict = FALSE)
  expect_equal(nrow(peaks.spct), 1)

  peaks.spct <- peaks(sun.spct, span = NULL, strict = TRUE, ignore_threshold = 0.9)
  expect_equal(nrow(peaks.spct), 1)

  peaks.spct <- peaks(sun.spct, span = NULL, strict = FALSE, ignore_threshold = -0.1)
  expect_equal(nrow(peaks.spct), 0)

  peaks.spct <- peaks(sun.spct)

  expect_equal(nrow(peaks.spct), 76)
  expect_equal(names(peaks.spct), c("w.length", "s.e.irrad"))
  expect_is(peaks.spct, "source_spct")

  peaks.spct <- peaks(sun.spct, span = NULL, unit.out = "photon")
  expect_equal(nrow(peaks.spct), 1)

  peaks.spct <- peaks(my.spct, unit.out = "photon")

  expect_equal(nrow(peaks.spct), 77)
  expect_equal(names(peaks.spct), c("w.length", "s.q.irrad"))
  expect_is(peaks.spct, "source_spct")

})

test_that("source_mspct", {

  spct.l <- list(A = sun.spct, B = sun.spct)
  my.mspct <- source_mspct(spct.l)

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
  expect_equal(nrow(valleys.spct), 1)

  valleys.spct <- valleys(sun.spct, span = NULL, strict = FALSE)
  expect_equal(nrow(valleys.spct), 14)

  valleys.spct <- valleys(sun.spct, span = NULL, ignore_threshold = -0.1)
  expect_equal(nrow(valleys.spct), 0)

  valleys.spct <- valleys(sun.spct, span = NULL, ignore_threshold = 0.1)
  expect_equal(nrow(valleys.spct), 1)

  valleys.spct <- valleys(sun.spct)

  expect_equal(nrow(valleys.spct), 75)
  expect_equal(names(valleys.spct), c("w.length", "s.e.irrad"))
  expect_is(valleys.spct, "source_spct")

  valleys.spct <- valleys(sun.spct, span = NULL,
                          unit.out = "photon")
  expect_equal(nrow(valleys.spct), 1)

  valleys.spct <- valleys(my.spct, unit.out = "photon")

  expect_equal(nrow(valleys.spct), 77)
  expect_equal(names(valleys.spct), c("w.length", "s.q.irrad"))
  expect_is(valleys.spct, "source_spct")

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

  spct.l <- list(A = sun.spct, B = sun.spct)
  my.mspct <- source_mspct(spct.l)

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



