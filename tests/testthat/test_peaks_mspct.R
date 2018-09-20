library("photobiology")

context("peaks")

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
  expect_equal(names(wls.mspct[[1]]), c("w.length", "s.e.irrad"))
  expect_is(wls.mspct[[1]], "source_spct")
  expect_is(wls.mspct, "source_mspct")

  wls.mspct <- wls_at_target(my.mspct, target = 1e-6, unit.out = "photon")

  expect_equal(length(wls.mspct), length(my.mspct))
  expect_equal(names(wls.mspct[[1]]), c("w.length", "s.q.irrad"))
  expect_is(wls.mspct[[1]], "source_spct")
  expect_is(wls.mspct, "source_mspct")

})


