library("photobiology")

context("peaks")

test_that("source_mspct", {

  spct.l <- list(A = sun.spct, B = sun.spct)
  my.mspct <- source_mspct(spct.l)

  peaks.mspct <- peaks(my.mspct)

  expect_equal(length(peaks.mspct), length(my.mspct))
  expect_equal(names(peaks.mspct[[1]]), c("w.length", "s.e.irrad"))

  peaks.mspct <- peaks(my.mspct, unit.out = "photon")

  expect_equal(length(peaks.mspct), length(my.mspct))
  expect_equal(names(peaks.mspct[[1]]), c("w.length", "s.q.irrad"))

})

