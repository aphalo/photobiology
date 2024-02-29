# commented-out tests lead to infinite recursion, but the tested statements
# work as expected at the R command line!
#
library("photobiology")

context("fscale_mspct")

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)

  scaled.mspct <- fscale(my.mspct)

  expect_equal(scaled.mspct[[1]][["s.e.irrad"]], scaled.mspct[[2]][["s.e.irrad"]])
  expect_equal(scaled.mspct[[1]][["s.e.irrad"]], scaled.mspct[[3]][["s.e.irrad"]])
  expect_equal(scaled.mspct[[1]][["s.e.irrad"]], scaled.mspct[[4]][["s.e.irrad"]])
  expect_equal(scaled.mspct[[1]][["s.e.irrad"]], scaled.mspct[[5]][["s.e.irrad"]])

  expect_equal(class(my.mspct), class(scaled.mspct))

  expect_equal(getScaled(scaled.mspct[[1]])$multiplier, 1)
  expect_equal(getScaled(scaled.mspct[[2]])$multiplier, 1/2)
  expect_equal(getScaled(scaled.mspct[[3]])$multiplier, 1/3)
  expect_equal(getScaled(scaled.mspct[[4]])$multiplier, 1/4)
  expect_equal(getScaled(scaled.mspct[[5]])$multiplier, 1/5)

  expect_equal(getScaled(scaled.mspct[[1]])$f, "mean")
  expect_equal(getScaled(scaled.mspct[[1]])$target, 1)
  expect_equal(getScaled(scaled.mspct[[1]])$range, c(400, 410))

  expect_equal(getScaled(scaled.mspct[[2]])$f, "mean")
  expect_equal(getScaled(scaled.mspct[[2]])$target, 1)
  expect_equal(getScaled(scaled.mspct[[2]])$range, c(400, 410))

  expect_equal(getScaled(scaled.mspct[[3]])$f, "mean")
  expect_equal(getScaled(scaled.mspct[[3]])$target, 1)
  expect_equal(getScaled(scaled.mspct[[3]])$range, c(400, 410))

  expect_equal(getScaled(scaled.mspct[[4]])$f, "mean")
  expect_equal(getScaled(scaled.mspct[[4]])$target, 1)
  expect_equal(getScaled(scaled.mspct[[4]])$range, c(400, 410))

  expect_equal(getScaled(scaled.mspct[[5]])$f, "mean")
  expect_equal(getScaled(scaled.mspct[[5]])$target, 1)
  expect_equal(getScaled(scaled.mspct[[5]])$range, c(400, 410))

})

test_that("scaling of long source_spct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)
  my.spct <- rbindspct(my.mspct, idfactor = "test.id")

  scaled.spct <- fscale(my.spct)

  expect_equal(getIdFactor(my.spct), getIdFactor(scaled.spct))
  expect_equal(colnames(my.spct), colnames(scaled.spct))
  expect_equal(getMultipleWl(my.spct), getMultipleWl(scaled.spct))
  expect_equal(nrow(my.spct), nrow(scaled.spct))
  expect_equal(my.spct$w.length, scaled.spct$w.length)
  expect_equal(class(my.spct), class(scaled.spct))

})
