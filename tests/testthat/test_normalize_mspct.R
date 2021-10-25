# commented-out tests lead to infinite recursion, but the tested statements
# work as expected at the R command line!
#
library("photobiology")

context("normalize_mspct")

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)

  normalized.mspct <- normalize(my.mspct)

  expect_lt(abs(max(normalized.mspct[[1]]$s.e.irrad) - 1), 1e-10)
  expect_lt(abs(max(normalized.mspct[[2]]$s.e.irrad) - 1), 1e-10)
  expect_lt(abs(max(normalized.mspct[[3]]$s.e.irrad) - 1), 1e-10)
  expect_lt(abs(max(normalized.mspct[[4]]$s.e.irrad) - 1), 1e-10)
  expect_lt(abs(max(normalized.mspct[[5]]$s.e.irrad) - 1), 1e-10)

  expect_equal(getNormalization(normalized.mspct[[1]])$norm.factor, 1)
  expect_equal(getNormalization(normalized.mspct[[2]])$norm.factor, 1/2)
  expect_equal(getNormalization(normalized.mspct[[3]])$norm.factor, 1/3)
  expect_equal(getNormalization(normalized.mspct[[4]])$norm.factor, 1/4)
  expect_equal(getNormalization(normalized.mspct[[5]])$norm.factor, 1/5)

  expect_equal(class(my.mspct), class(normalized.mspct))

})
