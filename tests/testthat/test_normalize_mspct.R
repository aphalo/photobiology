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

  expect_equal(normalized.mspct[[1]], normalized.mspct[[2]])
  expect_equal(normalized.mspct[[1]], normalized.mspct[[3]])
  expect_equal(normalized.mspct[[1]], normalized.mspct[[4]])
  expect_equal(normalized.mspct[[1]], normalized.mspct[[5]])

  expect_equal(class(my.mspct), class(normalized.mspct))

})
