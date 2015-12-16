library("photobiology")

context("interpolate_mspct")

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(A = my1.spct, B = my2.spct, C = my3.spct, D = my4.spct, E = my5.spct)
  my.mspct <- source_mspct(spct.l)

  interpolated.mspct <- interpolate_mspct(my.mspct, length.out = 21)

  expect_equal(length(interpolated.mspct), length(my.mspct))
  expect_equal(min(interpolated.mspct[[1]]), 400)
  expect_equal(max(interpolated.mspct[[1]]), 410)
  expect_equal(min(interpolated.mspct[[5]]), 400)
  expect_equal(max(interpolated.mspct[[5]]), 410)
  expect_equal(min(interpolated.mspct[["A"]]), 400)
  expect_equal(max(interpolated.mspct[["A"]]), 410)
  expect_equal(min(interpolated.mspct[["E"]]), 400)
  expect_equal(max(interpolated.mspct[["E"]]), 410)
  expect_equal(min(interpolated.mspct[[1]][["s.e.irrad"]]), 1)
  expect_equal(max(interpolated.mspct[[1]][["s.e.irrad"]]), 1)
  expect_equal(min(interpolated.mspct[[5]][["s.e.irrad"]]), 5)
  expect_equal(max(interpolated.mspct[[5]][["s.e.irrad"]]), 5)

  interpolated.mspct <- interpolate_mspct(my.mspct, w.length.out = 398:409, fill = 0)

  expect_equal(length(interpolated.mspct), length(my.mspct))
  expect_equal(min(interpolated.mspct[[1]]), 398)
  expect_equal(max(interpolated.mspct[[1]]), 409)
  expect_equal(min(interpolated.mspct[[5]]), 398)
  expect_equal(max(interpolated.mspct[[5]]), 409)
  expect_equal(min(interpolated.mspct[["A"]]), 398)
  expect_equal(max(interpolated.mspct[["A"]]), 409)
  expect_equal(min(interpolated.mspct[["E"]]), 398)
  expect_equal(max(interpolated.mspct[["E"]]), 409)
  expect_equal(min(interpolated.mspct[[1]][["s.e.irrad"]]), 0)
  expect_equal(max(interpolated.mspct[[1]][["s.e.irrad"]]), 1)
  expect_equal(min(interpolated.mspct[[5]][["s.e.irrad"]]), 0)
  expect_equal(max(interpolated.mspct[[5]][["s.e.irrad"]]), 5)

})

