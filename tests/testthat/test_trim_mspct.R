library("photobiology")

context("trim_mspct")

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(A = my1.spct, B = my2.spct, C = my3.spct, D = my4.spct, E = my5.spct)
  my.mspct <- source_mspct(spct.l)

  trimmed.mspct <- trim_mspct(my.mspct, range = c(402:410))

  expect_equal(length(trimmed.mspct), length(my.mspct))
  expect_equal(min(trimmed.mspct[[1]]), 402)
  expect_equal(max(trimmed.mspct[[1]]), 410)
  expect_equal(min(trimmed.mspct[[5]]), 402)
  expect_equal(max(trimmed.mspct[[5]]), 410)
  expect_equal(min(trimmed.mspct[["A"]]), 402)
  expect_equal(max(trimmed.mspct[["A"]]), 410)
  expect_equal(min(trimmed.mspct[["E"]]), 402)
  expect_equal(max(trimmed.mspct[["E"]]), 410)
  expect_equal(min(trimmed.mspct[[1]][["s.e.irrad"]]), 1)
  expect_equal(max(trimmed.mspct[[1]][["s.e.irrad"]]), 1)
  expect_equal(min(trimmed.mspct[[5]][["s.e.irrad"]]), 5)
  expect_equal(max(trimmed.mspct[[5]][["s.e.irrad"]]), 5)

  trimmed.mspct <- trim_mspct(my.mspct, range = c(402:410), fill = 0)

  expect_equal(min(trimmed.mspct[[1]][["s.e.irrad"]]), 0)
  expect_equal(max(trimmed.mspct[[1]][["s.e.irrad"]]), 1)
  expect_equal(min(trimmed.mspct[[5]][["s.e.irrad"]]), 0)
  expect_equal(max(trimmed.mspct[[5]][["s.e.irrad"]]), 5)

  trim_mspct(my.mspct, range = c(402:410), fill = 0, byref = TRUE)

  expect_equal(min(my.mspct[[1]][["s.e.irrad"]]), 0)
  expect_equal(max(my.mspct[[1]][["s.e.irrad"]]), 1)
  expect_equal(min(my.mspct[[5]][["s.e.irrad"]]), 0)
  expect_equal(max(my.mspct[[5]][["s.e.irrad"]]), 5)

})

