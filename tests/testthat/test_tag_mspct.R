library("photobiology")

context("tag_mspct")

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(A = my1.spct, B = my2.spct, C = my3.spct, D = my4.spct, E = my5.spct)
  my.mspct <- source_mspct(spct.l)

  tagged.mspct <- tag(my.mspct)

  expect_equal(length(tagged.mspct), length(my.mspct))
  expect_equal(names(tagged.mspct[[1]]), c("w.length", "s.e.irrad", "wl.color"))
  expect_equal(names(tagged.mspct[[5]]), c("w.length", "s.e.irrad", "wl.color"))
  expect_equal(min(tagged.mspct[[1]]), 400)
  expect_equal(max(tagged.mspct[[1]]), 410)
  expect_equal(min(tagged.mspct[[5]]), 400)
  expect_equal(max(tagged.mspct[[5]]), 410)
  expect_equal(min(tagged.mspct[["A"]]), 400)
  expect_equal(max(tagged.mspct[["A"]]), 410)
  expect_equal(min(tagged.mspct[["E"]]), 400)
  expect_equal(max(tagged.mspct[["E"]]), 410)
  expect_equal(min(tagged.mspct[[1]][["s.e.irrad"]]), 1)
  expect_equal(max(tagged.mspct[[1]][["s.e.irrad"]]), 1)
  expect_equal(min(tagged.mspct[[5]][["s.e.irrad"]]), 5)
  expect_equal(max(tagged.mspct[[5]][["s.e.irrad"]]), 5)

  tagged.mspct <- tag(my.mspct, w.band = waveband(c(402:410)))

  expect_equal(names(tagged.mspct[[1]]), c("w.length", "s.e.irrad", "wl.color", "wb.f"))
  expect_equal(names(tagged.mspct[[5]]), c("w.length", "s.e.irrad", "wl.color", "wb.f"))
  expect_equal(min(tagged.mspct[[1]]), 400)
  expect_equal(max(tagged.mspct[[1]]), 410)
  expect_equal(min(tagged.mspct[[5]]), 400)
  expect_equal(max(tagged.mspct[[5]]), 410)
  expect_equal(min(tagged.mspct[["A"]]), 400)
  expect_equal(max(tagged.mspct[["A"]]), 410)
  expect_equal(min(tagged.mspct[["E"]]), 400)
  expect_equal(max(tagged.mspct[["E"]]), 410)
  expect_equal(min(tagged.mspct[[1]][["s.e.irrad"]]), 1)
  expect_equal(max(tagged.mspct[[1]][["s.e.irrad"]]), 1)
  expect_equal(min(tagged.mspct[[5]][["s.e.irrad"]]), 5)
  expect_equal(max(tagged.mspct[[5]][["s.e.irrad"]]), 5)

  tag(my.mspct, w.band = waveband(c(402:410)), byref = TRUE)

  expect_equal(names(my.mspct[[1]]), c("w.length", "s.e.irrad", "wl.color", "wb.f"))
  expect_equal(names(my.mspct[[5]]), c("w.length", "s.e.irrad", "wl.color", "wb.f"))
  expect_equal(min(my.mspct[[1]]), 400)
  expect_equal(max(my.mspct[[1]]), 410)
  expect_equal(min(my.mspct[[5]]), 400)
  expect_equal(max(my.mspct[[5]]), 410)
  expect_equal(min(my.mspct[["A"]]), 400)
  expect_equal(max(my.mspct[["A"]]), 410)
  expect_equal(min(my.mspct[["E"]]), 400)
  expect_equal(max(my.mspct[["E"]]), 410)
  expect_equal(min(my.mspct[[1]][["s.e.irrad"]]), 1)
  expect_equal(max(my.mspct[[1]][["s.e.irrad"]]), 1)
  expect_equal(min(my.mspct[[5]][["s.e.irrad"]]), 5)
  expect_equal(max(my.mspct[[5]][["s.e.irrad"]]), 5)

  untagged.mspct <- untag(tagged.mspct)

#  expect_equal(untagged.mspct, my.mspct)

#  expect_equal(untagged.mspct[[1]], my.mspct[[1]])

#  expect_true(identical(untagged.mspct[[1]], my.mspct[[1]]))

  untag(tagged.mspct, byref = TRUE)

#  expect_equal(tagged.mspct, my.mspct)
  expect_equal(tagged.mspct, untagged.mspct)

  })

