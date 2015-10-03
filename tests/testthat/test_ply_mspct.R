# commented-out tests lead to infinite recursion, but the tested statements
# work as expected at the R command line!
#
library("photobiology")

context("msmsply")

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)

  result.mspct <- msmsply(my.mspct, `*`, e2 = 1)

  expect_equal(class(my.mspct), class(result.mspct))
  expect_equal(my.mspct, q2e(result.mspct, action = "replace"))

})

context("mslply")

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)

  result.lst <- mslply(my.mspct, `*`, e2 = 1)

  expect_equal("list", class(result.lst)[1])
  expect_equal(length(my.mspct), length(result.lst))
  expect_equal(names(my.mspct), names(result.lst))

})

context("msdply")

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)

  result.df <- msdply(my.mspct, `irrad`)

  expect_equal("tbl_df", class(result.df)[1])
  expect_equal(length(my.mspct), nrow(result.df))
  expect_equal(2, ncol(result.df))
  expect_equal(names(my.mspct), levels(result.df$spct.idx))
  expect_equal(1:5, result.df$irrad_Total / 10)

})

context("msaply")

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)

  result.ary <- msaply(my.mspct, min)

  expect_equal("numeric", class(result.ary)[1])
  expect_equal(length(my.mspct), length(result.ary))
  expect_true(is.numeric(result.ary))
  comment(result.ary) <- NULL
  expect_equal(rep(400, 5), result.ary)

  result.ary <- msaply(my.mspct, min, .drop = FALSE)

  expect_equal("matrix", class(result.ary)[1])
  expect_equal(length(my.mspct), nrow(result.ary))
  expect_true(is.numeric(result.ary))
  expect_equal(rep(400, 5), result.ary[ , 1])

  result.ary <- msaply(my.mspct, irrad)

  expect_equal("numeric", class(result.ary)[1])
  expect_equal(length(my.mspct), length(result.ary))
  expect_true(is.numeric(result.ary))
  expect_true(max(abs((1:5) * 10 -  result.ary)) < 5e12)

  result.ary <- msaply(my.mspct, irrad, .drop = FALSE)

  expect_equal("matrix", class(result.ary)[1])
  expect_equal(length(my.mspct), nrow(result.ary))
  expect_true(is.numeric(result.ary))
  expect_true(max(abs((1:5) * 10 -  result.ary)) < 5e12)

  result.ary <- msaply(my.mspct, range)

  expect_equal("matrix", class(result.ary)[1])
  expect_equal(length(my.mspct), nrow(result.ary))
  expect_true(is.numeric(result.ary))
  expect_equal(rep(400, 5), result.ary[ , 1])
  expect_equal(rep(410, 5), result.ary[ , 2])

  result.ary <- msaply(my.mspct, `range`, .drop = FALSE)

  expect_equal("matrix", class(result.ary)[1])
  expect_equal(length(my.mspct), nrow(result.ary))
  expect_true(is.numeric(result.ary))
  expect_equal(rep(400, 5), result.ary[ , 1])
  expect_equal(rep(410, 5), result.ary[ , 2])

  wb.lst <- list(a = waveband(c(400,405)), b = waveband(c(405,410)))

  result.ary <- msaply(my.mspct, irrad, wb.lst)

  expect_equal("matrix", class(result.ary)[1])
  expect_equal(length(my.mspct), nrow(result.ary))
  expect_true(is.numeric(result.ary))
  expect_equal((1:5) * 5, result.ary[ , 1])
  expect_equal((1:5) * 5, result.ary[ , 2])
  expect_true(max(abs((1:5) * 5 - result.ary[ , 1])) < 5e12)
  expect_true(max(abs((1:5) * 5 - result.ary[ , 2])) < 5e12)

})
