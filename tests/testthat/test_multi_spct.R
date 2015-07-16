library(photobiology)
library(data.table)

context("multi_spct")

test_that("source_multi_spct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  my.mspct <- source_multi_spct(list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct))

  expect_equal(class(my.mspct)[1:2], c("source_multi_spct", "generic_multi_spct") )
  expect_equal(attr(my.mspct, "mspct.version", exact = TRUE), 1)
  expect_equal(attr(my.mspct, "ncol", exact = TRUE), 1)
  expect_equal(attr(my.mspct, "byrow", exact = TRUE), FALSE)

  my_named.mspct <- source_multi_spct(list(one = my1.spct,
                                           two = my2.spct,
                                           three = my3.spct,
                                           four = my4.spct,
                                           five = my5.spct))
  expect_equal(class(my_named.mspct)[1:2], c("source_multi_spct", "generic_multi_spct") )
  expect_equal(attr(my_named.mspct, "mspct.version", exact = TRUE), 1)
  expect_equal(attr(my_named.mspct, "ncol", exact = TRUE), 1)
  expect_equal(attr(my_named.mspct, "byrow", exact = TRUE), FALSE)

  expect_equal(length(my.mspct), 5)
  expect_equal(length(my_named.mspct), 5)
  expect_null(names(my.mspct))
  expect_equal(names(my_named.mspct), c("one", "two", "three", "four", "five"))
  expect_true(is.data.frame(irrad(my.mspct)))
  expect_true(is.data.frame(irrad(my_named.mspct)))
  expect_equal(irrad(my.mspct)[["spct.idx"]], factor(1:5))
  expect_equal(irrad(my_named.mspct)[["spct.idx"]],
               factor(c("one", "two", "three", "four", "five")))
  expect_equal(round(irrad(my.mspct)[["irrad_Total"]], 3), 1:5 * 10)
  expect_equal(round(irrad(my_named.mspct)[["irrad_Total"]], 3), 1:5 * 10)

  expect_equal(round(irrad(my.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
                           )[["irrad_A"]], 3), 1:5 * 5)
  expect_equal(round(irrad(my.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
                           )[["irrad_B"]], 3), 1:5 * 5)

  expect_equal(round(irrad(my_named.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
  )[["irrad_A"]], 3), 1:5 * 5)
  expect_equal(round(irrad(my_named.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
  )[["irrad_B"]], 3), 1:5 * 5)

  expect_equal(min(my.mspct)[["spct.idx"]], factor(1:5))
  expect_equal(max(my_named.mspct)[["spct.idx"]],
               factor(c("one", "two", "three", "four", "five")))
  expect_equal(min(my.mspct)[["min"]], rep(400, 5))
  expect_equal(min(my_named.mspct)[["min"]], rep(400, 5))
  expect_equal(max(my.mspct)[["max"]], rep(410, 5))
  expect_equal(max(my_named.mspct)[["max"]], rep(410, 5))

  expect_equal(range(my.mspct)[["range_1"]], rep(400, 5))
  expect_equal(range(my_named.mspct)[["range_1"]], rep(400, 5))
  expect_equal(range(my.mspct)[["range_2"]], rep(410, 5))
  expect_equal(range(my_named.mspct)[["range_2"]], rep(410, 5))

})

## need to add similar tests for other classes
## need to add additional methods
