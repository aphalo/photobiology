library(photobiology)

context("multi_spct")

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  my.mspct <- source_mspct(list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct))

  expect_equal(class(my.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(attr(my.mspct, "mspct.version", exact = TRUE), 1)
  expect_equal(ncol(my.mspct), 1)
  expect_equal(attr(my.mspct, "byrow", exact = TRUE), FALSE)

  my_named.mspct <- source_mspct(list(one = my1.spct,
                                           two = my2.spct,
                                           three = my3.spct,
                                           four = my4.spct,
                                           five = my5.spct))
  expect_equal(class(my_named.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(attr(my_named.mspct, "mspct.version", exact = TRUE), 1)
  expect_equal(ncol(my_named.mspct), 1)
  expect_equal(attr(my_named.mspct, "byrow", exact = TRUE), FALSE)

  expect_equal(length(my.mspct), 5)
  expect_equal(names(my.mspct), paste("spct", 1:5, sep = "_"))
  expect_equal(length(my_named.mspct), 5)
  expect_equal(names(my_named.mspct),
               c("one", "two", "three", "four", "five"))
  expect_true(is.data.frame(irrad(my.mspct)))
  expect_true(is.data.frame(irrad(my_named.mspct)))
  expect_equal(irrad(my.mspct)[["spct.idx"]],
               factor(paste("spct", 1:5, sep = "_")))
  expect_equal(levels(irrad(my_named.mspct)[["spct.idx"]]),
               c("one", "two", "three", "four", "five"))
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

  expect_equal(min(my.mspct)[["spct.idx"]],
               factor(paste("spct", 1:5, sep = "_")))
  expect_equal(levels(max(my_named.mspct)[["spct.idx"]]),
               c("one", "two", "three", "four", "five"))
  expect_equal(min(my.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(min(my_named.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(max(my.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(max(my_named.mspct)[["max.wl"]], rep(410, 5))

  expect_equal(range(my.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(range(my_named.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(range(my.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(range(my_named.mspct)[["max.wl"]], rep(410, 5))

})

## need to add similar tests for other classes
## need to add additional methods
