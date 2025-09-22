library("photobiology")
context("raw_spct")

test_that("constructor", {

  empty.spct <- raw_spct()
  expect_true(is.raw_spct(empty.spct))
  expect_true(is.any_spct(empty.spct))
  expect_named(empty.spct, c("w.length", "counts"))
  expect_equal(nrow(empty.spct), 0L)

  my.spct <- raw_spct(w.length = 400:409, counts = 64000)

  expect_equal(my.spct[["counts"]], rep(64000, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "counts"))
  expect_true(is.raw_spct(my.spct))
  expect_true(is.any_spct(my.spct))
  expect_false(is.cps_spct(my.spct))
  expect_false(is.source_spct(my.spct))
  expect_false(is.filter_spct(my.spct))
  expect_false(is.reflector_spct(my.spct))
  expect_false(is.object_spct(my.spct))
  expect_false(is.response_spct(my.spct))
  expect_false(is.chroma_spct(my.spct))

  my.df <- data.frame(w.length = 400:409, counts = 64000)
  my.spct <- as.raw_spct(my.df)

  expect_equal(my.spct[["counts"]], rep(64000, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "counts"))
  expect_true(is.raw_spct(my.spct))
  expect_true(is.any_spct(my.spct))
})


test_that("oper", {

  my.spct <- raw_spct(w.length = 400:409, counts = 1)
  my.2.spct <- raw_spct(w.length = 400:409, counts = 2)

  expect_equal(class(my.spct)[1:2], c("raw_spct", "generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 3)

  expect_equal(my.spct + my.spct,  my.2.spct)
  expect_equal(my.spct * 2, my.2.spct)
  expect_equal(my.spct * 2L, my.2.spct)
  expect_equal(my.2.spct / 2, my.spct)
  expect_equal(my.2.spct / 2L, my.spct)
  expect_equal(-my.spct * -2, my.2.spct)
  expect_equal(-my.spct * -2L, my.2.spct)
  expect_equal(-my.2.spct / -2, my.spct)
  expect_equal(-my.2.spct / -2L, my.spct)
  expect_equal( 2 * my.spct, my.2.spct)
  expect_equal( 1 / (2 / my.2.spct), my.spct)
  expect_equal( 1 / my.spct, my.spct^-1)
  expect_equal(my.2.spct %/% 2L, my.spct)
  expect_equal(my.2.spct %% 2L, my.spct %% 1L)
})


test_that("math", {

  my.spct <- raw_spct(w.length = 400:409, counts = 1)

  expect_equal(log10(my.spct)[["counts"]],  rep(log10(1), length.out = 10))
  expect_equal(log(my.spct)[["counts"]],  rep(log(1), length.out = 10))
  expect_equal(log(my.spct, 2)[["counts"]],  rep(log(1, 2), length.out = 10))
  expect_equal(exp(my.spct)[["counts"]],  rep(exp(1), length.out = 10))
  expect_equal(sqrt(my.spct)[["counts"]],  rep(sqrt(1), length.out = 10))

})

