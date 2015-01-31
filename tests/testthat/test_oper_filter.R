library("photobiology")
context("filter.spct")

test_that("constructor T fraction", {

  my.spct <- filter.spct(w.length = 400:409, Tfr = 0.1)

  expect_warning(filter.spct(w.length = 400:409, Tfr = -0.1))
  expect_warning(filter.spct(w.length = 400:409, Tfr = 1.1))
  expect_equal(my.spct[["Tfr"]], rep(0.1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "Tfr"))
  expect_null(attr(my.spct, "time.unit", exact = TRUE))
})

test_that("constructor T percent", {

  my.spct <- filter.spct(w.length = 400:409, Tpc = 10)

  expect_warning(filter.spct(w.length = 400:409, Tpc = -0.1))
  expect_warning(filter.spct(w.length = 400:409, Tpc = 100.01))
  expect_equal(my.spct[["Tfr"]], rep(0.1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "Tfr"))
  expect_null(attr(my.spct, "time.unit", exact = TRUE))
})

test_that("constructor absorbance", {

  my.spct <- filter.spct(w.length = 400:409, A = 1)

  expect_warning(filter.spct(w.length = 400:409, A = -0.1))
  expect_equal(my.spct[["A"]], rep(1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "A"))
  expect_null(attr(my.spct, "time.unit", exact = TRUE))
})

test_that("oper default", {

  my.e.spct <- filter.spct(w.length = 400:409, Tfr = 0.1)
  my.2e.spct <- filter.spct(w.length = 400:409, Tfr = 0.2)

  options(photobiology.filter.qty = NULL)

  expect_warning(my.e.spct + my.e.spct)
  expect_equal(my.e.spct + my.e.spct, NA)
  expect_equal(my.e.spct * 2, my.2e.spct)
  expect_equal(my.e.spct * 2L, my.2e.spct)
  expect_equal(my.2e.spct / 2, my.e.spct)
  expect_equal(my.2e.spct / 2L, my.e.spct)
  expect_warning(-my.e.spct)
  expect_equal( 2 * my.e.spct, my.2e.spct)
  expect_equal( 1 / (2 / my.2e.spct), my.e.spct)
  expect_equal( 1 / my.e.spct, my.e.spct^-1)

})

test_that("oper transmittance", {

  my.e.spct <- filter.spct(w.length = 400:409, Tfr = 0.1)
  my.2e.spct <- filter.spct(w.length = 400:409, Tfr = 0.2)

  options(photobiology.filter.qty = "transmittance")

  expect_warning(my.e.spct + my.e.spct)
  expect_equal(my.e.spct + my.e.spct, NA)
  expect_equal(my.e.spct * 2, my.2e.spct)
  expect_equal(my.e.spct * 2L, my.2e.spct)
  expect_equal(my.2e.spct / 2, my.e.spct)
  expect_equal(my.2e.spct / 2L, my.e.spct)
  expect_warning(-my.e.spct)
  expect_equal( 2 * my.e.spct, my.2e.spct)
  expect_equal( 1 / (2 / my.2e.spct), my.e.spct)
  expect_equal( 1 / my.e.spct, my.e.spct^-1)

  options(photobiology.filter.qty = NULL)

})

test_that("oper absorbance", {

  my.e.spct <- filter.spct(w.length = 400:409, A = 1)
  my.2e.spct <- filter.spct(w.length = 400:409, A = 2)

  options(photobiology.filter.qty = "absorbance")

  expect_equal(my.e.spct + my.e.spct, my.2e.spct)
  expect_warning(my.e.spct * my.e.spct)
  expect_equal(my.e.spct * my.e.spct, NA)
  expect_equal(my.e.spct * 2, my.2e.spct)
  expect_equal(my.e.spct * 2L, my.2e.spct)
  expect_equal(my.2e.spct / 2, my.e.spct)
  expect_equal(my.2e.spct / 2L, my.e.spct)
  expect_equal(-my.e.spct, -1 * my.e.spct)
  expect_equal(-my.e.spct * -1, my.e.spct)
  expect_equal( 2 * my.e.spct, my.2e.spct)
  expect_equal( 1 / (2 / my.2e.spct), my.e.spct)
  expect_equal( 1 / my.e.spct, my.e.spct^-1)

  options(photobiology.filter.qty = NULL)

})

test_that("math default", {

  my.e.spct <- filter.spct(w.length = 400:409, Tfr = 0.1)
  my.2e.spct <- filter.spct(w.length = 400:409, Tfr = 0.2)

  expect_equal(log10(my.e.spct)[["Tfr"]],  rep(log10(0.1), length.out = 10))
  expect_equal(log(my.e.spct)[["Tfr"]],  rep(log(0.1), length.out = 10))
  expect_equal(log(my.e.spct, 2)[["Tfr"]],  rep(log(0.1, 2), length.out = 10))
  expect_equal(exp(my.e.spct)[["Tfr"]],  rep(exp(0.1), length.out = 10))
  expect_equal(sqrt(my.e.spct)[["Tfr"]],  rep(sqrt(0.1), length.out = 10))

})

test_that("math absorbance", {

  my.e.spct <- filter.spct(w.length = 400:409, Tfr = 0.1)
  my.2e.spct <- filter.spct(w.length = 400:409, Tfr = 0.2)

  options(photobiology.filter.qty = "absorbance")

  expect_equal(log10(my.e.spct)[["Tfr"]],  rep(log10(0.1), length.out = 10))
  expect_equal(log(my.e.spct)[["Tfr"]],  rep(log(0.1), length.out = 10))
  expect_equal(log(my.e.spct, 2)[["Tfr"]],  rep(log(0.1, 2), length.out = 10))
  expect_equal(exp(my.e.spct)[["Tfr"]],  rep(exp(0.1), length.out = 10))
  expect_equal(sqrt(my.e.spct)[["Tfr"]],  rep(sqrt(0.1), length.out = 10))

  options(photobiology.filter.qty = NULL)

})
