library("photobiology")
context("reflector.spct")

test_that("constructor fraction", {

  my.spct <- reflector.spct(w.length = 400:409, Rfr = 0.1)

  expect_error(reflector.spct(w.length = 400:409, Rfr = -0.1))
  expect_error(reflector.spct(w.length = 400:409, Rfr = 1.1))
  expect_equal(my.spct[["Rfr"]], rep(0.1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "Rfr"))
  expect_null(attr(my.spct, "time.unit", exact = TRUE))
})

test_that("constructor percent", {

  my.spct <- reflector.spct(w.length = 400:409, Rpc = 10)

  expect_error(reflector.spct(w.length = 400:409, Rpc = -0.1))
  expect_error(reflector.spct(w.length = 400:409, Rpc = 100.01))
  expect_equal(my.spct[["Rfr"]], rep(0.1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "Rfr"))
  expect_null(attr(my.spct, "time.unit", exact = TRUE))
})

test_that("oper", {

  my.e.spct <- reflector.spct(w.length = 400:409, Rfr = 0.1)
  my.2e.spct <- reflector.spct(w.length = 400:409, Rfr = 0.2)

  expect_equal(my.e.spct + my.e.spct,  my.2e.spct)
  expect_equal(my.e.spct * 2, my.2e.spct)
  expect_equal(my.e.spct * 2L, my.2e.spct)
  expect_equal(my.2e.spct / 2, my.e.spct)
  expect_equal(my.2e.spct / 2L, my.e.spct)
  expect_warning(-my.e.spct)
  expect_equal(my.e.spct + 0.1, my.2e.spct)
  expect_equal(-my.2e.spct / -2, my.e.spct)
  expect_equal(-my.2e.spct / -2L, my.e.spct)
  expect_equal( 2 * my.e.spct, my.2e.spct)
  expect_equal( 1 / (2 / my.2e.spct), my.e.spct)
  expect_equal( 1 / my.e.spct, my.e.spct^-1)

})


test_that("math", {

  my.e.spct <- reflector.spct(w.length = 400:409, Rfr = 0.1)
  my.2e.spct <- reflector.spct(w.length = 400:409, Rfr = 0.2)

  expect_equal(log10(my.e.spct)[["Rfr"]],  rep(log10(0.1), length.out = 10))
  expect_equal(log(my.e.spct)[["Rfr"]],  rep(log(0.1), length.out = 10))
  expect_equal(log(my.e.spct, 2)[["Rfr"]],  rep(log(0.1, 2), length.out = 10))
  expect_equal(exp(my.e.spct)[["Rfr"]],  rep(exp(0.1), length.out = 10))
  expect_equal(sqrt(my.e.spct)[["Rfr"]],  rep(sqrt(0.1), length.out = 10))

})

