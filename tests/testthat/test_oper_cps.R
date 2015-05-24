library("photobiology")
context("cps_spct")

test_that("constructor", {

  my.spct <- cps_spct(w.length = 400:409, cps = 1)

  expect_equal(my.spct[["cps"]], rep(1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "cps"))
})


test_that("oper", {

  my.spct <- cps_spct(w.length = 400:409, cps = 1)
  my.2.spct <- cps_spct(w.length = 400:409, cps= 2)

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
})


test_that("math", {

  my.spct <- cps_spct(w.length = 400:409, cps = 1)

  expect_equal(log10(my.spct)[["cps"]],  rep(log10(1), length.out = 10))
  expect_equal(log(my.spct)[["cps"]],  rep(log(1), length.out = 10))
  expect_equal(log(my.spct, 2)[["cps"]],  rep(log(1, 2), length.out = 10))
  expect_equal(exp(my.spct)[["cps"]],  rep(exp(1), length.out = 10))
  expect_equal(sqrt(my.spct)[["cps"]],  rep(sqrt(1), length.out = 10))

})

