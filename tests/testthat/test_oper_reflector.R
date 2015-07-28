library("photobiology")
context("reflector_spct")

test_that("constructor fraction", {

  my.spct <- reflector_spct(w.length = 400:409, Rfr = 0.1)
  expect_equal(class(my.spct)[1:2], c("reflector_spct", "generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 2)

  expect_warning(reflector_spct(w.length = 400:409, Rfr = -0.1))
  expect_warning(reflector_spct(w.length = 400:409, Rfr = 1.1))
  expect_equal(my.spct[["Rfr"]], rep(0.1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "Rfr"))
  expect_null(attr(my.spct, "time.unit", exact = TRUE))
})

test_that("constructor percent", {

  my.spct <- reflector_spct(w.length = 400:409, Rpc = 10)
  expect_equal(class(my.spct)[1:2], c("reflector_spct", "generic_spct") )

  expect_warning(reflector_spct(w.length = 400:409, Rpc = -0.1))
  expect_warning(reflector_spct(w.length = 400:409, Rpc = 100.01))
  expect_equal(my.spct[["Rfr"]], rep(0.1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "Rfr"))
  expect_null(attr(my.spct, "time.unit", exact = TRUE))
})

test_that("oper", {

  my.e.spct <- reflector_spct(w.length = 400:409, Rfr = 0.1)
  my.2e.spct <- reflector_spct(w.length = 400:409, Rfr = 0.2)

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

  my.e.spct <- reflector_spct(w.length = 400:409, Rfr = 0.1)
  my.2e.spct <- reflector_spct(w.length = 400:409, Rfr = 0.2)

  expect_equal(log10(my.e.spct)[["Rfr"]],  rep(log10(0.1), length.out = 10))
  expect_equal(log(my.e.spct)[["Rfr"]],  rep(log(0.1), length.out = 10))
  expect_equal(log(my.e.spct, 2)[["Rfr"]],  rep(log(0.1, 2), length.out = 10))
  expect_equal(exp(my.e.spct)[["Rfr"]],  rep(exp(0.1), length.out = 10))
  expect_equal(sqrt(my.e.spct)[["Rfr"]],  rep(sqrt(0.1), length.out = 10))

})

test_that("reflectance", {
  my.spct <- reflector_spct(w.length = 300:700, Rfr = 1)

  reflectance.result <- 400
  expect_equal(as.numeric(reflectance(my.spct)), 1, tolerance = 1e-6)
  expect_equal(as.numeric(reflectance(my.spct, quantity = "total")), reflectance.result, tolerance = 1e-6)
  expect_equal(as.numeric(reflectance(my.spct, quantity = "average")), 1, tolerance = 1e-6)
  expect_equal(as.numeric(reflectance(my.spct, quantity = "mean")), 1, tolerance = 1e-6)
  expect_equal(sum(as.numeric(reflectance(my.spct, quantity = "total",
                                          w.band = split_bands(my.spct, length.out = 3)))),
               reflectance.result)
  expect_equal(sum(as.numeric(reflectance(my.spct, quantity = "average",
                                          w.band = split_bands(my.spct, length.out = 3)))), 3)
  expect_equal(sum(as.numeric(reflectance(my.spct, quantity = "average",
                                          w.band = split_bands(my.spct, length.out = 5)))), 5)

  expect_equal(sum(as.numeric(reflectance(my.spct, quantity = "relative",
                                          w.band = split_bands(my.spct, length.out = 3)))), 1)
  expect_equal(sum(as.numeric(reflectance(my.spct, quantity = "relative",
                                          w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(reflectance(my.spct, quantity = "contribution",
                                          w.band = split_bands(my.spct, length.out = 3)))), 1)
  expect_less_than(sum(as.numeric(reflectance(my.spct, quantity = "contribution",
                                              w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(reflectance(trim_spct(my.spct, range = c(400, 600)),
                                          quantity = "contribution",
                                          w.band = split_bands(c(400, 600), length.out = 3)))), 1)


})
