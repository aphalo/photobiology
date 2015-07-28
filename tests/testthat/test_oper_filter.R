library("photobiology")
context("filter_spct")

test_that("constructor T fraction", {

  my.spct <- filter_spct(w.length = 400:409, Tfr = 0.1)
  expect_equal(class(my.spct)[1:2], c("filter_spct", "generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 2)

  expect_error(filter_spct(w.length = 400:409, Tfr = -0.1, strict.range = TRUE))
  expect_error(filter_spct(w.length = 400:409, Tfr = 1.1, strict.range = TRUE))
  expect_warning(filter_spct(w.length = 400:409, Tfr = -0.1, strict.range = FALSE))
  expect_warning(filter_spct(w.length = 400:409, Tfr = -0.1))
  expect_warning(T2A(filter_spct(w.length = 400:409, Tfr = 0)))
  expect_equal(my.spct[["Tfr"]], rep(0.1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "Tfr"))
  expect_null(attr(my.spct, "time.unit", exact = TRUE))
})

test_that("constructor T percent", {

  my.spct <- filter_spct(w.length = 400:409, Tpc = 10)
  expect_equal(class(my.spct)[1:2], c("filter_spct", "generic_spct") )

  expect_warning(filter_spct(w.length = 400:409, Tpc = -0.1))
  expect_warning(filter_spct(w.length = 400:409, Tpc = 100.01))
  expect_equal(my.spct[["Tfr"]], rep(0.1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "Tfr"))
  expect_null(attr(my.spct, "time.unit", exact = TRUE))
})

test_that("constructor absorbance", {

  my.spct <- filter_spct(w.length = 400:409, A = 1)
  expect_equal(class(my.spct)[1:2], c("filter_spct", "generic_spct") )

  expect_warning(filter_spct(w.length = 400:409, A = -0.1))
  expect_equal(my.spct[["A"]], rep(1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "A"))
  expect_null(attr(my.spct, "time.unit", exact = TRUE))
})

test_that("oper default", {

  my.e.spct <- filter_spct(w.length = 400:409, Tfr = 0.1)
  my.2e.spct <- filter_spct(w.length = 400:409, Tfr = 0.2)

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

  my.e.spct <- filter_spct(w.length = 400:409, Tfr = 0.1)
  my.2e.spct <- filter_spct(w.length = 400:409, Tfr = 0.2)

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

  my.e.spct <- filter_spct(w.length = 400:409, A = 1)
  my.2e.spct <- filter_spct(w.length = 400:409, A = 2)

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

  my.e.spct <- filter_spct(w.length = 400:409, Tfr = 0.1)

  expect_equal(log10(my.e.spct)[["Tfr"]],  rep(log10(0.1), length.out = 10))
  expect_equal(log(my.e.spct)[["Tfr"]],  rep(log(0.1), length.out = 10))
  expect_equal(log(my.e.spct, 2)[["Tfr"]],  rep(log(0.1, 2), length.out = 10))
  expect_equal(exp(my.e.spct)[["Tfr"]],  rep(exp(0.1), length.out = 10))
  expect_equal(sqrt(my.e.spct)[["Tfr"]],  rep(sqrt(0.1), length.out = 10))

})

test_that("math absorbance", {

  my.e.spct <- filter_spct(w.length = 400:409, A = 1)

  options(photobiology.filter.qty = "absorbance")

  expect_equal(log10(my.e.spct)[["A"]],  rep(log10(1), length.out = 10))
  expect_equal(log(my.e.spct)[["A"]],  rep(log(1), length.out = 10))
  expect_equal(log(my.e.spct, 2)[["A"]],  rep(log(1, 2), length.out = 10))
  expect_equal(exp(my.e.spct)[["A"]],  rep(exp(1), length.out = 10))
  expect_equal(sqrt(my.e.spct)[["A"]],  rep(sqrt(1), length.out = 10))

  options(photobiology.filter.qty = NULL)

})

test_that("transmittance", {
  my.spct <- filter_spct(w.length = 300:700, Tfr = 1)

  transmittance.result <- 400
  expect_equal(as.numeric(transmittance(my.spct)), 1, tolerance = 1e-6)
  expect_equal(as.numeric(transmittance(my.spct, quantity = "total")),
               transmittance.result, tolerance = 1e-6)
  expect_equal(as.numeric(transmittance(my.spct, quantity = "average")), 1, tolerance = 1e-6)
  expect_equal(as.numeric(transmittance(my.spct, quantity = "mean")), 1, tolerance = 1e-6)
  expect_equal(sum(as.numeric(transmittance(my.spct, quantity = "total",
                                       w.band = split_bands(my.spct, length.out = 3)))),
               transmittance.result)
  expect_equal(sum(as.numeric(transmittance(my.spct, quantity = "average",
                                       w.band = split_bands(my.spct, length.out = 3)))), 3)
  expect_equal(sum(as.numeric(transmittance(my.spct, quantity = "average",
                                       w.band = split_bands(my.spct, length.out = 5)))), 5)

  expect_equal(sum(as.numeric(transmittance(my.spct, quantity = "relative",
                                       w.band = split_bands(my.spct, length.out = 3)))), 1)
  expect_equal(sum(as.numeric(transmittance(my.spct, quantity = "relative",
                                       w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(transmittance(my.spct, quantity = "contribution",
                                       w.band = split_bands(my.spct, length.out = 3)))), 1)
  expect_less_than(sum(as.numeric(transmittance(my.spct, quantity = "contribution",
                                           w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(transmittance(trim_spct(my.spct, range = c(400, 600)),
                                       quantity = "contribution",
                                       w.band = split_bands(c(400, 600), length.out = 3)))), 1)


})

test_that("absorptance", {
  my.spct <- filter_spct(w.length = 300:700, Tfr = 0.5, Tfr.type = "internal")

  absorptance.result <- 0.5 * 400
  expect_equal(as.numeric(absorptance(my.spct)), 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(absorptance(my.spct, quantity = "total")),
               absorptance.result, tolerance = 1e-6)
  expect_equal(as.numeric(absorptance(my.spct, quantity = "average")), 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(absorptance(my.spct, quantity = "mean")), 0.5, tolerance = 1e-6)
  expect_equal(sum(as.numeric(absorptance(my.spct, quantity = "total",
                                            w.band = split_bands(my.spct, length.out = 3)))),
               absorptance.result)
  expect_equal(sum(as.numeric(absorptance(my.spct, quantity = "average",
                                            w.band = split_bands(my.spct, length.out = 3)))), 3 * 0.5)
  expect_equal(sum(as.numeric(absorptance(my.spct, quantity = "average",
                                            w.band = split_bands(my.spct, length.out = 5)))), 5 * 0.5)

  expect_equal(sum(as.numeric(absorptance(my.spct, quantity = "relative",
                                            w.band = split_bands(my.spct, length.out = 3)))), 1)
  expect_equal(sum(as.numeric(absorptance(my.spct, quantity = "relative",
                                            w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(absorptance(my.spct, quantity = "contribution",
                                            w.band = split_bands(my.spct, length.out = 3)))), 1)
  expect_less_than(sum(as.numeric(absorptance(my.spct, quantity = "contribution",
                                                w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(absorptance(trim_spct(my.spct, range = c(400, 600)),
                                            quantity = "contribution",
                                            w.band = split_bands(c(400, 600), length.out = 3)))), 1)


})

test_that("absorbance", {
  my.spct <- filter_spct(w.length = 300:700, Tfr = 0.3162278, Tfr.type = "internal") # A = 0.5

  absorbance.result <- 0.5 * 400
  expect_equal(as.numeric(absorbance(my.spct)), 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(absorbance(my.spct, quantity = "total")),
               absorbance.result, tolerance = 1e-6)
  expect_equal(as.numeric(absorbance(my.spct, quantity = "average")), 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(absorbance(my.spct, quantity = "mean")), 0.5, tolerance = 1e-6)
  expect_equal(sum(as.numeric(absorbance(my.spct, quantity = "total",
                                          w.band = split_bands(my.spct, length.out = 3)))),
               absorbance.result, tolerance = 1e-6)
  expect_equal(sum(as.numeric(absorbance(my.spct, quantity = "average",
                                          w.band = split_bands(my.spct, length.out = 3)))),
               3 * 0.5, tolerance = 1e-6)
  expect_equal(sum(as.numeric(absorbance(my.spct, quantity = "average",
                                          w.band = split_bands(my.spct, length.out = 5)))),
               5 * 0.5, tolerance = 1e-6)

  expect_equal(sum(as.numeric(absorbance(my.spct, quantity = "relative",
                                          w.band = split_bands(my.spct, length.out = 3)))),
               1, tolerance = 1e-6)
  expect_equal(sum(as.numeric(absorbance(my.spct, quantity = "relative",
                                          w.band = split_bands(c(400, 600), length.out = 3)))),
               1, tolerance = 1e-6)
  expect_equal(sum(as.numeric(absorbance(my.spct, quantity = "contribution",
                                          w.band = split_bands(my.spct, length.out = 3)))),
               1, tolerance = 1e-6)
  expect_less_than(sum(as.numeric(absorbance(my.spct, quantity = "contribution",
                                              w.band = split_bands(c(400, 600), length.out = 3)))),
                   1, tolerance = 1e-6)
  expect_equal(sum(as.numeric(absorbance(trim_spct(my.spct, range = c(400, 600)),
                                          quantity = "contribution",
                                          w.band = split_bands(c(400, 600), length.out = 3)))),
               1, tolerance = 1e-6)


})
