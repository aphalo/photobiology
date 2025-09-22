library("photobiology")
context("calibration_spct")

test_that("constructor", {

  empty.spct <- calibration_spct()
  expect_true(is.calibration_spct(empty.spct))
  expect_true(is.any_spct(empty.spct))
  expect_named(empty.spct, c("w.length", "irrad.mult"))
  expect_equal(nrow(empty.spct), 0L)

  my.spct <- calibration_spct(w.length = 400:409, irrad.mult = 1e-2)

  expect_equal(my.spct[["irrad.mult"]], rep(1e-2, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "irrad.mult"))
  expect_true(is.calibration_spct(my.spct))
  expect_true(is.any_spct(my.spct))
  expect_false(is.cps_spct(my.spct))
  expect_false(is.source_spct(my.spct))
  expect_false(is.filter_spct(my.spct))
  expect_false(is.reflector_spct(my.spct))
  expect_false(is.object_spct(my.spct))
  expect_false(is.response_spct(my.spct))
  expect_false(is.chroma_spct(my.spct))

  my.df <- data.frame(w.length = 400:409, irrad.mult = 1e-2)
  my.spct <- as.calibration_spct(my.df)

  expect_equal(my.spct[["irrad.mult"]], rep(1e-2, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "irrad.mult"))
  expect_true(is.calibration_spct(my.spct))
  expect_true(is.any_spct(my.spct))
})


test_that("oper", {

  my.spct <- calibration_spct(w.length = 400:409, irrad.mult = 1)
  my.2.spct <- calibration_spct(w.length = 400:409, irrad.mult = 2)

  expect_equal(class(my.spct)[1:2], c("calibration_spct", "generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 3)

  expect_equal(my.spct * 2, my.2.spct)
  expect_equal(my.spct * 2L, my.2.spct)
  expect_equal(my.2.spct / 2, my.spct)
  expect_equal(my.2.spct / 2L, my.spct)
  expect_warning(-my.spct * -2)
  expect_warning(-my.spct * -2L)
  expect_warning(-my.2.spct / -2)
  expect_warning(-my.2.spct / -2L)
  expect_equal(suppressWarnings(-my.spct * -2), my.2.spct)
  expect_equal(suppressWarnings(-my.spct * -2L), my.2.spct)
  expect_equal(suppressWarnings(-my.2.spct / -2), my.spct)
  expect_equal(suppressWarnings(-my.2.spct / -2L), my.spct)
  expect_equal( 2 * my.spct, my.2.spct)
  expect_equal( 1 / (2 / my.2.spct), my.spct)
  expect_equal( 1 / my.2.spct, my.2.spct^-1)
  expect_equal(my.2.spct %/% 2L, my.spct)
  expect_equal(my.2.spct %% 2L, my.spct %% 1L)

  my.cps.spct <- cps_spct(w.length = 400:409, cps = 100)
  my.source.spct <- source_spct(w.length = 400:409, s.e.irrad = 200)

  expect_equal(class(my.cps.spct * my.2.spct)[1:2],
               c("source_spct", "generic_spct") )
#  order of attributes is different
#  expect_equal(my.cps.spct * my.2.spct, my.source.spct)
#  expect_equal(my.2.spct * my.cps.spct, my.source.spct)

})


test_that("math", {

  my.spct <- calibration_spct(w.length = 400:409, irrad.mult = 1)

  expect_equal(log10(my.spct)[["irrad.mult"]],  rep(log10(1), length.out = 10))
  expect_equal(log(my.spct)[["irrad.mult"]],  rep(log(1), length.out = 10))
  expect_equal(log(my.spct, 2)[["irrad.mult"]],  rep(log(1, 2), length.out = 10))
  expect_equal(exp(my.spct)[["irrad.mult"]],  rep(exp(1), length.out = 10))
  expect_equal(sqrt(my.spct)[["irrad.mult"]],  rep(sqrt(1), length.out = 10))

})

