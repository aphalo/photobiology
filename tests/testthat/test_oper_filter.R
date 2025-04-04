library("photobiology")
context("filter_spct")

test_that("constructor T fraction", {

  empty.spct <- filter_spct()
  expect_true(is.filter_spct(empty.spct))
  expect_true(is.any_spct(empty.spct))
  expect_named(empty.spct, c("w.length", "Tfr"))
  expect_equal(nrow(empty.spct), 0L)

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

  expect_true(is.filter_spct(my.spct))
  expect_true(is.any_spct(my.spct))
  expect_false(is.raw_spct(my.spct))
  expect_false(is.source_spct(my.spct))
  expect_false(is.cps_spct(my.spct))
  expect_false(is.reflector_spct(my.spct))
  expect_false(is.object_spct(my.spct))
  expect_false(is.response_spct(my.spct))
  expect_false(is.chroma_spct(my.spct))

  my.df <- data.frame(w.length = 400:409, Tfr = 0.1)
  my.spct <- as.filter_spct(my.df)

  expect_equal(class(my.spct)[1:2], c("filter_spct", "generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 2)
  expect_named(my.spct, c("w.length", "Tfr"))
  expect_true(is.filter_spct(my.spct))
  expect_true(is.any_spct(my.spct))

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

  expect_true(is.filter_spct(my.spct))
  expect_true(is.any_spct(my.spct))

})

test_that("constructor absorbance", {

  my.spct <- filter_spct(w.length = 400:409, A = 1)
  expect_equal(class(my.spct)[1:2], c("filter_spct", "generic_spct") )

  expect_warning(filter_spct(w.length = 400:409, A = -0.1))
  expect_equal(my.spct[["A"]], rep(1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "A"))
  expect_null(attr(my.spct, "time.unit", exact = TRUE))

  expect_true(is.filter_spct(my.spct))
  expect_true(is.any_spct(my.spct))
})

test_that("oper default", {

  my.e.spct <- filter_spct(w.length = 400:409, Tfr = 0.1)
  my.2e.spct <- filter_spct(w.length = 400:409, Tfr = 0.2)

  options(photobiology.filter.qty = NULL)

  expect_error(my.e.spct + my.e.spct)
#  expect_equal(suppressWarnings(my.e.spct + my.e.spct), NA)
  expect_equal(my.e.spct * 2, my.2e.spct)
  expect_equal(my.e.spct * 2L, my.2e.spct)
  expect_equal(my.2e.spct / 2, my.e.spct)
  expect_equal(my.2e.spct / 2L, my.e.spct)
  expect_warning(-my.e.spct)
  expect_equal( 2 * my.e.spct, my.2e.spct)
  expect_equal(suppressWarnings( 1 / (2 / my.2e.spct)), my.e.spct)
  expect_equal(suppressWarnings( 1 / my.e.spct),
               suppressWarnings(my.e.spct^-1))

})

test_that("oper transmittance", {

  my.e.spct <- filter_spct(w.length = 400:409, Tfr = 0.1)
  my.2e.spct <- filter_spct(w.length = 400:409, Tfr = 0.2)

  options(photobiology.filter.qty = "transmittance")

  expect_error(my.e.spct + my.e.spct)
#  expect_equal(suppressWarnings(my.e.spct + my.e.spct), NA)
  expect_equal(my.e.spct * 2, my.2e.spct)
  expect_equal(my.e.spct * 2L, my.2e.spct)
  expect_equal(my.2e.spct / 2, my.e.spct)
  expect_equal(my.2e.spct / 2L, my.e.spct)
  expect_warning(-my.e.spct)
  expect_equal( 2 * my.e.spct, my.2e.spct)
  expect_equal(suppressWarnings( 1 / (2 / my.2e.spct)), my.e.spct)
  expect_equal(suppressWarnings( 1 / my.e.spct),
               suppressWarnings( my.e.spct^-1))

  options(photobiology.filter.qty = NULL)

})

test_that("oper absorptance", {

  my.e.spct <- filter_spct(w.length = 400:409, Afr = 0.1, Tfr.type = "internal")
  my.2e.spct <- filter_spct(w.length = 400:409, Afr = 0.2, Tfr.type = "internal")

  options(photobiology.filter.qty = "absorptance")

  expect_error(my.e.spct + my.e.spct)
#  expect_equal(suppressWarnings(my.e.spct + my.e.spct), NA)
  expect_equal(my.e.spct * 2, my.2e.spct)
  expect_equal(my.e.spct * 2L, my.2e.spct)
  expect_equal(my.2e.spct / 2, my.e.spct)
  expect_equal(my.2e.spct / 2L, my.e.spct)
  expect_warning(-my.e.spct)
  expect_equal( 2 * my.e.spct, my.2e.spct)
  expect_equal(suppressWarnings( 1 / (2 / my.2e.spct)), my.e.spct)
  expect_equal(suppressWarnings( 1 / my.e.spct),
               suppressWarnings( my.e.spct^-1))

  options(photobiology.filter.qty = NULL)

})

test_that("oper absorbance", {

  my.e.spct <- filter_spct(w.length = 400:409, A = 1, Tfr.type = "internal")
  my.2e.spct <- filter_spct(w.length = 400:409, A = 2, Tfr.type = "internal")

  options(photobiology.filter.qty = "absorbance")

#  expect_equal(my.e.spct + my.e.spct, my.2e.spct) different attributes
  expect_error(my.e.spct * my.e.spct)
#  expect_equal(suppressWarnings(my.e.spct * my.e.spct), NA)
  expect_equal(my.e.spct * 2, my.2e.spct)
  expect_equal(my.e.spct * 2L, my.2e.spct)
  expect_equal(my.2e.spct / 2, my.e.spct)
  expect_equal(my.2e.spct / 2L, my.e.spct)
  expect_warning(-my.e.spct)
  expect_warning(-1 * my.e.spct)
  expect_warning(my.e.spct * -1)
  expect_equal(suppressWarnings(-my.e.spct),
               suppressWarnings(-1 * my.e.spct))
  expect_equal(suppressWarnings(-my.e.spct * -1), my.e.spct)
  expect_equal( 2 * my.e.spct, my.2e.spct)
  expect_equal( 1 / (2 / my.2e.spct), my.e.spct)
  expect_equal( 1 / my.e.spct, my.e.spct^-1)

  options(photobiology.filter.qty = NULL)

})

test_that("math default", {

  options(photobiology.filter.qty = NULL)

  my.e.spct <- filter_spct(w.length = 400:409, Tfr = 0.1, Tfr.type = "internal")

  expect_warning(log10(my.e.spct))
  expect_equal(suppressWarnings(log10(my.e.spct)[["Tfr"]]),
               rep(log10(0.1), length.out = 10))
  expect_warning(log(my.e.spct))
  expect_equal(suppressWarnings(log(my.e.spct)[["Tfr"]]),
               rep(log(0.1), length.out = 10))
  expect_warning(log(my.e.spct, 2))
  expect_equal(suppressWarnings(log(my.e.spct, 2)[["Tfr"]]),
               rep(log(0.1, 2), length.out = 10))
  expect_warning(exp(my.e.spct))
  expect_equal(suppressWarnings(exp(my.e.spct)[["Tfr"]]),
               rep(exp(0.1), length.out = 10))
  expect_equal(sqrt(my.e.spct)[["Tfr"]],  rep(sqrt(0.1), length.out = 10))

})

test_that("math absorbance", {

  my.e.spct <- filter_spct(w.length = 400:409, A = 1, Tfr.type = "internal")

  options(photobiology.filter.qty = "absorbance")

  expect_equal(log10(my.e.spct)[["A"]],  rep(log10(1), length.out = 10))
  expect_equal(log(my.e.spct)[["A"]],  rep(log(1), length.out = 10))
  expect_equal(log(my.e.spct, 2)[["A"]],  rep(log(1, 2), length.out = 10))
  expect_equal(exp(my.e.spct)[["A"]],  rep(exp(1), length.out = 10))
  expect_equal(sqrt(my.e.spct)[["A"]],  rep(sqrt(1), length.out = 10))

  options(photobiology.filter.qty = NULL)

})

test_that("transmittance", {
  my.spct <- filter_spct(w.length = 300:700, Tfr = 1, Tfr.type = "internal")

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
  expect_lt(sum(as.numeric(transmittance(my.spct, quantity = "contribution",
                                           w.band = split_bands(c(400, 600), length.out = 3)))), 0.5)
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
  expect_equal(as.numeric(absorptance(my.spct, quantity = "average")), 0.5)
  expect_equal(as.numeric(absorptance(my.spct, quantity = "mean")), 0.5)
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
  expect_lt(sum(as.numeric(absorptance(my.spct, quantity = "contribution",
                                                w.band = split_bands(c(400, 600), length.out = 3)))), 0.5)
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
  expect_equal(sum(as.numeric(absorbance(my.spct, quantity = "contribution",
                                              w.band = split_bands(c(400, 600), length.out = 3)))),
                   0.5, tolerance = 1e-6)
  expect_equal(sum(as.numeric(absorbance(trim_spct(my.spct, range = c(400, 600)),
                                          quantity = "contribution",
                                          w.band = split_bands(c(400, 600), length.out = 3)))),
               1, tolerance = 1e-6)


})

test_that("Tfr_ratio", {
  uvb.wb <- waveband(c(280,315), wb.name = "UVB")
  blue.wb <- waveband(c(400,500), wb.name = "B")

  Tfr.ratio.result <- 0.01053392
  expect_equal(
    as.numeric(Tfr_ratio(polyester.spct, uvb.wb, blue.wb)),
    Tfr.ratio.result, tolerance = 1e-6)
  expect_equal(
    as.numeric(Tfr_ratio(polyester.spct, uvb.wb, blue.wb, quantity = "mean")),
    Tfr.ratio.result, tolerance = 1e-6)
  expect_equal(
    as.numeric(Tfr_ratio(polyester.spct, uvb.wb, blue.wb, quantity = "average")),
    Tfr.ratio.result, tolerance = 1e-6)

  Tfr.ratio.result <- 0.003686873
  expect_equal(
    as.numeric(Tfr_ratio(polyester.spct, uvb.wb, blue.wb, quantity = "total")),
    Tfr.ratio.result, tolerance = 1e-6)

  expect_error(Tfr_ratio(polyester.spct, uvb.wb, blue.wb, quantity = "bad argument"))
  expect_warning(Tfr_ratio(sun.spct, uvb.wb, blue.wb))

  expect_named(
    Tfr_ratio(polyester.spct, uvb.wb, blue.wb),
    "UVB:B[Tfr(wl):Tfr(wl)]")

  Tfr.fraction.result <- 0.01042412
  expect_equal(
    as.numeric(Tfr_fraction(polyester.spct, uvb.wb, blue.wb)),
    Tfr.fraction.result, tolerance = 1e-6)
  expect_equal(
    as.numeric(Tfr_fraction(polyester.spct, uvb.wb, blue.wb, quantity = "mean")),
    Tfr.fraction.result, tolerance = 1e-6)
  expect_equal(
    as.numeric(Tfr_fraction(polyester.spct, uvb.wb, blue.wb, quantity = "average")),
    Tfr.fraction.result, tolerance = 1e-6)

  Tfr.fraction.result <- 0.00367333
  expect_equal(
    as.numeric(Tfr_fraction(polyester.spct, uvb.wb, blue.wb, quantity = "total")),
    Tfr.fraction.result, tolerance = 1e-6)

  expect_error(Tfr_fraction(polyester.spct, uvb.wb, blue.wb, quantity = "bad argument"))
  expect_warning(Tfr_fraction(sun.spct, uvb.wb, blue.wb))

  expect_named(
    Tfr_fraction(polyester.spct, uvb.wb, blue.wb),
    "UVB:(UVB+B)[Tfr(wl):Tfr(wl)]")

  Tfr.normdiff.result <- 0.9791518
  expect_equal(
    as.numeric(Tfr_normdiff(polyester.spct, blue.wb, uvb.wb)),
    Tfr.normdiff.result, tolerance = 1e-6)
  expect_equal(
    as.numeric(Tfr_normdiff(polyester.spct, blue.wb, uvb.wb, quantity = "mean")),
    Tfr.normdiff.result, tolerance = 1e-6)
  expect_equal(
    as.numeric(Tfr_normdiff(polyester.spct, blue.wb, uvb.wb, quantity = "average")),
    Tfr.normdiff.result, tolerance = 1e-6)

  Tfr.normdiff.result <- 0.9926533
  expect_equal(
    as.numeric(Tfr_normdiff(polyester.spct, blue.wb, uvb.wb, quantity = "total")),
    Tfr.normdiff.result, tolerance = 1e-6)

  expect_error(Tfr_normdiff(polyester.spct, blue.wb, uvb.wb, quantity = "bad argument"))
  expect_warning(Tfr_normdiff(sun.spct, blue.wb, uvb.wb))

  expect_named(
    Tfr_normdiff(polyester.spct, blue.wb, uvb.wb),
    "(B-UVB):(B+UVB)[Tfr(wl):Tfr(wl)]")
})

