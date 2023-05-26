library("photobiology")
context("reflector_spct")

test_that("constructor fraction", {

  empty.spct <- reflector_spct()
  expect_true(is.reflector_spct(empty.spct))
  expect_true(is.any_spct(empty.spct))
  expect_named(empty.spct, c("w.length", "Rfr"))
  expect_equal(nrow(empty.spct), 0L)

  my.spct <- reflector_spct(w.length = 400:409, Rfr = 0.1)
  expect_equal(class(my.spct)[1:2], c("reflector_spct", "generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 2)

  expect_warning(reflector_spct(w.length = 400:409, Rfr = -0.1))
  expect_warning(reflector_spct(w.length = 400:409, Rfr = 1.1))
  expect_equal(my.spct[["Rfr"]], rep(0.1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "Rfr"))
  expect_null(attr(my.spct, "time.unit", exact = TRUE))

  expect_true(is.reflector_spct(my.spct))
  expect_true(is.any_spct(my.spct))
  expect_false(is.raw_spct(my.spct))
  expect_false(is.source_spct(my.spct))
  expect_false(is.cps_spct(my.spct))
  expect_false(is.filter_spct(my.spct))
  expect_false(is.object_spct(my.spct))
  expect_false(is.response_spct(my.spct))
  expect_false(is.chroma_spct(my.spct))

  my.df <- data.frame(w.length = 400:409, Rfr = 0.1)
  my.spct <- as.reflector_spct(my.df)

  expect_equal(class(my.spct)[1:2], c("reflector_spct", "generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 2)
  expect_named(my.spct, c("w.length", "Rfr"))
  expect_true(is.reflector_spct(my.spct))
  expect_true(is.any_spct(my.spct))

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

  expect_true(is.reflector_spct(my.spct))
  expect_true(is.any_spct(my.spct))

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
  expect_equal(suppressWarnings(-my.2e.spct / -2), my.e.spct)
  expect_equal(suppressWarnings(-my.2e.spct / -2L), my.e.spct)
  expect_equal(2 * my.e.spct, my.2e.spct)
  expect_equal(suppressWarnings( 1 / (2 / my.2e.spct)), my.e.spct)
  expect_warning(1 / my.e.spct)
  expect_warning(my.e.spct^-1)
  expect_equal(suppressWarnings( 1 / my.e.spct),
               suppressWarnings(my.e.spct^-1))
  expect_equal(my.2e.spct %/% 2L, my.e.spct %/% 1L)
  expect_equal(my.2e.spct %% 2L / 2, my.e.spct %% 1L)

})


test_that("math", {

  my.e.spct <- reflector_spct(w.length = 400:409, Rfr = 0.1)
  my.2e.spct <- reflector_spct(w.length = 400:409, Rfr = 0.2)

  expect_warning(log10(my.e.spct))
  expect_equal(suppressWarnings(log10(my.e.spct)[["Rfr"]]),
               rep(log10(0.1), length.out = 10))
  expect_warning(log(my.e.spct))
  expect_equal(suppressWarnings(log(my.e.spct)[["Rfr"]]),
               rep(log(0.1), length.out = 10))
  expect_warning(log(my.e.spct, 2))
  expect_equal(suppressWarnings(log(my.e.spct, 2)[["Rfr"]]),
               rep(log(0.1, 2), length.out = 10))
  expect_warning(exp(my.e.spct))
  expect_equal(suppressWarnings(exp(my.e.spct)[["Rfr"]]),
               rep(exp(0.1), length.out = 10))
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
  expect_lt(sum(as.numeric(reflectance(my.spct, quantity = "contribution",
                                              w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(reflectance(trim_spct(my.spct, range = c(400, 600)),
                                          quantity = "contribution",
                                          w.band = split_bands(c(400, 600), length.out = 3)))), 1)


})

test_that("Rfr_ratio", {
  uvb.wb <- waveband(c(280,315), wb.name = "UVB")
  blue.wb <- waveband(c(400,500), wb.name = "B")

  Rfr.ratio.result <- 0.827905455600525
  expect_equal(
    as.numeric(Rfr_ratio(Ler_leaf_rflt.spct, uvb.wb, blue.wb)),
    Rfr.ratio.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(Rfr_ratio(Ler_leaf_rflt.spct, uvb.wb, blue.wb, scale.factor = 10)),
    Rfr.ratio.result * 10, tolerance = 1e-12)
  expect_equal(
    as.numeric(Rfr_ratio(Ler_leaf_rflt.spct, uvb.wb, blue.wb, quantity = "mean")),
    Rfr.ratio.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(Rfr_ratio(Ler_leaf_rflt.spct, uvb.wb, blue.wb, quantity = "average")),
    Rfr.ratio.result, tolerance = 1e-12)

  Rfr.ratio.result <- 0.289766909460184
  expect_equal(
    as.numeric(Rfr_ratio(Ler_leaf_rflt.spct, uvb.wb, blue.wb, quantity = "total")),
    Rfr.ratio.result, tolerance = 1e-12)

  expect_error(Rfr_ratio(Ler_leaf_rflt.spct, uvb.wb, blue.wb, quantity = "bad argument"))
  expect_warning(Rfr_ratio(sun.spct, uvb.wb, blue.wb))

  expect_named(
    Rfr_ratio(Ler_leaf_rflt.spct, uvb.wb, blue.wb),
    "UVB:B[Rfr(wl):Rfr(wl)]")

  Rfr.fraction.result <- 0.45292575338834
  expect_equal(
    as.numeric(Rfr_fraction(Ler_leaf_rflt.spct, uvb.wb, blue.wb)),
    Rfr.fraction.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(Rfr_fraction(Ler_leaf_rflt.spct, uvb.wb, blue.wb, scale.factor = 10)),
    Rfr.fraction.result * 10, tolerance = 1e-12)
  expect_equal(
    as.numeric(Rfr_fraction(Ler_leaf_rflt.spct, uvb.wb, blue.wb, quantity = "mean")),
    Rfr.fraction.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(Rfr_fraction(Ler_leaf_rflt.spct, uvb.wb, blue.wb, quantity = "average")),
    Rfr.fraction.result, tolerance = 1e-12)

  Rfr.fraction.result <- 0.224666106204773
  expect_equal(
    as.numeric(Rfr_fraction(Ler_leaf_rflt.spct, uvb.wb, blue.wb, quantity = "total")),
    Rfr.fraction.result, tolerance = 1e-12)

  expect_error(Rfr_fraction(Ler_leaf_rflt.spct, uvb.wb, blue.wb, quantity = "bad argument"))
  expect_warning(Rfr_fraction(sun.spct, uvb.wb, blue.wb))

  expect_named(
    Rfr_fraction(Ler_leaf_rflt.spct, uvb.wb, blue.wb),
    "UVB:(UVB+B)[Rfr(wl):Rfr(wl)]")

  Rfr.normdiff.result <- 0.0941484932233198
  expect_equal(
    as.numeric(Rfr_normdiff(Ler_leaf_rflt.spct, blue.wb, uvb.wb)),
    Rfr.normdiff.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(Rfr_normdiff(Ler_leaf_rflt.spct, blue.wb, uvb.wb, scale.factor = 10)),
    Rfr.normdiff.result * 10, tolerance = 1e-12)
  expect_equal(
    as.numeric(Rfr_normdiff(Ler_leaf_rflt.spct, blue.wb, uvb.wb, quantity = "mean")),
    Rfr.normdiff.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(Rfr_normdiff(Ler_leaf_rflt.spct, blue.wb, uvb.wb, quantity = "average")),
    Rfr.normdiff.result, tolerance = 1e-12)

  Rfr.normdiff.result <- 0.550667787590454
  expect_equal(
    as.numeric(Rfr_normdiff(Ler_leaf_rflt.spct, blue.wb, uvb.wb, quantity = "total")),
    Rfr.normdiff.result, tolerance = 1e-12)

  expect_error(Rfr_normdiff(Ler_leaf_rflt.spct, blue.wb, uvb.wb, quantity = "bad argument"))
  expect_warning(Rfr_normdiff(sun.spct, blue.wb, uvb.wb))

  expect_named(
    Rfr_normdiff(Ler_leaf_rflt.spct, blue.wb, uvb.wb),
    "(B-UVB):(B+UVB)[Rfr(wl):Rfr(wl)]")
})

