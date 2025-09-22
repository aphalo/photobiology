library("lubridate")

context("source_spct")

test_that("constructor energy", {

  empty.spct <- source_spct()
  expect_true(is.source_spct(empty.spct))
  expect_true(is.any_spct(empty.spct))
  expect_named(empty.spct, c("w.length", "s.e.irrad"))
  expect_equal(nrow(empty.spct), 0L)

  my.spct <- source_spct(w.length = 400:409, s.e.irrad = 1)
  expect_equal(class(my.spct)[1:2], c("source_spct", "generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 3)
  expect_named(my.spct, c("w.length", "s.e.irrad"))

  expect_true(is.source_spct(my.spct))
  expect_true(is.any_spct(my.spct))
  expect_false(is.cps_spct(my.spct))
  expect_false(is.response_spct(my.spct))
  expect_false(is.filter_spct(my.spct))
  expect_false(is.reflector_spct(my.spct))
  expect_false(is.object_spct(my.spct))
  expect_false(is.raw_spct(my.spct))
  expect_false(is.chroma_spct(my.spct))

  my.df <- data.frame(w.length = 400:409, s.e.irrad = 1)
  my.spct <- as.source_spct(my.df)

  expect_equal(class(my.spct)[1:2], c("source_spct", "generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 3)
  expect_equal(my.spct[["s.e.irrad"]], rep(1, length.out = 10))
  expect_named(my.spct, c("w.length", "s.e.irrad"))
  expect_true(is.source_spct(my.spct))
  expect_true(is.any_spct(my.spct))

  my.s.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "second")
  my.h.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "hour")
  my.d.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "day")
  my.e.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "exposure")
  expect_warning(my.b.spct <- source_spct(w.length = 400:409, s.e.irrad = 1,
                                          time.unit = "zzz"))
  my.ds.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = lubridate::duration(1, "seconds"))
  my.dh.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = lubridate::duration(1, "hours"))

  expect_warning(my.b.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "zzz"))
  expect_equal(my.spct[["s.e.irrad"]], rep(1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "s.e.irrad"))
  expect_named(my.s.spct, c("w.length", "s.e.irrad"))
  expect_named(my.d.spct, c("w.length", "s.e.irrad"))
  expect_named(my.e.spct, c("w.length", "s.e.irrad"))

  expect_equal(getTimeUnit(my.spct), "second")
  expect_equal(getTimeUnit(my.s.spct), "second")
  expect_equal(getTimeUnit(my.h.spct), "hour")
  expect_equal(getTimeUnit(my.d.spct), "day")
  expect_equal(getTimeUnit(my.e.spct), "exposure")
  expect_equal(getTimeUnit(my.b.spct), "unknown")
  expect_equal(getTimeUnit(my.ds.spct), duration(1, "seconds"))
  expect_equal(getTimeUnit(my.dh.spct), duration(1, "hours"))

  expect_equal(getTimeUnit(my.spct, force.duration = TRUE), lubridate::duration(1, "seconds"))
  expect_equal(getTimeUnit(my.s.spct, force.duration = TRUE), lubridate::duration(1, "seconds"))
  expect_equal(getTimeUnit(my.h.spct, force.duration = TRUE), lubridate::duration(1, "hour"))
  expect_equal(getTimeUnit(my.d.spct, force.duration = TRUE), lubridate::duration(1, "day"))
  expect_equal(getTimeUnit(my.e.spct, force.duration = TRUE), lubridate::duration(NA_real_))
  expect_equal(getTimeUnit(my.b.spct, force.duration = TRUE), lubridate::duration(NA_real_))
  expect_equal(getTimeUnit(my.ds.spct, force.duration = TRUE), lubridate::duration(1, "seconds"))
  expect_equal(getTimeUnit(my.dh.spct, force.duration = TRUE), lubridate::duration(1, "hours"))
})

test_that("constructor photon", {

  my.spct <- source_spct(w.length = 400:409, s.q.irrad = 1)
  expect_equal(class(my.spct)[1:2], c("source_spct", "generic_spct") )

  my.s.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "second")
  my.h.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "hour")
  my.d.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "day")
  my.e.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "exposure")
  expect_warning(my.b.spct <- source_spct(w.length = 400:409, s.q.irrad = 1,
                                          time.unit = "zzz"))
  my.ds.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = duration(1, "seconds"))
  my.dh.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = duration(1, "hours"))

  expect_warning(my.b.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "zzz"))
  expect_equal(my.spct[["s.q.irrad"]], rep(1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "s.q.irrad"))
  expect_named(my.s.spct, c("w.length", "s.q.irrad"))
  expect_named(my.d.spct, c("w.length", "s.q.irrad"))
  expect_named(my.e.spct, c("w.length", "s.q.irrad"))
  expect_equal(getTimeUnit(my.spct), "second")
  expect_equal(getTimeUnit(my.s.spct), "second")
  expect_equal(getTimeUnit(my.h.spct), "hour")
  expect_equal(getTimeUnit(my.d.spct), "day")
  expect_equal(getTimeUnit(my.e.spct), "exposure")
  expect_equal(getTimeUnit(my.b.spct), "unknown")
  expect_equal(getTimeUnit(my.ds.spct), duration(1, "seconds"))
  expect_equal(getTimeUnit(my.dh.spct), duration(1, "hours"))
})

test_that("oper energy energy", {

  my.e.spct <- source_spct(w.length = 400:409, s.e.irrad = 1)
  my.2e.spct <- source_spct(w.length = 400:409, s.e.irrad = 2)

  options(photobiology.radiation.unit = "energy")

  expect_equal(my.e.spct + my.e.spct,  my.2e.spct)
  expect_equal(my.e.spct * 2, my.2e.spct)
  expect_equal(my.e.spct * 2L, my.2e.spct)
  expect_equal(my.2e.spct / 2, my.e.spct)
  expect_equal(my.2e.spct / 2L, my.e.spct)
  expect_silent(-my.e.spct * -2)
  expect_silent(-my.e.spct * -2L)
  expect_silent(-my.2e.spct / -2)
  expect_silent(-my.2e.spct / -2L)
  expect_equal(suppressWarnings(-my.e.spct * -2), my.2e.spct)
  expect_equal(suppressWarnings(-my.e.spct * -2L), my.2e.spct)
  expect_equal(suppressWarnings(-my.2e.spct / -2), my.e.spct)
  expect_equal(suppressWarnings(-my.2e.spct / -2L), my.e.spct)
  expect_equal( 2 * my.e.spct, my.2e.spct)
  expect_equal( 1 / (2 / my.2e.spct), my.e.spct)
  expect_equal( 1 / my.e.spct, my.e.spct^-1)
  expect_equal(my.2e.spct %/% 2L, my.e.spct)
  expect_equal(my.2e.spct %% 2L, my.e.spct %% 1L)

  options(photobiology.radiation.unit = NULL)
})

test_that("oper energy energy", {

  my.e.spct <- source_spct(w.length = 400:409, s.e.irrad = 1)
  my.2e.spct <- source_spct(w.length = 400:409, s.e.irrad = 2)

  options(photobiology.radiation.unit = NULL)

  expect_equal(my.e.spct + my.e.spct,  my.2e.spct)
  expect_equal(my.e.spct * 2, my.2e.spct)
  expect_equal(my.e.spct * 2L, my.2e.spct)
  expect_equal(my.2e.spct / 2, my.e.spct)
  expect_equal(my.2e.spct / 2L, my.e.spct)
  expect_silent(-my.e.spct * -2)
  expect_silent(-my.e.spct * -2L)
  expect_silent(-my.2e.spct / -2)
  expect_silent(-my.2e.spct / -2L)
  expect_equal(suppressWarnings(-my.e.spct * -2), my.2e.spct)
  expect_equal(suppressWarnings(-my.e.spct * -2L), my.2e.spct)
  expect_equal(suppressWarnings(-my.2e.spct / -2), my.e.spct)
  expect_equal(suppressWarnings(-my.2e.spct / -2L), my.e.spct)
  expect_equal( 2 * my.e.spct, my.2e.spct)
  expect_equal( 1 / (2 / my.2e.spct), my.e.spct)
  expect_equal( 1 / my.e.spct, my.e.spct^-1)

})

test_that("oper photon energy", {

  my.q.spct <- source_spct(w.length = 400:409, s.q.irrad = 1)
  my.2q.spct <- source_spct(w.length = 400:409, s.q.irrad = 2)

  options(photobiology.radiation.unit = "energy")

  expect_equal(my.q.spct + my.q.spct,  +my.2q.spct)
  expect_equal(my.q.spct * 2, +my.2q.spct)
  expect_equal(my.q.spct * 2L, +my.2q.spct)
  expect_equal(my.2q.spct / 2, +my.q.spct)
  expect_equal(my.2q.spct / 2L, +my.q.spct)
  expect_silent(-my.q.spct * -2)
  expect_silent(-my.q.spct * -2L)
  expect_silent(-my.2q.spct / -2)
  expect_silent(-my.2q.spct / -2L)
  expect_equal(suppressWarnings(-my.q.spct * -2), +my.2q.spct)
  expect_equal(suppressWarnings(-my.q.spct * -2L), +my.2q.spct)
  expect_equal(suppressWarnings(-my.2q.spct / -2), +my.q.spct)
  expect_equal(suppressWarnings(-my.2q.spct / -2L), +my.q.spct)
  expect_equal( 2 * my.q.spct, +my.2q.spct)
  expect_equal( 1 / (2 / my.2q.spct), +my.q.spct)
  expect_equal( 1 / my.q.spct, my.q.spct^-1)

  options(photobiology.radiation.unit = NULL)
})

test_that("oper photon photon", {

  my.q.spct <- source_spct(w.length = 400:409, s.q.irrad = 1)
  my.2q.spct <- source_spct(w.length = 400:409, s.q.irrad = 2)

  options(photobiology.radiation.unit = "photon")

  expect_equal(my.q.spct + my.q.spct,  my.2q.spct)
  expect_equal(my.q.spct * 2, my.2q.spct)
  expect_equal(my.q.spct * 2L, my.2q.spct)
  expect_equal(my.2q.spct / 2, my.q.spct)
  expect_equal(my.2q.spct / 2L, my.q.spct)
  expect_silent(-my.q.spct * -2)
  expect_silent(-my.q.spct * -2L)
  expect_silent(-my.2q.spct / -2)
  expect_silent(-my.2q.spct / -2L)
  expect_equal(suppressWarnings(-my.q.spct * -2), my.2q.spct)
  expect_equal(suppressWarnings(-my.q.spct * -2L), my.2q.spct)
  expect_equal(suppressWarnings(-my.2q.spct / -2), my.q.spct)
  expect_equal(suppressWarnings(-my.2q.spct / -2L), my.q.spct)
  expect_equal( 2 * my.q.spct, my.2q.spct)
  expect_equal( 1 / (2 / my.2q.spct), my.q.spct)
  expect_equal( 1 / my.q.spct, my.q.spct^-1)
  expect_equal(my.2q.spct %/% 2L, my.q.spct)
  expect_equal(my.2q.spct %% 2L, my.q.spct %% 1L)

  options(photobiology.radiation.unit = NULL)
})

test_that("oper energy photon", {

  my.e.spct <- source_spct(w.length = 400:409, s.e.irrad = 1)
  my.2e.spct <- source_spct(w.length = 400:409, s.e.irrad = 2)

  options(photobiology.radiation.unit = "photon")

  expect_equal(my.e.spct + my.e.spct,  +my.2e.spct)
  expect_equal(my.e.spct * 2, +my.2e.spct)
  expect_equal(my.e.spct * 2L, +my.2e.spct)
  expect_equal(my.2e.spct / 2, +my.e.spct)
  expect_equal(my.2e.spct / 2L, +my.e.spct)
  expect_silent(-my.e.spct * -2)
  expect_silent(-my.e.spct * -2L)
  expect_silent(-my.2e.spct / -2)
  expect_silent(-my.2e.spct / -2L)
  expect_equal(suppressWarnings(-my.e.spct * -2), +my.2e.spct)
  expect_equal(suppressWarnings(-my.e.spct * -2L), +my.2e.spct)
  expect_equal(suppressWarnings(-my.2e.spct / -2), +my.e.spct)
  expect_equal(suppressWarnings(-my.2e.spct / -2L), +my.e.spct)
  expect_equal( 2 * my.e.spct, +my.2e.spct)
  expect_equal(sum((1 / (2 / my.2e.spct) - my.e.spct)[["s.q.irrad"]]), 0)
  expect_equal( 1 / my.e.spct, my.e.spct^-1)

  options(photobiology.radiation.unit = NULL)
})

test_that("math energy energy", {

  my.e.spct <- source_spct(w.length = 400:409, s.e.irrad = 1)
  my.2e.spct <- source_spct(w.length = 400:409, s.e.irrad = 2)

  options(photobiology.radiation.unit = "energy")

  expect_equal(log10(my.e.spct)[["s.e.irrad"]],  rep(log10(1), length.out = 10))
  expect_equal(log(my.e.spct)[["s.e.irrad"]],  rep(log(1), length.out = 10))
  expect_equal(log(my.e.spct, 2)[["s.e.irrad"]],  rep(log(1, 2), length.out = 10))
  expect_equal(exp(my.e.spct)[["s.e.irrad"]],  rep(exp(1), length.out = 10))
  expect_equal(sqrt(my.e.spct)[["s.e.irrad"]],  rep(sqrt(1), length.out = 10))
  expect_equal(abs(my.e.spct)[["s.e.irrad"]],  rep(abs(1), length.out = 10))
  expect_equal(sign(my.e.spct)[["s.e.irrad"]],  rep(sign(1), length.out = 10))
  expect_equal(cos(my.e.spct)[["s.e.irrad"]],  rep(cos(1), length.out = 10))
  expect_equal(sin(my.e.spct)[["s.e.irrad"]],  rep(sin(1), length.out = 10))
  expect_equal(tan(my.e.spct)[["s.e.irrad"]],  rep(tan(1), length.out = 10))
  expect_equal(acos(my.e.spct)[["s.e.irrad"]],  rep(acos(1), length.out = 10))
  expect_equal(asin(my.e.spct)[["s.e.irrad"]],  rep(asin(1), length.out = 10))
  expect_equal(atan(my.e.spct)[["s.e.irrad"]],  rep(atan(1), length.out = 10))

  options(photobiology.radiation.unit = NULL)
})

test_that("math photon photon", {

  my.q.spct <- source_spct(w.length = 400:409, s.q.irrad = 1)
  my.2q.spct <- source_spct(w.length = 400:409, s.q.irrad = 2)

  options(photobiology.radiation.unit = "photon")

  expect_equal(log10(my.q.spct)[["s.q.irrad"]],  rep(log10(1), length.out = 10))
  expect_equal(log(my.q.spct)[["s.q.irrad"]],  rep(log(1), length.out = 10))
  expect_equal(log(my.q.spct, 2)[["s.q.irrad"]],  rep(log(1, 2), length.out = 10))
  expect_equal(exp(my.q.spct)[["s.q.irrad"]],  rep(exp(1), length.out = 10))
  expect_equal(sqrt(my.q.spct)[["s.q.irrad"]],  rep(sqrt(1), length.out = 10))
  expect_equal(abs(my.q.spct)[["s.q.irrad"]],  rep(abs(1), length.out = 10))
  expect_equal(sign(my.q.spct)[["s.q.irrad"]],  rep(sign(1), length.out = 10))
  expect_equal(cos(my.q.spct)[["s.q.irrad"]],  rep(cos(1), length.out = 10))
  expect_equal(sin(my.q.spct)[["s.q.irrad"]],  rep(sin(1), length.out = 10))
  expect_equal(tan(my.q.spct)[["s.q.irrad"]],  rep(tan(1), length.out = 10))
  expect_equal(acos(my.q.spct)[["s.q.irrad"]],  rep(acos(1), length.out = 10))
  expect_equal(asin(my.q.spct)[["s.q.irrad"]],  rep(asin(1), length.out = 10))
  expect_equal(atan(my.q.spct)[["s.q.irrad"]],  rep(atan(1), length.out = 10))

  options(photobiology.radiation.unit = NULL)
})

test_that("irrad e_irrad q_irrad", {
  irrad.result <- 269.12490033493787678
  expect_equal(as.numeric(irrad(sun.spct)), irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, quantity = "total")), irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, quantity = "average")),
               irrad.result / expanse(sun.spct), tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, quantity = "mean")),
               irrad.result / expanse(sun.spct), tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = "second")),
               irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = "hour")),
               irrad.result * 3600, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = duration(1))),
               irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = duration(0.5))),
               irrad.result * 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = duration(1, "minutes"))),
               irrad.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = minutes(1))),
               irrad.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = hms("00:01:00"))),
               irrad.result * 60, tolerance = 1e-6)
  expect_equal(sum(as.numeric(irrad(sun.spct, quantity = "relative",
                                    w.band = split_bands(sun.spct, length.out = 3)))), 1)
  expect_equal(sum(as.numeric(irrad(sun.spct, quantity = "relative",
                                    w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(irrad(sun.spct, quantity = "contribution",
                                    w.band = split_bands(sun.spct, length.out = 3)))), 1)
  expect_lt(sum(as.numeric(irrad(sun.spct, quantity = "contribution",
                                        w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(irrad(trim_spct(sun.spct, range = c(400, 600)),
                                    quantity = "contribution",
                                    w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_error(irrad(sun.spct, quantity = "bad input",
                     w.band = split_bands(sun.spct, length.out = 3)))

  expect_equal(as.numeric(e_irrad(sun.spct)), irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(e_irrad(sun.spct, time.unit = "second")),
               irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(e_irrad(sun.spct, time.unit = "hour")),
               irrad.result * 3600, tolerance = 1e-6)
  expect_equal(as.numeric(e_irrad(sun.spct, time.unit = duration(1))),
               irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(e_irrad(sun.spct, time.unit = duration(0.5))),
               irrad.result * 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(e_irrad(sun.spct, time.unit = duration(1, "minutes"))),
               irrad.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(e_irrad(sun.spct, time.unit = minutes(1))),
               irrad.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(e_irrad(sun.spct, time.unit = hms("00:01:00"))),
               irrad.result * 60, tolerance = 1e-6)
  expect_equal(sum(as.numeric(e_irrad(sun.spct, quantity = "relative",
                                      w.band = split_bands(sun.spct, length.out = 3)))), 1)
  expect_equal(sum(as.numeric(e_irrad(sun.spct, quantity = "relative",
                                      w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(e_irrad(sun.spct, quantity = "contribution",
                                      w.band = split_bands(sun.spct, length.out = 3)))), 1)
  expect_lt(sum(as.numeric(e_irrad(sun.spct, quantity = "contribution",
                                          w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(e_irrad(trim_spct(sun.spct, range = c(400, 600)),
                                      quantity = "contribution",
                                      w.band = split_bands(c(400, 600), length.out = 3)))), 1)


  irrad.result <- 0.001255353974944903089

  expect_equal(as.numeric(q_irrad(sun.spct)), irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(q_irrad(sun.spct, time.unit = "second")),
               irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(q_irrad(sun.spct, time.unit = "hour")),
               irrad.result * 3600, tolerance = 1e-6)
  expect_equal(as.numeric(q_irrad(sun.spct, time.unit = duration(1))),
               irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(q_irrad(sun.spct, time.unit = duration(0.5))),
               irrad.result * 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(q_irrad(sun.spct, time.unit = duration(1, "minutes"))),
               irrad.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(q_irrad(sun.spct, time.unit = minutes(1))),
               irrad.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(q_irrad(sun.spct, time.unit = hms("00:01:00"))),
               irrad.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(q_irrad(sun.spct, quantity = "relative")), 1)
  expect_equal(as.numeric(q_irrad(sun.spct, quantity = "relative.pc")), 100)
  expect_equal(as.numeric(q_irrad(sun.spct, quantity = "contribution")), 1)
  expect_equal(as.numeric(q_irrad(sun.spct, quantity = "contribution.pc")), 100)
  expect_equal(sum(as.numeric(q_irrad(sun.spct, quantity = "relative",
                                      w.band = split_bands(sun.spct, length.out = 3)))), 1)
  expect_equal(sum(as.numeric(q_irrad(sun.spct, quantity = "relative",
                                      w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(q_irrad(sun.spct, quantity = "contribution",
                                      w.band = split_bands(sun.spct, length.out = 3)))), 1)
  expect_lt(sum(as.numeric(q_irrad(sun.spct, quantity = "contribution",
                                          w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(q_irrad(trim_spct(sun.spct, range = c(400, 600)),
                                      quantity = "contribution",
                                      w.band = split_bands(c(400, 600), length.out = 3)))), 1)

  # cached multipliers
  irrad.result <- 269.1249

  expect_equal(as.numeric(irrad(sun.spct, use.cached.mult = TRUE)),
               irrad.result, tolerance = 1e-6)

  use_cached_mult_as_default(TRUE)
  expect_true(getOption("photobiology.use.cached.mult"))

  expect_equal(as.numeric(irrad(sun.spct)), irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, quantity = "total")), irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, quantity = "average")),
               irrad.result / expanse(sun.spct), tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, quantity = "mean")),
               irrad.result / expanse(sun.spct), tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = "second")),
               irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = "hour")),
               irrad.result * 3600, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = duration(1))),
               irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = duration(0.5))),
               irrad.result * 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = duration(1, "minutes"))),
               irrad.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = minutes(1))),
               irrad.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, time.unit = hms("00:01:00"))),
               irrad.result * 60, tolerance = 1e-6)
  expect_equal(sum(as.numeric(irrad(sun.spct, quantity = "relative",
                                    w.band = split_bands(sun.spct, length.out = 3)))), 1)
  expect_equal(sum(as.numeric(irrad(sun.spct, quantity = "relative",
                                    w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(irrad(sun.spct, quantity = "contribution",
                                    w.band = split_bands(sun.spct, length.out = 3)))), 1)
  expect_lt(sum(as.numeric(irrad(sun.spct, quantity = "contribution",
                                 w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(irrad(trim_spct(sun.spct, range = c(400, 600)),
                                    quantity = "contribution",
                                    w.band = split_bands(c(400, 600), length.out = 3)))), 1)

  use_cached_mult_as_default(FALSE)
  expect_true(!getOption("photobiology.use.cached.mult"))
})

test_that("fluence e_fluence q_fluence", {
  fluence.result <- 269.1249
  expect_error(fluence(sun.spct))
  expect_error(fluence(sun.spct, exposure.time = "second"))
  expect_error(fluence(sun.spct, exposure.time = "hour"))
  expect_equal(as.numeric(fluence(sun.spct, exposure.time = duration(1))),
               fluence.result, tolerance = 1e-6)
  expect_equal(as.numeric(fluence(sun.spct, exposure.time = duration(0.5))),
               fluence.result * 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(fluence(sun.spct, exposure.time = duration(1, "minutes"))),
               fluence.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(fluence(sun.spct, exposure.time = minutes(1))),
               fluence.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(fluence(sun.spct, exposure.time = hms("00:01:00"))),
               fluence.result * 60, tolerance = 1e-6)
  expect_error(e_fluence(sun.spct))
  expect_error(e_fluence(sun.spct, exposure.time = "second"))
  expect_error(e_fluence(sun.spct, exposure.time = "hour"))
  expect_equal(as.numeric(e_fluence(sun.spct, exposure.time = duration(1))),
               fluence.result, tolerance = 1e-6)
  expect_equal(as.numeric(e_fluence(sun.spct, exposure.time = duration(0.5))),
               fluence.result * 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(e_fluence(sun.spct, exposure.time = duration(1, "minutes"))),
               fluence.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(e_fluence(sun.spct, exposure.time = minutes(1))),
               fluence.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(e_fluence(sun.spct, exposure.time = hms("00:01:00"))),
               fluence.result * 60, tolerance = 1e-6)
  fluence.result <- 0.001255353974944903089
  expect_error(q_fluence(sun.spct))
  expect_error(q_fluence(sun.spct, exposure.time = "second"))
  expect_error(q_fluence(sun.spct, exposure.time = "hour"))
  expect_equal(as.numeric(q_fluence(sun.spct, exposure.time = duration(1))),
               fluence.result, tolerance = 1e-6)
  expect_equal(as.numeric(q_fluence(sun.spct, exposure.time = duration(0.5))),
               fluence.result * 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(q_fluence(sun.spct, exposure.time = duration(1, "minutes"))),
               fluence.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(q_fluence(sun.spct, exposure.time = minutes(1))),
               fluence.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(q_fluence(sun.spct, exposure.time = hms("00:01:00"))),
               fluence.result * 60, tolerance = 1e-6)

  # cached multipliers
  fluence.result <- 269.12490033493787678

  expect_equal(as.numeric(fluence(sun.spct,
                                  exposure.time = duration(1),
                                  use.cached.mult = TRUE)),
               fluence.result, tolerance = 1e-6)

  use_cached_mult_as_default(TRUE)
  expect_true(getOption("photobiology.use.cached.mult"))

  expect_error(fluence(sun.spct))
  expect_error(fluence(sun.spct, exposure.time = "second"))
  expect_error(fluence(sun.spct, exposure.time = "hour"))
  expect_equal(as.numeric(fluence(sun.spct, exposure.time = duration(1))),
               fluence.result, tolerance = 1e-6)
  expect_equal(as.numeric(fluence(sun.spct, exposure.time = duration(0.5))),
               fluence.result * 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(fluence(sun.spct, exposure.time = duration(1, "minutes"))),
               fluence.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(fluence(sun.spct, exposure.time = minutes(1))),
               fluence.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(fluence(sun.spct, exposure.time = hms("00:01:00"))),
               fluence.result * 60, tolerance = 1e-6)
  expect_error(e_fluence(sun.spct))
  expect_error(e_fluence(sun.spct, exposure.time = "second"))
  expect_error(e_fluence(sun.spct, exposure.time = "hour"))
  expect_equal(as.numeric(e_fluence(sun.spct, exposure.time = duration(1))),
               fluence.result, tolerance = 1e-6)
  expect_equal(as.numeric(e_fluence(sun.spct, exposure.time = duration(0.5))),
               fluence.result * 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(e_fluence(sun.spct, exposure.time = duration(1, "minutes"))),
               fluence.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(e_fluence(sun.spct, exposure.time = minutes(1))),
               fluence.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(e_fluence(sun.spct, exposure.time = hms("00:01:00"))),
               fluence.result * 60, tolerance = 1e-6)
  fluence.result <- 0.001255353974944903089
  expect_error(q_fluence(sun.spct))
  expect_error(q_fluence(sun.spct, exposure.time = "second"))
  expect_error(q_fluence(sun.spct, exposure.time = "hour"))
  expect_equal(as.numeric(q_fluence(sun.spct, exposure.time = duration(1))),
               fluence.result, tolerance = 1e-6)
  expect_equal(as.numeric(q_fluence(sun.spct, exposure.time = duration(0.5))),
               fluence.result * 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(q_fluence(sun.spct, exposure.time = duration(1, "minutes"))),
               fluence.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(q_fluence(sun.spct, exposure.time = minutes(1))),
               fluence.result * 60, tolerance = 1e-6)
  expect_equal(as.numeric(q_fluence(sun.spct, exposure.time = hms("00:01:00"))),
               fluence.result * 60, tolerance = 1e-6)

  use_cached_mult_as_default(FALSE)
  expect_true(!getOption("photobiology.use.cached.mult"))
})

test_that("ratio q_ratio e_ratio", {
  red.wb <- waveband(c(600,700), wb.name = "R")
  blue.wb <- waveband(c(400,500), wb.name = "B")
  narrow.wb <- waveband(c(460,480), wb.name = "B.narrow")

  q.ratio.result <- 0.83169703201068901
  expect_equal(
    as.numeric(q_ratio(sun.spct, blue.wb, red.wb)),
    q.ratio.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(q_ratio(sun.spct, blue.wb, red.wb, quantity = "mean")),
    q.ratio.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(q_ratio(sun.spct, blue.wb, red.wb, quantity = "average")),
    q.ratio.result, tolerance = 1e-12)
  expect_error(q_ratio(sun.spct, blue.wb, red.wb, quantity = "bad argument"))

  e.ratio.result <- 1.1922198464802185
  expect_equal(
    as.numeric(e_ratio(sun.spct, blue.wb, red.wb)),
    e.ratio.result, tolerance = 1e-112)
  expect_equal(
    as.numeric(e_ratio(sun.spct, blue.wb, red.wb, quantity = "mean")),
    e.ratio.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(e_ratio(sun.spct, blue.wb, red.wb, quantity = "average")),
    e.ratio.result, tolerance = 1e-12)
  expect_error(e_ratio(sun.spct, blue.wb, red.wb, quantity = "bad argument"))

  eq.ratio.result <- 264628.05810932186432
  expect_equal(
    as.numeric(eq_ratio(sun.spct, blue.wb)),
    eq.ratio.result, tolerance = 1e-12)
  qe.ratio.result <- 1 / eq.ratio.result
  expect_equal(
    as.numeric(qe_ratio(sun.spct, blue.wb)),
    qe.ratio.result, tolerance = 1e-12)

  q.ratio.result <- 0.94342725228671586724
  expect_equal(
    as.numeric(q_ratio(sun.spct, narrow.wb, red.wb, quantity = "mean")),
    q.ratio.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(q_ratio(sun.spct, narrow.wb, red.wb, quantity = "average")),
    q.ratio.result, tolerance = 1e-12)

  q.ratio.result <- 0.1886854504573431679
  expect_equal(
    as.numeric(q_ratio(sun.spct, narrow.wb, red.wb, quantity = "total")),
               q.ratio.result, tolerance = 1e-12)

  expect_named(
    q_ratio(sun.spct, blue.wb, red.wb, quantity = "total"),
    "B:R[q:q]")

  expect_named(
    q_ratio(sun.spct, blue.wb, red.wb, quantity = "mean"),
    "B:R[q(wl):q(wl)]")

  e.ratio.result <- 1.3007094880894196631
  expect_equal(
    as.numeric(e_ratio(sun.spct, narrow.wb, red.wb, quantity = "mean")),
    e.ratio.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(e_ratio(sun.spct, narrow.wb, red.wb, quantity = "average")),
    e.ratio.result, tolerance = 1e-12)

  e.ratio.result <- 0.26014189761788392152
  expect_equal(
    as.numeric(e_ratio(sun.spct, narrow.wb, red.wb, quantity = "total")),
    e.ratio.result, tolerance = 1e-12)

  expect_named(
    e_ratio(sun.spct, blue.wb, red.wb, quantity = "total"),
    "B:R[e:e]")

  expect_named(
    e_ratio(sun.spct, blue.wb, red.wb, quantity = "mean"),
    "B:R[e(wl):e(wl)]")
})

test_that("fraction q_fraction e_fraction", {
  red.wb <- waveband(c(600,700), wb.name = "R")
  blue.wb <- waveband(c(400,500), wb.name = "B")
  narrow.wb <- waveband(c(460,480), wb.name = "B.narrow")

  q.fraction.result <- 0.454058186193444
  expect_equal(
    as.numeric(q_fraction(sun.spct, blue.wb, red.wb)),
    q.fraction.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(q_fraction(sun.spct, blue.wb, red.wb, quantity = "mean")),
    q.fraction.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(q_fraction(sun.spct, blue.wb, red.wb, quantity = "average")),
    q.fraction.result, tolerance = 1e-12)
  expect_error(q_fraction(sun.spct, blue.wb, red.wb, quantity = "bad argument"))

  e.fraction.result <- 0.543841370834418
  expect_equal(
    as.numeric(e_fraction(sun.spct, blue.wb, red.wb)),
    e.fraction.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(e_fraction(sun.spct, blue.wb, red.wb, quantity = "mean")),
    e.fraction.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(e_fraction(sun.spct, blue.wb, red.wb, quantity = "average")),
    e.fraction.result, tolerance = 1e-12)
  expect_error(q_fraction(sun.spct, blue.wb, red.wb, quantity = "bad argument"))

  q.fraction.result <- 0.485445107953818
  expect_equal(
    as.numeric(q_fraction(sun.spct, narrow.wb, red.wb, quantity = "mean")),
    q.fraction.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(q_fraction(sun.spct, narrow.wb, red.wb, quantity = "average")),
    q.fraction.result, tolerance = 1e-12)

  q.fraction.result <- 0.158734550325948
  expect_equal(
    as.numeric(q_fraction(sun.spct, narrow.wb, red.wb, quantity = "total")),
    q.fraction.result, tolerance = 1e-12)

  expect_named(
    q_fraction(sun.spct, blue.wb, red.wb, quantity = "total"),
    "B:(B+R)[q:q]")

  expect_named(
    q_fraction(sun.spct, blue.wb, red.wb, quantity = "mean"),
    "B:(B+R)[q(wl):q(wl)]")

  e.fraction.result <- 0.565351468676547
  expect_equal(
    as.numeric(e_fraction(sun.spct, narrow.wb, red.wb, quantity = "mean")),
    e.fraction.result, tolerance = 1e-12)
  expect_equal(
    as.numeric(e_fraction(sun.spct, narrow.wb, red.wb, quantity = "average")),
    e.fraction.result, tolerance = 1e-12)

  e.fraction.result <- 0.206438574980837
  expect_equal(
    as.numeric(e_fraction(sun.spct, narrow.wb, red.wb, quantity = "total")),
    e.fraction.result, tolerance = 1e-12)

  expect_named(
    e_fraction(sun.spct, blue.wb, red.wb, quantity = "total"),
    "B:(B+R)[e:e]")

  expect_named(
    e_fraction(sun.spct, blue.wb, red.wb, quantity = "mean"),
    "B:(B+R)[e(wl):e(wl)]")
})
