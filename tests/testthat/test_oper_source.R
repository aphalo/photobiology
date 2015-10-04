library("photobiology")
library("lubridate")

context("source_spct")

test_that("constructor energy", {

  my.spct <- source_spct(w.length = 400:409, s.e.irrad = 1)
  expect_equal(class(my.spct)[1:2], c("source_spct", "generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 2)

  my.s.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "second")
  my.h.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "hour")
  my.d.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "day")
  my.e.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "exposure")
  my.b.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "zzz")
  my.ds.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = duration(1, "seconds"))
  my.dh.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = duration(1, "hours"))

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
})

test_that("constructor photon", {

  my.spct <- source_spct(w.length = 400:409, s.q.irrad = 1)
  expect_equal(class(my.spct)[1:2], c("source_spct", "generic_spct") )

  my.s.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "second")
  my.h.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "hour")
  my.d.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "day")
  my.e.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "exposure")
  my.b.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "zzz")
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
  expect_equal(-my.e.spct * -2, my.2e.spct)
  expect_equal(-my.e.spct * -2L, my.2e.spct)
  expect_equal(-my.2e.spct / -2, my.e.spct)
  expect_equal(-my.2e.spct / -2L, my.e.spct)
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
  expect_equal(-my.e.spct * -2, my.2e.spct)
  expect_equal(-my.e.spct * -2L, my.2e.spct)
  expect_equal(-my.2e.spct / -2, my.e.spct)
  expect_equal(-my.2e.spct / -2L, my.e.spct)
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
  expect_equal(-my.q.spct * -2, +my.2q.spct)
  expect_equal(-my.q.spct * -2L, +my.2q.spct)
  expect_equal(-my.2q.spct / -2, +my.q.spct)
  expect_equal(-my.2q.spct / -2L, +my.q.spct)
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
  expect_equal(-my.q.spct * -2, my.2q.spct)
  expect_equal(-my.q.spct * -2L, my.2q.spct)
  expect_equal(-my.2q.spct / -2, my.q.spct)
  expect_equal(-my.2q.spct / -2L, my.q.spct)
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
  expect_equal(-my.e.spct * -2, +my.2e.spct)
  expect_equal(-my.e.spct * -2L, +my.2e.spct)
  expect_equal(-my.2e.spct / -2, +my.e.spct)
  expect_equal(-my.2e.spct / -2L, +my.e.spct)
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

  options(photobiology.radiation.unit = NULL)
})

test_that("irrad e_irrad q_irrad", {
  irrad.result <- 269.1249
  expect_equal(as.numeric(irrad(sun.spct)), irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, quantity = "total")), irrad.result, tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, quantity = "average")),
               irrad.result / spread(sun.spct), tolerance = 1e-6)
  expect_equal(as.numeric(irrad(sun.spct, quantity = "mean")),
               irrad.result / spread(sun.spct), tolerance = 1e-6)
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
  expect_less_than(sum(as.numeric(irrad(sun.spct, quantity = "contribution",
                                        w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(irrad(trim_spct(sun.spct, range = c(400, 600)),
                                    quantity = "contribution",
                                    w.band = split_bands(c(400, 600), length.out = 3)))), 1)

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
  expect_less_than(sum(as.numeric(e_irrad(sun.spct, quantity = "contribution",
                                          w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(e_irrad(trim_spct(sun.spct, range = c(400, 600)),
                                      quantity = "contribution",
                                      w.band = split_bands(c(400, 600), length.out = 3)))), 1)


  irrad.result <- 0.001255336

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
  expect_less_than(sum(as.numeric(q_irrad(sun.spct, quantity = "contribution",
                                          w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(q_irrad(trim_spct(sun.spct, range = c(400, 600)),
                                      quantity = "contribution",
                                      w.band = split_bands(c(400, 600), length.out = 3)))), 1)
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
  fluence.result <- 0.001255336
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
})

