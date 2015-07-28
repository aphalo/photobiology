library("photobiology")
library("lubridate")
context("response_spct")

test_that("constructor energy", {

  my.spct <- response_spct(w.length = 400:409, s.e.response = 1)
  expect_equal(class(my.spct)[1:2], c("response_spct", "generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 2)

  my.s.spct <- response_spct(w.length = 400:409, s.e.response = 1, time.unit = "second")
  my.h.spct <- response_spct(w.length = 400:409, s.e.response = 1, time.unit = "hour")
  my.d.spct <- response_spct(w.length = 400:409, s.e.response = 1, time.unit = "day")
  my.e.spct <- response_spct(w.length = 400:409, s.e.response = 1, time.unit = "exposure")
  my.b.spct <- response_spct(w.length = 400:409, s.e.response = 1, time.unit = "zzz")
  my.ds.spct <- response_spct(w.length = 400:409, s.e.response = 1, time.unit = duration(1, "seconds"))
  my.dh.spct <- response_spct(w.length = 400:409, s.e.response = 1, time.unit = duration(1, "hours"))

  expect_warning(my.b.spct <- response_spct(w.length = 400:409, s.e.response = 1, time.unit = "zzz"))
  expect_equal(my.spct[["s.e.response"]], rep(1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "s.e.response"))
  expect_named(my.s.spct, c("w.length", "s.e.response"))
  expect_named(my.d.spct, c("w.length", "s.e.response"))
  expect_named(my.e.spct, c("w.length", "s.e.response"))
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

  my.spct <- response_spct(w.length = 400:409, s.q.response = 1)
  expect_equal(class(my.spct)[1:2], c("response_spct", "generic_spct") )

  my.s.spct <- response_spct(w.length = 400:409, s.q.response = 1, time.unit = "second")
  my.h.spct <- response_spct(w.length = 400:409, s.q.response = 1, time.unit = "hour")
  my.d.spct <- response_spct(w.length = 400:409, s.q.response = 1, time.unit = "day")
  my.e.spct <- response_spct(w.length = 400:409, s.q.response = 1, time.unit = "exposure")
  my.b.spct <- response_spct(w.length = 400:409, s.q.response = 1, time.unit = "zzz")
  my.ds.spct <- response_spct(w.length = 400:409, s.q.response = 1, time.unit = duration(1, "seconds"))
  my.dh.spct <- response_spct(w.length = 400:409, s.q.response = 1, time.unit = duration(1, "hours"))

  expect_warning(my.b.spct <- response_spct(w.length = 400:409, s.q.response = 1, time.unit = "zzz"))
  expect_equal(my.spct[["s.q.response"]], rep(1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "s.q.response"))
  expect_named(my.s.spct, c("w.length", "s.q.response"))
  expect_named(my.d.spct, c("w.length", "s.q.response"))
  expect_named(my.e.spct, c("w.length", "s.q.response"))
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

  my.e.spct <- response_spct(w.length = 400:409, s.e.response = 1)
  my.2e.spct <- response_spct(w.length = 400:409, s.e.response = 2)

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

  options(photobiology.radiation.unit = NULL)
})

test_that("oper energy energy", {

  my.e.spct <- response_spct(w.length = 400:409, s.e.response = 1)
  my.2e.spct <- response_spct(w.length = 400:409, s.e.response = 2)

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

  my.q.spct <- response_spct(w.length = 400:409, s.q.response = 1)
  my.2q.spct <- response_spct(w.length = 400:409, s.q.response = 2)

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

  my.q.spct <- response_spct(w.length = 400:409, s.q.response = 1)
  my.2q.spct <- response_spct(w.length = 400:409, s.q.response = 2)

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

  options(photobiology.radiation.unit = NULL)
})

test_that("oper energy photon", {

  my.e.spct <- response_spct(w.length = 400:409, s.e.response = 1)
  my.2e.spct <- response_spct(w.length = 400:409, s.e.response = 2)

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
  expect_equal( 1 / (2 / my.2e.spct), +my.e.spct)
  expect_equal( 1 / my.e.spct, my.e.spct^-1)

  options(photobiology.radiation.unit = NULL)
})

test_that("math energy energy", {

  my.e.spct <- response_spct(w.length = 400:409, s.e.response = 1)
  my.2e.spct <- response_spct(w.length = 400:409, s.e.response = 2)

  options(photobiology.radiation.unit = "energy")

  expect_equal(log10(my.e.spct)[["s.e.response"]],  rep(log10(1), length.out = 10))
  expect_equal(log(my.e.spct)[["s.e.response"]],  rep(log(1), length.out = 10))
  expect_equal(log(my.e.spct, 2)[["s.e.response"]],  rep(log(1, 2), length.out = 10))
  expect_equal(exp(my.e.spct)[["s.e.response"]],  rep(exp(1), length.out = 10))
  expect_equal(sqrt(my.e.spct)[["s.e.response"]],  rep(sqrt(1), length.out = 10))

  options(photobiology.radiation.unit = NULL)
})

test_that("math photon photon", {

  my.q.spct <- response_spct(w.length = 400:409, s.q.response = 1)
  my.2q.spct <- response_spct(w.length = 400:409, s.q.response = 2)

  options(photobiology.radiation.unit = "photon")

  expect_equal(log10(my.q.spct)[["s.q.response"]],  rep(log10(1), length.out = 10))
  expect_equal(log(my.q.spct)[["s.q.response"]],  rep(log(1), length.out = 10))
  expect_equal(log(my.q.spct, 2)[["s.q.response"]],  rep(log(1, 2), length.out = 10))
  expect_equal(exp(my.q.spct)[["s.q.response"]],  rep(exp(1), length.out = 10))
  expect_equal(sqrt(my.q.spct)[["s.q.response"]],  rep(sqrt(1), length.out = 10))

  options(photobiology.radiation.unit = NULL)
})

test_that("response e_response q_response", {
  my.spct <- response_spct(w.length = 300:700, s.e.response = 1)

  response.result <- 400
  expect_equal(as.numeric(response(my.spct)), response.result, tolerance = 1e-6)
  expect_equal(as.numeric(response(my.spct, quantity = "total")), response.result, tolerance = 1e-6)
  expect_equal(as.numeric(response(my.spct, quantity = "average")), 1, tolerance = 1e-6)
  expect_equal(as.numeric(response(my.spct, quantity = "mean")), 1, tolerance = 1e-6)
  expect_equal(as.numeric(response(my.spct, time.unit = "second")),
               response.result, tolerance = 1e-6)
  expect_equal(as.numeric(response(my.spct, time.unit = "hour")),
               response.result * 3600, tolerance = 1e-6)
  expect_equal(as.numeric(response(my.spct, time.unit = duration(1))),
               response.result, tolerance = 1e-6)
  expect_equal(as.numeric(response(my.spct, time.unit = duration(0.5))),
               response.result * 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(response(my.spct, time.unit = duration(1, "minutes"))),
               response.result * 60, tolerance = 1e-6)
  expect_equal(sum(as.numeric(response(my.spct, quantity = "total",
                                       w.band = split_bands(my.spct, length.out = 3)))),
               response.result)
  expect_equal(sum(as.numeric(response(my.spct, quantity = "average",
                                       w.band = split_bands(my.spct, length.out = 3)))), 3)
  expect_equal(sum(as.numeric(response(my.spct, quantity = "average",
                                       w.band = split_bands(my.spct, length.out = 5)))), 5)

  expect_equal(sum(as.numeric(response(my.spct, quantity = "relative",
                                      w.band = split_bands(my.spct, length.out = 3)))), 1)
  expect_equal(sum(as.numeric(response(my.spct, quantity = "relative",
                                      w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(response(my.spct, quantity = "contribution",
                                      w.band = split_bands(my.spct, length.out = 3)))), 1)
  expect_less_than(sum(as.numeric(response(my.spct, quantity = "contribution",
                                          w.band = split_bands(c(400, 600), length.out = 3)))), 1)
  expect_equal(sum(as.numeric(response(trim_spct(my.spct, range = c(400, 600)),
                                      quantity = "contribution",
                                      w.band = split_bands(c(400, 600), length.out = 3)))), 1)

  expect_equal(as.numeric(e_response(my.spct)), response.result, tolerance = 1e-6)
  expect_equal(as.numeric(e_response(my.spct, time.unit = "second")),
               response.result, tolerance = 1e-6)
  expect_equal(as.numeric(e_response(my.spct, time.unit = "hour")),
               response.result * 3600, tolerance = 1e-6)
  expect_equal(as.numeric(e_response(my.spct, time.unit = duration(1))),
               response.result, tolerance = 1e-6)
  expect_equal(as.numeric(e_response(my.spct, time.unit = duration(0.5))),
               response.result * 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(e_response(my.spct, time.unit = duration(1, "minutes"))),
               response.result * 60, tolerance = 1e-6)
  response.result <- 101360909
  expect_equal(as.numeric(q_response(my.spct)), response.result, tolerance = 1e-6)
  expect_equal(as.numeric(q_response(my.spct, time.unit = "second")),
               response.result, tolerance = 1e-6)
  expect_equal(as.numeric(q_response(my.spct, time.unit = "hour")),
               response.result * 3600, tolerance = 1e-6)
  expect_equal(as.numeric(q_response(my.spct, time.unit = duration(1))),
               response.result, tolerance = 1e-6)
  expect_equal(as.numeric(q_response(my.spct, time.unit = duration(0.5))),
               response.result * 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(q_response(my.spct, time.unit = duration(1, "minutes"))),
               response.result * 60, tolerance = 1e-6)
})
