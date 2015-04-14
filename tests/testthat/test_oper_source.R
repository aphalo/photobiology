library("photobiology")
context("source_spct")

test_that("constructor energy", {

  my.spct <- source_spct(w.length = 400:409, s.e.irrad = 1)
  my.s.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "second")
  my.h.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "hour")
  my.d.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "day")
  my.b.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "zzz")

  expect_warning(my.b.spct <- source_spct(w.length = 400:409, s.e.irrad = 1, time.unit = "zzz"))
  expect_equal(my.spct[["s.e.irrad"]], rep(1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "s.e.irrad"))
  expect_named(my.s.spct, c("w.length", "s.e.irrad"))
  expect_named(my.d.spct, c("w.length", "s.e.irrad"))
  expect_equal(getTimeUnit(my.spct), "second")
  expect_equal(getTimeUnit(my.s.spct), "second")
  expect_equal(getTimeUnit(my.h.spct), "hour")
  expect_equal(getTimeUnit(my.d.spct), "day")
  expect_equal(getTimeUnit(my.b.spct), "unknown")
})

test_that("constructor photon", {

  my.spct <- source_spct(w.length = 400:409, s.q.irrad = 1)
  my.s.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "second")
  my.h.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "hour")
  my.d.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "day")
  my.b.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "zzz")

  expect_warning(my.b.spct <- source_spct(w.length = 400:409, s.q.irrad = 1, time.unit = "zzz"))
  expect_equal(my.spct[["s.q.irrad"]], rep(1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "s.q.irrad"))
  expect_named(my.s.spct, c("w.length", "s.q.irrad"))
  expect_named(my.d.spct, c("w.length", "s.q.irrad"))
  expect_equal(getTimeUnit(my.spct), "second")
  expect_equal(getTimeUnit(my.s.spct), "second")
  expect_equal(getTimeUnit(my.h.spct), "hour")
  expect_equal(getTimeUnit(my.d.spct), "day")
  expect_equal(getTimeUnit(my.b.spct), "unknown")
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
  expect_equal( 1 / (2 / my.2e.spct), +my.e.spct)
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

  options(photobiology.radiation.unit = NULL)
})
