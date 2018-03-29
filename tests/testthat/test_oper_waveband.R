library(photobiology)
library(lubridate)

context("waveband")

test_that("construction", {

  wb1.fun <- function(x){x * 2}
  wb1 <- waveband(c(200, 350), SWF.q.fun = wb1.fun,
                  wb.name = "wb1", norm = 300, weight = "BSWF", SWF.norm = 300)

  expect_true(is.waveband(wb1))
  expect_true(is_effective(wb1))
  expect_equal(range(wb1), c(200, 350))
  expect_equal(midpoint(wb1), midpoint(c(200, 350)))
  expect_equal(expanse(wb1), expanse(c(200, 350)))
  expect_equal(min(wb1), 200)
  expect_equal(max(wb1), 350)

  wb2.fun <- function(x){x * 2}
  wb2 <- waveband(c(200, 350), SWF.e.fun = wb2.fun,
                  wb.name = "wb2", norm = 300, weight = "BSWF", SWF.norm = 300)

  expect_true(is.waveband(wb2))
  expect_true(is_effective(wb2))
  expect_equal(range(wb2), c(200, 350))
  expect_equal(midpoint(wb2), midpoint(c(200, 350)))
  expect_equal(expanse(wb2), expanse(c(200, 350)))
  expect_equal(min(wb2), 200)
  expect_equal(max(wb2), 350)

  wb3 <- waveband(300:310)

  expect_true(is.waveband(wb3))
  expect_false(is_effective(wb3))
  expect_equal(range(wb3), c(300, 310))
  expect_equal(midpoint(wb3), midpoint(c(300, 310)))
  expect_equal(expanse(wb3), expanse(c(300, 310)))
  expect_equal(min(wb3), 300)
  expect_equal(max(wb3), 310)

  wb3 <- waveband()

  expect_true(is.waveband(wb3))
  expect_false(is_effective(wb3))
  expect_true(all(is.na(range(wb3))))
  expect_true(all(is.na(midpoint(wb3))))
  expect_true(all(is.na(expanse(wb3))))
  expect_true(all(is.na(min(wb3))))
  expect_true(all(is.na(max(wb3))))

})

test_that("product energy", {

  wb1 <- waveband(c(400,700), wb.name = "wb1")
  wb2.fun <- function(x){x * 2}
  wb2 <- waveband(c(400,700), SWF.q.fun = wb2.fun, SWF.e.fun = wb2.fun,
                  wb.name = "wb2", norm = 300, weight = "BSWF", SWF.norm = 300)

  expect_true(is.waveband(wb1))
  expect_true(is.waveband(wb2))

  expect_false(is_effective(wb1))
  expect_true(is_effective(wb2))

  my.spct <- source_spct(w.length = 200:420, s.e.irrad = 1)

  my.s.spct <- source_spct(w.length = 200:420, s.e.irrad = 1, time.unit = "second")
  my.h.spct <- source_spct(w.length = 200:420, s.e.irrad = 1, time.unit = "hour")
  my.d.spct <- source_spct(w.length = 200:420, s.e.irrad = 1, time.unit = "day")
  my.e.spct <- source_spct(w.length = 200:420, s.e.irrad = 1, time.unit = "exposure")
  expect_warning(my.b.spct <- source_spct(w.length = 200:420, s.e.irrad = 1,
                                          time.unit = "zzz"))
  my.ds.spct <- source_spct(w.length = 200:420, s.e.irrad = 1, time.unit = duration(1, "seconds"))
  my.dh.spct <- source_spct(w.length = 200:420, s.e.irrad = 1, time.unit = duration(1, "hours"))

  expect_equal(getTimeUnit(my.spct * wb1), "second")
  expect_equal(getTimeUnit(my.spct * wb2), "second")
  expect_equal(getTimeUnit(my.h.spct * wb1), "hour")
  expect_equal(getTimeUnit(my.d.spct * wb1), "day")
  expect_equal(getTimeUnit(my.e.spct * wb1), "exposure")
  expect_equal(getTimeUnit(my.b.spct * wb1), "unknown")
  expect_equal(getTimeUnit(my.ds.spct * wb1), duration(1, "seconds"))
  expect_equal(getTimeUnit(my.dh.spct * wb1), duration(1, "hours"))
  expect_true(!is_effective(my.spct * wb1))
  expect_equal(getBSWFUsed(my.spct * wb1), "none")
  expect_true(is_effective(my.spct * wb2))
  expect_equal(getBSWFUsed(my.spct * wb2), as.character(labels(wb2)[2]))
})

test_that("product photon", {

  wb1 <- waveband(c(400,700), wb.name = "wb1")
  wb2.fun <- function(x){x * 2}
  wb2 <- waveband(c(400,700), SWF.q.fun = wb2.fun, SWF.e.fun = wb2.fun,
                  wb.name = "wb2", norm = 300, weight = "BSWF", SWF.norm = 300)

  expect_true(is.waveband(wb1))
  expect_true(is.waveband(wb2))

  expect_false(is_effective(wb1))
  expect_true(is_effective(wb2))

  my.spct <- source_spct(w.length = 200:420, s.q.irrad = 1)

  my.s.spct <- source_spct(w.length = 200:420, s.q.irrad = 1, time.unit = "second")
  my.h.spct <- source_spct(w.length = 200:420, s.q.irrad = 1, time.unit = "hour")
  my.d.spct <- source_spct(w.length = 200:420, s.q.irrad = 1, time.unit = "day")
  my.e.spct <- source_spct(w.length = 200:420, s.q.irrad = 1, time.unit = "exposure")
  expect_warning(my.b.spct <- source_spct(w.length = 200:420, s.q.irrad = 1,
                                          time.unit = "zzz"))
  my.ds.spct <- source_spct(w.length = 200:420, s.q.irrad = 1, time.unit = duration(1, "seconds"))
  my.dh.spct <- source_spct(w.length = 200:420, s.q.irrad = 1, time.unit = duration(1, "hours"))

  expect_equal(getTimeUnit(my.spct * wb2), "second")
  expect_equal(getTimeUnit(my.s.spct * wb1), "second")
  expect_equal(getTimeUnit(my.h.spct * wb1), "hour")
  expect_equal(getTimeUnit(my.d.spct * wb1), "day")
  expect_equal(getTimeUnit(my.e.spct * wb1), "exposure")
  expect_equal(getTimeUnit(my.b.spct * wb1), "unknown")
  expect_equal(getTimeUnit(my.ds.spct * wb1), duration(1, "seconds"))
  expect_equal(getTimeUnit(my.dh.spct * wb1), duration(1, "hours"))
  expect_true(!is_effective(my.spct * wb1))
  expect_equal(getBSWFUsed(my.spct * wb1), "none")
  expect_true(is_effective(my.spct * wb2))
  expect_equal(getBSWFUsed(my.spct * wb2), as.character(labels(wb2)[2]))
})

