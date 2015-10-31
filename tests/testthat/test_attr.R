library("photobiology")
library("lubridate")

context("set_get")

test_that("any_spct", {

  my.spct <- object_spct(w.length = 400:450, Tfr = 0.5, Rfr = 0.5)
  tested.time <- ymd_hms("2015-12-31 23:59:59")

  setWhenMeasured(my.spct, tested.time)
  expect_equal(getWhenMeasured(my.spct), tested.time)
  setWhenMeasured(my.spct, NULL)
  expect_true(is.na(getWhenMeasured(my.spct)))
  setWhenMeasured(my.spct, tested.time)
  expect_equal(getWhenMeasured(my.spct), tested.time)

  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)

  setWhereMeasured(my.spct, tested.location)
  expect_equal(getWhereMeasured(my.spct), tested.location)
  setWhereMeasured(my.spct, NULL)
  expect_true(is.data.frame(getWhereMeasured(my.spct)))
  expect_true(all(is.na(getWhereMeasured(my.spct))))

  tested.location <- data.frame(lon = 24.93545, lat = 60.16952,
                                address = "Helsinki")

  setWhereMeasured(my.spct, tested.location)
  expect_equal(getWhereMeasured(my.spct), tested.location)

  tested.location <- data.frame(lon = 1, lat = 2)

  setWhereMeasured(my.spct, lon = 1, lat = 2)
  expect_equal(getWhereMeasured(my.spct), tested.location)

  expect_error(setWhereMeasured(my.spct, 1L))
  expect_error(setWhereMeasured(my.spct, "here"))

  expect_equal(getSpctVersion(my.spct), 2L)
})

context("set_get_mspct")

test_that("any_mspct", {

  my.spct <- filter_spct(w.length = 400:450, Tfr = 0.5)
  tested.time1 <- ymd_hms("2015-12-31 23:59:59")
  tested.time2 <- ymd_hms("2015-12-30 23:59:59")
  my.mspct <- filter_mspct(list(A = my.spct, B = my.spct))
  my.mspct[["A"]] <- setWhenMeasured(my.mspct[["A"]], tested.time1)
  my.mspct[["B"]] <- setWhenMeasured(my.mspct[["B"]], tested.time2)
  expect_true(is.data.frame(getWhenMeasured(my.mspct)))
#  expect_true(is.data.frame(getWhenMeasured(my.mspct["A"])))
#  expect_true(is.data.frame(getWhenMeasured(my.mspct["B"])))
  expect_true(is.POSIXct(getWhenMeasured(my.mspct[["A"]])))
  expect_true(is.POSIXct(getWhenMeasured(my.mspct[["B"]])))
  expect_true(all(is.POSIXct(getWhenMeasured(my.mspct)[["when.measured"]])))
  expect_equal(getWhenMeasured(my.mspct)[["when.measured"]][1], tested.time1)
  expect_equal(getWhenMeasured(my.mspct)[["when.measured"]][2], tested.time2)
  expect_equal(getWhenMeasured(my.mspct[["A"]]), tested.time1)
  expect_equal(getWhenMeasured(my.mspct[["B"]]), tested.time2)

  tested.location1 <- data.frame(lon = 10, lat = 20)
  tested.location2 <- data.frame(lon = 15, lat = 25)
  my.mspct[["A"]] <- setWhereMeasured(my.mspct[["A"]], tested.location1)
  my.mspct[["B"]] <- setWhereMeasured(my.mspct[["B"]], tested.location2)
  expect_true(is.data.frame(getWhereMeasured(my.mspct[["A"]])))
  expect_true(is.data.frame(getWhereMeasured(my.mspct[["B"]])))
  expect_true(all(is.numeric(getWhereMeasured(my.mspct)[["lon"]])))
  expect_true(all(is.numeric(getWhereMeasured(my.mspct)[["lat"]])))
  expect_equal(getWhereMeasured(my.mspct)[["lon"]][1], tested.location1[[1]])
  expect_equal(getWhereMeasured(my.mspct)[["lat"]][1], tested.location1[[2]])
  expect_equal(getWhereMeasured(my.mspct)[["lon"]][2], tested.location2[[1]])
  expect_equal(getWhereMeasured(my.mspct)[["lat"]][2], tested.location2[[2]])
  expect_equal(getWhereMeasured(my.mspct[["A"]]), tested.location1)
  expect_equal(getWhereMeasured(my.mspct[["B"]]), tested.location2)

})

context("conversions")

test_that("object_spct", {

  my.spct <- object_spct(w.length = 400:450, Tfr = 0.5, Rfr = 0.5)
  setTfrType(my.spct, "internal")
  setRfrType(my.spct, "specular")

  expect_equal(getTfrType(my.spct), "internal")
  expect_equal(getTfrType(as.filter_spct(my.spct)), "internal")

  expect_equal(getRfrType(my.spct), "specular")
  expect_equal(getRfrType(as.reflector_spct(my.spct)), "specular")

  expect_equal(getTfrType(as.object_spct(my.spct)), "internal")
  expect_equal(getRfrType(as.object_spct(my.spct)), "specular")

  expect_equal(getTfrType(as.object_spct(as.filter_spct(my.spct))), "internal")
  expect_equal(getRfrType(as.object_spct(as.reflector_spct(my.spct))), "specular")

})

context("fscale_spct")

test_that("fscale", {

  my.spct <- q2e(sun.spct, action = "replace")

  expect_less_than(abs(integrate_spct(fscale(my.spct, f = "total")) - 1), 1e-6)
  expect_less_than(abs(average_spct(fscale(my.spct, f = "mean")) - 1), 1e-6)
  expect_warning(irrad(fscale(my.spct, f = "mean")))
  expect_equal(irrad(fscale(my.spct, f = "mean")), NA)
  expect_named(fscale(my.spct), names(my.spct))
  expect_equal(class(fscale(my.spct)), class(my.spct))
  expect_error(fscale(my.spct, range = 100))
  expect_error(normalize(my.spct, range = c(100, 100)))
  expect_true(is_scaled(fscale(my.spct)))
  expect_false(is_scaled(my.spct))
  expect_equal(is.source_spct(fscale(my.spct)), is.source_spct(my.spct))
  expect_equal(is.filter_spct(fscale(my.spct)), is.filter_spct(my.spct))
  expect_equal(getTimeUnit(fscale(my.spct)), getTimeUnit(my.spct))
  expect_equal(comment(fscale(my.spct)), comment(my.spct))
})

context("integrate_spct average_spct")

test_that("integrate_spct", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)

  expect_equivalent(integrate_spct(my.spct), sum(my.spct$s.e.irrad) - 1)
  expect_equivalent(average_spct(my.spct), 1)
  expect_equivalent(average_spct(my.spct * 2), 2)
  expect_named(average_spct(my.spct), "e.irrad")


  my.spct <- source_spct(w.length=seq(from=1000, to=2000, by=10), s.e.irrad = 1)

  expect_equivalent(integrate_spct(my.spct), (sum(my.spct$s.e.irrad) - 1) * 10)
  expect_equivalent(average_spct(my.spct), 1)
  expect_equivalent(average_spct(my.spct * 2), 2)

  e2q(my.spct, byref = TRUE)

  expect_equivalent(average_spct(my.spct), c(1, 1.2538837047156523583e-05))
  expect_named(average_spct(my.spct), c("e.irrad", "q.irrad"))

  e2q(my.spct, action="replace", byref = TRUE)

  expect_equivalent(average_spct(my.spct), 1.2538837047156523583e-05)
  expect_named(average_spct(my.spct), "q.irrad")

})
