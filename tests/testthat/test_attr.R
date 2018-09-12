library("photobiology")
library("lubridate")

context("set_get")

test_that("any_spct", {

  my.spct <- object_spct(w.length = 400:450, Tfr = 0.5, Rfr = 0.5)
  tested.time <- ymd_hms("2015-12-31 23:59:59", tz = "UTC")

  setWhenMeasured(my.spct, tested.time)
  expect_equal(getWhenMeasured(my.spct), tested.time)
  expect_is(getWhenMeasured(my.spct), "POSIXct")

  setWhenMeasured(my.spct, NULL)
  expect_true(is.na(getWhenMeasured(my.spct)))
  expect_is(getWhenMeasured(my.spct), "POSIXct")

  setWhenMeasured(my.spct, tested.time)
  expect_equal(getWhenMeasured(my.spct), tested.time)
  expect_is(getWhenMeasured(my.spct), "POSIXct")

  tested.date <- ymd("2015-12-30", tz = "UTC")
  target <- lubridate::as_datetime(tested.date, tz = "UTC")
  setWhenMeasured(my.spct, tested.date)
  expect_equal(getWhenMeasured(my.spct), target)
  expect_is(getWhenMeasured(my.spct), "POSIXct")

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

  tested.locationz <- data.frame(lat = 2, lon = 1)
  setWhereMeasured(my.spct, tested.location)
  getWhenMeasured(my.spct)

  tested.what <- "user message"

  setWhatMeasured(my.spct, tested.what)
  expect_equal(getWhatMeasured(my.spct), tested.what)

  setMultipleWl(my.spct, 2)
  expect_equal(getMultipleWl(my.spct), 2)

  setMultipleWl(my.spct, 1)
  expect_equal(getMultipleWl(my.spct), 1)

  my.descriptor <- list(spectrometer.name = "fake",
                        spectrometer.sn = "12345AB c",
                        bench.grating = "fake 123",
                        bench.slit = "10um")
  class(my.descriptor) <- c("instr_desc", class(my.descriptor))
  expect_true(!isValidInstrDesc(my.spct))
  expect_equal(length(getInstrDesc(my.spct)), 4)
  expect_is(getInstrDesc(my.spct), "instr_desc")
  setInstrDesc(my.spct, my.descriptor)
  expect_true(isValidInstrDesc(my.spct))
  expect_equal(getInstrDesc(my.spct), my.descriptor)
  expect_is(getInstrDesc(my.spct), "instr_desc")
  expect_equal(length(getInstrDesc(my.spct)), 4)
  expect_identical(trimInstrDesc(my.spct, c("*")), my.spct)
  expect_equal(length(getInstrDesc(my.spct)), 4)
  expect_identical(trimInstrDesc(my.spct, names(my.descriptor)), my.spct)
  expect_equal(length(getInstrDesc(my.spct)), 4)
  trimInstrDesc(my.spct, c("-", "bench.slit"))
  expect_true(isValidInstrDesc(my.spct))
  expect_is(getInstrDesc(my.spct), "instr_desc")
  expect_equal(length(getInstrDesc(my.spct)), 3)
  expect_equal(names(getInstrDesc(my.spct)), setdiff(names(my.descriptor), "bench.slit"))

  my.settings <- list(integ.time = 321,
                      tot.time = 1000,
                      num.scans = 50,
                      rel.signal = 0.86)
  class(my.settings) <- c("instr_settings", class(my.settings))
  expect_equal(length(getInstrSettings(my.spct)), 4)
  expect_is(getInstrSettings(my.spct), "instr_settings")
  setInstrSettings(my.spct, my.settings)
  expect_equal(getInstrSettings(my.spct), my.settings)
  expect_equal(length(getInstrSettings(my.spct)), 4)
  expect_identical(trimInstrSettings(my.spct, c("*")), my.spct)
  expect_is(getInstrSettings(my.spct), "instr_settings")
  expect_identical(trimInstrSettings(my.spct, names(my.settings)), my.spct)
  expect_is(getInstrSettings(my.spct), "instr_settings")
  trimInstrSettings(my.spct, c("-", "num.scans"))
  expect_true(isValidInstrSettings(my.spct))
  expect_is(getInstrSettings(my.spct), "instr_settings")
  expect_equal(length(getInstrSettings(my.spct)), 3)
  expect_equal(names(getInstrSettings(my.spct)), setdiff(names(my.settings), "num.scans"))

  expect_equal(getSpctVersion(my.spct), 2L)
})

context("set_get_mspct")

test_that("any_mspct", {

  my.spct <- filter_spct(w.length = 400:450, Tfr = 0.5)
  tested.time1 <- ymd_hms("2015-12-31 23:59:59 UTC")
  tested.time2 <- ymd_hms("2015-12-30 23:59:59 UTC")
  my.mspct <- filter_mspct(list(A = my.spct, B = my.spct))
  my.mspct[["A"]] <- setWhenMeasured(my.mspct[["A"]], tested.time1)
  my.mspct[["B"]] <- setWhenMeasured(my.mspct[["B"]], tested.time2)
  expect_is(getWhenMeasured(my.mspct), "data.frame")
  expect_is(getWhenMeasured(my.mspct["A"]), "data.frame")
  expect_is(getWhenMeasured(my.mspct["B"]), "data.frame")
  expect_is(getWhenMeasured(my.mspct[c("A", "B")]), "data.frame")
  expect_is(getWhenMeasured(my.mspct[c("B", "A")]), "data.frame")
  expect_is(getWhenMeasured(my.mspct[["A"]]), "POSIXct")
  expect_is(getWhenMeasured(my.mspct[["B"]]), "POSIXct")
  expect_is(getWhenMeasured(my.mspct)[["when.measured"]], "POSIXct")
  expect_equal(getWhenMeasured(my.mspct)[["when.measured"]][1], tested.time1)
  expect_equal(getWhenMeasured(my.mspct)[["when.measured"]][2], tested.time2)
  expect_equal(getWhenMeasured(my.mspct[["A"]]), tested.time1)
  expect_equal(getWhenMeasured(my.mspct[["B"]]), tested.time2)

  expect_error(setWhenMeasured(my.mspct, "A"))
  expect_error(setWhenMeasured(my.mspct, 100))
  expect_error(setWhenMeasured(my.mspct, c(tested.time1, tested.time2)))

  setWhenMeasured(my.mspct, list(tested.time1, tested.time2))
  expect_true(is.POSIXct(getWhenMeasured(my.mspct[["A"]])))
  expect_true(is.POSIXct(getWhenMeasured(my.mspct[["B"]])))
  expect_true(all(is.POSIXct(getWhenMeasured(my.mspct)[["when.measured"]])))
  expect_equal(getWhenMeasured(my.mspct)[["when.measured"]][1], tested.time1)
  expect_equal(getWhenMeasured(my.mspct)[["when.measured"]][2], tested.time2)
  expect_equal(getWhenMeasured(my.mspct[["A"]]), tested.time1)
  expect_equal(getWhenMeasured(my.mspct[["B"]]), tested.time2)

  setWhenMeasured(my.mspct, tested.time1)
  expect_true(is.POSIXct(getWhenMeasured(my.mspct[["A"]])))
  expect_true(is.POSIXct(getWhenMeasured(my.mspct[["B"]])))
  expect_true(all(is.POSIXct(getWhenMeasured(my.mspct)[["when.measured"]])))
  expect_equal(getWhenMeasured(my.mspct)[["when.measured"]][1], tested.time1)
  expect_equal(getWhenMeasured(my.mspct)[["when.measured"]][2], tested.time1)
  expect_equal(getWhenMeasured(my.mspct[["A"]]), tested.time1)
  expect_equal(getWhenMeasured(my.mspct[["B"]]), tested.time1)

  expect_equal(ncol(getWhenMeasured(my.mspct, idx = F)), 1L)
  expect_equal(ncol(getWhenMeasured(my.mspct, idx = T)), 2L)
  expect_equal(ncol(getWhenMeasured(my.mspct, idx = NULL)), 2L)
  expect_equal(ncol(getWhenMeasured(my.mspct, idx = "abc")), 2L)

  tested.location1 <- data.frame(lon = 10, lat = 20)
  tested.location2 <- data.frame(lon = 15, lat = 25)
  my.mspct[["A"]] <- setWhereMeasured(my.mspct[["A"]], tested.location1)
  my.mspct[["B"]] <- setWhereMeasured(my.mspct[["B"]], tested.location2)
  expect_true(is.data.frame(getWhereMeasured(my.mspct[["A"]])))
  expect_true(is.data.frame(getWhereMeasured(my.mspct[["B"]])))
  expect_true(all(is.numeric(getWhereMeasured(my.mspct)[["lon"]])))
  expect_true(all(is.numeric(getWhereMeasured(my.mspct)[["lat"]])))
  expect_equal(getWhereMeasured(my.mspct)[["lon"]][1], tested.location1[["lon"]])
  expect_equal(getWhereMeasured(my.mspct)[["lat"]][1], tested.location1[["lat"]])
  expect_equal(getWhereMeasured(my.mspct)[["lon"]][2], tested.location2[["lon"]])
  expect_equal(getWhereMeasured(my.mspct)[["lat"]][2], tested.location2[["lat"]])
  expect_equal(getWhereMeasured(my.mspct[["A"]]), tested.location1)
  expect_equal(getWhereMeasured(my.mspct[["B"]]), tested.location2)

  expect_error(setWhereMeasured(my.mspct, "A"))
  expect_error(setWhereMeasured(my.mspct, 100))
  expect_error(setWhereMeasured(my.mspct, c(tested.location1, tested.location2)))

  setWhereMeasured(my.mspct, NULL)
  expect_true(is.data.frame(getWhereMeasured(my.mspct[["A"]])))
  expect_true(is.data.frame(getWhereMeasured(my.mspct[["B"]])))
  expect_true(all(is.numeric(getWhereMeasured(my.mspct)[["lon"]])))
  expect_true(all(is.numeric(getWhereMeasured(my.mspct)[["lat"]])))

  setWhereMeasured(my.mspct, tested.location1)
  expect_true(is.data.frame(getWhereMeasured(my.mspct[["A"]])))
  expect_true(is.data.frame(getWhereMeasured(my.mspct[["B"]])))
  expect_true(all(is.numeric(getWhereMeasured(my.mspct)[["lon"]])))
  expect_true(all(is.numeric(getWhereMeasured(my.mspct)[["lat"]])))
  expect_equal(getWhereMeasured(my.mspct[["A"]]), tested.location1)
  expect_equal(getWhereMeasured(my.mspct[["B"]]), tested.location1)

  tested.locations <- rbind(tested.location1, tested.location2)
  setWhereMeasured(my.mspct, tested.locations)
  expect_true(is.data.frame(getWhereMeasured(my.mspct[["A"]])))
  expect_true(is.data.frame(getWhereMeasured(my.mspct[["B"]])))
  expect_true(all(is.numeric(getWhereMeasured(my.mspct)[["lon"]])))
  expect_true(all(is.numeric(getWhereMeasured(my.mspct)[["lat"]])))
  expect_equal(getWhereMeasured(my.mspct)[["lon"]][1], tested.location1[["lon"]])
  expect_equal(getWhereMeasured(my.mspct)[["lat"]][1], tested.location1[["lat"]])
  expect_equal(getWhereMeasured(my.mspct)[["lon"]][2], tested.location2[["lon"]])
  expect_equal(getWhereMeasured(my.mspct)[["lat"]][2], tested.location2[["lat"]])
  expect_equal(getWhereMeasured(my.mspct[["A"]]), tested.locations[1, ]) # row names match
  expect_equal(getWhereMeasured(my.mspct[["B"]]), tested.locations[2, ]) # row names match

  tested.locations <- list(tested.location1, tested.location2)
  setWhereMeasured(my.mspct, tested.locations)
  expect_true(is.data.frame(getWhereMeasured(my.mspct[["A"]])))
  expect_true(is.data.frame(getWhereMeasured(my.mspct[["B"]])))
  expect_true(all(is.numeric(getWhereMeasured(my.mspct)[["lon"]])))
  expect_true(all(is.numeric(getWhereMeasured(my.mspct)[["lat"]])))
  expect_equal(getWhereMeasured(my.mspct)[["lon"]][1], tested.location1[["lon"]])
  expect_equal(getWhereMeasured(my.mspct)[["lat"]][1], tested.location1[["lat"]])
  expect_equal(getWhereMeasured(my.mspct)[["lon"]][2], tested.location2[["lon"]])
  expect_equal(getWhereMeasured(my.mspct)[["lat"]][2], tested.location2[["lat"]])
  expect_equal(getWhereMeasured(my.mspct[["A"]]), tested.locations[[1]])
  expect_equal(getWhereMeasured(my.mspct[["B"]]), tested.locations[[2]])

  tested.location1z <- data.frame(lat = 20, lon = 10)
  tested.location2z <- data.frame(lat = 25, lon = 15)
  tested.locationsz <- rbind(tested.location1z, tested.location2z)
  setWhereMeasured(my.mspct, tested.locationsz)
  expect_equal(getWhereMeasured(my.mspct)[["lon"]][1], tested.location1[["lon"]])
  expect_equal(getWhereMeasured(my.mspct)[["lat"]][1], tested.location1[["lat"]])
  expect_equal(getWhereMeasured(my.mspct)[["lon"]][2], tested.location2[["lon"]])
  expect_equal(getWhereMeasured(my.mspct)[["lat"]][2], tested.location2[["lat"]])
  expect_equal(getWhereMeasured(my.mspct[["A"]])[["lon"]], tested.location1[["lon"]])
  expect_equal(getWhereMeasured(my.mspct[["A"]])[["lat"]], tested.location1[["lat"]])
  expect_equal(getWhereMeasured(my.mspct)[["lon"]][2], tested.location2[["lon"]])
  expect_equal(getWhereMeasured(my.mspct)[["lat"]][2], tested.location2[["lat"]])

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

  expect_lt(abs(integrate_spct(fscale(my.spct, f = "total")) - 1), 1e-6)
  expect_lt(abs(average_spct(fscale(my.spct, f = "mean")) - 1), 1e-6)
  expect_warning(irrad(fscale(my.spct, f = "mean")))
  expect_equivalent(suppressWarnings(irrad(fscale(my.spct, f = "mean"))), 520)
  expect_named(fscale(my.spct), names(my.spct))
  expect_equal(class(fscale(my.spct)), class(my.spct))
  expect_warning(fscale(my.spct, range = 100))
  expect_error(fscale(my.spct, range = 281))
  expect_warning(fscale(my.spct, range = c(100, 100)))
  expect_error(fscale(my.spct, range = c(281, 281)))
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

context("return same attributes")


test_that("trim_wl attr", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct, range = 110:200)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct, range = 110:200,
                                                fill = 0)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct, range = 100:190)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct, range = 100:190,
                                                fill = 0)))),
               character(0) )

    expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct, range = 100:210)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct, range = 100:210,
                                                fill = 0)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct)))),
               character(0) )

})

test_that("clip_wl attr", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clip_wl(my.spct, range = 110:200)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clip_wl(my.spct, range = 100:190)))),
               character(0) )


  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clip_wl(my.spct, range = 100:210)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clip_wl(my.spct)))),
               character(0) )

})

test_that("clean attr", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clean(my.spct, range = 110:200)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clean(my.spct, range = 100:190)))),
               character(0) )


  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clean(my.spct, range = 100:210)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clean(my.spct)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clean(my.spct, range.s.data = c(0,0.9))))),
               character(0) )
})

test_that("fshift attr", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 100:105)))),
               character(0) )


  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 195:200)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 145:150)))),
               character(0) )

})

test_that("fshift attr2", {

  my.spct <- raw_spct(w.length=100:200, counts = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 100:105)))),
               character(0) )


  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 195:200)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 145:150)))),
               character(0) )

  my.spct <- raw_spct(w.length=100:200, counts1 = 1, counts2 = 2)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 100:105)))),
               character(0) )


  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 195:200)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 145:150)))),
               character(0) )
})

test_that("fscale attr", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(fscale(my.spct))),
               names(attributes(my.spct))),
               "scaled" )

  expect_equal(setdiff(names(attributes(fscale(my.spct, range = 100:110))),
                       names(attributes(my.spct))),
               "scaled" )

  expect_equal(setdiff(names(attributes(fscale(my.spct, f = "total"))),
                       names(attributes(my.spct))),
               "scaled" )

})

test_that("normalize attr", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(normalize(my.spct))),
                       names(attributes(my.spct))),
               "normalized" )

  expect_equal(setdiff(names(attributes(normalize(my.spct, range = 100:110))),
                       names(attributes(my.spct))),
               "normalized" )

  expect_equal(setdiff(names(attributes(normalize(my.spct, norm = "max"))),
                       names(attributes(my.spct))),
               "normalized" )

  expect_equal(setdiff(names(attributes(normalize(my.spct, norm = 150))),
                       names(attributes(my.spct))),
               "normalized" )
})

test_that("peaks attr", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(peaks(my.spct))),
                       names(attributes(my.spct))),
               character(0)  )

  expect_equal(setdiff(names(attributes(valleys(my.spct))),
                       names(attributes(my.spct))),
               character(0)  )

})

test_that("smooth_spct attr", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(setdiff(names(attributes(smooth_spct(my.spct))),
                       names(attributes(my.spct))), "comment"),
               character(0) )

})

test_that("extract attr", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(my.spct[2:50, ]))),
               character(0) )

})

test_that("merge attr in operators", {
  expect_equal((white_led.source_spct + white_led.source_spct) / 2, white_led.source_spct)
  expect_equal(white_led.source_spct + white_led.source_spct, white_led.source_spct * 2)
  expect_equal(white_led.source_spct + 0, white_led.source_spct)
  expect_equal(white_led.source_spct - white_led.source_spct, white_led.source_spct * 0)
})
