library("lubridate")

context("set_get_attr_spct")

test_that("single-spct: attrs are good", {

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

  tested.location <-
    validate_geocode(data.frame(lon = 24.93545, lat = 60.16952))

  setWhereMeasured(my.spct, tested.location)
  expect_equal(getWhereMeasured(my.spct), tested.location)
  setWhereMeasured(my.spct, NULL)
  expect_true(is.data.frame(getWhereMeasured(my.spct)))
  expect_true(all(is.na(getWhereMeasured(my.spct))))

  tested.location <-
    validate_geocode(data.frame(lon = 24.93545, lat = 60.16952,
                                address = "Helsinki",
                                stringsAsFactors = FALSE))

  setWhereMeasured(my.spct, tested.location)
  expect_equal(getWhereMeasured(my.spct), tested.location)

  tested.location <- validate_geocode(data.frame(lon = 1, lat = 2))
  expected.location <-
    validate_geocode(data.frame(lon = 1, lat = 2,
                                address = NA_character_,
                                stringsAsFactors = FALSE))

  setWhereMeasured(my.spct, lon = 1, lat = 2)
  expect_equal(getWhereMeasured(my.spct), expected.location)

  expect_error(setWhereMeasured(my.spct, 1L))
  expect_error(setWhereMeasured(my.spct, "here"))

  tested.locationz <- validate_geocode(data.frame(lat = 2, lon = 1))
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
  expect_equal(length(getInstrDesc(my.spct)), 4)
  expect_true(is.na(getInstrDesc(my.spct)$bench.slit))

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

  expect_equal(getSpctVersion(my.spct), 3L)
})

context("get_attr_long_spct")

test_that("long-spct: attrs are good", {

  my.spct <- sun_evening.spct

  expect_is(getWhenMeasured(my.spct), "list")
  expect_named(getWhenMeasured(my.spct), levels(my.spct$spct.idx))

  expect_is(getWhereMeasured(my.spct), "tbl_df")
  expect_named(getWhereMeasured(my.spct), c("lat", "lon", "address"))

  expect_is(getMultipleWl(my.spct), "integer")

  expect_true(isValidInstrDesc(my.spct))
  expect_equal(length(getInstrDesc(my.spct)), 6)
  expect_is(getInstrDesc(my.spct), "instr_desc")

  expect_equal(length(getInstrSettings(my.spct)), 16)
  expect_is(getInstrSettings(my.spct), "instr_settings")

  expect_equal(getSpctVersion(my.spct), 3L)

  expect_is(getIdFactor(my.spct), "character")
  expect_equal(length(getIdFactor(my.spct)), 1L)
  expect_equal(getIdFactor(my.spct), "spct.idx")

  setIdFactor(my.spct, "time")
  expect_is(getIdFactor(my.spct), "character")
  expect_equal(length(getIdFactor(my.spct)), 1L)
  expect_equal(getIdFactor(my.spct), "time")
  expect_true("time" %in% colnames(my.spct))

  setIdFactor(my.spct, NULL)
  expect_is(getIdFactor(my.spct), "character")
  expect_equal(length(getIdFactor(my.spct)), 1L)
  expect_true(is.na(getIdFactor(my.spct)))
  expect_true("time" %in% colnames(my.spct))

  setIdFactor(my.spct, "time")
  expect_is(getIdFactor(my.spct), "character")
  expect_equal(length(getIdFactor(my.spct)), 1L)
  expect_equal(getIdFactor(my.spct), "time")
  expect_true("time" %in% colnames(my.spct))

  setIdFactor(my.spct, NULL)
  expect_error(setIdFactor(my.spct, "inexistent.variable"))
})


context("set_get_attr_mspct")

test_that("mspct: attrs are good", {

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

  tested.location1 <- validate_geocode(data.frame(lon = 10, lat = 20))
  tested.location2 <- validate_geocode(data.frame(lon = 15, lat = 25))
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

  tested.location1z <- validate_geocode(data.frame(lat = 20, lon = 10))
  tested.location2z <- validate_geocode(data.frame(lat = 25, lon = 15))
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

context("RfrType and TfrType")

test_that("object_spct: Tfr and Rfr are good", {

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
  expect_message(as.object_spct(as.reflector_spct(my.spct)))
  expect_named(as.object_spct(as.reflector_spct(my.spct)),
               c("w.length", "Rfr", "Tfr"), ignore.order = TRUE)

  setTfrType(my.spct, "total")
  setRfrType(my.spct, "total")
  expect_equal(getTfrType(my.spct), "total")
  expect_equal(getRfrType(my.spct), "total")
  expect_silent(as.object_spct(as.reflector_spct(my.spct)))
  expect_named(as.object_spct(as.reflector_spct(my.spct)),
               c("w.length", "Rfr", "Tfr"), ignore.order = TRUE)
  expect_equal(getTfrType(as.object_spct(as.filter_spct(my.spct))), "total")
  expect_equal(getRfrType(as.object_spct(as.reflector_spct(my.spct))), "total")

})

