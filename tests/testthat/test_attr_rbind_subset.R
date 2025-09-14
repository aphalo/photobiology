library("photobiology")
library("lubridate")

context("wide_long_wide_plain")

test_that("raw_spct", {

  my.spct <- white_led.raw_spct
  tested.location <- validate_geocode(data.frame(lon = 24.93545, lat = 60.16952))
  setWhereMeasured(my.spct, tested.location)
  tested.time <- getWhenMeasured(white_led.cps_spct)
  tested.what <- getWhatMeasured(my.spct)
  tested.instr.desc <- getInstrDesc(my.spct)
  tested.instr.settings <- getInstrSettings(my.spct)
  normalization.flag <- is_normalized(my.spct)
  scaling.flag <- is_scaled(my.spct)

  spct.l <- list(one = my.spct, two = my.spct)

  expect_named(spct.l, c("one", "two"))
  expect_equal(getWhenMeasured(spct.l[["one"]]), tested.time)
  expect_equal(getWhenMeasured(spct.l[["two"]]), tested.time)
  expect_equal(getWhereMeasured(spct.l[["one"]]), tested.location)
  expect_equal(getWhereMeasured(spct.l[["two"]]), tested.location)
  expect_equal(getWhatMeasured(spct.l[["one"]]), tested.what)
  expect_equal(getWhatMeasured(spct.l[["two"]]), tested.what)
  expect_equal(getInstrDesc(spct.l[["one"]]), tested.instr.desc)
  expect_equal(getInstrDesc(spct.l[["two"]]), tested.instr.desc)
  expect_equal(getInstrSettings(spct.l[["one"]]), tested.instr.settings)
  expect_equal(getInstrSettings(spct.l[["two"]]), tested.instr.settings)
  expect_equal(is_normalized(spct.l[["one"]]), normalization.flag)
  expect_equal(is_normalized(spct.l[["two"]]), normalization.flag)
  expect_equal(is_scaled(spct.l[["one"]]), scaling.flag)
  expect_equal(is_scaled(spct.l[["two"]]), scaling.flag)

  ## collection
  my.mspct <- raw_mspct(spct.l)

  expect_named(my.mspct, c("one", "two"))
  expect_equal(getWhenMeasured(my.mspct[["one"]]), tested.time)
  expect_equal(getWhenMeasured(my.mspct[["two"]]), tested.time)
  expect_equal(getWhereMeasured(my.mspct[["one"]]), tested.location)
  expect_equal(getWhereMeasured(my.mspct[["two"]]), tested.location)
  expect_equal(getWhatMeasured(my.mspct[["one"]]), tested.what)
  expect_equal(getWhatMeasured(my.mspct[["two"]]), tested.what)
  expect_equal(getInstrDesc(my.mspct[["one"]]), tested.instr.desc)
  expect_equal(getInstrDesc(my.mspct[["two"]]), tested.instr.desc)
  expect_equal(getInstrSettings(my.mspct[["one"]]), tested.instr.settings)
  expect_equal(getInstrSettings(my.mspct[["two"]]), tested.instr.settings)
  expect_equal(is_normalized(my.mspct[["one"]]), normalization.flag)
  expect_equal(is_normalized(my.mspct[["two"]]), normalization.flag)
  expect_equal(is_scaled(my.mspct[["one"]]), scaling.flag)
  expect_equal(is_scaled(my.mspct[["two"]]), scaling.flag)

  expect_is(getInstrDesc(my.mspct[["one"]]), "instr_desc")
  expect_is(getInstrDesc(my.mspct[["two"]]), "instr_desc")
  expect_is(getInstrSettings(my.mspct[["one"]]), "instr_settings")
  expect_is(getInstrSettings(my.mspct[["two"]]), "instr_settings")

  ## long form
  long.spct <- rbindspct(my.mspct)

  expect_named(getWhenMeasured(long.spct), c("one", "two"))
  expect_named(getWhatMeasured(long.spct), c("one", "two"))
  expect_equal(getWhereMeasured(long.spct)$spct.idx, c("one", "two"))
  expect_named(getInstrDesc(long.spct), c("one", "two"))
  expect_named(getInstrSettings(long.spct), c("one", "two"))

  expect_equal(getWhenMeasured(long.spct)[["one"]], tested.time)
  expect_equal(getWhenMeasured(long.spct)[["two"]], tested.time)
  expect_equal(getWhereMeasured(long.spct)[1, -1], tested.location)
  expect_equal(getWhereMeasured(long.spct)[2, -1], tested.location)
  expect_equal(getWhatMeasured(long.spct)[["one"]], tested.what)
  expect_equal(getWhatMeasured(long.spct)[["two"]], tested.what)
  expect_equal(getInstrDesc(long.spct)[["one"]], tested.instr.desc)
  expect_equal(getInstrDesc(long.spct)[["two"]], tested.instr.desc)
  expect_equal(getInstrSettings(long.spct)[["one"]], tested.instr.settings)
  expect_equal(getInstrSettings(long.spct)[["two"]], tested.instr.settings)
  # expect_equal(is_normalized(long.spct)[["one"]], normalization.flag)
  # expect_equal(is_normalized(long.spct)[["two"]], normalization.flag)
  # expect_equal(is_scaled(long.spct)[["one"]], scaling.flag)
  # expect_equal(is_scaled(long.spct)[["two"]], scaling.flag)

  expect_true(!inherits(getInstrDesc(long.spct), "instr_desc"))
  expect_is(getInstrDesc(long.spct)[["one"]], "instr_desc")
  expect_is(getInstrDesc(long.spct)[["two"]], "instr_desc")
  expect_true(!inherits(getInstrSettings(long.spct), "instr_settings"))
  expect_is(getInstrSettings(long.spct)[["one"]], "instr_settings")
  expect_is(getInstrSettings(long.spct)[["two"]], "instr_settings")

  ## recovered collection
  recovered.mspct <- subset2mspct(long.spct)

  expect_named(recovered.mspct, c("one", "two"))

  expect_equal(getWhenMeasured(recovered.mspct[["one"]]), tested.time)
  expect_equal(getWhenMeasured(recovered.mspct[["two"]]), tested.time)
  expect_equal(getWhereMeasured(recovered.mspct[["one"]]), tested.location)
  expect_equal(getWhereMeasured(recovered.mspct[["two"]]), tested.location)
  expect_equal(getWhatMeasured(recovered.mspct[["one"]]), tested.what)
  expect_equal(getWhatMeasured(recovered.mspct[["two"]]), tested.what)
  expect_equal(getInstrDesc(recovered.mspct[["one"]]), tested.instr.desc)
  expect_equal(getInstrDesc(recovered.mspct[["two"]]), tested.instr.desc)
  expect_equal(getInstrSettings(recovered.mspct[["one"]]), tested.instr.settings)
  expect_equal(getInstrSettings(recovered.mspct[["two"]]), tested.instr.settings)
  expect_equal(is_normalized(recovered.mspct[["one"]]), normalization.flag)
  expect_equal(is_normalized(recovered.mspct[["two"]]), normalization.flag)
  expect_equal(is_scaled(recovered.mspct[["one"]]), scaling.flag)
  expect_equal(is_scaled(recovered.mspct[["two"]]), scaling.flag)

  expect_is(getInstrDesc(recovered.mspct[["one"]]), "instr_desc")
  expect_is(getInstrDesc(recovered.mspct[["two"]]), "instr_desc")
  expect_is(getInstrSettings(recovered.mspct[["one"]]), "instr_settings")
  expect_is(getInstrSettings(recovered.mspct[["two"]]), "instr_settings")
})

test_that("cps_spct", {

  my.spct <- white_led.cps_spct
  tested.location <- validate_geocode(data.frame(lon = 24.93545, lat = 60.16952))
  setWhereMeasured(my.spct, tested.location)
  tested.time <- getWhenMeasured(white_led.cps_spct)
  tested.what <- getWhatMeasured(my.spct)
  tested.instr.desc <- getInstrDesc(my.spct)
  tested.instr.settings <- getInstrSettings(my.spct)
  normalization.flag <- is_normalized(my.spct)
  scaling.flag <- is_scaled(my.spct)

  spct.l <- list(one = my.spct, two = my.spct)

  expect_named(spct.l, c("one", "two"))
  expect_equal(getWhenMeasured(spct.l[["one"]]), tested.time)
  expect_equal(getWhenMeasured(spct.l[["two"]]), tested.time)
  expect_equal(getWhereMeasured(spct.l[["one"]]), tested.location)
  expect_equal(getWhereMeasured(spct.l[["two"]]), tested.location)
  expect_equal(getWhatMeasured(spct.l[["one"]]), tested.what)
  expect_equal(getWhatMeasured(spct.l[["two"]]), tested.what)
  expect_equal(getInstrDesc(spct.l[["one"]]), tested.instr.desc)
  expect_equal(getInstrDesc(spct.l[["two"]]), tested.instr.desc)
  expect_equal(getInstrSettings(spct.l[["one"]]), tested.instr.settings)
  expect_equal(getInstrSettings(spct.l[["two"]]), tested.instr.settings)
  expect_equal(is_normalized(spct.l[["one"]]), normalization.flag)
  expect_equal(is_normalized(spct.l[["two"]]), normalization.flag)
  expect_equal(is_scaled(spct.l[["one"]]), scaling.flag)
  expect_equal(is_scaled(spct.l[["two"]]), scaling.flag)

  ## collection
  my.mspct <- cps_mspct(spct.l)

  expect_named(my.mspct, c("one", "two"))
  expect_equal(getWhenMeasured(my.mspct[["one"]]), tested.time)
  expect_equal(getWhenMeasured(my.mspct[["two"]]), tested.time)
  expect_equal(getWhereMeasured(my.mspct[["one"]]), tested.location)
  expect_equal(getWhereMeasured(my.mspct[["two"]]), tested.location)
  expect_equal(getWhatMeasured(my.mspct[["one"]]), tested.what)
  expect_equal(getWhatMeasured(my.mspct[["two"]]), tested.what)
  expect_equal(getInstrDesc(my.mspct[["one"]]), tested.instr.desc)
  expect_equal(getInstrDesc(my.mspct[["two"]]), tested.instr.desc)
  expect_equal(getInstrSettings(my.mspct[["one"]]), tested.instr.settings)
  expect_equal(getInstrSettings(my.mspct[["two"]]), tested.instr.settings)
  expect_equal(is_normalized(my.mspct[["one"]]), normalization.flag)
  expect_equal(is_normalized(my.mspct[["two"]]), normalization.flag)
  expect_equal(is_scaled(my.mspct[["one"]]), scaling.flag)
  expect_equal(is_scaled(my.mspct[["two"]]), scaling.flag)

  expect_is(getInstrDesc(my.mspct[["one"]]), "instr_desc")
  expect_is(getInstrDesc(my.mspct[["two"]]), "instr_desc")
  expect_is(getInstrSettings(my.mspct[["one"]]), "instr_settings")
  expect_is(getInstrSettings(my.mspct[["two"]]), "instr_settings")

  ## long form
  long.spct <- rbindspct(my.mspct)

  expect_named(getWhenMeasured(long.spct), c("one", "two"))
  expect_named(getWhatMeasured(long.spct), c("one", "two"))
  expect_equal(getWhereMeasured(long.spct)[["spct.idx"]], c("one", "two"))
  expect_named(getInstrDesc(long.spct), c("one", "two"))
  expect_named(getInstrSettings(long.spct), c("one", "two"))

  expect_equal(getWhenMeasured(long.spct)[["one"]], tested.time)
  expect_equal(getWhenMeasured(long.spct)[["two"]], tested.time)
  expect_equal(getWhereMeasured(long.spct)[1, -1], tested.location)
  expect_equal(getWhereMeasured(long.spct)[2, -1], tested.location)
  expect_equal(getWhatMeasured(long.spct)[["one"]], tested.what)
  expect_equal(getWhatMeasured(long.spct)[["two"]], tested.what)
  expect_equal(getInstrDesc(long.spct)[["one"]], tested.instr.desc)
  expect_equal(getInstrDesc(long.spct)[["two"]], tested.instr.desc)
  expect_equal(getInstrSettings(long.spct)[["one"]], tested.instr.settings)
  expect_equal(getInstrSettings(long.spct)[["two"]], tested.instr.settings)
  # expect_equal(is_normalized(long.spct)[["one"]], normalization.flag)
  # expect_equal(is_normalized(long.spct)[["two"]], normalization.flag)
  # expect_equal(is_scaled(long.spct)[["one"]], scaling.flag)
  # expect_equal(is_scaled(long.spct)[["two"]], scaling.flag)

  expect_true(!inherits(getInstrDesc(long.spct), "instr_desc"))
  expect_is(getInstrDesc(long.spct)[["one"]], "instr_desc")
  expect_is(getInstrDesc(long.spct)[["two"]], "instr_desc")
  expect_true(!inherits(getInstrSettings(long.spct), "instr_settings"))
  expect_is(getInstrSettings(long.spct)[["one"]], "instr_settings")
  expect_is(getInstrSettings(long.spct)[["two"]], "instr_settings")

  ## recovered collection
  recovered.mspct <- subset2mspct(long.spct)

  expect_named(recovered.mspct, c("one", "two"))

  expect_equal(getWhenMeasured(recovered.mspct[["one"]]), tested.time)
  expect_equal(getWhenMeasured(recovered.mspct[["two"]]), tested.time)
  expect_equal(getWhereMeasured(recovered.mspct[["one"]]), tested.location)
  expect_equal(getWhereMeasured(recovered.mspct[["two"]]), tested.location)
  expect_equal(getWhatMeasured(recovered.mspct[["one"]]), tested.what)
  expect_equal(getWhatMeasured(recovered.mspct[["two"]]), tested.what)
  expect_equal(getInstrDesc(recovered.mspct[["one"]]), tested.instr.desc)
  expect_equal(getInstrDesc(recovered.mspct[["two"]]), tested.instr.desc)
  expect_equal(getInstrSettings(recovered.mspct[["one"]]), tested.instr.settings)
  expect_equal(getInstrSettings(recovered.mspct[["two"]]), tested.instr.settings)
  expect_equal(is_normalized(recovered.mspct[["one"]]), normalization.flag)
  expect_equal(is_normalized(recovered.mspct[["two"]]), normalization.flag)
  expect_equal(is_scaled(recovered.mspct[["one"]]), scaling.flag)
  expect_equal(is_scaled(recovered.mspct[["two"]]), scaling.flag)

  expect_is(getInstrDesc(recovered.mspct[["one"]]), "instr_desc")
  expect_is(getInstrDesc(recovered.mspct[["two"]]), "instr_desc")
  expect_is(getInstrSettings(recovered.mspct[["one"]]), "instr_settings")
  expect_is(getInstrSettings(recovered.mspct[["two"]]), "instr_settings")
})

test_that("source_spct", {

  my.spct <- white_led.source_spct
  tested.location <- validate_geocode(data.frame(lon = 24.93545, lat = 60.16952))
  setWhereMeasured(my.spct, tested.location)
  tested.time <- getWhenMeasured(white_led.cps_spct)
  tested.what <- getWhatMeasured(my.spct)
  tested.instr.desc <- getInstrDesc(my.spct)
  tested.instr.settings <- getInstrSettings(my.spct)
  normalization.flag <- is_normalized(my.spct)
  scaling.flag <- is_scaled(my.spct)
  time.unit <- getTimeUnit(my.spct)

  spct.l <- list(one = my.spct, two = my.spct)

  expect_named(spct.l, c("one", "two"))
  expect_equal(getTimeUnit(spct.l[["one"]]), time.unit)
  expect_equal(getTimeUnit(spct.l[["two"]]), time.unit)
  expect_equal(getWhenMeasured(spct.l[["one"]]), tested.time)
  expect_equal(getWhenMeasured(spct.l[["two"]]), tested.time)
  expect_equal(getWhereMeasured(spct.l[["one"]]), tested.location)
  expect_equal(getWhereMeasured(spct.l[["two"]]), tested.location)
  expect_equal(getWhatMeasured(spct.l[["one"]]), tested.what)
  expect_equal(getWhatMeasured(spct.l[["two"]]), tested.what)
  expect_equal(getInstrDesc(spct.l[["one"]]), tested.instr.desc)
  expect_equal(getInstrDesc(spct.l[["two"]]), tested.instr.desc)
  expect_equal(getInstrSettings(spct.l[["one"]]), tested.instr.settings)
  expect_equal(getInstrSettings(spct.l[["two"]]), tested.instr.settings)
  expect_equal(is_normalized(spct.l[["one"]]), normalization.flag)
  expect_equal(is_normalized(spct.l[["two"]]), normalization.flag)
  expect_equal(is_scaled(spct.l[["one"]]), scaling.flag)
  expect_equal(is_scaled(spct.l[["two"]]), scaling.flag)

  ## collection
  my.mspct <- source_mspct(spct.l)

  expect_named(my.mspct, c("one", "two"))
  expect_equal(getTimeUnit(my.mspct[["one"]]), time.unit)
  expect_equal(getTimeUnit(my.mspct[["two"]]), time.unit)
  expect_equal(getWhenMeasured(my.mspct[["one"]]), tested.time)
  expect_equal(getWhenMeasured(my.mspct[["two"]]), tested.time)
  expect_equal(getWhereMeasured(my.mspct[["one"]]), tested.location)
  expect_equal(getWhereMeasured(my.mspct[["two"]]), tested.location)
  expect_equal(getWhatMeasured(my.mspct[["one"]]), tested.what)
  expect_equal(getWhatMeasured(my.mspct[["two"]]), tested.what)
  expect_equal(getInstrDesc(my.mspct[["one"]]), tested.instr.desc)
  expect_equal(getInstrDesc(my.mspct[["two"]]), tested.instr.desc)
  expect_equal(getInstrSettings(my.mspct[["one"]]), tested.instr.settings)
  expect_equal(getInstrSettings(my.mspct[["two"]]), tested.instr.settings)
  expect_equal(is_normalized(my.mspct[["one"]]), normalization.flag)
  expect_equal(is_normalized(my.mspct[["two"]]), normalization.flag)
  expect_equal(is_scaled(my.mspct[["one"]]), scaling.flag)
  expect_equal(is_scaled(my.mspct[["two"]]), scaling.flag)

  expect_is(getInstrDesc(my.mspct[["one"]]), "instr_desc")
  expect_is(getInstrDesc(my.mspct[["two"]]), "instr_desc")
  expect_is(getInstrSettings(my.mspct[["one"]]), "instr_settings")
  expect_is(getInstrSettings(my.mspct[["two"]]), "instr_settings")

  ## long form
  long.spct <- rbindspct(my.mspct)

  expect_named(getWhenMeasured(long.spct), c("one", "two"))
  expect_named(getWhatMeasured(long.spct), c("one", "two"))
  expect_equal(getWhereMeasured(long.spct)[["spct.idx"]], c("one", "two"))
  expect_named(getInstrDesc(long.spct), c("one", "two"))
  expect_named(getInstrSettings(long.spct), c("one", "two"))

  expect_equal(getTimeUnit(long.spct), time.unit)
  expect_equal(getWhenMeasured(long.spct)[["one"]], tested.time)
  expect_equal(getWhenMeasured(long.spct)[["two"]], tested.time)
  expect_equal(getWhereMeasured(long.spct)[1, -1], tested.location)
  expect_equal(getWhereMeasured(long.spct)[2, -1], tested.location)
  expect_equal(getWhatMeasured(long.spct)[["one"]], tested.what)
  expect_equal(getWhatMeasured(long.spct)[["two"]], tested.what)
  expect_equal(getInstrDesc(long.spct)[["one"]], tested.instr.desc)
  expect_equal(getInstrDesc(long.spct)[["two"]], tested.instr.desc)
  expect_equal(getInstrSettings(long.spct)[["one"]], tested.instr.settings)
  expect_equal(getInstrSettings(long.spct)[["two"]], tested.instr.settings)
  # expect_equal(is_normalized(long.spct)[["one"]], normalization.flag)
  # expect_equal(is_normalized(long.spct)[["two"]], normalization.flag)
  # expect_equal(is_scaled(long.spct)[["one"]], scaling.flag)
  # expect_equal(is_scaled(long.spct)[["two"]], scaling.flag)

  expect_true(!inherits(getInstrDesc(long.spct), "instr_desc"))
  expect_is(getInstrDesc(long.spct)[["one"]], "instr_desc")
  expect_is(getInstrDesc(long.spct)[["two"]], "instr_desc")
  expect_true(!inherits(getInstrSettings(long.spct), "instr_settings"))
  expect_is(getInstrSettings(long.spct)[["one"]], "instr_settings")
  expect_is(getInstrSettings(long.spct)[["two"]], "instr_settings")

  ## attributes
  zz <- rbindspct(sun_evening.mspct, attrs.simplify = TRUE)
  expect_equal(length(comment(zz)), 1L)
  expect_equal(nrow(where_measured(zz)), 1L)
  expect_equal(length(what_measured(zz)), 1L)
  expect_equal(length(how_measured(zz)), 1L)
  expect_equal(length(when_measured(zz)), length(sun_evening.mspct))
  zz <- rbindspct(sun_evening.mspct, attrs.simplify = FALSE)
  expect_equal(length(comment(zz)), 1L)
  expect_equal(nrow(where_measured(zz)), length(sun_evening.mspct))
  expect_equal(length(what_measured(zz)), length(sun_evening.mspct))
  expect_equal(length(how_measured(zz)), length(sun_evening.mspct))
  expect_equal(length(when_measured(zz)), length(sun_evening.mspct))

  ## idfactor
  zz <- rbindspct(sun_evening.mspct, idfactor = FALSE, attrs.simplify = FALSE)
  expect_equal(colnames(zz), c( "w.length", "s.e.irrad"))
  expect_true(setequal(colnames(where_measured(zz)), c("spct.idx", "lat", "lon", "address")))

  zz <- rbindspct(sun_evening.mspct, idfactor = TRUE, attrs.simplify = FALSE)
  expect_equal(colnames(zz), c( "w.length", "s.e.irrad", "spct.idx"))
  expect_true(setequal(colnames(where_measured(zz)), c("spct.idx", "lat", "lon", "address", "spct.idx")))

  ## recovered collection
  recovered.mspct <- subset2mspct(long.spct)

  expect_named(recovered.mspct, c("one", "two"))

  expect_equal(getTimeUnit(recovered.mspct[["one"]]), time.unit)
  expect_equal(getTimeUnit(recovered.mspct[["two"]]), time.unit)
  expect_equal(getWhenMeasured(recovered.mspct[["one"]]), tested.time)
  expect_equal(getWhenMeasured(recovered.mspct[["two"]]), tested.time)
  expect_equal(getWhereMeasured(recovered.mspct[["one"]]), tested.location)
  expect_equal(getWhereMeasured(recovered.mspct[["two"]]), tested.location)
  expect_equal(getWhatMeasured(recovered.mspct[["one"]]), tested.what)
  expect_equal(getWhatMeasured(recovered.mspct[["two"]]), tested.what)
  expect_equal(getInstrDesc(recovered.mspct[["one"]]), tested.instr.desc)
  expect_equal(getInstrDesc(recovered.mspct[["two"]]), tested.instr.desc)
  expect_equal(getInstrSettings(recovered.mspct[["one"]]), tested.instr.settings)
  expect_equal(getInstrSettings(recovered.mspct[["two"]]), tested.instr.settings)
  expect_equal(is_normalized(recovered.mspct[["one"]]), normalization.flag)
  expect_equal(is_normalized(recovered.mspct[["two"]]), normalization.flag)
  expect_equal(is_scaled(recovered.mspct[["one"]]), scaling.flag)
  expect_equal(is_scaled(recovered.mspct[["two"]]), scaling.flag)

  expect_is(getInstrDesc(recovered.mspct[["one"]]), "instr_desc")
  expect_is(getInstrDesc(recovered.mspct[["two"]]), "instr_desc")
  expect_is(getInstrSettings(recovered.mspct[["one"]]), "instr_settings")
  expect_is(getInstrSettings(recovered.mspct[["two"]]), "instr_settings")

  # long spct in a collection
  mixed.mspct <- source_mspct(list(a = long.spct, b = my.spct))
  expanded.mspct <- subset2mspct(mixed.mspct)

  expect_named(expanded.mspct, c("a.one", "a.two", "b"))

  expect_equal(getTimeUnit(expanded.mspct[["a.one"]]), time.unit)
  expect_equal(getTimeUnit(expanded.mspct[["a.two"]]), time.unit)
  expect_equal(getTimeUnit(expanded.mspct[["b"]]), time.unit)
  expect_equal(getWhenMeasured(expanded.mspct[["a.one"]]), tested.time)
  expect_equal(getWhenMeasured(expanded.mspct[["a.two"]]), tested.time)
  expect_equal(getWhenMeasured(expanded.mspct[["b"]]), tested.time)
  expect_equal(getWhereMeasured(expanded.mspct[["a.one"]]), tested.location)
  expect_equal(getWhereMeasured(expanded.mspct[["a.two"]]), tested.location)
  expect_equal(getWhereMeasured(expanded.mspct[["b"]]), tested.location)
  expect_equal(getWhatMeasured(expanded.mspct[["a.one"]]), tested.what)
  expect_equal(getWhatMeasured(expanded.mspct[["a.two"]]), tested.what)
  expect_equal(getWhatMeasured(expanded.mspct[["b"]]), tested.what)
  expect_equal(getInstrDesc(expanded.mspct[["a.one"]]), tested.instr.desc)
  expect_equal(getInstrDesc(expanded.mspct[["a.two"]]), tested.instr.desc)
  expect_equal(getInstrDesc(expanded.mspct[["b"]]), tested.instr.desc)
  expect_equal(getInstrSettings(expanded.mspct[["a.one"]]), tested.instr.settings)
  expect_equal(getInstrSettings(expanded.mspct[["a.two"]]), tested.instr.settings)
  expect_equal(getInstrSettings(expanded.mspct[["b"]]), tested.instr.settings)
  expect_equal(is_normalized(expanded.mspct[["a.one"]]), normalization.flag)
  expect_equal(is_normalized(expanded.mspct[["a.two"]]), normalization.flag)
  expect_equal(is_normalized(expanded.mspct[["b"]]), normalization.flag)
  expect_equal(is_scaled(expanded.mspct[["a.one"]]), scaling.flag)
  expect_equal(is_scaled(expanded.mspct[["a.two"]]), scaling.flag)
  expect_equal(is_scaled(expanded.mspct[["b"]]), scaling.flag)

  expect_is(getInstrDesc(expanded.mspct[["a.one"]]), "instr_desc")
  expect_is(getInstrDesc(expanded.mspct[["a.two"]]), "instr_desc")
  expect_is(getInstrDesc(expanded.mspct[["b"]]), "instr_desc")
  expect_is(getInstrSettings(expanded.mspct[["a.one"]]), "instr_settings")
  expect_is(getInstrSettings(expanded.mspct[["a.two"]]), "instr_settings")
  expect_is(getInstrSettings(expanded.mspct[["b"]]), "instr_settings")

  ## build long spctrum from mixed collection
  long3.spct <- rbindspct(mixed.mspct)

  spct.idx.rle <- rle(as.character(long3.spct$spct.idx))
  expect_equal(spct.idx.rle$lengths, rep(1421L, 3L))
  expect_equal(spct.idx.rle$values, c("a.one", "a.two", "b"))

  expect_equal(getTimeUnit(long3.spct), time.unit)
  expect_named(getWhenMeasured(long3.spct), c("a.one", "a.two", "b"))
  expect_equal(getWhenMeasured(long3.spct)[["a.one"]], tested.time)
  expect_equal(getWhenMeasured(long3.spct)[["a.two"]], tested.time)
  expect_equal(getWhenMeasured(long3.spct)[["b"]], tested.time)
  expect_equal(nrow(getWhereMeasured(long3.spct)), 3L)
  expect_equal(getWhereMeasured(long3.spct)[["spct.idx"]], c("a.one", "a.two", "b"))
  expect_named(getWhatMeasured(long3.spct), c("a.one", "a.two", "b"))
  expect_equal(getWhatMeasured(long3.spct)[["a.one"]], tested.what)
  expect_equal(getWhatMeasured(long3.spct)[["a.two"]], tested.what)
  expect_equal(getWhatMeasured(long3.spct)[["b"]], tested.what)
  expect_named(getInstrDesc(long3.spct), c("a.one", "a.two", "b"))
  expect_equal(getInstrDesc(long3.spct)[["a.one"]], tested.instr.desc)
  expect_equal(getInstrDesc(long3.spct)[["a.two"]], tested.instr.desc)
  expect_equal(getInstrDesc(long3.spct)[["b"]], tested.instr.desc)
  expect_named(getInstrSettings(long3.spct), c("a.one", "a.two", "b"))
  expect_equal(getInstrSettings(long3.spct)[["a.one"]], tested.instr.settings)
  expect_equal(getInstrSettings(long3.spct)[["a.two"]], tested.instr.settings)
  expect_equal(getInstrSettings(long3.spct)[["b"]], tested.instr.settings)
  expect_equal(is_normalized(long3.spct), normalization.flag)
  expect_equal(is_scaled(long3.spct), scaling.flag)

  expect_is(getInstrDesc(long3.spct)[["a.one"]], "instr_desc")
  expect_is(getInstrDesc(long3.spct)[["a.two"]], "instr_desc")
  expect_is(getInstrDesc(long3.spct)[["b"]], "instr_desc")
  expect_is(getInstrSettings(long3.spct)[["a.one"]], "instr_settings")
  expect_is(getInstrSettings(long3.spct)[["a.two"]], "instr_settings")
  expect_is(getInstrSettings(long3.spct)[["b"]], "instr_settings")
})

test_that("object_spct", {
  # We test only the specific attributes as all the rest of the code is
  # shared and so tested avobe.

  my.spct <- white_body.spct

  Tfr.type.test <- getTfrType(my.spct)
  Rfr.type.test <- getRfrType(my.spct)

  spct.l <- list(one = my.spct, two = my.spct)

  expect_named(spct.l, c("one", "two"))
  expect_equal(getTfrType(spct.l[["one"]]), Tfr.type.test)
  expect_equal(getTfrType(spct.l[["two"]]), Tfr.type.test)
  expect_equal(getRfrType(spct.l[["one"]]), Rfr.type.test)
  expect_equal(getRfrType(spct.l[["two"]]), Rfr.type.test)

  ## collection
  my.mspct <- object_mspct(spct.l)

  expect_named(my.mspct, c("one", "two"))
  expect_equal(getTfrType(my.mspct[["one"]]), Tfr.type.test)
  expect_equal(getTfrType(my.mspct[["two"]]), Tfr.type.test)
  expect_equal(getRfrType(my.mspct[["one"]]), Rfr.type.test)
  expect_equal(getRfrType(my.mspct[["two"]]), Rfr.type.test)

  ## long form
  long.spct <- rbindspct(my.mspct)

  expect_equal(getTfrType(long.spct), Tfr.type.test)
  expect_equal(getRfrType(long.spct), Rfr.type.test)

  ## recovered collection
  recovered.mspct <- subset2mspct(long.spct)

  expect_named(recovered.mspct, c("one", "two"))
  expect_equal(getTfrType(recovered.mspct[["one"]]), Tfr.type.test)
  expect_equal(getTfrType(recovered.mspct[["two"]]), Tfr.type.test)
  expect_equal(getRfrType(recovered.mspct[["one"]]), Rfr.type.test)
  expect_equal(getRfrType(recovered.mspct[["two"]]), Rfr.type.test)
})

test_that("filter_spct", {
  # We test only the specific attributes as all the rest of the code is
  # shared and so tested avobe.

  my.spct <- yellow_gel.spct

  Tfr.type.test <- getTfrType(my.spct)

  spct.l <- list(one = my.spct, two = my.spct)

  expect_named(spct.l, c("one", "two"))
  expect_equal(getTfrType(spct.l[["one"]]), Tfr.type.test)
  expect_equal(getTfrType(spct.l[["two"]]), Tfr.type.test)

  ## collection
  my.mspct <- filter_mspct(spct.l)

  expect_named(my.mspct, c("one", "two"))
  expect_equal(getTfrType(my.mspct[["one"]]), Tfr.type.test)
  expect_equal(getTfrType(my.mspct[["two"]]), Tfr.type.test)

  ## long form
  long.spct <- rbindspct(my.mspct)

  expect_equal(getTfrType(long.spct), Tfr.type.test)

  ## recovered collection
  recovered.mspct <- subset2mspct(long.spct)

  expect_named(recovered.mspct, c("one", "two"))
  expect_equal(getTfrType(recovered.mspct[["one"]]), Tfr.type.test)
  expect_equal(getTfrType(recovered.mspct[["two"]]), Tfr.type.test)
})

test_that("reflector_spct", {
  # We test only the specific attributes as all the rest of the code is
  # shared and so tested avobe.

  my.spct <- as.reflector_spct(white_body.spct)

  # temporary fix until data object is made valid
  setTfrType(my.spct, "internal")
  setRfrType(my.spct, "total")

  Tfr.type.test <- getTfrType(my.spct)
  Rfr.type.test <- getRfrType(my.spct)

  spct.l <- list(one = my.spct, two = my.spct)

  expect_named(spct.l, c("one", "two"))
  expect_equal(getRfrType(spct.l[["one"]]), Rfr.type.test)
  expect_equal(getRfrType(spct.l[["two"]]), Rfr.type.test)

  ## collection
  my.mspct <- reflector_mspct(spct.l)

  expect_named(my.mspct, c("one", "two"))
  expect_equal(getRfrType(my.mspct[["one"]]), Rfr.type.test)
  expect_equal(getRfrType(my.mspct[["two"]]), Rfr.type.test)

  ## long form
  long.spct <- rbindspct(my.mspct)

   expect_equal(getRfrType(long.spct), Rfr.type.test)

  ## recovered collection
  recovered.mspct <- subset2mspct(long.spct)

  expect_named(recovered.mspct, c("one", "two"))
  expect_equal(getRfrType(recovered.mspct[["one"]]), Rfr.type.test)
  expect_equal(getRfrType(recovered.mspct[["two"]]), Rfr.type.test)
})

