library("photobiology")
library("lubridate")

context("wide_long_wide")

test_that("cps_spct", {

  my.spct <- white_led.cps_spct
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
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

  expect_is(getInstrDesc(my.mspct[["one"]]), "instr_desc")
  expect_is(getInstrDesc(my.mspct[["two"]]), "instr_desc")
  expect_is(getInstrSettings(my.mspct[["one"]]), "instr_settings")
  expect_is(getInstrSettings(my.mspct[["two"]]), "instr_settings")

  ## long form
  long.spct <- rbindspct(my.mspct)

  expect_named(getWhenMeasured(long.spct), c("one", "two"))
  expect_named(getWhatMeasured(long.spct), c("one", "two"))
  expect_named(getWhereMeasured(long.spct), c("one", "two"))
  expect_named(getInstrDesc(long.spct), c("one", "two"))
  expect_named(getInstrSettings(long.spct), c("one", "two"))

  expect_equal(getWhenMeasured(long.spct)[["one"]], tested.time)
  expect_equal(getWhenMeasured(long.spct)[["two"]], tested.time)
  expect_equal(getWhereMeasured(long.spct)[["one"]], tested.location)
  expect_equal(getWhereMeasured(long.spct)[["two"]], tested.location)
  expect_equal(getWhatMeasured(long.spct)[["one"]], tested.what)
  expect_equal(getWhatMeasured(long.spct)[["two"]], tested.what)
  expect_equal(getInstrDesc(long.spct)[["one"]], tested.instr.desc)
  expect_equal(getInstrDesc(long.spct)[["two"]], tested.instr.desc)
  expect_equal(getInstrSettings(long.spct)[["one"]], tested.instr.settings)
  expect_equal(getInstrSettings(long.spct)[["two"]], tested.instr.settings)

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

})
