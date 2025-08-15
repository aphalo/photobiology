library("lubridate")

context("attributes2tb source_spct")

test_that("tbNULL", {

  my.spct <- white_led.source_spct
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.time <- getWhenMeasured(white_led.cps_spct)
  tested.what <- getWhatMeasured(my.spct)
  setHowMeasured(my.spct, "Array spectrometer")
  tested.how <- getHowMeasured(my.spct)
  tested.instr.desc <- getInstrDesc(my.spct)
  tested.instr.settings <- getInstrSettings(my.spct)
  normalization.flag <- is_normalized(my.spct)
  scaling.flag <- is_scaled(my.spct)
  tested.comment <- comment(my.spct)

  spct.l <- list(one = my.spct, two = my.spct)
  my.mspct <- source_mspct(spct.l)

  attributes.tb <- lon2tb(my.mspct)
  attributes.tb <- lat2tb(my.mspct, attributes.tb)
  attributes.tb <- geocode2tb(my.mspct, attributes.tb)
  attributes.tb <- when_measured2tb(my.mspct, attributes.tb)
  attributes.tb <- what_measured2tb(my.mspct, attributes.tb)
  attributes.tb <- how_measured2tb(my.mspct, attributes.tb)
  attributes.tb <- instr_desc2tb(my.mspct, attributes.tb)
  attributes.tb <- instr_settings2tb(my.mspct, attributes.tb)
  attributes.tb <- normalized2tb(my.mspct, attributes.tb)
  attributes.tb <- scaled2tb(my.mspct, attributes.tb)
  attributes.tb <- time_unit2tb(my.mspct, attributes.tb)
  attributes.tb <- multiple_wl2tb(my.mspct, attributes.tb)
  attributes.tb <- BSWF_used2tb(my.mspct, attributes.tb)

  attributes.tb <- comment2tb(my.mspct, attributes.tb)

  expect_named(my.mspct, c("one", "two"))
  expect_is(attributes.tb[["spct.idx"]], "factor")
  expect_equal(names(my.mspct), as.character(attributes.tb[["spct.idx"]]))

  expect_equal(getWhenMeasured(my.mspct[["one"]])[1],
               attributes.tb[["when.measured"]][1])
  expect_equal(getWhenMeasured(my.mspct[["two"]])[1],
               attributes.tb[["when.measured"]][2])

  expect_equal(getWhatMeasured(my.mspct[["one"]])[1],
               attributes.tb[["what.measured"]][1])
  expect_equal(getWhatMeasured(my.mspct[["two"]])[1],
               attributes.tb[["what.measured"]][2])

  expect_equal(getWhereMeasured(my.mspct[["one"]])[["lon"]],
               attributes.tb[["lon"]][1])
  expect_equal(getWhereMeasured(my.mspct[["one"]])[["lat"]],
               attributes.tb[["lat"]][1])
  expect_equal(getWhereMeasured(my.mspct[["one"]]),
               attributes.tb[["geocode"]][[1]])

  expect_equal(getHowMeasured(my.mspct[["one"]]),
               attributes.tb[["how.measured"]][[1]])
  expect_equal(getHowMeasured(my.mspct[["two"]]),
               attributes.tb[["how.measured"]][[2]])

  expect_equivalent(getInstrDesc(my.mspct[["one"]]),
               attributes.tb[["instr.desc"]][["one"]])
  expect_equivalent(getInstrDesc(my.mspct[["two"]]),
               attributes.tb[["instr.desc"]][["two"]])

  expect_equivalent(getInstrSettings(my.mspct[["one"]]),
               attributes.tb[["instr.settings"]][["one"]])
  expect_equivalent(getInstrSettings(my.mspct[["two"]]),
               attributes.tb[["instr.settings"]][["two"]])

  expect_equal(is_normalized(my.mspct[["one"]]),
               as.logical(attributes.tb[["normalized"]][[1]]))
  expect_equal(is_normalized(my.mspct[["two"]]),
               as.logical(attributes.tb[["normalized"]][[2]]))

  expect_equal(is_scaled(my.mspct[["one"]]) + 1,
               attributes.tb[["scaled"]][["one"]]$multiplier)
  expect_equal(is_scaled(my.mspct[["two"]]) + 1,
               attributes.tb[["scaled"]][["two"]]$multiplier)

  expect_equal(getTimeUnit(my.mspct[["one"]], force.duration = TRUE),
               attributes.tb[["time.unit"]][[1]])
  expect_equal(getTimeUnit(my.mspct[["two"]], force.duration = TRUE),
               attributes.tb[["time.unit"]][[2]])

  expect_equal(getMultipleWl(my.mspct[["one"]]),
               attributes.tb[["multiple.wl"]][[1]])
  expect_equal(getMultipleWl(my.mspct[["two"]]),
               attributes.tb[["multiple.wl"]][[2]])

  expect_equal(getBSWFUsed(my.mspct[["one"]]),
               attributes.tb[["bswf.used"]][[1]])
  expect_equal(getBSWFUsed(my.mspct[["two"]]),
               attributes.tb[["bswf.used"]][[2]])

  expect_equal(comment(my.mspct[["one"]]),
               attributes.tb[["comment"]][["one"]])
  expect_equal(comment(my.mspct[["two"]]),
               attributes.tb[["comment"]][["two"]])

})

test_that("spct_attributes", {

  my.spct <- white_led.source_spct
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.time <- getWhenMeasured(white_led.cps_spct)
  tested.what <- getWhatMeasured(my.spct)
  setHowMeasured(my.spct, "Array spectrometer")
  tested.how <- getHowMeasured(my.spct)
  tested.instr.desc <- getInstrDesc(my.spct)
  tested.instr.settings <- getInstrSettings(my.spct)
  normalization.flag <- is_normalized(my.spct)
  scaling.flag <- is_scaled(my.spct)
  tested.comment <- comment(my.spct)

  spct.l <- list(one = my.spct, two = my.spct)
  my.mspct <- source_mspct(spct.l)

  attributes.tb <- spct_metadata(my.mspct)

  expect_named(my.mspct, c("one", "two"))
  expect_is(attributes.tb[["spct.idx"]], "factor")
  expect_equal(names(my.mspct), as.character(attributes.tb[["spct.idx"]]))

  expect_equal(getWhenMeasured(my.mspct[["one"]])[1],
               attributes.tb[["when.measured"]][1])
  expect_equal(getWhenMeasured(my.mspct[["two"]])[1],
               attributes.tb[["when.measured"]][2])

  expect_equal(getWhatMeasured(my.mspct[["one"]])[1],
               attributes.tb[["what.measured"]][1])
  expect_equal(getWhatMeasured(my.mspct[["two"]])[1],
               attributes.tb[["what.measured"]][2])

  expect_equal(getWhereMeasured(my.mspct[["one"]])[["lon"]],
               attributes.tb[["lon"]][1])
  expect_equal(getWhereMeasured(my.mspct[["one"]])[["lat"]],
               attributes.tb[["lat"]][1])

  expect_equal(getWhereMeasured(my.mspct[["two"]])[["lon"]],
               attributes.tb[["lon"]][2])
  expect_equal(getWhereMeasured(my.mspct[["two"]])[["lat"]],
               attributes.tb[["lat"]][2])

  expect_equal(is_normalized(my.mspct[["one"]]),
               as.logical(attributes.tb[["normalized"]][[1]]))
  expect_equal(is_normalized(my.mspct[["two"]]),
               as.logical(attributes.tb[["normalized"]][[2]]))

  expect_equal(is_scaled(my.mspct[["one"]]) + 1,
               attributes.tb[["multiplier"]][[1]])
  expect_equal(is_scaled(my.mspct[["two"]]) + 1,
               attributes.tb[["multiplier"]][[2]])

  expect_equal(getTimeUnit(my.mspct[["one"]], force.duration = TRUE),
               attributes.tb[["time.unit"]][[1]])
  expect_equal(getTimeUnit(my.mspct[["two"]], force.duration = TRUE),
               attributes.tb[["time.unit"]][[2]])

  expect_equal(getBSWFUsed(my.mspct[["one"]]),
               attributes.tb[["bswf.used"]][[1]])
  expect_equal(getBSWFUsed(my.mspct[["two"]]),
               attributes.tb[["bswf.used"]][[2]])

})

test_that("tb_irrad", {

  my.spct <- white_led.source_spct
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.time <- getWhenMeasured(white_led.cps_spct)
  tested.what <- getWhatMeasured(my.spct)
  # tested.instr.desc <- getInstrDesc(my.spct)
  # tested.instr.settings <- getInstrSettings(my.spct)
  # normalization.flag <- is_normalized(my.spct)
  # scaling.flag <- is_scaled(my.spct)

  spct.l <- list(one = my.spct, two = my.spct)
  my.mspct <- source_mspct(spct.l)

  attributes.tb <- q_irrad(my.mspct)

  attributes.tb <- lon2tb(my.mspct, attributes.tb)
  attributes.tb <- lat2tb(my.mspct, attributes.tb)
  attributes.tb <- geocode2tb(my.mspct, attributes.tb)
  attributes.tb <- when_measured2tb(my.mspct, attributes.tb)
  attributes.tb <- what_measured2tb(my.mspct, attributes.tb)

  expect_named(my.mspct, c("one", "two"))
  expect_is(attributes.tb[["spct.idx"]], "factor")
  expect_equal(names(my.mspct), as.character(attributes.tb[["spct.idx"]]))

  expect_equal(getWhenMeasured(my.mspct[["one"]])[1],
               attributes.tb[["when.measured"]][1])
  expect_equal(getWhenMeasured(my.mspct[["two"]])[1],
               attributes.tb[["when.measured"]][2])

  expect_equal(getWhatMeasured(my.mspct[["one"]])[1],
               attributes.tb[["what.measured"]][1])
  expect_equal(getWhatMeasured(my.mspct[["two"]])[1],
               attributes.tb[["what.measured"]][2])

  expect_equal(getWhereMeasured(my.mspct[["one"]])[["lon"]],
               attributes.tb[["lon"]][1])
  expect_equal(getWhereMeasured(my.mspct[["one"]])[["lat"]],
               attributes.tb[["lat"]][1])
  expect_equal(getWhereMeasured(my.mspct[["one"]]),
               attributes.tb[["geocode"]][[1]])

  expect_equal(getWhereMeasured(my.mspct[["two"]])[["lon"]],
               attributes.tb[["lon"]][2])
  expect_equal(getWhereMeasured(my.mspct[["two"]])[["lat"]],
               attributes.tb[["lat"]][2])
  expect_equal(getWhereMeasured(my.mspct[["two"]]),
               attributes.tb[["geocode"]][[2]])

})

test_that("add_attr2tb", {

  my.spct <- white_led.source_spct
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.time <- getWhenMeasured(white_led.cps_spct)
  tested.what <- getWhatMeasured(my.spct)
  setHowMeasured(my.spct, "Array spectrometer")
  tested.how <- getHowMeasured(my.spct)
  tested.instr.desc <- getInstrDesc(my.spct)
  tested.instr.settings <- getInstrSettings(my.spct)
  normalization.flag <- is_normalized(my.spct)
  scaling.flag <- is_scaled(my.spct)
  tested.comment <- comment(my.spct)

  spct.l <- list(one = my.spct, two = my.spct)
  my.mspct <- source_mspct(spct.l)

  attributes.tb <- q_irrad(my.mspct)

  attributes.tb <- add_attr2tb(tb = attributes.tb,
                               mspct = my.mspct,
                               col.names = c("lon",
                                             "lat",
                                             "geocode",
                                             "when.measured",
                                             "what.measured",
                                             "how.measured"))

  expect_named(my.mspct, c("one", "two"))
  expect_is(attributes.tb[["spct.idx"]], "factor")
  expect_equal(names(my.mspct), as.character(attributes.tb[["spct.idx"]]))

  expect_equal(getWhenMeasured(my.mspct[["one"]])[1],
               attributes.tb[["when.measured"]][1])
  expect_equal(getWhenMeasured(my.mspct[["two"]])[1],
               attributes.tb[["when.measured"]][2])

  expect_equal(getWhatMeasured(my.mspct[["one"]])[1],
               attributes.tb[["what.measured"]][1])
  expect_equal(getWhatMeasured(my.mspct[["two"]])[1],
               attributes.tb[["what.measured"]][2])

  expect_equal(getWhereMeasured(my.mspct[["one"]])[["lon"]],
               attributes.tb[["lon"]][1])
  expect_equal(getWhereMeasured(my.mspct[["one"]])[["lat"]],
               attributes.tb[["lat"]][1])
  expect_equal(getWhereMeasured(my.mspct[["one"]]),
               attributes.tb[["geocode"]][[1]])

  expect_equal(getWhereMeasured(my.mspct[["two"]])[["lon"]],
               attributes.tb[["lon"]][2])
  expect_equal(getWhereMeasured(my.mspct[["two"]])[["lat"]],
               attributes.tb[["lat"]][2])
  expect_equal(getWhereMeasured(my.mspct[["two"]]),
               attributes.tb[["geocode"]][[2]])

  attributes.tb <- q_irrad(my.mspct)

  expect_warning(add_attr2tb(tb = attributes.tb,
                               mspct = my.mspct,
                               col.names = c("test1",
                                             "lon",
                                             "lat",
                                             "test2",
                                             "geocode",
                                             "when.measured",
                                             "what.measured")))

  attributes.tb <- q_irrad(my.mspct)

  attributes.tb <- add_attr2tb(tb = attributes.tb,
                               mspct = my.mspct,
                               col.names = c("where.measured",
                                             "when.measured",
                                             "what.measured",
                                             "how.measured",
                                             "normalized",
                                             "scaled",
                                             "comment",
                                             "time.unit",
                                             "multiple.wl",
                                             "bswf.used",
                                             "instr.sn",
                                             "instr.desc",
                                             "instr.settings"))

  expect_named(my.mspct, c("one", "two"))
  expect_is(attributes.tb[["spct.idx"]], "factor")
  expect_equal(names(my.mspct), as.character(attributes.tb[["spct.idx"]]))

  expect_equal(getWhenMeasured(my.mspct[["one"]])[1],
               attributes.tb[["when.measured"]][1])
  expect_equal(getWhenMeasured(my.mspct[["two"]])[1],
               attributes.tb[["when.measured"]][2])

  expect_equal(getWhatMeasured(my.mspct[["one"]])[1],
               attributes.tb[["what.measured"]][1])
  expect_equal(getWhatMeasured(my.mspct[["two"]])[1],
               attributes.tb[["what.measured"]][2])

  expect_equal(getWhereMeasured(my.mspct[["one"]]),
               attributes.tb[["where.measured"]][[1]])
  expect_equal(getWhereMeasured(my.mspct[["one"]]),
               attributes.tb[["where.measured"]][[1]])

  expect_equal(is_normalized(my.mspct[["one"]]),
               as.logical(attributes.tb[["normalized"]][[1]]))
  expect_equal(is_normalized(my.mspct[["two"]]),
               as.logical(attributes.tb[["normalized"]][[2]]))

  expect_equal(is_scaled(my.mspct[["one"]]) + 1,
               attributes.tb[["scaled"]][[1]]$multiplier)
  expect_equal(is_scaled(my.mspct[["two"]]) + 1,
               attributes.tb[["scaled"]][[2]]$multiplier)

  expect_equivalent(getInstrDesc(my.mspct[["one"]]),
                    attributes.tb[["instr.desc"]][["one"]])
  expect_equivalent(getInstrDesc(my.mspct[["two"]]),
                    attributes.tb[["instr.desc"]][["two"]])

  expect_equivalent(getInstrDesc(my.mspct[["one"]])["spectrometer.sn"],
                    attributes.tb[["instr.sn"]][[1]])
  expect_equivalent(getInstrDesc(my.mspct[["two"]])["spectrometer.sn"],
                    attributes.tb[["instr.sn"]][[1]])

  expect_equivalent(getInstrSettings(my.mspct[["one"]]),
                    attributes.tb[["instr.settings"]][["one"]])
  expect_equivalent(getInstrSettings(my.mspct[["two"]]),
                    attributes.tb[["instr.settings"]][["two"]])

  expect_equal(getTimeUnit(my.mspct[["one"]], force.duration = TRUE),
               attributes.tb[["time.unit"]][[1]])
  expect_equal(getTimeUnit(my.mspct[["two"]], force.duration = TRUE),
               attributes.tb[["time.unit"]][[2]])

  expect_equal(getMultipleWl(my.mspct[["one"]]),
               attributes.tb[["multiple.wl"]][[1]])
  expect_equal(getMultipleWl(my.mspct[["two"]]),
               attributes.tb[["multiple.wl"]][[2]])

  expect_equal(getBSWFUsed(my.mspct[["one"]]),
               attributes.tb[["bswf.used"]][[1]])
  expect_equal(getBSWFUsed(my.mspct[["two"]]),
               attributes.tb[["bswf.used"]][[2]])

  attributes.tb <- q_irrad(my.mspct)

  attributes.tb <- add_attr2tb(tb = attributes.tb,
                               mspct = my.mspct,
                               col.names = c("instr.desc"),
                               unnest = TRUE)

  expect_true(setequal(colnames(attributes.tb),
                       c("spct.idx", "Q_Total", "time", "spectrometer.name",
                         "spectrometer.sn", "bench.grating", "bench.slit")))

  attributes.tb <- q_irrad(my.mspct)

  attributes.tb <- add_attr2tb(tb = attributes.tb,
                               mspct = my.mspct,
                               col.names = c("instr.settings"),
                               unnest = TRUE)

  expect_true(setequal(colnames(attributes.tb),
                       c("spct.idx", "Q_Total", "pix.selector", "HDR.mult",
                         "target.margin", "max.integ.time", "min.integ.time",
                         "tot.time.range", "integ.time", "num.scans",
                         "corr.elect.dark", "corr.sensor.nl", "boxcar.width",
                         "linearized", "tot.time", "rel.signal")))

  attributes.tb <- ori.tb <- q_irrad(my.mspct)

  # bad column names
  expect_warning(attributes.tb <-
                   add_attr2tb(tb = attributes.tb,
                               mspct = my.mspct,
                               col.names = c("test1",
                                             "test2")))
  # no columns added
  expect_true(setequal(colnames(attributes.tb),
                       c("spct.idx", "Q_Total")))

  expect_equal(attributes.tb, ori.tb)

  attributes.tb <- q_irrad(my.mspct)

  expect_warning(attributes.tb <-
                   add_attr2tb(tb = attributes.tb,
                             mspct = my.mspct,
                             col.names = c("test1",
                                           "lon",
                                           "lat",
                                           "test2",
                                           "geocode",
                                           "when.measured",
                                           "what.measured")))

  # only good columns added
  expect_true(setequal(colnames(attributes.tb),
                       c("spct.idx", "Q_Total", "lon", "lat", "geocode",
                         "when.measured", "what.measured")))

})
