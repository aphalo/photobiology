library("lubridate")

context("attributes2tb")

test_that("tbNULL", {

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

  attributes.tb <- lon2tb(my.mspct)
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
  # tested.instr.desc <- getInstrDesc(my.spct)
  # tested.instr.settings <- getInstrSettings(my.spct)
  # normalization.flag <- is_normalized(my.spct)
  # scaling.flag <- is_scaled(my.spct)

  spct.l <- list(one = my.spct, two = my.spct)
  my.mspct <- source_mspct(spct.l)

  attributes.tb <- q_irrad(my.mspct)

  attributes.tb <- add_attr2tb(tb = attributes.tb,
                               mspct = my.mspct,
                               col.names = c("lon",
                                              "lat",
                                              "geocode",
                                              "when.measured",
                                              "what.measured"))

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

  # attributes.tb <- add_attr2tb(tb = attributes.tb,
  #                              mspct = my.mspct,
  #                              col.names = c("test1",
  #                                            "lon",
  #                                            "lat",
  #                                            "test2",
  #                                            "geocode",
  #                                            "when.measured",
  #                                            "what.measured"))
  # expect_named(my.mspct, c("one", "two"))
  # expect_is(attributes.tb[["spct.idx"]], "factor")
  # expect_equal(names(my.mspct), as.character(attributes.tb[["spct.idx"]]))
  #
  # expect_equal(getWhenMeasured(my.mspct[["one"]])[1],
  #              attributes.tb[["when.measured"]][1])
  # expect_equal(getWhenMeasured(my.mspct[["two"]])[1],
  #              attributes.tb[["when.measured"]][2])
  #
  # expect_equal(getWhatMeasured(my.mspct[["one"]])[1],
  #              attributes.tb[["what.measured"]][1])
  # expect_equal(getWhatMeasured(my.mspct[["two"]])[1],
  #              attributes.tb[["what.measured"]][2])
  #
  # expect_equal(getWhereMeasured(my.mspct[["one"]])[["lon"]],
  #              attributes.tb[["lon"]][1])
  # expect_equal(getWhereMeasured(my.mspct[["one"]])[["lat"]],
  #              attributes.tb[["lat"]][1])
  # expect_equal(getWhereMeasured(my.mspct[["one"]]),
  #              attributes.tb[["geocode"]][[1]])
  #
  # expect_equal(getWhereMeasured(my.mspct[["two"]])[["lon"]],
  #              attributes.tb[["lon"]][2])
  # expect_equal(getWhereMeasured(my.mspct[["two"]])[["lat"]],
  #              attributes.tb[["lat"]][2])
  # expect_equal(getWhereMeasured(my.mspct[["two"]]),
  #              attributes.tb[["geocode"]][[2]])

  attributes.tb <- q_irrad(my.mspct)

  expect_warning(add_attr2tb(tb = attributes.tb,
                             mspct = my.mspct,
                             col.names = c("test1",
                                           "test2")))

  # expect_equivalent(attributes.tb,
  #                   add_attr2tb(tb = attributes.tb,
  #                               mspct = my.mspct,
  #                               col.names = c("test1",
  #                                             "test2")))
})
