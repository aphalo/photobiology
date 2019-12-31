context("wl_thin")

test_that("source_spct", {

  my.spct <- sun.spct

  energy_as_default()
  expect_known_value(wl_thin(my.spct), "./data/wl-thin-default-value-e")
  expect_known_value(wl_thin(my.spct, max.wl.step = 40), "./data/wl-thin-default-value-step-e")
  expect_known_value(wl_thin(my.spct, max.slope.delta = 0.003), "./data/wl-thin-default-value-slope-e")

  photon_as_default()
  expect_known_value(wl_thin(my.spct), "./data/wl-thin-default-value-q")
  expect_known_value(wl_thin(my.spct, max.wl.step = 40), "./data/wl-thin-default-value-step-q")
  expect_known_value(wl_thin(my.spct, max.slope.delta = 0.003), "./data/wl-thin-default-value-slope-q")

  unset_radiation_unit_default()
})

test_that("response_spct", {

  my.spct <- ccd.spct

  energy_as_default()
  expect_known_value(wl_thin(my.spct), "./data/wl-thin-default-value-re")
  expect_known_value(wl_thin(my.spct, max.wl.step = 40), "./data/wl-thin-default-value-step-re")
  expect_known_value(wl_thin(my.spct, max.slope.delta = 0.003), "./data/wl-thin-default-value-slope-re")

  photon_as_default()
  expect_known_value(wl_thin(my.spct), "./data/wl-thin-default-value-rq")
  expect_known_value(wl_thin(my.spct, max.wl.step = 40), "./data/wl-thin-default-value-step-rq")
  expect_known_value(wl_thin(my.spct, max.slope.delta = 0.003), "./data/wl-thin-default-value-slope-rq")

  unset_radiation_unit_default()
})

test_that("filter_spct", {

  my.spct <- yellow_gel.spct

  Tfr_as_default()
  expect_known_value(wl_thin(my.spct), "./data/wl-thin-default-value-tfr")
  expect_known_value(wl_thin(my.spct, max.wl.step = 40), "./data/wl-thin-default-value-step-tfr")
  expect_known_value(wl_thin(my.spct, max.slope.delta = 0.003), "./data/wl-thin-default-value-slope-tfr")

  unset_filter_qty_default()
})

test_that("reflector_spct", {

  my.spct <- Ler_leaf_rflt.spct

  expect_known_value(wl_thin(my.spct), "./data/wl-thin-default-value-rfr")
  expect_known_value(wl_thin(my.spct, max.wl.step = 40), "./data/wl-thin-default-value-step-rfr")
  expect_known_value(wl_thin(my.spct, max.slope.delta = 0.003), "./data/wl-thin-default-value-slope-rfr")
})

test_that("object_spct", {

  my.spct <- Ler_leaf.spct

  expect_error(wl_thin(my.spct))
  expect_known_value(wl_thin(my.spct, col.names = "Tfr"), "./data/wl-thin-default-value-otfr")
  expect_known_value(wl_thin(my.spct, col.names = "Rfr"), "./data/wl-thin-default-value-orfr")
})

test_that("chroma_spct", {

  my.spct <- ciexyzCC10.spct

  expect_warning(wl_thin(my.spct))
  expect_equal(suppressWarnings(wl_thin(my.spct)), my.spct)
})

test_that("calibration_spct", {

  my.spct <- calibration_spct(w.length = 400:450, irrad.mult = 1)

  expect_warning(wl_thin(my.spct))
  expect_equal(suppressWarnings(wl_thin(my.spct)), my.spct)
})
