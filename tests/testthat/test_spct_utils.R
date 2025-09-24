context("thin_wl")

update_all <- FALSE

test_that("source_spct", {

  my.spct <- sun.spct

  energy_as_default()
  expect_equal(nrow(thin_wl(my.spct)), 508L)
  expect_known_value(thin_wl(my.spct), "./data/wl-thin-default-value-e", update = update_all)
  expect_known_value(thin_wl(my.spct, max.wl.step = 40), "./data/wl-thin-default-value-step-e", update = update_all)
  expect_known_value(thin_wl(my.spct, max.slope.delta = 0.003), "./data/wl-thin-default-value-slope-e", update = update_all)

  photon_as_default()
  expect_equal(nrow(thin_wl(my.spct)), 507L)
  expect_known_value(thin_wl(my.spct), "./data/wl-thin-default-value-q", update = update_all)
  expect_known_value(thin_wl(my.spct, max.wl.step = 40), "./data/wl-thin-default-value-step-q", update = update_all)
  expect_known_value(thin_wl(my.spct, max.slope.delta = 0.003), "./data/wl-thin-default-value-slope-q", update = update_all)

  unset_radiation_unit_default()
})

test_that("response_spct", {

  my.spct <- ccd.spct

  energy_as_default()
  expect_equal(nrow(thin_wl(my.spct)), 159L)
  expect_known_value(thin_wl(my.spct), "./data/wl-thin-default-value-re", update = update_all)
  expect_known_value(thin_wl(my.spct, max.wl.step = 40), "./data/wl-thin-default-value-step-re", update = update_all)
  expect_known_value(thin_wl(my.spct, max.slope.delta = 0.003), "./data/wl-thin-default-value-slope-re", update = update_all)

  photon_as_default()
  expect_equal(nrow(thin_wl(my.spct)), 161L)
  expect_known_value(thin_wl(my.spct), "./data/wl-thin-default-value-rq", update = update_all)
  expect_known_value(thin_wl(my.spct, max.wl.step = 40), "./data/wl-thin-default-value-step-rq", update = update_all)
  expect_known_value(thin_wl(my.spct, max.slope.delta = 0.003), "./data/wl-thin-default-value-slope-rq", update = update_all)

  unset_radiation_unit_default()
})

test_that("filter_spct", {

  my.spct <- yellow_gel.spct

  Tfr_as_default()
  expect_equal(nrow(thin_wl(my.spct)), 344L)
  expect_known_value(thin_wl(my.spct), "./data/wl-thin-default-value-tfr", update = update_all)
  expect_known_value(thin_wl(my.spct, max.wl.step = 40), "./data/wl-thin-default-value-step-tfr", update = update_all)
  expect_known_value(thin_wl(my.spct, max.slope.delta = 0.003), "./data/wl-thin-default-value-slope-tfr", update = update_all)

  unset_filter_qty_default()
})

test_that("reflector_spct", {

  my.spct <- Ler_leaf_rflt.spct

  expect_equal(nrow(thin_wl(my.spct)), 580L)
  expect_known_value(thin_wl(my.spct), "./data/wl-thin-default-value-rfr", update = update_all)
  expect_known_value(thin_wl(my.spct, max.wl.step = 40), "./data/wl-thin-default-value-step-rfr", update = update_all)
  expect_known_value(thin_wl(my.spct, max.slope.delta = 0.003), "./data/wl-thin-default-value-slope-rfr", update = update_all)
})

test_that("object_spct", {

  my.spct <- Ler_leaf.spct

  expect_error(thin_wl(my.spct))
  expect_equal(nrow(thin_wl(my.spct, col.names = "Tfr")), 641L)
  expect_known_value(thin_wl(my.spct, col.names = "Tfr"), "./data/wl-thin-default-value-otfr", update = update_all)
  expect_known_value(thin_wl(my.spct, col.names = "Rfr"), "./data/wl-thin-default-value-orfr", update = update_all)
})

test_that("chroma_spct", {

  my.spct <- ciexyzCC10.spct

  expect_warning(thin_wl(my.spct))
  expect_equal(suppressWarnings(thin_wl(my.spct)), my.spct)
})

test_that("calibration_spct", {

  my.spct <- calibration_spct(w.length = 400:450, irrad.mult = 1)

  expect_warning(thin_wl(my.spct))
  expect_equal(suppressWarnings(thin_wl(my.spct)), my.spct)
})

test_that("generic_mspct", {

  my.mspct <- source_mspct(list(sun.spct))

  energy_as_default()
  expect_equal(nrow(thin_wl(my.mspct)[[1]]), 508L)
  expect_known_value(thin_wl(my.mspct), "./data/wl-thin-default-value-mspct-e", update = update_all)

  photon_as_default()
  expect_equal(nrow(thin_wl(my.mspct)[[1]]), 507L)
  expect_known_value(thin_wl(my.mspct), "./data/wl-thin-default-value-mspct-q", update = update_all)

  unset_radiation_unit_default()
})
