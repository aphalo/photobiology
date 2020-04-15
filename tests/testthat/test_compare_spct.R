context("compare_spct")

test_that("source_mspct", {

  my.mspct <- source_mspct(list(sun1 = sun.spct, sun2 = sun.spct * 2))
  my.mspct <- clip_wl(my.mspct, range = c(400, 450)) # make tests faster

  energy_as_default()
  expect_known_value(compare_spct(my.mspct), "./data/compare-spct-default-value-e")
  expect_named(compare_spct(my.mspct),
               c("w.length", "wl.min", "wl.max", "sun1.irrad", "sun2.irrad",
                 "comparison.result"))
  expect_equal(unique(na.omit(compare_spct(my.mspct)[["comparison.result"]])), 2)
  expect_is(unique(na.omit(compare_spct(my.mspct)[["comparison.result"]])), "numeric")
  expect_is(unique(na.omit(compare_spct(my.mspct, .comparison.fun = `>`)[["comparison.result"]])), "logical")

  default.spct <- compare_spct(my.mspct)
  expect_equal(nrow(default.spct), 5)
  expect_known_value(compare_spct(my.mspct, returned.value = "data.frame"), "./data/compare-spct-df-value-e")
  expect_known_value(compare_spct(my.mspct, returned.value = "spectrum"), "./data/compare-spct-spct-value-e")
  expect_known_value(compare_spct(my.mspct, returned.value = "tagged.spectrum"), "./data/compare-spct-tag-value-e")
  expect_warning(compare_spct(my.mspct, returned.value = "zzz"))
  expect_equal(suppressWarnings(compare_spct(my.mspct, returned.value = "zzz")),
               compare_spct(my.mspct, returned.value = "data.frame"))

  photon_as_default()
  expect_known_value(compare_spct(my.mspct), "./data/compare-spct-default-value-q")
  expect_named(compare_spct(my.mspct),
               c("w.length", "wl.min", "wl.max", "sun1.irrad", "sun2.irrad",
                 "comparison.result"))
  expect_equal(unique(na.omit(compare_spct(my.mspct)[["comparison.result"]])), 2)
  expect_is(unique(na.omit(compare_spct(my.mspct)[["comparison.result"]])), "numeric")
  expect_is(unique(na.omit(compare_spct(my.mspct, .comparison.fun = `>`)[["comparison.result"]])), "logical")

  default.spct <- compare_spct(my.mspct)
  expect_equal(nrow(default.spct), 5)
  expect_known_value(compare_spct(my.mspct, returned.value = "data.frame"), "./data/compare-spct-df-value-e2")
  expect_known_value(compare_spct(my.mspct, returned.value = "spectrum"), "./data/compare-spct-spct-value-e2")
  expect_known_value(compare_spct(my.mspct, returned.value = "tagged.spectrum"), "./data/compare-spct-tag-value-e2")
  expect_warning(compare_spct(my.mspct, returned.value = "zzz"))
  expect_equal(suppressWarnings(compare_spct(my.mspct, returned.value = "zzz")),
               compare_spct(my.mspct, returned.value = "data.frame"))

  unset_radiation_unit_default()

})

test_that("response_mspct", {

  energy_as_default()
  my.mspct <- response_mspct(list(ccd1 = ccd.spct, ccd2 = ccd.spct * 2))
  my.mspct <- clip_wl(my.mspct, range = c(400, 450)) # make tests faster

  expect_known_value(compare_spct(my.mspct), "./data/compare-spct-default-value-re")
  expect_named(compare_spct(my.mspct),
               c("w.length", "wl.min", "wl.max", "ccd1.response", "ccd2.response",
                 "comparison.result"))
  expect_equal(unique(na.omit(compare_spct(my.mspct)[["comparison.result"]])), 2)
  expect_is(unique(na.omit(compare_spct(my.mspct)[["comparison.result"]])), "numeric")
  expect_is(unique(na.omit(compare_spct(my.mspct, .comparison.fun = `>`)[["comparison.result"]])), "logical")

  photon_as_default()
  my.mspct <- response_mspct(list(ccd1 = ccd.spct, ccd2 = ccd.spct * 2))
  my.mspct <- clip_wl(my.mspct, range = c(400, 450)) # make tests faster

  expect_known_value(compare_spct(my.mspct), "./data/compare-spct-default-value-rq")
  expect_named(compare_spct(my.mspct),
               c("w.length", "wl.min", "wl.max", "ccd1.response", "ccd2.response",
                 "comparison.result"))
  expect_equal(unique(na.omit(compare_spct(my.mspct)[["comparison.result"]])), 2)
  expect_is(unique(na.omit(compare_spct(my.mspct)[["comparison.result"]])), "numeric")
  expect_is(unique(na.omit(compare_spct(my.mspct, .comparison.fun = `>`)[["comparison.result"]])), "logical")

  unset_radiation_unit_default()

})

test_that("filter_mspct", {

  Tfr_as_default()
  my.mspct <- filter_mspct(list(pet1 = clean(polyester.spct), pet2 = clean(polyester.spct) / 2))
  my.mspct <- clip_wl(my.mspct, range = c(400, 450)) # make tests faster

  expect_known_value(compare_spct(my.mspct), "./data/compare-spct-default-value-tfr")
  expect_named(compare_spct(my.mspct),
               c("w.length", "wl.min", "wl.max", "pet1.transmittance", "pet2.transmittance",
                 "comparison.result"))
  expect_equal(unique(na.omit(compare_spct(my.mspct)[["comparison.result"]])), 1/2)
  expect_is(unique(na.omit(compare_spct(my.mspct)[["comparison.result"]])), "numeric")
  expect_is(unique(na.omit(compare_spct(my.mspct, .comparison.fun = `>`)[["comparison.result"]])), "logical")

  A_as_default()
  my.mspct <- filter_mspct(list(pet1 = clean(polyester.spct), pet2 = clean(polyester.spct) * 2))
  my.mspct <- clip_wl(my.mspct, range = c(400, 450)) # make tests faster

  expect_known_value(compare_spct(my.mspct, .summary.fun = absorbance), "./data/compare-spct-default-value-A")
  expect_named(compare_spct(my.mspct, .summary.fun = absorbance),
               c("w.length", "wl.min", "wl.max", "pet1.absorbance", "pet2.absorbance",
                 "comparison.result"))
  expect_equal(unique(na.omit(compare_spct(my.mspct, .summary.fun = absorbance)[["comparison.result"]])), 2)
  expect_is(unique(na.omit(compare_spct(my.mspct, .summary.fun = absorbance)[["comparison.result"]])), "numeric")
  expect_is(unique(na.omit(compare_spct(my.mspct, .summary.fun = absorbance, .comparison.fun = `>`)[["comparison.result"]])), "logical")

  unset_filter_qty_default()

})


test_that("reflector_mspct", {

  my.mspct <- reflector_mspct(list(ler1 = Ler_leaf_rflt.spct, ler2 = Ler_leaf_rflt.spct / 2))
  my.mspct <- clip_wl(my.mspct, range = c(400, 450)) # make tests faster

  expect_known_value(compare_spct(my.mspct), "./data/compare-spct-default-value-Rfr")
  expect_named(compare_spct(my.mspct),
               c("w.length", "wl.min", "wl.max", "ler1.reflectance", "ler2.reflectance",
                 "comparison.result"))
  expect_equal(unique(na.omit(compare_spct(my.mspct)[["comparison.result"]])), 1/2)
  expect_is(unique(na.omit(compare_spct(my.mspct)[["comparison.result"]])), "numeric")
  expect_is(unique(na.omit(compare_spct(my.mspct, .comparison.fun = `>`)[["comparison.result"]])), "logical")

})



