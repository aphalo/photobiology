context("summary")

update_all <- FALSE

test_that("source_spct", {

  energy_as_default()

  my.spct <- sun.spct
  expect_known_value(summary(my.spct),
                     "./data/summary-sun-default-e", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"),
                     "./data/summary-sun-none-e", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"),
                     "./data/summary-sun-collection-e", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"),
                     "./data/summary-sun-each-e", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection", which.metadata = "none"),
                     "./data/summary-sun-none-e-xcoll", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection", which.metadata = "all"),
                     "./data/summary-sun-all-e-xcoll", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection", which.metadata = "what.measured"),
                     "./data/summary-sun-what-e-xcoll", update = update_all)

  my.spct <- sun_evening.spct
  expect_known_value(summary(my.spct),
                     "./data/summary-eve-default-e", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"),
                     "./data/summary-eve-none-e", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"),
                     "./data/summary-eve-collection-e", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"),
                     "./data/summary-eve-each-e", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection", which.metadata = "none"),
                     "./data/summary-eve-none-e-xcoll", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection", which.metadata = "all"),
                     "./data/summary-eve-all-e-xcoll", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection", which.metadata = "what.measured"),
                     "./data/summary-eve-what-e-xcoll", update = update_all)

  photon_as_default()

  my.spct <- sun.spct
  expect_known_value(summary(my.spct), "./data/summary-sun-default-q", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-sun-none-q", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-sun-collection-q", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-sun-each-q", update = update_all)

  my.spct <- sun_evening.spct
  expect_known_value(summary(my.spct), "./data/summary-eve-default-q", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-eve-none-q", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-eve-collection-q", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-eve-each-q", update = update_all)

  unset_radiation_unit_default()
})

test_that("response_spct", {

  my.spct <- q2e(ccd.spct)

  energy_as_default()
  expect_known_value(summary(my.spct), "./data/summary-default-re", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-default-none-re", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-collection-re", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-each-re", update = update_all)

  photon_as_default()
  expect_known_value(summary(my.spct), "./data/summary-default-rq", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-default-none-rq", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-collection-rq", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-each-rq", update = update_all)

  unset_radiation_unit_default()
})

test_that("filter_spct", {

  my.spct <- yellow_gel.spct

  Tfr_as_default()
  expect_known_value(summary(my.spct), "./data/summary-y-default-tfr", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-y-sun-none-tfr", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-y-collection-tfr", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-y-each-tfr", update = update_all)

  Afr_as_default()
  expect_known_value(summary(my.spct), "./data/summary-y-default-afr", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-y-sun-none-afr", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-y-collection-afr", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-y-each-afr", update = update_all)

  A_as_default()
  expect_known_value(summary(my.spct), "./data/summary-y-default-a", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-y-sun-none-a", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-y-collection-a", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-y-each-a", update = update_all)

  unset_filter_qty_default()

  my.spct <- two_filters.spct

  Tfr_as_default()
  expect_known_value(summary(my.spct),
                     "./data/summary-2f-default-tfr", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"),
                     "./data/summary-2f-sun-none-tfr", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"),
                     "./data/summary-2f-collection-tfr", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"),
                     "./data/summary-2f-each-tfr", update = update_all)
  expect_known_value(summary(my.spct, expand = "each", which.metadtata = "none"),
                     "./data/summary-2f-each-none-tfr", update = update_all)
  expect_known_value(summary(my.spct, expand = "each", which.metadtata = "all"),
                     "./data/summary-2f-each-all-tfr", update = update_all)
  expect_known_value(summary(my.spct, expand = "each", which.metadtata = "what.measured"),
                     "./data/summary-2f-each-what-tfr", update = update_all)

  Afr_as_default()
  expect_known_value(summary(my.spct), "./data/summary-2f-default-afr", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-2f-sun-none-afr", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-2f-collection-afr", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-2f-each-afr", update = update_all)

  A_as_default()
  expect_known_value(summary(my.spct), "./data/summary-2f-default-a", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-2f-sun-none-a", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-2f-collection-a", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-2f-each-a", update = update_all)

  unset_filter_qty_default()

 })

test_that("reflector_spct", {

  my.spct <- Ler_leaf_rflt.spct

  expect_known_value(summary(my.spct), "./data/summary-default-rfr", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-none-rfr", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-collection-rfr", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-each-rfr", update = update_all)
})

test_that("object_spct", {

  my.spct <- Ler_leaf.spct

  expect_known_value(summary(my.spct), "./data/summary-default-objt", update = update_all)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-none-objt", update = update_all)
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-collection-objt", update = update_all)
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-each-objt", update = update_all)
})
