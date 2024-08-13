context("summary")

test_that("source_spct", {

  energy_as_default()

  my.spct <- sun.spct
  expect_known_value(summary(my.spct), "./data/summary-sun-default-e", update = FALSE)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-sun-none-e")
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-sun-collection-e")
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-sun-each-e")

  my.spct <- sun_evening.spct
  expect_known_value(summary(my.spct), "./data/summary-eve-default-e", update = FALSE)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-eve-none-e")
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-eve-collection-e")
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-eve-each-e")

  # photon_as_default()
  #
  # my.spct <- sun.spct
  # expect_known_value(summary(my.spct), "./data/summary-sun-default-q", update = FALSE)
  # expect_known_value(summary(my.spct, expand = "none"), "./data/summary-sun-none-q")
  # expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-sun-collection-q")
  # expect_known_value(summary(my.spct, expand = "each"), "./data/summary-sun-each-q")
  #
  # my.spct <- sun_evening.spct
  # expect_known_value(summary(my.spct), "./data/summary-eve-default-q", update = FALSE)
  # expect_known_value(summary(my.spct, expand = "none"), "./data/summary-eve-none-q")
  # expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-eve-collection-q")
  # expect_known_value(summary(my.spct, expand = "each"), "./data/summary-eve-each-q")

  unset_radiation_unit_default()
})

test_that("response_spct", {

  my.spct <- ccd.spct

  energy_as_default()
  expect_known_value(summary(my.spct), "./data/summary-default-re")
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-default-none-re")
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-collection-re")
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-each-re")

  # photon_as_default()
  # expect_known_value(summary(my.spct), "./data/summary-default-rq")
  # expect_known_value(summary(my.spct, expand = "none"), "./data/summary-default-none-rq")
  # expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-collection-rq")
  # expect_known_value(summary(my.spct, expand = "each"), "./data/summary-each-rq")

  unset_radiation_unit_default()
})

test_that("filter_spct", {

  my.spct <- yellow_gel.spct

  Tfr_as_default()
  expect_known_value(summary(my.spct), "./data/summary-y-default-tfr", update = FALSE)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-y-sun-none-tfr")
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-y-collection-tfr")
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-y-each-tfr")

  # Afr_as_default()
  # expect_known_value(summary(my.spct), "./data/summary-y-default-afr", update = FALSE)
  # expect_known_value(summary(my.spct, expand = "none"), "./data/summary-y-sun-none-afr")
  # expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-y-collection-afr")
  # expect_known_value(summary(my.spct, expand = "each"), "./data/summary-y-each-afr")
  #
  # A_as_default()
  # expect_known_value(summary(my.spct), "./data/summary-y-default-a", update = FALSE)
  # expect_known_value(summary(my.spct, expand = "none"), "./data/summary-y-sun-none-a")
  # expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-y-collection-a")
  # expect_known_value(summary(my.spct, expand = "each"), "./data/summary-y-each-a")
  #
  unset_filter_qty_default()

  my.spct <- two_filters.spct

  Tfr_as_default()
  expect_known_value(summary(my.spct), "./data/summary-2f-default-tfr", update = FALSE)
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-2f-sun-none-tfr")
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-2f-collection-tfr")
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-2f-each-tfr")

  # Afr_as_default()
  # expect_known_value(summary(my.spct), "./data/summary-2f-default-afr", update = FALSE)
  # expect_known_value(summary(my.spct, expand = "none"), "./data/summary-2f-sun-none-afr")
  # expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-2f-collection-afr")
  # expect_known_value(summary(my.spct, expand = "each"), "./data/summary-2f-each-afr")
  #
  # A_as_default()
  # expect_known_value(summary(my.spct), "./data/summary-2f-default-a", update = FALSE)
  # expect_known_value(summary(my.spct, expand = "none"), "./data/summary-2f-sun-none-a")
  # expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-2f-collection-a")
  # expect_known_value(summary(my.spct, expand = "each"), "./data/summary-2f-each-a")
  #
  unset_filter_qty_default()

 })

test_that("reflector_spct", {

  my.spct <- Ler_leaf_rflt.spct

  expect_known_value(summary(my.spct), "./data/summary-default-rfr")
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-none-rfr")
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-collection-rfr")
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-each-rfr")
})

test_that("object_spct", {

  my.spct <- Ler_leaf.spct

  expect_known_value(summary(my.spct), "./data/summary-default-objt")
  expect_known_value(summary(my.spct, expand = "none"), "./data/summary-none-objt")
  expect_known_value(summary(my.spct, expand = "collection"), "./data/summary-collection-objt")
  expect_known_value(summary(my.spct, expand = "each"), "./data/summary-each-objt")
})
