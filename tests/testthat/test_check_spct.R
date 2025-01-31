context("source_spct")

options(photobiology.strict.range = TRUE)

test_that("check_spct works for source_spct", {
  # tests in test_data.r test for false positives implicitly
  # here we test for failing tests

  empty.spct <- source_spct()
  expect_silent(check_spct(empty.spct))
  bad.spct <- reference.spct <- e2q(sun.spct, action = "add")
  expect_error(bad.spct[1, "w.length"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 0)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 1000)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[length(bad.spct), "w.length"] <- 100)
  expect_equal(bad.spct, reference.spct)
  # negative fluxes are valid as they depend on direction taken as reference
  expect_silent(bad.spct[1, "s.e.irrad"] <- -1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "s.e.irrad"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "s.q.irrad"] <- -1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "s.q.irrad"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
})

context("response_spct")

test_that("check_spct works for response_spct", {
  # tests in test_data.r test for false positives implicitly
  # here we test for failing tests

  empty.spct <- response_spct()
  expect_silent(check_spct(empty.spct))
  bad.spct <- reference.spct <- e2q(photodiode.spct, action = "add")
  expect_error(bad.spct[1, "w.length"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 0)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 1000)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[length(bad.spct), "w.length"] <- 100)
  expect_equal(bad.spct, reference.spct)
  # negative responses are valid as they depend on what is taken as reference
  expect_silent(bad.spct[1, "s.e.response"] <- -1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "s.e.response"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "s.q.response"] <- -1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "s.q.response"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
})

context("filter_spct")

test_that("check_spct works for filter_spct", {
  # tests in test_data.r test for false positives implicitly
  # here we test for failing tests

  empty.spct <- filter_spct()
  expect_silent(check_spct(empty.spct))

  bad.spct <- reference.spct <- polyester.spct
  expect_error(bad.spct[1, "w.length"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 0)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 1000)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[length(bad.spct), "w.length"] <- 100)
  expect_equal(bad.spct, reference.spct)
  # negative transmittances are invalid, valid range is 0..1
  expect_error(bad.spct[1, "Tfr"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1:(nrow(bad.spct)/250 + 2), "Tfr"] <- -2e-4)
  expect_equal(bad.spct, reference.spct)
  expect_silent(bad.spct[1, "Tfr"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "Tfr"] <- 1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])

  bad.spct <- reference.spct <- any2Afr(polyester.spct, action = "replace")
  expect_error(bad.spct[1, "w.length"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 0)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 1000)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[length(bad.spct), "w.length"] <- 100)
  expect_equal(bad.spct, reference.spct)
  # negative transmittances are invalid, valid range is 0..1
  expect_error(bad.spct[1, "Afr"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1:(nrow(bad.spct)/250 + 2), "Afr"] <- -2e-4)
  expect_equal(bad.spct, reference.spct)
  expect_silent(bad.spct[1, "Afr"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "Afr"] <- 1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])

  bad.spct <- reference.spct <- any2A(polyester.spct, action = "replace")
  expect_error(bad.spct[1, "w.length"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 0)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 1000)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[length(bad.spct), "w.length"] <- 100)
  expect_equal(bad.spct, reference.spct)
  # negative absorbances are valid
  expect_error(bad.spct[1, "A"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1:(nrow(bad.spct)/250 + 2), "A"] <- -2e-7)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1:(nrow(bad.spct)/250 + 2), "A"] <- 21)
  expect_equal(bad.spct, reference.spct)
  expect_silent(bad.spct[1, "A"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "A"] <- 1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
})

context("reflector_spct")

test_that("check_spct works for reflector_spct", {
  # tests in test_data.r test for false positives implicitly
  # here we test for failing tests

  empty.spct <- reflector_spct()
  expect_silent(check_spct(empty.spct))

  bad.spct <- reference.spct <- green_leaf.spct
  expect_error(bad.spct[1, "w.length"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 0)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 1000)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[length(bad.spct), "w.length"] <- 100)
  expect_equal(bad.spct, reference.spct)
  # negative reflectances are invalid, valid range is 0..1
  expect_error(bad.spct[1, "Rfr"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1:(nrow(bad.spct)/250 + 2), "Rfr"] <- -2e-4)
  expect_equal(bad.spct, reference.spct)
  expect_silent(bad.spct[1, "Rfr"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "Rfr"] <- 1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
})

context("object_spct")

test_that("check_spct works for object_spct", {
  # tests in test_data.r test for false positives implicitly
  # here we test for failing tests

  empty.spct <- object_spct()
  expect_silent(check_spct(empty.spct))

  bad.spct <- reference.spct <- Ler_leaf.spct
  expect_error(bad.spct[1, "w.length"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 0)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 1000)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[length(bad.spct), "w.length"] <- 100)
  expect_equal(bad.spct, reference.spct)
  # negative transmittances are invalid, valid range is 0..1
  expect_error(bad.spct[1, "Rfr"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_silent(bad.spct[1, "Rfr"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "Rfr"] <- 1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  # negative reflectances are invalid, valid range is 0..1
  bad.spct <- reference.spct <- Ler_leaf.spct
  expect_error(bad.spct[1, "Rfr"] <- -1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  bad.spct <- reference.spct <- Ler_leaf.spct
  expect_error(bad.spct[1:(nrow(bad.spct)/250 + 2), "Rfr"] <- -2e-4)
  expect_equal(bad.spct, reference.spct)
  expect_silent(bad.spct[1, "Rfr"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "Rfr"] <- 1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])

  expect_error(bad.spct[1, "w.length"] <- -1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_error(bad.spct[1, "w.length"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_error(bad.spct[1, "w.length"] <- 1000)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_error(bad.spct[length(bad.spct), "w.length"] <- 100)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  # negative transmittances are invalid, valid range is 0..1
  bad.spct <- reference.spct <- Ler_leaf.spct
  expect_error(bad.spct[1, "Tfr"] <- -1)
  expect_equal(bad.spct, reference.spct)
  bad.spct <- reference.spct <- Ler_leaf.spct
  expect_error(bad.spct[1:(nrow(bad.spct)/250 + 2), "Tfr"] <- -2e-4)
  expect_equal(bad.spct, reference.spct)
  expect_silent(bad.spct[1, "Tfr"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "Tfr"] <- 1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
})

context("calibration_spct")

test_that("check_spct works for reflector_spct", {
  # tests in test_data.r test for false positives implicitly
  # here we test for failing tests

  empty.spct <- calibration_spct()
  expect_silent(check_spct(empty.spct))

  bad.spct <- reference.spct <- calibration_spct(w.length = 300:400, irrad.mult = 1)
  expect_error(bad.spct[1, "w.length"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 0)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 1000)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[length(bad.spct), "w.length"] <- 100)
  expect_equal(bad.spct, reference.spct)
  # negative multipliers are invalid, valid range is 0..1
  expect_error(bad.spct[1, "irrad.mult"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_silent(bad.spct[1, "irrad.mult"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "irrad.mult"] <- 1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
})

context("solute_spct")

test_that("check_spct works for reflector_spct", {
  # tests in test_data.r test for false positives implicitly
  # here we test for failing tests

  empty.spct <- reflector_spct()
  expect_silent(check_spct(empty.spct))

  bad.spct <- reference.spct <- green_leaf.spct
  expect_error(bad.spct[1, "w.length"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 0)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 1000)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[length(bad.spct), "w.length"] <- 100)
  expect_equal(bad.spct, reference.spct)
  # negative reflectances are invalid, valid range is 0..1
  expect_error(bad.spct[1, "Rfr"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1:(nrow(bad.spct)/250 + 2), "Rfr"] <- -2e-4)
  expect_equal(bad.spct, reference.spct)
  expect_silent(bad.spct[1, "Rfr"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "Rfr"] <- 1)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
})

context("solute_spct")

test_that("check_spct works for reflector_spct", {
  # tests in test_data.r test for false positives implicitly
  # here we test for failing tests

  empty.spct <- solute_spct()
  expect_silent(check_spct(empty.spct))

  # bad.spct <- reference.spct <- green_leaf.spct
  # expect_error(bad.spct[1, "w.length"] <- -1)
  # expect_equal(bad.spct, reference.spct)
  # expect_error(bad.spct[1, "w.length"] <- 0)
  # expect_equal(bad.spct, reference.spct)
  # expect_error(bad.spct[1, "w.length"] <- 1000)
  # expect_equal(bad.spct, reference.spct)
  # expect_error(bad.spct[length(bad.spct), "w.length"] <- 100)
  # expect_equal(bad.spct, reference.spct)
  # # negative reflectances are invalid, valid range is 0..1
  # expect_error(bad.spct[1, "Rfr"] <- -1)
  # expect_equal(bad.spct, reference.spct)
  # expect_silent(bad.spct[1, "Rfr"] <- 0)
  # expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  # expect_silent(bad.spct[1, "Rfr"] <- 1)
  # expect_equal(bad.spct[-1, ], reference.spct[-1, ])
})


context("cps_spct")

test_that("check_spct works for cps_spct", {
  # tests in test_data.r test for false positives implicitly
  # here we test for failing tests

  empty.spct <- cps_spct()
  expect_silent(check_spct(empty.spct))

  bad.spct <- reference.spct <- white_led.cps_spct
  expect_error(bad.spct[1, "w.length"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 0)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 1000)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[length(bad.spct), "w.length"] <- 100)
  expect_equal(bad.spct, reference.spct)
  # negative reflectances are invalid, valid range is 0..1
  expect_error(bad.spct[1, "cps"] <- -300000)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1:(nrow(bad.spct)/100 + 2), "cps"] <- -30000)
  expect_equal(bad.spct, reference.spct)
  expect_silent(bad.spct[1, "cps"] <- 0)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
  expect_silent(bad.spct[1, "cps"] <- Inf)
  expect_equal(bad.spct[-1, ], reference.spct[-1, ])
})

context("raw_spct")

test_that("check_spct works for raw_spct", {
  # tests in test_data.r test for false positives implicitly
  # here we test for failing tests

  empty.spct <- raw_spct()
  expect_silent(check_spct(empty.spct))

  bad.spct <- reference.spct <- white_led.raw_spct
  expect_error(bad.spct[1, "w.length"] <- -1)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 0)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[1, "w.length"] <- 1000)
  expect_equal(bad.spct, reference.spct)
  expect_error(bad.spct[length(bad.spct), "w.length"] <- 100)
  expect_equal(bad.spct, reference.spct)
  # negative counts are possible during computations
  # no range check implemented
})
