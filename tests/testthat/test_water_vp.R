library("photobiology")
library("lubridate")

context("water_calc")

test_that("water_reference_values", {

  expect_equal(round(water_vp_sat(20), 3), 2338.023)
  expect_equal(round(water_vp_sat(20, method = "tetens"), 3), 2338.023)
  expect_equal(round(water_vp_sat(20, method = "Tetens"), 3), 2338.023)
  expect_equal(round(water_vp_sat(20, method = "magnus"), 3), 2333.441)
  expect_equal(round(water_vp_sat(20, method = "Magnus"), 3), 2333.441)

  expect_equal(round(water_vp_sat(-50), 3), 6.078)
  expect_equal(round(water_vp_sat(-50, method = "tetens"), 3), 6.078)
  expect_equal(round(water_vp_sat(-50, method = "Tetens"), 3), 6.078)
  expect_equal(round(water_vp_sat(-50, method = "magnus"), 3), 6.359)
  expect_equal(round(water_vp_sat(-50, method = "Magnus"), 3), 6.359)

})

test_that("water_vp_sat", {
  #  test.path <- tempfile()
  test.path <- "water-vapour-test-values"

  expect_known_value(
    water_vp_sat(-50:100),
    file = test.path
  )

  expect_warning(water_vp_sat(-50.001))
#  expect_is(water_vp_sat(-50.001), "numeric")
#  expect_true(is.na(water_vp_sat(-50.001)))

})

test_that("water_vp_sat_magnus", {
  #  test.path <- tempfile()
  test.path <- "water-vapour-magnus-test-values"

  expect_known_value(
    water_vp_sat(-50:100, method = "magnus"),
    file = test.path
  )

  expect_warning(water_vp_sat(-50.001, method = "magnus"))
#  expect_is(water_vp_sat(-50.001, method = "magnus"), "numeric")
#  expect_true(is.na(water_vp_sat(-50.001, method = "magnus")))

})

test_that("water_vp_sat_ice", {
  #  test.path <- tempfile()
  test.path <- "water-vapour-ice-test-values"

  expect_known_value(
    signif(water_vp_sat(-30:0, over.ice = TRUE), 12),
    file = test.path
  )

  expect_warning(water_vp_sat(-30:10, over.ice = TRUE))

  expect_warning(water_vp_sat(1, over.ice = TRUE))

})

test_that("water_vp_sat_ice_magnus", {
  #  test.path <- tempfile()
  test.path <- "water-vapour-ice-magnus-test-values"

  expect_known_value(
    signif(water_vp_sat(-30:0, over.ice = TRUE, method = "magnus"), 12),
    file = test.path
  )

  expect_warning(water_vp_sat(-30:10, over.ice = TRUE, method = "magnus"))

  expect_warning(water_vp_sat(1, over.ice = TRUE, method = "magnus"))

})

test_that("water_vapour_dp", {
  #  test.path <- tempfile()
  test.path <- "dew-point-test-values"

  expect_known_value(
    round(water_dp((20:100) * 10), 6),
    file = test.path
  )

})

test_that("water_vapour_dp_ice", {
  #  test.path <- tempfile()
  test.path <- "dew-point-ice-test-values"

  expect_known_value(
    water_dp(50:500, over.ice = TRUE),
    file = test.path
  )

})

test_that("water_vapour_mvc", {
  #  test.path <- tempfile()
  test.path <- "mvc-test-values"

  expect_known_value(
    water_vp2mvc((0:200) * 10, 20),
    file = test.path
  )

})

test_that("water_vapour_mvc_ice", {
  #  test.path <- tempfile()
  test.path <- "mvc-ice-test-values"

  expect_known_value(
    water_dp(50:500, over.ice = TRUE),
    file = test.path
  )

})

test_that("water_consistency", {

  expect_equal(
    round(water_dp(water_vp_sat(-50:100)), 6),
    -50:100
  )

  expect_equal(
    round(water_mvc2vp(water_vp2mvc((0:200) * 10, 100), 100), 6),
    (0:200) * 10
  )

})
