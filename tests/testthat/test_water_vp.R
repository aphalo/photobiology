library("photobiology")
library("lubridate")

context("water_calc")

test_that("water_reference_values", {

  expect_equal(round(water_vp_sat(20), 3), 2338.023)
  expect_equal(round(water_vp_sat(20, method = "tetens"), 3), 2338.023)
  expect_equal(round(water_vp_sat(20, method = "Tetens"), 3), 2338.023)
  expect_equal(round(water_vp_sat(20, method = "magnus"), 3), 2333.441)
  expect_equal(round(water_vp_sat(20, method = "Magnus"), 3), 2333.441)
  expect_equal(round(water_vp_sat(20, method = "wexler"), 3), 2339.262)
  expect_equal(round(water_vp_sat(20, method = "Wexler"), 3), 2339.262)
  expect_equal(round(water_vp_sat(20, method = "goff.gratch"), 3), 2335.856)
  expect_equal(round(water_vp_sat(20, method = "Goff.Gratch"), 3), 2335.856)

  expect_equal(round(water_vp_sat(-30), 3), 50.179)
  expect_equal(round(water_vp_sat(-30, method = "tetens"), 3), 50.179)
  expect_equal(round(water_vp_sat(-30, method = "Tetens"), 3), 50.179)
  expect_equal(round(water_vp_sat(-30, method = "magnus"), 3), 51.064)
  expect_equal(round(water_vp_sat(-30, method = "Magnus"), 3), 51.064)
  expect_equal(round(water_vp_sat(-50, method = "magnus"), 3), 6.359)
  expect_equal(round(water_vp_sat(-50, method = "Magnus"), 3), 6.359)
  expect_equal(round(water_vp_sat(-50, method = "wexler"), 3), 6.438)
  expect_equal(round(water_vp_sat(-50, method = "Wexler"), 3), 6.438)
  expect_equal(round(water_vp_sat(-50, method = "goff.gratch"), 3), 6.349)
  expect_equal(round(water_vp_sat(-50, method = "Goff.Gratch"), 3), 6.349)
  expect_equal(round(water_vp_sat(-100, method = "wexler"), 7), 0.0036174)
  expect_equal(round(water_vp_sat(-100, method = "Wexler"), 7), 0.0036174)
  expect_equal(round(water_vp_sat(100, method = "wexler"), 0), 101418)
  expect_equal(round(water_vp_sat(100, method = "Wexler"), 0), 101418)

})

test_that("water_vp_sat", {
  #  test.path <- tempfile()
  test.path <- "./data/water-vapour-test-values"

  expect_known_value(
    water_vp_sat(-40:50),
    file = test.path
  )

  expect_warning(water_vp_sat(-40.001))
#  expect_is(water_vp_sat(-50.001), "numeric")
#  expect_true(is.na(water_vp_sat(-50.001)))

})

test_that("water_vp_sat_magnus", {
  #  test.path <- tempfile()
  test.path <- "./data/water-vapour-magnus-test-values"

  expect_known_value(
    water_vp_sat(-80:50, method = "magnus"),
    file = test.path
  )

  expect_warning(water_vp_sat(-80.001, method = "magnus"))
#  expect_is(water_vp_sat(-50.001, method = "magnus"), "numeric")
#  expect_true(is.na(water_vp_sat(-50.001, method = "magnus")))

})

test_that("water_vp_sat_wexler", {
  #  test.path <- tempfile()
  test.path <- "./data/water-vapour-wexler-test-values"

  expect_known_value(
    water_vp_sat(-100:100, method = "wexler"),
    file = test.path
  )

  expect_warning(water_vp_sat(-100.001, method = "wexler"))
  #  expect_is(water_vp_sat(-50.001, method = "magnus"), "numeric")
  #  expect_true(is.na(water_vp_sat(-50.001, method = "magnus")))

})

test_that("water_vp_sat_ice", {
  #  test.path <- tempfile()
  test.path <- "./data/water-vapour-ice-test-values"

  expect_known_value(
    signif(water_vp_sat(-30:0, over.ice = TRUE), 12),
    file = test.path
  )

  expect_warning(water_vp_sat(-30:10, over.ice = TRUE))

  expect_warning(water_vp_sat(1, over.ice = TRUE))

})

test_that("water_vp_sat_ice_magnus", {
  #  test.path <- tempfile()
  test.path <- "./data/water-vapour-ice-magnus-test-values"

  expect_known_value(
    signif(water_vp_sat(-80:0, over.ice = TRUE, method = "magnus"), 12),
    file = test.path
  )

  expect_warning(water_vp_sat(-30:10, over.ice = TRUE, method = "magnus"))

  expect_warning(water_vp_sat(1, over.ice = TRUE, method = "magnus"))

})

test_that("water_vp_sat_ice_wexler", {
  #  test.path <- tempfile()
  test.path <- "./data/water-vapour-ice-wexler-test-values"

  expect_known_value(
    signif(water_vp_sat(-100:0, over.ice = TRUE, method = "wexler"), 12),
    file = test.path
  )

  expect_warning(water_vp_sat(-30:10, over.ice = TRUE, method = "wexler"))

  expect_warning(water_vp_sat(1, over.ice = TRUE, method = "wexler"))

})

test_that("water_vp_sat_ice_goff.gratch", {
  #  test.path <- tempfile()
  test.path <- "./data/water-vapour-ice-goff.gratch-test-values"

  expect_known_value(
    signif(water_vp_sat(-50:0, over.ice = TRUE, method = "goff.gratch"), 12),
    file = test.path
  )

  expect_warning(water_vp_sat(-50:10, over.ice = TRUE, method = "goff.gratch"))

  expect_warning(water_vp_sat(1, over.ice = TRUE, method = "goff.gratch"))

})

test_that("water_vapour_dp", {
  #  test.path <- tempfile()
  test.path <- "./data/dew-point-test-values"

  expect_known_value(
    round(water_dp((6:110) * 100), 6),
    file = test.path
  )

})

test_that("water_vapour_dp_magnus", {
  #  test.path <- tempfile()
  test.path <- "./data/dew-point-magnus-test-values"

  expect_known_value(
    round(water_dp((6:110) * 100, method = "magnus"), 6),
    file = test.path
  )

})

test_that("water_vapour_dp_wexler", {
  #  test.path <- tempfile()
  test.path <- "./data/dew-point-wexler-test-values"

  expect_known_value(
    round(water_dp((6:900) * 100, method = "wexler"), 6),
    file = test.path
  )

})

test_that("water_vapour_dp_ice", {
  #  test.path <- tempfile()
  test.path <- "./data/dew-point-ice-test-values"

  expect_known_value(
    water_dp(40:600, over.ice = TRUE),
    file = test.path
  )

  expect_equal(
    water_dp((4:60) * 10, over.ice = TRUE),
    water_fp((4:60) * 10)
  )

  expect_equal(
    water_fp((61:500) * 10, over.ice = FALSE),
    water_dp((61:500) * 10)
  )

})

test_that("water_vapour_mvc", {
  #  test.path <- tempfile()
  test.path <- "./data/mvc-test-values"

  expect_known_value(
    water_vp2mvc((0:200) * 10, 20),
    file = test.path
  )

})

test_that("water_vapour_mvc_ice", {
  #  test.path <- tempfile()
  test.path <- "./data/mvc-ice-test-values"

  expect_known_value(
    water_dp(50:500, over.ice = TRUE),
    file = test.path
  )

})

test_that("water_vapour_RH", {
  #  test.path <- tempfile()
  test.path <- "./data/rh-test-values"

  expect_known_value(
    round(water_vp2RH((5:100) * 100, 50), 6),
    file = test.path
  )

  test.path <- "./data/rh-tetens-test-values"
  expect_known_value(
    round(water_vp2RH((5:100) * 100, 50, method = "tetens"), 6),
    file = test.path
  )

  test.path <- "./data/rh-magnus-test-values"
  expect_known_value(
    round(water_vp2RH((5:100) * 100, 50, method = "magnus"), 6),
    file = test.path
  )

  test.path <- "./data/rh-wexler-test-values"
  expect_known_value(
    round(water_vp2RH((5:100) * 100, 50, method = "wexler"), 6),
    file = test.path
  )

  test.path <- "./data/rh-goff.gratch-test-values"
  expect_known_value(
    round(water_vp2RH((5:100) * 100, 50, method = "goff.gratch"), 6),
    file = test.path
  )

})

test_that("water_consistency", {

  expect_equal(
    round(water_dp(water_vp_sat(-30:50)), 6),
    -30:50
  )

  expect_equal(
    round(water_dp(water_vp_sat(-80:49, method = "magnus"), method = "magnus"), 6),
    -80:49
  )

  # using ITS90 coefs
  expect_equal(
    round(water_dp(water_vp_sat(-99:40, method = "wexler"), method = "wexler"), 4),
    -99:40
  )

  expect_equal(
    round(water_dp(water_vp_sat(41:99, method = "wexler"), method = "wexler"), 3),
    41:99
  )

  expect_equal(
    round(water_mvc2vp(water_vp2mvc((0:200) * 10, 100), 100), 6),
    (0:200) * 10
  )

})


test_that("ET_ref", {

  expect_equal(
    round(ET_ref((1:10) * 10, 800, 5, 200), 5),
    c(0.34555, 0.34560, 0.34563, 0.34565, 0.34566, 0.34567, 0.34568, 0.34568, 0.34569, 0.34569)
   )

  expect_equal(
    round(ET_ref(50, (0:10) * 100, 5, 200), 5),
    c(0.34567, 0.34567, 0.34567, 0.34567, 0.34566, 0.34566, 0.34566, 0.34566, 0.34566, 0.34566, 0.34566)
  )

  expect_equal(
    round(ET_ref(20, 2000, 0:10, 200), 5),
    c(0.34558, 0.34558, 0.34558, 0.34557, 0.34557, 0.34557, 0.34556, 0.34556, 0.34556, 0.34555, 0.34555)
  )

  expect_equal(
    round(ET_ref(20, 2000, 5, (1:10) * 50), 5),
    c(0.08640, 0.17279, 0.25918, 0.34557, 0.43196, 0.51834, 0.60473, 0.69112, 0.77751, 0.86390)
  )

})
