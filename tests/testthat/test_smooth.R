library("photobiology")

context("smooth_spct")

test_that("smooth_default", {
  #  test.path <- tempfile()
  test.path <- "smooth-default-test-value"

  testthat::expect_known_value(
    smooth_spct(green_leaf.spct),
    file = test.path
  )
}
)

test_that("smooth_lowess", {
  #  test.path <- tempfile()
  test.path <- "smooth-lowess-test-value"

  testthat::expect_known_value(
    smooth_spct(green_leaf.spct,
                method = "lowess"),
    file = test.path
  )
}
)

test_that("smooth_supsmu", {
  #  test.path <- tempfile()
  test.path <- "smooth-supsmu-test-value"

  testthat::expect_known_value(
    smooth_spct(green_leaf.spct,
                method = "supsmu"),
    file = test.path
  )
}
)

x <- sun.spct
x$s.e.irrad[1:100] <- NA_real_
expect_warning(smooth_spct(x))
