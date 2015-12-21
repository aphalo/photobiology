library("photobiology")

context("interpolate_wl")

test_that("source_spct", {

  my.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)

  expect_equal(interpolate_wl(my.spct), my.spct)
  expect_equal(interpolate_wl(my.spct, w.length.out = NULL, length.out = NULL), my.spct)

  interpolated.spct <- interpolate_wl(my.spct, length.out = 21)

  expect_equal(nrow(interpolated.spct), 21)
  expect_equal(min(interpolated.spct), 400)
  expect_equal(max(interpolated.spct), 410)
  expect_named(interpolated.spct, names(my.spct))
  expect_equal(max(interpolated.spct[["s.e.irrad"]]), 1)
  expect_equal(min(interpolated.spct[["s.e.irrad"]]), 1)

  my1.spct <- interpolate_wl(my.spct, length.out = 0)
  expect_equal(nrow(my1.spct), 0)
  expect_named(my1.spct, names(my.spct))

  my1.spct <- interpolate_wl(my.spct, length.out = 1)
  expect_equal(nrow(my1.spct), 1)
  expect_named(my1.spct, names(my.spct))

  my1.spct <- interpolate_wl(my.spct, length.out = 2)
  expect_equal(nrow(my1.spct), 2)
  expect_named(my1.spct, names(my.spct))

  my1.spct <- interpolate_wl(my.spct, length.out = 10)
  expect_equal(nrow(my1.spct), 10)
  expect_named(my1.spct, names(my.spct))

  my1.spct <- interpolate_wl(my.spct, length.out = NA)
  expect_equal(nrow(my1.spct), 0)
  expect_named(my1.spct, names(my.spct))

  my1.spct <- interpolate_wl(my.spct, length.out = numeric())
  expect_equal(nrow(my1.spct), 0)
  expect_named(my1.spct, names(my.spct))

  expect_warning(interpolated.spct <-
                   interpolate_wl(my.spct, w.length.out = c(398:409, NA, 250), fill = 0))

  expect_equal(nrow(interpolated.spct), length(c(398:409, 250)))
  expect_equal(min(interpolated.spct), 250)
  expect_equal(max(interpolated.spct), 409)
  expect_equal(min(interpolated.spct$s.e.irrad), 0)
  expect_equal(max(interpolated.spct$s.e.irrad), 1)
  expect_named(interpolated.spct, names(my.spct))

  expect_warning(interpolated.spct <-
                   interpolate_wl(my.spct, w.length.out = c(398:409, NA, 250), fill = NA))

  expect_equal(nrow(interpolated.spct), length(c(398:409, 250)))
  expect_equal(min(interpolated.spct), 250)
  expect_equal(max(interpolated.spct), 409)
  expect_equal(min(interpolated.spct$s.e.irrad, na.rm = TRUE), 1)
  expect_equal(max(interpolated.spct$s.e.irrad, na.rm = TRUE), 1)
  expect_named(interpolated.spct, names(my.spct))

  expect_warning(interpolated.spct <-
                   interpolate_wl(my.spct, w.length.out = c(398:409, NA, 250), fill = NULL))

  expect_equal(nrow(interpolated.spct), 10) # because of trimming
  expect_equal(min(interpolated.spct), 400)
  expect_equal(max(interpolated.spct), 409)
  expect_equal(min(interpolated.spct$s.e.irrad), 1)
  expect_equal(max(interpolated.spct$s.e.irrad), 1)
  expect_named(interpolated.spct, names(my.spct))

})

