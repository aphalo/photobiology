library("photobiology")

context("generic_spct")

test_that("constructor", {
  my.spct <- data.frame(w.length = 400:409, x = 1)
  setGenericSpct(my.spct)
  expect_equal(class(my.spct)[1], c("generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 1)

  expect_equal(names(my.spct), c("w.length", "x"))
  expect_is(my.spct, "generic_spct")
  expect_is(my.spct, "data.frame")
  expect_equal(min(my.spct), 400)
  expect_equal(max(my.spct), 409)
  expect_equal(range(my.spct), c(400, 409))
  expect_equal(spread(my.spct), 9)
  expect_equal(midpoint(my.spct), (400 + 409) / 2)
  my.b.spct <- data.frame(w.length = 0:101, x = 1)
  expect_error(setGenericSpct(my.b.spct))
  my.b.spct <- data.frame(w.length = 9999:10001, x = 1)
  expect_error(setGenericSpct(my.b.spct))
})
