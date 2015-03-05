library("photobiology")
context("normalize.spct")

test_that("normalize", {

  my.spct <- copy(sun.spct)

  expect_equal(max(normalize(my.spct)$s.e.irrad), 1)
  expect_equal(max(normalize(my.spct, norm = "max")$s.e.irrad), 1)
  expect_warning(irrad(normalize(my.spct, norm = "max")))
  expect_equal(irrad(normalize(my.spct, norm = "max")), NA)
  expect_named(normalize(my.spct), setdiff(names(my.spct), "s.q.irrad"))
  expect_equal(class(normalize(my.spct)), class(my.spct))
  expect_error(normalize(my.spct, range = 100))
  expect_error(normalize(my.spct, range = c(100, 100)))
  expect_true(is.normalized(normalize(my.spct)))
  expect_false(is.normalized(my.spct))
  expect_equal(is.source.spct(normalize(my.spct)), is.source.spct(my.spct))
  expect_equal(is.filter.spct(normalize(my.spct)), is.filter.spct(my.spct))
  expect_equal(getTimeUnit(normalize(my.spct)), getTimeUnit(my.spct))
  expect_equal(comment(normalize(my.spct)), comment(my.spct))
})

