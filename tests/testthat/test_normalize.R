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
  expect_true(is_normalized(normalize(my.spct)))
  expect_false(is_normalized(my.spct))
  expect_equal(is.source_spct(normalize(my.spct)), is.source_spct(my.spct))
  expect_equal(is.filter_spct(normalize(my.spct)), is.filter_spct(my.spct))
  expect_equal(getTimeUnit(normalize(my.spct)), getTimeUnit(my.spct))
  expect_equal(comment(normalize(my.spct)), comment(my.spct))
})

context("scale.spct")

test_that("scale", {

  my.spct <- q2e(sun.spct, action = "replace")

  expect_equivalent(integrate_spct(scale(my.spct, f = "total")), 1)
  expect_less_than(integrate_spct(scale(my.spct, f = "mean")) * average_spct(my.spct) - irrad(my.spct), 0.25)
  expect_warning(irrad(scale(my.spct, f = "mean")))
  expect_equal(irrad(scale(my.spct, f = "mean")), NA)
  expect_named(scale(my.spct), setdiff(names(my.spct), "s.q.irrad"))
  expect_equal(class(scale(my.spct)), class(my.spct))
  expect_error(scale(my.spct, range = 100))
  expect_error(normalize(my.spct, range = c(100, 100)))
  expect_true(is_scaled(scale(my.spct)))
  expect_false(is_scaled(my.spct))
  expect_equal(is.source_spct(scale(my.spct)), is.source_spct(my.spct))
  expect_equal(is.filter_spct(scale(my.spct)), is.filter_spct(my.spct))
  expect_equal(getTimeUnit(scale(my.spct)), getTimeUnit(my.spct))
  expect_equal(comment(scale(my.spct)), comment(my.spct))
})

context("integrate_spct average_spct")

test_that("integrate_spct", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)

  expect_equivalent(integrate_spct(my.spct), sum(my.spct$s.e.irrad) - 1)
  expect_equivalent(average_spct(my.spct), 1)
  expect_equivalent(average_spct(my.spct * 2), 2)
  expect_named(average_spct(my.spct), "e.irrad")


  my.spct <- source_spct(w.length=seq(from=1000, to=2000, by=10), s.e.irrad = 1)

  expect_equivalent(integrate_spct(my.spct), (sum(my.spct$s.e.irrad) - 1) * 10)
  expect_equivalent(average_spct(my.spct), 1)
  expect_equivalent(average_spct(my.spct * 2), 2)

  e2q(my.spct, byref = TRUE)

  expect_equivalent(average_spct(my.spct), c(1, 1.2538837047156523583e-05))
  expect_named(average_spct(my.spct), c("e.irrad", "q.irrad"))

  e2q(my.spct, action="replace", byref = TRUE)

  expect_equivalent(average_spct(my.spct), 1.2538837047156523583e-05)
  expect_named(average_spct(my.spct), "q.irrad")


})
