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

context("rescale.spct")

test_that("Rescale", {

  my.spct <- q2e(sun.spct, action = "replace")

  expect_equivalent(integrate_spct(Rescale(my.spct, f = "total")), 1)
  expect_less_than(integrate_spct(Rescale(my.spct, f = "mean")) * average_spct(my.spct) - irrad(my.spct), 0.25)
  expect_warning(irrad(Rescale(my.spct, f = "mean")))
  expect_equal(irrad(Rescale(my.spct, f = "mean")), NA)
  expect_named(Rescale(my.spct), setdiff(names(my.spct), "s.q.irrad"))
  expect_equal(class(Rescale(my.spct)), class(my.spct))
  expect_error(Rescale(my.spct, range = 100))
  expect_error(normalize(my.spct, range = c(100, 100)))
  expect_true(is.rescaled(Rescale(my.spct)))
  expect_false(is.rescaled(my.spct))
  expect_equal(is.source.spct(Rescale(my.spct)), is.source.spct(my.spct))
  expect_equal(is.filter.spct(Rescale(my.spct)), is.filter.spct(my.spct))
  expect_equal(getTimeUnit(Rescale(my.spct)), getTimeUnit(my.spct))
  expect_equal(comment(Rescale(my.spct)), comment(my.spct))
})

context("integrate_spct average_spct")

test_that("integrate_spct", {

  my.spct <- source.spct(w.length=100:200, s.e.irrad = 1)

  expect_equivalent(integrate_spct(my.spct), sum(my.spct$s.e.irrad) - 1)
  expect_equivalent(average_spct(my.spct), 1)
  expect_equivalent(average_spct(my.spct * 2), 2)
  expect_named(average_spct(my.spct), "e.irrad")


  my.spct <- source.spct(w.length=seq(from=1000, to=2000, by=10), s.e.irrad = 1)

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
