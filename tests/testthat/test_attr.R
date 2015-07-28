library("photobiology")

context("conversions")

test_that("object_spct", {

  my.spct <- object_spct(w.length = 400:450, Tfr = 0.5, Rfr = 0.5)
  setTfrType(my.spct, "internal")
  setRfrType(my.spct, "specular")

  expect_equal(getTfrType(my.spct), "internal")
  expect_equal(getTfrType(as.filter_spct(my.spct)), "internal")

  expect_equal(getRfrType(my.spct), "specular")
  expect_equal(getRfrType(as.reflector_spct(my.spct)), "specular")

  expect_equal(getTfrType(as.object_spct(my.spct)), "internal")
  expect_equal(getRfrType(as.object_spct(my.spct)), "specular")

  expect_equal(getTfrType(as.object_spct(as.filter_spct(my.spct))), "internal")
  expect_equal(getRfrType(as.object_spct(as.reflector_spct(my.spct))), "specular")

})

context("fscale_spct")

test_that("fscale", {

  my.spct <- q2e(sun.spct, action = "replace")

  expect_equivalent(integrate_spct(fscale(my.spct, f = "total")), 1)
  expect_less_than(average_spct(fscale(my.spct, f = "mean")), 1)
  expect_warning(irrad(fscale(my.spct, f = "mean")))
  expect_equal(irrad(fscale(my.spct, f = "mean")), NA)
  expect_named(fscale(my.spct), setdiff(names(my.spct), "s.q.irrad"))
  expect_equal(class(fscale(my.spct)), class(my.spct))
  expect_error(fscale(my.spct, range = 100))
  expect_error(normalize(my.spct, range = c(100, 100)))
  expect_true(is_scaled(fscale(my.spct)))
  expect_false(is_scaled(my.spct))
  expect_equal(is.source_spct(fscale(my.spct)), is.source_spct(my.spct))
  expect_equal(is.filter_spct(fscale(my.spct)), is.filter_spct(my.spct))
  expect_equal(getTimeUnit(fscale(my.spct)), getTimeUnit(my.spct))
  expect_equal(comment(fscale(my.spct)), comment(my.spct))
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
