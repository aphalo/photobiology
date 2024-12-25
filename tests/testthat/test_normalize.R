library("photobiology")

context("setNormalized")

test_that("setNormalized source_spct", {
  energy_as_default()

  my.spct <- sun.spct
  my.spct <- q2e(my.spct, action = "replace")
  my_norm.spct <- normalize(my.spct)
  my_unset_norm.spct <- my_norm.spct
  setNormalized(my_unset_norm.spct, norm = FALSE)
  my_defaultset_norm.spct <- my_norm.spct
  setNormalized(my_defaultset_norm.spct)

  expect_false(is_normalized(my.spct))
  expect_true(is_normalized(my_norm.spct))
  expect_false(is_normalized(my_unset_norm.spct))
  expect_false(is_normalized(my_defaultset_norm.spct))

  # check query function and method consistency
  expect_equal(normalization(my_norm.spct),
               normalization(my_norm.spct))

  expect_true(all(is.na(unlist(normalization(my.spct)))))
  expect_true(!any(is.na(unlist(normalization(my_norm.spct)))))
  expect_true(all(is.na(unlist(normalization(my_unset_norm.spct)))))
  expect_true(all(is.na(unlist(normalization(my_defaultset_norm.spct)))))

  expect_named(normalization(my.spct),
               c("norm.type", "norm.wl", "norm.factors", "norm.cols", "norm.range"))
  expect_named(normalization(my_norm.spct),
               c("norm.type", "norm.wl", "norm.factors", "norm.cols", "norm.range"))
  expect_named(normalization(my_unset_norm.spct),
               c("norm.type", "norm.wl", "norm.factors", "norm.cols", "norm.range"))
  expect_named(normalization(my_defaultset_norm.spct),
               c("norm.type", "norm.wl", "norm.factors", "norm.cols", "norm.range"))

  my.spct <- sun.spct
  expect_silent(setNormalized(as.data.frame(my.spct)))
  my.spct <- sun.spct
  expect_silent(setNormalized(as.data.frame(my.spct), FALSE))
  my.spct <- sun.spct
  expect_warning(setNormalized(as.data.frame(my.spct), TRUE))

  temp.spct <- sun.spct
  setNormalized(temp.spct, norm = 450)
  expect_true(is_normalized(temp.spct))
  expect_equal(getNormalized(temp.spct), 450) # backwards compatibility
  expect_equal(sum(is.na(unlist(normalization(temp.spct)))), 5L)
  expect_named(normalization(temp.spct),
               c("norm.type", "norm.wl", "norm.factors", "norm.cols", "norm.range"))
  expect_equal(normalization(temp.spct)[["norm.wl"]], 450)

  temp.spct <- sun.spct
  setNormalized(temp.spct, norm = TRUE, norm.type = "max")
  expect_true(is_normalized(temp.spct))
  expect_equal(getNormalized(temp.spct), TRUE) # backwards compatibility
  expect_equal(sum(is.na(unlist(normalization(temp.spct)))), 5L)
  expect_named(normalization(temp.spct),
               c("norm.type", "norm.wl", "norm.factors", "norm.cols", "norm.range"))
  expect_equal(normalization(temp.spct)[["norm.type"]], "max")
  expect_true(is.na(normalization(temp.spct)[["norm.wl"]]))

  temp.spct <- sun.spct
  expect_silent(setNormalized(temp.spct, norm = TRUE, norm.cols = c("s.e.irrad", "s.q.irrad")))
  temp.spct <- sun.spct
  expect_silent(setNormalized(temp.spct, norm = TRUE, norm.cols = "s.q.irrad"))
  temp.spct <- sun.spct
  expect_silent(setNormalized(temp.spct, norm = TRUE, norm.cols = "s.e.irrad"))
  temp.spct <- sun.spct
  expect_error(setNormalized(temp.spct, norm = TRUE, norm.cols = "zzzzzz"))
  temp.spct <- sun.spct
  expect_silent(setNormalized(temp.spct, norm = TRUE, norm.cols = c("s.q.irrad", "zzzzzz")))
  temp.spct <- sun.spct
  expect_silent(setNormalized(temp.spct, norm = TRUE, norm.cols = c("s.e.irrad", "zzzzzz")))
  temp.spct <- sun.spct
  expect_silent(setNormalized(temp.spct, norm = TRUE, norm.cols = c("s.q.irrad", NA)))
  temp.spct <- sun.spct
  expect_error(setNormalized(temp.spct, norm = TRUE, norm.cols = ""))
  temp.spct <- sun.spct
  expect_silent(setNormalized(temp.spct, norm = TRUE, norm.cols = character()))

  temp.spct <- sun.spct
  setNormalized(temp.spct, norm = TRUE)
  expect_true(is_normalized(temp.spct))
  expect_equal(getNormalized(temp.spct), TRUE) # backwards compatibility
  expect_equal(sum(is.na(unlist(normalization(temp.spct)))), 6L)
  expect_named(normalization(temp.spct),
               c("norm.type", "norm.wl", "norm.factors", "norm.cols", "norm.range"))

})

context("normalize.spct")

test_that("normalize source_spct", {
  energy_as_default()

  my.spct <- q2e(sun.spct, action = "replace")
  my_norm.spct <- normalize(my.spct)

  # expect_equal(my.spct, setNormalized(my_norm.spct))
  # expect_equal(my.spct, setNormalized(my_norm.spct, norm = FALSE))

  # check query function and method consistency
  expect_equal(normalization(my_norm.spct),
               getNormalisation(my_norm.spct))

  # check norm = "skip" is a no-op
  expect_equal(my.spct, normalize(my.spct, norm = "skip"))

  # check norm = "undo" reverts the normalizatiom
  expect_equal(my.spct[["w.length"]],
               normalize(my_norm.spct, norm = "undo")[["w.length"]])
  expect_equal(my.spct[["s.e.irrad"]],
               normalize(my_norm.spct, norm = "undo")[["s.e.irrad"]])
  expect_contains(attributes(normalize(my_norm.spct, norm = "undo")),
                  attributes(my.spct))

  # check that default is norm = "max"
  my_norm_max.spct <- normalize(my.spct, norm = "max")
  expect_equal(getNormalization(my_norm.spct),
               getNormalization(my_norm_max.spct))

  # check that unnecessary update does not change the result
  my_norm_umax.spct <- normalize(my_norm_max.spct, norm = "update")
  expect_equal(getNormalization(my_norm_umax.spct),
               getNormalization(my_norm_max.spct))

  # check that old style normalization is handled correctly
  my_norm.spct <- normalize(my.spct)

  my_old_style_norm.spct <- my_norm.spct
  attr(my_old_style_norm.spct, "normalization") <- NULL
  expect_true(is_normalised(my_old_style_norm.spct))
  expect_no_warning(normalization(my_old_style_norm.spct))
  expect_no_error(normalization(my_old_style_norm.spct))
  expect_no_message(normalization(my_old_style_norm.spct))
  expect_false(all(is.na(unlist(normalization(my_old_style_norm.spct)))))
  expect_equal(normalization(my_old_style_norm.spct)[["norm.wl"]], 451)
  expect_warning(normalize(my_old_style_norm.spct, norm = "update"))

  # check that old style normalization is handled correctly
  my_vold_style_norm.spct <- my_old_style_norm.spct
  attr(my_vold_style_norm.spct, "normalized") <- TRUE
  expect_true(is_normalised(my_vold_style_norm.spct))
  expect_no_warning(normalization(my_vold_style_norm.spct))
  expect_no_error(normalization(my_vold_style_norm.spct))
  expect_no_message(normalization(my_vold_style_norm.spct))
  expect_true(all(is.na(unlist(normalization(my_vold_style_norm.spct)))))
  expect_true(is.na(normalization(my_vold_style_norm.spct)[["norm.wl"]]))
  expect_warning(normalize(my_vold_style_norm.spct, norm = "update"))

  my_norm_500.spct <- normalize(my.spct, norm = 500)
  my_norm_u500.spct <- normalize(my_norm_500.spct, norm = "update")
  expect_equal(getNormalization(my_norm_500.spct),
               getNormalization(my_norm_u500.spct))

  my.q.spct <- e2q(sun.spct, action = "replace")
  my_norm_qmax.spct <- normalize(my.q.spct, norm = "max")
  my_norm_qumax.spct <- normalize(my_norm_qmax.spct, norm = "update")
  expect_equal(getNormalization(my_norm_qmax.spct),
               getNormalization(my_norm_qumax.spct))

  my_norm_q500.spct <- normalize(my.q.spct, norm = 500)
  my_norm_qu500.spct <- normalize(my_norm_q500.spct, norm = "update")
  expect_equal(getNormalization(my_norm_q500.spct),
               getNormalization(my_norm_qu500.spct))

  my.eq.spct <- e2q(sun.spct, action = "add")
  my_norm_eqmax.spct <- normalize(my.eq.spct, norm = "max")
  my_norm_equmax.spct <- normalize(my_norm_eqmax.spct, norm = "update")
  expect_equal(getNormalization(my_norm_eqmax.spct),
               getNormalization(my_norm_equmax.spct))

  my_norm_eq500.spct <- normalize(my.eq.spct, norm = 500)
  my_norm_equ500.spct <- normalize(my_norm_eq500.spct, norm = "update")
  expect_equal(getNormalization(my_norm_eq500.spct),
               getNormalization(my_norm_equ500.spct))

  my_scaled.spct <- fscale(my.spct)
  my_norm_smax.spct <- normalize(my_scaled.spct, norm = "max")
  my_norm_ssmax.spct <- normalize(my_scaled.spct, norm = "max", keep.scaling = TRUE)

  expect_equal(normalize(my.spct, norm = "skip"), my.spct)
  expect_equal(normalize(my.spct, norm = "update"), my.spct)
  # expect_equal(normalize(my_norm_qmax.spct, norm = "update"),
  #              my_norm_max.spct)
  # expect_equal(normalize(my_norm_emax.spct, norm = "update"),
  #              normalize(my_norm_qmax.spct, norm = "update"))
  expect_equal(normalize(my_norm_qmax.spct, norm = "update"),
               normalize(my_norm_qmax.spct, norm = "max"))
  expect_equal(getNormalized(my_norm_max.spct),
               getNormalization(my_norm_max.spct)[["norm.wl"]])
  expect_equal(getNormalized(my_norm_qmax.spct),
               getNormalization(my_norm_qmax.spct)[["norm.wl"]])
  expect_equal(getNormalized(my_norm_qmax.spct),
               getNormalized(my_norm_qumax.spct))
  expect_true(all(is.na(unlist(getNormalization(my.spct)))))
  expect_false(all(is.na(unlist(getNormalization(my_norm_max.spct)))))
  expect_equal(getNormalization(my_norm_max.spct)[["norm.type"]], "max")
  expect_equal(getNormalization(my_norm_max.spct)[["norm.wl"]], 451, tolerance = 1e-4)
  expect_equal(getNormalization(my_norm_max.spct)[["norm.factors"]], 1.218823, tolerance = 1e-5)
  expect_equal(getNormalization(my_norm_max.spct)[["norm.cols"]], "s.e.irrad")
  expect_equal(max(my_norm_max.spct$s.e.irrad), 1, tolerance = 1e-5)
  expect_equal(max(my_norm_max.spct$s.e.irrad), 1, tolerance = 1e-5)
  expect_warning(irrad(normalize(my.spct, norm = "max")))
  expect_equal(suppressWarnings(irrad(normalize(my.spct, norm = "max"))), NA_real_)
  expect_named(my_norm_max.spct, names(my.spct))
  expect_equal(class(my_norm_max.spct), class(my.spct))
  expect_error(normalize(my.spct, range = 100))
  expect_error(normalize(my.spct, range = c(100, 100)))
  expect_true(is_normalized(my_norm_max.spct))
  expect_false(is_normalized(my.spct))
  expect_equal(is.source_spct(my_norm_max.spct), is.source_spct(my.spct))
  expect_equal(is.response_spct(my_norm_max.spct), is.response_spct(my.spct))
  expect_equal(is.filter_spct(my_norm_max.spct), is.filter_spct(my.spct))
  expect_equal(getTimeUnit(my_norm_max.spct), getTimeUnit(my.spct))
  expect_equal(comment(my_norm_max.spct), comment(my.spct))
  expect_equal(getNormalized(normalize(my.spct, norm = 400)), 400)
  expect_equal(getNormalized(normalize(my.spct, norm = 400.2)), 400.2)
  expect_equal(getNormalized(my_norm_max.spct), 451)
})

test_that("normalize response_spct", {

  my.spct <- q2e(ccd.spct, action = "replace")
  my_norm_max.spct <- normalize(my.spct, norm = "max")
  getNormalization(my_norm_max.spct)
  my_norm_500.spct <- normalize(my.spct, norm = 500)
  getNormalization(my_norm_500.spct)
  my_norm_umax.spct <- normalize(my_norm_max.spct, norm = "update")
  getNormalization(my_norm_umax.spct)

  my.q.spct <- e2q(ccd.spct, action = "replace")
  my_norm_qmax.spct <- normalize(my.q.spct, norm = "max")
  getNormalization(my_norm_qmax.spct)
  my_norm_q500.spct <- normalize(my.q.spct, norm = 500)
  getNormalization(my_norm_q500.spct)
  my_norm_qumax.spct <- normalize(my_norm_qmax.spct, norm = "update")
  getNormalization(my_norm_qumax.spct)

  my.eq.spct <- e2q(ccd.spct, action = "add")
  my_norm_eqmax.spct <- normalize(my.eq.spct, norm = "max")
  getNormalization(my_norm_eqmax.spct)
  my_norm_eq500.spct <- normalize(my.eq.spct, norm = 500)
  getNormalization(my_norm_eq500.spct)
  my_norm_equmax.spct <- normalize(my_norm_eqmax.spct, norm = "update")
  getNormalization(my_norm_equmax.spct)

  expect_equal(normalize(my.spct, norm = "skip"), my.spct)
  expect_equal(normalize(my.spct, norm = "update"), my.spct)
  # expect_equal(normalize(my_norm_qmax.spct, norm = "update"),
  #              my_norm_max.spct)
  # expect_equal(normalize(my_norm_emax.spct, norm = "update"),
  #              normalize(my_norm_qmax.spct, norm = "update"))
  expect_equal(normalize(my_norm_qmax.spct, norm = "update"),
               normalize(my_norm_qmax.spct, norm = "max"))
  expect_equal(getNormalized(my_norm_max.spct),
               getNormalization(my_norm_max.spct)[["norm.wl"]])
  expect_equal(getNormalized(my_norm_qmax.spct),
               getNormalization(my_norm_qmax.spct)[["norm.wl"]])
  expect_equal(getNormalized(my_norm_qmax.spct),
               getNormalized(my_norm_qumax.spct))
  expect_true(all(is.na(unlist(getNormalization(my.spct)))))
  expect_false(all(is.na(unlist(getNormalization(my_norm_max.spct)))))
  expect_equal(getNormalization(my_norm_max.spct)[["norm.type"]], "max")
  expect_equal(getNormalization(my_norm_max.spct)[["norm.wl"]], 742.6704, tolerance = 1e-4)
  expect_equal(getNormalization(my_norm_max.spct)[["norm.factors"]], 228802.2, tolerance = 1e-5)
  expect_equal(getNormalization(my_norm_max.spct)[["norm.cols"]], "s.e.response")
  expect_equal(max(normalize(my.spct)$s.e.response), 1, tolerance = 1e-5)
  expect_equal(max(my_norm_max.spct$s.e.response), 1, tolerance = 1e-5)
  expect_warning(response(normalize(my.spct, norm = "max")))
  expect_equal(suppressWarnings(response(normalize(my.spct, norm = "max"))), NA_real_)
  expect_named(my_norm_max.spct, names(my.spct))
  expect_equal(class(my_norm_max.spct), class(my.spct))
  expect_error(normalize(my.spct, range = 100))
  expect_error(normalize(my.spct, range = c(100, 100)))
  expect_true(is_normalized(my_norm_max.spct))
  expect_false(is_normalized(my.spct))
  expect_equal(is.source_spct(my_norm_max.spct), is.source_spct(my.spct))
  expect_equal(is.response_spct(my_norm_max.spct), is.response_spct(my.spct))
  expect_equal(is.filter_spct(my_norm_max.spct), is.filter_spct(my.spct))
  expect_equal(getTimeUnit(my_norm_max.spct), getTimeUnit(my.spct))
  expect_equal(comment(my_norm_max.spct), comment(my.spct))
  expect_equal(getNormalized(normalize(my.spct, norm = 400)), 400)
  expect_equal(getNormalized(normalize(my.spct, norm = 400.2)), 400.2)
  expect_equal(getNormalized(my_norm_max.spct), 742.6704, tolerance = 1e-4)
})

context("fscale.spct")

test_that("fscale", {

  my.spct <- q2e(sun.spct, action = "replace")

  expect_lt(abs(integrate_spct(fscale(my.spct, f = "total")) - 1), 1e6)
  expect_lt(abs(average_spct(fscale(my.spct, f = "mean")) - 1), 1e6)
  expect_warning(irrad(fscale(my.spct, f = "mean")))
  expect_false(suppressWarnings(is.na(irrad(fscale(my.spct, f = "mean")))))
  expect_named(fscale(my.spct), names(my.spct))
  expect_equal(class(fscale(my.spct)), class(my.spct))
  expect_error(fscale(my.spct, range = 281))
  expect_true(is_scaled(fscale(my.spct)))
  expect_false(is_scaled(my.spct))
  expect_equal(is.source_spct(fscale(my.spct)), is.source_spct(my.spct))
  expect_equal(is.filter_spct(fscale(my.spct)), is.filter_spct(my.spct))
  expect_equal(getTimeUnit(fscale(my.spct)), getTimeUnit(my.spct))
  expect_equal(comment(fscale(my.spct)), comment(my.spct))
  expect_true(is_scaled(fscale(my.spct)))
  expect_false(is_scaled(my.spct))
})

context("fshift.spct")

test_that("fshift", {

  my.spct <- trim_spct(sun.spct, range = c(200, 700), fill = 0)
  my.spct <- q2e(my.spct, action = "replace")

  expect_equal(suppressWarnings(irrad(fshift(my.spct - 1), f = "min")),
               irrad(my.spct))
  expect_equal(irrad(fshift(my.spct + 1, f = "min")), irrad(my.spct))
  expect_named(fshift(my.spct), names(my.spct))
  expect_equal(class(fshift(my.spct)), class(my.spct))
  expect_equal(range(fshift(my.spct)), range(my.spct))
  expect_equal(is_scaled(fshift(my.spct)), is_scaled(my.spct))
  expect_equal(is.source_spct(fshift(my.spct)), is.source_spct(my.spct))
  expect_equal(is.filter_spct(fshift(my.spct)), is.filter_spct(my.spct))
  expect_equal(getTimeUnit(fshift(my.spct)), getTimeUnit(my.spct))
  expect_equal(comment(fshift(my.spct)), comment(my.spct))
})

context("clean.spct")

test_that("clean", {

  my.spct <- q2e(sun.spct, action = "replace")
  expect_warning(my.spct[1, "s.e.irrad"] <- -1)

  expect_equal(irrad(clean(my.spct)), irrad(sun.spct))
  expect_equal(clean(my.spct), q2e(sun.spct, action = "replace"))
  expect_named(clean(my.spct), names(my.spct))
  expect_equal(class(clean(my.spct)), class(my.spct))
  expect_equal(range(clean(my.spct)), range(my.spct))
  expect_equal(is_scaled(clean(my.spct)), is_scaled(my.spct))
  expect_equal(is.source_spct(clean(my.spct)), is.source_spct(my.spct))
  expect_equal(is.filter_spct(clean(my.spct)), is.filter_spct(my.spct))
  expect_equal(getTimeUnit(clean(my.spct)), getTimeUnit(my.spct))
  expect_equal(comment(clean(my.spct)), comment(my.spct))
})

context("integrate_spct average_spct")

test_that("integrate_spct", {

  my.spct <- source_spct(w.length = 100:200, s.e.irrad = 1)

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

