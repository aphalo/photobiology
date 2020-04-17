context("conversions")

test_that("T2A", {
  f.spct <- filter_spct(w.length = 300:320, Tfr = 0.1)
  expect_silent(T2A(f.spct))
  expect_silent(T2A(f.spct, clean = TRUE))
  Tfr_as_default()
  expect_silent(T2A(f.spct))
  expect_silent(T2A(f.spct, clean = TRUE))
  Afr_as_default()
  expect_silent(T2A(f.spct))
  expect_silent(T2A(f.spct, clean = TRUE))
  A_as_default()
  expect_silent(T2A(f.spct))
  expect_silent(T2A(f.spct, clean = TRUE))
  unset_user_defaults()

  expect_true(all(T2A(f.spct)[["A"]] == 1))

  # add "ignore.order = TRUE" if needed!
  expect_named(T2A(f.spct), c("w.length", "Tfr", "A"))
  expect_named(T2A(f.spct, action = "add"), c("w.length", "Tfr", "A"))
  expect_named(T2A(f.spct, action = "replace"), c("w.length", "A"))

  })

test_that("A2T", {
  f.spct <- filter_spct(w.length = 300:320, A = 1)
  expect_silent(A2T(f.spct))
  expect_silent(A2T(f.spct, clean = TRUE))
  Tfr_as_default()
  expect_silent(A2T(f.spct))
  expect_silent(A2T(f.spct, clean = TRUE))
  Afr_as_default()
  expect_silent(A2T(f.spct))
  expect_silent(A2T(f.spct, clean = TRUE))
  A_as_default()
  expect_silent(A2T(f.spct))
  expect_silent(A2T(f.spct, clean = TRUE))
  unset_user_defaults()

  expect_true(all(A2T(f.spct)[["Tfr"]] == 0.1))

  # add "ignore.order = TRUE" if needed!
  expect_named(A2T(f.spct), c("w.length", "A", "Tfr"))
  expect_named(A2T(f.spct, action = "add"), c("w.length", "A", "Tfr"))
  expect_named(A2T(f.spct, action = "replace"), c("w.length", "Tfr"))
})

test_that("T and A values", {
  fa.spct <- filter_spct(w.length = 300:320, A = 1)
  ft.spct <- filter_spct(w.length = 300:320, Tfr = 0.1)
  expect_equal(fa.spct, T2A(ft.spct, action = "replace"))
  expect_equal(A2T(fa.spct, action = "replace"), ft.spct)

  expect_equal(ft.spct, A2T(ft.spct))
  expect_equal(fa.spct, T2A(fa.spct))
})

test_that("T2Afr", {
  f.spct <- filter_spct(w.length = 300:320, Tfr = 0.1, Rfr = 0.2, Tfr.type = "total")
  expect_silent(T2Afr(f.spct))
  expect_silent(T2Afr(f.spct, clean = TRUE))
  expect_silent(T2Afr(f.spct, clean = FALSE))
  Tfr_as_default()
  expect_silent(T2Afr(f.spct))
  expect_silent(T2Afr(f.spct, clean = TRUE))
  Afr_as_default()
  expect_silent(T2Afr(f.spct))
  expect_silent(T2Afr(f.spct, clean = TRUE))
  A_as_default()
  expect_silent(T2Afr(f.spct))
  expect_silent(T2Afr(f.spct, clean = TRUE))
  unset_user_defaults()

  expect_true(all(T2Afr(f.spct)[["Afr"]] == 0.7))

  # add "ignore.order = TRUE" if needed!
  expect_named(T2Afr(f.spct), c("w.length", "Tfr", "Rfr", "Afr"))
  expect_named(T2Afr(f.spct, action = "add"), c("w.length", "Tfr", "Rfr", "Afr"))
  expect_named(T2Afr(f.spct, action = "replace"), c("w.length", "Rfr", "Afr"))

})

test_that("T2Afr", {
  f.spct <- filter_spct(w.length = 300:320, Tfr = 0.1, Rfr = 0.2, Tfr.type = "total")
  expect_silent(T2Afr(f.spct))
  expect_silent(T2Afr(f.spct, clean = TRUE))
  expect_silent(T2Afr(f.spct, clean = FALSE))
  Tfr_as_default()
  expect_silent(T2Afr(f.spct))
  expect_silent(T2Afr(f.spct, clean = TRUE))
  Afr_as_default()
  expect_silent(T2Afr(f.spct))
  expect_silent(T2Afr(f.spct, clean = TRUE))
  A_as_default()
  expect_silent(T2Afr(f.spct))
  expect_silent(T2Afr(f.spct, clean = TRUE))
  unset_user_defaults()

  expect_true(all(T2Afr(f.spct)[["Afr"]] == 0.7))

  # add "ignore.order = TRUE" if needed!
  expect_named(T2Afr(f.spct), c("w.length", "Tfr", "Rfr", "Afr"))
  expect_named(T2Afr(f.spct, action = "add"), c("w.length", "Tfr", "Rfr", "Afr"))
  expect_named(T2Afr(f.spct, action = "replace"), c("w.length", "Rfr", "Afr"))

})

test_that("Afr2T", {
  f.spct <- filter_spct(w.length = 300:320, Afr = 0.7, Rfr = 0.2, Tfr.type = "total")
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  Tfr_as_default()
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  Afr_as_default()
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  A_as_default()
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  unset_user_defaults()

  expect_true(all((Afr2T(f.spct)[["Tfr"]] - 0.1) < 0.0001))

  # add "ignore.order = TRUE" if needed!
  expect_named(Afr2T(f.spct), c("w.length", "Afr", "Rfr", "Tfr"))
  expect_named(Afr2T(f.spct, action = "add"), c("w.length", "Afr", "Rfr", "Tfr"))
  expect_named(Afr2T(f.spct, action = "replace"), c("w.length", "Rfr", "Tfr"))
})

test_that("Afr2T_internal", {
  f.spct <- filter_spct(w.length = 300:320, Afr = 0.9, Tfr.type = "internal")
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  Tfr_as_default()
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  Afr_as_default()
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  A_as_default()
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  unset_user_defaults()

  expect_true(all((Afr2T(f.spct)[["Tfr"]] - 0.1) < 0.0001))

  # add "ignore.order = TRUE" if needed!
  expect_named(Afr2T(f.spct), c("w.length", "Afr", "Tfr"))
  expect_named(Afr2T(f.spct, action = "add"), c("w.length", "Afr", "Tfr"))
  expect_named(Afr2T(f.spct, action = "replace"), c("w.length", "Tfr"))
})

test_that("Afr2T_properties", {
  f.spct <- filter_spct(w.length = 300:320, Afr = 0.7, Tfr.type = "total")
  filter_properties(f.spct) <- list(Rfr.constant = 0.2, thickness = NA_real_, attenuation.mode = NA_character_)
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  Tfr_as_default()
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  Afr_as_default()
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  A_as_default()
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  unset_user_defaults()

  expect_true(all((Afr2T(f.spct)[["Tfr"]] - 0.1) < 0.0001))

  # add "ignore.order = TRUE" if needed!
  expect_named(Afr2T(f.spct), c("w.length", "Afr", "Tfr"))
  expect_named(Afr2T(f.spct, action = "add"), c("w.length", "Afr", "Tfr"))
  expect_named(Afr2T(f.spct, action = "replace"), c("w.length", "Tfr"))
})

test_that("Afr2T_object_spct", {
  f.spct <- object_spct(w.length = 300:320, Afr = 0.7, Rfr = 0.2, Tfr.type = "total")
  expect_is(Afr2T(f.spct), "object_spct")
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  Tfr_as_default()
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  Afr_as_default()
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  A_as_default()
  expect_silent(Afr2T(f.spct))
  expect_silent(Afr2T(f.spct, clean = TRUE))
  unset_user_defaults()

  expect_true(all((Afr2T(f.spct)[["Tfr"]] - 0.1) < 0.0001))

  # add "ignore.order = TRUE" if needed!
  expect_named(Afr2T(f.spct), c("w.length", "Rfr", "Afr", "Tfr"))
  expect_named(Afr2T(f.spct, action = "add"), c("w.length", "Rfr", "Afr", "Tfr"))
  expect_named(Afr2T(f.spct, action = "replace"), c("w.length", "Rfr", "Tfr"))
})

test_that("any2Afr", {
  f.spct <- filter_spct(w.length = 300:320, Tfr = 0.1, Rfr = 0.2, Tfr.type = "total")
  expect_silent(any2Afr(f.spct))
  expect_silent(any2Afr(f.spct, clean = TRUE))
  expect_silent(any2Afr(f.spct, clean = FALSE))

  expect_true(all(any2Afr(f.spct)[["Afr"]] == 0.7))

  # add "ignore.order = TRUE" if needed!
  expect_named(any2Afr(f.spct), c("w.length", "Tfr", "Rfr", "Afr"))
  expect_named(any2Afr(f.spct, action = "add"), c("w.length", "Tfr", "Rfr", "Afr"))
  expect_named(any2Afr(f.spct, action = "replace"), c("w.length", "Rfr", "Afr"))

})

test_that("any2T_Afr", {
  f.spct <- filter_spct(w.length = 300:320, Afr = 0.7, Rfr = 0.2, Tfr.type = "total")
  expect_silent(any2T(f.spct))
  expect_silent(any2T(f.spct, clean = TRUE))

  expect_true(all((Afr2T(f.spct)[["Tfr"]] - 0.1) < 0.0001))

  # add "ignore.order = TRUE" if needed!
  expect_named(any2T(f.spct), c("w.length", "Afr", "Rfr", "Tfr"))
  expect_named(any2T(f.spct, action = "add"), c("w.length", "Afr", "Rfr", "Tfr"))
  expect_named(any2T(f.spct, action = "replace"), c("w.length", "Rfr", "Tfr"))
})

test_that("any2T_Tfr", {
  f.spct <- filter_spct(w.length = 300:320, Tfr = 0.1, Rfr = 0.2, Tfr.type = "total")
  expect_silent(any2T(f.spct))
  expect_silent(any2T(f.spct, clean = TRUE))

  expect_true(all(any2T(f.spct)[["Tfr"]] == 0.1))

  # add "ignore.order = TRUE" if needed!
  expect_named(any2T(f.spct), c("w.length", "Tfr", "Rfr"))
  expect_named(any2T(f.spct, action = "add"), c("w.length", "Tfr", "Rfr"))
  expect_named(any2T(f.spct, action = "replace"), c("w.length", "Tfr", "Rfr"))
})


test_that("any2T_A", {
  f.spct <- filter_spct(w.length = 300:320, A = 1)
  expect_silent(any2T(f.spct))
  expect_silent(any2T(f.spct, clean = TRUE))
  expect_true(all(any2T(f.spct)[["Tfr"]] == 0.1))

  # add "ignore.order = TRUE" if needed!
  expect_named(any2T(f.spct), c("w.length", "A", "Tfr"))
  expect_named(any2T(f.spct, action = "add"), c("w.length", "A", "Tfr"))
  expect_named(any2T(f.spct, action = "replace"), c("w.length", "Tfr"))
})

test_that("convertTfrType", {
  f.spct <- filter_spct(w.length = 300:400, Tfr = 0.5, Rfr = 0.5, Tfr.type = "total")
  expect_equal(convertTfrType(f.spct), f.spct)
  expect_silent(convertTfrType(f.spct, Tfr.type = "total"))
  expect_true(all(convertTfrType(f.spct, Tfr.type = "total")[["Tfr"]] == 0.5))
  expect_true(all(convertTfrType(f.spct, Tfr.type = "total")[["Rfr"]] == 0.5))
  # add "ignore.order = TRUE" if needed!
  expect_named(convertTfrType(f.spct), c("w.length", "Tfr", "Rfr"))

  expect_silent(convertTfrType(f.spct, Tfr.type = "internal"))
  expect_true(all(convertTfrType(f.spct, Tfr.type = "internal")[["Tfr"]] == 1))
  expect_true(all(convertTfrType(f.spct, Tfr.type = "internal")[["Rfr"]] == 0.5))
  # add "ignore.order = TRUE" if needed!
  expect_named(convertTfrType(f.spct, Tfr.type = "internal"), c("w.length", "Tfr", "Rfr"))
})

test_that("convertTfrType_properties", {
  f.spct <- filter_spct(w.length = 300:400, Tfr = 0.5, Tfr.type = "total")
  filter_properties(f.spct) <- list(Rfr.constant = 0.5, thickness = NA_real_, attenuation.mode = NA_character_)
  expect_equal(convertTfrType(f.spct), f.spct)
  expect_silent(convertTfrType(f.spct, Tfr.type = "total"))
  expect_true(all(convertTfrType(f.spct, Tfr.type = "total")[["Tfr"]] == 0.5))
  expect_true(all(convertTfrType(f.spct, Tfr.type = "total")[["Rfr"]] == 0.5))
  # add "ignore.order = TRUE" if needed!
  expect_named(convertTfrType(f.spct), c("w.length", "Tfr"))

  expect_silent(convertTfrType(f.spct, Tfr.type = "internal"))
  expect_true(all(convertTfrType(f.spct, Tfr.type = "internal")[["Tfr"]] == 1))
  expect_true(all(convertTfrType(f.spct, Tfr.type = "internal")[["Rfr"]] == 0.5))
  # add "ignore.order = TRUE" if needed!
  expect_named(convertTfrType(f.spct, Tfr.type = "internal"), c("w.length", "Tfr"))
})

test_that("convertTfrType", {
  f.spct <- filter_spct(w.length = 300:400, Tfr = 1, Rfr = 0.5, Tfr.type = "internal")
  expect_equal(convertTfrType(f.spct), f.spct)
  expect_silent(convertTfrType(f.spct, Tfr.type = "total"))
  expect_true(all(convertTfrType(f.spct, Tfr.type = "total")[["Tfr"]] == 0.5))
  expect_true(all(convertTfrType(f.spct, Tfr.type = "total")[["Rfr"]] == 0.5))
  # add "ignore.order = TRUE" if needed!
  expect_named(convertTfrType(f.spct), c("w.length", "Tfr", "Rfr"))

  expect_silent(convertTfrType(f.spct, Tfr.type = "internal"))
  expect_true(all(convertTfrType(f.spct, Tfr.type = "internal")[["Tfr"]] == 1))
  expect_true(all(convertTfrType(f.spct, Tfr.type = "internal")[["Rfr"]] == 0.5))
  # add "ignore.order = TRUE" if needed!
  expect_named(convertTfrType(f.spct, Tfr.type = "internal"), c("w.length", "Tfr", "Rfr"))
})


test_that("e2q", {
  s.spct <- source_spct(w.length = 300:400, s.e.irrad = 1)
  expect_silent(e2q(s.spct))
  photon_as_default()
  expect_silent(e2q(s.spct))
  energy_as_default()
  expect_silent(e2q(s.spct))
  unset_user_defaults()

  expect_equal_to_reference(e2q(s.spct), "e2q.rds")

  # add "ignore.order = TRUE" if needed!
  expect_named(e2q(s.spct), c("w.length", "s.e.irrad", "s.q.irrad"))
  expect_named(e2q(s.spct, action = "add"), c("w.length", "s.e.irrad", "s.q.irrad"))
  expect_named(e2q(s.spct, action = "replace"), c("w.length", "s.q.irrad"))

})

test_that("q2e", {
  s.spct <- source_spct(w.length = 300:400, s.q.irrad = 1)
  expect_silent(q2e(s.spct))
  photon_as_default()
  expect_silent(q2e(s.spct))
  energy_as_default()
  expect_silent(q2e(s.spct))
  unset_user_defaults()

  expect_equal_to_reference(q2e(s.spct), "q2e.rds")

  # add "ignore.order = TRUE" if needed!
  expect_named(q2e(s.spct), c("w.length", "s.q.irrad", "s.e.irrad"))
  expect_named(q2e(s.spct, action = "add"), c("w.length", "s.q.irrad", "s.e.irrad"))
  expect_named(q2e(s.spct, action = "replace"), c("w.length", "s.e.irrad"))

})

test_that("e and q values", {
  sq.spct <- source_spct(w.length = 300:400, s.q.irrad = 1)
  se.spct <- source_spct(w.length = 300:400, s.e.irrad = 1)

  expect_equal(sq.spct, e2q(sq.spct))
  expect_equal(se.spct, q2e(se.spct))
})

test_that("convertThickness_properties", {
  f.spct <- filter_spct(w.length = 300:400, Tfr = 0.5, Tfr.type = "internal")
  filter_properties(f.spct) <- list(Rfr.constant = 0.1, thickness = 1e-3, attenuation.mode = "absorption")
  expect_equal(convertThickness(f.spct), f.spct)
  expect_silent(convertThickness(f.spct, thickness = 2e-3))
  expect_true(all(convertThickness(f.spct, thickness = 2e-3)[["Tfr"]] == 0.25))
  # add "ignore.order = TRUE" if needed!
  expect_named(convertThickness(f.spct, thickness = 2e-3), c("w.length", "Tfr"))

  f.spct <- filter_spct(w.length = 300:400, Tfr = 0.4, Tfr.type = "total")
  filter_properties(f.spct) <- list(Rfr.constant = 0.1, thickness = 1e-3, attenuation.mode = "absorption")
  expect_silent(convertThickness(f.spct, thickness = 2e-3))
  expect_true(all(round(convertThickness(f.spct, thickness = 2e-3)[["Tfr"]], 6) == 0.177778))
  # add "ignore.order = TRUE" if needed!
  expect_named(convertThickness(f.spct, thickness = 2e-3), c("w.length", "Tfr"))
})
