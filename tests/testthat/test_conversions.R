library("photobiology")

context("conversions")

test_that("T2A", {
  f.spct <- filter_spct(w.length = 300:400, Tfr = 0.1)
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
  f.spct <- filter_spct(w.length = 300:400, A = 1)
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
  fa.spct <- filter_spct(w.length = 300:400, A = 1)
  ft.spct <- filter_spct(w.length = 300:400, Tfr = 0.1)
  expect_equal(fa.spct, T2A(ft.spct, action = "replace"))
  expect_equal(A2T(fa.spct, action = "replace"), ft.spct)

  expect_equal(ft.spct, A2T(ft.spct))
  expect_equal(fa.spct, T2A(fa.spct))
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

