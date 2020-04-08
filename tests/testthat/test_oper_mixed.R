library("photobiology")
context("source.spct and reflector.spct")

test_that("oper", {

  my.rf.spct <- reflector_spct(w.length = 400:409, Rfr = 0.5)
  my.insc.spct <- source_spct(w.length = 400:409, s.e.irrad = 2)

  expect_equivalent(my.insc.spct * my.rf.spct,
                    source_spct(w.length = 400:409, s.e.irrad = 1))
  expect_equivalent(my.rf.spct * my.insc.spct,
                    source_spct(w.length = 400:409, s.e.irrad = 1))
  expect_warning(my.rf.spct / my.insc.spct)
  expect_equivalent(suppressWarnings(my.rf.spct / my.insc.spct),
                    source_spct(w.length = 400:409, s.e.irrad = 0.25))
  expect_equivalent(my.insc.spct / my.rf.spct,
                    source_spct(w.length = 400:409, s.e.irrad = 4))
  expect_error(my.rf.spct + my.insc.spct)
#  expect_true(is.na(suppressWarnings(my.rf.spct + my.insc.spct)))
  expect_error(my.rf.spct - my.insc.spct)
#  expect_true(is.na(suppressWarnings(my.rf.spct - my.insc.spct)))
  expect_error(my.insc.spct + my.rf.spct)
#  expect_true(is.na(suppressWarnings(my.insc.spct + my.rf.spct)))
  expect_error(my.insc.spct - my.rf.spct)
#  expect_true(is.na(suppressWarnings(my.insc.spct - my.rf.spct)))
})

context("source.spct and filter.spct")

test_that("oper", {

  my.ft.spct <- filter_spct(w.length = 400:409, Tfr = 0.5)
  my.insc.spct <- source_spct(w.length = 400:409, s.e.irrad = 2)

  expect_equivalent(my.insc.spct * my.ft.spct,
                    source_spct(w.length = 400:409, s.e.irrad = 1))
  expect_equivalent(my.ft.spct * my.insc.spct,
                    source_spct(w.length = 400:409, s.e.irrad = 1))
  expect_warning(my.ft.spct / my.insc.spct)
  expect_equivalent(suppressWarnings(my.ft.spct / my.insc.spct),
                    source_spct(w.length = 400:409, s.e.irrad = 0.25))
  expect_equivalent(my.insc.spct / my.ft.spct,
                    source_spct(w.length = 400:409, s.e.irrad = 4))
  expect_error(my.ft.spct + my.insc.spct)
#  expect_true(is.na(suppressWarnings(my.ft.spct + my.insc.spct)))
  expect_error(my.ft.spct - my.insc.spct)
#  expect_true(is.na(suppressWarnings(my.ft.spct - my.insc.spct)))
  expect_error(my.insc.spct + my.ft.spct)
#  expect_true(is.na(suppressWarnings(my.insc.spct + my.ft.spct)))
  expect_error(my.insc.spct - my.ft.spct)
#  expect_true(is.na(suppressWarnings(my.insc.spct - my.ft.spct)))
})

