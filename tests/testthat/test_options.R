library("photobiology")

context("options")

test_that("verbose", {
  expect_silent(unset_user_defaults())
  expect_silent(old <- verbose_as_default())
  expect_equal(old$photobiology.verbose, getOption("verbose"))
  expect_equal(getOption("photobiology.verbose"), TRUE)
  expect_silent(old <- verbose_as_default(FALSE))
  expect_true(old$photobiology.verbose)
  expect_equal(getOption("photobiology.verbose"), FALSE)
  expect_silent(unset_user_defaults())
})

test_that("strict.range", {
  expect_silent(unset_user_defaults())
  expect_silent(old <- strict_range_as_default())
  expect_null(old$photobiology.strict.range)
  expect_equal(getOption("photobiology.strict.range"), TRUE)
  expect_silent(old <- strict_range_as_default(FALSE))
  expect_true(old$photobiology.strict.range)
  expect_equal(getOption("photobiology.strict.range"), FALSE)
  expect_silent(old <- strict_range_as_default(NA_integer_))
  expect_false(old$photobiology.strict.range)
  expect_true(is.na(getOption("photobiology.strict.range")))
  expect_silent(unset_user_defaults())
})

test_that("filter.qty", {
  expect_silent(unset_user_defaults())
  expect_silent(old <- Tfr_as_default())
  expect_null(old$photobiology.filter.qty)
  expect_equal(getOption("photobiology.filter.qty"), "transmittance")
  expect_silent(old <- Afr_as_default())
  expect_equal(old$photobiology.filter.qty, "transmittance")
  expect_equal(getOption("photobiology.filter.qty"), "absorptance")
  expect_silent(old <- A_as_default())
  expect_equal(old$photobiology.filter.qty, "absorptance")
  expect_equal(getOption("photobiology.filter.qty"), "absorbance")
  expect_silent(unset_user_defaults())

  expect_equal(using_Tfr(getOption("photobiology.filter.qty")), "transmittance")
  expect_equal(using_Afr(getOption("photobiology.filter.qty")), "absorptance")
  expect_equal(using_A(getOption("photobiology.filter.qty")), "absorbance")
})


test_that("radiation.unit", {
  expect_silent(unset_user_defaults())
  expect_silent(old <- photon_as_default())
  expect_null(old$photobiology.radiation.unit)
  expect_equal(getOption("photobiology.radiation.unit"), "photon")
  expect_silent(old <- energy_as_default())
  expect_equal(old$photobiology.radiation.unit, "photon")
  expect_equal(getOption("photobiology.radiation.unit"), "energy")
  expect_silent(old <- quantum_as_default())
  expect_equal(old$photobiology.radiation.unit, "energy")
  expect_equal(getOption("photobiology.radiation.unit"), "photon")
  expect_silent(unset_user_defaults())

  expect_equal(using_energy(getOption("photobiology.radiation.unit")), "energy")
  expect_equal(using_photon(getOption("photobiology.radiation.unit")), "photon")
  expect_equal(using_quantum(getOption("photobiology.radiation.unit")), "photon")
})

test_that("waveband.trim", {
  expect_silent(unset_user_defaults())
  expect_silent(old <- wb_trim_as_default())
  expect_null(old$photobiology.waveband.trim)
  expect_equal(getOption("photobiology.waveband.trim"), TRUE)
  expect_silent(old <- wb_trim_as_default(FALSE))
  expect_true(old$photobiology.waveband.trim)
  expect_equal(getOption("photobiology.waveband.trim"), FALSE)
  expect_silent(unset_user_defaults())
})

test_that("use.cached.mult", {
  expect_silent(unset_user_defaults())
  expect_silent(old <- use_cached_mult_as_default())
  expect_null(old$photobiology.use.cached.mult)
  expect_equal(getOption("photobiology.use.cached.mult"), TRUE)
  expect_silent(old <- use_cached_mult_as_default(FALSE))
  expect_true(old$photobiology.use.cached.mult)
  expect_equal(getOption("photobiology.use.cached.mult"), FALSE)
  expect_silent(unset_user_defaults())
})


