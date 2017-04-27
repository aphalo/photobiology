library("photobiology")

context("trim_wl")

test_that("waveband", {

  my.spct <- source_spct(w.length = 400:450, s.e.irrad = 0.5, time.unit = "second")
  my.wb <- waveband(c(300,500))
  expect_equal(spread(trim_wl(my.wb, my.spct)), spread(my.spct))
  expect_equal(range(trim_wl(my.wb, my.spct)), range(my.spct))

  my.wb <- waveband(c(100, 150))
  expect_equal(trim_wl(my.wb, my.spct), waveband())

})

