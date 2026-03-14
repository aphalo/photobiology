library("lubridate")

context("fscale_spct")

test_that("fscale copies attrs correctly", {
  energy_as_default()
  my.spct <- q2e(sun.spct, action = "replace")

  expect_lt(abs(integrate_spct(fscale(my.spct, f = "total")) - 1), 1e-6)
  expect_lt(abs(average_spct(fscale(my.spct, f = "mean")) - 1), 1e-6)
  expect_warning(irrad(fscale(my.spct, f = "mean")))
#  behaviour changed in 0.13.3
#  expect_equivalent(suppressWarnings(irrad(fscale(my.spct, f = "mean"))), 520)
  expect_true(is.na(suppressWarnings(irrad(fscale(my.spct, f = "mean")))))
  expect_no_warning(irrad(fscale(my.spct, f = "mean"), allow.scaled = TRUE))
  expect_equivalent(irrad(fscale(my.spct, f = "mean"),
                          allow.scaled = TRUE), 520)
  expect_named(fscale(my.spct), names(my.spct))
  expect_equal(class(fscale(my.spct)), class(my.spct))
  expect_warning(fscale(my.spct, range = 100))
  expect_error(fscale(my.spct, range = 281))
  expect_warning(fscale(my.spct, range = c(100, 100)))
  expect_error(fscale(my.spct, range = c(281, 281)))
  expect_true(is_scaled(fscale(my.spct)))
  expect_false(is_scaled(my.spct))
  expect_equal(is.source_spct(fscale(my.spct)), is.source_spct(my.spct))
  expect_equal(is.filter_spct(fscale(my.spct)), is.filter_spct(my.spct))
  expect_equal(getTimeUnit(fscale(my.spct)), getTimeUnit(my.spct))
  expect_equal(comment(fscale(my.spct)), comment(my.spct))
})

context("integrate_spct average_spct")

test_that("integrate_spct handles attrs correctly", {
  energy_as_default()
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

context("invariant attributes after operation")

test_that("trim_wl handles attrs correctly", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct, range = 110:200)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct, range = 110:200,
                                                fill = 0)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct, range = 100:190)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct, range = 100:190,
                                                fill = 0)))),
               character(0) )

    expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct, range = 100:210)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct, range = 100:210,
                                                fill = 0)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(trim_wl(my.spct)))),
               character(0) )

})

test_that("clip_wl attr", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clip_wl(my.spct, range = 110:200)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clip_wl(my.spct, range = 100:190)))),
               character(0) )


  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clip_wl(my.spct, range = 100:210)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clip_wl(my.spct)))),
               character(0) )

})

test_that("clean handles attrs correctly", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clean(my.spct, range = 110:200)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clean(my.spct, range = 100:190)))),
               character(0) )


  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clean(my.spct, range = 100:210)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clean(my.spct)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(clean(my.spct, range.s.data = c(0,0.9))))),
               character(0) )
})

test_that("fshift handles attrs correctly", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 100:105)))),
               character(0) )


  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 195:200)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 145:150)))),
               character(0) )

})

test_that("fshift x2 handles attrs correctly", {

  my.spct <- raw_spct(w.length=100:200, counts = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 100:105)))),
               character(0) )


  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 195:200)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 145:150)))),
               character(0) )

  my.spct <- raw_spct(w.length=100:200, counts1 = 1, counts2 = 2)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 100:105)))),
               character(0) )


  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 195:200)))),
               character(0) )

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(fshift(my.spct, range = 145:150)))),
               character(0) )
})

test_that("fscale handles attrs correctly", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(fscale(my.spct))),
               names(attributes(my.spct))),
               "scaled" )

  expect_equal(setdiff(names(attributes(fscale(my.spct, range = 100:110))),
                       names(attributes(my.spct))),
               "scaled" )

  expect_equal(setdiff(names(attributes(fscale(my.spct, f = "total"))),
                       names(attributes(my.spct))),
               "scaled" )

})

test_that("normalize handles attrs correctly", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- validate_geocode(data.frame(lon = 24.93545, lat = 60.16952))
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(normalize(my.spct))),
                       c(names(attributes(my.spct)))),
               c("normalized", "normalization"))

  expect_equal(setdiff(names(attributes(normalize(my.spct, range = c(100, 150)))),
                       names(attributes(my.spct))),
               c("normalized", "normalization"))

  expect_equal(setdiff(names(attributes(normalize(my.spct, norm = "max"))),
                       names(attributes(my.spct))),
               c("normalized", "normalization"))

  expect_equal(setdiff(names(attributes(normalize(my.spct, norm = 130))),
                       names(attributes(my.spct))),
               c("normalized", "normalization"))
})

test_that("peaks handles attrs correctly", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(peaks(my.spct))),
                       names(attributes(my.spct))),
               character(0)  )

  expect_equal(setdiff(names(attributes(valleys(my.spct))),
                       names(attributes(my.spct))),
               character(0)  )

})

test_that("wls_at_target handles attrs correctly", {

  my.spct <- sun.spct[200:300, ]
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(wls_at_target(my.spct, 0.7))),
                       c(names(attributes(my.spct)),
                       "idfactor", "instr.desc", "instr.settings", "how.measured")),
               character(0)  )

  expect_equal(setdiff(names(attributes(wls_at_target(my.spct))),
                       c(names(attributes(my.spct)),
                         "idfactor", "instr.desc", "instr.settings", "how.measured")),
               character(0)  )

})

test_that("smooth_spct handles attrs correctly", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(setdiff(names(attributes(smooth_spct(my.spct))),
                       names(attributes(my.spct))), "comment"),
               character(0) )

})

test_that("extract handles attrs correctly", {

  my.spct <- source_spct(w.length=100:200, s.e.irrad = 1)
  tested.time <- ymd_hms("2015-12-31 23:59:59")
  setWhenMeasured(my.spct, tested.time)
  tested.location <- data.frame(lon = 24.93545, lat = 60.16952)
  setWhereMeasured(my.spct, tested.location)
  tested.what <- "user message"
  setWhatMeasured(my.spct, tested.what)

  expect_equal(setdiff(names(attributes(my.spct)),
                       names(attributes(my.spct[2:50, ]))),
               character(0) )

})

