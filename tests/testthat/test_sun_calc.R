library("photobiology")
library("lubridate")

context("sun_calc")

test_that("daylength", {

  expect_equal(day_length(ymd("2014-12-21"), lat = 85), 0)
  expect_equal(day_length(ymd("2014-06-21"), lat = 85), 24)
  expect_equal(night_length(ymd("2014-12-21"), lat = 85), 24)
  expect_equal(night_length(ymd("2014-06-21"), lat = 85), 0)
  expect_equal(round(day_length(ymd("2014-12-21"), lat = 0), 2), 12)
  expect_equal(round(day_length(ymd("2014-06-21"), lat = 0), 2), 12)
  expect_equal(round(night_length(ymd("2014-12-21"), lat = 0), 2), 12)
  expect_equal(round(night_length(ymd("2014-06-21"), lat = 0), 2), 12)

  expect_equal(round(day_length(ymd("2014-03-20"), lat = 45), 1), 12)
  expect_equal(round(day_length(ymd("2014-09-22"), lat = 45), 1), 12)
  expect_equal(round(night_length(ymd("2014-03-20"), lat = 45), 1), 12)
  expect_equal(round(night_length(ymd("2014-09-22"), lat = 45), 1), 12)

  expect_gt(round(day_length(ymd("2014-03-21"), lat = 45, twilight = "civil"), 2), 12)
  expect_gt(round(day_length(ymd("2014-09-21"), lat = 45, twilight = "civil"), 2), 12)
  expect_lt(round(night_length(ymd("2014-03-21"), lat = 45, twilight = "civil"), 2), 12)
  expect_lt(round(night_length(ymd("2014-09-21"), lat = 45, twilight = "civil"), 2), 12)

  expect_gt(round(day_length(ymd("2014-03-21"), lat = 45, twilight = "nautical"), 2), 12)
  expect_gt(round(day_length(ymd("2014-09-21"), lat = 45, twilight = "nautical"), 2), 12)
  expect_lt(round(night_length(ymd("2014-03-21"), lat = 45, twilight = "nautical"), 2), 12)
  expect_lt(round(night_length(ymd("2014-09-21"), lat = 45, twilight = "nautical"), 2), 12)

  expect_gt(round(day_length(ymd("2014-03-21"), lat = 45, twilight = "astronomical"), 2), 12)
  expect_gt(round(day_length(ymd("2014-09-21"), lat = 45, twilight = "astronomical"), 2), 12)
  expect_lt(round(night_length(ymd("2014-03-21"), lat = 45, twilight = "astronomical"), 2), 12)
  expect_lt(round(night_length(ymd("2014-09-21"), lat = 45, twilight = "astronomical"), 2), 12)

  expect_gt(round(day_length(ymd("2014-03-21"), lat = 45, twilight = -1), 2), 12)
  expect_gt(round(day_length(ymd("2014-09-21"), lat = 45, twilight = -1), 2), 12)
  expect_lt(round(night_length(ymd("2014-03-21"), lat = 45, twilight = -1), 2), 12)
  expect_lt(round(night_length(ymd("2014-09-21"), lat = 45, twilight = -1), 2), 12)

  expect_lt(round(day_length(ymd("2014-03-21"), lat = 45, twilight = +1), 2), 12)
  expect_lt(round(day_length(ymd("2014-09-21"), lat = 45, twilight = +1), 2), 12)
  expect_gt(round(night_length(ymd("2014-03-21"), lat = 45, twilight = +1), 2), 12)
  expect_gt(round(night_length(ymd("2014-09-21"), lat = 45, twilight = +1), 2), 12)

  expect_error(day_length(ymd("2014-03-21"), lat = 45, twilight = "bad"))
  expect_error(day_length(ymd("2014-03-21"), lat = 45, twilight = +91))
  expect_error(day_length(ymd("2014-03-21"), lat = 45, twilight = -91))

})

test_that("daylength.geocode", {

  expect_equal(day_length(ymd("2014-12-21"), lat = 85, lon = 0),
               day_length(ymd("2014-12-21"), geocode = data.frame(lon = 0, lat = 85)))

  expect_equal(day_length(ymd("2014-12-21"), geocode = data.frame(lon = 0, lat = 85)), 0)
  expect_equal(day_length(ymd("2014-06-21"), geocode = data.frame(lon = 0, lat = 85)), 24)
  expect_equal(night_length(ymd("2014-12-21"), geocode = data.frame(lon = 0, lat = 85)), 24)
  expect_equal(night_length(ymd("2014-06-21"), geocode = data.frame(lon = 0, lat = 85)), 0)
  expect_equal(round(day_length(ymd("2014-12-21"),
                                geocode = data.frame(lon = 0, lat = 0)), 2), 12)
  expect_equal(round(day_length(ymd("2014-06-21"),
                                geocode = data.frame(lon = 0, lat = 0)), 2), 12)
  expect_equal(round(night_length(ymd("2014-12-21"),
                                  geocode = data.frame(lon = 0, lat = 0)), 2), 12)
  expect_equal(round(night_length(ymd("2014-06-21"),
                                  geocode = data.frame(lon = 0, lat = 0)), 2), 12)

  expect_equal(round(day_length(ymd("2014-03-20"),
                                geocode = data.frame(lon = 0, lat = 45)), 1), 12)

  expect_error(day_length(ymd("2014-03-20"), geocode = "zzz"))
  expect_error(day_length(ymd("2014-03-20"), geocode = "45"))
#  expect_error(day_length(ymd("2014-03-20"), lat = 91))
#  expect_error(day_length(ymd("2014-03-20"), lon = 180))
})
