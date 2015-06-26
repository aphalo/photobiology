library("photobiology")
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

  expect_more_than(round(day_length(ymd("2014-03-21"), lat = 45, twilight = "civil"), 2), 12)
  expect_more_than(round(day_length(ymd("2014-09-21"), lat = 45, twilight = "civil"), 2), 12)
  expect_less_than(round(night_length(ymd("2014-03-21"), lat = 45, twilight = "civil"), 2), 12)
  expect_less_than(round(night_length(ymd("2014-09-21"), lat = 45, twilight = "civil"), 2), 12)

  expect_more_than(round(day_length(ymd("2014-03-21"), lat = 45, twilight = "nautical"), 2), 12)
  expect_more_than(round(day_length(ymd("2014-09-21"), lat = 45, twilight = "nautical"), 2), 12)
  expect_less_than(round(night_length(ymd("2014-03-21"), lat = 45, twilight = "nautical"), 2), 12)
  expect_less_than(round(night_length(ymd("2014-09-21"), lat = 45, twilight = "nautical"), 2), 12)

  expect_more_than(round(day_length(ymd("2014-03-21"), lat = 45, twilight = "astronomical"), 2), 12)
  expect_more_than(round(day_length(ymd("2014-09-21"), lat = 45, twilight = "astronomical"), 2), 12)
  expect_less_than(round(night_length(ymd("2014-03-21"), lat = 45, twilight = "astronomical"), 2), 12)
  expect_less_than(round(night_length(ymd("2014-09-21"), lat = 45, twilight = "astronomical"), 2), 12)

  expect_more_than(round(day_length(ymd("2014-03-21"), lat = 45, twilight = -1), 2), 12)
  expect_more_than(round(day_length(ymd("2014-09-21"), lat = 45, twilight = -1), 2), 12)
  expect_less_than(round(night_length(ymd("2014-03-21"), lat = 45, twilight = -1), 2), 12)
  expect_less_than(round(night_length(ymd("2014-09-21"), lat = 45, twilight = -1), 2), 12)

  expect_less_than(round(day_length(ymd("2014-03-21"), lat = 45, twilight = +1), 2), 12)
  expect_less_than(round(day_length(ymd("2014-09-21"), lat = 45, twilight = +1), 2), 12)
  expect_more_than(round(night_length(ymd("2014-03-21"), lat = 45, twilight = +1), 2), 12)
  expect_more_than(round(night_length(ymd("2014-09-21"), lat = 45, twilight = +1), 2), 12)

  expect_error(day_length(ymd("2014-03-21"), lat = 45, twilight = "bad"))
  expect_error(day_length(ymd("2014-03-21"), lat = 45, twilight = +91))
  expect_error(day_length(ymd("2014-03-21"), lat = 45, twilight = -91))

})


