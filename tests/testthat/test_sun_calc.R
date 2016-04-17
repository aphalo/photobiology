library("photobiology")
library("lubridate")
library("ggmap")

context("sun_calc")

test_that("sun_angles", {
  sun.angles <- sun_angles(ymd_hms("2012-12-22 12:00:00", tz = "UTC"),
                           lon = rep(0, 4), lat = c(-41, -3, 3, 41),
                           use_refraction = FALSE)
  expect_equal(sun.angles$longitude, c(0, 0, 0, 0))
  expect_equal(sun.angles$latitude, c(-41, -3, 3, 41))
  expect_true(all(abs(sun.angles$azimuth -
                        c(359.0787, 180.7965, 180.6247, 180.3083)) < 1e-4))
  expect_true(all(abs(sun.angles$elevation -
                        c(72.43112, 69.56493, 63.56539, 25.56642)) < 1e-4))
  expect_true(all(abs(sun.angles$diameter -
                  c(0.542063, 0.542063, 0.542063, 0.542063)) < 1e-4))
  expect_true(all(abs(sun.angles$distance -
                    c(0.9836494, 0.9836494, 0.9836494, 0.9836494)) < 1e-4))

  sun.angles <- sun_angles(ymd_hms("2012-12-22 12:00:00", tz = "UTC"),
                           lon = rep(0, 4), lat = c(-41, -3, 3, 41),
                           use_refraction = TRUE)
  expect_equal(sun.angles$longitude, c(0, 0, 0, 0))
  expect_equal(sun.angles$latitude, c(-41, -3, 3, 41))
  expect_true(all(abs(sun.angles$azimuth -
                        c(359.0787, 180.7965, 180.6247, 180.3083)) < 1e-4))
  expect_true(all(abs(sun.angles$elevation -
                        c(72.43616, 69.57086, 63.57329, 25.59966)) < 1e-4))
  expect_true(all(abs(sun.angles$diameter -
                        c(0.542063, 0.542063, 0.542063, 0.542063)) < 1e-4))
  expect_true(all(abs(sun.angles$distance -
                        c(0.9836494, 0.9836494, 0.9836494, 0.9836494)) < 1e-4))
}
)

# test_that("sunrise_time", {
#   expect_lt(
#     sunrise_time(ymd("2016-04-17"),
#                  geocode = geocode("Helsinki, Finland"),
#                  twilight = "none") -
#     ymd_hms("2016-04-17 03:03:37", tz = "UTC"), duration(1, "seconds")
#   )
#   expect_lt(
#     sunrise_time(ymd("2016-04-17"),
#                geocode = geocode("Oulu, Finland"),
#                twilight = "none") -
#     ymd_hms("2016-04-17 02:43:04", tz = "UTC"), duration(1, "seconds")
#   )
#   expect_lt(
#       sunrise_time(ymd("2016-04-17"),
#                geocode = geocode("Utsjoki, Finland"),
#                twilight = "none") -
#       ymd_hms("2016-04-17 02:08:33", tz = "UTC"), duration(1, "seconds")
#   )
# })

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
  expect_error(day_length(ymd("2014-03-20"), lat = 91))
  expect_error(day_length(ymd("2014-03-20"), lon = 181))
  expect_error(day_length(ymd("2014-03-20"), lat = -91))
  expect_error(day_length(ymd("2014-03-20"), lon = -181))
})
