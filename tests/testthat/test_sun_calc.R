library("photobiology")
library("lubridate")
# library("ggmap")

context("sun_calc")

test_that("sun_angles", {
  sun.angles <-
    sun_angles(ymd_hms("2012-12-22 12:00:00", tz = "UTC"),
               geocode = data.frame(lon = 0, lat = c(-41, -3, 3, 41)),
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

  sun.angles <-
    sun_angles(ymd_hms("2012-12-22 12:00:00", tz = "UTC"),
               geocode = data.frame(lon = 0, lat = c(-41, -3, 3, 41)),
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

test_that("sunrise_time", {
  expect_lt(
    abs(as.numeric(sunrise_time(ymd("2016-04-17", tz = "UTC"),
                            geocode = data.frame(lon = 24.93838, lat = 60.16986),
#                           geocode = geocode("Helsinki, Finland"),
                            twilight = "none") -
                 ymd_hms("2016-04-17 03:01:42", tz = "UTC"))), 1
  )
  expect_lt(
    abs(as.numeric(sunrise_time(ymd("2016-04-17", tz = "UTC"),
                            geocode = data.frame(lon = 25.46508, lat = 65.01209),
#                            geocode = geocode("Oulu, Finland"),
                            twilight = "none") -
                 ymd_hms("2016-04-17 02:40:36", tz = "UTC"))), 1
  )
  expect_lt(
    abs(as.numeric(sunrise_time(ymd("2016-04-17", tz = "UTC"),
                            geocode = data.frame(lon = 27.02853, lat = 69.90905),
#                           geocode = geocode("Utsjoki, Finland"),
                            twilight = "none") -
                 ymd_hms("2016-04-17 02:05:09", tz = "UTC"))), 1
  )
  expect_lt(
    abs(as.numeric(sunrise_time(ymd("2016-04-17", tz = "UTC"),
                            geocode = data.frame(lon = -68.30295, lat = -54.80191),
#                           geocode = geocode("Ushuaia, Argentina"),
                            twilight = "none") -
                 ymd_hms("2016-04-17 11:35:32", tz = "UTC"))), 1
  )
  expect_lt(
    as.numeric(sunrise_time(ymd("2016-04-17", tz = "UTC"),
                            geocode = data.frame(lon = 23.67027, lat = 77.5536),
                            #
                            twilight = "none") -
                 ymd_hms("2016-04-17 00:24:26", tz = "UTC")), 1
  )
  expect_lt(
    as.numeric(sunrise_time(ymd("2016-04-21", tz = "UTC"),
                            geocode = data.frame(lon = 23.67027, lat = 77.5536),
                            #
                            twilight = "none") -
                 ymd_hms("2016-04-20 23:10:19", tz = "UTC")), 1
  )
  expect_true(
    is.na(sunrise_time(ymd("2016-04-22", tz = "UTC"),
                            geocode = data.frame(lon = 23.67027, lat = 77.5536),
                            #
                            twilight = "none"))
  )
})

test_that("sunset_time", {
  expect_lt(
    as.numeric(sunset_time(ymd("2016-04-17", tz = "UTC"),
                           geocode = data.frame(lon = 24.93838, lat = 60.16986),
#                            geocode = geocode("Helsinki, Finland"),
                            twilight = "none") -
                 ymd_hms("2016-04-17 17:37:35", tz = "UTC")), 1
  )
  expect_lt(
    as.numeric(sunset_time(ymd("2016-04-17", tz = "UTC"),
                           geocode = data.frame(lon = 25.46508, lat = 65.01209),
#                           geocode = geocode("Oulu, Finland"),
                            twilight = "none") -
                 ymd_hms("2016-04-17 17:54:27", tz = "UTC")), 1
  )
  expect_lt(
    as.numeric(sunset_time(ymd("2016-04-17", tz = "UTC"),
                           geocode = data.frame(lon = 27.02853, lat = 69.90905),
#                          geocode = geocode("Utsjoki, Finland"),
                            twilight = "none") -
                 ymd_hms("2016-04-17 18:17:24", tz = "UTC")), 1
  )
})

test_that("sunrise_time_vectorized", {
  expect_equal(
    length(sunrise_time(ymd("2016-04-17", tz = "UTC") + days(0:5),
                        geocode = data.frame(lon = 24.93838, lat = 60.16986),
#                       geocode = geocode("Helsinki, Finland"),
                            twilight = "none")), 6)
  expect_equal(
    sunrise_time(ymd("2016-04-17", tz = "UTC") + days(0:5),
                        geocode = data.frame(lon = 24.93838, lat = 60.16986),
                        #                       geocode = geocode("Helsinki, Finland"),
                        twilight = "none")[1],
    sunrise_time(ymd("2016-04-17", tz = "UTC"),
                 geocode = data.frame(lon = 24.93838, lat = 60.16986),
                 #                       geocode = geocode("Helsinki, Finland"),
                 twilight = "none"))
  expect_equal(
    length(sunrise_time(ymd("2016-04-20", tz = "UTC") + days(0:5),
                            geocode = data.frame(lon = 23.67027, lat = 77.5536),
                            #                           geocode = geocode("Ushuaia, Argentina"),
                            twilight = "none")), 6)
  expect_equal(
    sunrise_time(ymd("2016-04-17", tz = "UTC") + days(0:5),
                 geocode = data.frame(lon = 23.67027, lat = 77.5536),
                 #                           geocode = geocode("Ushuaia, Argentina"),
                 twilight = "none")[1],
    sunrise_time(ymd("2016-04-17", tz = "UTC"),
                 geocode = data.frame(lon = 23.67027, lat = 77.5536),
                 #                           geocode = geocode("Ushuaia, Argentina"),
                 twilight = "none"))
  expect_equal(
    sunrise_time(ymd("2016-04-20", tz = "UTC") + days(0:5),
                 geocode = data.frame(lon = 23.67027, lat = 77.5536),
                 #                           geocode = geocode("Ushuaia, Argentina"),
                 twilight = "none")[3],
    sunrise_time(ymd("2016-04-22", tz = "UTC"),
                 geocode = data.frame(lon = 23.67027, lat = 77.5536),
                 #                           geocode = geocode("Ushuaia, Argentina"),
                 twilight = "none"))
})

test_that("sunset_time_vectorized", {
  expect_equal(
    length(sunset_time(ymd("2016-04-17", tz = "UTC") + days(0:5),
                       geocode = data.frame(lon = 24.93838, lat = 60.16986),
#                     geocode = geocode("Helsinki, Finland"),
                        twilight = "none")), 6)
})

test_that("daylength", {

  expect_equal(day_length(ymd("2014-12-21"),
                          geocode = data.frame(lat = 85, lon = 0)),
               0)
  expect_equal(day_length(ymd("2014-06-21"),
                          geocode = data.frame(lat = 85, lon = 0)),
               24)
  expect_equal(night_length(ymd("2014-12-21"),
                            geocode = data.frame(lat = 85, lon = 0)),
               24)
  expect_equal(night_length(ymd("2014-06-21"),
                            geocode = data.frame(lat = 85, lon = 0)),
               0)
  expect_equal(round(day_length(ymd("2014-12-21"),
                                geocode = data.frame(lat = 0, lon = 0)), 3),
               12.121)
  expect_equal(round(day_length(ymd("2014-06-21"),
                                geocode = data.frame(lat = 0, lon = 0)), 3),
               12.121)
  expect_equal(round(night_length(ymd("2014-12-21"),
                                  geocode = data.frame(lat = 0, lon = 0)), 3),
               11.879)
  expect_equal(round(night_length(ymd("2014-06-21"),
                                  geocode = data.frame(lat = 0, lon = 0)), 3),
               11.879)

  expect_equal(round(day_length(ymd("2014-03-20"),
                                geocode = data.frame(lat = 45, lon = 0)), 3),
               12.162)
  expect_equal(round(day_length(ymd("2014-09-22"),
                                geocode = data.frame(lat = 45, lon = 0)), 3),
               12.173)
  expect_equal(round(night_length(ymd("2014-03-20"),
                                  geocode = data.frame(lat = 45, lon = 0)), 3),
               11.838)
  expect_equal(round(night_length(ymd("2014-09-22"),
                                  geocode = data.frame(lat = 45, lon = 0)), 3),
               11.827)

  expect_gt(round(day_length(ymd("2014-03-21"),
                             geocode = data.frame(lat = 45, lon = 0),
                             twilight = "civil"), 2), 12)
  expect_gt(round(day_length(ymd("2014-09-21"),
                             geocode = data.frame(lat = 45, lon = 0),
                             twilight = "civil"), 2), 12)
  expect_lt(round(night_length(ymd("2014-03-21"),
                               geocode = data.frame(lat = 45, lon = 0),
                               twilight = "civil"), 2), 12)
  expect_lt(round(night_length(ymd("2014-09-21"),
                               geocode = data.frame(lat = 45, lon = 0),
                               twilight = "civil"), 2), 12)

  expect_gt(round(day_length(ymd("2014-03-21"),
                             geocode = data.frame(lat = 45, lon = 0),
                             twilight = "nautical"), 2), 12)
  expect_gt(round(day_length(ymd("2014-09-21"),
                             geocode = data.frame(lat = 45, lon = 0),
                             twilight = "nautical"), 2), 12)
  expect_lt(round(night_length(ymd("2014-03-21"),
                               geocode = data.frame(lat = 45, lon = 0),
                               twilight = "nautical"), 2), 12)
  expect_lt(round(night_length(ymd("2014-09-21"),
                               geocode = data.frame(lat = 45, lon = 0),
                               twilight = "nautical"), 2), 12)

  expect_gt(round(day_length(ymd("2014-03-21"),
                             geocode = data.frame(lat = 45, lon = 0),
                             twilight = "astronomical"), 2), 12)
  expect_gt(round(day_length(ymd("2014-09-21"),
                             geocode = data.frame(lat = 45, lon = 0),
                             twilight = "astronomical"), 2), 12)
  expect_lt(round(night_length(ymd("2014-03-21"),
                               geocode = data.frame(lat = 45, lon = 0),
                               twilight = "astronomical"), 2), 12)
  expect_lt(round(night_length(ymd("2014-09-21"),
                               geocode = data.frame(lat = 45, lon = 0),
                               twilight = "astronomical"), 2), 12)

  expect_gt(round(day_length(ymd("2014-03-21"),
                             geocode = data.frame(lat = 45, lon = 0),
                             twilight = -1), 2), 12)
  expect_gt(round(day_length(ymd("2014-09-21"),
                             geocode = data.frame(lat = 45, lon = 0),
                             twilight = -1), 2), 12)
  expect_lt(round(night_length(ymd("2014-03-21"),
                               geocode = data.frame(lat = 45, lon = 0),
                               twilight = -1), 2), 12)
  expect_lt(round(night_length(ymd("2014-09-21"),
                               geocode = data.frame(lat = 45, lon = 0),
                               twilight = -1), 2), 12)

  expect_lt(round(day_length(ymd("2014-03-21"),
                             geocode = data.frame(lat = 45, lon = 0),
                             twilight = +1), 2), 12)
  expect_lt(round(day_length(ymd("2014-09-21"),
                             geocode = data.frame(lat = 45, lon = 0),
                             twilight = +1), 2), 12)
  expect_gt(round(night_length(ymd("2014-03-21"),
                               geocode = data.frame(lat = 45, lon = 0),
                               twilight = +1), 2), 12)
  expect_gt(round(night_length(ymd("2014-09-21"),
                               geocode = data.frame(lat = 45, lon = 0),
                               twilight = +1), 2), 12)

  expect_error(day_length(ymd("2014-03-21"),
                          geocode = data.frame(lat = 45, lon = 0),
                          twilight = "bad"))
  expect_error(day_length(ymd("2014-03-21"),
                          geocode = data.frame(lat = 45, lon = 0),
                          twilight = +91))
  expect_error(day_length(ymd("2014-03-21"),
                          geocode = data.frame(lat = 45, lon = 0),
                          twilight = -91))

})

