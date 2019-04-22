library("photobiology")
library("lubridate")

context("sun_calc")

test_that("sun_angles_24h", {
  #  test.path <- tempfile()
  test.path <- "./data/sun-angles-test-value"

  testthat::expect_known_value(
    sun_angles(time = ymd_hms("2012-10-22 12:00:00", tz = "UTC") + hours(0:24),
               geocode = data.frame(lon = 0, lat = c(89, 60, 45, 30, 0),
                                    address = "test"),
               use.refraction = FALSE),
    file = test.path
  )
}
)

test_that("sun_angles_geocode_vectorized", {
  sun.angles <-
    sun_angles(ymd_hms("2012-12-22 12:00:00", tz = "UTC"),
               geocode = data.frame(lon = 0, lat = c(-41, -3, 3, 41)),
               use.refraction = FALSE)
  expect_equal(sun.angles$longitude, c(0, 0, 0, 0))
  expect_equal(sun.angles$latitude, c(-41, -3, 3, 41))
  expect_true(all(abs(sun.angles$azimuth -
                        c(359.0886, 180.7880, 180.6180, 180.3050)) < 1e-4))
  expect_true(all(abs(sun.angles$elevation -
                        c(72.43013, 69.56602, 63.56646, 25.56747)) < 1e-4))

  sun.angles <-
    sun_angles(ymd_hms("2012-12-22 12:00:00", tz = "UTC"),
               geocode = data.frame(lon = 0, lat = c(-41, -3, 3, 41)),
               use.refraction = TRUE)
  expect_equal(sun.angles$longitude, c(0, 0, 0, 0))
  expect_equal(sun.angles$latitude, c(-41, -3, 3, 41))
  expect_true(all(abs(sun.angles$azimuth -
                        c(359.0886, 180.7880, 180.6180, 180.3050)) < 1e-4))
  expect_true(all(abs(sun.angles$elevation -
                        c(72.43523, 69.57203, 63.57448, 25.60103)) < 1e-4))
}
)

test_that("sun_angles_times_vectorized", {
  sun.angles <-
    sun_angles(ymd_hms("2012-10-22 12:00:00", tz = "UTC") + weeks(0:3),
               geocode = data.frame(lon = 0, lat = 41),
               use.refraction = FALSE)
  expect_equal(sun.angles$longitude, c(0, 0, 0, 0))
  expect_equal(sun.angles$latitude, c(41, 41, 41, 41))
  expect_true(all(abs(sun.angles$azimuth -
                        c(184.8340, 184.8604, 184.7117, 184.3951)) < 1e-4))
  expect_true(all(abs(sun.angles$elevation -
                        c(37.58375, 35.19321, 32.98920, 31.01003)) < 1e-4))

  sun.angles <-
    sun_angles(ymd_hms("2012-10-22 12:00:00", tz = "UTC") + weeks(0:3),
               geocode = data.frame(lon = 0, lat = 41),
               use.refraction = TRUE)
  expect_equal(sun.angles$longitude, c(0, 0, 0, 0))
  expect_equal(sun.angles$latitude, c(41, 41, 41, 41))
  expect_true(all(abs(sun.angles$azimuth -
                        c(184.8340, 184.8604, 184.7117, 184.3951)) < 1e-4))
  expect_true(all(abs(sun.angles$elevation -
                        c(37.60468, 35.21603, 33.01399, 31.03679)) < 1e-4))


  sun.angles <-
    sun_angles(ymd_hms("2012-10-22 16:30:00", tz = "UTC") + weeks(0:3),
               geocode = data.frame(lon = 0, lat = 41),
               use.refraction = FALSE)
  expect_equal(sun.angles$longitude, c(0, 0, 0, 0))
  expect_equal(sun.angles$latitude, c(41, 41, 41, 41))
  expect_true(all(abs(sun.angles$azimuth -
                        c(249.1641, 247.5674, 245.9957, 244.4684)) < 1e-4))
  expect_true(all(abs(sun.angles$elevation -
                        c(6.123102, 4.343891, 2.798342, 1.519157)) < 1e-4))

  sun.angles <-
    sun_angles(ymd_hms("2012-10-22 16:30:00", tz = "UTC") + weeks(0:3),
               geocode = data.frame(lon = 0, lat = 41),
               use.refraction = TRUE)
  expect_equal(sun.angles$longitude, c(0, 0, 0, 0))
  expect_equal(sun.angles$latitude, c(41, 41, 41, 41))
  expect_true(all(abs(sun.angles$azimuth -
                        c(249.1641, 247.5674, 245.9957, 244.4684)) < 1e-4))
  expect_true(all(abs(sun.angles$elevation -
                        c(6.259474, 4.521639, 3.036655, 1.837319)) < 1e-4))
}
)

test_that("sunrise_time", {
  expect_lt(
    abs(as.numeric(sunrise_time(ymd("2016-04-17", tz = "UTC"), tz = "UTC",
                                geocode = data.frame(lon = 24.93838, lat = 60.16986,
                                                     address = "Helsinki, Finland"),
                                twilight = "none") -
                     ymd_hms("2016-04-17 03:02:32", tz = "UTC"))), 1
  )
  expect_lt(
    abs(as.numeric(sunrise_time(ymd("2016-04-17", tz = "UTC"), tz = "UTC",
                                geocode = data.frame(lon = 25.46508, lat = 65.01209,
                                                     address = "Oulu, Finland"),
                                twilight = "none") -
                     ymd_hms("2016-04-17 02:41:40", tz = "UTC"))), 1
  )
  expect_lt(
    abs(as.numeric(sunrise_time(ymd("2016-04-17", tz = "UTC"), tz = "UTC",
                                geocode = data.frame(lon = 27.02853, lat = 69.90905,
                                                     adress = "Utsjoki, Finland"),
                                twilight = "none") -
                     ymd_hms("2016-04-17 02:06:34", tz = "UTC"))), 1
  )
  expect_lt(
    abs(as.numeric(sunrise_time(ymd("2016-04-17", tz = "UTC"), tz = "UTC",
                                geocode = data.frame(lon = -68.30295, lat = -54.80191,
                                                     address = "Ushuaia, Argentina"),
                                twilight = "none") -
                     ymd_hms("2016-04-17 11:34:59", tz = "UTC"))), 1
  )
  expect_lt(
    as.numeric(sunrise_time(ymd("2016-04-17", tz = "UTC"), tz = "UTC",
                            geocode = data.frame(lon = 23.67027, lat = 77.5536),
                            twilight = "none") -
                 ymd_hms("2016-04-17 00:28:16", tz = "UTC")), 1
  )
  expect_lt(
    as.numeric(sunrise_time(ymd("2016-04-21", tz = "UTC"), tz = "UTC",
                            geocode = data.frame(lon = 23.67027, lat = 77.5536),
                            twilight = "none") -
                 ymd_hms("2016-04-20 23:18:52", tz = "UTC")), 1
  )
  expect_true(
    is.na(sunrise_time(ymd("2016-04-23", tz = "UTC"), tz = "UTC",
                       geocode = data.frame(lon = 23.67027, lat = 77.5536),
                       twilight = "none"))
  )
})

test_that("sunset_time", {
  expect_lt(
    as.numeric(sunset_time(ymd("2016-04-17", tz = "UTC"), tz = "UTC",
                           geocode = data.frame(lon = 24.93838, lat = 60.16986),
                           #                            geocode = geocode("Helsinki, Finland"),
                           twilight = "none") -
                 ymd_hms("2016-04-17 17:36:52", tz = "UTC")), 1
  )
  expect_lt(
    as.numeric(sunset_time(ymd("2016-04-17", tz = "UTC"), tz = "UTC",
                           geocode = data.frame(lon = 25.46508, lat = 65.01209),
                           #                           geocode = geocode("Oulu, Finland"),
                           twilight = "none") -
                 ymd_hms("2016-04-17 17:53:31", tz = "UTC")), 1
  )
  expect_lt(
    as.numeric(sunset_time(ymd("2016-04-17", tz = "UTC"), tz = "UTC",
                           geocode = data.frame(lon = 27.02853, lat = 69.90905),
                           #                          geocode = geocode("Utsjoki, Finland"),
                           twilight = "none") -
                 ymd_hms("2016-04-17 18:16:06", tz = "UTC")), 1
  )
})

test_that("sunrise_time_vectorized", {
  expect_equal(
    length(sunrise_time(ymd("2016-04-17", tz = "UTC") + days(0:5),
                        geocode = data.frame(lon = 24.93838, lat = 60.16986),
                        #                       geocode = geocode("Helsinki, Finland"),
                        twilight = "none")), 6)
  expect_equal(
    sunrise_time(ymd("2016-04-17", tz = "UTC") + days(0:5), tz = "UTC",
                 geocode = data.frame(lon = 24.93838, lat = 60.16986),
                 #                       geocode = geocode("Helsinki, Finland"),
                 twilight = "none")[1],
    sunrise_time(ymd("2016-04-17", tz = "UTC"),
                 geocode = data.frame(lon = 24.93838, lat = 60.16986), tz = "UTC",
                 #                       geocode = geocode("Helsinki, Finland"),
                 twilight = "none"))
  expect_equal(
    length(sunrise_time(ymd("2016-04-20", tz = "UTC") + days(0:5),
                        geocode = data.frame(lon = 23.67027, lat = 77.5536),
                        #                           geocode = geocode("Ushuaia, Argentina"),
                        twilight = "none")), 6)
  expect_equal(
    sunrise_time(ymd("2016-04-17", tz = "UTC") + days(0:5), tz = "UTC",
                 geocode = data.frame(lon = 23.67027, lat = 77.5536),
                 #                           geocode = geocode("Ushuaia, Argentina"),
                 twilight = "none")[1],
    sunrise_time(ymd("2016-04-17", tz = "UTC"), tz = "UTC",
                 geocode = data.frame(lon = 23.67027, lat = 77.5536),
                 #                           geocode = geocode("Ushuaia, Argentina"),
                 twilight = "none"))
  expect_equal(
    sunrise_time(ymd("2016-04-20", tz = "UTC") + days(0:5), tz = "UTC",
                 geocode = data.frame(lon = 23.67027, lat = 77.5536),
                 #                           geocode = geocode("Ushuaia, Argentina"),
                 twilight = "none")[3],
    sunrise_time(ymd("2016-04-22", tz = "UTC"), tz = "UTC",
                 geocode = data.frame(lon = 23.67027, lat = 77.5536),
                 #                           geocode = geocode("Ushuaia, Argentina"),
                 twilight = "none"))
})

test_that("sunset_time_vectorized", {
  expect_equal(
    length(sunset_time(ymd("2016-04-17", tz = "UTC") + days(0:5), tz = "UTC",
                       geocode = data.frame(lon = 24.93838, lat = 60.16986),
                       #                     geocode = geocode("Helsinki, Finland"),
                       twilight = "none")), 6)
})

test_that("daylength", {

  expect_equal(day_length(now(tzone = "UTC")), day_length(today(tzone = "UTC"))) # is conversion o.k.?
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
               12.147)
  expect_equal(round(day_length(ymd("2014-09-22"),
                                geocode = data.frame(lat = 45, lon = 0)), 3),
               12.188)
  expect_equal(round(night_length(ymd("2014-03-20"),
                                  geocode = data.frame(lat = 45, lon = 0)), 3),
               11.853)
  expect_equal(round(night_length(ymd("2014-09-22"),
                                  geocode = data.frame(lat = 45, lon = 0)), 3),
               11.812)

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

