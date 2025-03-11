library("photobiology")

context("extract")

test_that("source_spct", {

  my.spct <- source_spct(w.length = 400:450, s.e.irrad = 0.5, time.unit = "second")
  my_z.spct <- my.spct
  my_z.spct$z <- 1

  expect_equal(class(my.spct)[1], "source_spct")
  expect_equal(class(my.spct[ 1, ])[1], "source_spct")
  expect_equal(class(my.spct[ FALSE, ])[1], "source_spct")
  expect_equal(class(my.spct[ , -1])[1], "numeric")
  expect_equal(class(my.spct[ , 1])[1], "integer")
  expect_equal(class(my.spct[ FALSE, 1])[1], "integer")

  expect_equal(class(my_z.spct)[1], "source_spct")
  expect_equal(class(my_z.spct[ 1, ])[1], "source_spct")
  expect_equal(class(my_z.spct[ FALSE, ])[1], "source_spct")
  expect_equal(class(my_z.spct[ , 1])[1], "integer")
  expect_equal(class(my_z.spct[ , -1])[1], "tbl_df")
  expect_warning(my_z.spct[ , -2])
  expect_equal(suppressWarnings(class(my_z.spct[ , -2])[1]), "source_spct")
  expect_warning(my_z.spct[ 1 , -2])
  expect_equal(suppressWarnings(class(my_z.spct[ 1, -2])[1]), "source_spct")

  expect_equal(getTimeUnit(my_z.spct), getTimeUnit(my.spct))
  expect_equal(getTimeUnit(my.spct[1:10, ]), getTimeUnit(my.spct))
  expect_equal(getTimeUnit(my.spct[FALSE, ]), getTimeUnit(my.spct))
  expect_equal(getTimeUnit(my_z.spct[ , -3]), getTimeUnit(my.spct))
  expect_equal(getTimeUnit(my_z.spct[-(1:5) , -3]), getTimeUnit(my.spct))
  expect_equal(getTimeUnit(my_z.spct[ , -3]), getTimeUnit(my.spct))

  expect_equal(nrow(my_z.spct[ , -3]), nrow(my.spct))
  expect_equal(ncol(my_z.spct[ , -3]), ncol(my.spct))
  expect_equal(ncol(my.spct[FALSE , ]), ncol(my.spct))
  expect_equal(nrow(my.spct[FALSE , ]), 0)
})


context("replace")

test_that("source_spct", {

  my.spct <- source_spct(w.length = 400:450, s.e.irrad = 0.5, time.unit = "second")
  my_z.spct <- my.spct
  my_z.spct[1, 2] <- 1
  expect_equal(class(my_z.spct), class(my.spct))
  my_z.spct[ , 2] <- 1
  expect_equal(class(my_z.spct), class(my.spct))
  expect_silent(my_z.spct[ , 2] <- -1) # negative irradiances are valid!
  expect_equal(class(my_z.spct), class(my.spct))
  expect_silent(my_z.spct[ , 2] <- -(my.spct[ , 2])) # negative irradiances are valid!
  expect_equal(class(my_z.spct), class(my.spct))
  my_z.spct <- my.spct
  expect_error(my_z.spct[ , 1] <- -(my.spct[ , 1]))
  my_z.spct <- my.spct
  expect_error(my_z.spct[ , 1] <- 0)
  my_z.spct <- my.spct
  expect_error(my_z.spct$w.length <- 0)
  my_z.spct <- my.spct
  expect_error(my_z.spct[1 , 1] <- my.spct[ 2, 1])

  my_z.spct <- my.spct
  my_z.spct[1, 2] <- 1
  expect_equal(getTimeUnit(my_z.spct), getTimeUnit(my.spct))
  my_z.spct <- my.spct
  my_z.spct[ , 2] <- 1
  expect_equal(getTimeUnit(my_z.spct), getTimeUnit(my.spct))
})

