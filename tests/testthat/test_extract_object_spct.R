library("photobiology")

context("extract")

test_that("object_spct", {

  my.spct <- object_spct(w.length = 400:450,
                         Tfr = 0.5, Tfr.type = "total",
                         Rfr = 0.5, Rfr.type = "total")
  my_z.spct <- my.spct
  my_z.spct$z <- 1

  expect_equal(class(my.spct)[1], "object_spct")
  expect_equal(class(my.spct[ 1, ])[1], "object_spct")
  expect_equal(class(my.spct[ FALSE, ])[1], "object_spct")
  expect_equal(class(my.spct[ , -1])[1], "tbl_df")
  expect_equal(class(my.spct[ , 1])[1], "integer")
  expect_equal(class(my.spct[ FALSE, 1])[1], "integer")

  expect_equal(class(my_z.spct)[1], "object_spct")
  expect_equal(class(my_z.spct[ 1, ])[1], "object_spct")
  expect_equal(class(my_z.spct[ FALSE, ])[1], "object_spct")
  expect_equal(class(my_z.spct[ , 1])[1], "integer")
  expect_equal(class(my_z.spct[ , -1])[1], "tbl_df")
  expect_equal(suppressWarnings(class(my_z.spct[ , -2])[1]), "object_spct")
  expect_warning(my_z.spct[ , -2])
#  expect_equal(class(my_z.spct[ FALSE, -2])[1], "object_spct")
#  expect_warning(my_z.spct[ FALSE , -2])

  expect_equal(getTfrType(my_z.spct), getTfrType(my.spct))
  expect_equal(getTfrType(my.spct[1:10, ]), getTfrType(my.spct))
  expect_equal(getTfrType(my.spct[FALSE, ]), getTfrType(my.spct))
  expect_equal(getTfrType(my_z.spct[ , -4]), getTfrType(my.spct))
  expect_equal(getTfrType(my_z.spct[-(1:5) , -4]), getTfrType(my.spct))
  expect_equal(getTfrType(my_z.spct[ , -4]), getTfrType(my.spct))

  expect_equal(getRfrType(my_z.spct), getRfrType(my.spct))
  expect_equal(getRfrType(my.spct[1:10, ]), getRfrType(my.spct))
  expect_equal(getRfrType(my.spct[FALSE, ]), getRfrType(my.spct))
  expect_equal(getRfrType(my_z.spct[ , -4]), getRfrType(my.spct))
  expect_equal(getRfrType(my_z.spct[-(1:5) , -4]), getRfrType(my.spct))
  expect_equal(getRfrType(my_z.spct[ , -4]), getRfrType(my.spct))

  expect_equal(nrow(my_z.spct[ , -4]), nrow(my.spct))
  expect_equal(ncol(my_z.spct[ , -4]), ncol(my.spct))
  expect_equal(ncol(my.spct[FALSE , ]), ncol(my.spct))
  expect_equal(nrow(my.spct[FALSE , ]), 0)
})


context("replace")

test_that("object_spct", {

  my.spct <- object_spct(w.length = 400:450,
                         Tfr = 0.5, Tfr.type = "total",
                         Rfr = 0.5, Rfr.type = "total")
  my_z.spct <- my.spct
  my_z.spct[1, 2] <- 1
  expect_equal(class(my_z.spct), class(my.spct))
  my_z.spct[ , 2] <- 1
  expect_equal(class(my_z.spct), class(my.spct))
  expect_warning(my_z.spct[ , 2] <- -1)
  expect_equal(class(my_z.spct), class(my.spct))
  expect_warning(my_z.spct[ , 2] <- -(my.spct[ , 2]))
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
  expect_equal(getTfrType(my_z.spct), getTfrType(my.spct))
  expect_equal(getRfrType(my_z.spct), getRfrType(my.spct))
  my_z.spct <- my.spct
  my_z.spct[ , 2] <- 1
  expect_equal(getTfrType(my_z.spct), getTfrType(my.spct))
  expect_equal(getRfrType(my_z.spct), getRfrType(my.spct))
})

