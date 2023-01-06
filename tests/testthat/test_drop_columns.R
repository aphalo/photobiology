library("photobiology")

context("drop_user_cols")

test_that("source_spct", {

  sun.spct <- sun.spct[1:20, ]
  my.spct <- q2e(sun.spct, action = "add")
  my.spct$A <- "A"
  my.spct$two <- 2

  dropped.spct <- drop_user_cols(my.spct)

  expect_named(my.spct, c("w.length", "s.e.irrad", "s.q.irrad", "A", "two"))
  expect_named(dropped.spct, c("w.length", "s.e.irrad", "s.q.irrad"))
  expect_is(dropped.spct, "source_spct")
  expect_equal(nrow(my.spct), nrow(dropped.spct))
  expect_equal(my.spct$w.length, dropped.spct$w.length)
  expect_equal(my.spct$s.e.irrad, dropped.spct$s.e.irrad)
  expect_equal(my.spct$s.q.irrad, dropped.spct$s.q.irrad)
  # column 1 is "spct.idx"
  expect_equal(spct_metadata(my.spct)[ , -1],
               spct_metadata(dropped.spct)[ , -1])

  # energy
  my.spct <- q2e(sun.spct, action = "replace")
  my.spct$A <- "A"
  my.spct$two <- 2

  dropped.spct <- drop_user_cols(my.spct)

  expect_named(my.spct, c("w.length", "s.e.irrad", "A", "two"))
  expect_named(dropped.spct, c("w.length", "s.e.irrad"))
  expect_is(dropped.spct, "source_spct")
  expect_equal(nrow(my.spct), nrow(dropped.spct))
  expect_equal(my.spct$w.length, dropped.spct$w.length)
  expect_equal(my.spct$s.e.irrad, dropped.spct$s.e.irrad)
  # column 1 is "spct.idx"
  expect_equal(spct_metadata(my.spct)[ , -1],
               spct_metadata(dropped.spct)[ , -1])

  # photon
  my.spct <- e2q(sun.spct, action = "replace")
  my.spct$A <- "A"
  my.spct$two <- 2

  dropped.spct <- drop_user_cols(my.spct)

  expect_named(my.spct, c("w.length", "s.q.irrad", "A", "two"))
  expect_named(dropped.spct, c("w.length", "s.q.irrad"))
  expect_is(dropped.spct, "source_spct")
  expect_equal(nrow(my.spct), nrow(dropped.spct))
  expect_equal(my.spct$w.length, dropped.spct$w.length)
  expect_equal(my.spct$s.q.irrad, dropped.spct$s.q.irrad)
  # column 1 is "spct.idx"
  expect_equal(spct_metadata(my.spct)[ , -1],
               spct_metadata(dropped.spct)[ , -1])

})

test_that("response_spct", {

  ccd.spct <- ccd.spct[1:20, ]
  my.spct <- q2e(ccd.spct, action = "add")
  my.spct$A <- "A"
  my.spct$two <- 2

  dropped.spct <- drop_user_cols(my.spct)

  expect_named(my.spct, c("w.length", "s.q.response", "s.e.response", "A", "two"))
  expect_named(dropped.spct, c("w.length", "s.q.response", "s.e.response"))
  expect_is(dropped.spct, "response_spct")
  expect_equal(nrow(my.spct), nrow(dropped.spct))
  expect_equal(my.spct$w.length, dropped.spct$w.length)
  expect_equal(my.spct$s.e.response, dropped.spct$s.e.response)
  expect_equal(my.spct$s.q.response, dropped.spct$s.q.response)
  # column 1 is "spct.idx"
  expect_equal(spct_metadata(my.spct)[ , -1],
               spct_metadata(dropped.spct)[ , -1])

  # energy
  my.spct <- q2e(ccd.spct, action = "replace")
  my.spct$A <- "A"
  my.spct$two <- 2

  dropped.spct <- drop_user_cols(my.spct)

  expect_named(my.spct, c("w.length", "s.e.response", "A", "two"))
  expect_named(dropped.spct, c("w.length", "s.e.response"))
  expect_is(dropped.spct, "response_spct")
  expect_equal(nrow(my.spct), nrow(dropped.spct))
  expect_equal(my.spct$w.length, dropped.spct$w.length)
  expect_equal(my.spct$s.e.response, dropped.spct$s.e.response)
  # column 1 is "spct.idx"
  expect_equal(spct_metadata(my.spct)[ , -1],
               spct_metadata(dropped.spct)[ , -1])

  # photon
  my.spct <- e2q(ccd.spct, action = "replace")
  my.spct$A <- "A"
  my.spct$two <- 2

  dropped.spct <- drop_user_cols(my.spct)

  expect_named(my.spct, c("w.length", "s.q.response", "A", "two"))
  expect_named(dropped.spct, c("w.length", "s.q.response"))
  expect_is(dropped.spct, "response_spct")
  expect_equal(nrow(my.spct), nrow(dropped.spct))
  expect_equal(my.spct$w.length, dropped.spct$w.length)
  expect_equal(my.spct$s.q.response, dropped.spct$s.q.response)
  # column 1 is "spct.idx"
  expect_equal(spct_metadata(my.spct)[ , -1],
               spct_metadata(dropped.spct)[ , -1])

})

test_that("filter_spct", {
  polyester.spct <- polyester.spct[1:20, ]
  my.spct <- T2A(polyester.spct, action = "add")
  my.spct$Z <- "Z"
  my.spct$two <- 2

  dropped.spct <- drop_user_cols(my.spct)

  expect_named(my.spct, c("w.length", "Tfr", "A", "Z", "two"))
  expect_named(dropped.spct, c("w.length", "Tfr", "A"))
  expect_is(dropped.spct, "filter_spct")
  expect_equal(nrow(my.spct), nrow(dropped.spct))
  expect_equal(my.spct$w.length, dropped.spct$w.length)
  expect_equal(my.spct$Tfr, dropped.spct$Tfr)
  expect_equal(my.spct$A, dropped.spct$A)
  # column 1 is "spct.idx"
  expect_equal(spct_metadata(my.spct)[ , -1],
               spct_metadata(dropped.spct)[ , -1])

  my.spct <- T2A(polyester.spct, action = "replace")
  my.spct$Z <- "Z"
  my.spct$two <- 2

  dropped.spct <- drop_user_cols(my.spct)

  expect_named(my.spct, c("w.length", "A", "Z", "two"))
  expect_named(dropped.spct, c("w.length", "A"))
  expect_is(dropped.spct, "filter_spct")
  expect_equal(nrow(my.spct), nrow(dropped.spct))
  expect_equal(my.spct$w.length, dropped.spct$w.length)
  expect_equal(my.spct$A, dropped.spct$A)
  # column 1 is "spct.idx"
  expect_equal(spct_metadata(my.spct)[ , -1],
               spct_metadata(dropped.spct)[ , -1])

  my.spct <- A2T(polyester.spct, action = "replace")
  my.spct$Z <- "Z"
  my.spct$two <- 2

  dropped.spct <- drop_user_cols(my.spct)

  expect_named(my.spct, c("w.length", "Tfr", "Z", "two"))
  expect_named(dropped.spct, c("w.length", "Tfr"))
  expect_is(dropped.spct, "filter_spct")
  expect_equal(nrow(my.spct), nrow(dropped.spct))
  expect_equal(my.spct$w.length, dropped.spct$w.length)
  expect_equal(my.spct$Tfr, dropped.spct$Tfr)
  # column 1 is "spct.idx"
  expect_equal(spct_metadata(my.spct)[ , -1],
               spct_metadata(dropped.spct)[ , -1])

})

