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
  expect_equal(class(my_z.spct[ , -2])[1], "source_spct")
  expect_warning(my_z.spct[ , -2])
#  expect_equal(class(my_z.spct[ FALSE, -2])[1], "source_spct")
#  expect_warning(my_z.spct[ FALSE , -2])

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

context("extract_mspct")

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)

  expect_equal(class(my.mspct)[1], "source_mspct")
  expect_equal(class(my.mspct[1:3])[1], "source_mspct")
  expect_equal(class(my.mspct[-1])[1], "source_mspct")
  expect_equal(class(my.mspct[2])[1], "source_mspct")
  expect_equal(class(my.mspct[[2]])[1], "source_spct")
  expect_equal(class(my.mspct[["spct_2"]])[1], "source_spct")
  expect_equal(class(my.mspct$spct_2)[1], "source_spct")

  expect_equal(my.mspct[[1]], my1.spct)
  expect_equal(my.mspct[[5]], my5.spct)
  expect_equal(my.mspct[["spct_1"]], my1.spct)
  expect_equal(my.mspct[["spct_5"]], my5.spct)

  expect_equal(names(my.mspct["spct_1"]), "spct_1")
  expect_equal(names(my.mspct[1:2]), c("spct_1", "spct_2"))
  expect_equal(names(my.mspct[c(1,5)]), c("spct_1", "spct_5"))

  expect_equal(length(my.mspct[-1]), length(my.mspct) - 1)
  expect_equal(length(my.mspct[-5]), length(my.mspct) - 1)

})

context("replace")

test_that("source_spct", {

  my.spct <- source_spct(w.length = 400:450, s.e.irrad = 0.5, time.unit = "second")
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
  expect_equal(getTimeUnit(my_z.spct), getTimeUnit(my.spct))
  my_z.spct <- my.spct
  my_z.spct[ , 2] <- 1
  expect_equal(getTimeUnit(my_z.spct), getTimeUnit(my.spct))
})

context("replace_mspct")

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)
  my_z.mspct <- my.mspct

  expect_error(my_z.mspct[1] <- 1)
  expect_error(my_z.mspct[1] <- "a")
  expect_error(my_z.mspct[1] <- TRUE)
  expect_error(my_z.mspct[1] <- as.generic_spct(my1.spct))

  expect_error(my_z.mspct[1] <- my1.spct)
  expect_equal(my_z.mspct[1:2] <- my_z.mspct[1:2], my_z.mspct)

  my_z.mspct <- my.mspct

  expect_error(my_z.mspct["spct_1"] <- 1)
  expect_error(my_z.mspct["spct_1"] <- "a")
  expect_error(my_z.mspct["spct_1"] <- TRUE)
  expect_error(my_z.mspct["spct_1"] <- as.generic_spct(my1.spct))

  my_z.mspct <- my.mspct

  expect_error(my_z.mspct[[1]] <- 1)
  expect_error(my_z.mspct[[1]] <- "a")
  expect_error(my_z.mspct[[1]] <- TRUE)
  expect_error(my_z.mspct[[1]] <- as.generic_spct(my1.spct))
  my_z.mspct[["spct_6"]] <- my5.spct
  expect_equal(length(my_z.mspct), length(my.mspct) + 1)
  my_z.mspct[[6]] <- NULL
  expect_equal(my_z.mspct, my.mspct)

  my_z.mspct <- my.mspct

  expect_error(my_z.mspct[["spct_1"]] <- 1)
  expect_error(my_z.mspct[["spct_1"]] <- "a")
  expect_error(my_z.mspct[["spct_1"]] <- TRUE)
  expect_error(my_z.mspct[["spct_1"]] <- as.generic_spct(my1.spct))

})
