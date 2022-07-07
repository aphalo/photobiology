# commented-out tests lead to infinite recursion, but the tested statements
# work as expected at the R command line!
#
library("photobiology")

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
  expect_error(suppressWarnings(my_z.mspct[1] <- as.generic_spct(my1.spct)))

  expect_error(suppressWarnings(my_z.mspct[1] <- my1.spct))
  my_z.mspct[1:2] <- my.mspct[1:2]
  expect_equal(my_z.mspct, my.mspct)

  my_z.mspct <- my.mspct

  expect_error(my_z.mspct["spct_1"] <- 1)
  expect_error(my_z.mspct["spct_1"] <- "a")
  expect_error(my_z.mspct["spct_1"] <- TRUE)
  expect_error(suppressWarnings(my_z.mspct["spct_1"] <- as.generic_spct(my1.spct)))

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

context("combine_mspct")

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)

  spct.l1 <- list(my1.spct, my2.spct)
  spct.l2 <- list(my3.spct, my4.spct, my5.spct)
  my1.mspct <- source_mspct(spct.l1)
  my2.mspct <- source_mspct(spct.l2)
  my12.mspct <- c(my1.mspct, my2.mspct)

  expect_equal(length(my12.mspct), length(my1.mspct) + length(my2.mspct))

  expect_error(c(my1.mspct, my1.spct))

})
