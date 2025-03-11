library("photobiology")

update_all <- FALSE

context("s_mean")

test_that("source_mspct", {

  load("./data/test-lamps-mspct.Rda")

  row100 <- c(lamps.mspct[[1]][[2]][100], lamps.mspct[[2]][[2]][100], lamps.mspct[[3]][[2]][100])
  expect_equal(s_mean(lamps.mspct)[[2]][100], mean(row100))
  expect_equal(s_median(lamps.mspct)[[2]][100], median(row100))
  expect_equal(s_var(lamps.mspct)[[2]][100], var(row100))
  expect_equal(s_sd(lamps.mspct)[[2]][100], sd(row100))
  expect_equal(s_sum(lamps.mspct)[[2]][100], sum(row100))
#  expect_equal(s_prod(lamps.mspct)[[2]][100], prod(row100))

  print.flag <- FALSE # TRUE leads to warnings about encoding
  expect_known_value(s_mean(lamps.mspct), "data/s-mean-value", update = update_all)
  expect_known_value(s_median(lamps.mspct), "data/s-median-value", update = update_all)
  expect_known_value(s_mean_se(lamps.mspct), "data/s-mean-se-value", update = update_all)
  expect_known_value(s_sd(lamps.mspct), "data/s-sd-value", update = update_all)
  expect_known_value(s_se(lamps.mspct), "data/s-se-value", update = update_all)
  expect_known_value(s_var(lamps.mspct), "data/s-var-value", update = update_all)
  expect_known_value(s_sum(lamps.mspct), "data/s-sum-value", update = update_all)
  expect_warning(s_prod(lamps.mspct))

  expect_known_value(s_mean(lamps.mspct, na.rm = TRUE), "data/s-mean-na-value", update = update_all)
  expect_known_value(s_median(lamps.mspct, na.rm = TRUE), "data/s-median-na-value", update = update_all)
  expect_known_value(s_mean_se(lamps.mspct, na.rm = TRUE), "data/s-mean-se-na-value", update = update_all)
  expect_known_value(s_sd(lamps.mspct, na.rm = TRUE), "data/s-sd-na-value", update = update_all)
  expect_known_value(s_se(lamps.mspct, na.rm = TRUE), "data/s-se-na-value", update = update_all)
  expect_known_value(s_var(lamps.mspct, na.rm = TRUE), "data/s-var-na-value", update = update_all)
  expect_known_value(s_sum(lamps.mspct, na.rm = TRUE), "data/s-sum-na-value", update = update_all)
  expect_warning(s_prod(lamps.mspct, na.rm = TRUE))
})

