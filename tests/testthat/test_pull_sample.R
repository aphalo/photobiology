library("photobiology")

context("pull sample")

test_that("source_mspct", {

  set.seed(987654321)
  sampled.mspct <- pull_sample(sun_evening.mspct, size = 1)
  expect_is(sampled.mspct, "source_mspct")
  expect_equal(length(sampled.mspct), 1)
  expect_named(sampled.mspct, "time.05")
  expect_equal(sampled.mspct, sun_evening.mspct["time.05"])

  set.seed(987654321)
  sampled.spct <- pull_sample(sun_evening.mspct, size = 1, simplify = TRUE)
  expect_is(sampled.spct, "source_spct")
  expect_equal(nrow(sampled.spct), nrow(sun_evening.mspct[["time.05"]]))
  expect_equal(sampled.spct, sun_evening.mspct[["time.05"]])

  set.seed(987654321)
  sampled.mspct <- pull_sample(sun_evening.mspct, size = 2)
  expect_is(sampled.mspct, "source_mspct")
  expect_equal(length(sampled.mspct), 2)
  expect_named(sampled.mspct, c("time.04", "time.05"))
  expect_equal(sampled.mspct, sun_evening.mspct[c("time.04", "time.05")])

})

test_that("source_spct", {

  set.seed(987654321)
  sampled.spct <- pull_sample(sun_evening.spct, size = 1)
  expect_is(sampled.spct, "source_spct")
  expect_equal(getMultipleWl(sampled.spct), 1)
  expect_equal(as.character(unique(sampled.spct[["spct.idx"]])), "time.05")

  set.seed(987654321)
  sampled.spct <- pull_sample(sun_evening.spct, size = 2)
  expect_is(sampled.spct, "source_spct")
  expect_equal(getMultipleWl(sampled.spct), 2)
  expect_true(all(as.character(unique(sampled.spct[["spct.idx"]])) %in% c("time.04", "time.05")))

})

test_that("filter_mspct", {

  set.seed(987654321)
  sampled.mspct <- pull_sample(two_filters.mspct, size = 1)
  expect_is(sampled.mspct, "filter_mspct")
  expect_equal(length(sampled.mspct), 1)
  expect_equal(sampled.mspct, two_filters.mspct[1])

  set.seed(987654321)
  sampled.spct <- pull_sample(two_filters.mspct, size = 1, simplify = TRUE)
  expect_is(sampled.spct, "filter_spct")
  expect_equal(sampled.spct, two_filters.mspct[[1]])

  set.seed(987654321)
  sampled.mspct <- pull_sample(two_filters.mspct, size = 2)
  expect_is(sampled.mspct, "filter_mspct")
  expect_equal(length(sampled.mspct), 2)
  expect_equal(sampled.mspct, two_filters.mspct[1:2])

})

test_that("filter_spct", {

  set.seed(987654321)
  sampled.spct <- pull_sample(two_filters.spct, size = 1)
  expect_is(sampled.spct, "filter_spct")
  expect_equal(getMultipleWl(sampled.spct), 1)

  set.seed(987654321)
  sampled.spct <- pull_sample(two_filters.spct, size = 2)
  expect_is(sampled.spct, "filter_spct")
  expect_equal(getMultipleWl(sampled.spct), 2)

})

