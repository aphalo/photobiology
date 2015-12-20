library("photobiology")

context("cpp")

test_that("integrate", {

  x <- 100:110
  y <- rep(1, times = length(x))

  expect_equal(integrate_xy(x, y), 10)

})


test_that("insert", {

  x <- 100:110
  y <- rep(1, times = length(x))

  expect_equal(insert_hinges(x, y, numeric()), data.frame(x, y))
  expect_equal(insert_hinges(x, y, 100.1), data.frame(x = sort(c(x, 100.1)), y = 1))

  expect_equal(put_hinges(x, y, numeric()), y)
  expect_equal(put_hinges(x, y, 100.1), rep(1, times = length(x) + 1))
})

