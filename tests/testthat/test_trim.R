library("photobiology")

context("trim_spct")

test_that("source_spct", {

  my.spct <- source_spct(w.length = 400:450, s.e.irrad = 0.5, time.unit = "second")

  expect_equal(length(trim_spct(my.spct)), length(my.spct))
  expect_equal(min(trim_spct(my.spct)), min(my.spct))
  expect_equal(max(trim_spct(my.spct)), max(my.spct))
  expect_equal(length(trim_spct(my.spct)[["w.length"]]), length(my.spct[["w.length"]]))
  expect_equal(length(trim_spct(my.spct, range = c(NA, NA))), length(my.spct))

  expect_equal(max(trim_spct(my.spct, low.limit = 401)), max(my.spct))
  expect_equal(min(trim_spct(my.spct, high.limit = 449)), min(my.spct))
  expect_equal(max(trim_spct(my.spct, range = c(401, NA))), max(my.spct))
  expect_equal(min(trim_spct(my.spct, range = c(NA, 449))), min(my.spct))

  expect_equal(max(trim_spct(my.spct, low.limit = 401, fill = NA)), max(my.spct))
  expect_equal(min(trim_spct(my.spct, high.limit = 449, fill = NA)), min(my.spct))
  expect_equal(max(trim_spct(my.spct, low.limit = 401, fill = NULL)), max(my.spct))
  expect_equal(min(trim_spct(my.spct, high.limit = 449, fill = NULL)), min(my.spct))

  expect_equal(max(trim_spct(my.spct, low.limit = 401.1, fill = NA)), max(my.spct))
  expect_equal(min(trim_spct(my.spct, high.limit = 449.1, fill = NA)), min(my.spct))
  expect_equal(max(trim_spct(my.spct, low.limit = 401.1, fill = NULL)), max(my.spct))
  expect_equal(min(trim_spct(my.spct, high.limit = 449.1, fill = NULL)), min(my.spct))

  expect_equal(min(trim_spct(my.spct, low.limit = 401.1, fill = NA)), min(my.spct))
  expect_equal(max(trim_spct(my.spct, high.limit = 449.1, fill = NA)), max(my.spct))
  expect_equal(min(trim_spct(my.spct, low.limit = 401.1, fill = NULL)), 401.1)
  expect_equal(max(trim_spct(my.spct, high.limit = 449.1, fill = NULL)), 449.1)

  my_z.spct <- trim_spct(my.spct, range = c(min(my.spct) + 0.2, max(my.spct) - 0.2))

  expect_equal(class(my_z.spct)[1], "source_spct")
  expect_equal(getTimeUnit(my_z.spct), getTimeUnit(my.spct))
  expect_equal(getBSWFUsed(my_z.spct), getBSWFUsed(my.spct))
  expect_equal(is_effective(my_z.spct), is_effective(my.spct))
  expect_equal(is_normalized(my_z.spct), is_normalized(my.spct))
  expect_equal(getNormalized(my_z.spct), getNormalized(my.spct))
  expect_equal(getScaled(my_z.spct), getScaled(my.spct))

  my_z.spct <- trim_spct(normalize(my.spct), range = c(min(my.spct) + 0.2, max(my.spct) - 0.2))

  expect_equal(class(my_z.spct)[1], "source_spct")
  expect_equal(getTimeUnit(my_z.spct), getTimeUnit(my.spct))
  expect_equal(getBSWFUsed(my_z.spct), getBSWFUsed(my.spct))
  expect_equal(is_effective(my_z.spct), is_effective(my.spct))
  expect_equal(is_normalized(my_z.spct), is_normalized(my.spct))
  expect_equal(getNormalized(my_z.spct), getNormalized(my.spct))
  expect_equal(getScaled(my_z.spct), getScaled(my.spct))

  my_z.spct <- trim_spct(fscale(my.spct), range = c(min(my.spct) + 0.2, max(my.spct) - 0.2))

  expect_equal(class(my_z.spct)[1], "source_spct")
  expect_equal(getTimeUnit(my_z.spct), getTimeUnit(my.spct))
  expect_equal(getBSWFUsed(my_z.spct), getBSWFUsed(my.spct))
  expect_equal(is_effective(my_z.spct), is_effective(my.spct))
  expect_equal(is_normalized(my_z.spct), is_normalized(my.spct))
  expect_equal(getNormalized(my_z.spct), getNormalized(my.spct))
  expect_equal(getScaled(my_z.spct), getScaled(my.spct))

  my_z.spct <- trim_spct(setTimeUnit(my.spct, time.unit = "day"), range = c(min(my.spct) + 0.2, max(my.spct) - 0.2))

})

