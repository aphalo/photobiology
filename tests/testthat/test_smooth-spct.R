library("photobiology")

context("smooth_spct")

test_that("source_spct", {

  expect_equal(getTimeUnit(sun.spct), getTimeUnit(smooth_spct(sun.spct)))
  expect_equal(getWhenMeasured(sun.spct), getWhenMeasured(smooth_spct(sun.spct)))
  expect_equal(getWhatMeasured(sun.spct), getWhatMeasured(smooth_spct(sun.spct)))
  expect_equal(getWhereMeasured(sun.spct), getWhereMeasured(smooth_spct(sun.spct)))

  expect_known_output(comment(smooth_spct(sun.spct)), "./data/smooth-comment-output", print = TRUE)

  expect_known_output(smooth_spct(sun.spct), "./data/smooth-default-output", print = TRUE)
  expect_known_output(smooth_spct(sun.spct, method = "custom"), "./data/smooth-custom-output", print = TRUE)
  expect_known_output(smooth_spct(sun.spct, method = "lowess"), "./data/smooth-lowess-output", print = TRUE)
  expect_known_output(smooth_spct(sun.spct, method = "supsmu"), "./data/smooth-supsmu-output", print = TRUE)

  expect_known_output(smooth_spct(sun.spct, strength = 0), "./data/smooth-default-0-output", print = TRUE)
  expect_known_output(smooth_spct(sun.spct, method = "custom", strength = 0), "./data/smooth-custom-0-output", print = TRUE)
  expect_known_output(smooth_spct(sun.spct, method = "lowess", strength = 0.1), "./data/smooth-lowess-0-output", print = TRUE)
  expect_known_output(smooth_spct(sun.spct, method = "supsmu", strength = 0), "./data/smooth-supsmu-0-output", print = TRUE)

  expect_known_output(smooth_spct(sun.spct, strength = 1), "./data/smooth-default-1-output", print = TRUE)
  expect_known_output(smooth_spct(sun.spct, method = "custom", strength = 1), "./data/smooth-custom-1-output", print = TRUE)
  expect_known_output(smooth_spct(sun.spct, method = "lowess", strength = 1), "./data/smooth-lowess-1-output", print = TRUE)
  expect_known_output(smooth_spct(sun.spct, method = "supsmu", strength = 1), "./data/smooth-supsmu-1-output", print = TRUE)

  expect_known_output(smooth_spct(sun.spct, strength = 2), "./data/smooth-default-2-output", print = TRUE)
  expect_known_output(smooth_spct(sun.spct, method = "custom", strength = 2), "./data/smooth-custom-2-output", print = TRUE)
  expect_known_output(smooth_spct(sun.spct, method = "lowess", strength = 2), "./data/smooth-lowess-2-output", print = TRUE)
  expect_known_output(smooth_spct(sun.spct, method = "supsmu", strength = 2), "./data/smooth-supsmu-2-output", print = TRUE)

})

test_that("filter_spct", {

  expect_equal(getTimeUnit(yellow_gel.spct), getTimeUnit(smooth_spct(yellow_gel.spct)))
  expect_equal(getWhenMeasured(yellow_gel.spct), getWhenMeasured(smooth_spct(yellow_gel.spct)))
  expect_equal(getWhatMeasured(yellow_gel.spct), getWhatMeasured(smooth_spct(yellow_gel.spct)))
  expect_equal(getWhereMeasured(yellow_gel.spct), getWhereMeasured(smooth_spct(yellow_gel.spct)))

  expect_known_output(comment(smooth_spct(yellow_gel.spct)), "./data/smooth-flt-comment-output", print = TRUE)

  expect_known_output(smooth_spct(yellow_gel.spct), "./data/smooth-flt-default-output", print = TRUE)
  expect_known_output(smooth_spct(yellow_gel.spct, method = "custom"), "./data/smooth-flt-custom-output", print = TRUE)
  expect_known_output(smooth_spct(yellow_gel.spct, method = "lowess"), "./data/smooth-flt-lowess-output", print = TRUE)
  expect_known_output(smooth_spct(yellow_gel.spct, method = "supsmu"), "./data/smooth-flt-supsmu-output", print = TRUE)

  expect_known_output(smooth_spct(yellow_gel.spct, strength = 0), "./data/smooth-flt-default-0-output", print = TRUE)
  expect_known_output(smooth_spct(yellow_gel.spct, method = "custom", strength = 0), "./data/smooth-flt-custom-0-output", print = TRUE)
  expect_known_output(smooth_spct(yellow_gel.spct, method = "lowess", strength = 0.1), "./data/smooth-flt-lowess-0-output", print = TRUE)
  expect_known_output(smooth_spct(yellow_gel.spct, method = "supsmu", strength = 0), "./data/smooth-flt-supsmu-0-output", print = TRUE)

  expect_known_output(smooth_spct(yellow_gel.spct, strength = 1), "./data/smooth-flt-default-1-output", print = TRUE)
  expect_known_output(smooth_spct(yellow_gel.spct, method = "custom", strength = 1), "./data/smooth-flt-custom-1-output", print = TRUE)
  expect_known_output(smooth_spct(yellow_gel.spct, method = "lowess", strength = 1), "./data/smooth-flt-lowess-1-output", print = TRUE)
  expect_known_output(smooth_spct(yellow_gel.spct, method = "supsmu", strength = 1), "./data/smooth-flt-supsmu-1-output", print = TRUE)

  expect_known_output(smooth_spct(yellow_gel.spct, strength = 2), "./data/smooth-flt-default-2-output", print = TRUE)
  expect_known_output(smooth_spct(yellow_gel.spct, method = "custom", strength = 2), "./data/smooth-flt-custom-2-output", print = TRUE)
  expect_known_output(smooth_spct(yellow_gel.spct, method = "lowess", strength = 2), "./data/smooth-flt-lowess-2-output", print = TRUE)
  expect_known_output(smooth_spct(yellow_gel.spct, method = "supsmu", strength = 2), "./data/smooth-flt-supsmu-2-output", print = TRUE)

})

test_that("reflector_spct", {

  expect_known_output(smooth_spct(green_leaf.spct), "./data/smooth-rflt-default-output", print = TRUE)
  expect_known_output(smooth_spct(green_leaf.spct, method = "custom"), "./data/smooth-rflt-custom-output", print = TRUE)
  expect_known_output(smooth_spct(green_leaf.spct, method = "lowess"), "./data/smooth-rflt-lowess-output", print = TRUE)
  expect_known_output(smooth_spct(green_leaf.spct, method = "supsmu"), "./data/smooth-rflt-supsmu-output", print = TRUE)

  expect_known_output(smooth_spct(green_leaf.spct, strength = 0), "./data/smooth-rflt-default-0-output", print = TRUE)
  expect_known_output(smooth_spct(green_leaf.spct, method = "custom", strength = 0), "./data/smooth-rflt-custom-0-output", print = TRUE)
  expect_known_output(smooth_spct(green_leaf.spct, method = "lowess", strength = 0.1), "./data/smooth-rflt-lowess-0-output", print = TRUE)
  expect_known_output(smooth_spct(green_leaf.spct, method = "supsmu", strength = 0), "./data/smooth-rflt-supsmu-0-output", print = TRUE)

  expect_known_output(smooth_spct(green_leaf.spct, strength = 1), "./data/smooth-rflt-default-1-output", print = TRUE)
  expect_known_output(smooth_spct(green_leaf.spct, method = "custom", strength = 1), "./data/smooth-rflt-custom-1-output", print = TRUE)
  expect_known_output(smooth_spct(green_leaf.spct, method = "lowess", strength = 1), "./data/smooth-rflt-lowess-1-output", print = TRUE)
  expect_known_output(smooth_spct(green_leaf.spct, method = "supsmu", strength = 1), "./data/smooth-rflt-supsmu-1-output", print = TRUE)

  expect_known_output(smooth_spct(green_leaf.spct, strength = 2), "./data/smooth-rflt-default-2-output", print = TRUE)
  expect_known_output(smooth_spct(green_leaf.spct, method = "custom", strength = 2), "./data/smooth-rflt-custom-2-output", print = TRUE)
  expect_known_output(smooth_spct(green_leaf.spct, method = "lowess", strength = 2), "./data/smooth-rflt-lowess-2-output", print = TRUE)
  expect_known_output(smooth_spct(green_leaf.spct, method = "supsmu", strength = 2), "./data/smooth-rflt-supsmu-2-output", print = TRUE)

})

test_that("response_spct", {

  expect_known_output(smooth_spct(ccd.spct), "./data/smooth-rsp-default-output", print = TRUE)
  expect_known_output(smooth_spct(ccd.spct, method = "custom"), "./data/smooth-rsp-custom-output", print = TRUE)
  expect_known_output(smooth_spct(ccd.spct, method = "lowess"), "./data/smooth-rsp-lowess-output", print = TRUE)
  expect_known_output(smooth_spct(ccd.spct, method = "supsmu"), "./data/smooth-rsp-supsmu-output", print = TRUE)

  expect_known_output(smooth_spct(ccd.spct, strength = 0), "./data/smooth-rsp-default-0-output", print = TRUE)
  expect_known_output(smooth_spct(ccd.spct, method = "custom", strength = 0), "./data/smooth-rsp-custom-0-output", print = TRUE)
  expect_known_output(smooth_spct(ccd.spct, method = "lowess", strength = 0.1), "./data/smooth-rsp-lowess-0-output", print = TRUE)
  expect_known_output(smooth_spct(ccd.spct, method = "supsmu", strength = 0), "./data/smooth-rsp-supsmu-0-output", print = TRUE)

  expect_known_output(smooth_spct(ccd.spct, strength = 1), "./data/smooth-rsp-default-1-output", print = TRUE)
  expect_known_output(smooth_spct(ccd.spct, method = "custom", strength = 1), "./data/smooth-rsp-custom-1-output", print = TRUE)
  expect_known_output(smooth_spct(ccd.spct, method = "lowess", strength = 1), "./data/smooth-rsp-lowess-1-output", print = TRUE)
  expect_known_output(smooth_spct(ccd.spct, method = "supsmu", strength = 1), "./data/smooth-rsp-supsmu-1-output", print = TRUE)

  expect_known_output(smooth_spct(ccd.spct, strength = 2), "./data/smooth-rsp-default-2-output", print = TRUE)
  expect_known_output(smooth_spct(ccd.spct, method = "custom", strength = 2), "./data/smooth-rsp-custom-2-output", print = TRUE)
  expect_known_output(smooth_spct(ccd.spct, method = "lowess", strength = 2), "./data/smooth-rsp-lowess-2-output", print = TRUE)
  expect_known_output(smooth_spct(ccd.spct, method = "supsmu", strength = 2), "./data/smooth-rsp-supsmu-2-output", print = TRUE)

})

