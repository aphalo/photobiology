library("photobiology")
# need to ad cases for solute_spct

context("smooth_spct")

update_all <- FALSE

test_that("smoothing of long source_spct", {

  my1.spct <- my2.spct <- my3.spct <- my4.spct <- my5.spct <- sun.spct[200:300]

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)
  my.spct <- rbindspct(my.mspct, idfactor = "test.id")

  smoothed.spct <- smooth_spct(my.spct)

  expect_equal(getIdFactor(my.spct), getIdFactor(smoothed.spct))
  expect_equal(sort(colnames(my.spct)), sort(colnames(smoothed.spct)))
  expect_equal(getMultipleWl(my.spct), getMultipleWl(smoothed.spct))
  expect_equal(nrow(my.spct), nrow(smoothed.spct))
  expect_equal(my.spct$w.length, smoothed.spct$w.length)
  expect_equal(class(my.spct), class(smoothed.spct))

})


test_that("source_spct", {

  x <- sun.spct[200:300]
  x$s.e.irrad[10:20] <- NA_real_
  expect_error(smooth_spct(x))

  expect_equal(smooth_spct(sun.spct, method = "skip"), sun.spct)

  expect_known_value(smooth_spct(x, na.rm = TRUE), "./data/smooth-na-value", update = update_all)

  expect_equal(getTimeUnit(sun.spct), getTimeUnit(smooth_spct(sun.spct)))
  expect_equal(getWhenMeasured(sun.spct), getWhenMeasured(smooth_spct(sun.spct)))
  expect_equal(getWhatMeasured(sun.spct), getWhatMeasured(smooth_spct(sun.spct)))
  expect_equal(getWhereMeasured(sun.spct), getWhereMeasured(smooth_spct(sun.spct)))

  expect_known_value(comment(smooth_spct(sun.spct)), "./data/smooth-comment-value", update = update_all)

  expect_known_value(smooth_spct(sun.spct), "./data/smooth-default-value", update = update_all)
  expect_known_value(smooth_spct(sun.spct, method = "custom"), "./data/smooth-custom-value", update = update_all)
  expect_known_value(smooth_spct(sun.spct, method = "lowess"), "./data/smooth-lowess-value", update = update_all)
  expect_known_value(smooth_spct(sun.spct, method = "supsmu"), "./data/smooth-supsmu-value", update = update_all)

  expect_known_value(smooth_spct(sun.spct, strength = 0), "./data/smooth-default-0-value", update = update_all)
  expect_known_value(smooth_spct(sun.spct, method = "custom", strength = 0), "./data/smooth-custom-0-value", update = update_all)
  expect_known_value(smooth_spct(sun.spct, method = "lowess", strength = 0), "./data/smooth-lowess-0-value", update = update_all)
  expect_known_value(smooth_spct(sun.spct, method = "supsmu", strength = 0), "./data/smooth-supsmu-0-value", update = update_all)

  expect_known_value(smooth_spct(sun.spct, strength = 1), "./data/smooth-default-1-value", update = update_all)
  expect_known_value(smooth_spct(sun.spct, method = "custom", strength = 1), "./data/smooth-custom-1-value", update = update_all)
  expect_known_value(smooth_spct(sun.spct, method = "lowess", strength = 1), "./data/smooth-lowess-1-value", update = update_all)
  expect_known_value(smooth_spct(sun.spct, method = "supsmu", strength = 1), "./data/smooth-supsmu-1-value", update = update_all)

  expect_known_value(smooth_spct(sun.spct, strength = 2), "./data/smooth-default-2-value", update = update_all)
  expect_known_value(smooth_spct(sun.spct, method = "custom", strength = 2), "./data/smooth-custom-2-value", update = update_all)
  expect_known_value(smooth_spct(sun.spct, method = "lowess", strength = 2), "./data/smooth-lowess-2-value", update = update_all)
  expect_known_value(smooth_spct(sun.spct, method = "supsmu", strength = 2), "./data/smooth-supsmu-2-value", update = update_all)

})

test_that("filter_spct", {

  x <- yellow_gel.spct[200:400]
  x$Tfr[10:20] <- NA_real_
  expect_error(smooth_spct(x))

  expect_known_value(smooth_spct(x, na.rm = TRUE), "./data/smooth-flt-na-value", update = update_all)

  expect_equal(getTimeUnit(yellow_gel.spct), getTimeUnit(smooth_spct(yellow_gel.spct)))
  expect_equal(getWhenMeasured(yellow_gel.spct), getWhenMeasured(smooth_spct(yellow_gel.spct)))
  expect_equal(getWhatMeasured(yellow_gel.spct), getWhatMeasured(smooth_spct(yellow_gel.spct)))
  expect_equal(getWhereMeasured(yellow_gel.spct), getWhereMeasured(smooth_spct(yellow_gel.spct)))

  expect_known_value(comment(smooth_spct(yellow_gel.spct)), "./data/smooth-flt-comment-value", update = update_all)

  expect_known_value(smooth_spct(yellow_gel.spct), "./data/smooth-flt-default-value", update = update_all)
  expect_known_value(smooth_spct(yellow_gel.spct, method = "custom"), "./data/smooth-flt-custom-value", update = update_all)
  expect_known_value(smooth_spct(yellow_gel.spct, method = "lowess"), "./data/smooth-flt-lowess-value", update = update_all)
  expect_known_value(smooth_spct(yellow_gel.spct, method = "supsmu"), "./data/smooth-flt-supsmu-value", update = update_all)

  expect_known_value(smooth_spct(yellow_gel.spct, strength = 0), "./data/smooth-flt-default-0-value", update = update_all)
  expect_known_value(smooth_spct(yellow_gel.spct, method = "custom", strength = 0), "./data/smooth-flt-custom-0-value", update = update_all)
  expect_known_value(smooth_spct(yellow_gel.spct, method = "lowess", strength = 0.1), "./data/smooth-flt-lowess-0-value", update = update_all)
  expect_known_value(smooth_spct(yellow_gel.spct, method = "supsmu", strength = 0), "./data/smooth-flt-supsmu-0-value", update = update_all)

  expect_known_value(smooth_spct(yellow_gel.spct, strength = 1), "./data/smooth-flt-default-1-value", update = update_all)
  expect_known_value(smooth_spct(yellow_gel.spct, method = "custom", strength = 1), "./data/smooth-flt-custom-1-value", update = update_all)
  expect_known_value(smooth_spct(yellow_gel.spct, method = "lowess", strength = 1), "./data/smooth-flt-lowess-1-value", update = update_all)
  expect_known_value(smooth_spct(yellow_gel.spct, method = "supsmu", strength = 1), "./data/smooth-flt-supsmu-1-value", update = update_all)

  expect_known_value(smooth_spct(yellow_gel.spct, strength = 2), "./data/smooth-flt-default-2-value", update = update_all)
  expect_known_value(smooth_spct(yellow_gel.spct, method = "custom", strength = 2), "./data/smooth-flt-custom-2-value", update = update_all)
  expect_known_value(smooth_spct(yellow_gel.spct, method = "lowess", strength = 2), "./data/smooth-flt-lowess-2-value", update = update_all)
  expect_known_value(smooth_spct(yellow_gel.spct, method = "supsmu", strength = 2), "./data/smooth-flt-supsmu-2-value", update = update_all)

})

test_that("reflector_spct", {

  x <- green_leaf.spct
  x$Rfr[10:20] <- NA_real_
  expect_error(smooth_spct(x))

  expect_known_value(smooth_spct(x, na.rm = TRUE), "./data/smooth-rflt-na-value", update = update_all)

  expect_equal(getTimeUnit(green_leaf.spct), getTimeUnit(smooth_spct(green_leaf.spct)))
  expect_equal(getWhenMeasured(green_leaf.spct), getWhenMeasured(smooth_spct(green_leaf.spct)))
  expect_equal(getWhatMeasured(green_leaf.spct), getWhatMeasured(smooth_spct(green_leaf.spct)))
  expect_equal(getWhereMeasured(green_leaf.spct), getWhereMeasured(smooth_spct(green_leaf.spct)))

  expect_known_value(comment(smooth_spct(green_leaf.spct)), "./data/smooth-rflt-comment-value", update = update_all)

  expect_known_value(smooth_spct(green_leaf.spct), "./data/smooth-rflt-default-value", update = update_all)
  expect_known_value(smooth_spct(green_leaf.spct, method = "custom"), "./data/smooth-rflt-custom-value", update = update_all)
  expect_known_value(smooth_spct(green_leaf.spct, method = "lowess"), "./data/smooth-rflt-lowess-value", update = update_all)
  expect_known_value(smooth_spct(green_leaf.spct, method = "supsmu"), "./data/smooth-rflt-supsmu-value", update = update_all)

  expect_known_value(smooth_spct(green_leaf.spct, strength = 0), "./data/smooth-rflt-default-0-value", update = update_all)
  expect_known_value(smooth_spct(green_leaf.spct, method = "custom", strength = 0), "./data/smooth-rflt-custom-0-value", update = update_all)
  expect_known_value(smooth_spct(green_leaf.spct, method = "lowess", strength = 0.1), "./data/smooth-rflt-lowess-0-value", update = update_all)
  expect_known_value(smooth_spct(green_leaf.spct, method = "supsmu", strength = 0), "./data/smooth-rflt-supsmu-0-value", update = update_all)

  expect_known_value(smooth_spct(green_leaf.spct, strength = 1), "./data/smooth-rflt-default-1-value", update = update_all)
  expect_known_value(smooth_spct(green_leaf.spct, method = "custom", strength = 1), "./data/smooth-rflt-custom-1-value", update = update_all)
  expect_known_value(smooth_spct(green_leaf.spct, method = "lowess", strength = 1), "./data/smooth-rflt-lowess-1-value", update = update_all)
  expect_known_value(smooth_spct(green_leaf.spct, method = "supsmu", strength = 1), "./data/smooth-rflt-supsmu-1-value", update = update_all)

  expect_known_value(smooth_spct(green_leaf.spct, strength = 2), "./data/smooth-rflt-default-2-value", update = update_all)
  expect_known_value(smooth_spct(green_leaf.spct, method = "custom", strength = 2), "./data/smooth-rflt-custom-2-value", update = update_all)
  expect_known_value(smooth_spct(green_leaf.spct, method = "lowess", strength = 2), "./data/smooth-rflt-lowess-2-value", update = update_all)
  expect_known_value(smooth_spct(green_leaf.spct, method = "supsmu", strength = 2), "./data/smooth-rflt-supsmu-2-value", update = update_all)

})

test_that("response_spct", {

  x <- ccd.spct
  x$s.q.response[10:20] <- NA_real_
  expect_error(smooth_spct(x))

  expect_known_value(smooth_spct(x, na.rm = TRUE), "./data/smooth-rsp-na-value", update = update_all)

  expect_equal(getTimeUnit(ccd.spct), getTimeUnit(smooth_spct(ccd.spct)))
  expect_equal(getWhenMeasured(ccd.spct), getWhenMeasured(smooth_spct(ccd.spct)))
  expect_equal(getWhatMeasured(ccd.spct), getWhatMeasured(smooth_spct(ccd.spct)))
  expect_equal(getWhereMeasured(ccd.spct), getWhereMeasured(smooth_spct(ccd.spct)))

  expect_known_value(comment(smooth_spct(ccd.spct)), "./data/smooth-rsp-comment-value", update = update_all)

  expect_known_value(smooth_spct(ccd.spct), "./data/smooth-rsp-default-value", update = update_all)
  expect_known_value(smooth_spct(ccd.spct, method = "custom"), "./data/smooth-rsp-custom-value", update = update_all)
  expect_known_value(smooth_spct(ccd.spct, method = "lowess"), "./data/smooth-rsp-lowess-value", update = update_all)
  expect_known_value(smooth_spct(ccd.spct, method = "supsmu"), "./data/smooth-rsp-supsmu-value", update = update_all)

  expect_known_value(smooth_spct(ccd.spct, strength = 0), "./data/smooth-rsp-default-0-value", update = update_all)
  expect_known_value(smooth_spct(ccd.spct, method = "custom", strength = 0), "./data/smooth-rsp-custom-0-value", update = update_all)
  expect_known_value(smooth_spct(ccd.spct, method = "lowess", strength = 0.1), "./data/smooth-rsp-lowess-0-value", update = update_all)
  expect_known_value(smooth_spct(ccd.spct, method = "supsmu", strength = 0), "./data/smooth-rsp-supsmu-0-value", update = update_all)

  expect_known_value(smooth_spct(ccd.spct, strength = 1), "./data/smooth-rsp-default-1-value", update = update_all)
  expect_known_value(smooth_spct(ccd.spct, method = "custom", strength = 1), "./data/smooth-rsp-custom-1-value", update = update_all)
  expect_known_value(smooth_spct(ccd.spct, method = "lowess", strength = 1), "./data/smooth-rsp-lowess-1-value", update = update_all)
  expect_known_value(smooth_spct(ccd.spct, method = "supsmu", strength = 1), "./data/smooth-rsp-supsmu-1-value", update = update_all)

  expect_known_value(smooth_spct(ccd.spct, strength = 2), "./data/smooth-rsp-default-2-value", update = update_all)
  expect_known_value(smooth_spct(ccd.spct, method = "custom", strength = 2), "./data/smooth-rsp-custom-2-value", update = update_all)
  expect_known_value(smooth_spct(ccd.spct, method = "lowess", strength = 2), "./data/smooth-rsp-lowess-2-value", update = update_all)
  expect_known_value(smooth_spct(ccd.spct, method = "supsmu", strength = 2), "./data/smooth-rsp-supsmu-2-value", update = update_all)

})

