library("photobiology")

context("merge2object_spct")

test_that("default-wl", {

  rfr.spct <- Ler_leaf_rflt.spct
  tfr.spct <- Ler_leaf_trns.spct
  # io
  expect_silent(merge2object_spct(rfr.spct, tfr.spct))
  expect_error(merge2object_spct(tfr.spct, data.frame(a = 1:20, b = 20:1)))
  expect_error(merge2object_spct(rfr.spct, data.frame(a = 1:20, b = 20:1)))
  expect_error(merge2object_spct(data.frame(a = 1:20, b = 20:1), tfr.spct))
  expect_error(merge2object_spct(data.frame(a = 1:20, b = 20:1), rfr.spct))
  # value
  expect_s3_class(merge2object_spct(rfr.spct, tfr.spct), "object_spct")
  expect_false(nrow(merge2object_spct(tfr.spct, rfr.spct)) == nrow(merge2object_spct(rfr.spct, tfr.spct)))
  expect_equal(tfr.spct[["Tfr"]], merge2object_spct(tfr.spct, rfr.spct)[["Tfr"]])
  expect_equal(rfr.spct[["Rfr"]], merge2object_spct(rfr.spct, tfr.spct)[["Rfr"]])
  expect_equal(tfr.spct[["w.length"]], merge2object_spct(tfr.spct, rfr.spct)[["w.length"]])
  expect_equal(rfr.spct[["w.length"]], merge2object_spct(rfr.spct, tfr.spct)[["w.length"]])
  expect_equal(getTfrType(tfr.spct), getTfrType(merge2object_spct(tfr.spct, rfr.spct)))
  expect_equal(getRfrType(rfr.spct), getRfrType(merge2object_spct(tfr.spct, rfr.spct)))
})

test_that("argument-wl", {

  rfr.spct <- Ler_leaf_rflt.spct
  tfr.spct <- Ler_leaf_trns.spct
  # io
  expect_silent(merge2object_spct(rfr.spct, tfr.spct, w.length.out = 380:800))
  expect_error(merge2object_spct(rfr.spct, tfr.spct, w.length.out = 800:300))
  expect_error(merge2object_spct(rfr.spct, tfr.spct, w.length.out =rep(700, 200)))
  expect_error(merge2object_spct(tfr.spct, data.frame(a = 1:20, b = 20:1), w.length.out = 380:800))
  expect_error(merge2object_spct(rfr.spct, data.frame(a = 1:20, b = 20:1), w.length.out = 380:800))
  expect_error(merge2object_spct(data.frame(a = 1:20, b = 20:1), tfr.spct, w.length.out = 380:800))
  expect_error(merge2object_spct(data.frame(a = 1:20, b = 20:1), rfr.spct, w.length.out = 380:800))
  # value
  expect_s3_class(merge2object_spct(rfr.spct, tfr.spct, w.length.out = 380:800), "object_spct")
  # fails because attributes depend on the order of the arguments
  # expect_equal(merge2object_spct(tfr.spct, rfr.spct, w.length.out = 380:800),
  #              merge2object_spct(rfr.spct, tfr.spct, w.length.out = 380:800))
  expect_equal(380:800, merge2object_spct(tfr.spct, rfr.spct, w.length.out = 380:800)[["w.length"]])
  expect_equal(getTfrType(tfr.spct), getTfrType(merge2object_spct(tfr.spct, rfr.spct, w.length.out = 380:800)))
  expect_equal(getRfrType(rfr.spct), getRfrType(merge2object_spct(tfr.spct, rfr.spct, w.length.out = 380:800)))
})


