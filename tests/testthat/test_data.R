library("photobiology")

context("example_data_classes")

test_that("any_spct", {
  expect_silent(check_spct(D65.illuminant.spct))
  expect_s3_class(D65.illuminant.spct,
                  c("source_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(A.illuminant.spct))
  expect_s3_class(A.illuminant.spct,
                  c("source_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(sun.spct))
  expect_s3_class(sun.spct,
                  c("source_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(sun.daily.spct))
  expect_s3_class(sun.daily.spct,
                  c("source_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(white_led.cps_spct))
  expect_s3_class(white_led.cps_spct,
                  c("cps_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(white_led.raw_spct))
  expect_s3_class(white_led.raw_spct,
                  c("raw_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(white_led.source_spct))
  expect_s3_class(white_led.source_spct,
                  c("source_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(polyester.spct))
  expect_s3_class(polyester.spct,
                  c("filter_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(yellow_gel.spct))
  expect_s3_class(yellow_gel.spct,
                  c("filter_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(green_leaf.spct))
  expect_s3_class(green_leaf.spct,
                  c("reflector_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(Ler_leaf.spct))
  expect_s3_class(Ler_leaf.spct,
                  c("object_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(Ler_leaf_trns.spct))
  expect_s3_class(Ler_leaf_trns.spct,
                  c("filter_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(Ler_leaf_trns_i.spct))
  expect_s3_class(Ler_leaf_trns_i.spct,
                  c("filter_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(Ler_leaf_rflt.spct))
  expect_s3_class(Ler_leaf_rflt.spct,
                  c("reflector_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(clear.spct))
  expect_s3_class(clear.spct,
                  c("filter_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(opaque.spct))
  expect_s3_class(opaque.spct,
                  c("filter_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(white_body.spct))
  expect_s3_class(white_body.spct,
                  c("object_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(black_body.spct))
  expect_s3_class(black_body.spct,
                  c("object_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(clear_body.spct))
  expect_s3_class(clear_body.spct,
                  c("object_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(ciev10.spct))
  expect_s3_class(ciev10.spct,
                  c("chroma_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(ciev2.spct))
  expect_s3_class(ciev2.spct,
                  c("chroma_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(ciexyzCC10.spct))
  expect_s3_class(ciexyzCC10.spct,
                  c("chroma_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(ciexyzCC2.spct))
  expect_s3_class(ciexyzCC2.spct,
                  c("chroma_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(ciexyzCMF10.spct))
  expect_s3_class(ciexyzCMF10.spct,
                  c("chroma_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(ciexyzCMF2.spct))
  expect_s3_class(ciexyzCMF2.spct,
                  c("chroma_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(beesxyzCMF.spct))
  expect_s3_class(beesxyzCMF.spct,
                  c("chroma_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(ccd.spct))
  expect_s3_class(ccd.spct,
                  c("response_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_silent(check_spct(photodiode.spct))
  expect_s3_class(photodiode.spct,
                  c("response_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
})

context("example_data_numbers")

test_that("source_spct_data", {
  expect_equal(round(sun.spct[c(1, 100, 200, 300), "s.e.irrad"], 7),
               c(0.0000000, 0.4969714, 0.7869773, 0.6573545))
  expect_equal(sun.spct[c(1, 100, 200, 300), "w.length"],
               c(280, 378, 478, 578))
  expect_equal(round(sun.daily.spct[c(1, 100, 200, 300), "s.e.irrad"], 2),
               c(0.00, 19111.46, 31504.73, 26651.81))
  expect_equal(sun.daily.spct[c(1, 100, 200, 300), "w.length"],
               c(280, 378, 478, 578))
  expect_equal(round(white_led.source_spct[c(1, 100, 200, 300, 500), "s.e.irrad"], 7),
               c(0, 0, 0, 0, 0.1140026))
  expect_equal(white_led.source_spct[c(1, 100, 200, 300, 500), "w.length"],
               c(251.16, 298.05, 345.19, 392.09, 485.15))
})

test_that("response_spct_data", {
  expect_equal(round(ccd.spct[c(1, 50, 100, 150), "s.q.response"], 7),
               c(0.6228370, 0.6164968, 0.7302462, 0.3238695))
  expect_equal(round(ccd.spct[c(1, 50, 100, 150), "w.length"], 4),
               c(205.7957, 356.4036, 709.5212, 958.3967))
  expect_equal(round(photodiode.spct[c(1, 25, 50, 75), "s.e.response"], 8),
               c(0.02064991, 0.11857684, 0.18065186, 0.06727758))
  expect_equal(round(photodiode.spct[c(1, 25, 50, 75), "w.length"], 4),
               c(300.0000, 389.5422, 480.6899, 501.9237))
})

test_that("filter_spct_data", {
  expect_equal(round(polyester.spct[c(1, 100, 200, 300), "Tfr"], 4),
               c(0.0048, 0.8547, 0.9247, 0.9180))
  expect_equal(polyester.spct[c(1, 100, 200, 300), "w.length"],
               c(240, 360, 481, 604))
  expect_equal(round(yellow_gel.spct[c(1, 100, 200, 300, 400), "Tfr"], 5),
               c(0.00271, 0.00001, 0.00001, 0.88236, 0.90093))
  expect_equal(yellow_gel.spct[c(1, 100, 200, 300, 400), "w.length"],
               c(190, 291, 443, 566, 761))
})

test_that("object_spct_data", {
  expect_equal(round(Ler_leaf.spct[c(1, 500, 1000, 1500, 2000), "Tfr"], 8),
               c(0.00000000, 0.01057526, 0.03689410, 0.07965146, 0.41962118))
  expect_equal(round(Ler_leaf.spct[c(1, 500, 1000, 1500, 2000), "Rfr"], 8),
               c(0.04668778, 0.03785808, 0.06234386, 0.08139910, 0.56057470))
  expect_equal(Ler_leaf.spct[c(1, 500, 1000, 1500, 2000), "w.length"],
               c(250.00, 374.75, 499.75, 624.75, 749.75))
})
