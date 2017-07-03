library("photobiology")

context("example_data")

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
  expect_warning(check_spct(Ler_leaf.spct))
  expect_s3_class(Ler_leaf.spct,
                  c("object_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_warning(check_spct(Ler_leaf_trns.spct))
  expect_s3_class(Ler_leaf_trns.spct,
                  c("filter_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_warning(check_spct(Ler_leaf_trns_i.spct))
  expect_s3_class(Ler_leaf_trns_i.spct,
                  c("filter_spct", "generic_spct", "tbl_df", "tbl", "data.frame"))
  expect_warning(check_spct(Ler_leaf_rflt.spct))
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

