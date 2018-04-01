library(photobiology)

context("multi_spct.empty")

test_that("constructors", {
  empty.mspct <- response_mspct()
expect_true(is.response_mspct(empty.mspct))
expect_true(is.any_mspct(empty.mspct))
expect_true(is.null(names(empty.mspct)))
expect_length(empty.mspct, 0L)

empty.mspct <- reflector_mspct()
expect_true(is.reflector_mspct(empty.mspct))
expect_true(is.any_mspct(empty.mspct))
expect_true(is.null(names(empty.mspct)))

empty.mspct <- object_mspct()
expect_true(is.object_mspct(empty.mspct))
expect_true(is.any_mspct(empty.mspct))
expect_true(is.null(names(empty.mspct)))
expect_length(empty.mspct, 0L)

empty.mspct <- cps_mspct()
expect_true(is.cps_mspct(empty.mspct))
expect_true(is.any_mspct(empty.mspct))
expect_true(is.null(names(empty.mspct)))

empty.mspct <- raw_mspct()
expect_true(is.raw_mspct(empty.mspct))
expect_true(is.any_mspct(empty.mspct))
expect_true(is.null(names(empty.mspct)))
expect_length(empty.mspct, 0L)

empty.mspct <- calibration_mspct()
expect_true(is.calibration_mspct(empty.mspct))
expect_true(is.any_mspct(empty.mspct))
expect_true(is.null(names(empty.mspct)))
expect_length(empty.mspct, 0L)

empty.mspct <- chroma_mspct()
expect_true(is.chroma_mspct(empty.mspct))
expect_true(is.any_mspct(empty.mspct))
expect_true(is.null(names(empty.mspct)))
expect_length(empty.mspct, 0L)

})
