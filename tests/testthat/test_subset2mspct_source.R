## need to add similar tests for additional classes
## Some attributes get simplified
## The attributes can change when rebuilding the data objects!!

context("subsetting source_spct to collection")

test_that("subset2mspct works with no normalization and IDs", {
  # test shared code only once
  my.spct <- sun_evening.spct

  expect_silent(
    subset2mspct(subset(my.spct,
                        spct.idx == "time.01"),
                 drop.idx = FALSE)
    )
  expect_silent(
    subset2mspct(subset(my.spct,
                        spct.idx %in% c("time.01", "time.02")),
                 drop.idx = FALSE)
  )
  expect_silent(
    subset2mspct(subset(my.spct,
                        spct.idx %in% c("time.01", "time.02", "time.03")),
                 drop.idx = FALSE)
    )
  expect_silent(
    subset2mspct(subset(my.spct,
                        spct.idx %in% c("time.01", "time.02", "time.03", "time.04")),
                 drop.idx = FALSE)
  )

  expect_silent(
    collection.mspct <- subset2mspct(my.spct, drop.idx = FALSE)
  )

  expect_equal(class(my.spct), class(collection.mspct[[1]]))
  expect_equal(class(my.spct), class(collection.mspct[[2]]))
  expect_equal(class(my.spct), class(collection.mspct[[3]]))
  expect_equal(class(my.spct), class(collection.mspct[[4]]))
  expect_equal(class(my.spct), class(collection.mspct[[5]]))

  expect_named(collection.mspct[[1]], colnames(my.spct))
  expect_true(setequal(names(collection.mspct),
                       levels(my.spct[[getIdFactor(my.spct)]])))
  expect_equal(length(collection.mspct), getMultipleWl(my.spct))

  expect_equal(getMultipleWl(collection.mspct[[1]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[2]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[3]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[4]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[5]]), 1)

  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[1]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[2]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[3]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[4]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[5]]))

  expect_equal(names(when_measured(my.spct)), names(collection.mspct))
  expect_equal(when_measured(my.spct)[[1]], when_measured(collection.mspct[[1]]))
  expect_equal(when_measured(my.spct)[[2]], when_measured(collection.mspct[[2]]))
  expect_equal(when_measured(my.spct)[[3]], when_measured(collection.mspct[[3]]))
  expect_equal(when_measured(my.spct)[[4]], when_measured(collection.mspct[[4]]))
  expect_equal(when_measured(my.spct)[[5]], when_measured(collection.mspct[[5]]))

  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[1]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[2]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[3]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[4]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[5]]))

  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[1]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[2]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[3]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[4]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[5]]))

  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[1]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[2]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[3]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[4]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[5]]))

  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[1]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[2]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[3]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[4]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[5]]))

  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[1]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[2]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[3]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[4]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[5]]))

})

test_that("subset2mspct works with no normalization and no IDs", {
  # test shared code only once
  my.spct <- sun_evening.spct

  expect_silent(collection.mspct <- subset2mspct(my.spct, drop.idx = TRUE))

  expect_equal(class(my.spct), class(collection.mspct[[1]]))
  expect_equal(class(my.spct), class(collection.mspct[[2]]))
  expect_equal(class(my.spct), class(collection.mspct[[3]]))
  expect_equal(class(my.spct), class(collection.mspct[[4]]))
  expect_equal(class(my.spct), class(collection.mspct[[5]]))

  expect_named(collection.mspct[[1]], setdiff(colnames(my.spct),
                                              getIdFactor(my.spct)))
  expect_true(setequal(names(collection.mspct),
                       setdiff(levels(my.spct[[getIdFactor(my.spct)]]),
                               getIdFactor(my.spct))))
  expect_equal(length(collection.mspct), getMultipleWl(my.spct))

  expect_equal(getMultipleWl(collection.mspct[[1]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[2]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[3]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[4]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[5]]), 1)

  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[1]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[2]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[3]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[4]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[5]]))

  expect_equal(names(when_measured(my.spct)), names(collection.mspct))
  expect_equal(when_measured(my.spct)[[1]], when_measured(collection.mspct[[1]]))
  expect_equal(when_measured(my.spct)[[2]], when_measured(collection.mspct[[2]]))
  expect_equal(when_measured(my.spct)[[3]], when_measured(collection.mspct[[3]]))
  expect_equal(when_measured(my.spct)[[4]], when_measured(collection.mspct[[4]]))
  expect_equal(when_measured(my.spct)[[5]], when_measured(collection.mspct[[5]]))

  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[1]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[2]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[3]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[4]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[5]]))

  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[1]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[2]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[3]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[4]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[5]]))

  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[1]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[2]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[3]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[4]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[5]]))

  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[1]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[2]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[3]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[4]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[5]]))

  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[1]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[2]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[3]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[4]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[5]]))

  expect_equal(names(getNormalized(my.spct)), names(collection.mspct))
  expect_equal(getNormalized(my.spct)[[1]], getNormalized(collection.mspct[[1]]))
  expect_equal(getNormalized(my.spct)[[2]], getNormalized(collection.mspct[[2]]))
  expect_equal(getNormalized(my.spct)[[3]], getNormalized(collection.mspct[[3]]))
  expect_equal(getNormalized(my.spct)[[4]], getNormalized(collection.mspct[[4]]))
  expect_equal(getNormalized(my.spct)[[5]], getNormalized(collection.mspct[[5]]))

  expect_equal(names(getNormalized(my.spct)), names(collection.mspct))
  expect_equal(getNormalized(my.spct)[[1]], getNormalized(collection.mspct[[1]]))
  expect_equal(getNormalized(my.spct)[[2]], getNormalized(collection.mspct[[2]]))
  expect_equal(getNormalized(my.spct)[[3]], getNormalized(collection.mspct[[3]]))
  expect_equal(getNormalized(my.spct)[[4]], getNormalized(collection.mspct[[4]]))
  expect_equal(getNormalized(my.spct)[[5]], getNormalized(collection.mspct[[5]]))

})

test_that("subset2mspct works with normalization and IDs", {
  # test shared code only once
  my.spct <- normalize(sun_evening.spct)

  expect_silent(collection.mspct <- subset2mspct(my.spct, drop.idx = FALSE))
  expect_named(collection.mspct[[1]], colnames(my.spct))
  expect_true(setequal(names(collection.mspct),
                       levels(my.spct[[getIdFactor(my.spct)]])))
  expect_equal(length(collection.mspct), getMultipleWl(my.spct))

  expect_equal(getMultipleWl(collection.mspct[[1]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[2]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[3]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[4]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[5]]), 1)

  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[1]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[2]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[3]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[4]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[5]]))

  expect_equal(when_measured(my.spct)[[1]], when_measured(collection.mspct[[1]]))
  expect_equal(when_measured(my.spct)[[2]], when_measured(collection.mspct[[2]]))
  expect_equal(when_measured(my.spct)[[3]], when_measured(collection.mspct[[3]]))
  expect_equal(when_measured(my.spct)[[4]], when_measured(collection.mspct[[4]]))
  expect_equal(when_measured(my.spct)[[5]], when_measured(collection.mspct[[5]]))

  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[1]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[2]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[3]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[4]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[5]]))

  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[1]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[2]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[3]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[4]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[5]]))

  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[1]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[2]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[3]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[4]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[5]]))

  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[1]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[2]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[3]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[4]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[5]]))

  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[1]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[2]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[3]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[4]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[5]]))

  expect_equal(names(getNormalized(my.spct)), names(collection.mspct))
  expect_equal(getNormalized(my.spct)[[1]], getNormalized(collection.mspct[[1]]))
  expect_equal(getNormalized(my.spct)[[2]], getNormalized(collection.mspct[[2]]))
  expect_equal(getNormalized(my.spct)[[3]], getNormalized(collection.mspct[[3]]))
  expect_equal(getNormalized(my.spct)[[4]], getNormalized(collection.mspct[[4]]))
  expect_equal(getNormalized(my.spct)[[5]], getNormalized(collection.mspct[[5]]))

  expect_equal(names(getNormalization(my.spct)), names(collection.mspct))
  expect_equal(getNormalization(my.spct)[[1]], getNormalization(collection.mspct[[1]]))
  expect_equal(getNormalization(my.spct)[[2]], getNormalization(collection.mspct[[2]]))
  expect_equal(getNormalization(my.spct)[[3]], getNormalization(collection.mspct[[3]]))
  expect_equal(getNormalization(my.spct)[[4]], getNormalization(collection.mspct[[4]]))
  expect_equal(getNormalization(my.spct)[[5]], getNormalization(collection.mspct[[5]]))

})

test_that("subset2mspct works with normalization and no IDs", {
  # test shared code only once
  my.spct <- normalize(sun_evening.spct)

  expect_silent(collection.mspct <- subset2mspct(my.spct, drop.idx = TRUE))

  expect_equal(class(my.spct), class(collection.mspct[[1]]))
  expect_equal(class(my.spct), class(collection.mspct[[2]]))
  expect_equal(class(my.spct), class(collection.mspct[[3]]))
  expect_equal(class(my.spct), class(collection.mspct[[4]]))
  expect_equal(class(my.spct), class(collection.mspct[[5]]))

  expect_named(collection.mspct[[1]], setdiff(colnames(my.spct),
                                              getIdFactor(my.spct)))
  expect_true(setequal(names(collection.mspct),
                       setdiff(levels(my.spct[[getIdFactor(my.spct)]]),
                               getIdFactor(my.spct))))
  expect_equal(length(collection.mspct), getMultipleWl(my.spct))

  expect_equal(getMultipleWl(collection.mspct[[1]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[2]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[3]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[4]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[5]]), 1)

  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[1]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[2]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[3]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[4]]))
  expect_equal(where_measured(my.spct), where_measured(collection.mspct[[5]]))

  expect_equal(names(when_measured(my.spct)), names(collection.mspct))
  expect_equal(when_measured(my.spct)[[1]], when_measured(collection.mspct[[1]]))
  expect_equal(when_measured(my.spct)[[2]], when_measured(collection.mspct[[2]]))
  expect_equal(when_measured(my.spct)[[3]], when_measured(collection.mspct[[3]]))
  expect_equal(when_measured(my.spct)[[4]], when_measured(collection.mspct[[4]]))
  expect_equal(when_measured(my.spct)[[5]], when_measured(collection.mspct[[5]]))

  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[1]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[2]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[3]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[4]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[5]]))

  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[1]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[2]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[3]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[4]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[5]]))

  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[1]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[2]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[3]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[4]]))
  expect_equal(instr_descriptor(my.spct), instr_descriptor(collection.mspct[[5]]))

  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[1]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[2]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[3]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[4]]))
  expect_equal(instr_settings(my.spct), instr_settings(collection.mspct[[5]]))

  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[1]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[2]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[3]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[4]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[5]]))

  expect_equal(names(getNormalized(my.spct)), names(collection.mspct))
  expect_equal(getNormalized(my.spct)[[1]], getNormalized(collection.mspct[[1]]))
  expect_equal(getNormalized(my.spct)[[2]], getNormalized(collection.mspct[[2]]))
  expect_equal(getNormalized(my.spct)[[3]], getNormalized(collection.mspct[[3]]))
  expect_equal(getNormalized(my.spct)[[4]], getNormalized(collection.mspct[[4]]))
  expect_equal(getNormalized(my.spct)[[5]], getNormalized(collection.mspct[[5]]))

  expect_equal(names(getNormalization(my.spct)), names(collection.mspct))
  expect_equal(getNormalization(my.spct)[[1]], getNormalization(collection.mspct[[1]]))
  expect_equal(getNormalization(my.spct)[[2]], getNormalization(collection.mspct[[2]]))
  expect_equal(getNormalization(my.spct)[[3]], getNormalization(collection.mspct[[3]]))
  expect_equal(getNormalization(my.spct)[[4]], getNormalization(collection.mspct[[4]]))
  expect_equal(getNormalization(my.spct)[[5]], getNormalization(collection.mspct[[5]]))

})

