## need to add similar tests for additional classes
## Some attributes get simplified
## The attributes can change when rebuilding the data objects!!

context("subsetting filter_spct to collection")

test_that("subset2mspct works with no normalization and IDs", {
  # test shared code only once
  my.spct <- two_filters.spct

  expect_silent(collection.mspct <- subset2mspct(my.spct, drop.idx = FALSE))

  expect_equal(class(my.spct), class(collection.mspct[[1]]))
  expect_equal(class(my.spct), class(collection.mspct[[2]]))

  expect_named(collection.mspct[[1]], colnames(my.spct))
  expect_true(setequal(names(collection.mspct),
                       levels(my.spct[[getIdFactor(my.spct)]])))
  expect_equal(length(collection.mspct), getMultipleWl(my.spct))

  expect_equal(getMultipleWl(collection.mspct[[1]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[2]]), 1)

  expect_equal(names(where_measured(my.spct, .bind.geocodes = FALSE)),
               names(collection.mspct))
  expect_equal(where_measured(my.spct, .bind.geocodes = FALSE)[[1]],
               where_measured(collection.mspct[[1]]))
  expect_equal(where_measured(my.spct, .bind.geocodes = FALSE)[[2]],
               where_measured(collection.mspct[[2]]))

  expect_equal(names(when_measured(my.spct)), names(collection.mspct))
  expect_equal(when_measured(my.spct)[[1]], when_measured(collection.mspct[[1]]))
  expect_equal(when_measured(my.spct)[[2]], when_measured(collection.mspct[[2]]))

  expect_equal(names(what_measured(my.spct)), names(collection.mspct))
  expect_equal(what_measured(my.spct)[[1]], what_measured(collection.mspct[[1]]))
  expect_equal(what_measured(my.spct)[[2]], what_measured(collection.mspct[[2]]))

  expect_equal(names(how_measured(my.spct)), names(collection.mspct))
  expect_equal(how_measured(my.spct)[[1]], how_measured(collection.mspct[[1]]))
  expect_equal(how_measured(my.spct)[[2]], how_measured(collection.mspct[[2]]))

  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[1]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[2]]))

})

test_that("subset2mspct works with no normalization and no IDs", {
  # test shared code only once
  my.spct <- two_filters.spct

  expect_silent(collection.mspct <- subset2mspct(my.spct, drop.idx = TRUE))

  expect_equal(class(my.spct), class(collection.mspct[[1]]))
  expect_equal(class(my.spct), class(collection.mspct[[2]]))

  expect_named(collection.mspct[[1]], setdiff(colnames(my.spct),
                                              getIdFactor(my.spct)))
  expect_true(setequal(names(collection.mspct),
                       setdiff(levels(my.spct[[getIdFactor(my.spct)]]),
                               getIdFactor(my.spct))))
  expect_equal(length(collection.mspct), getMultipleWl(my.spct))

  expect_equal(getMultipleWl(collection.mspct[[1]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[2]]), 1)

  expect_equal(names(where_measured(my.spct, .bind.geocodes = FALSE)),
               names(collection.mspct))
  expect_equal(where_measured(my.spct, .bind.geocodes = FALSE)[[1]],
               where_measured(collection.mspct[[1]]))
  expect_equal(where_measured(my.spct, .bind.geocodes = FALSE)[[2]],
               where_measured(collection.mspct[[2]]))

  expect_equal(names(when_measured(my.spct)), names(collection.mspct))
  expect_equal(when_measured(my.spct)[[1]], when_measured(collection.mspct[[1]]))
  expect_equal(when_measured(my.spct)[[2]], when_measured(collection.mspct[[2]]))

  expect_equal(names(what_measured(my.spct)), names(collection.mspct))
  expect_equal(what_measured(my.spct)[[1]], what_measured(collection.mspct[[1]]))
  expect_equal(what_measured(my.spct)[[2]], what_measured(collection.mspct[[2]]))

  expect_equal(names(how_measured(my.spct)), names(collection.mspct))
  expect_equal(how_measured(my.spct)[[1]], how_measured(collection.mspct[[1]]))
  expect_equal(how_measured(my.spct)[[2]], how_measured(collection.mspct[[2]]))

  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[1]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[2]]))
})

test_that("subset2mspct works with normalization and IDs", {
  # test shared code only once
  my.spct <- normalize(two_filters.spct)

  expect_silent(collection.mspct <- subset2mspct(my.spct, drop.idx = FALSE))

  expect_equal(class(my.spct), class(collection.mspct[[1]]))
  expect_equal(class(my.spct), class(collection.mspct[[2]]))

  expect_named(collection.mspct[[1]], colnames(my.spct))
  expect_true(setequal(names(collection.mspct),
                       levels(my.spct[[getIdFactor(my.spct)]])))
  expect_equal(length(collection.mspct), getMultipleWl(my.spct))

  expect_equal(getMultipleWl(collection.mspct[[1]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[2]]), 1)

  expect_equal(where_measured(my.spct),
               where_measured(collection.mspct[[1]]))
  expect_equal(where_measured(my.spct),
               where_measured(collection.mspct[[2]]))

  expect_equal(names(when_measured(my.spct)), names(collection.mspct))
  expect_equal(when_measured(my.spct)[[1]], when_measured(collection.mspct[[1]]))
  expect_equal(when_measured(my.spct)[[2]], when_measured(collection.mspct[[2]]))

  expect_equal(names(what_measured(my.spct)), names(collection.mspct))
  expect_equal(what_measured(my.spct)[[1]], what_measured(collection.mspct[[1]]))
  expect_equal(what_measured(my.spct)[[2]], what_measured(collection.mspct[[2]]))

  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[1]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[2]]))

  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[1]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[2]]))

  expect_equal(names(getNormalized(my.spct)), names(collection.mspct))
  expect_equal(getNormalized(my.spct)[[1]], getNormalized(collection.mspct[[1]]))
  expect_equal(getNormalized(my.spct)[[2]], getNormalized(collection.mspct[[2]]))

  expect_equal(names(getNormalization(my.spct)), names(collection.mspct))
  expect_equal(getNormalization(my.spct)[[1]], getNormalization(collection.mspct[[1]]))
  expect_equal(getNormalization(my.spct)[[2]], getNormalization(collection.mspct[[2]]))

})

test_that("subset2mspct works with normalization and no IDs", {
  # test shared code only once
  my.spct <- normalize(sun_evening.spct)

  expect_silent(collection.mspct <- subset2mspct(my.spct, drop.idx = TRUE))

  expect_equal(class(my.spct), class(collection.mspct[[1]]))
  expect_equal(class(my.spct), class(collection.mspct[[2]]))

  expect_named(collection.mspct[[1]], setdiff(colnames(my.spct),
                                              getIdFactor(my.spct)))
  expect_true(setequal(names(collection.mspct),
                       setdiff(levels(my.spct[[getIdFactor(my.spct)]]),
                               getIdFactor(my.spct))))
  expect_equal(length(collection.mspct), getMultipleWl(my.spct))

  expect_equal(getMultipleWl(collection.mspct[[1]]), 1)
  expect_equal(getMultipleWl(collection.mspct[[2]]), 1)

  expect_equal(where_measured(my.spct),
               where_measured(collection.mspct[[1]]))
  expect_equal(where_measured(my.spct),
               where_measured(collection.mspct[[2]]))

  expect_equal(names(when_measured(my.spct)), names(collection.mspct))
  expect_equal(when_measured(my.spct)[[1]], when_measured(collection.mspct[[1]]))
  expect_equal(when_measured(my.spct)[[2]], when_measured(collection.mspct[[2]]))

  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[1]]))
  expect_equal(what_measured(my.spct), what_measured(collection.mspct[[2]]))

  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[1]]))
  expect_equal(how_measured(my.spct), how_measured(collection.mspct[[2]]))

  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[1]]))
  expect_equal(is_normalized(my.spct), is_normalized(collection.mspct[[2]]))

  expect_equal(names(getNormalized(my.spct)), names(collection.mspct))
  expect_equal(getNormalized(my.spct)[[1]], getNormalized(collection.mspct[[1]]))
  expect_equal(getNormalized(my.spct)[[2]], getNormalized(collection.mspct[[2]]))

  expect_equal(names(getNormalization(my.spct)), names(collection.mspct))
  expect_equal(getNormalization(my.spct)[[1]], getNormalization(collection.mspct[[1]]))
  expect_equal(getNormalization(my.spct)[[2]], getNormalization(collection.mspct[[2]]))

})

