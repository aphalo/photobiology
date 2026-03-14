context("attributes under math opers")

## the order in which attributes are stored is altered
## the order of the columns in the dataframe can also change
## this is irrelevant as long as names are used to access them
test_that("merge attr in operations with s.e.irrad", {
  energy_as_default()
  reference.spct <- q2e(white_led.source_spct, action = "replace")

  result.spct <- (reference.spct + reference.spct) / 2

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$s.e.irrad, reference.spct$s.e.irrad)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               setdiff(names(attributes(reference.spct)),
                       c("linearized", "merged.cps.cols")),
               ignore.order = TRUE)

  result.spct <- reference.spct + 0

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$s.e.irrad, reference.spct$s.e.irrad)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               names(attributes(reference.spct)),
               ignore.order = TRUE)

  result1.spct <- reference.spct - reference.spct
  result2.spct <- reference.spct * 0

  expect_named(result1.spct, names(result2.spct), ignore.order = TRUE)
  expect_equal(result1.spct$w.length, result2.spct$w.length)
  expect_equal(result1.spct$s.e.irrad, result2.spct$s.e.irrad)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result1.spct),
               setdiff(names(attributes(result2.spct)),
                       c("linearized", "merged.cps.cols")),
               ignore.order = TRUE)
})

## the order in which attributes are stored is altered
## the order of the columns in the dataframe can also change
## this is irrelevant as long as names are used to access them
test_that("merge attr in operations with s.q.irrad", {
  photon_as_default()
  reference.spct <- e2q(white_led.source_spct, action = "replace")

  result.spct <- (reference.spct + reference.spct) / 2

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$s.q.irrad, reference.spct$s.q.irrad)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               setdiff(names(attributes(reference.spct)),
                       c("linearized", "merged.cps.cols")),
               ignore.order = TRUE)

  result.spct <- reference.spct + 0

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$s.q.irrad, reference.spct$s.q.irrad)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               names(attributes(reference.spct)),
               ignore.order = TRUE)

  result1.spct <- reference.spct - reference.spct
  result2.spct <- reference.spct * 0

  expect_named(result1.spct, names(result2.spct), ignore.order = TRUE)
  expect_equal(result1.spct$w.length, result2.spct$w.length)
  expect_equal(result1.spct$s.q.irrad, result2.spct$s.q.irrad)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result1.spct),
               setdiff(names(attributes(result2.spct)),
                       c("linearized", "merged.cps.cols")),
               ignore.order = TRUE)
})

## the order in which attributes are stored is altered
## the order of the columns in the dataframe can also change
## this is irrelevant as long as names are used to access them
test_that("merge attr in operations with s.e.response", {
  energy_as_default()
  reference.spct <- q2e(ccd.spct, action = "replace")

  result.spct <- (reference.spct + reference.spct) / 2

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$s.e.response, reference.spct$s.e.response)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               setdiff(names(attributes(reference.spct)),
                       c("linearized", "merged.cps.cols")),
               ignore.order = TRUE)

  result.spct <- reference.spct + 0

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$s.e.response, reference.spct$s.e.response)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               names(attributes(reference.spct)),
               ignore.order = TRUE)

  result1.spct <- reference.spct - reference.spct
  result2.spct <- reference.spct * 0

  expect_named(result1.spct, names(result2.spct), ignore.order = TRUE)
  expect_equal(result1.spct$w.length, result2.spct$w.length)
  expect_equal(result1.spct$s.e.response, result2.spct$s.e.response)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result1.spct),
               setdiff(names(attributes(result2.spct)),
                       c("linearized", "merged.cps.cols")),
               ignore.order = TRUE)
})

## the order in which attributes are stored is altered
## the order of the columns in the dataframe can also change
## this is irrelevant as long as names are used to access them
test_that("merge attr in operations with s.q.response", {
  photon_as_default()
  reference.spct <- e2q(ccd.spct, action = "replace")

  result.spct <- (reference.spct + reference.spct) / 2

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$s.q.response, reference.spct$s.q.response)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               setdiff(names(attributes(reference.spct)),
                       c("linearized", "merged.cps.cols")),
               ignore.order = TRUE)

  result.spct <- reference.spct + 0

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$s.q.response, reference.spct$s.q.response)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               names(attributes(reference.spct)),
               ignore.order = TRUE)

  result1.spct <- reference.spct - reference.spct
  result2.spct <- reference.spct * 0

  expect_named(result1.spct, names(result2.spct), ignore.order = TRUE)
  expect_equal(result1.spct$w.length, result2.spct$w.length)
  expect_equal(result1.spct$s.q.response, result2.spct$s.q.response)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result1.spct),
               setdiff(names(attributes(result2.spct)),
                       c("linearized", "merged.cps.cols")),
               ignore.order = TRUE)
})

## the order in which attributes are stored is altered
## the order of the columns in the dataframe can also change
## this is irrelevant as long as names are used to access them
test_that("merge attr in operations with Tfr", {
  Tfr_as_default()
  reference.spct <- polyester.spct

  expect_error(reference.spct + reference.spct)
  expect_error(reference.spct - reference.spct)

  result.spct <- sqrt((reference.spct * reference.spct))

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$Tfr, reference.spct$Tfr)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               setdiff(names(attributes(reference.spct)),
                       c("linearized", "merged.cps.cols")),
               ignore.order = TRUE)

  result.spct <- reference.spct + 0

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$Tfr, reference.spct$Tfr)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               names(attributes(reference.spct)),
               ignore.order = TRUE)

  result.spct <- reference.spct * 1
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$Tfr, reference.spct$Tfr)
  expect_named(attributes(result.spct),
               names(attributes(reference.spct)),
               ignore.order = TRUE)

  expect_equal(polyester.spct / Inf, polyester.spct * 0)
})

## the order in which attributes are stored is altered
## the order of the columns in the dataframe can also change
## this is irrelevant as long as names are used to access them
test_that("merge attr in operations with Afr", {
  Afr_as_default()
  reference.spct <- any2Afr(polyester.spct, action = "replace")

  expect_error(reference.spct + reference.spct)
  expect_error(reference.spct - reference.spct)

  result.spct <- sqrt((reference.spct * reference.spct))

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$Afr, reference.spct$Afr)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               setdiff(names(attributes(reference.spct)),
                       c("linearized", "merged.cps.cols")),
               ignore.order = TRUE)

  result.spct <- reference.spct + 0

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$Afr, reference.spct$Afr)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               names(attributes(reference.spct)),
               ignore.order = TRUE)

  result.spct <- reference.spct * 1
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$Afr, reference.spct$Afr)
  expect_named(attributes(result.spct),
               names(attributes(reference.spct)),
               ignore.order = TRUE)

  expect_equal(polyester.spct / Inf, polyester.spct * 0)
})

## the order in which attributes are stored is altered
## the order of the columns in the dataframe can also change
## this is irrelevant as long as names are used to access them
test_that("merge attr in operations with A", {
  A_as_default()
  reference.spct <- any2A(polyester.spct, action = "replace")

  expect_error(reference.spct * reference.spct)
  expect_error(reference.spct / reference.spct)

  result.spct <- (reference.spct + reference.spct) / 2

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$A, reference.spct$A)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               setdiff(names(attributes(reference.spct)),
                       c("linearized", "merged.cps.cols")),
               ignore.order = TRUE)

  result.spct <- reference.spct + 0

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$A, reference.spct$A)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               names(attributes(reference.spct)),
               ignore.order = TRUE)

  result.spct <- reference.spct * 1
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$A, reference.spct$A)
  expect_named(attributes(result.spct),
               names(attributes(reference.spct)),
               ignore.order = TRUE)

  expect_equal(reference.spct / Inf, reference.spct * 0)
})

## the order in which attributes are stored is altered
## the order of the columns in the dataframe can also change
## this is irrelevant as long as names are used to access them
test_that("merge attr in operations with Rfr", {
  reference.spct <- green_leaf.spct

  #  expect_error(reference.spct * reference.spct)
  #  expect_error(reference.spct / reference.spct)

  result.spct <- sqrt((reference.spct * reference.spct))

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$Rfr, reference.spct$Rfr)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               setdiff(names(attributes(reference.spct)),
                       c("linearized", "merged.cps.cols")),
               ignore.order = TRUE)

  result.spct <- reference.spct + 0

  expect_named(result.spct, names(reference.spct), ignore.order = TRUE)
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$Rfr, reference.spct$Rfr)
  # these missing attributes are irrelevant after maths operations
  # additional ones could need to be dropped!
  expect_named(attributes(result.spct),
               names(attributes(reference.spct)),
               ignore.order = TRUE)

  result.spct <- reference.spct * 1
  expect_equal(result.spct$w.length, reference.spct$w.length)
  expect_equal(result.spct$Rfr, reference.spct$Rfr)
  expect_named(attributes(result.spct),
               names(attributes(reference.spct)),
               ignore.order = TRUE)

  expect_equal(reference.spct / Inf, reference.spct * 0)
})

# tests for object_spct

# tests for solute_spct
