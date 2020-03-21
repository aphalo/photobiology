library("photobiology")

context("na.exclude")

test_that("source_spct", {

  my.spct <- my.ref.spct <- e2q(sun.spct, action = "add")
  my.spct[c(4,8), c("s.e.irrad", "s.q.irrad")] <- NA
  my.ref.spct <- my.ref.spct[-c(4 ,8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 3)

  my.spct <- my.ref.spct <- q2e(sun.spct, action = "replace")
  my.spct[c(4,8), "s.e.irrad"] <- NA
  my.ref.spct <- my.ref.spct[-c(4 ,8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 2)

  my.spct <- my.ref.spct <- e2q(sun.spct, action = "replace")
  my.spct[c(4,8), "s.q.irrad"] <- NA
  my.ref.spct <- my.ref.spct[-c(4 ,8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 2)

})

test_that("response_spct", {

  my.spct <- my.ref.spct <- ccd.spct
  my.spct[c(4,8), "s.q.response"] <- NA
  my.ref.spct <- my.ref.spct[-c(4 ,8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 2)

  my.spct <- my.ref.spct <- q2e(ccd.spct, action = "add")
  my.spct[c(4,8), c("s.e.response", "s.q.response")] <- NA
  my.ref.spct <- my.ref.spct[-c(4 ,8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 3)

  my.spct <- my.ref.spct <- q2e(ccd.spct, action = "replace")
  my.spct[c(4,8), "s.e.response"] <- NA
  my.ref.spct <- my.ref.spct[-c(4 ,8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 2)

})

test_that("filter_spct", {

  my.spct <- my.ref.spct <- polyester.spct
  my.spct[c(4,8), "Tfr"] <- NA
  my.ref.spct <- my.ref.spct[-c(4 ,8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 2)

  my.spct <- my.ref.spct <- T2A(polyester.spct, action = "add")
  my.spct[c(4,8), c("Tfr", "A")] <- NA
  my.ref.spct <- my.ref.spct[-c(4 ,8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 3)

  my.spct <- my.ref.spct <- T2A(polyester.spct, action = "replace")
  my.spct[c(4,8), "A"] <- NA
  my.ref.spct <- my.ref.spct[-c(4 ,8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 2)

  })

test_that("response_spct", {

  my.spct <- my.ref.spct <- green_leaf.spct
  my.spct[c(4,8), "Rfr"] <- NA
  my.ref.spct <- my.ref.spct[-c(4 ,8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 2)

})

test_that("cps_spct", {

  my.spct <- my.ref.spct <- white_led.cps_spct
  my.spct[c(4,8), "cps"] <- NA
  my.ref.spct <- my.ref.spct[-c(4, 8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 2)

})

test_that("chroma_spct", {

  my.spct <- my.ref.spct <- beesxyzCMF.spct
  my.spct[c(4, 8), c("x", "y", "z")] <- NA
  my.ref.spct <- my.ref.spct[-c(4 ,8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 4)

  my.spct <- my.ref.spct <- beesxyzCMF.spct
  my.spct[c(4, 8), c("x", "y")] <- NA
  my.ref.spct <- my.ref.spct[-c(4 ,8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 4)

  my.spct <- my.ref.spct <- beesxyzCMF.spct
  my.spct[c(4, 8), "z"] <- NA
  my.ref.spct <- my.ref.spct[-c(4 ,8), ]
#  expect_equal(na.exclude(my.spct), my.ref.spct)
  expect_equal(as.numeric(na.action(na.exclude(my.spct))), c(4, 8))
  expect_is(na.action(na.exclude(my.spct)), "exclude")
  expect_equal(ncol(na.exclude(my.spct)), 4)

})

