library("photobiology")
context("solute_spct")

test_that("constructor fraction", {

  empty.spct <- solute_spct()
  expect_true(is.solute_spct(empty.spct))
  expect_true(is.any_spct(empty.spct))
  expect_named(empty.spct, c("w.length", "K.mole"))
  expect_equal(nrow(empty.spct), 0L)

  my.spct <- solute_spct(w.length = 400:409, K.mole = 0.1)
  expect_equal(class(my.spct)[1:2], c("solute_spct", "generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 2)

  expect_warning(solute_spct(w.length = 400:409, K.mole = -0.1))
  expect_equal(my.spct[["K.mole"]], rep(0.1, length.out = 10))
  expect_equal(my.spct[["w.length"]], 400:409)
  expect_named(my.spct, c("w.length", "K.mole"))
  expect_null(attr(my.spct, "time.unit", exact = TRUE))

  expect_true(is.solute_spct(my.spct))
  expect_true(is.any_spct(my.spct))
  expect_false(is.raw_spct(my.spct))
  expect_false(is.source_spct(my.spct))
  expect_false(is.cps_spct(my.spct))
  expect_false(is.filter_spct(my.spct))
  expect_false(is.object_spct(my.spct))
  expect_false(is.response_spct(my.spct))
  expect_false(is.chroma_spct(my.spct))

  my.df <- data.frame(w.length = 400:409, K.mole = 0.1)
  my.spct <- as.solute_spct(my.df)

  expect_equal(class(my.spct)[1:2], c("solute_spct", "generic_spct") )
  expect_equal(attr(my.spct, "spct.version", exact = TRUE), 2)
  expect_named(my.spct, c("w.length", "K.mole"))
  expect_true(is.solute_spct(my.spct))
  expect_true(is.any_spct(my.spct))

})

test_that("oper", {

  my.e.spct <- solute_spct(w.length = 400:409, K.mole = 0.1)
  my.2e.spct <- solute_spct(w.length = 400:409, K.mole = 0.2)

  expect_equal(my.e.spct + my.e.spct,  my.2e.spct)
  expect_equal(my.e.spct * 2, my.2e.spct)
  expect_equal(my.e.spct * 2L, my.2e.spct)
  expect_equal(my.2e.spct / 2, my.e.spct)
  expect_equal(my.2e.spct / 2L, my.e.spct)
  expect_warning(-my.e.spct)
  expect_equal(my.e.spct + 0.1, my.2e.spct)
  expect_equal(suppressWarnings(-my.2e.spct / -2), my.e.spct)
  expect_equal(suppressWarnings(-my.2e.spct / -2L), my.e.spct)
  expect_equal(2 * my.e.spct, my.2e.spct)
  expect_equal(suppressWarnings( 1 / (2 / my.2e.spct)), my.e.spct)
  expect_equal(suppressWarnings( 1 / my.e.spct),
               suppressWarnings(my.e.spct^-1))
  expect_equal(my.2e.spct %/% 2L, my.e.spct %/% 1L)
  expect_equal(my.2e.spct %% 2L / 2, my.e.spct %% 1L)

})


test_that("math", {

  my.e.spct <- solute_spct(w.length = 400:409, K.mole = 0.1)
  my.2e.spct <- solute_spct(w.length = 400:409, K.mole = 0.2)

  expect_warning(log10(my.e.spct))
  expect_equal(suppressWarnings(log10(my.e.spct)[["K.mole"]]),
               rep(log10(0.1), length.out = 10))
  expect_warning(log(my.e.spct))
  expect_equal(suppressWarnings(log(my.e.spct)[["K.mole"]]),
               rep(log(0.1), length.out = 10))
  expect_warning(log(my.e.spct, 2))
  expect_equal(suppressWarnings(log(my.e.spct, 2)[["K.mole"]]),
               rep(log(0.1, 2), length.out = 10))
  expect_equal(suppressWarnings(exp(my.e.spct)[["K.mole"]]),
               rep(exp(0.1), length.out = 10))
  expect_equal(sqrt(my.e.spct)[["K.mole"]],  rep(sqrt(0.1), length.out = 10))

})
