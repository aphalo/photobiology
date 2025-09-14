context("normalize methods for multiple spectra")

test_that("normalization of source_mspct works", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)

  normalized.mspct <- normalize(my.mspct)

  expect_lt(abs(max(normalized.mspct[[1]]$s.e.irrad) - 1), 1e-10)
  expect_lt(abs(max(normalized.mspct[[2]]$s.e.irrad) - 1), 1e-10)
  expect_lt(abs(max(normalized.mspct[[3]]$s.e.irrad) - 1), 1e-10)
  expect_lt(abs(max(normalized.mspct[[4]]$s.e.irrad) - 1), 1e-10)
  expect_lt(abs(max(normalized.mspct[[5]]$s.e.irrad) - 1), 1e-10)

  expect_equal(getNormalization(normalized.mspct[[1]])$norm.factor, 1)
  expect_equal(getNormalization(normalized.mspct[[2]])$norm.factor, 1/2)
  expect_equal(getNormalization(normalized.mspct[[3]])$norm.factor, 1/3)
  expect_equal(getNormalization(normalized.mspct[[4]])$norm.factor, 1/4)
  expect_equal(getNormalization(normalized.mspct[[5]])$norm.factor, 1/5)

  expect_equal(getNormalization(normalized.mspct[[1]])$norm.type, "max")
  expect_equal(getNormalization(normalized.mspct[[2]])$norm.type, "max")
  expect_equal(getNormalization(normalized.mspct[[3]])$norm.type, "max")
  expect_equal(getNormalization(normalized.mspct[[4]])$norm.type, "max")
  expect_equal(getNormalization(normalized.mspct[[5]])$norm.type, "max")

  expect_equal(getNormalization(normalized.mspct[[1]])$norm.wl, 400)
  expect_equal(getNormalization(normalized.mspct[[2]])$norm.wl, 400)
  expect_equal(getNormalization(normalized.mspct[[3]])$norm.wl, 400)
  expect_equal(getNormalization(normalized.mspct[[4]])$norm.wl, 400)
  expect_equal(getNormalization(normalized.mspct[[5]])$norm.wl, 400)

  expect_equal(getNormalization(normalized.mspct[[1]])$norm.range, c(400, 410))
  expect_equal(getNormalization(normalized.mspct[[2]])$norm.range, c(400, 410))
  expect_equal(getNormalization(normalized.mspct[[3]])$norm.range, c(400, 410))
  expect_equal(getNormalization(normalized.mspct[[4]])$norm.range, c(400, 410))
  expect_equal(getNormalization(normalized.mspct[[5]])$norm.range, c(400, 410))

  expect_equal(class(my.mspct), class(normalized.mspct))

})

test_that("normalization of long source_spct works", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)
  my.spct <- rbindspct(my.mspct, idfactor = "test.id")

  expect_false(is_normalised(my.spct))
  normalized.spct <- normalize(my.spct)
  expect_true(is_normalised(normalized.spct))

  expect_equal(getIdFactor(my.spct), getIdFactor(normalized.spct))
  expect_equal(colnames(my.spct), colnames(normalized.spct))
  expect_equal(getMultipleWl(my.spct), getMultipleWl(normalized.spct))
  expect_equal(nrow(my.spct), nrow(normalized.spct))
  expect_equal(my.spct$w.length, normalized.spct$w.length)
  expect_equal(class(my.spct), class(normalized.spct))

  expect_equal(names(getNormalized(normalized.spct)),
               levels(my.spct[[getIdFactor(my.spct)]]))
  expect_equal(names(getNormalization(normalized.spct)),
               levels(my.spct[[getIdFactor(my.spct)]]))
  expect_true(all(unlist(is_normalized(normalized.spct), use.names = FALSE)))
  expect_false(anyNA(unlist(getNormalisation(normalized.spct), use.names = FALSE)))

  expect_equal(getNormalization(normalized.spct)[[1]]$norm.factor, 1)
  expect_equal(getNormalization(normalized.spct)[[2]]$norm.factor, 1/2)
  expect_equal(getNormalization(normalized.spct)[[3]]$norm.factor, 1/3)
  expect_equal(getNormalization(normalized.spct)[[4]]$norm.factor, 1/4)
  expect_equal(getNormalization(normalized.spct)[[5]]$norm.factor, 1/5)

  expect_equal(getNormalization(normalized.spct)[[1]]$norm.type, "max")
  expect_equal(getNormalization(normalized.spct)[[2]]$norm.type, "max")
  expect_equal(getNormalization(normalized.spct)[[3]]$norm.type, "max")
  expect_equal(getNormalization(normalized.spct)[[4]]$norm.type, "max")
  expect_equal(getNormalization(normalized.spct)[[5]]$norm.type, "max")

  expect_equal(getNormalization(normalized.spct)[[1]]$norm.wl, 400)
  expect_equal(getNormalization(normalized.spct)[[2]]$norm.wl, 400)
  expect_equal(getNormalization(normalized.spct)[[3]]$norm.wl, 400)
  expect_equal(getNormalization(normalized.spct)[[4]]$norm.wl, 400)
  expect_equal(getNormalization(normalized.spct)[[5]]$norm.wl, 400)

  expect_equal(getNormalization(normalized.spct)[[1]]$norm.range, c(400, 410))
  expect_equal(getNormalization(normalized.spct)[[2]]$norm.range, c(400, 410))
  expect_equal(getNormalization(normalized.spct)[[3]]$norm.range, c(400, 410))
  expect_equal(getNormalization(normalized.spct)[[4]]$norm.range, c(400, 410))
  expect_equal(getNormalization(normalized.spct)[[5]]$norm.range, c(400, 410))


  norm.l <-
    lapply(list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct),
           normalize)
  norm.mspct <- source_mspct(norm.l)
  expect_true(all(unlist(is_normalized(norm.mspct), use.names = FALSE)))
  expect_false(anyNA(unlist(getNormalisation(norm.mspct), use.names = FALSE)))

  norm.spct <- rbindspct(norm.mspct, idfactor = "test.id")
  expect_true(all(unlist(is_normalized(norm.mspct), use.names = FALSE)))
  expect_false(anyNA(unlist(getNormalisation(norm.mspct), use.names = FALSE)))

  expect_equal(names(getNormalized(norm.spct)),
               levels(norm.spct[[getIdFactor(norm.spct)]]))
  expect_equal(names(getNormalization(norm.mspct)),
               levels(norm.spct[[getIdFactor(norm.spct)]]))

})
