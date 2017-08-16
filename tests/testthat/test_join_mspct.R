library("photobiology")

context("join_mspct")

test_that("source_mspct", {
  my.mspct <- source_mspct(list(sun1 = sun.spct, sun2 = sun.spct * 2))
  expect_is(my.mspct, "source_mspct")

  expect_silent(my.df <- join_mspct(my.mspct))
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "sun1", "sun2"))
  expect_equal(my.df[["w.length"]], my.mspct[["sun1"]][["w.length"]])
  expect_equal(my.df[["w.length"]], my.mspct[["sun2"]][["w.length"]])
  expect_equal(my.df[["sun1"]], my.mspct[["sun1"]][["s.e.irrad"]])
  expect_equal(my.df[["sun1"]], sun.spct[["s.e.irrad"]])
  expect_equal(my.df[["sun2"]], my.mspct[["sun2"]][["s.e.irrad"]])

  expect_silent(my.df <- join_mspct(my.mspct, unit.out = "energy"))
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "sun1", "sun2"))
  expect_equal(my.df[["w.length"]], my.mspct[["sun1"]][["w.length"]])
  expect_equal(my.df[["w.length"]], my.mspct[["sun2"]][["w.length"]])
  expect_equal(my.df[["sun1"]], my.mspct[["sun1"]][["s.e.irrad"]])
  expect_equal(my.df[["sun1"]], sun.spct[["s.e.irrad"]])
  expect_equal(my.df[["sun2"]], my.mspct[["sun2"]][["s.e.irrad"]])

  expect_silent(my.df <- join_mspct(my.mspct, unit.out = "photon"))
  expect_silent(myq.df <- join_mspct(my.mspct, unit.out = "quantum"))
  expect_equal(my.df, myq.df)
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "sun1", "sun2"))
  expect_equal(my.df[["w.length"]], my.mspct[["sun1"]][["w.length"]])
  expect_equal(my.df[["w.length"]], my.mspct[["sun2"]][["w.length"]])
  expect_equal(my.df[["sun1"]], my.mspct[["sun1"]][["s.q.irrad"]])
  expect_equal(my.df[["sun1"]], sun.spct[["s.q.irrad"]])
  expect_equal(my.df[["sun2"]], e2q(my.mspct[["sun2"]])[["s.q.irrad"]])

  # boundary cases

  my1.mspct <- source_mspct(list(sun1 = sun.spct))
  expect_is(my1.mspct, "source_mspct")

  expect_silent(my.df <- join_mspct(my1.mspct))
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "sun1"))
  expect_equal(ncol(my.df), 2L)
  expect_equal(my.df[["sun1"]], my1.mspct[["sun1"]][["s.e.irrad"]])
  expect_equal(my.df[["sun1"]], sun.spct[["s.e.irrad"]])

  my0.mspct <- source_mspct()
  expect_is(my0.mspct, "source_mspct")

  expect_silent(my.df <- join_mspct(my0.mspct))
  expect_is(my.df, "data.frame")
  expect_equal(nrow(my.df), 0L)
  expect_equal(ncol(my.df), 0L)

})

test_that("response_mspct", {
  my.mspct <- response_mspct(list(ccd1 = ccd.spct, ccd2 = ccd.spct * 2))
  expect_is(my.mspct, "response_mspct")

  expect_silent(my.df <- join_mspct(my.mspct))
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "ccd1", "ccd2"))
  expect_equal(my.df[["w.length"]], my.mspct[["ccd1"]][["w.length"]])
  expect_equal(my.df[["w.length"]], my.mspct[["ccd2"]][["w.length"]])
  expect_equal(my.df[["ccd1"]], q2e(my.mspct[["ccd1"]])[["s.e.response"]])
  expect_equal(my.df[["ccd1"]], q2e(ccd.spct)[["s.e.response"]])
  expect_equal(my.df[["ccd2"]], q2e(my.mspct[["ccd2"]])[["s.e.response"]])

  expect_silent(my.df <- join_mspct(my.mspct, unit.out = "energy"))
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "ccd1", "ccd2"))
  expect_equal(my.df[["w.length"]], my.mspct[["ccd1"]][["w.length"]])
  expect_equal(my.df[["w.length"]], my.mspct[["ccd2"]][["w.length"]])
  expect_equal(my.df[["ccd1"]], q2e(my.mspct[["ccd1"]])[["s.e.response"]])
  expect_equal(my.df[["ccd1"]], q2e(ccd.spct)[["s.e.response"]])
  expect_equal(my.df[["ccd2"]], q2e(my.mspct)[["ccd2"]][["s.e.response"]])

  expect_silent(my.df <- join_mspct(my.mspct, unit.out = "photon"))
  expect_silent(myx.df <- join_mspct(my.mspct, unit.out = "quantum"))
  expect_equal(my.df, myx.df)
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "ccd1", "ccd2"))
  expect_equal(my.df[["w.length"]], my.mspct[["ccd1"]][["w.length"]])
  expect_equal(my.df[["w.length"]], my.mspct[["ccd2"]][["w.length"]])
  expect_equal(my.df[["ccd1"]], my.mspct[["ccd1"]][["s.q.response"]])
  expect_equal(my.df[["ccd1"]], ccd.spct[["s.q.response"]])
  expect_equal(my.df[["ccd2"]], e2q(my.mspct[["ccd2"]])[["s.q.response"]])

  # boundary cases

  my1.mspct <- response_mspct(list(ccd1 = ccd.spct))
  expect_is(my1.mspct, "response_mspct")

  expect_silent(my.df <- join_mspct(my1.mspct))
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "ccd1"))
  expect_equal(ncol(my.df), 2L)
  expect_equal(my.df[["ccd1"]], q2e(my1.mspct[["ccd1"]])[["s.e.response"]])
  expect_equal(my.df[["ccd1"]], q2e(ccd.spct)[["s.e.response"]])

  my0.mspct <- response_mspct()
  expect_is(my0.mspct, "response_mspct")

  expect_silent(my.df <- join_mspct(my0.mspct))
  expect_is(my.df, "data.frame")
  expect_equal(nrow(my.df), 0L)
  expect_equal(ncol(my.df), 0L)

})

test_that("filter_mspct", {
  my.mspct <- filter_mspct(list(pet1 = polyester.spct, pet2 = polyester.spct / 2))
  expect_is(my.mspct, "filter_mspct")

  expect_silent(my.df <- join_mspct(my.mspct))
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "pet1", "pet2"))
  expect_equal(my.df[["w.length"]], my.mspct[["pet1"]][["w.length"]])
  expect_equal(my.df[["w.length"]], my.mspct[["pet2"]][["w.length"]])
  expect_equal(my.df[["pet1"]], my.mspct[["pet1"]][["Tfr"]])
  expect_equal(my.df[["pet1"]], polyester.spct[["Tfr"]])
  expect_equal(my.df[["pet2"]], my.mspct[["pet2"]][["Tfr"]])

  expect_silent(my.df <- join_mspct(my.mspct, qty.out = "transmittance"))
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "pet1", "pet2"))
  expect_equal(my.df[["w.length"]], my.mspct[["pet1"]][["w.length"]])
  expect_equal(my.df[["w.length"]], my.mspct[["pet2"]][["w.length"]])
  expect_equal(my.df[["pet1"]], my.mspct[["pet1"]][["Tfr"]])
  expect_equal(my.df[["pet1"]], polyester.spct[["Tfr"]])
  expect_equal(my.df[["pet2"]], my.mspct[["pet2"]][["Tfr"]])

  expect_silent(my.df <- join_mspct(my.mspct, qty.out = "absorbance"))
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "pet1", "pet2"))
  expect_equal(my.df[["w.length"]], my.mspct[["pet1"]][["w.length"]])
  expect_equal(my.df[["w.length"]], my.mspct[["pet2"]][["w.length"]])
  expect_equal(my.df[["pet1"]], T2A(my.mspct[["pet1"]])[["A"]])
  expect_equal(my.df[["pet1"]], T2A(polyester.spct)[["A"]])
  expect_equal(my.df[["pet2"]], T2A(my.mspct[["pet2"]])[["A"]])

  # boundary cases

  my1.mspct <- filter_mspct(list(pet1 = polyester.spct))
  expect_is(my1.mspct, "filter_mspct")

  expect_silent(my.df <- join_mspct(my1.mspct))
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "pet1"))
  expect_equal(ncol(my.df), 2L)
  expect_equal(my.df[["pet1"]], my1.mspct[["pet1"]][["Tfr"]])
  expect_equal(my.df[["pet1"]], polyester.spct[["Tfr"]])

  my0.mspct <- filter_mspct()
  expect_is(my0.mspct, "filter_mspct")

  expect_silent(my.df <- join_mspct(my0.mspct))
  expect_is(my.df, "data.frame")
  expect_equal(nrow(my.df), 0L)
  expect_equal(ncol(my.df), 0L)

})

test_that("reflector_mspct", {
  my.mspct <- reflector_mspct(list(ler1 = Ler_leaf_rflt.spct, ler2 = Ler_leaf_rflt.spct / 2))
  expect_is(my.mspct, "reflector_mspct")

  expect_silent(my.df <- join_mspct(my.mspct))
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "ler1", "ler2"))
  expect_equal(my.df[["w.length"]], my.mspct[["ler1"]][["w.length"]])
  expect_equal(my.df[["w.length"]], my.mspct[["ler2"]][["w.length"]])
  expect_equal(my.df[["ler1"]], my.mspct[["ler1"]][["Rfr"]])
  expect_equal(my.df[["ler1"]], Ler_leaf_rflt.spct[["Rfr"]])
  expect_equal(my.df[["ler2"]], my.mspct[["ler2"]][["Rfr"]])

  # boundary cases

  my1.mspct <- reflector_mspct(list(ler1 = Ler_leaf_rflt.spct))
  expect_is(my1.mspct, "reflector_mspct")

  expect_silent(my.df <- join_mspct(my1.mspct))
  expect_is(my.df, "data.frame")
  expect_named(my.df, c("w.length", "ler1"))
  expect_equal(ncol(my.df), 2L)
  expect_equal(my.df[["ler1"]], my1.mspct[["ler1"]][["Rfr"]])
  expect_equal(my.df[["ler1"]], Ler_leaf_rflt.spct[["Rfr"]])

  my0.mspct <- reflector_mspct()
  expect_is(my0.mspct, "reflector_mspct")

  expect_silent(my.df <- join_mspct(my0.mspct))
  expect_is(my.df, "data.frame")
  expect_equal(nrow(my.df), 0L)
  expect_equal(ncol(my.df), 0L)

})

test_that("generic_mspct", {
  my.mspct <- generic_mspct(list(ler = Ler_leaf_rflt.spct, pet = polyester.spct))
  expect_is(my.mspct, "generic_mspct")

  expect_error(my.df <- join_mspct(my.mspct))

})

test_that("object_mspct", {
  my.mspct <- object_mspct(list(ler1 = Ler_leaf.spct, ler2 = Ler_leaf.spct))
  expect_is(my.mspct, "object_mspct")

  expect_error(my.df <- join_mspct(my.mspct))

})

test_that("chroma_mspct", {
  my.mspct <- chroma_mspct(list(cie10 = ciexyzCMF10.spct, cie2 = ciexyzCMF2.spct))
  expect_is(my.mspct, "chroma_mspct")

  expect_error(my.df <- join_mspct(my.mspct))

})
