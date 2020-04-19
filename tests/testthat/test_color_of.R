context("set_get")

test_that("default", {
  expect_warning(color_of("abc"))
  expect_warning(color_of(NA_character_))
  expect_warning(color_of(polyester.spct))
})

test_that("numeric", {
  expect_equal(color_of(numeric()), character())
  expect_silent(color_of(numeric()))
  expect_is(color_of(NA_real_), "character")
  expect_true(is.na(color_of(NA_real_)))
  expect_silent(color_of(NA_real_))
  expect_is(color_of(NA_integer_), "character")
  expect_true(is.na(color_of(NA_integer_)))
  expect_silent(color_of(NA_integer_))
  expect_error(color_of(-1))
  expect_warning(color_of(300, chroma.type = "zz"))
  expect_warning(color_of(300, chroma.type = 2))
  expect_warning(color_of(300, chroma.type = sun.spct))
  expect_equal(unname(color_of(200)), "#000000")
  expect_equal(names(color_of(200)), "wl.200.nm.CMF")
  expect_equal(unname(color_of(430)), "#1F00FF")
  expect_equal(names(color_of(430)), "wl.430.nm.CMF")
  expect_equal(unname(color_of(430, chroma.type = "CMF")), "#1F00FF")
  expect_equal(names(color_of(430, chroma.type = "CMF")), "wl.430.nm.CMF")
  expect_equal(unname(color_of(430, chroma.type = "CC")),  "#1000E0")
  expect_equal(names(color_of(430, chroma.type = "CC")), "wl.430.nm.CC")
  expect_equal(unname(color_of(430, chroma.type = beesxyzCMF.spct)),  "#00FF00")
  expect_equal(names(color_of(430, chroma.type = beesxyzCMF.spct)), "wl.430.nm.chroma")
  expect_equal(unname(color_of(c(460, 500, 600))),  c("#0000FF", "#00A92E", "#FF2E00"))
  expect_equal(names(color_of(c(460, 500, 600))), c("wl.460.nm.CMF", "wl.500.nm.CMF" ,"wl.600.nm.CMF"))
  expect_equal(unname(color_of(c(460, 500, 600), chroma.type = "CMF")),
               c("#0000FF", "#00A92E", "#FF2E00"))
  expect_equal(names(color_of(c(460, 500, 600), chroma.type = "CMF")),
               c("wl.460.nm.CMF", "wl.500.nm.CMF" ,"wl.600.nm.CMF"))
  expect_equal(unname(color_of(c(460, 500, 600), chroma.type = "CC")),
               c( "#0000DD", "#00FF4E", "#FF1900"))
  expect_equal(names(color_of(c(460, 500, 600), chroma.type = "CC")),
               c("wl.460.nm.CC", "wl.500.nm.CC" ,"wl.600.nm.CC"))
  expect_equal(unname(color_of(c(460, 500, 600), chroma.type = beesxyzCMF.spct)),
               c( "#00FB00", "#FF0005", "#E00004"))
  expect_equal(names(color_of(c(460, 500, 600), chroma.type = beesxyzCMF.spct)),
               c("wl.460.nm.chroma", "wl.500.nm.chroma" ,"wl.600.nm.chroma"))
})

test_that("source_spct", {
  my.spct <- sun.spct[300:400]
#  wl_range(my.spct)
  expect_warning(color_of(my.spct, chroma.type = "zz"))
  expect_warning(color_of(my.spct, chroma.type = 2))
  expect_warning(color_of(my.spct, chroma.type = sun.spct))
  expect_equal(unname(color_of(my.spct)), "#FF1400")
  expect_equal(names(color_of(my.spct)), "source.CMF")
  expect_equal(unname(color_of(my.spct, chroma.type = "CMF")), "#FF1400")
  expect_equal(names(color_of(my.spct, chroma.type = "CMF")), "source.CMF")
  expect_equal(unname(color_of(my.spct, chroma.type = "CC")),  "#FF0000")
  expect_equal(names(color_of(my.spct, chroma.type = "CC")), "source.CC")
  expect_equal(unname(color_of(my.spct, chroma.type = beesxyzCMF.spct)),  "#860002")
  expect_equal(names(color_of(my.spct, chroma.type = beesxyzCMF.spct)), "source.chroma")
})

test_that("source_mspct", {
  my.mspct <- source_mspct(list(sun1 = sun.spct[300:400], sun2 = sun.spct[400:500]))
#  wl_range(my.mspct)
  expect_equal(color_of(my.mspct)[["color"]], c( "#FF1400", "#060000"))
  expect_equal(as.character(color_of(my.mspct)[["spct.idx"]]), names(my.mspct))
  expect_equal(color_of(my.mspct, chroma.type = "CMF")[["color"]], c( "#FF1400", "#060000"))
  expect_equal(as.character(color_of(my.mspct, chroma.type = "CMF")[["spct.idx"]]), names(my.mspct))
  expect_equal(color_of(my.mspct, chroma.type = "CC")[["color"]], c( "#FF0000", "#FF0000"))
  expect_equal(as.character(color_of(my.mspct, chroma.type = "CC")[["spct.idx"]]), names(my.mspct))
})

test_that("list", {
  my.ls <- list(waveband(c(400, 500), wb.name = "wb1"),
                waveband(c(500, 600), wb.name = "wb2"))
  expect_equal(unname(color_of(my.ls)), c("#000EFF", "#5FFF00"))
  expect_equal(names(color_of(my.ls)), c("wb1.CMF", "wb2.CMF"))
})
