library(photobiology)

context("multi_spct")

test_that("constructors", {

  my.mspct <- source_mspct(list(sun1 = sun.spct, sun2 = sun.spct))
  expect_true(is.source_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("sun1", "sun2"))
  expect_true(is.source_spct(my.mspct[["sun1"]]))
  expect_true(is.source_spct(my.mspct[["sun2"]]))

  my.mspct <- as.generic_mspct(my.mspct)
  expect_false(is.source_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("sun1", "sun2"))
  expect_true(is.source_spct(my.mspct[["sun1"]]))
  expect_true(is.source_spct(my.mspct[["sun2"]]))
  expect_true(is.any_spct(my.mspct[["sun1"]]))
  expect_true(is.any_spct(my.mspct[["sun2"]]))

  my.mspct <- as.generic_mspct(my.mspct, force.spct.class = TRUE)
  expect_false(is.source_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("sun1", "sun2"))
  expect_false(is.source_spct(my.mspct[["sun1"]]))
  expect_false(is.source_spct(my.mspct[["sun2"]]))
  expect_true(is.any_spct(my.mspct[["sun1"]]))
  expect_true(is.any_spct(my.mspct[["sun2"]]))

  my.mspct <- as.source_mspct(my.mspct)
  expect_true(is.source_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("sun1", "sun2"))
  expect_true(is.source_spct(my.mspct[["sun1"]]))
  expect_true(is.source_spct(my.mspct[["sun2"]]))

  expect_error(as.filter_spct(my.mspct))
  expect_error(as.reflector_spct(my.mspct))
  expect_error(as.object_spct(my.mspct))
  expect_error(as.response_spct(my.mspct))
  expect_error(as.cps_spct(my.mspct))
  expect_error(as.raw_spct(my.mspct))

  empty.mspct <- source_mspct()
  expect_true(is.source_mspct(empty.mspct))
  expect_true(is.any_mspct(empty.mspct))
  expect_true(is.null(names(empty.mspct)))

  empty.mspct <- response_mspct()
  expect_true(is.response_mspct(empty.mspct))
  expect_true(is.any_mspct(empty.mspct))
  expect_true(is.null(names(empty.mspct)))

  empty.mspct <- filter_mspct()
  expect_true(is.filter_mspct(empty.mspct))
  expect_true(is.any_mspct(empty.mspct))
  expect_true(is.null(names(empty.mspct)))

  empty.mspct <- reflector_mspct()
  expect_true(is.reflector_mspct(empty.mspct))
  expect_true(is.any_mspct(empty.mspct))
  expect_true(is.null(names(empty.mspct)))

  empty.mspct <- object_mspct()
  expect_true(is.object_mspct(empty.mspct))
  expect_true(is.any_mspct(empty.mspct))
  expect_true(is.null(names(empty.mspct)))

  empty.mspct <- cps_mspct()
  expect_true(is.cps_mspct(empty.mspct))
  expect_true(is.any_mspct(empty.mspct))
  expect_true(is.null(names(empty.mspct)))

  empty.mspct <- raw_mspct()
  expect_true(is.raw_mspct(empty.mspct))
  expect_true(is.any_mspct(empty.mspct))
  expect_true(is.null(names(empty.mspct)))

  empty.mspct <- chroma_mspct()
  expect_true(is.chroma_mspct(empty.mspct))
  expect_true(is.any_mspct(empty.mspct))
  expect_true(is.null(names(empty.mspct)))

  })

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)

  expect_equal(paste("spct", 1:length(spct.l), sep = "_"), names(my.mspct))

  expect_equal(class(my.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(attr(my.mspct, "mspct.version", exact = TRUE), 2)
  expect_equal(ncol(my.mspct), 1)
  expect_equal(nrow(my.mspct), length(spct.l))
  expect_equal(attr(my.mspct, "mspct.byrow", exact = TRUE), FALSE)

  named_spct.l <- list(one = my1.spct,
                       two = my2.spct,
                       three = my3.spct,
                       four = my4.spct,
                       five = my5.spct)

  my_named.mspct <- source_mspct(named_spct.l)

  expect_equal(names(named_spct.l), names(my_named.mspct))
  expect_equal(class(my_named.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(attr(my_named.mspct, "mspct.version", exact = TRUE), 2)
  expect_equal(ncol(my_named.mspct), 1)
  expect_equal(nrow(my.mspct), length(named_spct.l))
  expect_equal(attr(my_named.mspct, "mspct.byrow", exact = TRUE), FALSE)

  expect_equal(length(my.mspct), 5)
  expect_equal(length(my_named.mspct), 5)

  expect_true(is.data.frame(irrad(my.mspct)))
  expect_true(is.data.frame(irrad(my_named.mspct)))

  expect_equal(irrad(my.mspct)[["spct.idx"]],
               factor(paste("spct", 1:5, sep = "_")))
  expect_equal(levels(irrad(my_named.mspct)[["spct.idx"]]),
               c("one", "two", "three", "four", "five"))

  expect_equal(round(irrad(my.mspct)[["irrad_Total"]], 3), 1:5 * 10)
  expect_equal(round(irrad(my_named.mspct)[["irrad_Total"]], 3), 1:5 * 10)

# irrad -------------------------------------------------------------------

  expect_equal(round(irrad(my.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
                           )[["irrad_A"]], 3), 1:5 * 5)
  expect_equal(round(irrad(my.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
                           )[["irrad_B"]], 3), 1:5 * 5)

  expect_equal(round(irrad(my_named.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
  )[["irrad_A"]], 3), 1:5 * 5)
  expect_equal(round(irrad(my_named.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
  )[["irrad_B"]], 3), 1:5 * 5)

# min ---------------------------------------------------------------------

  expect_equal(min(my.mspct)[["spct.idx"]],
               factor(paste("spct", 1:5, sep = "_")))
  expect_equal(levels(max(my_named.mspct)[["spct.idx"]]),
               c("one", "two", "three", "four", "five"))
  expect_equal(min(my.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(min(my_named.mspct)[["min.wl"]], rep(400, 5))

  # max ---------------------------------------------------------------------

  expect_equal(max(my.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(max(my_named.mspct)[["max.wl"]], rep(410, 5))

  # range -------------------------------------------------------------------

  expect_equal(range(my.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(range(my_named.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(range(my.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(range(my_named.mspct)[["max.wl"]], rep(410, 5))

  # constructor methods for 'wide' data frames ------------------------------

  my_wide.df <- data.frame(w.length = 300:400, A = 1, B = 2, C = 3)
  my1_df.mspct <- split2source_mspct(my_wide.df)
  expect_equal(c("A", "B", "C"), names(my1_df.mspct))
  expect_equal(nrow(my_wide.df), nrow(my1_df.mspct[[1]]))
  expect_equal(class(my1_df.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(class(my1_df.mspct[[1]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my1_df.mspct[[2]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my1_df.mspct[[3]])[1:2], c("source_spct", "generic_spct") )

# constructor methods for 'long' data frames --------------------------------

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           s.e.irrad = c(rep(1, 11), rep(2, 11), rep(3, 11)),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "source_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("source_spct", "generic_spct") )

  # constructor methods for 'long' spct objects -----------------------------

  my_long.spct <- rbindspct(spct.l)
  my3_df.mspct <- subset2mspct(my_long.spct)

  expect_equal(paste("spct", 1:5, sep = "_"), names(my3_df.mspct))
  expect_equal(class(my3_df.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(class(my3_df.mspct[[1]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[2]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[3]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[4]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[5]])[1:2], c("source_spct", "generic_spct") )

  # print -------------------------------------------------------------------

  expect_equal(print(my.mspct), my.mspct)

  # clean -------------------------------------------------------------------

  expect_equal(clean(my.mspct), my.mspct)
})



## need to add similar tests for other classes
## need to add additional methods
