library(photobiology)

test.print <- FALSE

context("multi_spct")

test_that("constructors", {

  my.mspct <- solute_mspct(list(water = water.spct, pha = phenylalanine.spct))
  expect_true(is.solute_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("water", "pha"))
  expect_true(is.solute_spct(my.mspct[["water"]]))
  expect_true(is.solute_spct(my.mspct[["pha"]]))

  my.mspct <- as.generic_mspct(my.mspct)
  expect_false(is.solute_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("water", "pha"))
  expect_true(is.solute_spct(my.mspct[["water"]]))
  expect_true(is.solute_spct(my.mspct[["pha"]]))
  expect_true(is.any_spct(my.mspct[["water"]]))
  expect_true(is.any_spct(my.mspct[["pha"]]))

  my.mspct <- as.generic_mspct(my.mspct, force.spct.class = TRUE)
  expect_false(is.solute_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("water", "pha"))
  expect_false(is.filter_spct(my.mspct[["water"]]))
  expect_false(is.filter_spct(my.mspct[["pha"]]))
  expect_true(is.any_spct(my.mspct[["water"]]))
  expect_true(is.any_spct(my.mspct[["pha"]]))

  my.mspct <- as.solute_mspct(my.mspct)
  expect_true(is.solute_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("water", "pha"))
  expect_true(is.solute_spct(my.mspct[["water"]]))
  expect_true(is.solute_spct(my.mspct[["pha"]]))

  my.mspct <- as.solute_mspct(water.spct)
  expect_true(is.solute_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("spct_1"))
  expect_true(is.solute_spct(my.mspct[["spct_1"]]))

  my.mspct <- as.solute_mspct(data.frame(w.length = 400:500, K.mole = 0.5))
  expect_true(is.solute_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("spct_1"))
  expect_true(is.solute_spct(my.mspct[["spct_1"]]))

  my.mspct <- as.solute_mspct(tibble::tibble(w.length = 400:500, K.mole = 0.5))
  expect_true(is.solute_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("spct_1"))
  expect_true(is.solute_spct(my.mspct[["spct_1"]]))

  my.mspct <- as.solute_mspct(list(data.frame(w.length = 400:500, K.mole = 0.1),
                                   data.frame(w.length = 400:500, K.mole = 0.2)))
  expect_true(is.solute_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("spct_1", "spct_2"))
  expect_true(is.solute_spct(my.mspct[["spct_1"]]))
  expect_true(is.solute_spct(my.mspct[["spct_1"]]))

  my.mspct <- as.solute_mspct(list(A = data.frame(w.length = 400:500, K.mole = 0.1),
                                   B = data.frame(w.length = 400:500, K.mole = 0.2)))
  expect_true(is.solute_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("A", "B"))
  expect_true(is.solute_spct(my.mspct[["A"]]))
  expect_true(is.solute_spct(my.mspct[["B"]]))

  expect_message(as.solute_mspct(1))
  expect_message(as.solute_mspct("abc"))
  expect_message(as.solute_mspct(TRUE))
  expect_error(as.solute_mspct(list(w.length = 400:500,
                                    K.mole = rep(0.3, 101))))

  expect_error(as.solute_spct(my.mspct))
  expect_error(as.filter_spct(my.mspct))
  expect_error(as.reflector_spct(my.mspct))
  expect_error(as.object_spct(my.mspct))
  expect_error(as.response_spct(my.mspct))
  expect_error(as.cps_spct(my.mspct))
  expect_error(as.raw_spct(my.mspct))

  empty.mspct <- solute_mspct()
  expect_true(is.solute_mspct(empty.mspct))
  expect_true(is.any_mspct(empty.mspct))
  expect_true(is.null(names(empty.mspct)))

  })

test_that("solute_mspct", {

  my1.spct <- solute_spct(w.length = 400:410, K.mole = 1)
  my2.spct <- solute_spct(w.length = 400:410, K.mole = 0.2)
  my3.spct <- solute_spct(w.length = 400:410, K.mole = 0.3)
  my4.spct <- solute_spct(w.length = 400:410, K.mole = 0.4)
  my5.spct <- solute_spct(w.length = 400:410, K.mole = 0.5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- solute_mspct(spct.l)

  expect_equal(paste("spct", seq_len(length(spct.l)), sep = "_"), names(my.mspct))

  expect_equal(class(my.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(attr(my.mspct, "mspct.version", exact = TRUE), 3)
  expect_equal(ncol(my.mspct), 1)
  expect_equal(nrow(my.mspct), length(spct.l))
  expect_equal(attr(my.mspct, "mspct.byrow", exact = TRUE), FALSE)

  named_spct.l <- list(one = my1.spct,
                       two = my2.spct,
                       three = my3.spct,
                       four = my4.spct,
                       five = my5.spct)

  my_named.mspct <- solute_mspct(named_spct.l)

  expect_equal(names(named_spct.l), names(my_named.mspct))
  expect_equal(class(my_named.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(attr(my_named.mspct, "mspct.version", exact = TRUE), 3)
  expect_equal(ncol(my_named.mspct), 1)
  expect_equal(nrow(my.mspct), length(named_spct.l))
  expect_equal(attr(my_named.mspct, "mspct.byrow", exact = TRUE), FALSE)

  expect_equal(length(my.mspct), 5)
  expect_equal(length(my_named.mspct), 5)

  # min ---------------------------------------------------------------------

  expect_equal(wl_min(my.mspct)[["spct.idx"]],
               factor(paste("spct", 1:5, sep = "_")))
  expect_equal(levels(max(my_named.mspct)[["spct.idx"]]),
               c("one", "two", "three", "four", "five"))
  expect_equal(wl_min(my.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(wl_min(my_named.mspct)[["min.wl"]], rep(400, 5))

  # max ---------------------------------------------------------------------

  expect_equal(wl_max(my.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(wl_max(my_named.mspct)[["max.wl"]], rep(410, 5))

  # range -------------------------------------------------------------------

  expect_equal(wl_range(my.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(wl_range(my_named.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(wl_range(my.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(wl_range(my_named.mspct)[["max.wl"]], rep(410, 5))

  # constructor methods for 'wide' data frames ------------------------------

  my_wide.df <- data.frame(w.length = 300:400, A = 1, B = 0.5, C = 0.1)
  my1_df.mspct <- split2solute_mspct(my_wide.df)
  expect_equal(c("A", "B", "C"), names(my1_df.mspct))
  expect_equal(nrow(my_wide.df), nrow(my1_df.mspct[[1]]))
  expect_equal(class(my1_df.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(class(my1_df.mspct[[1]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my1_df.mspct[[2]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my1_df.mspct[[3]])[1:2], c("solute_spct", "generic_spct") )

# constructor methods for 'long' data frames --------------------------------

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           K.mole = c(rep(1, 11), rep(0.5, 11), rep(0.1, 11)),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "solute_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("solute_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 2L )
  expect_equal(ncol(my2_df.mspct[[2]]), 2L )
  expect_equal(ncol(my2_df.mspct[[3]]), 2L )

  expect_named(my2_df.mspct[[1]], c("w.length", "K.mole") )
  expect_named(my2_df.mspct[[2]], c("w.length", "K.mole") )
  expect_named(my2_df.mspct[[3]], c("w.length", "K.mole") )

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           K.mole = c(rep(1, 11), rep(0.5, 11), rep(0.1, 11)),
                           other = c(rep("one", 11), rep("two", 11), rep("three", 11) ),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "solute_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("solute_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 3L )
  expect_equal(ncol(my2_df.mspct[[2]]), 3L )
  expect_equal(ncol(my2_df.mspct[[3]]), 3L )

  expect_named(my2_df.mspct[[1]], c("w.length", "K.mole", "other") )
  expect_named(my2_df.mspct[[2]], c("w.length", "K.mole", "other") )
  expect_named(my2_df.mspct[[3]], c("w.length", "K.mole", "other") )

# constructor methods for 'long' spct objects -----------------------------

  my_long.spct <- rbindspct(spct.l)
  my3_df.mspct <- subset2mspct(my_long.spct)

  expect_equal(paste("spct", 1:5, sep = "_"), names(my3_df.mspct))
  expect_equal(class(my3_df.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(class(my3_df.mspct[[1]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[2]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[3]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[4]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[5]])[1:2], c("solute_spct", "generic_spct") )

  expect_equal(ncol(my3_df.mspct[[1]]), 2L )
  expect_equal(ncol(my3_df.mspct[[2]]), 2L )
  expect_equal(ncol(my3_df.mspct[[3]]), 2L )

  expect_named(my3_df.mspct[[1]], c("w.length", "K.mole") )
  expect_named(my3_df.mspct[[2]], c("w.length", "K.mole") )
  expect_named(my3_df.mspct[[3]], c("w.length", "K.mole") )

  # print -------------------------------------------------------------------

  if (test.print) expect_equal(print(my.mspct), my.mspct)

  # clean -------------------------------------------------------------------

  expect_equal(clean(my.mspct), my.mspct)
})


test_that("solute_mspct_absorption", {

  my1.spct <- solute_spct(w.length = 400:410, K.mole = 1, K.type = "absorption")
  my2.spct <- solute_spct(w.length = 400:410, K.mole = 0.2, K.type = "absorption")
  my3.spct <- solute_spct(w.length = 400:410, K.mole = 0.3, K.type = "absorption")
  my4.spct <- solute_spct(w.length = 400:410, K.mole = 0.4, K.type = "absorption")
  my5.spct <- solute_spct(w.length = 400:410, K.mole = 0.5, K.type = "absorption")

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- solute_mspct(spct.l)

  expect_equal(paste("spct", seq_len(length(spct.l)), sep = "_"), names(my.mspct))

  expect_equal(class(my.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(attr(my.mspct, "mspct.version", exact = TRUE), 3)
  expect_equal(ncol(my.mspct), 1)
  expect_equal(nrow(my.mspct), length(spct.l))
  expect_equal(attr(my.mspct, "mspct.byrow", exact = TRUE), FALSE)

  named_spct.l <- list(one = my1.spct,
                       two = my2.spct,
                       three = my3.spct,
                       four = my4.spct,
                       five = my5.spct)

  my_named.mspct <- solute_mspct(named_spct.l)

  expect_equal(names(named_spct.l), names(my_named.mspct))
  expect_equal(class(my_named.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(attr(my_named.mspct, "mspct.version", exact = TRUE), 3)
  expect_equal(ncol(my_named.mspct), 1)
  expect_equal(nrow(my.mspct), length(named_spct.l))
  expect_equal(attr(my_named.mspct, "mspct.byrow", exact = TRUE), FALSE)

  expect_equal(length(my.mspct), 5)
  expect_equal(length(my_named.mspct), 5)



  # min ---------------------------------------------------------------------

  expect_equal(wl_min(my.mspct)[["spct.idx"]],
               factor(paste("spct", 1:5, sep = "_")))
  expect_equal(levels(max(my_named.mspct)[["spct.idx"]]),
               c("one", "two", "three", "four", "five"))
  expect_equal(wl_min(my.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(wl_min(my_named.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(wl_min(my.mspct), wl_min(my.mspct))
  expect_equal(wl_min(my_named.mspct), wl_min(my_named.mspct))

  # max ---------------------------------------------------------------------

  expect_equal(wl_max(my.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(wl_max(my_named.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(wl_max(my.mspct), wl_max(my.mspct))
  expect_equal(wl_max(my_named.mspct), wl_max(my_named.mspct))

  # range -------------------------------------------------------------------

  expect_equal(wl_range(my.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(wl_range(my_named.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(wl_range(my.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(wl_range(my_named.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(wl_range(my.mspct), wl_range(my.mspct))
  expect_equal(wl_range(my_named.mspct), wl_range(my_named.mspct))

  # midpoint ----------------------------------------------------------------

  expect_equal(wl_midpoint(my.mspct)[["midpoint.wl"]], rep(405, 5))
  expect_equal(wl_midpoint(my_named.mspct)[["midpoint.wl"]], rep(405, 5))
  expect_equal(wl_midpoint(my.mspct), wl_midpoint(my.mspct))
  expect_equal(wl_midpoint(my_named.mspct), wl_midpoint(my_named.mspct))

  # expanse -----------------------------------------------------------------

  expect_equal(wl_expanse(my.mspct)[["expanse_V1"]], rep(10, 5))
  expect_equal(wl_expanse(my_named.mspct)[["expanse_V1"]], rep(10, 5))
  expect_equal(wl_expanse(my.mspct), wl_expanse(my.mspct))
  expect_equal(wl_expanse(my_named.mspct), wl_expanse(my_named.mspct))
  expect_message(spread(my.mspct[[1]]))

  # constructor methods for 'wide' data frames ------------------------------

  my_wide.df <- data.frame(w.length = 300:400, A = 1, B = 0.5, C = 0.1)
  my1_df.mspct <- split2solute_mspct(my_wide.df)
  expect_equal(c("A", "B", "C"), names(my1_df.mspct))
  expect_equal(nrow(my_wide.df), nrow(my1_df.mspct[[1]]))
  expect_equal(class(my1_df.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(class(my1_df.mspct[[1]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my1_df.mspct[[2]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my1_df.mspct[[3]])[1:2], c("solute_spct", "generic_spct") )

  # constructor methods for 'long' data frames --------------------------------

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           K.mole = c(rep(1, 11), rep(0.5, 11), rep(0.1, 11)),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "solute_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("solute_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 2L )
  expect_equal(ncol(my2_df.mspct[[2]]), 2L )
  expect_equal(ncol(my2_df.mspct[[3]]), 2L )

  expect_named(my2_df.mspct[[1]], c("w.length", "K.mole") )
  expect_named(my2_df.mspct[[2]], c("w.length", "K.mole") )
  expect_named(my2_df.mspct[[3]], c("w.length", "K.mole") )

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           K.mole = c(rep(1, 11), rep(0.5, 11), rep(0.1, 11)),
                           other = c(rep("one", 11), rep("two", 11), rep("three", 11) ),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "solute_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("solute_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 3L )
  expect_equal(ncol(my2_df.mspct[[2]]), 3L )
  expect_equal(ncol(my2_df.mspct[[3]]), 3L )

  expect_named(my2_df.mspct[[1]], c("w.length", "K.mole", "other") )
  expect_named(my2_df.mspct[[2]], c("w.length", "K.mole", "other") )
  expect_named(my2_df.mspct[[3]], c("w.length", "K.mole", "other") )

  # constructor methods for 'long' spct objects -----------------------------

  my_long.spct <- rbindspct(spct.l)
  my3_df.mspct <- subset2mspct(my_long.spct)

  expect_equal(paste("spct", 1:5, sep = "_"), names(my3_df.mspct))
  expect_equal(class(my3_df.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(class(my3_df.mspct[[1]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[2]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[3]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[4]])[1:2], c("solute_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[5]])[1:2], c("solute_spct", "generic_spct") )

  expect_equal(ncol(my3_df.mspct[[1]]), 2L )
  expect_equal(ncol(my3_df.mspct[[2]]), 2L )
  expect_equal(ncol(my3_df.mspct[[3]]), 2L )

  expect_named(my3_df.mspct[[1]], c("w.length", "K.mole") )
  expect_named(my3_df.mspct[[2]], c("w.length", "K.mole") )
  expect_named(my3_df.mspct[[3]], c("w.length", "K.mole") )

  # print -------------------------------------------------------------------

  if (test.print) expect_equal(print(my.mspct), my.mspct)

  # clean -------------------------------------------------------------------

  expect_equal(clean(my.mspct), my.mspct)
})

test_that("solute_mspct_attr", {

  my1.spct <- solute_spct(w.length = 400:410, K.mole = 0.5, K.type = "absorption")
  setWhatMeasured(my1.spct, "first spectrum")
  setWhenMeasured(my1.spct, lubridate::ymd("2018-03-03", tz = "UTC"))
  setWhereMeasured(my1.spct, data.frame(lat = -30, lon = +80))
  my2.spct <- solute_spct(w.length = 400:410, K.mole = 0.5, K.type = "absorption")
  setWhatMeasured(my2.spct, "second spectrum")
  setWhenMeasured(my2.spct, lubridate::ymd_hm("2018-03-03 12:30", tz = "UTC"))
  setWhereMeasured(my2.spct, data.frame(lat = 5, lon = 20))

  spct.l <- list(my1.spct, my2.spct)
  my.mspct <- solute_mspct(spct.l)

  expect_equal(paste("spct", seq_len(length(spct.l)), sep = "_"), names(my.mspct))

  expect_equal(class(my.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(attr(my.mspct, "mspct.version", exact = TRUE), 3)
  expect_equal(ncol(my.mspct), 1)
  expect_equal(nrow(my.mspct), length(spct.l))
  expect_equal(attr(my.mspct, "mspct.byrow", exact = TRUE), FALSE)

  named_spct.l <- list(one = my1.spct,
                       two = my2.spct)

  my_named.mspct <- solute_mspct(named_spct.l)

  expect_equal(names(named_spct.l), names(my_named.mspct))
  expect_equal(class(my_named.mspct)[1:2], c("solute_mspct", "generic_mspct") )
  expect_equal(attr(my_named.mspct, "mspct.version", exact = TRUE), 3)
  expect_equal(ncol(my_named.mspct), 1)
  expect_equal(nrow(my.mspct), length(named_spct.l))
  expect_equal(attr(my_named.mspct, "mspct.byrow", exact = TRUE), FALSE)

  expect_equal(length(my.mspct), 2)
  expect_equal(length(my_named.mspct), 2)

})


