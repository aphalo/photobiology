library(photobiology)

context("multi_spct")

test_that("constructors", {

  spct.l <- as.list(photobiology::filter_cps.mspct)
  my.mspct <- cps_mspct(spct.l)
  expect_true(is.cps_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("filter", "no.filter"))
  expect_true(is.cps_spct(my.mspct[["filter"]]))
  expect_true(is.cps_spct(my.mspct[["no.filter"]]))

  my.mspct <- as.generic_mspct(my.mspct)
  expect_false(is.cps_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("filter", "no.filter"))
  expect_true(is.cps_spct(my.mspct[["filter"]]))
  expect_true(is.cps_spct(my.mspct[["no.filter"]]))

  my.mspct <- as.generic_mspct(my.mspct, force.spct.class = TRUE)
  expect_false(is.cps_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("filter", "no.filter"))
  expect_false(is.cps_spct(my.mspct[["filter"]]))
  expect_false(is.cps_spct(my.mspct[["no.filter"]]))
  expect_true(is.any_spct(my.mspct[["filter"]]))
  expect_true(is.any_spct(my.mspct[["no.filter"]]))

  my.mspct <- as.cps_mspct(my.mspct)
  expect_true(is.cps_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("filter", "no.filter"))
  expect_true(is.cps_spct(my.mspct[["filter"]]))
  expect_true(is.cps_spct(my.mspct[["no.filter"]]))

  expect_error(as.filter_spct(my.mspct))
  expect_error(as.reflector_spct(my.mspct))
  expect_error(as.object_spct(my.mspct))
  expect_error(as.response_spct(my.mspct))
  expect_error(as.source_spct(my.mspct))
  expect_error(as.cps_spct(my.mspct))
  expect_error(as.raw_spct(my.mspct))

  empty.mspct <- cps_mspct()
  expect_true(is.cps_mspct(empty.mspct))
  expect_true(is.any_mspct(empty.mspct))
  expect_true(is.null(names(empty.mspct)))

 })

test_that("cps_mspct", {

  my1.spct <- setCpsSpct(data.frame(w.length = 400:410, cps_1 = 0.1, cps_2 = 0.1))
  my2.spct <- setCpsSpct(data.frame(w.length = 400:410, cps_1 = 0.2, cps_2 = 0.2))
  my3.spct <- setCpsSpct(data.frame(w.length = 400:410, cps_1 = 0.3, cps_2 = 0.3))
  my4.spct <- setCpsSpct(data.frame(w.length = 400:410, cps_1 = 0.4, cps_2 = 0.4))
  my5.spct <- setCpsSpct(data.frame(w.length = 400:410, cps_1 = 0.5, cps_2 = 0.5))

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- cps_mspct(spct.l)

  expect_equal(paste("spct", 1:length(spct.l), sep = "_"), names(my.mspct))

  expect_equal(class(my.mspct)[1:2], c("cps_mspct", "generic_mspct") )
  expect_equal(attr(my.mspct, "mspct.version", exact = TRUE), 2)
  expect_equal(ncol(my.mspct), 1)
  expect_equal(nrow(my.mspct), length(spct.l))
  expect_equal(attr(my.mspct, "mspct.byrow", exact = TRUE), FALSE)

  named_spct.l <- list(one = my1.spct,
                       two = my2.spct,
                       three = my3.spct,
                       four = my4.spct,
                       five = my5.spct)

  my_named.mspct <- cps_mspct(named_spct.l)

  expect_equal(names(named_spct.l), names(my_named.mspct))
  expect_equal(class(my_named.mspct)[1:2], c("cps_mspct", "generic_mspct") )
  expect_equal(attr(my_named.mspct, "mspct.version", exact = TRUE), 2)
  expect_equal(ncol(my_named.mspct), 1)
  expect_equal(nrow(my.mspct), length(named_spct.l))
  expect_equal(attr(my_named.mspct, "mspct.byrow", exact = TRUE), FALSE)

  expect_equal(length(my.mspct), 5)
  expect_equal(length(my_named.mspct), 5)

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

  # my_wide.df <- data.frame(w.length = 300:400, A = 1, B = 2, C = 3)
  # my1_df.mspct <- split2filter_mspct(my_wide.df)
  # expect_equal(c("A", "B", "C"), names(my1_df.mspct))
  # expect_equal(nrow(my_wide.df), nrow(my1_df.mspct[[1]]))
  # expect_equal(class(my1_df.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  # expect_equal(class(my1_df.mspct[[1]])[1:2], c("filter_spct", "generic_spct") )
  # expect_equal(class(my1_df.mspct[[2]])[1:2], c("filter_spct", "generic_spct") )
  # expect_equal(class(my1_df.mspct[[3]])[1:2], c("filter_spct", "generic_spct") )

# constructor methods for 'long' data frames --------------------------------

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           cps = c(rep(1, 11), rep(2, 11), rep(3, 11)),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "cps_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("cps_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("cps_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("cps_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("cps_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 2L )
  expect_equal(ncol(my2_df.mspct[[2]]), 2L )
  expect_equal(ncol(my2_df.mspct[[3]]), 2L )

  expect_named(my2_df.mspct[[1]], c("w.length", "cps") )
  expect_named(my2_df.mspct[[2]], c("w.length", "cps") )
  expect_named(my2_df.mspct[[3]], c("w.length", "cps") )

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           cps_1 = c(rep(1, 11), rep(2, 11), rep(3, 11)),
                           cps_2 = c(rep(11, 11), rep(12, 11), rep(13, 11)),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "cps_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("cps_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("cps_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("cps_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("cps_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 3L )
  expect_equal(ncol(my2_df.mspct[[2]]), 3L )
  expect_equal(ncol(my2_df.mspct[[3]]), 3L )

  expect_named(my2_df.mspct[[1]], c("w.length", "cps_1", "cps_2") )
  expect_named(my2_df.mspct[[2]], c("w.length", "cps_1", "cps_2") )
  expect_named(my2_df.mspct[[3]], c("w.length", "cps_1", "cps_2") )

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           cps = c(rep(1, 11), rep(2, 11), rep(3, 11)),
                           other = c(rep("one", 11), rep("two", 11), rep("three", 11)),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "cps_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("cps_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("cps_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("cps_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("cps_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 3L )
  expect_equal(ncol(my2_df.mspct[[2]]), 3L )
  expect_equal(ncol(my2_df.mspct[[3]]), 3L )

  expect_named(my2_df.mspct[[1]], c("w.length", "cps", "other") )
  expect_named(my2_df.mspct[[2]], c("w.length", "cps", "other") )
  expect_named(my2_df.mspct[[3]], c("w.length", "cps", "other") )


  # constructor methods for 'long' spct objects -----------------------------

  expect_warning(my_long.spct <- rbindspct(spct.l))
  my3_df.mspct <- subset2mspct(my_long.spct)

  expect_equal(paste("spct", 1:5, sep = "_"), names(my3_df.mspct))
  expect_equal(class(my3_df.mspct)[1:2], c("cps_mspct", "generic_mspct") )
  expect_equal(class(my3_df.mspct[[1]])[1:2], c("cps_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[2]])[1:2], c("cps_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[3]])[1:2], c("cps_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[4]])[1:2], c("cps_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[5]])[1:2], c("cps_spct", "generic_spct") )

  # print -------------------------------------------------------------------

  expect_equal(print(my.mspct), my.mspct)

  # clean -------------------------------------------------------------------

  expect_equal(clean(my.mspct), my.mspct)
})


test_that("raw_mspct", {

  my1.spct <- setRawSpct(data.frame(w.length = 400:410, counts_1 = 0.1, counts_2 = 0.1))
  my2.spct <- setRawSpct(data.frame(w.length = 400:410, counts_1 = 0.2, counts_2 = 0.2))
  my3.spct <- setRawSpct(data.frame(w.length = 400:410, counts_1 = 0.3, counts_2 = 0.3))
  my4.spct <- setRawSpct(data.frame(w.length = 400:410, counts_1 = 0.4, counts_2 = 0.4))
  my5.spct <- setRawSpct(data.frame(w.length = 400:410, counts_1 = 0.5, counts_2 = 0.5))

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- raw_mspct(spct.l)

  expect_equal(paste("spct", 1:length(spct.l), sep = "_"), names(my.mspct))

  expect_equal(class(my.mspct)[1:2], c("raw_mspct", "generic_mspct") )
  expect_equal(attr(my.mspct, "mspct.version", exact = TRUE), 2)
  expect_equal(ncol(my.mspct), 1)
  expect_equal(nrow(my.mspct), length(spct.l))
  expect_equal(attr(my.mspct, "mspct.byrow", exact = TRUE), FALSE)

  named_spct.l <- list(one = my1.spct,
                       two = my2.spct,
                       three = my3.spct,
                       four = my4.spct,
                       five = my5.spct)

  my_named.mspct <- raw_mspct(named_spct.l)

  expect_equal(names(named_spct.l), names(my_named.mspct))
  expect_equal(class(my_named.mspct)[1:2], c("raw_mspct", "generic_mspct") )
  expect_equal(attr(my_named.mspct, "mspct.version", exact = TRUE), 2)
  expect_equal(ncol(my_named.mspct), 1)
  expect_equal(nrow(my.mspct), length(named_spct.l))
  expect_equal(attr(my_named.mspct, "mspct.byrow", exact = TRUE), FALSE)

  expect_equal(length(my.mspct), 5)
  expect_equal(length(my_named.mspct), 5)

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

  # my_wide.df <- data.frame(w.length = 300:400, A = 1, B = 2, C = 3)
  # my1_df.mspct <- split2filter_mspct(my_wide.df)
  # expect_equal(c("A", "B", "C"), names(my1_df.mspct))
  # expect_equal(nrow(my_wide.df), nrow(my1_df.mspct[[1]]))
  # expect_equal(class(my1_df.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  # expect_equal(class(my1_df.mspct[[1]])[1:2], c("filter_spct", "generic_spct") )
  # expect_equal(class(my1_df.mspct[[2]])[1:2], c("filter_spct", "generic_spct") )
  # expect_equal(class(my1_df.mspct[[3]])[1:2], c("filter_spct", "generic_spct") )

  # constructor methods for 'long' data frames --------------------------------

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           cps = c(rep(1, 11), rep(2, 11), rep(3, 11)),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "cps_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("cps_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("cps_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("cps_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("cps_spct", "generic_spct") )

  # constructor methods for 'long' spct objects -----------------------------

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           counts = c(rep(1, 11), rep(2, 11), rep(3, 11)),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "raw_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("raw_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("raw_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("raw_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("raw_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 2L )
  expect_equal(ncol(my2_df.mspct[[2]]), 2L )
  expect_equal(ncol(my2_df.mspct[[3]]), 2L )

  expect_named(my2_df.mspct[[1]], c("w.length", "counts") )
  expect_named(my2_df.mspct[[2]], c("w.length", "counts") )
  expect_named(my2_df.mspct[[3]], c("w.length", "counts") )

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           counts_1 = c(rep(1, 11), rep(2, 11), rep(3, 11)),
                           counts_2 = c(rep(11, 11), rep(12, 11), rep(13, 11)),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "raw_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("raw_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("raw_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("raw_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("raw_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 3L )
  expect_equal(ncol(my2_df.mspct[[2]]), 3L )
  expect_equal(ncol(my2_df.mspct[[3]]), 3L )

  expect_named(my2_df.mspct[[1]], c("w.length", "counts_1", "counts_2") )
  expect_named(my2_df.mspct[[2]], c("w.length", "counts_1", "counts_2") )
  expect_named(my2_df.mspct[[3]], c("w.length", "counts_1", "counts_2") )

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           counts = c(rep(1, 11), rep(2, 11), rep(3, 11)),
                           other = c(rep("one", 11), rep("two", 11), rep("three", 11)),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "raw_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("raw_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("raw_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("raw_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("raw_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 3L )
  expect_equal(ncol(my2_df.mspct[[2]]), 3L )
  expect_equal(ncol(my2_df.mspct[[3]]), 3L )

  expect_named(my2_df.mspct[[1]], c("w.length", "counts", "other") )
  expect_named(my2_df.mspct[[2]], c("w.length", "counts", "other") )
  expect_named(my2_df.mspct[[3]], c("w.length", "counts", "other") )

  # print -------------------------------------------------------------------

  expect_equal(print(my.mspct), my.mspct)

  # clean -------------------------------------------------------------------

  expect_equal(clean(my.mspct), my.mspct)
})


## need to add similar tests for other classes
## need to add additional methods
