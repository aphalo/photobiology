library(photobiology)

test.print <- FALSE

context("multi_spct")

test_that("constructors", {

  my.mspct <- filter_mspct(list(yellow = yellow_gel.spct, pet = polyester.spct))
  expect_true(is.filter_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("yellow", "pet"))
  expect_true(is.filter_spct(my.mspct[["yellow"]]))
  expect_true(is.filter_spct(my.mspct[["pet"]]))

  my.mspct <- as.generic_mspct(my.mspct)
  expect_false(is.filter_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("yellow", "pet"))
  expect_true(is.filter_spct(my.mspct[["yellow"]]))
  expect_true(is.filter_spct(my.mspct[["pet"]]))
  expect_true(is.any_spct(my.mspct[["yellow"]]))
  expect_true(is.any_spct(my.mspct[["pet"]]))

  my.mspct <- as.generic_mspct(my.mspct, force.spct.class = TRUE)
  expect_false(is.filter_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("yellow", "pet"))
  expect_false(is.filter_spct(my.mspct[["yellow"]]))
  expect_false(is.filter_spct(my.mspct[["pet"]]))
  expect_true(is.any_spct(my.mspct[["yellow"]]))
  expect_true(is.any_spct(my.mspct[["pet"]]))

  my.mspct <- as.filter_mspct(my.mspct)
  expect_true(is.filter_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("yellow", "pet"))
  expect_true(is.filter_spct(my.mspct[["yellow"]]))
  expect_true(is.filter_spct(my.mspct[["pet"]]))

  my.mspct <- as.filter_mspct(polyester.spct)
  expect_true(is.filter_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("spct_1"))
  expect_true(is.filter_spct(my.mspct[["spct_1"]]))

  my.mspct <- as.filter_mspct(data.frame(w.length = 400:500, Tfr = 0.5))
  expect_true(is.filter_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("spct_1"))
  expect_true(is.filter_spct(my.mspct[["spct_1"]]))

  my.mspct <- as.filter_mspct(tibble::tibble(w.length = 400:500, Tfr = 0.5))
  expect_true(is.filter_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("spct_1"))
  expect_true(is.filter_spct(my.mspct[["spct_1"]]))

  my.mspct <- as.filter_mspct(list(data.frame(w.length = 400:500, Tfr = 0.1),
                                   data.frame(w.length = 400:500, Tfr = 0.2)))
  expect_true(is.filter_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("spct_1", "spct_2"))
  expect_true(is.filter_spct(my.mspct[["spct_1"]]))
  expect_true(is.filter_spct(my.mspct[["spct_1"]]))

  my.mspct <- as.filter_mspct(list(A = data.frame(w.length = 400:500, Tfr = 0.1),
                                   B = data.frame(w.length = 400:500, Tfr = 0.2)))
  expect_true(is.filter_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("A", "B"))
  expect_true(is.filter_spct(my.mspct[["A"]]))
  expect_true(is.filter_spct(my.mspct[["B"]]))

  expect_message(as.filter_mspct(1))
  expect_message(as.filter_mspct("abc"))
  expect_message(as.filter_mspct(TRUE))
  expect_error(as.filter_mspct(list(w.length = 400:500,
                                    Tfr = rep(0.3, 101))))

  expect_error(as.filter_spct(my.mspct))
  expect_error(as.reflector_spct(my.mspct))
  expect_error(as.object_spct(my.mspct))
  expect_error(as.response_spct(my.mspct))
  expect_error(as.cps_spct(my.mspct))
  expect_error(as.raw_spct(my.mspct))

  empty.mspct <- filter_mspct()
  expect_true(is.filter_mspct(empty.mspct))
  expect_true(is.any_mspct(empty.mspct))
  expect_true(is.null(names(empty.mspct)))

  })

test_that("filter_mspct", {

  my1.spct <- filter_spct(w.length = 400:410, Tfr = 1)
  my2.spct <- filter_spct(w.length = 400:410, Tfr = 0.2)
  my3.spct <- filter_spct(w.length = 400:410, Tfr = 0.3)
  my4.spct <- filter_spct(w.length = 400:410, Tfr = 0.4)
  my5.spct <- filter_spct(w.length = 400:410, Tfr = 0.5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- filter_mspct(spct.l)

  expect_equal(paste("spct", seq_len(length(spct.l)), sep = "_"), names(my.mspct))

  expect_equal(class(my.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(attr(my.mspct, "mspct.version", exact = TRUE), 3)
  expect_equal(ncol(my.mspct), 1)
  expect_equal(nrow(my.mspct), length(spct.l))
  expect_equal(attr(my.mspct, "mspct.byrow", exact = TRUE), FALSE)

  named_spct.l <- list(one = my1.spct,
                       two = my2.spct,
                       three = my3.spct,
                       four = my4.spct,
                       five = my5.spct)

  my_named.mspct <- filter_mspct(named_spct.l)

  expect_equal(names(named_spct.l), names(my_named.mspct))
  expect_equal(class(my_named.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(attr(my_named.mspct, "mspct.version", exact = TRUE), 3)
  expect_equal(ncol(my_named.mspct), 1)
  expect_equal(nrow(my.mspct), length(named_spct.l))
  expect_equal(attr(my_named.mspct, "mspct.byrow", exact = TRUE), FALSE)

  expect_equal(length(my.mspct), 5)
  expect_equal(length(my_named.mspct), 5)

  expect_true(is.data.frame(transmittance(my.mspct)))
  expect_true(is.data.frame(transmittance(my_named.mspct)))

  expect_equal(transmittance(my.mspct)[["spct.idx"]],
               factor(paste("spct", 1:5, sep = "_")))
  expect_equal(levels(transmittance(my_named.mspct)[["spct.idx"]]),
               c("one", "two", "three", "four", "five"))

#  expect_equal(round(transmittance(my.mspct)[["transmittance_Total"]], 3), 1:5 * 10)
#  expect_equal(round(transmittance(my_named.mspct)[["transmittance_Total"]], 3), 1:5 * 10)

# transmittance -------------------------------------------------------------------

  expect_equal(round(transmittance(my.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
                           )[["Tfr(wl)_A"]], 4), c(1, 2:5 * 0.1))
  expect_equal(round(transmittance(my.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
                           )[["Tfr(wl)_B"]], 4), c(1, 2:5 * 0.1))

  expect_equal(round(transmittance(my_named.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
  )[["Tfr(wl)_A"]], 4), c(1, 2:5 * 0.1))
  expect_equal(round(transmittance(my_named.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
  )[["Tfr(wl)_B"]], 4), c(1, 2:5 * 0.1))

  # relative

  expect_equal(round(transmittance(my.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "relative"
  )[["Tfr/Tfrsum_A"]], 4), rep(0.5, 5))
  expect_equal(round(transmittance(my.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "relative"
  )[["Tfr/Tfrsum_A"]], 4), rep(0.5, 5))

  expect_equal(round(transmittance(my_named.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "relative"
  )[["Tfr/Tfrsum_A"]], 4), rep(0.5, 5))
  expect_equal(round(transmittance(my_named.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "relative"
  )[["Tfr/Tfrsum_A"]], 4), rep(0.5, 5))

  # contribution

  expect_equal(round(transmittance(my.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "contribution"
  )[["Tfr/Tfrtot_A"]], 4), rep(0.5, 5))
  expect_equal(round(transmittance(my.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "contribution"
  )[["Tfr/Tfrtot_A"]], 4), rep(0.5, 5))

  expect_equal(round(transmittance(my_named.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "contribution"
  )[["Tfr/Tfrtot_A"]], 4), rep(0.5, 5))
  expect_equal(round(transmittance(my_named.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "contribution"
  )[["Tfr/Tfrtot_A"]], 4), rep(0.5, 5))

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

  my_wide.df <- data.frame(w.length = 300:400, A = 1, B = 0.5, C = 0.1)
  my1_df.mspct <- split2filter_mspct(my_wide.df)
  expect_equal(c("A", "B", "C"), names(my1_df.mspct))
  expect_equal(nrow(my_wide.df), nrow(my1_df.mspct[[1]]))
  expect_equal(class(my1_df.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(class(my1_df.mspct[[1]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my1_df.mspct[[2]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my1_df.mspct[[3]])[1:2], c("filter_spct", "generic_spct") )

# constructor methods for 'long' data frames --------------------------------

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           Tfr = c(rep(1, 11), rep(0.5, 11), rep(0.1, 11)),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "filter_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("filter_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 2L )
  expect_equal(ncol(my2_df.mspct[[2]]), 2L )
  expect_equal(ncol(my2_df.mspct[[3]]), 2L )

  expect_named(my2_df.mspct[[1]], c("w.length", "Tfr") )
  expect_named(my2_df.mspct[[2]], c("w.length", "Tfr") )
  expect_named(my2_df.mspct[[3]], c("w.length", "Tfr") )

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           Tfr = c(rep(1, 11), rep(0.5, 11), rep(0.1, 11)),
                           other = c(rep("one", 11), rep("two", 11), rep("three", 11) ),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "filter_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("filter_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 3L )
  expect_equal(ncol(my2_df.mspct[[2]]), 3L )
  expect_equal(ncol(my2_df.mspct[[3]]), 3L )

  expect_named(my2_df.mspct[[1]], c("w.length", "Tfr", "other") )
  expect_named(my2_df.mspct[[2]], c("w.length", "Tfr", "other") )
  expect_named(my2_df.mspct[[3]], c("w.length", "Tfr", "other") )

# constructor methods for 'long' spct objects -----------------------------

  my_long.spct <- rbindspct(spct.l)
  my3_df.mspct <- subset2mspct(my_long.spct)

  expect_equal(paste("spct", 1:5, sep = "_"), names(my3_df.mspct))
  expect_equal(class(my3_df.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(class(my3_df.mspct[[1]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[2]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[3]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[4]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[5]])[1:2], c("filter_spct", "generic_spct") )

  expect_equal(ncol(my3_df.mspct[[1]]), 2L )
  expect_equal(ncol(my3_df.mspct[[2]]), 2L )
  expect_equal(ncol(my3_df.mspct[[3]]), 2L )

  expect_named(my3_df.mspct[[1]], c("w.length", "Tfr") )
  expect_named(my3_df.mspct[[2]], c("w.length", "Tfr") )
  expect_named(my3_df.mspct[[3]], c("w.length", "Tfr") )

  # print -------------------------------------------------------------------

  if (test.print) expect_equal(print(my.mspct), my.mspct)

  # clean -------------------------------------------------------------------

  expect_equal(clean(my.mspct), my.mspct)
})


test_that("filter_mspct_internal", {

  my1.spct <- filter_spct(w.length = 400:410, Tfr = 1, Tfr.type = "internal")
  my2.spct <- filter_spct(w.length = 400:410, Tfr = 0.2, Tfr.type = "internal")
  my3.spct <- filter_spct(w.length = 400:410, Tfr = 0.3, Tfr.type = "internal")
  my4.spct <- filter_spct(w.length = 400:410, Tfr = 0.4, Tfr.type = "internal")
  my5.spct <- filter_spct(w.length = 400:410, Tfr = 0.5, Tfr.type = "internal")

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- filter_mspct(spct.l)

  expect_equal(paste("spct", seq_len(length(spct.l)), sep = "_"), names(my.mspct))

  expect_equal(class(my.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(attr(my.mspct, "mspct.version", exact = TRUE), 3)
  expect_equal(ncol(my.mspct), 1)
  expect_equal(nrow(my.mspct), length(spct.l))
  expect_equal(attr(my.mspct, "mspct.byrow", exact = TRUE), FALSE)

  named_spct.l <- list(one = my1.spct,
                       two = my2.spct,
                       three = my3.spct,
                       four = my4.spct,
                       five = my5.spct)

  my_named.mspct <- filter_mspct(named_spct.l)

  expect_equal(names(named_spct.l), names(my_named.mspct))
  expect_equal(class(my_named.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(attr(my_named.mspct, "mspct.version", exact = TRUE), 3)
  expect_equal(ncol(my_named.mspct), 1)
  expect_equal(nrow(my.mspct), length(named_spct.l))
  expect_equal(attr(my_named.mspct, "mspct.byrow", exact = TRUE), FALSE)

  expect_equal(length(my.mspct), 5)
  expect_equal(length(my_named.mspct), 5)

  expect_true(is.data.frame(transmittance(my.mspct)))
  expect_true(is.data.frame(transmittance(my_named.mspct)))

  expect_equal(transmittance(my.mspct)[["spct.idx"]],
               factor(paste("spct", 1:5, sep = "_")))
  expect_equal(levels(transmittance(my_named.mspct)[["spct.idx"]]),
               c("one", "two", "three", "four", "five"))

  #  expect_equal(round(transmittance(my.mspct)[["transmittance_Total"]], 3), 1:5 * 10)
  #  expect_equal(round(transmittance(my_named.mspct)[["transmittance_Total"]], 3), 1:5 * 10)

  # transmittance -------------------------------------------------------------------

  expect_equal(round(transmittance(my.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B"))
  )[["Tfr(wl)_A"]], 4), c(1, 2:5 * 0.1))
  expect_equal(round(transmittance(my.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B"))
  )[["Tfr(wl)_B"]], 4), c(1, 2:5 * 0.1))

  expect_equal(round(transmittance(my_named.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B"))
  )[["Tfr(wl)_A"]], 4), c(1, 2:5 * 0.1))
  expect_equal(round(transmittance(my_named.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B"))
  )[["Tfr(wl)_B"]], 4), c(1, 2:5 * 0.1))

  # relative

  expect_equal(round(transmittance(my.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "relative"
  )[["Tfr/Tfrsum_A"]], 4), rep(0.5, 5))
  expect_equal(round(transmittance(my.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "relative"
  )[["Tfr/Tfrsum_A"]], 4), rep(0.5, 5))

  expect_equal(round(transmittance(my_named.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "relative"
  )[["Tfr/Tfrsum_A"]], 4), rep(0.5, 5))
  expect_equal(round(transmittance(my_named.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "relative"
  )[["Tfr/Tfrsum_A"]], 4), rep(0.5, 5))

  # absorptance -------------------------------------------------------------------

  expect_equal(round(absorptance(my.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B"))
  )[["Afr(wl)_A"]], 4), 1 - c(1, 2:5 * 0.1))
  expect_equal(round(absorptance(my.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B"))
  )[["Afr(wl)_B"]], 4), 1 - c(1, 2:5 * 0.1))

  expect_equal(round(absorptance(my_named.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B"))
  )[["Afr(wl)_A"]], 4), 1 - c(1, 2:5 * 0.1))
  expect_equal(round(absorptance(my_named.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B"))
  )[["Afr(wl)_B"]], 4), 1 - c(1, 2:5 * 0.1))

  # total

  # relative

  expect_equal(round(absorptance(my.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "relative"
  )[["Afr/Afrsum_A"]], 4), c(NaN, rep(0.5, 4)))
  expect_equal(round(absorptance(my.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "relative"
  )[["Afr/Afrsum_B"]], 4), c(NaN, rep(0.5, 4)))

  expect_equal(round(absorptance(my_named.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "relative"
  )[["Afr/Afrsum_A"]], 4), c(NaN, rep(0.5, 4)))
  expect_equal(round(absorptance(my_named.mspct,
                                   list(waveband(c(400,405), wb.name = "A"),
                                        waveband(c(405,410), wb.name = "B")),
                                   quantity = "relative"
  )[["Afr/Afrsum_B"]], 4), c(NaN, rep(0.5, 4)))

  # contribution

  expect_equal(round(absorptance(my.mspct,
                                 list(waveband(c(400,405), wb.name = "A"),
                                      waveband(c(405,410), wb.name = "B")),
                                 quantity = "contribution"
  )[["Afr/Afrtot_A"]], 4), c(NaN, rep(0.5, 4)))
  expect_equal(round(absorptance(my.mspct,
                                 list(waveband(c(400,405), wb.name = "A"),
                                      waveband(c(405,410), wb.name = "B")),
                                 quantity = "contribution"
  )[["Afr/Afrtot_B"]], 4), c(NaN, rep(0.5, 4)))

  expect_equal(round(absorptance(my_named.mspct,
                                 list(waveband(c(400,405), wb.name = "A"),
                                      waveband(c(405,410), wb.name = "B")),
                                 quantity = "contribution"
  )[["Afr/Afrtot_A"]], 4), c(NaN, rep(0.5, 4)))
  expect_equal(round(absorptance(my_named.mspct,
                                 list(waveband(c(400,405), wb.name = "A"),
                                      waveband(c(405,410), wb.name = "B")),
                                 quantity = "contribution"
  )[["Afr/Afrtot_B"]], 4), c(NaN, rep(0.5, 4)))

  # absorbance -------------------------------------------------------------------

  expect_equal(round(absorbance(my.mspct,
                                 list(waveband(c(400,405), wb.name = "A"),
                                      waveband(c(405,410), wb.name = "B"))
  )[["A(wl)_A"]], 3), c(0, 0.699, 0.523, 0.398, 0.301))
  expect_equal(round(absorbance(my.mspct,
                                 list(waveband(c(400,405), wb.name = "A"),
                                      waveband(c(405,410), wb.name = "B"))
  )[["A(wl)_B"]], 3),  c(0, 0.699, 0.523, 0.398, 0.301))

  expect_equal(round(absorbance(my_named.mspct,
                                 list(waveband(c(400,405), wb.name = "A"),
                                      waveband(c(405,410), wb.name = "B"))
  )[["A(wl)_A"]], 3),  c(0, 0.699, 0.523, 0.398, 0.301))
  expect_equal(round(absorbance(my_named.mspct,
                                 list(waveband(c(400,405), wb.name = "A"),
                                      waveband(c(405,410), wb.name = "B"))
  )[["A(wl)_B"]], 3),  c(0, 0.699, 0.523, 0.398, 0.301))

  # relative

  expect_equal(round(absorbance(my.mspct,
                                 list(waveband(c(400,405), wb.name = "A"),
                                      waveband(c(405,410), wb.name = "B")),
                                 quantity = "relative"
  )[["A/Asum_A"]], 4), c(NaN, rep(0.5, 4)))
  expect_equal(round(absorbance(my.mspct,
                                 list(waveband(c(400,405), wb.name = "A"),
                                      waveband(c(405,410), wb.name = "B")),
                                 quantity = "relative"
  )[["A/Asum_A"]], 4), c(NaN, rep(0.5, 4)))

  expect_equal(round(absorbance(my_named.mspct,
                                 list(waveband(c(400,405), wb.name = "A"),
                                      waveband(c(405,410), wb.name = "B")),
                                 quantity = "relative"
  )[["A/Asum_A"]], 4), c(NaN, rep(0.5, 4)))
  expect_equal(round(absorbance(my_named.mspct,
                                 list(waveband(c(400,405), wb.name = "A"),
                                      waveband(c(405,410), wb.name = "B")),
                                 quantity = "relative"
  )[["A/Asum_A"]], 4), c(NaN, rep(0.5, 4)))

  # contribution

  expect_equal(round(absorbance(my.mspct,
                                list(waveband(c(400,405), wb.name = "A"),
                                     waveband(c(405,410), wb.name = "B")),
                                quantity = "contribution"
  )[["A/Atot_A"]], 4), c(NaN, rep(0.5, 4)))
  expect_equal(round(absorbance(my.mspct,
                                list(waveband(c(400,405), wb.name = "A"),
                                     waveband(c(405,410), wb.name = "B")),
                                quantity = "contribution"
  )[["A/Atot_A"]], 4), c(NaN, rep(0.5, 4)))

  expect_equal(round(absorbance(my_named.mspct,
                                list(waveband(c(400,405), wb.name = "A"),
                                     waveband(c(405,410), wb.name = "B")),
                                quantity = "contribution"
  )[["A/Atot_A"]], 4), c(NaN, rep(0.5, 4)))
  expect_equal(round(absorbance(my_named.mspct,
                                list(waveband(c(400,405), wb.name = "A"),
                                     waveband(c(405,410), wb.name = "B")),
                                quantity = "contribution"
  )[["A/Atot_A"]], 4), c(NaN, rep(0.5, 4)))

  # min ---------------------------------------------------------------------

  expect_equal(min(my.mspct)[["spct.idx"]],
               factor(paste("spct", 1:5, sep = "_")))
  expect_equal(levels(max(my_named.mspct)[["spct.idx"]]),
               c("one", "two", "three", "four", "five"))
  expect_equal(min(my.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(min(my_named.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(min(my.mspct), wl_min(my.mspct))
  expect_equal(min(my_named.mspct), wl_min(my_named.mspct))

  # max ---------------------------------------------------------------------

  expect_equal(max(my.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(max(my_named.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(max(my.mspct), wl_max(my.mspct))
  expect_equal(max(my_named.mspct), wl_max(my_named.mspct))

  # range -------------------------------------------------------------------

  expect_equal(range(my.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(range(my_named.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(range(my.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(range(my_named.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(range(my.mspct), wl_range(my.mspct))
  expect_equal(range(my_named.mspct), wl_range(my_named.mspct))

  # midpoint ----------------------------------------------------------------

  expect_equal(midpoint(my.mspct)[["midpoint.wl"]], rep(405, 5))
  expect_equal(midpoint(my_named.mspct)[["midpoint.wl"]], rep(405, 5))
  expect_equal(midpoint(my.mspct), wl_midpoint(my.mspct))
  expect_equal(midpoint(my_named.mspct), wl_midpoint(my_named.mspct))

  # expanse -----------------------------------------------------------------

  expect_equal(expanse(my.mspct)[["expanse_V1"]], rep(10, 5))
  expect_equal(expanse(my_named.mspct)[["expanse_V1"]], rep(10, 5))
  expect_equal(expanse(my.mspct), wl_expanse(my.mspct))
  expect_equal(expanse(my_named.mspct), wl_expanse(my_named.mspct))
  expect_message(spread(my.mspct[[1]]))

  # constructor methods for 'wide' data frames ------------------------------

  my_wide.df <- data.frame(w.length = 300:400, A = 1, B = 0.5, C = 0.1)
  my1_df.mspct <- split2filter_mspct(my_wide.df)
  expect_equal(c("A", "B", "C"), names(my1_df.mspct))
  expect_equal(nrow(my_wide.df), nrow(my1_df.mspct[[1]]))
  expect_equal(class(my1_df.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(class(my1_df.mspct[[1]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my1_df.mspct[[2]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my1_df.mspct[[3]])[1:2], c("filter_spct", "generic_spct") )

  # constructor methods for 'long' data frames --------------------------------

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           Tfr = c(rep(1, 11), rep(0.5, 11), rep(0.1, 11)),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "filter_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("filter_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 2L )
  expect_equal(ncol(my2_df.mspct[[2]]), 2L )
  expect_equal(ncol(my2_df.mspct[[3]]), 2L )

  expect_named(my2_df.mspct[[1]], c("w.length", "Tfr") )
  expect_named(my2_df.mspct[[2]], c("w.length", "Tfr") )
  expect_named(my2_df.mspct[[3]], c("w.length", "Tfr") )

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           Tfr = c(rep(1, 11), rep(0.5, 11), rep(0.1, 11)),
                           other = c(rep("one", 11), rep("two", 11), rep("three", 11) ),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "filter_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("filter_spct", "generic_spct") )

  expect_equal(ncol(my2_df.mspct[[1]]), 3L )
  expect_equal(ncol(my2_df.mspct[[2]]), 3L )
  expect_equal(ncol(my2_df.mspct[[3]]), 3L )

  expect_named(my2_df.mspct[[1]], c("w.length", "Tfr", "other") )
  expect_named(my2_df.mspct[[2]], c("w.length", "Tfr", "other") )
  expect_named(my2_df.mspct[[3]], c("w.length", "Tfr", "other") )

  # constructor methods for 'long' spct objects -----------------------------

  my_long.spct <- rbindspct(spct.l)
  my3_df.mspct <- subset2mspct(my_long.spct)

  expect_equal(paste("spct", 1:5, sep = "_"), names(my3_df.mspct))
  expect_equal(class(my3_df.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(class(my3_df.mspct[[1]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[2]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[3]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[4]])[1:2], c("filter_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[5]])[1:2], c("filter_spct", "generic_spct") )

  expect_equal(ncol(my3_df.mspct[[1]]), 2L )
  expect_equal(ncol(my3_df.mspct[[2]]), 2L )
  expect_equal(ncol(my3_df.mspct[[3]]), 2L )

  expect_named(my3_df.mspct[[1]], c("w.length", "Tfr") )
  expect_named(my3_df.mspct[[2]], c("w.length", "Tfr") )
  expect_named(my3_df.mspct[[3]], c("w.length", "Tfr") )

  # print -------------------------------------------------------------------

  if (test.print) expect_equal(print(my.mspct), my.mspct)

  # clean -------------------------------------------------------------------

  expect_equal(clean(my.mspct), my.mspct)
})

test_that("filter_mspct_attr", {

  my1.spct <- filter_spct(w.length = 400:410, Tfr = 0.5, Tfr.type = "internal")
  setWhatMeasured(my1.spct, "first spectrum")
  setWhenMeasured(my1.spct, lubridate::ymd("2018-03-03", tz = "UTC"))
  setWhereMeasured(my1.spct, data.frame(lat = -30, lon = +80))
  my2.spct <- filter_spct(w.length = 400:410, Tfr = 0.5, Tfr.type = "internal")
  setWhatMeasured(my2.spct, "second spectrum")
  setWhenMeasured(my2.spct, lubridate::ymd_hm("2018-03-03 12:30", tz = "UTC"))
  setWhereMeasured(my2.spct, data.frame(lat = 5, lon = 20))

  spct.l <- list(my1.spct, my2.spct)
  my.mspct <- filter_mspct(spct.l)

  expect_equal(paste("spct", seq_len(length(spct.l)), sep = "_"), names(my.mspct))

  expect_equal(class(my.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(attr(my.mspct, "mspct.version", exact = TRUE), 3)
  expect_equal(ncol(my.mspct), 1)
  expect_equal(nrow(my.mspct), length(spct.l))
  expect_equal(attr(my.mspct, "mspct.byrow", exact = TRUE), FALSE)

  named_spct.l <- list(one = my1.spct,
                       two = my2.spct)

  my_named.mspct <- filter_mspct(named_spct.l)

  expect_equal(names(named_spct.l), names(my_named.mspct))
  expect_equal(class(my_named.mspct)[1:2], c("filter_mspct", "generic_mspct") )
  expect_equal(attr(my_named.mspct, "mspct.version", exact = TRUE), 3)
  expect_equal(ncol(my_named.mspct), 1)
  expect_equal(nrow(my.mspct), length(named_spct.l))
  expect_equal(attr(my_named.mspct, "mspct.byrow", exact = TRUE), FALSE)

  expect_equal(length(my.mspct), 2)
  expect_equal(length(my_named.mspct), 2)

  expect_true(is.data.frame(transmittance(my.mspct)))
  expect_true(is.data.frame(transmittance(my_named.mspct)))

  expect_equal(transmittance(my.mspct)[["spct.idx"]],
               factor(paste("spct", 1:2, sep = "_")))
  expect_equal(levels(transmittance(my_named.mspct)[["spct.idx"]]),
               c("one", "two"))

  # transmittance -------------------------------------------------------------------

  expect_equal(
    names(
      transmittance(my.mspct,
            list(waveband(c(400,405), wb.name = "A"),
                 waveband(c(405,410), wb.name = "B")),
            attr2tb = c(when.measured = "when",
                        what.measured = "what",
                        lat = "latitude",
                        lon = "longitude")
      )), c("spct.idx", "Tfr(wl)_A", "Tfr(wl)_B",
            "when", "what", "latitude", "longitude")
  )

  expect_equal(
    names(
      transmittance(my.mspct,
            list(waveband(c(400,405), wb.name = "A"),
                 waveband(c(405,410), wb.name = "B")),
            attr2tb = c(when.measured = "when",
                        what.measured = "what",
                        lat = "latitude",
                        lon = "longitude"),
            idx = "spectrum"
      )), c("spectrum", "Tfr(wl)_A", "Tfr(wl)_B",
            "when", "what", "latitude", "longitude")
  )

  expect_equal(
    transmittance(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
  )
  expect_equal(
    transmittance(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["what"]], c("first spectrum", "second spectrum")
  )
  expect_equal(
    transmittance(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["latitude"]], c(-30, 5)
  )
  expect_equal(
    transmittance(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["longitude"]], c(80, 20)
  )

  expect_equal(
    transmittance(my_named.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
  )
  expect_equal(
    transmittance(my_named.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["what"]], c("first spectrum", "second spectrum")
  )
  expect_equal(
    transmittance(my_named.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["latitude"]], c(-30, 5)
  )
  expect_equal(
    transmittance(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["longitude"]], c(80, 20)
  )

  # relative

  expect_equal(
    transmittance(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
  )
  expect_equal(
    transmittance(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["what"]], c("first spectrum", "second spectrum")
  )
  expect_equal(
    transmittance(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["latitude"]], c(-30, 5)
  )
  expect_equal(
    transmittance(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["longitude"]], c(80, 20)
  )

  expect_equal(
    transmittance(my_named.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
  )
  expect_equal(
    transmittance(my_named.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["what"]], c("first spectrum", "second spectrum")
  )
  expect_equal(
    transmittance(my_named.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["latitude"]], c(-30, 5)
  )
  expect_equal(
    transmittance(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["longitude"]], c(80, 20)
  )


# absorptance -------------------------------------------------------------------

expect_equal(
  names(
    absorptance(my.mspct,
                  list(waveband(c(400,405), wb.name = "A"),
                       waveband(c(405,410), wb.name = "B")),
                  attr2tb = c(when.measured = "when",
                              what.measured = "what",
                              lat = "latitude",
                              lon = "longitude")
    )), c("spct.idx", "Afr(wl)_A", "Afr(wl)_B",
          "when", "what", "latitude", "longitude")
)

expect_equal(
  names(
    absorptance(my.mspct,
                  list(waveband(c(400,405), wb.name = "A"),
                       waveband(c(405,410), wb.name = "B")),
                  attr2tb = c(when.measured = "when",
                              what.measured = "what",
                              lat = "latitude",
                              lon = "longitude"),
                  idx = "spectrum"
    )), c("spectrum", "Afr(wl)_A", "Afr(wl)_B",
          "when", "what", "latitude", "longitude")
)

expect_equal(
  absorptance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
)
expect_equal(
  absorptance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["what"]], c("first spectrum", "second spectrum")
)
expect_equal(
  absorptance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["latitude"]], c(-30, 5)
)
expect_equal(
  absorptance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["longitude"]], c(80, 20)
)

expect_equal(
  absorptance(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
)
expect_equal(
  absorptance(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["what"]], c("first spectrum", "second spectrum")
)
expect_equal(
  absorptance(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["latitude"]], c(-30, 5)
)
expect_equal(
  absorptance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["longitude"]], c(80, 20)
)

# relative

expect_equal(
  absorptance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
)
expect_equal(
  absorptance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["what"]], c("first spectrum", "second spectrum")
)
expect_equal(
  absorptance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["latitude"]], c(-30, 5)
)
expect_equal(
  absorptance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["longitude"]], c(80, 20)
)

expect_equal(
  absorptance(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
)
expect_equal(
  absorptance(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["what"]], c("first spectrum", "second spectrum")
)
expect_equal(
  absorptance(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["latitude"]], c(-30, 5)
)
expect_equal(
  absorptance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["longitude"]], c(80, 20)
)


# absorbance -------------------------------------------------------------------

expect_equal(
  names(
    absorbance(my.mspct,
                  list(waveband(c(400,405), wb.name = "A"),
                       waveband(c(405,410), wb.name = "B")),
                  attr2tb = c(when.measured = "when",
                              what.measured = "what",
                              lat = "latitude",
                              lon = "longitude")
    )), c("spct.idx", "A(wl)_A", "A(wl)_B",
          "when", "what", "latitude", "longitude")
)

expect_equal(
  names(
    absorbance(my.mspct,
                  list(waveband(c(400,405), wb.name = "A"),
                       waveband(c(405,410), wb.name = "B")),
                  attr2tb = c(when.measured = "when",
                              what.measured = "what",
                              lat = "latitude",
                              lon = "longitude"),
                  idx = "spectrum"
    )), c("spectrum", "A(wl)_A", "A(wl)_B",
          "when", "what", "latitude", "longitude")
)

expect_equal(
  absorbance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
)
expect_equal(
  absorbance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["what"]], c("first spectrum", "second spectrum")
)
expect_equal(
  absorbance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["latitude"]], c(-30, 5)
)
expect_equal(
  absorbance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["longitude"]], c(80, 20)
)

expect_equal(
  absorbance(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
)
expect_equal(
  absorbance(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["what"]], c("first spectrum", "second spectrum")
)
expect_equal(
  absorbance(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["latitude"]], c(-30, 5)
)
expect_equal(
  absorbance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
  )[["longitude"]], c(80, 20)
)

# relative

expect_equal(
  absorbance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
)
expect_equal(
  absorbance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["what"]], c("first spectrum", "second spectrum")
)
expect_equal(
  absorbance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["latitude"]], c(-30, 5)
)
expect_equal(
  absorbance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["longitude"]], c(80, 20)
)

expect_equal(
  absorbance(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
)
expect_equal(
  absorbance(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["what"]], c("first spectrum", "second spectrum")
)
expect_equal(
  absorbance(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["latitude"]], c(-30, 5)
)
expect_equal(
  absorbance(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude"),
                quantity = "relative"
  )[["longitude"]], c(80, 20)
)

})


